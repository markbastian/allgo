(ns allgo.search.adversarial
  "Adversarial search: Russell and Norvig's chapter 5, as steps.

    minimax           figure 5.3: every line of play to the end, MAX
                      taking the best of its moves and MIN the worst
    alpha-beta        figure 5.7: the same value, without the branches
                      that cannot change it
    heuristic         section 5.3: either of them cut off at `:depth`,
                      an evaluation function standing in for the rest of
                      the game; with `:order`, moves tried best-looking
                      first, and with `:table?`, a transposition table
    deepening         iterative deepening on top: depth 1, 2, 3 ..., each
                      pass trying first the move the last one chose, until
                      the game is searched to its end or the depth runs out
    expectiminimax    section 5.5: chance nodes valued at the expectation
                      of their outcomes, for games with dice
    Monte Carlo       figure 5.11: no evaluation function at all -- each
                      step one playout, down the tree by UCB1 to a new
                      node and then at random to the end of the game, the
                      result counted all the way back up

  Minimax, alpha-beta and expectiminimax are recursive; they run on
  `allgo.search/call-seq`, a step each node visited or valued, so that
  a search can be watched, or spread across animation frames, or stopped.
  Each machine state's `:context` has the count of `:nodes` visited and
  `:pruned` (moves never looked at), and at the root, the `:best` move so
  far, the `:scores` of the moves valued, and their `:bounds` -- `:upper`
  for a move alpha-beta cut short, whose score is only a bound. Monte Carlo's states are its
  tree.

  A game is a map:

    :initial     the state play begins in
    :to-move     (fn [state]) -> the player to move, or :chance
    :actions     (fn [state]) -> the moves there
    :result      (fn [state move]) -> the state it leads to
    :terminal?   (fn [state]) -> whether the game is over
    :utility     (fn [state player]) -> what the end is worth to `player`
    :evaluate    (fn [state player]) -> an estimate of the same, for a
                 state the search is cut off at; it should stay inside the
                 utilities, so that a real win always outranks a guess
    :chances     (fn [state]) -> `[[outcome probability] ...]` at a chance
                 node, each outcome passed to `:result`

  Minimax and alpha-beta assume two players and a zero sum: the value is
  the utility to the player to move at the root, and the other one
  minimizes it."
  (:require [allgo.random :as random]
            [allgo.search :as search]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Minimax, alpha-beta and expectiminimax as frames

(defn- table-lookup
  "A transposition table's answer for `state` searched `remaining` plies
  deeper: a value to return, or the window narrowed by a bound."
  [table state remaining alpha beta]
  (when-let [{:keys [value depth flag]} (get table state)]
    ;; Deep enough: a full search takes only another full search's word.
    (when (if (nil? remaining) (nil? depth) (and depth (>= (long depth) (long remaining))))
      (case flag
        :exact {:value value}
        :lower (let [a (max alpha value)] (if (>= a beta) {:value value} {:alpha a}))
        :upper (let [b (min beta value)] (if (<= b alpha) {:value value} {:beta b}))))))

(defn- game-frame
  "One node of minimax, alpha-beta or expectiminimax, as a frame. A
  decision frame tries its moves in turn, each a child frame, keeping the
  best value -- most for MAX, least for MIN -- and with pruning returns as
  soon as the value is outside the window (alpha, beta) the frames above
  it can still use. A chance frame sums its outcomes' values weighted by
  their probabilities."
  [{:keys [to-move actions result terminal? utility evaluate chances order]}
   {:keys [root-player prune? depth table?]}]
  (fn [{:keys [state ply alpha beta returned] :as frame} ctx]
    (if-not (contains? frame :moves)
      ;; First visit.
      (let [ctx (update ctx :nodes inc)
            remaining (when depth (- (long depth) (long ply)))
            ;; Not at the root, which must find its best move itself.
            hit (when (and table? (pos? (long ply))) (table-lookup (:table ctx) state remaining alpha beta))
            player (to-move state)]
        (cond
          (terminal? state) [:return (utility state root-player) nil ctx]
          (and depth (>= (long ply) (long depth)))
          [:return (evaluate state root-player) nil (update ctx :cutoffs inc)]
          (contains? hit :value) [:return (:value hit) nil ctx]

          (= :chance player)
          (let [outs (chances state)]
            [:call {:state (result state (ffirst outs)) :ply (inc (long ply)) :alpha ##-Inf :beta ##Inf}
             (assoc frame :moves (rest outs) :move (first outs) :chance? true :v 0.0) ctx])

          :else
          (let [alpha (or (:alpha hit) alpha) beta (or (:beta hit) beta)
                moves (cond->> (actions state) order (order state))
                max? (= player root-player)
                frame (assoc frame :alpha alpha :beta beta :alpha0 alpha :beta0 beta
                             :max? max? :v (if max? ##-Inf ##Inf) :moves moves)]
            [:call {:state (result state (first moves)) :ply (inc (long ply)) :alpha alpha :beta beta}
             (assoc frame :moves (rest moves) :move (first moves)) ctx])))

      (let [{:keys [moves move v max? chance?]} frame
            r returned
            frame (dissoc frame :returned)]
        (if chance?
          (let [v (+ v (* (second move) r))]
            (if (seq moves)
              (let [o (first moves)]
                [:call {:state (result state (first o)) :ply (inc (long ply)) :alpha ##-Inf :beta ##Inf}
                 (assoc frame :moves (rest moves) :move o :v v) ctx])
              [:return v nil ctx]))
          (let [better? (if max? (> r v) (< r v))
                v (if better? r v)
                ctx (cond-> ctx
                      (zero? (long ply)) (update :scores assoc move r)
                      ;; A move that could not beat the best so far was cut
                      ;; short: its score is only an upper bound.
                      (zero? (long ply)) (update :bounds assoc move
                                                 (if (and prune? (<= r (:alpha frame))) :upper :exact))
                      (and (zero? (long ply)) better?) (assoc :best move :value r))
                alpha (if max? (max alpha v) alpha)
                beta (if max? beta (min beta v))
                cut? (and prune? (if max? (>= v beta) (<= v alpha)))
                store (fn [ctx]
                        (if table?
                          (assoc-in ctx [:table state]
                                    {:value v :depth (when depth (- (long depth) (long ply)))
                                     :flag (cond (<= v (:alpha0 frame)) :upper
                                                 (>= v (:beta0 frame)) :lower
                                                 :else :exact)})
                          ctx))]
            (cond
              cut? [:return v nil (store (update ctx :pruned + (count moves)))]
              (seq moves)
              [:call {:state (result state (first moves)) :ply (inc (long ply)) :alpha alpha :beta beta}
               (assoc frame :moves (rest moves) :move (first moves) :v v :alpha alpha :beta beta) ctx]
              :else [:return v nil (store ctx)])))))))

(defn- game-seq [game state opts]
  (let [opts (assoc opts :root-player ((:to-move game) state))]
    (search/call-seq (game-frame game opts)
                     {:state state :ply 0 :alpha ##-Inf :beta ##Inf}
                     {:nodes 0 :pruned 0 :cutoffs 0 :scores {} :best nil
                      :table (or (:table opts) {})})))

(defn minimax-seq
  "Every machine state of a minimax search from `state`: with `:depth`, cut
  off there and valued by the game's `:evaluate`; at chance nodes, the
  expectation -- expectiminimax. The last state's `:value` is the value
  of `state` to the player to move, and its `:context` `:best` the move."
  ([game state] (minimax-seq game state {}))
  ([game state opts] (game-seq (cond-> game (not (:order? opts)) (dissoc :order))
                               state (assoc opts :prune? false))))

(defn alpha-beta-seq
  "Every machine state of an alpha-beta search from `state`. Options as
  for `minimax-seq`, and `:order?` to try moves in the game's `:order`,
  `:table?` to keep a transposition table, carried over from `:table`."
  ([game state] (alpha-beta-seq game state {}))
  ([game state opts]
   (game-seq (cond-> game (not (:order? opts)) (dissoc :order))
             state (assoc opts :prune? true))))

(defn- decision [states]
  (let [{:keys [value context]} (last states)]
    (assoc (select-keys context [:best :scores :bounds :nodes :pruned :cutoffs]) :value value)))

(defn minimax
  "`{:best :value :scores :nodes}` for the player to move in `state`."
  ([game state] (minimax game state {}))
  ([game state opts] (decision (minimax-seq game state opts))))

(defn alpha-beta
  "`{:best :value :scores :nodes :pruned}` for the player to move in
  `state`; the same best value as `minimax`, for less work."
  ([game state] (alpha-beta game state {}))
  ([game state opts] (decision (alpha-beta-seq game state opts))))

;; ---------------------------------------------------------------------------
;; Iterative deepening

(defn deepening-seq
  "Alpha-beta searches from `state` to depths 1, 2, ... `:max-depth`
  (default 64), each pass trying first the move the pass before chose and
  sharing its transposition table, until a pass meets no cutoff -- it
  searched to the end of the game -- or the depth runs out. The machine
  states of every pass, one after another, each with the `:depth` of its
  pass; a pass's last state, with `:value`, holds its decision."
  [{:keys [order] :as game} state {:keys [max-depth table?] :or {max-depth 64} :as opts}]
  (letfn [(pass [d best table]
            (lazy-seq
             (let [first-best (fn [s moves]
                                (let [moves (if order (order s moves) moves)]
                                  (if (and best (= s state) (some #{best} moves))
                                    (cons best (remove #{best} moves))
                                    moves)))
                   states (map #(assoc % :depth d)
                               (alpha-beta-seq (assoc game :order first-best) state
                                               (assoc opts :depth d :order? true :table table)))
                   end (last states)]
               (if (or (>= d (long max-depth)) (zero? (long (:cutoffs (:context end)))))
                 states
                 (concat states (pass (inc d) (:best (:context end))
                                      (when table? (:table (:context end)))))))))]
    (pass 1 nil {})))

;; ---------------------------------------------------------------------------
;; Monte Carlo tree search

(defn- playout-to-end
  "Moves at random from `state` to the end of the game: SIMULATION."
  [{:keys [to-move actions result terminal? chances]} rng state]
  (loop [s state]
    (cond
      (terminal? s) s
      (= :chance (to-move s)) (let [r (double (rng))
                                    outs (chances s)
                                    o (loop [[[o p] & more] outs acc 0.0]
                                        (if (or (empty? more) (< r (+ acc p))) o (recur more (+ acc p))))]
                                (recur (result s o)))
      :else (recur (result s (random/pick rng (vec (actions s))))))))

(defn- ucb
  "UCB1: the average utility of `child` to the player who moved into it,
  plus `c` times how little it has been tried beside its parent."
  [tree parent-n child c]
  (let [{:keys [n w]} (tree child)]
    (+ (/ w n) (* c (math/sqrt (/ (math/log parent-n) n))))))

(defn- playout-step
  "One PLAYOUT: down the tree by UCB1 while every move from a node has a
  node of its own; then one new node, a random game from it to the end,
  and the result added to every node on the way, each from the view of
  the player who chose to move into it."
  [{:keys [game tree root rng c] :as mcts}]
  (let [{:keys [to-move actions result terminal? utility chances]} game
        [path leaf]
        (loop [s root path [root]]
          (cond
            (terminal? s) [path s]
            (= :chance (to-move s))
            ;; Chance nodes are passed through by sampling, not choosing.
            (let [r (double (rng))
                  o (loop [[[o p] & more] (chances s) acc 0.0]
                      (if (or (empty? more) (< r (+ acc p))) o (recur more (+ acc p))))
                  nxt (result s o)]
              (if (contains? tree nxt) (recur nxt (conj path nxt)) [(conj path nxt) nxt]))
            :else
            (let [kids (map #(result s %) (actions s))
                  fresh (vec (remove #(contains? tree %) kids))]
              (if (seq fresh)
                (let [k (random/pick rng fresh)] [(conj path k) k])
                (let [n (:n (tree s))
                      best (apply max-key #(ucb tree n % c) kids)]
                  (recur best (conj path best)))))))
        end (playout-to-end game rng leaf)
        tree (first
              (reduce (fn [[t last-mover] [parent node]]
                        ;; A chance node's outcome is credited to whoever
                        ;; moved into the chance node.
                        (let [mover (if (= :chance (to-move parent)) last-mover (to-move parent))
                              u (if (some? mover) (utility end mover) 0.0)]
                          [(update t node (fn [{:keys [n w] :or {n 0 w 0.0}}] {:n (inc n) :w (+ w u)}))
                           mover]))
                      [(update tree root (fn [{:keys [n] :or {n 0}}] {:n (inc n) :w 0.0})) nil]
                      (partition 2 1 path)))]
    (assoc mcts :tree tree :playouts (inc (long (:playouts mcts))) :path path)))

(defn monte-carlo-seq
  "Every state of a Monte Carlo tree search from `state`, one per playout:
  `{:tree :playouts :path}`, the tree a map from each state in it to its
  playouts `:n` and their total utility `:w` to the player who moved into
  it. `:c` is UCB1's exploration constant (default sqrt 2), `:rng` the
  randomness (default seed 1). Unbounded: take as many as there is time
  for."
  ([game state] (monte-carlo-seq game state {}))
  ([game state {:keys [rng c]}]
   (iterate playout-step {:game game :root state :tree {state {:n 0 :w 0.0}} :playouts 0
                          :rng (or rng (random/rng 1)) :c (or c (math/sqrt 2.0))})))

(defn monte-carlo-choice
  "The move from the root of a Monte Carlo search's state with the most
  playouts -- the most robust measure of the best, since a move is only
  played often if it kept looking good -- with each move's
  `{:n :q}`, playouts and average utility to the player to move."
  [{:keys [game tree root]}]
  (let [{:keys [actions result]} game
        stats (into {} (for [a (actions root) :let [{:keys [n w]} (tree (result root a))]]
                         [a {:n (or n 0) :q (if (and n (pos? n)) (/ w n) 0.0)}]))]
    {:best (when (seq stats) (key (apply max-key (comp :n val) stats)))
     :moves stats}))

(defn monte-carlo
  "The move `playouts` playouts of Monte Carlo tree search choose."
  ([game state playouts] (monte-carlo game state playouts {}))
  ([game state playouts opts]
   (monte-carlo-choice (nth (monte-carlo-seq game state opts) playouts))))
