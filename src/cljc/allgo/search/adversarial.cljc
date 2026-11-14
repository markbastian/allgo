(ns allgo.search.adversarial
  "Adversarial search: Russell and Norvig's chapter 5, as steps.

    minimax           figure 5.3: every line of play to the end, MAX
                      taking the best of its moves and MIN the worst
    alpha-beta        figure 5.7: the same value, without the branches
                      that cannot change it
    heuristic         the chapter's cutoff: either of them cut off at `:depth`,
                      an evaluation function standing in for the rest of
                      the game; with `:order`, moves tried best-looking
                      first, and with `:table?`, a transposition table
    deepening         iterative deepening on top: depth 1, 2, 3 ..., each
                      pass trying first the move the last one chose, until
                      the game is searched to its end or the depth runs out
    expectiminimax    for games of chance: chance nodes valued at the expectation
                      of their outcomes, for games with dice
    Monte Carlo       figure 5.11: no evaluation function at all -- each
                      step one playout, down the tree by UCB1 to a new
                      node and then at random to the end of the game, the
                      result counted all the way back up
                      -- or, for a game of hidden information, ISMCTS
                      (Cowling, Powley and Whitehouse 2012): each playout
                      in a world dealt at random among those the player
                      cannot tell apart, all sharing one tree keyed by
                      actions, so that it searches what the player knows
                      rather than peeking at the true state

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
    :determinize (fn [state player rng]) -> optional, for Monte Carlo:
                 a state `player` cannot tell from `state`, its hidden
                 parts dealt at random

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
;;
;; One search for perfect and hidden information alike. What tells the
;; two apart is what names a node: the state reached (closed loop, the
;; book's figure 5.11) or the actions that reached it (open loop, which
;; information set MCTS needs). Everything else -- UCB1, expansion,
;; simulation, backing up -- is the same, and with no hidden information
;; the two agree but for transpositions.

(defn- choose
  "An action from `acts` for `state`: by `policy`'s weights if there is
  one, else uniformly."
  [rng policy state acts]
  (if policy
    (random/weighted rng (map (juxt identity #(policy state %)) acts))
    (random/pick rng acts)))

(defn- cut? [cutoff start state plies]
  (cond (nil? cutoff) false
        (number? cutoff) (>= (long plies) (long cutoff))
        :else (cutoff start state plies)))

(defn- simulate
  "SIMULATION from `start`: by the policy, through chance by sampling,
  until the game ends, no move is left, or the cutoff says stop."
  [{:keys [to-move actions result terminal? chances]} rng policy cutoff start]
  (loop [s start plies 0]
    (cond
      (terminal? s) s
      (cut? cutoff start s plies) s
      (= :chance (to-move s)) (recur (result s (random/weighted rng (chances s))) plies)
      :else (let [acts (vec (actions s))]
              (if (empty? acts) s (recur (result s (choose rng policy s acts)) (inc plies)))))))

(defn- ucb
  "UCB1, with exploration measured against how often the action was
  *available* -- under hidden information, in some sampled worlds it is
  not legal at all. With nothing hidden every child is available
  whenever its parent is visited, and this is the textbook UCB1."
  [{:keys [n w avail]} c]
  (+ (/ (double w) (long n)) (* (double c) (math/sqrt (/ (math/log (long (or avail n))) (long n))))))

(defn- child
  "The node `a` leads to from `s`, reached by `path`: `[key state]`. Keyed
  by state, the state is worked out to name it; keyed by path it is
  worked out only if the search goes there, so `state` may be nil."
  [result node-key path s a]
  (if (= node-key :state)
    (let [s' (result s a)] [s' s'])
    [(conj path a) nil]))

(defn- playout-step
  "One PLAYOUT. Deal a world the searching player cannot tell from the
  real one (just the real one, with nothing hidden); go down the tree by
  UCB1 among the actions legal in that world while each has a node; add
  one node, simulate from it, and credit the result to every node on the
  way, each from the view of the player who chose to move into it --
  counting, at every node passed, which of its children were available.
  Chance is sampled where it falls, never a node of its own."
  [{:keys [game tree root player rng c policy cutoff node-key] :as mcts}]
  (let [{:keys [to-move actions result terminal? utility evaluate chances determinize]} game
        world (if determinize (determinize root player rng) root)
        root-key (if (= node-key :state) root [])
        [ks movers seen leaf]
        (loop [s world path [] ks [root-key] movers [] seen []]
          (cond
            (terminal? s) [ks movers seen s]
            (= :chance (to-move s))
            (recur (result s (random/weighted rng (chances s))) path ks movers seen)
            :else
            (let [acts (vec (actions s))
                  kids (mapv (fn [a] (into [a] (child result node-key path s a))) acts)
                  have (filterv #(contains? tree (nth % 1)) kids)
                  fresh (filterv #(not (contains? tree (nth % 1))) kids)
                  seen (conj seen (mapv #(nth % 1) have))
                  movers (conj movers (to-move s))
                  step (fn [[a k s']] [(conj path a) (conj ks k) (or s' (result s a))])]
              (cond
                (empty? acts) [ks (pop movers) (pop seen) s]
                (seq fresh)
                (let [by-action (into {} (map (juxt first identity)) fresh)
                      [_ ks' s'] (step (by-action (choose rng policy s (mapv first fresh))))]
                  [ks' movers seen s'])
                :else
                (let [[path' ks' s'] (step (apply max-key #(ucb (tree (nth % 1)) c) have))]
                  (recur s' path' ks' movers seen))))))
        end (simulate game rng policy cutoff leaf)
        value (memoize (fn [p] (if (or (terminal? end) (nil? evaluate)) (utility end p) (evaluate end p))))
        tree (reduce (fn [t i]
                       (let [t (reduce (fn [t k] (update-in t [k :avail] (fnil inc 0))) t (seen i))]
                         (update t (ks (inc i))
                                 (fn [{:keys [n w avail]}]
                                   {:n (inc (long (or n 0)))
                                    :w (+ (double (or w 0.0)) (double (value (movers i))))
                                    :avail (or avail 1)}))))
                     (update-in tree [root-key :n] inc)
                     (range (count movers)))]
    (assoc mcts :tree tree :playouts (inc (long (:playouts mcts))) :path ks)))

(defn monte-carlo-seq
  "Every state of a Monte Carlo tree search for the player to move in
  `state`, one per playout. Unbounded: take as many as there is time for.

  Figure 5.11 as it stands is UCB1 down a tree keyed by state, random
  playouts to the end, the result counted back up. Under hidden
  information that search cheats -- it plays the state it was handed, as
  though it could see every hand -- so the game may have

    :determinize (fn [state player rng]) -> a state `player` cannot tell
                 from `state`, its hidden parts -- the cards in other
                 hands, the order of the deck -- dealt at random among
                 those consistent with what they have seen

  and then each playout first deals such a world, and the tree is keyed
  by the actions that reach a node, so that every world shares it and its
  statistics average over what the player cannot know: information set
  MCTS (Cowling, Powley and Whitehouse 2012, the single-observer form).
  An action legal in only some worlds is explored against how often it
  was available.

  Options:

    :node-key  :state, a node per state, so transpositions share their
               statistics; or :path, a node per line of play, which hidden
               information needs and which also suits games whose chance
               has too many outcomes to tell apart. Default :path if the
               game has :determinize, else :state.
    :policy    (fn [state action]) -> a positive weight; the tree's new
               moves and the playouts' moves are drawn in proportion.
               Default uniform. A little knowledge here makes playouts
               look like games rather than noise.
    :cutoff    stop a playout early and score it with the game's
               `:evaluate`: a number of plies, or (fn [start state plies])
               -> true to stop, `start` being where the playout began.
    :rng       the randomness (default seed 1); :c UCB1's constant
               (default sqrt 2).

  Each state is `{:tree :playouts :path ...}`, the tree a map from each
  node's key to `{:n :w :avail}` -- playouts through it, their total
  utility to the player who moved into it, and how often it could have
  been chosen -- and `:path` the keys the last playout went through."
  ([game state] (monte-carlo-seq game state {}))
  ([game state {:keys [rng c policy cutoff node-key]}]
   (let [node-key (or node-key (if (:determinize game) :path :state))]
     (iterate playout-step {:game game :root state :player ((:to-move game) state)
                            :tree {(if (= node-key :state) state []) {:n 0 :w 0.0}}
                            :playouts 0 :rng (or rng (random/rng 1)) :c (or c (math/sqrt 2.0))
                            :policy policy :cutoff cutoff :node-key node-key}))))

(defn monte-carlo-choice
  "The move from the root of a Monte Carlo search's state with the most
  playouts -- the most robust measure of the best, since a move is only
  played often if it kept looking good -- with each move's `{:n :q
  :avail}`: playouts, average utility to the player to move, and how
  often it was available."
  [{:keys [game tree root node-key]}]
  (let [{:keys [actions result]} game
        stats (into {} (for [a (actions root)
                             :let [{:keys [n w avail]} (tree (first (child result node-key [] root a)))]]
                         [a {:n (or n 0) :q (if (and n (pos? (long n))) (/ (double w) (long n)) 0.0)
                             :avail (or avail 0)}]))]
    {:best (when (seq stats) (key (apply max-key (comp :n val) stats)))
     :moves stats}))

(defn monte-carlo
  "The move `playouts` playouts of Monte Carlo tree search choose."
  ([game state playouts] (monte-carlo game state playouts {}))
  ([game state playouts opts]
   (monte-carlo-choice (nth (monte-carlo-seq game state opts) playouts))))
