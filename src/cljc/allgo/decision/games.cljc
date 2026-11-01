(ns allgo.decision.games
  "Multiagent decision making: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 18, which
  teaches its algorithms in the text rather than in figures.

    non-cooperative  normal-form games: best responses, dominant
    games            strategies and their iterated elimination (a step a
                     strategy removed), pure Nash equilibria, Pareto
                     optimality and social welfare; mixed equilibria of
                     two-player games by support enumeration; and the
                     zero-sum game's maximin strategy and value by linear
                     programming, von Neumann's theorem
    repeated games   strategies that see the history -- always defect,
                     tit-for-tat, grim trigger -- played round by round
    sequential       alternating-offers bargaining by backward induction,
    bargaining       a step a round from the end, toward Rubinstein's
                     split
    cooperative      the Shapley value of a characteristic function
    games
    social choice    plurality, Borda count, Condorcet winners, and
                     instant runoff, a step a round
    auctions         the English auction, a step a raise; first- and
                     second-price sealed bids; and VCG for identical items

  A normal-form game is `{:actions [[...] ...] :payoff f}`: each player's
  actions, and `(f profile)` -> each player's utility, a profile being a
  vector of one action index per player. `bimatrix` makes one from a
  table of payoff pairs. A mixed strategy is a vector of probabilities
  over a player's actions."
  (:require [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.optimize :as opt]))

;; ---------------------------------------------------------------------------
;; Normal-form games

(defn bimatrix
  "A two-player game from `table`, `(table i j)` the pair `[row-utility
  column-utility]` when the row player plays its i-th action and the
  column player its j-th."
  [row-actions column-actions table]
  {:actions [(vec row-actions) (vec column-actions)]
   :payoff (fn [[i j]] (get-in table [i j]))})

(defn zero-sum
  "A two-player zero-sum game from the row player's utilities `matrix`."
  [row-actions column-actions matrix]
  (bimatrix row-actions column-actions (mapv (fn [row] (mapv (fn [u] [u (- u)]) row)) matrix)))

(defn- player-count [{:keys [actions]}] (count actions))

(defn profiles
  "Every pure-strategy profile of the game, as vectors of action indices,
  over the actions still in `:alive` (default all)."
  [{:keys [actions alive]}]
  (let [alive (or alive (mapv #(vec (range (count %))) actions))]
    (reduce (fn [ps choices] (for [p ps c choices] (conj p c))) [[]] alive)))

(defn- alive [{:keys [actions alive]} i]
  (if alive (alive i) (range (count (actions i)))))

(defn utility
  "Player `i`'s utility at `profile`."
  [{:keys [payoff]} i profile]
  (nth (payoff profile) i))

(defn best-responses
  "Player `i`'s best actions when the others play as in `profile`."
  [game i profile]
  (let [us (map (fn [a] [a (utility game i (assoc profile i a))]) (alive game i))
        best (apply max (map second us))]
    (set (keep (fn [[a u]] (when (== u best) a)) us))))

(defn pure-nash
  "The pure-strategy Nash equilibria: profiles where every player's action
  is a best response to the others'."
  [game]
  (filter (fn [p] (every? #(contains? (best-responses game % p) (p %)) (range (player-count game))))
          (profiles game)))

(defn- dominates?
  "Whether player `i`'s action `a` strictly dominates `b`: better against
  every profile of the others'."
  [game i a b]
  (every? (fn [p] (> (utility game i (assoc p i a)) (utility game i (assoc p i b))))
          (profiles game)))

(defn dominant-strategy
  "Player `i`'s strictly dominant action, or nil."
  [game i]
  (first (filter (fn [a] (every? #(or (= a %) (dominates? game i a %)) (alive game i)))
                 (alive game i))))

(defn dominant-equilibrium
  "The profile of every player's dominant action, or nil if one has none."
  [game]
  (let [ds (mapv #(dominant-strategy game %) (range (player-count game)))]
    (when (every? some? ds) ds)))

(defn elimination-seq
  "Iterated elimination of strictly dominated strategies, a step a
  strategy: each state the game with `:alive`, the actions left to each
  player, and `:removed` `[player action]`, ending when none is dominated."
  [game]
  (let [step (fn [g]
               (first (for [i (range (player-count g))
                            b (alive g i)
                            :when (some #(and (not= % b) (dominates? g i % b)) (alive g i))]
                        (-> g
                            (update-in [:alive i] (fn [as] (vec (remove #{b} as))))
                            (assoc :removed [i b])))))
        g0 (assoc game :alive (mapv #(vec (range (count %))) (:actions game)) :removed nil)]
    (take-while some? (iterate step g0))))

(defn pareto-optimal?
  "Whether no profile makes some player better off and none worse."
  [game profile]
  (let [u (vec ((:payoff game) profile))]
    (not-any? (fn [p] (let [v ((:payoff game) p)]
                        (and (every? true? (map >= v u)) (some true? (map > v u)))))
              (profiles game))))

(defn utilitarian
  "The profiles maximizing the sum of utilities, the social welfare."
  [game]
  (let [w (fn [p] (reduce + ((:payoff game) p)))
        best (apply max (map w (profiles game)))]
    (filter #(== best (w %)) (profiles game))))

;; ---------------------------------------------------------------------------
;; Mixed strategies

(defn expected-utility
  "Player `i`'s expected utility in a two-player game when the row player
  mixes by `p` and the column player by `q`."
  [game i p q]
  (reduce + (for [r (range (count p)) c (range (count q))]
              (* (p r) (q c) (utility game i [r c])))))

(defn maximin
  "The row player's maximin mixed strategy in a two-player zero-sum game
  and the game's value, by linear programming: `{:strategy :value}`. With
  every utility made positive by adding a constant, the strategy is x /
  sum x for the x >= 0 minimizing sum x subject to x . column >= 1 for
  every column; the value is 1 / sum x, less the constant. By von
  Neumann's theorem the column player's minimax gives the same value."
  [game]
  (let [m (count (alive game 0)) n (count (alive game 1))
        A (mapv (fn [i] (mapv (fn [j] (utility game 0 [i j])) (range n))) (range m))
        shift (- 1.0 (apply min (flatten A)))
        cols (mapv (fn [j] (mapv (fn [i] (+ shift (get-in A [i j]))) (range m))) (range n))
        {:keys [x]} (opt/linear-program (vec (repeat m 1.0)) cols (vec (repeat n 1.0)) {:kinds (repeat n :>=)})
        total (reduce + x)]
    {:strategy (mapv #(/ % total) x) :value (- (/ 1.0 total) shift)}))

(defn- subsets [xs]
  (reduce (fn [ss x] (concat ss (map #(conj % x) ss))) [[]] xs))

(defn- indifferent
  "The opponent's mix over `their` that makes every action in `mine`
  equally good for the player whose utility `u` gives: solve u . q = v on
  `mine`, sum q = 1. `[q v]` with q over all `n` actions, or nil."
  [u mine their n]
  (let [k (count their)
        rows (concat (map (fn [a] (conj (mapv #(u a %) their) -1.0)) mine)
                     [(conj (vec (repeat k 1.0)) 0.0)])
        b (conj (vec (repeat (count mine) 0.0)) 1.0)]
    (when-let [sol (ls/gauss (vec rows) b)]
      (let [q (reduce (fn [q [j x]] (assoc q j x)) (vec (repeat n 0.0)) (map vector their sol))]
        [q (peek sol)]))))

(defn mixed-nash
  "The Nash equilibria of a two-player game, mixed or pure, by support
  enumeration: for each pair of equal-sized supports, the mixes that make
  each player indifferent among their own support; kept when they are
  probabilities and no action outside a support does better. Every
  equilibrium of a nondegenerate game is found. `[[p q] ...]`."
  [game]
  (let [m (count (get-in game [:actions 0])) n (count (get-in game [:actions 1]))
        u0 (fn [i j] (utility game 0 [i j]))
        u1 (fn [i j] (utility game 1 [i j]))
        eps 1e-9]
    (distinct
     (for [I (rest (subsets (range m))) J (rest (subsets (range n)))
           :when (= (count I) (count J))
           :let [[q v0] (indifferent u0 I J n)
                 [p v1] (indifferent (fn [j i] (u1 i j)) J I m)]
           :when (and q p (every? #(>= % (- eps)) q) (every? #(>= % (- eps)) p))
           :let [row-u (fn [i] (reduce + (map-indexed (fn [j qj] (* qj (u0 i j))) q)))
                 col-u (fn [j] (reduce + (map-indexed (fn [i pi] (* pi (u1 i j))) p)))]
           :when (and (every? #(<= (row-u %) (+ v0 eps)) (range m))
                      (every? #(<= (col-u %) (+ v1 eps)) (range n)))]
       [(mapv #(max 0.0 %) p) (mapv #(max 0.0 %) q)]))))

;; ---------------------------------------------------------------------------
;; Repeated games

(def always-defect (fn [_mine _theirs] :defect))
(def always-cooperate (fn [_mine _theirs] :cooperate))
(def tit-for-tat
  "Cooperate first, then do what the other did last."
  (fn [_mine theirs] (or (peek theirs) :cooperate)))
(def grim
  "Cooperate until the other defects, then defect for ever."
  (fn [_mine theirs] (if (some #{:defect} theirs) :defect :cooperate)))

(def prisoners-dilemma
  "The prisoner's dilemma as payoffs to `[row column]` for each pair of
  `:cooperate` (refuse to testify) and `:defect` (testify): years in
  prison, negated."
  {[:cooperate :cooperate] [-1 -1] [:cooperate :defect] [-10 0]
   [:defect :cooperate] [0 -10] [:defect :defect] [-5 -5]})

(defn play-seq
  "Two strategies -- each (fn [own-history other's-history]) -> action --
  playing the stage game `payoffs` round after round: each state
  `{:histories [h1 h2] :totals [t1 t2]}`."
  [payoffs s1 s2]
  (iterate (fn [{[h1 h2] :histories [t1 t2] :totals}]
             (let [a1 (s1 h1 h2) a2 (s2 h2 h1)
                   [u1 u2] (payoffs [a1 a2])]
               {:histories [(conj h1 a1) (conj h2 a2)] :totals [(+ t1 u1) (+ t2 u2)]}))
           {:histories [[] []] :totals [0 0]}))

(defn tournament
  "Every pair of `strategies` (a map of name to strategy), each against
  itself too, played for `rounds`: each name's total over its games."
  [payoffs strategies rounds]
  (let [names (keys strategies)]
    (reduce (fn [scores [a b]]
              (let [[ta tb] (:totals (nth (play-seq payoffs (strategies a) (strategies b)) rounds))]
                (if (= a b)
                  (update scores a + ta)
                  (-> scores (update a + ta) (update b + tb)))))
            (zipmap names (repeat 0))
            (for [a names b names :when (<= (compare (str a) (str b)) 0)] [a b]))))

;; ---------------------------------------------------------------------------
;; Bargaining

(defn bargaining-seq
  "Alternating-offers bargaining over a pie of 1 by backward induction,
  with discount factors `d1` and `d2`: the share the player who proposes
  keeps when there are 1, 2, 3 ... rounds left. In the last round the
  proposer keeps everything; earlier, it offers the other exactly what
  the other would keep as proposer next round, discounted. It converges
  to Rubinstein's split, (1 - d2) / (1 - d1 d2) for player 1."
  [d1 d2]
  (map first
       (iterate (fn [[keep1 keep2]]
                  ;; one more round: player 1 proposes, player 2 would
                  ;; propose next and keep keep2 of it
                  [(- 1.0 (* d2 keep2)) (- 1.0 (* d1 keep1))])
                [1.0 1.0])))

(defn rubinstein
  "The subgame-perfect share of the first proposer, `[share1 share2]`."
  [d1 d2]
  (let [s (/ (- 1.0 d2) (- 1.0 (* d1 d2)))] [s (- 1.0 s)]))

;; ---------------------------------------------------------------------------
;; Cooperative games

(defn- permutations [xs]
  (if (empty? xs) [[]] (for [x xs p (permutations (remove #{x} xs))] (cons x p))))

(defn shapley
  "The Shapley value of each of `players` for the characteristic function
  `v`, (fn [coalition-set]) -> its value: each player's marginal
  contribution averaged over every order in which the coalition could
  form. Exact, by enumerating the orders."
  [players v]
  (let [orders (permutations players)
        n (count orders)]
    (reduce (fn [acc order]
              (first (reduce (fn [[acc before] p]
                               [(update acc p + (/ (- (v (conj before p)) (v before)) n)) (conj before p)])
                             [acc #{}] order)))
            (zipmap players (repeat 0))
            orders)))

;; ---------------------------------------------------------------------------
;; Social choice

(defn plurality
  "The candidates with the most first preferences among `ballots`, each a
  ranking."
  [ballots]
  (let [f (frequencies (map first ballots)) best (apply max (vals f))]
    (set (keep (fn [[c k]] (when (= k best) c)) f))))

(defn borda
  "The Borda count: each ballot gives its last choice 0 points, the next
  1, and so on up. `{candidate points}`."
  [ballots]
  (reduce (fn [scores ballot]
            (reduce (fn [s [c k]] (update s c (fnil + 0) k)) scores
                    (map vector (reverse ballot) (range))))
          {} ballots))

(defn prefer-count
  "How many ballots rank `a` above `b`."
  [ballots a b]
  (count (filter (fn [ballot] (let [rank (zipmap ballot (range))] (< (rank a) (rank b)))) ballots)))

(defn condorcet-winner
  "The candidate who beats every other head to head, or nil -- there may
  be none, since majority preference can cycle."
  [ballots]
  (let [cs (distinct (apply concat ballots))]
    (first (filter (fn [a] (every? #(or (= a %) (> (prefer-count ballots a %) (prefer-count ballots % a))) cs)) cs))))

(defn runoff-seq
  "Instant runoff, a step a round: each state the candidates still in;
  each round the one with fewest first preferences among those in drops,
  until one holds a majority. Ties drop the candidate last in sort order."
  [ballots]
  (let [firsts (fn [in] (frequencies (keep #(first (filter in %)) ballots)))]
    (take-while some?
                (iterate (fn [in]
                           (when (and in (> (count in) 1))
                             (let [f (firsts in)
                                   total (reduce + (vals f))]
                               (when-not (some #(> (* 2 %) total) (vals f))
                                 (let [loser (last (sort-by (fn [c] [(- (get f c 0)) (str c)]) in))]
                                   (disj in loser))))))
                         (set (apply concat ballots))))))

(defn runoff-winner
  "The instant-runoff winner."
  [ballots]
  (let [in (last (runoff-seq ballots))
        f (frequencies (keep #(first (filter in %)) ballots))]
    (key (apply max-key val f))))

;; ---------------------------------------------------------------------------
;; Auctions

(defn english-seq
  "An English auction, a step a raise of `increment` from `reserve`:
  bidders with private `values` stay in while the price is at most their
  value. Each state `{:price :in}`, ending when one bidder is left -- who
  wins at about the second-highest value -- or none."
  [values reserve increment]
  (take-while some?
              (iterate (fn [{:keys [price in]}]
                         (when (> (count in) 1)
                           (let [price (+ price increment)]
                             {:price price :in (set (filter #(<= price (values %)) in))})))
                       {:price reserve :in (set (filter #(<= reserve (values %)) (range (count values))))})))

(defn sealed-bid
  "A sealed-bid auction of one item: `{:winner :price}`, the highest
  bidder paying its own bid (`:first`) or the second-highest bid
  (`:second`, Vickrey's, where bidding one's true value is a dominant
  strategy). Ties go to the lower index."
  [bids kind]
  (let [order (sort-by (fn [i] [(- (bids i)) i]) (range (count bids)))
        w (first order)]
    {:winner w :price (if (= kind :first) (bids w) (if (second order) (bids (second order)) 0))}))

(defn vcg
  "The Vickrey-Clarke-Groves mechanism for `k` identical items and bidders
  who each want one: the `k` highest bids win, and each pays the harm it
  does the others -- what they would have had without it, less what they
  have -- which comes to the (k+1)-th highest bid. `{:winners :prices}`, the
  prices a map from winner to payment."
  [bids k]
  (let [order (sort-by (fn [i] [(- (bids i)) i]) (range (count bids)))
        winners (take k order)
        welfare (fn [excluded] (reduce + (map bids (take k (remove #{excluded} order)))))
        others-with (fn [w] (- (reduce + (map bids winners)) (bids w)))]
    {:winners (vec winners)
     :prices (into {} (map (fn [w] [w (- (welfare w) (others-with w))]) winners))}))

(defn auction-utility
  "Bidder `i`'s utility, its `value` less what it pays if it wins, in a
  sealed-bid auction of `kind`."
  [bids i value kind]
  (let [{:keys [winner price]} (sealed-bid bids kind)]
    (if (= winner i) (- value price) 0)))
