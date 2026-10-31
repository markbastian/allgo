(ns allgo.search.games
  "Games for `allgo.search.adversarial`, in the shape it takes, for the
  demo.

    tic-tac-toe    small enough for minimax to search to the end, which
                   is how it is known that perfect play draws
    connect-four   seven columns of six; too large to search to the end
                   in a browser, so the searches are cut off and an
                   evaluation function guesses: each line of four cells
                   that one player could still complete scores for them by
                   how many of its cells they hold

  A state is `{:board :to-move ...}`, the board a vector of :x, :o and
  nil, and :x moves first. A win is worth more the sooner it comes --
  1 plus a hundredth for each empty cell left -- so that a search that
  sees two wins takes the nearer."
  (:require [clojure.math :as math]))

(defn other [p] (if (= p :x) :o :x))

(defn- win-utility [{:keys [winner board]} player]
  (if winner
    (let [v (+ 1.0 (* 0.01 (count (filter nil? board))))]
      (if (= winner player) v (- v)))
    0.0))

;; ---------------------------------------------------------------------------
;; Tic-tac-toe

(def ^:private ttt-lines
  [[0 1 2] [3 4 5] [6 7 8] [0 3 6] [1 4 7] [2 5 8] [0 4 8] [2 4 6]])

(defn- ttt-winner [board]
  (some (fn [[a b c]] (let [p (board a)] (when (and p (= p (board b) (board c))) p))) ttt-lines))

(def tic-tac-toe
  {:initial {:board (vec (repeat 9 nil)) :to-move :x}
   :to-move :to-move
   :actions (fn [{:keys [board]}] (filterv #(nil? (board %)) (range 9)))
   :result (fn [{:keys [board to-move]} k]
             (let [b (assoc board k to-move)]
               {:board b :to-move (other to-move) :winner (ttt-winner b)}))
   :terminal? (fn [{:keys [board winner]}] (or (some? winner) (every? some? board)))
   :utility win-utility
   ;; Center, then corners, then edges: the moves that sit on most lines.
   :order (fn [_ moves] (sort-by #(case % 4 0 (0 2 6 8) 1 2) moves))
   :evaluate (fn [_ _] 0.0)})

;; ---------------------------------------------------------------------------
;; Connect Four

(def columns 7)
(def rows 6)

(defn- cell [c r] (+ c (* columns r)))

(def ^:private windows
  "Every line of four cells: across, up, and both diagonals."
  (vec (for [c (range columns) r (range rows)
             [dc dr] [[1 0] [0 1] [1 1] [1 -1]]
             :let [cells (for [k (range 4)] [(+ c (* k dc)) (+ r (* k dr))])]
             :when (every? (fn [[x y]] (and (< -1 x columns) (< -1 y rows))) cells)]
         (mapv (fn [[x y]] (cell x y)) cells))))

(def ^:private windows-through
  "For each cell, the lines of four it is part of."
  (reduce (fn [m w] (reduce (fn [m k] (update m k (fnil conj []) w)) m w))
          {} windows))

(defn- four? [board k p]
  (some (fn [w] (every? #(= p (board %)) w)) (windows-through k)))

(defn- c4-evaluate
  "Lines of four still open to one player alone, weighted by the cells of
  them it holds -- 1, 4, 16 for one, two, three -- for `player` minus
  against, squashed by tanh to stay inside the utilities."
  [{:keys [board]} player]
  (let [score (reduce (fn [acc w]
                        (let [mine (count (filter #(= player (board %)) w))
                              theirs (count (filter #(= (other player) (board %)) w))]
                          (cond
                            (and (pos? mine) (zero? theirs)) (+ acc (math/pow 4.0 (dec mine)))
                            (and (pos? theirs) (zero? mine)) (- acc (math/pow 4.0 (dec theirs)))
                            :else acc)))
                      0.0 windows)]
    (* 0.9 (math/tanh (/ score 40.0)))))

(def connect-four
  {:initial {:board (vec (repeat (* columns rows) nil)) :heights (vec (repeat columns 0)) :to-move :x}
   :to-move :to-move
   :actions (fn [{:keys [heights]}] (filterv #(< (heights %) rows) (range columns)))
   :result (fn [{:keys [board heights to-move]} c]
             (let [k (cell c (heights c))
                   b (assoc board k to-move)]
               {:board b :heights (update heights c inc) :to-move (other to-move)
                :winner (when (four? b k to-move) to-move) :last c}))
   :terminal? (fn [{:keys [winner heights]}] (or (some? winner) (every? #(= rows %) heights)))
   :utility win-utility
   ;; The middle columns first: they sit on the most lines of four.
   :order (fn [_ moves] (sort-by #(abs (- % 3)) moves))
   :evaluate c4-evaluate})

(defn play
  "The state after `moves` from the start of `game`."
  [game moves]
  (reduce (:result game) (:initial game) moves))
