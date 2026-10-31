(ns allgo.search.problems
  "Problems to search, in the shapes `allgo.search` and
  `allgo.search.local` take, for the demos. The book's worked problems --
  Romania, the 8-puzzle, the vacuum worlds -- are in the tests that
  check its answers to them, `allgo.search-test`.

    n-queens      the complete-state formulation of section 4.1: a queen
                  per column, a state the row of each, a move one queen
                  to another row of its column, the value minus the
                  number of pairs that attack each other
    tsp           a tour of random cities, a move reversing a stretch of
                  it (2-opt)
    grid          a maze of walls and mud on a grid, for watching"
  (:require [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; n-queens

(defn attacking-pairs
  "The pairs of queens that attack each other, directly or not -- the
  book's h for n-queens."
  [rows]
  (let [n (count rows)]
    (count (for [i (range n) j (range (inc i) n)
                 :let [a (rows i) b (rows j)]
                 :when (or (= a b) (= (abs (- a b)) (- j i)))]
             [i j]))))

(defn queens
  "n-queens as a local search problem, from a random board."
  [n rng]
  (let [rows (vec (repeatedly n #(random/below rng n)))]
    {:initial rows
     :value (fn [s] (- (attacking-pairs s)))
     :goal? (fn [s] (zero? (attacking-pairs s)))
     ;; All n(n-1) successors: each queen to each other row of its column.
     :neighbors (fn [s] (for [c (range n) r (range n) :when (not= r (s c))] (assoc s c r)))
     :random-neighbor (fn [rng s]
                        (let [c (random/below rng n)
                              r (random/below rng (dec n))]
                          (assoc s c (if (>= r (s c)) (inc r) r))))
     :random-state (fn [rng] (vec (repeatedly n #(random/below rng n))))}))

(defn queens-genetic
  "n-queens for the genetic algorithm: an individual is the rows as a
  string of digits, its fitness the pairs that do not attack each
  other -- 28 for a solution of eight."
  [n]
  {:fitness (fn [rows] (- (quot (* n (dec n)) 2) (attacking-pairs rows)))
   :alphabet (vec (range n))
   :goal? (fn [rows] (zero? (attacking-pairs rows)))})

;; ---------------------------------------------------------------------------
;; The traveling salesperson

(defn tour-length [cities tour]
  (let [n (count tour)]
    (reduce + (for [k (range n)
                    :let [[x0 y0] (cities (tour k)) [x1 y1] (cities (tour (mod (inc k) n)))]]
                (math/sqrt (+ (* (- x1 x0) (- x1 x0)) (* (- y1 y0) (- y1 y0))))))))

(defn- reverse-stretch [tour i j]
  (into (subvec tour 0 i) (concat (reverse (subvec tour i j)) (subvec tour j))))

(defn tsp
  "A tour of `n` cities placed at random in the unit square, as a local
  search problem: the value is minus the tour's length, and a move
  reverses a stretch of the tour, which uncrosses two of its edges."
  [n rng]
  (let [cities (vec (repeatedly n (fn [] [(rng) (rng)])))
        tour (random/shuffle rng (vec (range n)))]
    {:cities cities
     :initial tour
     :value (fn [t] (- (tour-length cities t)))
     :neighbors (fn [t] (for [i (range n) j (range (+ i 2) (inc n))] (reverse-stretch t i j)))
     :random-neighbor (fn [rng t]
                        (let [i (random/below rng (dec n))
                              j (+ i 2 (random/below rng (- n i 1)))]
                          (reverse-stretch t i (min j n))))
     :random-state (fn [rng] (random/shuffle rng (vec (range n))))}))

;; ---------------------------------------------------------------------------
;; A grid, for watching

(defn grid
  "A `w` by `h` maze from `seed`: walls in about `density` of the cells,
  patches of mud that cost `mud` to enter, and a start and goal at
  opposite corners, kept clear. Cells are `[x y]`; `:walls` and `:mud`
  are sets. As a search, moves are to the four neighbors, and the
  heuristic is the Manhattan distance."
  [{:keys [w h seed density mud] :or {w 40 h 30 seed 1 density 0.28 mud 5}}]
  (let [rng (random/rng seed 7 0 0)
        start [1 (quot h 2)]
        goal [(- w 2) (quot h 2)]
        clear #{start goal}
        walls (set (for [x (range w) y (range h)
                         :when (and (not (clear [x y])) (< (double (rng)) density))]
                     [x y]))
        ;; Mud in blobs around a few centers.
        centers (vec (repeatedly 5 (fn [] [(random/below rng w) (random/below rng h)])))
        muddy (set (for [x (range w) y (range h) [cx cy] centers
                         :when (and (not (walls [x y]))
                                    (< (+ (math/pow (- x cx) 2) (math/pow (- y cy) 2)) 18.0))]
                     [x y]))
        inside? (fn [[x y]] (and (< -1 x w) (< -1 y h)))
        neighbors (fn [[x y]]
                    (filterv (fn [c] (and (inside? c) (not (walls c))))
                             [[(inc x) y] [x (inc y)] [(dec x) y] [x (dec y)]]))]
    {:w w :h h :walls walls :mud muddy
     :start start :goal goal
     :neighbors neighbors
     :cost (fn [_ to] (if (muddy to) mud 1))
     :heuristic (fn [[x y] [gx gy]] (+ (abs (- x gx)) (abs (- y gy))))}))
