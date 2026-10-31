(ns allgo.search-test
  (:require [allgo.random :as random]
            [allgo.search :as s]
            [allgo.search.local :as l]
            [allgo.search.nondeterministic :as nd]
            [allgo.search.online :as o]
            [allgo.search.problems :as p]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; ---------------------------------------------------------------------------
;; The worked problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapters 3 and 4
;;
;; Romania's roads and straight-line distances are as aima-python and
;; aima-java (both MIT) have them: the two agree on every road, and
;; aima-java alone has the straight-line table.

(def roads
  "The map of Romania: each road once, with its length in kilometers."
  {[:arad :zerind] 75 [:arad :sibiu] 140 [:arad :timisoara] 118
   [:zerind :oradea] 71 [:oradea :sibiu] 151 [:timisoara :lugoj] 111
   [:lugoj :mehadia] 70 [:mehadia :drobeta] 75 [:drobeta :craiova] 120
   [:sibiu :fagaras] 99 [:sibiu :rimnicu-vilcea] 80
   [:rimnicu-vilcea :pitesti] 97 [:rimnicu-vilcea :craiova] 146
   [:craiova :pitesti] 138 [:fagaras :bucharest] 211
   [:pitesti :bucharest] 101 [:bucharest :giurgiu] 90
   [:bucharest :urziceni] 85 [:urziceni :hirsova] 98
   [:hirsova :eforie] 86 [:urziceni :vaslui] 142 [:vaslui :iasi] 92
   [:iasi :neamt] 87})

(def straight-line-to-bucharest
  "The straight-line distance from each city to Bucharest."
  {:arad 366 :bucharest 0 :craiova 160 :drobeta 242 :eforie 161
   :fagaras 176 :giurgiu 77 :hirsova 151 :iasi 226 :lugoj 244
   :mehadia 241 :neamt 234 :oradea 380 :pitesti 100 :rimnicu-vilcea 193
   :sibiu 253 :timisoara 329 :urziceni 80 :vaslui 199 :zerind 374})

(def ^:private road-length
  (reduce (fn [m [[a b] d]] (assoc m [a b] d [b a] d)) {} roads))

(def ^:private road-neighbors
  (reduce (fn [m [a b]] (update m a (fnil conj []) b)) {} (keys road-length)))

(defn romania
  "A search from `start` to `goal` on the map; with `goal` Bucharest, the
  straight-line distance is its heuristic."
  ([start] (romania start :bucharest))
  ([start goal]
   (cond-> {:start start :goal goal
            :neighbors (fn [c] (sort (road-neighbors c)))
            :cost (fn [a b] (road-length [a b]))}
     (= goal :bucharest) (assoc :heuristic (fn [c _] (straight-line-to-bucharest c))))))

;; The 8-puzzle

(def eight-goal
  "The book's goal: the blank, 0, in the top left corner."
  [0 1 2 3 4 5 6 7 8])

(defn- tile-moves [n blank]
  (let [r (quot blank n) c (mod blank n)]
    (cond-> []
      (pos? r) (conj (- blank n))
      (< r (dec n)) (conj (+ blank n))
      (pos? c) (conj (dec blank))
      (< c (dec n)) (conj (inc blank)))))

(defn puzzle-neighbors
  "The boards one slide from `board`, a vector of n*n tiles with 0 for
  the blank."
  [board]
  (let [n (long (math/sqrt (count board)))
        blank (first (keep-indexed (fn [k t] (when (zero? t) k)) board))]
    (mapv (fn [k] (assoc board blank (board k) k 0)) (tile-moves n blank))))

(defn misplaced
  "h1: the tiles not where the goal has them, the blank not counted."
  [board goal]
  (count (filter (fn [k] (and (pos? (board k)) (not= (board k) (goal k)))) (range (count board)))))

(defn manhattan
  "h2: the sum over the tiles of how many rows and columns each is from
  where the goal has it."
  [board goal]
  (let [n (long (math/sqrt (count board)))
        where (into {} (map-indexed (fn [k t] [t k]) goal))]
    (reduce + (for [k (range (count board))
                    :let [t (board k)]
                    :when (pos? t)
                    :let [g (where t)]]
                (+ (abs (- (quot k n) (quot g n))) (abs (- (mod k n) (mod g n))))))))

(defn scramble
  "A board `moves` random slides from `goal`, never undoing the slide
  before: solvable, and usually about that far from solved."
  [goal moves rng]
  (loop [b goal prev nil k 0]
    (if (= k moves)
      b
      (let [options (vec (remove #(= % prev) (puzzle-neighbors b)))
            nxt (random/pick rng options)]
        (recur nxt b (inc k))))))

(defn eight-puzzle
  "A search from `board` to the goal, with `heuristic` (default h2)."
  ([board] (eight-puzzle board manhattan))
  ([board heuristic]
   {:start board :goal eight-goal :neighbors puzzle-neighbors :heuristic heuristic}))

;; The vacuum world

;; A state is [where left-dirty? right-dirty?], `where` :left or :right.

(def vacuum-actions [:suck :left :right])

(defn vacuum-result
  "The one state an action leads to in the ordinary vacuum world."
  [[where l r] action]
  (case action
    :left [:left l r]
    :right [:right l r]
    :suck (if (= where :left) [where false r] [where l false])))

(defn erratic-results
  "The states an action may lead to in the erratic vacuum world (section
  4.3): sucking a dirty square cleans it and sometimes the other square
  too; sucking a clean one sometimes dirties it."
  [[where l r :as s] action]
  (let [here-dirty? (if (= where :left) l r)]
    (if (= action :suck)
      (if here-dirty?
        (hash-set (vacuum-result s :suck) [where false false])
        (hash-set s (if (= where :left) [where true r] [where l true])))
      #{(vacuum-result s action)})))

(defn vacuum-clean? [[_ l r]] (not (or l r)))

(def vacuum-states
  (vec (for [w [:left :right] l [true false] r [true false]] [w l r])))

;; ---------------------------------------------------------------------------

(def ^:private best-route [:arad :sibiu :rimnicu-vilcea :pitesti :bucharest])

(defn- steps-to-solve
  "Steps a search's seq takes to solve, counted without holding the seq."
  [states]
  (count (take-while (complement :solved?) states)))

(deftest romania-test
  (let [r (romania :arad)]
    (testing "the optimal searches find the book's 418 km route"
      (doseq [f [s/dijkstra s/a-star s/ida-star s/recursive-best-first s/bidirectional s/beam]]
        (is (= best-route (f r)))
        (is (= 418 (s/path-cost r (f r))))))
    (testing "greedy best-first goes by Fagaras, 450 km, as the book shows"
      (is (= [:arad :sibiu :fagaras :bucharest] (s/greedy r)))
      (is (= 450 (s/path-cost r (s/greedy r)))))
    (testing "breadth-first and iterative deepening find the fewest roads, three"
      (is (= 4 (count (s/breadth-first r)) (count (s/iterative-deepening r)))))
    (testing "weighted A* stays within its weight of the best"
      (is (<= (s/path-cost r (s/weighted-a-star r)) (* 2 418))))
    (testing "depth-limited search: cut off at 2, found at 3"
      (is (= :cutoff (:failed (last (s/depth-limited-seq (assoc r :limit 2))))))
      (is (= 4 (count (s/depth-limited (assoc r :limit 3))))))
    (testing "depth-first finds a real route, if not a good one"
      (let [route (s/depth-first r)]
        (is (= [:arad :bucharest] [(first route) (last route)]))
        (is (every? (fn [[a b]] (some #{b} ((:neighbors r) a))) (partition 2 1 route)))))
    (testing "A* expands only states whose f is at most 418"
      (let [expanded (keep :current (take-while (complement :solved?) (s/a-star-seq r)))]
        (is (every? #(<= (+ (get-in (s/goal-state (s/a-star-seq r)) [:costs %])
                            (straight-line-to-bucharest %))
                         418)
                    expanded))))))

(deftest cheapest-not-first-test
  ;; The goal is reached first by its dear road, then cheaply; uniform
  ;; cost must not stop at the first.
  (let [g {:a {:g 10 :b 1} :b {:g 1} :g {}}
        t {:start :a :goal :g :neighbors (comp keys g) :cost (fn [x y] (get-in g [x y]))}]
    (doseq [f [s/dijkstra s/a-star s/ida-star s/recursive-best-first s/bidirectional]]
      (is (= [:a :b :g] (f t))))))

(deftest goal-predicate-test
  (testing "a goal can be any state passing a test"
    (let [r (assoc (romania :arad) :goal? #{:craiova :fagaras} :goal nil :heuristic nil)]
      (is (= [:arad :sibiu :fagaras] (s/breadth-first r)))
      (is (= [:arad :sibiu :fagaras] (s/iterative-deepening r)))
      (is (= :fagaras (last (s/dijkstra r))))
      (is (= 239 (s/path-cost r (s/dijkstra r)))))))

(deftest eight-puzzle-test
  (let [boards (map (fn [k] (scramble eight-goal 30 (random/rng k))) (range 1 9))]
    (testing "every optimal search agrees on the length, and h2 does less work than h1"
      (doseq [b boards]
        (let [h1 (eight-puzzle b misplaced)
              h2 (eight-puzzle b manhattan)
              moves (count (s/breadth-first h2))]
          (doseq [path (cond-> [(s/a-star h1) (s/a-star h2) (s/ida-star h2) (s/recursive-best-first h2)
                                (s/bidirectional (dissoc h2 :heuristic))]
                         ;; About three times the work a move deeper: 1.5 s
                         ;; at twenty, so only the shallow ones.
                         (<= moves 15) (conj (s/iterative-deepening h2)))]
            (is (= moves (count path)))
            (is (= eight-goal (last path))))
          (is (<= (count (s/weighted-a-star h2)) (* 2 moves)))
          (is (<= (steps-to-solve (s/a-star-seq h2)) (steps-to-solve (s/a-star-seq h1)))))))
    (testing "every move slides one tile into the blank"
      (let [path (s/a-star (eight-puzzle (first boards)))]
        (is (every? (fn [[a b]] (some #{b} (puzzle-neighbors a))) (partition 2 1 path)))))))

(deftest grid-test
  (let [g (p/grid {:w 30 :h 20 :seed 2})
        best (s/path-cost g (s/dijkstra g))]
    (testing "the optimal searches agree on the cost through walls and mud"
      (is (= best (s/path-cost g (s/a-star g)) (s/path-cost g (s/bidirectional g))
             (s/path-cost g (s/bidirectional (dissoc g :heuristic))))))
    (testing "A* does less work than uniform cost; weighted A* less again, for a worse route"
      (is (< (steps-to-solve (s/a-star-seq g)) (steps-to-solve (s/dijkstra-seq g))))
      (is (< (steps-to-solve (s/weighted-a-star-seq g)) (steps-to-solve (s/a-star-seq g))))
      (is (<= best (s/path-cost g (s/weighted-a-star g)) (* 2 best))))
    (testing "a beam too narrow can lose the way"
      (let [routes (map #(s/beam (assoc g :beam-width %)) [1 2 4 64])]
        (is (some nil? routes))
        (is (= best (s/path-cost g (last routes))))))))

(deftest call-seq-test
  (testing "a recursion run a frame at a time: the sum of a tree"
    (let [tree [1 [2 [3] [4]] [5]]
          frame (fn [{:keys [node kids total returned] :as f}]
                  (let [total (+ (or total (first node)) (or returned 0))
                        kids (if (contains? f :kids) kids (rest node))]
                    (if (seq kids)
                      [:call {:node (first kids)} (-> f (assoc :kids (rest kids) :total total) (dissoc :returned))]
                      [:return total])))]
      (is (= 15 (:value (last (s/call-seq frame {:node tree}))))))))

;; ---------------------------------------------------------------------------
;; Local search

(defn- queens-trials [n opts]
  (for [k (range n)]
    (let [rng (random/rng (inc k))
          q (p/queens 8 rng)]
      ((:goal? q) (:state (l/hill-climbing (merge q {:rng rng} opts)))))))

(deftest hill-climbing-test
  (testing "steepest ascent solves about one 8-queens board in seven; with sideways moves, nearly all"
    (let [plain (frequencies (queens-trials 500 {}))
          sideways (frequencies (queens-trials 500 {:sideways 100}))]
      (is (< 0.09 (/ (plain true) 500.0) 0.19))
      (is (< 0.9 (/ (sideways true) 500.0)))))
  (testing "every move is uphill, and it stops on a peak"
    (let [rng (random/rng 3)
          q (p/queens 8 rng)
          states (l/hill-climbing-seq (assoc q :rng rng))
          values (map :current-value states)]
      (is (apply < (distinct values)))
      (is (every? #(<= ((:value q) %) (:current-value (last states))) ((:neighbors q) (:current (last states)))))))
  (testing "random restarts always get there"
    (doseq [k (range 1 21)]
      (let [rng (random/rng k)
            q (p/queens 8 rng)]
        (is ((:goal? q) (:current (last (l/random-restart-seq l/hill-climbing-seq (assoc q :rng rng))))))))))

(deftest annealing-test
  (testing "simulated annealing shortens a random tour to well under half"
    (let [rng (random/rng 5)
          t (p/tsp 40 rng)
          r (l/simulated-annealing (assoc t :rng rng :schedule (l/exponential-schedule 0.5 0.999 1e-4)))]
      (is (< (- (:value r)) (* 0.5 (- ((:value t) (:initial t))))))))
  (testing "and solves most 8-queens boards"
    (let [solved (for [k (range 1 41)]
                   (let [rng (random/rng k) q (p/queens 8 rng)]
                     ((:goal? q) (:state (l/simulated-annealing (assoc q :rng rng))))))]
      (is (> (count (filter true? solved)) 28)))))

(deftest beam-and-genetic-test
  (testing "local beam ends no worse than the best state it began with"
    (let [rng (random/rng 2)
          q (p/queens 8 rng)
          states (l/local-beam-seq (assoc q :rng rng :k 6))]
      (is (>= (:current-value (last states)) (:current-value (first states))))))
  (testing "the genetic algorithm breeds an 8-queens solution"
    (doseq [k (range 1 6)]
      (let [rng (random/rng k)
            g (p/queens-genetic 8)
            r (l/genetic (merge g {:rng rng :population (l/random-population rng 100 8 (:alphabet g))}))]
        (is ((:goal? g) (:state r)))
        (is (= 28 (:value r)))))))

;; ---------------------------------------------------------------------------
;; Nondeterministic, sensorless and online

(def ^:private erratic
  {:initial [:left true true] :actions (constantly vacuum-actions)
   :results erratic-results :goal? vacuum-clean?})

(deftest and-or-test
  (testing "the erratic vacuum world's plan is the book's: Suck, and if the right is still dirty, Right, Suck"
    (is (= [:suck {[:left false false] []
                   [:left false true] [:right {[:right false true] [:suck {[:right false false] []}]}]}]
           (nd/and-or erratic))))
  (testing "it reaches a clean world from every start, however the world chooses"
    (doseq [s vacuum-states
            :let [prob (assoc erratic :initial s) plan (nd/and-or prob)]
            choose [first last]]
      (is (vacuum-clean? (last (nd/execute prob plan s choose))))))
  (testing "no plan when no action can be relied on"
    (is (nil? (nd/and-or (assoc erratic :actions (constantly [:left :right])))))))

(deftest sensorless-test
  (let [sp (nd/sensorless {:initial vacuum-states :actions vacuum-actions
                           :result vacuum-result :goal? vacuum-clean?})
        path (s/breadth-first sp)]
    (testing "blind, four actions clean both squares from anywhere"
      (is (= 4 (count (nd/actions sp path))))
      (is (every? vacuum-clean? (last path))))))

(defn- maze [seed]
  (let [g (p/grid {:w 30 :h 20 :seed seed :density 0.25 :mud 1})
        dirs {:e [1 0] :n [0 1] :w [-1 0] :s [0 -1]}
        walkable (fn [[x y]] (and (< -1 x 30) (< -1 y 20) (not ((:walls g) [x y]))))]
    {:grid g
     :start (:start g) :actions (constantly [:e :n :w :s])
     :result (fn [[x y] a] (let [[dx dy] (dirs a) c [(+ x dx) (+ y dy)]] (if (walkable c) c [x y])))
     :goal? #(= % (:goal g))
     :h (fn [c] ((:heuristic g) c (:goal g)))}))

(deftest online-test
  (testing "both agents find the goal of a maze they have never seen"
    (doseq [seed [1 2 3]
            agent [o/online-dfs-agent o/lrta-star-agent]]
      (let [m (maze seed)
            end (last (o/run-seq agent m 20000))]
        (is ((:goal? m) (:state end)) [seed])
        (is (= :stop (:action end))))))
  (testing "online DFS, walled in, explores all it can reach and stops"
    (let [m (assoc (maze 1) :goal? (constantly false))
          run (o/run-seq o/online-dfs-agent m 100000)
          ;; Every cell reachable from the start: what a breadth-first
          ;; search for a goal off the grid has seen when it gives up.
          reachable (set (keys (:came-from (last (s/breadth-first-seq
                                                  (assoc (:grid m) :goal [-1 -1]))))))]
      (is (= :stop (:action (last run))))
      (is (< (:steps (last run)) 100000))
      (is (= reachable (set (map :state run)))))))
