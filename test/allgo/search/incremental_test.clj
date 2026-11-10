(ns allgo.search.incremental-test
  (:require [allgo.search :as search]
            [allgo.search.incremental :as incremental]
            [clojure.test :refer [deftest is testing]]))

;; D* Lite, after LaValle, *Planning Algorithms*, Cambridge University
;; Press, 2006, section 12.3.2: a robot crossing a 20 x 20 grid whose walls
;; it sees only when beside them.

(def ^:private n 20)
(def ^:private walls (set (concat (for [y (range 0 15)] [8 y]) (for [x (range 3 18)] [x 15]) (for [y (range 5 18)] [14 y]))))
(defn- nbrs [[x y]] (for [[dx dy] [[1 0] [-1 0] [0 1] [0 -1]] :let [c [(+ x dx) (+ y dy)]] :when (and (< -1 (c 0) n) (< -1 (c 1) n))] c))
(defn- manhattan [[a b] [c d]] (+ (abs (- a c)) (abs (- b d))))
(defn- sense [p] (for [q (nbrs p) :when (walls q) r (nbrs q)] [q r ##Inf]))

(deftest d-star-lite-test
  (let [problem {:start [0 0] :goal [19 19] :neighbors nbrs :cost (fn [_ _] 1.0) :heuristic manhattan}
        states (vec (incremental/navigate-seq problem sense))
        path (mapv :start states)
        known (fn [i] (set (map first (mapcat sense (subvec path 0 (inc i))))))
        a-star (fn [i] (search/a-star-seq {:start (path i) :goal [19 19] :neighbors #(remove (known i) (nbrs %))
                                           :cost (fn [_ _] 1) :heuristic manhattan}))]
    (testing "believing the unseen grid clear, the robot discovers the walls as it goes and still reaches the goal, never entering one"
      (is (= [19 19] (peek path)))
      (is (not-any? walls path)))
    (testing "after every repair the plan is optimal for the map as then known: its cost is A*'s from scratch"
      (is (every? (fn [i] (= (double (dec (count (search/path (search/goal-state (a-star i))))))
                             (incremental/path-cost (states i))))
                  (range (count states)))))
    (testing "repairing rather than replanning: under a tenth of the expansions of an A* search from scratch at every move"
      (is (< (* 10 (:expanded (peek states))) (reduce + (map #(count (a-star %)) (range (count states)))))))))

(deftest planned-path-test
  (testing "the plan read out from the robot's position: start to goal, a move per unit of its cost"
    (let [st (incremental/plan {:start [0 0] :goal [19 19] :neighbors nbrs :cost (fn [_ _] 1.0) :heuristic manhattan})
          p (incremental/planned-path st)]
      (is (= [0 0] (first p)))
      (is (= [19 19] (peek p)))
      (is (= (incremental/path-cost st) (double (dec (count p))))))))
