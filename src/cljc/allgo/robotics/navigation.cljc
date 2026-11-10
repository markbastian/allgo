(ns allgo.robotics.navigation
  "Navigation functions: LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, section 8.2 -- a plan that is not a path but a
  feedback: from wherever the robot finds itself, which way to go.

  A navigation function assigns each free cell a value that is zero at
  the goal and has, everywhere else, a neighbor with a smaller one; so
  stepping always downhill reaches the goal, from every cell that can
  reach it at all, and there are no local minima to stall in, as a
  potential field's has. The wavefront (NF1, section 8.2.3) is the
  simplest: a breadth-first search outward from the goal, each cell's
  value its number of steps from it. With diagonal moves weighed by root
  two, Dijkstra's algorithm gives the optimal one, each value the
  shortest path's length.

  Cells are `[i j]`; `free?` says which the robot may enter."
  (:require [allgo.search :as search]
            [clojure.math :as math]))

(def ^:private four [[1 0] [-1 0] [0 1] [0 -1]])
(def ^:private eight [[1 0] [-1 0] [0 1] [0 -1] [1 1] [1 -1] [-1 1] [-1 -1]])

(defn- neighbors [free? moves [i j]]
  (for [[di dj] moves :let [c [(+ i di) (+ j dj)]] :when (free? c)] c))

(defn- wavefront-step
  "The next wave: the unlabeled free cells beside the last, one more."
  [free? {:keys [values wave]}]
  (let [d (inc (values (first wave)))
        next-wave (vec (distinct (for [c wave n (neighbors free? four c) :when (not (values n))] n)))]
    {:values (into values (map vector next-wave (repeat d))) :wave next-wave}))

(defn wavefront-seq
  "NF1, a wave a step: from the `goal` cell, each step labels the free
  cells next to the last wave -- four ways -- one more than it. Each state
  `{:values :wave}`, until the wave is empty: every cell that can reach
  the goal labeled with its distance in steps."
  [free? goal]
  (take-while (comp seq :wave) (iterate (partial wavefront-step free?) {:values {goal 0} :wave [goal]})))

(defn wavefront
  "The NF1 navigation function: a map from each cell that can reach
  `goal` to its distance in steps."
  [free? goal]
  (:values (last (wavefront-seq free? goal))))

(defn optimal-navigation
  "The optimal navigation function with diagonal moves: each cell's
  shortest-path length to `goal`, moves along an axis costing 1 and along
  a diagonal root 2 -- Dijkstra's search outward from the goal,
  `allgo.search`'s."
  [free? goal]
  (let [final (last (search/dijkstra-seq {:start goal :goal ::none
                                          :neighbors #(neighbors free? eight %)
                                          :cost (fn [[a b] [c d]] (if (and (not= a c) (not= b d)) (math/sqrt 2.0) 1.0))}))]
    (:costs final)))

(defn descend
  "The feedback plan run from `cell`: a step at a time to the neighbor of
  least value, eight ways, until the goal (value 0): the cells visited.
  Nil if `cell` cannot reach the goal."
  [nf free? cell]
  (when (contains? nf cell)
    (loop [c cell path [cell]]
      (if (zero? (nf c))
        path
        (let [n (apply min-key #(get nf % ##Inf) (neighbors free? eight c))]
          (recur n (conj path n)))))))
