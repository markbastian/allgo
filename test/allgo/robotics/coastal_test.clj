(ns allgo.robotics.coastal-test
  "Coastal navigation, after Russell and Norvig, *Artificial Intelligence:
  A Modern Approach*, 4th edition, section 26.6: a hall with landmarks
  along its top wall only, and a robot whose uncertainty grows a step at
  a time away from them."
  (:require [allgo.robotics.coastal :as coastal]
            [clojure.test :refer [deftest is testing]]))

(def ^:private hall
  {:size [21 9] :walls #{} :range 1 :cap 8
   :landmarks (set (for [x (range 0 21 3)] [x 8]))})

(defn- mean-spread [path] (/ (reduce + (map second path)) (double (count path))))

(deftest shortest-is-lost
  (testing "the shortest way along the bottom of the hall, twenty steps, arrives as uncertain as the robot can be"
    (let [path (coastal/shortest hall [0 0] [20 0])]
      (is (= 21 (count path)))
      (is (= 8 (second (peek path)))))))

(deftest docking-needs-the-coast
  (testing "a dock out on the floor that needs the spread at most six -- the nearest landmark six steps away: the plan goes up to the landmarks first and comes down at the end"
    (let [path (coastal/plan hall [0 0] [20 3] :dock 6)]
      (is (<= (second (peek path)) 6))
      (is (> (count path) (count (coastal/shortest hall [0 0] [20 3]))))
      (is (some (fn [[[_ y]]] (>= y 7)) path) "it visits the coast")))
  (testing "and one needing more precision than any approach can bring -- five -- has no plan"
    (is (nil? (coastal/plan hall [0 0] [20 3] :dock 5)))))

(deftest uncertainty-as-a-cost
  (testing "weighting the spread in the cost, the robot keeps near the landmarks: a longer way, known far better along it"
    (let [coastal (coastal/plan hall [0 0] [20 0] :weight 0.5)
          shortest (coastal/shortest hall [0 0] [20 0])]
      (is (> (count coastal) (count shortest)))
      (is (< (* 2 (mean-spread coastal)) (mean-spread shortest))))))
