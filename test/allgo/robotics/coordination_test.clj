(ns allgo.robotics.coordination-test
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.robotics.coordination :as coordination]
            [clojure.test :refer [deftest is testing]]))

;; Time-varying problems and multiple robots, after LaValle, *Planning
;; Algorithms*, Cambridge University Press, 2006, sections 7.1 and 7.2.
;; The scenes are made up here, small enough to check by hand.

(deftest velocity-tuning-test
  (let [path [[0.0 0.0] [10.0 0.0]]
        opts {:ds 0.5 :dt 0.25 :clearance 1.0}
        pts (coordination/samples path 0.5)]
    (testing "with nothing in the way, a sample a tick: 20 ticks of a quarter second"
      (is (= 21 (count pts)))
      (is (= [20 5.0] (peek (coordination/velocity-tuning path [] opts)))))
    (testing "a car crossing the path: the robot waits for it, three ticks, and is never within the clearance"
      (let [car (fn [t] [5.0 (- 6.0 (* 2.4 t))])
            s (coordination/velocity-tuning path [car] opts)]
        (is (= [20 5.75] (peek s)))
        (is (= 3 (count (filter (fn [[[i] [j]]] (= i j)) (partition 2 1 s)))))
        (is (every? (fn [[i t]] (> (polygon/distance (pts i) (car t)) 1.0)) s))
        (is (some (fn [[i t]] (<= (polygon/distance (pts i) (car t)) 1.0)) (coordination/velocity-tuning path [] opts)))))))

(deftest coordination-diagram-test
  (testing "two robots crossing at right angles: one lets the other through, a tick or two later than alone"
    (let [d (coordination/coordination-diagram [[0.0 0.0] [10.0 0.0]] [[5.0 -5.0] [5.0 5.0]] 0.5 1.0)
          s (coordination/coordinate d)]
      (is (= [21 21] (:size d)))
      (is (contains? (:blocked d) [10 10]))
      (is (= [0 0] (first s)))
      (is (= [20 20] (peek s)))
      (is (not-any? (:blocked d) s))
      (is (< 21 (count s) 26))
      (is (every? (fn [[[a b] [c e]]] (and (<= 0 (- c a) 1) (<= 0 (- e b) 1))) (partition 2 1 s))))))

(deftest prioritized-test
  (let [free? (fn [[i j]] (or (and (= j 0) (<= 0 i 6)) (= [2 1] [i j])))
        east [[0 0] [6 0]]
        west [[6 0] [0 0]]]
    (testing "a corridor with a bay near one end: planned first, the robot near the bay steps in to let the other by"
      (let [paths (coordination/prioritized free? [west east])]
        (is (every? some? paths))
        (is (= [] (coordination/conflicts paths)))
        (is (some #{[2 1]} (second paths)))))
    (testing "planned the other way round, the robot far from the bay cannot reach it in time: prioritized planning is incomplete"
      (let [paths (coordination/prioritized free? [east west])]
        (is (some? (first paths)))
        (is (nil? (second paths)))))))
