(ns allgo.robotics.navigation-test
  (:require [allgo.robotics.navigation :as nav]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Navigation functions, after LaValle, *Planning Algorithms*, Cambridge
;; University Press, 2006, section 8.2: a cup a potential field stalls in.

(def ^:private cup
  "A 20 by 14 room with a cup open to the left, its inside facing the
  start: walls on three sides of a pocket."
  (let [walls (set (concat (for [y (range 3 11)] [12 y]) (for [x (range 7 13)] [x 3]) (for [x (range 7 13)] [x 10])))]
    (fn [[x y]] (and (< -1 x 20) (< -1 y 14) (not (walls [x y]))))))

(deftest wavefront-test
  (let [nf (nav/wavefront cup [18 6])]
    (testing "NF1: the goal 0, its neighbors 1, each cell its distance in steps -- every free cell reached"
      (is (zero? (nf [18 6])))
      (is (= 1 (nf [17 6])))
      (is (= (count (for [x (range 20) y (range 14) :when (cup [x y])] [x y])) (count nf))))
    (testing "no local minima: every free cell but the goal has a neighbor lower than itself"
      (is (every? (fn [[c v]] (or (zero? v) (some #(< (get nf % ##Inf) v) (for [[dx dy] [[1 0] [-1 0] [0 1] [0 -1]]] [(+ (c 0) dx) (+ (c 1) dy)])))) nf)))
    (testing "so descending from deep in the cup, where a potential field stalls, reaches the goal"
      (let [path (nav/descend nf cup [10 6])]
        (is (= [18 6] (peek path)))
        (is (every? cup path))))))

(deftest optimal-navigation-test
  (testing "with diagonal moves: in open space each value is the octile distance -- a long axis and root-2 diagonals"
    (let [open (fn [[x y]] (and (< -1 x 10) (< -1 y 10)))
          nf (nav/optimal-navigation open [0 0])]
      (is (< (abs (- (+ 3.0 (* 4 (math/sqrt 2.0))) (nf [7 4]))) 1e-9))
      (is (= 7 (dec (count (nav/descend nf open [7 4]))))))))
