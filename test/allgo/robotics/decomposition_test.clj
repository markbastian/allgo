(ns allgo.robotics.decomposition-test
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.robotics.decomposition :as decomposition]
            [allgo.robotics.motion :as motion]
            [clojure.test :refer [deftest is testing]]))

;; Vertical cell decomposition and boustrophedon coverage, after LaValle,
;; *Planning Algorithms*, Cambridge University Press, 2006, sections 6.2.2
;; and 7.6. The checks are exact: the cells' areas against the free area.

(def ^:private square [[1.0 1.0] [3.0 1.0] [3.0 3.0] [1.0 3.0]])
(def ^:private triangle [[3.5 3.5] [4.5 3.8] [4.0 4.5]])
(def ^:private room {:bounds [[-1.0 -1.0] [5.0 5.0]] :obstacles [square triangle] :radius 0.0})

(deftest decomposition-test
  (let [cells (decomposition/vertical-decomposition room)
        doors (decomposition/adjacency cells)]
    (testing "the trapezoids tile the free space: their areas sum to the room's less the obstacles'"
      (is (< (abs (- (reduce + (map #(polygon/area (decomposition/corners %)) cells))
                     (- 36.0 4.0 (polygon/area triangle))))
             1e-9)))
    (testing "a vertical from each corner: below, beside and above the square, then round the triangle"
      (is (= 8 (count cells)))
      (is (= 18 (count doors)))
      (is (every? (fn [[[i j] p]] (= p (doors [j i]))) doors)))
    (testing "the path crosses each boundary at its door, straight within each cell, so it is free"
      (let [p (decomposition/path room [0.0 2.0] [4.0 2.0])]
        (is (= [0.0 2.0] (first p)))
        (is (= [4.0 2.0] (peek p)))
        (is (motion/path-free? room p))))
    (testing "no path into an obstacle"
      (is (nil? (decomposition/path room [0.0 2.0] [2.0 2.0]))))))

(deftest boustrophedon-test
  (testing "the sweep visits every cell once and passes within half a strip of every free point"
    (let [width 0.25
          {:keys [order path]} (decomposition/boustrophedon room [0.0 2.0] width)
          probes (for [x (range -0.95 5.0 0.1) y (range -0.95 5.0 0.1) :when (motion/free? room [x y])] [x y])
          gap (fn [q] (reduce min (map (fn [[a b]] (polygon/segment-distance q a b)) (partition 2 1 path))))]
      (is (= (range 8) (sort order)))
      (is (motion/path-free? room path))
      (is (<= (reduce max (map gap probes)) (+ (/ width 2) 1e-5))))))

(deftest walls-on-the-floor-test
  (let [world {:bounds [[0.0 0.0] [12.0 8.0]] :radius 0.0
               :obstacles [[[2.0 0.0] [2.6 0.0] [2.6 5.5] [2.0 5.5]]
                           [[4.5 2.5] [6.5 2.0] [7.0 4.0] [5.0 4.8]]
                           [[5.0 6.2] [9.0 6.2] [9.0 6.8] [5.0 6.8]]
                           [[8.5 0.8] [10.0 0.8] [10.0 4.8] [9.4 4.8] [9.4 1.4] [8.5 1.4]]
                           [[3.6 6.0] [4.2 6.0] [4.2 8.0] [3.6 8.0]]]}
        cells (decomposition/vertical-decomposition world)]
    (testing "obstacles standing on the floor or hanging from the ceiling leave no empty cells, and the rest still tile the free space"
      (is (every? (fn [c] (pos? (polygon/area (decomposition/corners c)))) cells))
      (is (< (abs (- (reduce + (map #(polygon/area (decomposition/corners %)) cells))
                     (- 96.0 (reduce + (map polygon/area (:obstacles world))))))
             1e-9)))
    (testing "where a cell's two doors lie on one side, the way between them goes by its center, clear of the corner between them"
      (is (motion/path-free? world (:path (decomposition/boustrophedon world [1.0 1.0] 0.3))))
      (is (motion/path-free? world (decomposition/path world [1.0 1.0] [11.0 7.2]))))))
