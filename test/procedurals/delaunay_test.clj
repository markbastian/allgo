(ns procedurals.delaunay-test
  (:require [procedurals.delaunay :as d]
            [clojure.test :refer [deftest is testing]])
  (:import (java.util Random)))

(def bounds [[0 0] [480 480]])
(def box-area (* 480.0 480.0))

(defn- rng [seed] (Random. seed))
(defn- points [^Random g n]
  (vec (repeatedly n #(vector (* 480 (.nextDouble g)) (* 480 (.nextDouble g))))))

(defn- area
  "Shoelace, unsigned."
  [poly]
  (abs (* 0.5 (reduce + (map (fn [[[x1 y1] [x2 y2]]] (- (* x1 y2) (* x2 y1)))
                             (partition 2 1 (conj (vec poly) (first poly))))))))

(deftest bounded-cells-tile-the-rectangle
  ;; Areas summing to the rectangle rules out gaps, overlaps and dropped
  ;; sites in one measurement.
  (let [g (rng 1234)]
    (doseq [n [1 5 40 120]]
      (let [pts   (points g n)
            cells (d/bounded-cells pts bounds)
            total (reduce + (map area (vals cells)))]
        (is (= n (count cells)) (str "every one of " n " sites gets a cell"))
        (is (< (abs (- total box-area)) 1e-3)
            (str "n=" n " cells fill the rectangle"))))))

(deftest bounded-cells-are-clipped-to-bounds
  (let [g     (rng 55)
        cells (d/bounded-cells (points g 40) bounds)]
    (doseq [[_ cell] cells, [x y] cell]
      (is (and (<= -1e-6 x 480.000001) (<= -1e-6 y 480.000001))))
    (doseq [[_ cell] cells]
      (is (>= (count cell) 3) "every cell is a polygon"))))

(deftest clip-polygon-behaviour
  (let [square [[0 0] [10 0] [10 10] [0 10]]]
    (testing "a polygon inside the window is unchanged in area"
      (is (< (abs (- 100.0 (area (d/clip-polygon square [[-1 -1] [11 11]])))) 1e-9)))
    (testing "clipping in half halves the area"
      (is (< (abs (- 50.0 (area (d/clip-polygon square [[0 0] [5 10]])))) 1e-9)))
    (testing "a polygon wholly outside clips away"
      (is (empty? (d/clip-polygon square [[20 20] [30 30]]))))))

(deftest triangulation-is-delaunay
  ;; The defining property: no site lies strictly inside any triangle's
  ;; circumcircle.
  (let [g    (rng 88)
        pts  (points g 40)
        tris (d/triangulate pts)]
    (is (pos? (count tris)))
    (doseq [{:keys [points circle]} tris
            :let [{:keys [center radius-sq]} circle]
            p pts
            :when (not (some #{p} points))]
      (is (> (+ (* (- (first p) (first center)) (- (first p) (first center)))
                (* (- (second p) (second center)) (- (second p) (second center))))
             (* radius-sq (- 1 1e-9)))
          (str p " lies inside the circumcircle of " points)))))

(deftest voronoi-cells-cover-every-site
  (let [g     (rng 3)
        pts   (points g 30)
        cells (d/voronoi-diagram pts)]
    (is (= (set pts) (set (keys cells))))))
