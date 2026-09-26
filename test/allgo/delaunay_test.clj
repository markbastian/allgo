(ns allgo.delaunay-test
  (:require [allgo.geometry.delaunay :as d]
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

(deftest clip-polygon-behavior
  (let [square [[0 0] [10 0] [10 10] [0 10]]]
    (testing "a polygon inside the window is unchanged in area"
      (is (< (abs (- 100.0 (area (d/clip-polygon square [[-1 -1] [11 11]])))) 1e-9)))
    (testing "clipping in half halves the area"
      (is (< (abs (- 50.0 (area (d/clip-polygon square [[0 0] [5 10]])))) 1e-9)))
    (testing "a polygon wholly outside clips away"
      (is (empty? (d/clip-polygon square [[20 20] [30 30]]))))))

(defn- naive-triangulate
  "Bowyer-Watson against a linear scan of every triangle: the
  implementation `triangulate` used to be, kept here to check the
  spatial index against.

  The index rests on one claim -- that a circumcircle containing `p` has
  a bounding box containing `p`, so looking in `p`'s grid cell finds
  every triangle the scan would have found. If that were ever wrong the
  result would still be a plausible triangulation, just not the Delaunay
  one, and the property test below would mostly still pass. This is the
  test that would not."
  [points]
  (let [sq (fn [x] (* x x))
        d2 (fn [[x1 y1] [x2 y2]] (+ (sq (- x2 x1)) (sq (- y2 y1))))
        circum (fn [[ax ay] [bx by] [cx cy]]
                 (let [dd (* 2.0 (+ (* ax (- by cy)) (* bx (- cy ay)) (* cx (- ay by))))]
                   (when-not (zero? dd)
                     (let [a2 (+ (sq ax) (sq ay)) b2 (+ (sq bx) (sq by)) c2 (+ (sq cx) (sq cy))
                           ux (/ (+ (* a2 (- by cy)) (* b2 (- cy ay)) (* c2 (- ay by))) dd)
                           uy (/ (+ (* a2 (- cx bx)) (* b2 (- ax cx)) (* c2 (- bx ax))) dd)]
                       {:center [ux uy] :radius-sq (d2 [ux uy] [ax ay])}))))
        mk (fn [p1 p2 p3] {:points [p1 p2 p3] :circle (circum p1 p2 p3)})
        inside? (fn [{:keys [circle]} p]
                  (and circle (<= (d2 (:center circle) p) (:radius-sq circle))))
        xs (map first points) ys (map second points)
        mx (/ (+ (apply min xs) (apply max xs)) 2.0)
        my (/ (+ (apply min ys) (apply max ys)) 2.0)
        dm (* 20 (max 1.0 (- (apply max xs) (apply min xs)) (- (apply max ys) (apply min ys))))
        st (mk [(- mx dm) (- my dm)] [mx (+ my dm)] [(+ mx dm) (- my dm)])
        sp (set (:points st))
        step (fn [tris p]
               (let [{bad true ok false} (group-by #(inside? % p) tris)
                     boundary (->> bad
                                   (mapcat (fn [{[a b c] :points}] [#{a b} #{b c} #{c a}]))
                                   frequencies
                                   (keep (fn [[e n]] (when (= n 1) e))))]
                 (into (vec ok)
                       (map (fn [e] (let [[p1 p2] (vec e)] (mk p1 p2 p))))
                       boundary)))]
    (remove #(some sp (:points %)) (reduce step [st] points))))

(defn- shape [tris] (set (map #(set (:points %)) tris)))

(deftest index-matches-the-scan
  (testing "the indexed triangulation is the one the linear scan gives"
    (let [g (rng 4242)]
      (doseq [n [4 9 30 120 260]]
        (let [pts (points g n)]
          (is (= (shape (naive-triangulate pts)) (shape (d/triangulate pts)))
              (str "n=" n))))))

  (testing "including the awkward inputs"
    (doseq [[label pts] [["collinear" [[0.0 0.0] [1.0 1.0] [2.0 2.0] [3.0 3.0]]]
                         ["repeated point" [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 0.0]]]
                         ["exact lattice" (vec (for [i (range 7) j (range 7)]
                                                 [(double i) (double j)]))]
                         ["one triangle" [[0.0 0.0] [1.0 0.0] [0.0 1.0]]]]]
      (is (= (shape (naive-triangulate pts)) (shape (d/triangulate pts))) label)))

  (testing "fewer than three points make no triangles, and do not throw"
    (is (empty? (d/triangulate [])))
    (is (empty? (d/triangulate [[0.0 0.0]])))
    (is (empty? (d/triangulate [[0.0 0.0] [1.0 0.0]])))))

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
