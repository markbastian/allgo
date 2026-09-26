(ns allgo.voronoi3d-test
  (:require [allgo.geometry.voronoi3d :as v3]
            [clojure.test :refer [deftest is testing]])
  (:import (java.util Random)))

(def bounds [[0 0 0] [100 100 100]])
(def box-volume (* 100.0 100.0 100.0))

(defn- rng [seed] (Random. seed))
(defn- sites [^Random g n]
  (vec (repeatedly n #(vector (* 100 (.nextDouble g)) (* 100 (.nextDouble g)) (* 100 (.nextDouble g))))))

(defn- tet-volume [r a b c]
  (let [[ax ay az] (mapv - a r) [bx by bz] (mapv - b r) [cx cy cz] (mapv - c r)]
    (/ (abs (+ (* ax (- (* by cz) (* bz cy)))
               (- (* ay (- (* bx cz) (* bz cx))))
               (* az (- (* bx cy) (* by cx)))))
       6.0)))

(defn- cell-volume
  "Fan every face back to the site. The site lies inside its own convex
  cell, so these tetrahedra partition it."
  [site {:keys [faces]}]
  (reduce + (for [f faces [a b] (partition 2 1 (rest f))]
              (tet-volume site (first f) a b))))

(deftest cells-tile-the-bounds-exactly
  ;; The single strongest invariant available: cells summing to the volume of
  ;; the box rules out gaps, overlaps and dropped sites all at once.
  (let [g (rng 42)]
    (doseq [n [2 5 20 60]]
      (let [pts (sites g n)
            d   (v3/diagram pts bounds)
            total (reduce + (map (fn [[s c]] (cell-volume s c)) d))]
        (is (= n (count d)) (str "every one of " n " sites gets a cell"))
        (is (< (abs (- total box-volume)) 1e-4)
            (str "n=" n " cells fill the box: " total " vs " box-volume))))))

(deftest cells-stay-within-bounds
  (let [g (rng 7)
        d (v3/diagram (sites g 40) bounds)]
    (doseq [[_ {:keys [faces]}] d, face faces, [x y z] face]
      (is (and (<= -1e-6 x 100.000001) (<= -1e-6 y 100.000001) (<= -1e-6 z 100.000001))))))

(deftest every-face-is-a-polygon
  (let [g (rng 11)
        d (v3/diagram (sites g 30) bounds)]
    (doseq [[_ {:keys [faces]}] d]
      (is (>= (count faces) 4) "a bounded convex polyhedron needs four faces")
      (doseq [f faces] (is (>= (count f) 3) "each face is a polygon")))))

(deftest adjacency-is-symmetric
  ;; Neighbors are read off surviving faces rather than recorded when a cut
  ;; happens, because a later cut can remove an earlier face entirely. If that
  ;; regressed, adjacency would become order-dependent and asymmetric.
  (let [g (rng 3)]
    (doseq [n [5 20 60]]
      (let [d (v3/diagram (sites g n) bounds)]
        (doseq [[site {:keys [neighbors]}] d, other neighbors]
          (is (contains? (:neighbors (get d other)) site)
              (str site " lists " other " but not the reverse")))))))

(deftest seeding-does-not-change-the-answer
  (testing "a frame-coherence seed may only save work, never alter the result"
    (let [g   (rng 99)
          pts (sites g 40)
          idx (v3/diagram pts bounds nil)
          seeded (v3/diagram pts bounds (v3/neighbor-hints idx pts))]
      (doseq [s pts]
        (is (= (:neighbors (get idx s)) (:neighbors (get seeded s)))
            "same neighbors with and without a seed")
        (is (< (abs (- (cell-volume s (get idx s)) (cell-volume s (get seeded s)))) 1e-6)
            "same cell volume with and without a seed"))))
  (testing "a deliberately wrong seed still yields the right answer"
    (let [g   (rng 5)
          pts (sites g 25)
          plain  (v3/diagram pts bounds)
          hints  (vec (repeat (count pts) #{0 1 2}))
          seeded (v3/diagram pts bounds hints)]
      (doseq [s pts]
        (is (< (abs (- (cell-volume s (get plain s)) (cell-volume s (get seeded s)))) 1e-6))))))

(deftest neighbor-indices-agree-with-neighbor-points
  (let [g   (rng 17)
        pts (sites g 30)
        d   (v3/diagram pts bounds nil)]
    (doseq [s pts
            :let [{:keys [neighbors neighbor-idx]} (get d s)]]
      (is (= neighbors (into #{} (map #(nth pts %)) neighbor-idx)))))
  (testing "the two-argument form skips indexing entirely"
    (let [g (rng 17)
          d (v3/diagram (sites g 10) bounds)]
      (is (every? #(nil? (:neighbor-idx %)) (vals d))))))

(deftest edges-are-undirected-pairs
  (let [g (rng 23)
        d (v3/diagram (sites g 20) bounds)
        e (v3/edges d)]
    (is (every? #(= 2 (count %)) e))
    (is (= e (into #{} (for [[s {:keys [neighbors]}] d, o neighbors] #{s o}))))))

(deftest separating-plane-predicate
  (let [faces (mapv (fn [pts] {:pts pts}) (v3/box-faces [[0 0 0] [10 10 10]]))]
    (is (v3/separating-plane? faces 1.0 0.0 0.0 20.0) "plane clear of the box separates")
    (is (v3/separating-plane? faces 1.0 0.0 0.0 10.0) "plane flush with the far face still separates")
    (is (not (v3/separating-plane? faces 1.0 0.0 0.0 5.0)) "plane through the box does not")
    (is (not (v3/separating-plane? faces 1.0 0.0 0.0 -1.0)) "plane past the near face does not")))
