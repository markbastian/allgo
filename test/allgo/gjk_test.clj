(ns allgo.gjk-test
  (:require [allgo.geometry.gjk :as gjk]
            [clojure.test :refer [deftest is testing]])
  (:import (java.util Random)))

(defn- rng [seed] (Random. seed))
(defn- r [^Random g lo hi] (+ lo (* (.nextDouble g) (- hi lo))))
(defn- mag [[x y z]] (Math/sqrt (+ (* x x) (* y y) (* z z))))
(defn- sub [a b] (mapv - a b))
(defn- about
  ([a b] (about a b 1e-6))
  ([a b tol] (< (Math/abs (double (- a b))) tol)))

(deftest distance-matches-analytic
  (testing "axis-separated boxes"
    (is (about 2.0 (:distance (gjk/distance (gjk/box [0 0 0] [1 1 1])
                                            (gjk/box [3 0 0] [4 1 1]))))))
  (testing "corner-to-corner boxes"
    (is (about (* 2 (Math/sqrt 3))
               (:distance (gjk/distance (gjk/box [0 0 0] [1 1 1])
                                        (gjk/box [3 3 3] [4 4 4]))))))
  (testing "spheres, distance is centre separation less both radii"
    (is (about 2.0 (:distance (gjk/distance (gjk/sphere [0 0 0] 1) (gjk/sphere [5 0 0] 2)))))
    (is (about 3.0 (:distance (gjk/distance (gjk/sphere [0 0 0] 1) (gjk/sphere [3 4 0] 1))))))
  (testing "a point against a box, nearest feature is a vertex"
    (is (about (Math/sqrt 3)
               (:distance (gjk/distance (gjk/point-cloud [[-1 -1 -1]])
                                        (gjk/box [0 0 0] [1 1 1]))))))
  (testing "a point cloud is the hull of its points"
    (is (about 3.0 (:distance (gjk/distance
                               (gjk/point-cloud (for [x [0 2] y [0 2] z [0 2]] [x y z]))
                               (gjk/box [5 0 0] [6 2 2]))))))
  (testing "shapes in contact are at distance zero"
    (is (about 0.0 (:distance (gjk/distance (gjk/box [0 0 0] [1 1 1])
                                            (gjk/box [1 0 0] [2 1 1])))))))

(deftest witness-points-lie-on-the-shapes
  (let [{:keys [on-a on-b distance]}
        (gjk/distance (gjk/box [0 0 0] [1 1 1]) (gjk/box [3 3 3] [4 4 4]))]
    (is (about distance (mag (sub on-b on-a)))
        "the witnesses are exactly `distance` apart")
    (is (every? #(about 1.0 %) on-a) "nearest point of A is its far corner")
    (is (every? #(about 3.0 %) on-b) "nearest point of B is its near corner")))

(deftest thin-slab-is-not-reported-as-overlapping
  ;; Boxes almost touching on one axis while overlapping broadly on the
  ;; others make the Minkowski difference a thin slab, which drives GJK's
  ;; simplex near-coplanar. Enclosure must not be inferred from a flat
  ;; tetrahedron: every same-side test there is ambiguous.
  (doseq [gap [1e-1 1e-2 1e-3 1e-4]]
    (let [res (gjk/distance (gjk/box [0 0 0] [2 2 2])
                            (gjk/box [(+ 2 gap) -5 -5] [10 5 5]))]
      (is (not (:overlap? res)) (str "gap " gap " must not read as overlap"))
      (is (about gap (:distance res) 1e-9) (str "gap " gap " measured exactly")))))

(deftest overlap-detection
  (is (gjk/intersects? (gjk/box [0 0 0] [2 2 2]) (gjk/box [1 1 1] [3 3 3])))
  (is (gjk/intersects? (gjk/sphere [0 0 0] 2) (gjk/sphere [1 0 0] 2)))
  (is (not (gjk/intersects? (gjk/box [0 0 0] [1 1 1]) (gjk/box [3 3 3] [4 4 4]))))
  (is (not (gjk/intersects? (gjk/sphere [0 0 0] 1) (gjk/sphere [5 0 0] 1)))))

(deftest penetration-depth-matches-analytic
  (testing "boxes overlapping on each axis in turn"
    (is (about 1.0 (:depth (gjk/penetration (gjk/box [0 0 0] [2 2 2]) (gjk/box [1 0 0] [3 2 2])))))
    (is (about 0.5 (:depth (gjk/penetration (gjk/box [0 0 0] [2 2 2]) (gjk/box [1.5 0 0] [3 2 2])))))
    (is (about 1.0 (:depth (gjk/penetration (gjk/box [0 0 0] [4 4 4]) (gjk/box [0 3 0] [4 7 4])))))
    (is (about 1.0 (:depth (gjk/penetration (gjk/box [0 0 0] [4 4 4]) (gjk/box [0 0 3] [4 4 7]))))))
  (testing "a box wholly inside another must travel far enough to escape it"
    (is (about 6.0 (:depth (gjk/penetration (gjk/box [0 0 0] [10 10 10])
                                            (gjk/box [4 4 4] [6 6 6]))))))
  (testing "spheres; EPA approaches a curved surface from inside, so it converges rather than terminating exactly"
    (is (about 1.0 (:depth (gjk/penetration (gjk/sphere [0 0 0] 2) (gjk/sphere [3 0 0] 2))) 1e-3))
    (is (about 0.5 (:depth (gjk/penetration (gjk/sphere [0 0 0] 2) (gjk/sphere [3 4 0] 3.5))) 1e-3))))

(deftest penetration-depth-exact-for-random-boxes
  ;; A - B of two axis-aligned boxes is itself a box, so the depth has a
  ;; closed form to check every sample against.
  (let [g (rng 20260910)]
    (doseq [_ (range 300)]
      (let [c    [(r g -1.5 1.5) (r g -1.5 1.5) (r g -1.5 1.5)]
            want (apply min (map #(- 2.0 (Math/abs (double %))) c))
            got  (:depth (gjk/penetration (gjk/box [0 0 0] [2 2 2])
                                          (gjk/box c (mapv + c [2 2 2]))))]
        (is (about want got 1e-9) (str "c=" c))))))

(deftest pushing-by-the-penetration-vector-separates
  (let [g (rng 7)]
    (doseq [_ (range 200)]
      (let [c [(r g -1.5 1.5) (r g -1.5 1.5) (r g -1.5 1.5)]
            a (gjk/box [0 0 0] [2 2 2])
            b (gjk/box c (mapv + c [2 2 2]))]
        (when (gjk/intersects? a b)
          (let [{:keys [depth normal]} (gjk/penetration a b)
                moved (gjk/translate b (mapv #(* % (+ depth 0.01)) normal))]
            (is (not (gjk/intersects? a moved))
                (str "still overlapping after push, c=" c))))))))
