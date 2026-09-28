(ns allgo.sphere-geometry-test
  "Rotations about the axes, great circles, lines against a sphere, and
  overlapping disks."
  (:require [allgo.geometry.disk :as disk]
            [allgo.geometry.rotation :as rot]
            [allgo.geometry.sphere :as sphere]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- near? [a b] (every? #(< (abs %) 1e-12) (map - (flatten a) (flatten b))))

(deftest rotations
  (doseq [[frame vector] [[rot/rx rot/rotate-x] [rot/ry rot/rotate-y] [rot/rz rot/rotate-z]]]
    (testing "turning the frame by t is turning the vector by -t, and the transpose"
      (is (near? (frame 0.4) (vector -0.4)))
      (is (near? (frame 0.4) (lin/transpose (vector 0.4))))))
  (testing "a vector along x turned a quarter about z points along y"
    (is (near? (lin/mat-vec (rot/rotate-z (/ math/PI 2)) [1.0 0.0 0.0]) [0.0 1.0 0.0])))
  (testing "a chain applies its last matrix first"
    (is (near? (rot/chain (rot/rz 0.3) (rot/rx 0.2)) (lin/mat-mul (rot/rz 0.3) (rot/rx 0.2))))))

(deftest great-circles
  (testing "a quarter of the equator, due east"
    (let [{:keys [angle azimuth]} (sphere/great-circle [0.0 0.0] [0.0 (/ math/PI 2)])]
      (is (< (abs (- angle (/ math/PI 2))) 1e-15))
      (is (< (abs (- azimuth (/ math/PI 2))) 1e-15))))
  (testing "setting out along the azimuth for the angle arrives"
    (let [p1 [0.7 -1.8] p2 [-0.3 2.4]
          {:keys [angle azimuth]} (sphere/great-circle p1 p2)]
      (is (near? (sphere/destination p1 angle azimuth) p2)))))

(deftest lines-and-spheres
  (testing "a segment passing just outside, and just inside"
    (is (sphere/segment-clears? [-5.0 1.001 0.0] [5.0 1.001 0.0] 1.0))
    (is (not (sphere/segment-clears? [-5.0 0.999 0.0] [5.0 0.999 0.0] 1.0))))
  (testing "the far intersection lands on the sphere"
    (let [o [0.2 -0.3 0.1] u (v3/normalize [1.0 2.0 -0.5])
          t (sphere/far-intersection o u 2.0)]
      (is (< (abs (- (v3/length (v3/add-scaled o u t)) 2.0)) 1e-12))
      (is (pos? t))))
  (testing "a line that never reaches the radius"
    (is (nil? (sphere/far-intersection [0.0 5.0 0.0] [1.0 0.0 0.0] 2.0))))
  (testing "the angle between vectors, clamped at the ends"
    (is (zero? (v3/angle [1.0 1.0 0.0] [2.0 2.0 0.0])))
    (is (< (abs (- (v3/angle [1.0 0.0 0.0] [0.0 0.0 3.0]) (/ math/PI 2))) 1e-15))))

(deftest overlapping-disks
  (is (zero? (disk/overlap-fraction 1.0 1.0 2.5)) "apart")
  (is (= 1.0 (disk/overlap-fraction 1.0 2.0 0.5)) "covered")
  (is (< (abs (- (disk/overlap-fraction 1.0 0.5 0.2) 0.25)) 1e-15) "a smaller disk inside hides its area's share")
  (testing "half-way: two equal disks a radius apart share 39 percent"
    ;; 2 acos(1/2) - (1/2) sqrt(3), over pi
    (is (< (abs (- (disk/overlap-fraction 1.0 1.0 1.0) (/ (- (* 2 (math/acos 0.5)) (* 0.5 (math/sqrt 3.0))) math/PI))) 1e-12))))
