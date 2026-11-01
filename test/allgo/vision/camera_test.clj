(ns allgo.vision.camera-test
  (:require [allgo.vision.camera :as camera]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Image formation and two-view geometry, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25,
;; checked by the geometry worked out here.

(defn- close? [a b tol] (every? #(< (abs %) tol) (map - a b)))

(def ^:private cam (camera/look-at [0.0 0.0 -5.0] [0.0 0.0 0.0] 500.0 [320.0 240.0]))

(deftest projection-test
  (testing "the target lands on the center; a unit up lands f / Z above it"
    (is (close? [320.0 240.0] (camera/project cam [0.0 0.0 0.0]) 1e-9))
    (is (close? [320.0 140.0] (camera/project cam [0.0 1.0 0.0]) 1e-9)))
  (testing "twice as far, half the size: perspective"
    (let [near (camera/project cam [0.0 1.0 0.0]) far (camera/project cam [0.0 1.0 5.0])]
      (is (close? [50.0] [(- 240.0 (far 1))] 1e-9))
      (is (close? [100.0] [(- 240.0 (near 1))] 1e-9))))
  (testing "parallel lines converge on a vanishing point -- for lines along the view, the center"
    (let [at (fn [x z] (camera/project cam [x -1.0 z]))
          gap (fn [z] (let [[a b] (at -2.0 z) [c d] (at 2.0 z)] (math/hypot (- a c) (- b d))))]
      (is (> (gap 10.0) (gap 100.0) (gap 1000.0)))
      (is (close? [320.0 240.0] (at 2.0 1e9) 1e-3))))
  (testing "nothing behind the camera is seen"
    (is (nil? (camera/project cam [0.0 0.0 -10.0]))))
  (testing "scaled orthographic projection approaches perspective as an object's depth shrinks relative to its distance"
    (let [err (fn [d] (let [p [1.0 1.0 0.5]
                            c (camera/look-at [0.0 0.0 (- d)] [0.0 0.0 0.0] 500.0 [0.0 0.0])
                            [a b] (camera/project c p) [o q] (camera/project-orthographic c d p)]
                        (math/hypot (- a o) (- b q))))]
      (is (> (err 5.0) (err 50.0) (err 500.0))))))

(deftest triangulation-test
  (testing "a point seen by two cameras is recovered from its two images"
    (let [cam2 (camera/look-at [1.0 0.5 -5.0] [0.0 0.0 0.0] 500.0 [320.0 240.0])]
      (doseq [p [[0.3 -0.4 0.7] [-1.0 1.0 2.0] [0.0 0.0 0.0]]]
        (is (close? p (camera/triangulate cam (camera/project cam p) cam2 (camera/project cam2 p)) 1e-9)))))
  (testing "a ray back through a point's image passes through the point"
    (let [p [0.5 -0.2 1.5]
          [o d] (camera/ray cam (camera/project cam p))
          t (/ (reduce + (map * (map - p o) d)) 1.0)]
      (is (close? p (map + o (map #(* t %) d)) 1e-9))))
  (testing "stereo depth: f b / d"
    (is (== 10.0 (camera/stereo-depth 500.0 0.2 10.0)))))

(deftest lambert-test
  (testing "brightest facing the light, by the cosine off it, and ambient only when turned away"
    (is (close? [0.8] [(camera/lambert 0.8 [0.0 0.0 1.0] [0.0 0.0 1.0] 1.0 0.0)] 1e-12))
    (is (close? [0.4] [(camera/lambert 0.8 [0.0 0.0 1.0] [(math/sin (/ math/PI 3)) 0.0 0.5] 1.0 0.0)] 1e-12))
    (is (close? [0.08] [(camera/lambert 0.8 [0.0 0.0 1.0] [0.0 0.0 -1.0] 1.0 0.1)] 1e-12))))
