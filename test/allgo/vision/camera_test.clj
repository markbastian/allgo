(ns allgo.vision.camera-test
  (:require [allgo.random :as random]
            [allgo.vision.camera :as camera]
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

(deftest stereopsis-test
  (testing "the chapter's figures for human eyes: a 6 cm baseline and a disparity of 5 seconds of arc"
    (let [arcsec (/ math/PI (* 180 3600))]
      (testing "resolve 0.4 mm at a meter"
        (is (< (abs (- 0.4 (* 10 (camera/depth-resolution 6.0 100.0 (* 5 arcsec))))) 0.005)))
      (testing "and 0.036 mm at 30 cm, enough to thread a needle"
        (is (< (abs (- 0.036 (* 10 (camera/depth-resolution 6.0 30.0 (* 5 arcsec))))) 0.0005)))
      (testing "the angular disparity is its inverse"
        (is (< (abs (- (* 5 arcsec) (camera/angular-disparity 6.0 100.0 (camera/depth-resolution 6.0 100.0 (* 5 arcsec))))) 1e-15))))))

(defn- unit [x] (let [n (math/sqrt (reduce + (map * x x)))] (mapv #(/ % n) x)))

(deftest shading-test
  (let [light (unit [0.5 -0.5 0.7071])
        matte {:albedo 0.6}
        shiny {:albedo 0.6 :specular 0.5 :shininess 30}]
    (testing "Lambert's law, plus the ambient: facing the light, rho (ambient + I0); turned away, rho ambient"
      (is (< (abs (- (* 0.6 1.1) (camera/shade matte light light [0.0 0.0 1.0] 1.0 0.1))) 1e-12))
      (is (< (abs (- (* 0.6 0.1) (camera/shade matte (mapv - light) light [0.0 0.0 1.0] 1.0 0.1))) 1e-12)))
    (testing "a point in shadow gets only the ambient"
      (is (< (abs (- (* 0.6 0.1) (camera/shade matte light light [0.0 0.0 1.0] 1.0 0.1 :shadowed true))) 1e-12)))
    (testing "the specularity is brightest where the light's mirror direction is the view: the normal halfway between them"
      (let [im (camera/sphere-image 81 81 36.0 shiny light 1.0 0.1 0.0)
            i (apply max-key (:px im) (range (count (:px im))))
            h (unit (mapv + light [0.0 0.0 1.0]))]
        (is (<= (abs (- (mod i 81) (+ 40 (* 36 (h 0))))) 1.0))
        (is (<= (abs (- (quot i 81) (+ 40 (* 36 (h 1))))) 1.0))))
    (testing "color a channel at a time: a red surface under white light, and a white one under red"
      (is (close? [0.6 0.0 0.0] (camera/shade {:albedo [0.6 0.0 0.0]} light light [0.0 0.0 1.0] 1.0 0.0) 1e-12))
      (is (close? [0.6 0.0 0.0] (camera/shade {:albedo 0.6} light light [0.0 0.0 1.0] [1.0 0.0 0.0] 0.0) 1e-12)))))

(deftest color-constancy-test
  (let [rng (random/rng 1)
        albedos (vec (repeatedly 40 (fn [] (vec (repeatedly 3 #(+ 0.1 (* 0.8 (rng))))))))
        light [1.0 0.8 0.5]
        under (fn [albedos light] {:w (count albedos) :h 1 :px (mapv #(mapv * % light) albedos)})
        white (/ (reduce + light) 3)]
    (testing "where the surfaces do average to gray, gray-world finds the light's color and undoes it exactly"
      (let [balanced (into albedos (map (fn [a] (mapv #(- 1.0 %) a)) albedos))
            {:keys [illuminant image]} (camera/gray-world (under balanced light))]
        (is (close? illuminant (mapv #(/ % white) light) 1e-12))
        (is (every? true? (map #(close? %1 %2 1e-12) (:px image) (:px (under balanced [white white white])))))))
    (testing "where they only roughly do, it comes close: each channel within 0.04 of the scene under white light"
      (let [{:keys [image]} (camera/gray-world (under albedos light))]
        (is (every? true? (map #(close? %1 %2 0.04) (:px image) (:px (under albedos [white white white])))))))))

(def ^:private house
  (mapv #(mapv double %) [[0 0 0] [2 0 0] [2 0 2] [0 0 2] [0 1.5 0] [2 1.5 0] [2 1.5 2] [0 1.5 2] [1 2.3 0] [1 2.3 2]]))

(deftest pose-test
  (doseq [[position target] [[[6.0 3.0 -5.0] [1.0 1.0 1.0]] [[-4.0 2.0 7.0] [1.0 0.8 1.0]] [[1.0 9.0 -2.0] [1.0 1.0 1.0]]]]
    (let [truth (camera/look-at position target [0.0 -1.0 0.0] 300.0 [160.0 120.0])
          pixels (mapv #(camera/project truth %) house)]
      (testing "from a known object's exact images, the camera where it was"
        (let [{:keys [camera rms]} (camera/pose house pixels 300.0 [160.0 120.0])]
          (is (close? position (:position camera) 1e-8))
          (is (< rms 1e-8))))
      (testing "and from images half a pixel out, within 0.15 of it, eight to ten away"
        (let [rng (java.util.Random. 1)
              noisy (mapv (fn [p] (mapv #(+ % (* 0.5 (.nextGaussian rng))) p)) pixels)]
          (is (close? position (:position (:camera (camera/pose house noisy 300.0 [160.0 120.0]))) 0.15)))))))
