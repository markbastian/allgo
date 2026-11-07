(ns allgo.vision.flow-test
  (:require [allgo.random :as random]
            [allgo.vision.flow :as flow]
            [allgo.vision.image :as img]
            [clojure.test :refer [deftest is testing]]))

;; Optical flow, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 25, on a textured image moved by
;; a known amount.

(def ^:private texture (img/pattern (random/rng 3) 6 8.0))
(defn- frame [dx dy] (img/image 40 40 (fn [x y] (texture (- x dx) (- y dy)))))

(deftest block-matching-test
  (testing "a whole-pixel move is found exactly, everywhere on a grid"
    (let [field (flow/flow-field (frame 0 0) (frame 2 1) 8 6 3 4 :block)]
      (is (= 16 (count field)))
      (is (every? #(= [2 1] %) (vals field))))))

(deftest lucas-kanade-test
  (testing "a whole-pixel move to a hair, and a subpixel one to a few hundredths"
    (let [[u v] (flow/lucas-kanade (frame 0 0) (frame 2 1) 20 20 4)]
      (is (< (abs (- u 2.0)) 1e-6))
      (is (< (abs (- v 1.0)) 1e-6)))
    (let [[u v] (flow/lucas-kanade (frame 0 0) (frame 0.6 -0.3) 20 20 4)]
      (is (< (abs (- u 0.6)) 0.03))
      (is (< (abs (- v -0.3)) 0.03))))
  (testing "a featureless window cannot say how it moved: the aperture problem"
    (is (nil? (flow/lucas-kanade-seq (img/image 20 20 (constantly 0.5)) (img/image 20 20 (constantly 0.5)) 10 10 3)))))

(deftest moving-camera-test
  (let [t [0.2 -0.1 0.5]
        points (for [x (range -0.5 0.51 0.25) y (range -0.4 0.41 0.2)] [x y (+ 3.0 (* 4 (abs (+ x y))))])
        field (into {} (map (fn [[x y z]] [[x y] (flow/translational-flow t x y z)]) points))]
    (testing "the flow of a translating camera radiates from the focus of expansion, (Tx / Tz, Ty / Tz)"
      (let [[x0 y0] (flow/focus-of-expansion field)]
        (is (< (abs (- x0 0.4)) 1e-12))
        (is (< (abs (- y0 -0.2)) 1e-12))))
    (testing "and every point's distance from it over its speed is Z / Tz, the time to contact, whatever the depth"
      (is (every? (fn [[x y z]] (< (abs (- (/ z 0.5) (flow/time-to-contact [0.4 -0.2] [x y] (field [x y])))) 1e-9)) points)))
    (testing "moving sideways, the inverse of the ratio of two flows' lengths is the ratio of their depths: motion parallax"
      (is (< (abs (- (/ 3.0 7.0) (flow/depth-ratio (flow/translational-flow [1.0 0.0 0.0] 0.1 0.2 3.0)
                                                   (flow/translational-flow [1.0 0.0 0.0] -0.3 0.0 7.0)))) 1e-12))))
  (testing "measured: Lucas-Kanade's flow as the camera nears a textured wall finds the focus at the center, and the time to contact"
    (let [pat (img/pattern (random/rng 3) 12 9.0)
          z 10.0 tz 0.5 s (/ z (- z tz))
          a (img/image 64 64 (fn [x y] (pat (+ x 0.5) (+ y 0.5))))
          b (img/image 64 64 (fn [x y] (pat (+ 32 (/ (- (+ x 0.5) 32) s)) (+ 32 (/ (- (+ y 0.5) 32) s)))))
          field (flow/flow-field a b 6 12 4 3 :lucas-kanade)
          [x0 y0 :as foe] (flow/focus-of-expansion field)
          ttcs (for [[p v] field :when (> (abs (- (p 0) x0)) 8)] (flow/time-to-contact foe p v))]
      ;; pixel x is centered at x + 0.5: the true focus is at 31.5
      (is (< (abs (- x0 31.5)) 0.5))
      (is (< (abs (- y0 31.5)) 0.5))
      ;; the second frame is 0.5 nearer: (Z - Tz) / Tz frames to go
      (is (< (abs (- 19.0 (/ (reduce + ttcs) (count ttcs)))) 0.5)))))
