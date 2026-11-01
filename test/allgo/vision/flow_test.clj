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
