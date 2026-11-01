(ns allgo.vision.stereo-test
  (:require [allgo.random :as random]
            [allgo.vision.image :as img]
            [allgo.vision.stereo :as stereo]
            [clojure.test :refer [deftest is testing]]))

;; Binocular stereo, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 25, on a synthetic pair whose
;; disparities are known: a textured backdrop, and nearer, a textured box.

(def ^:private layers
  [{:texture (img/pattern (random/rng 5) 8 5.0) :disparity 2}
   {:texture (img/pattern (random/rng 6) 8 4.0) :disparity 7 :box [14 10 34 26]}])

(deftest disparity-test
  (let [[left right] (stereo/stereo-pair layers 48 36)
        truth (stereo/stereo-truth layers 48 36)
        d (stereo/disparity-map left right 3 10)
        pixels (for [y (range 36) x (range 14 48)] [x y])]
    (testing "the disparity is right at nine pixels in ten, the misses at the box's borders where one eye sees what the other cannot"
      (is (> (count (filter (fn [[x y]] (= (img/at d x y) (img/at truth x y))) pixels)) (* 0.9 (count pixels)))))
    (testing "and the box comes out nearer than the backdrop"
      (let [depths (stereo/depth-map d 500.0 0.1)]
        (is (< (img/at depths 24 18) (img/at depths 40 4)))))))
