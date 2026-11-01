(ns allgo.vision.image-test
  (:require [allgo.random :as random]
            [allgo.vision.image :as img]
            [clojure.test :refer [deftest is testing]]))

;; Images and filtering, after Russell and Norvig, *Artificial Intelligence:
;; A Modern Approach*, 4th edition, chapter 25.

(defn- close? [a b tol] (< (abs (- a b)) tol))

(deftest filtering-test
  (testing "a Gaussian kernel sums to one and is symmetric"
    (let [k (img/gaussian-kernel 1.5)]
      (is (close? 1.0 (reduce + k) 1e-12))
      (is (= k (vec (reverse k))))))
  (testing "smoothing leaves a constant image alone and keeps a noisy one's mean"
    (let [flat (img/image 10 8 (constantly 0.4))
          noisy (img/add-noise flat 0.1 (random/rng 1))
          mean (fn [{:keys [px]}] (/ (reduce + px) (count px)))]
      (is (every? #(close? 0.4 % 1e-12) (:px (img/smooth flat 1.0))))
      (is (close? (mean noisy) (mean (img/smooth noisy 1.0)) 0.01))
      (is (< (apply max (:px (img/smooth noisy 2.0))) (apply max (:px noisy))) "and evens out its extremes")))
  (testing "the gradient of a ramp is its slope, away from the edges"
    (let [ramp (img/image 12 12 (fn [x y] (+ (* 0.3 x) (* -0.2 y))))
          [gx gy] (img/gradients ramp)]
      (is (close? 0.3 (img/at gx 5 5) 1e-12))
      (is (close? -0.2 (img/at gy 5 5) 1e-12))))
  (testing "a color image filters a channel at a time"
    (let [c (img/image 6 6 (fn [x _] [(* 0.1 x) 0.5 1.0]))]
      (is (every? #(close? 0.5 (% 1) 1e-12) (:px (img/smooth c 1.0)))))))

(deftest sampling-test
  (testing "bilinear sampling is exact on pixels and linear between them"
    (let [ramp (img/image 8 8 (fn [x y] (+ x (* 10 y))))]
      (is (close? 23.0 (img/sample ramp 3.0 2.0) 1e-12))
      (is (close? 28.25 (img/sample ramp 3.25 2.5) 1e-12)))))

(deftest render-test
  (testing "shapes drawn in order, later ones in front, sampled at pixel centers"
    (let [im (img/render 10 10 0.0 [{:kind :rect :box [2 2 8 8] :value 0.5}
                                    {:kind :disk :center [5 5] :radius 1.5 :value 1.0}])]
      (is (= 0.0 (img/at im 0 0)))
      (is (= 0.5 (img/at im 2 2)))
      (is (= 1.0 (img/at im 4 4))))))
