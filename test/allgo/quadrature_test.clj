(ns allgo.quadrature-test
  "The composite rules against integrals known in closed form."
  (:require [allgo.numerics.quadrature :as q]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest simpsons-rule
  (testing "exact for a cubic"
    (is (< (abs (- (q/simpson #(+ (* % % %) (* 2 %) 1.0) 0.0 2.0 4) 10.0)) 1e-12)))
  (testing "fourth-order: halving h cuts the error sixteenfold"
    (let [err #(abs (- (q/simpson math/exp 0.0 1.0 %) (- math/E 1.0)))]
      (is (< 15.0 (/ (err 8) (err 16)) 17.0)))))

(deftest trapezoid-rule
  (testing "second-order in general"
    (let [err #(abs (- (q/trapezoid math/exp 0.0 1.0 %) (- math/E 1.0)))]
      (is (< 3.9 (/ (err 8) (err 16)) 4.1))))
  (testing "but spectrally accurate over a period of a periodic function"
    (is (< (abs (- (q/trapezoid #(math/exp (math/cos %)) 0.0 (* 2 math/PI) 16)
                   (* 2 math/PI 1.2660658777520082)))
           1e-12))))
