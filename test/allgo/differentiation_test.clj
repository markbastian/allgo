(ns allgo.differentiation-test
  "Numerical differentiation against derivatives known in closed form."
  (:require [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest numerical-jacobian
  (testing "against the Jacobian of polar coordinates, known in closed form"
    (let [f (fn [[x y]] [(math/hypot x y) (math/atan2 y x)])
          [x y] [3.0 -4.0] r2 25.0 r 5.0
          J (diff/jacobian f [x y] {:angles #{1}})]
      (is (every? #(< (abs %) 1e-10)
                  (flatten (lin/mat-sub J [[(/ x r) (/ y r)] [(- (/ y r2)) (/ x r2)]]))))))
  (testing "and across the 2 pi seam, which the wrapping keeps smooth"
    (let [J (diff/jacobian (fn [[x y]] [(math/atan2 y x)]) [-1.0 1e-12] {:angles #{0}})]
      (is (< (abs (- (get-in J [0 1]) -1.0)) 1e-8)))))

(deftest plain-central-differences
  (testing "without Richardson's refinement, the h^2 error shows"
    (let [f (fn [[x]] [(math/exp x)])
          plain (get-in (diff/jacobian f [1.0] {:steps [1e-2] :richardson? false}) [0 0])
          refined (get-in (diff/jacobian f [1.0] {:steps [1e-2]}) [0 0])]
      ;; the plain difference is off by e h^2 / 6
      (is (< (abs (- (- plain math/E) (/ (* math/E 1e-4) 6.0))) 1e-9))
      (is (< (abs (- refined math/E)) 1e-9)))))
