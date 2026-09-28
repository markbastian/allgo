(ns allgo.interpolation-test
  "Derivatives of the quadratic through three points."
  (:require [allgo.numerics.interpolation :as interp]
            [clojure.test :refer [deftest is testing]]))

(deftest derivatives-at-middle
  (testing "exact for a quadratic, whatever the spacing"
    (let [f (fn [t] [(+ 1.0 (* 2.0 t) (* 3.0 t t)) (- 4.0 (* 0.5 t t))])
          [d1 d2] (interp/derivatives-at-middle [(f -0.7) (f 0.0) (f 1.3)] -0.7 1.3)]
      (is (every? #(< (abs %) 1e-12) (map - d1 [2.0 0.0])))
      (is (every? #(< (abs %) 1e-12) (map - d2 [6.0 -1.0]))))))
