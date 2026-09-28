(ns allgo.us76-test
  "The 1976 standard atmosphere's defining equations against its own Table
  I, and the exponential model built on it."
  (:require [allgo.astro.us76 :as us76]
            [clojure.test :refer [deftest is testing]]))

(defn- within? [a b tol] (< (abs (- a b)) (* tol (abs b))))

(deftest defining-equations
  (testing "below 86 km, the equations give Table I's densities to its printed digits"
    ;; Table I, geometric altitude, metric units
    (doseq [[z rho] [[0.0 1.2250] [30.0 1.8410e-2] [42.0 2.9948e-3] [46.0 1.7142e-3]
                     [50.0 1.0269e-3] [86.0 6.958e-6]]]
      (is (within? (:rho (us76/lower z)) rho 1e-4) (str z " km")))))

(deftest exponential-model
  (testing "exact at each base, from the equations below 86 km and Table I above"
    (doseq [[h rho] [[0.0 1.225] [100.0 5.604e-7] [400.0 2.803e-12] [1000.0 3.561e-15]]]
      (is (within? (us76/exponential-density h) rho 1e-3) (str h))))
  (testing "continuous across each base"
    (doseq [[h] (rest us76/exponential-table)]
      (is (within? (us76/exponential-density (- h 1e-9)) (us76/exponential-density h) 1e-8) (str h))))
  (testing "and between bases within a few percent of the table it stands in for"
    ;; Table I's values at 125 and 170 km, between bases
    (is (within? (us76/exponential-density 125.0) 1.291e-8 0.06))
    (is (within? (us76/exponential-density 170.0) 7.815e-10 0.07))
    (is (within? (us76/exponential-density 374.0) 4.478e-12 0.1)))
  (testing "falling all the way up, and on past 1000 km"
    (is (apply > (map us76/exponential-density (range 0.0 1500.0 50.0))))))
