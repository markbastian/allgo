(ns allgo.optimize-test
  (:require [allgo.numerics.optimize :as opt]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(< (abs (- %1 %2)) tol) a b)))

(deftest rosenbrock
  (testing "down Rosenbrock's curved valley from its classic start to (1, 1)"
    (let [f (fn [[x y]] (+ (* 100.0 (math/pow (- y (* x x)) 2)) (math/pow (- 1.0 x) 2)))
          {:keys [x f]} (opt/nelder-mead f [-1.2 1.0] {:max-iter 2000})]
      (is (close? x [1.0 1.0] 1e-6))
      (is (< f 1e-12)))))

(deftest quadratic
  (testing "a skewed quadratic bowl in five variables"
    (let [c [1.0 -2.0 3.0 0.5 -0.25]
          f (fn [x] (reduce + (map-indexed (fn [i [a b]] (* (inc i) (math/pow (- a b) 2)))
                                           (map vector x c))))
          {:keys [x]} (opt/nelder-mead f [0.0 0.0 0.0 0.0 0.0] {:step 1.0 :max-iter 5000})]
      (is (close? x c 1e-6)))))

(deftest undefined-regions
  (testing "nil where the function is undefined is uphill"
    (let [f (fn [[x]] (when (> x 0.5) (math/pow (- x 2.0) 2)))
          {:keys [x]} (opt/nelder-mead f [0.6] {:step 5.0})]
      (is (close? x [2.0] 1e-6)))))
