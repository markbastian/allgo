(ns allgo.chapra-optimization-test
  "Optimization (Chapra and Canale, *Numerical Methods for Engineers*,
  7th ed., chapters 13-15): the book's examples -- maxima, found as the
  minima of the negated functions -- Rosenbrock's valley for every
  multidimensional method, and linear and constrained programs with known
  answers."
  (:require [allgo.numerics.optimize :as opt]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

;; Chapter 13's function, 2 sin x - x^2/10, its maximum at 1.4276
(defn- g [x] (- (* 2.0 (math/sin x)) (/ (* x x) 10.0)))
(def ^:private x-max 1.4275517787645942)

(deftest one-variable
  (testing "Example 13.2: parabolic interpolation from 0, 1 and 4"
    (is (close? (opt/parabolic #(- (g %)) 0.0 1.0 4.0) x-max 1e-8)))
  (testing "Example 13.3: Newton's method from 2.5"
    (is (close? (opt/newton-1d #(- (* 2.0 (math/cos %)) (/ % 5.0)) #(- (* -2.0 (math/sin %)) 0.2) 2.5) x-max 1e-14)))
  (testing "Brent's method on [0, 4], to the square root of the precision a minimum allows"
    (let [{:keys [x f]} (opt/brent-minimize #(- (g %)) 0.0 4.0)]
      (is (close? x x-max 2e-8))
      (is (close? (- f) 1.7757 1e-4)))))

(defn- rosenbrock [[x y]] (+ (* 100.0 (math/pow (- y (* x x)) 2.0)) (math/pow (- 1.0 x) 2.0)))
(defn- rosenbrock-grad [[x y]] [(- (* -400.0 x (- y (* x x))) (* 2.0 (- 1.0 x))) (* 200.0 (- y (* x x)))])

(deftest several-variables
  (testing "Example 14.1: random search for the maximum of y - x - 2x^2 - 2xy - y^2,
            1.25 at (-1, 1.5) (the book prints its value as 1.5)"
    (let [h (fn [[x y]] (- (- y x (* 2.0 x x) (* 2.0 x y) (* y y))))
          {:keys [x f]} (opt/random-search h [-2.0 1.0] [2.0 3.0] 20000 {:rng (random/rng 7)})]
      (is (close? x [-1.0 1.5] 0.05))
      (is (close? (- f) 1.25 1e-3))))
  (testing "Example 14.4: steepest ascent on 2xy + 2x - x^2 - 2y^2 from (-1, 1) to (2, 1)"
    (let [h (fn [[x y]] (- (- (+ (* 2.0 x y) (* 2.0 x)) (* x x) (* 2.0 y y))))]
      (is (close? (:x (opt/steepest-descent h [-1.0 1.0])) [2.0 1.0] 1e-6))
      (testing "and every other method on it"
        (doseq [m [opt/univariate opt/powell opt/conjugate-gradient opt/newton opt/marquardt opt/quasi-newton]]
          (is (close? (:x (m h [-1.0 1.0])) [2.0 1.0] 1e-6))))))
  (testing "Rosenbrock's valley from (-1.2, 1), the hard case"
    (doseq [[label run] [["Powell" #(opt/powell rosenbrock [-1.2 1.0])]
                         ["Fletcher-Reeves" #(opt/conjugate-gradient rosenbrock rosenbrock-grad [-1.2 1.0])]
                         ["Newton" #(opt/newton rosenbrock [-1.2 1.0])]
                         ["Marquardt" #(opt/marquardt rosenbrock [-1.2 1.0])]
                         ["BFGS" #(opt/quasi-newton rosenbrock rosenbrock-grad [-1.2 1.0])]
                         ["DFP" #(opt/quasi-newton rosenbrock rosenbrock-grad [-1.2 1.0] {:update :dfp})]]]
      (testing label
        (is (close? (:x (run)) [1.0 1.0] 1e-4)))))
  (testing "a quadratic in n variables: conjugate gradients in n line searches"
    (let [q (fn [x] (reduce + (map-indexed (fn [i xi] (* (inc i) (- xi i) (- xi i))) x)))
          {:keys [x iterations]} (opt/conjugate-gradient q (vec (repeat 5 0.0)))]
      (is (close? x [0.0 1.0 2.0 3.0 4.0] 1e-6))
      (is (<= iterations 6)))))

(deftest linear-programs
  (testing "Examples 15.1-15.2: the gas-processing plant, Z = 150 x1 + 175 x2 at its most"
    (let [{:keys [x value]} (opt/linear-program [-150.0 -175.0]
                                                [[7.0 11.0] [10.0 8.0] [1.0 0.0] [0.0 1.0]]
                                                [77.0 80.0 9.0 6.0])]
      (is (close? x [(/ 44.0 9.0) (/ 35.0 9.0)] 1e-9))
      (is (close? (- value) 1413.8888888888889 1e-9))))
  (testing ">= and = constraints: minimize x + y with x + 2y >= 4, 3x + y >= 6, x - y = 0.5"
    (let [{:keys [x value]} (opt/linear-program [1.0 1.0] [[1.0 2.0] [3.0 1.0] [1.0 -1.0]] [4.0 6.0 0.5]
                                                {:kinds [:>= :>= :=]})]
      ;; along x - y = 0.5 the first constraint binds first: x = 5/3, y = 7/6
      (is (close? x [(/ 5.0 3.0) (/ 7.0 6.0)] 1e-12))
      (is (close? value (/ 17.0 6.0) 1e-12))))
  (testing "infeasible and unbounded programs are said to be"
    (is (= :infeasible (:status (opt/linear-program [1.0] [[1.0] [1.0]] [1.0 2.0] {:kinds [:<= :>=]}))))
    (is (= :unbounded (:status (opt/linear-program [-1.0 0.0] [[1.0 -1.0]] [1.0]))))))

(deftest constrained
  (testing "penalty functions: the nearest point to (2, 1) with x + y <= 2 is (1.5, 0.5)"
    (let [{:keys [x]} (opt/penalty (fn [[x y]] (+ (math/pow (- x 2.0) 2.0) (math/pow (- y 1.0) 2.0)))
                                   [(fn [[x y]] (- (+ x y) 2.0))] [] [0.0 0.0])]
      (is (close? x [1.5 0.5] 1e-6))))
  (testing "and with x = y, (1.5, 1.5)"
    (let [{:keys [x]} (opt/penalty (fn [[x y]] (+ (math/pow (- x 2.0) 2.0) (math/pow (- y 1.0) 2.0)))
                                   [] [(fn [[x y]] (- x y))] [0.0 0.0])]
      (is (close? x [1.5 1.5] 1e-6)))))
