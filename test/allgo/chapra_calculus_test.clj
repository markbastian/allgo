(ns allgo.chapra-calculus-test
  "Numerical integration and differentiation (Chapra and Canale,
  *Numerical Methods for Engineers*, 7th ed., chapters 21-23): the book's
  running example, f(x) = 0.2 + 25x - 200x^2 + 675x^3 - 900x^4 + 400x^5
  on [0, 0.8], by every rule, and its other worked examples."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.quadrature :as q]
            [allgo.numerics.regression :as reg]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(defn- f [x] (+ 0.2 (* 25.0 x) (* -200.0 x x) (* 675.0 x x x) (* -900.0 x x x x) (* 400.0 x x x x x)))
(def ^:private exact
  (let [x 0.8] (+ (* 0.2 x) (* 12.5 x x) (* (/ -200.0 3.0) x x x) (* 168.75 x x x x) (* -180.0 (math/pow x 5.0)) (* (/ 400.0 6.0) (math/pow x 6.0)))))

(deftest newton-cotes
  (testing "the exact integral, 1.640533"
    (is (close? exact 1.640533 1e-6)))
  (testing "Examples 21.1, 21.4, 21.6: single applications -- 0.1728, 1.367467, 1.519170"
    (is (close? (q/newton-cotes f 0.0 0.8 2) 0.1728 1e-12))
    (is (close? (q/newton-cotes f 0.0 0.8 3) 1.367467 1e-6))
    (is (close? (q/newton-cotes f 0.0 0.8 4) 1.519170 1e-6))
    (is (close? (q/simpson-3-8 f 0.0 0.8 3) 1.519170 1e-6)))
  (testing "Boole's and the six-point rule are exact for this quintic"
    (is (close? (q/newton-cotes f 0.0 0.8 5) exact 1e-12))
    (is (close? (q/newton-cotes f 0.0 0.8 6) exact 1e-12)))
  (testing "the open formulas are exact to their degree"
    (let [cubic #(+ 1.0 (* 2.0 %) (* -3.0 % %) (* 4.0 % % %))
          exact3 (- (+ 1.0 1.0 -1.0 1.0) 0.0)]
      (doseq [n [4 5 6]]
        (is (close? (q/newton-cotes cubic 0.0 1.0 n {:open? true}) exact3 1e-12)))))
  (testing "Examples 21.7-21.8: the unequally spaced data of Table 21.3"
    (let [xs [0.0 0.12 0.22 0.32 0.36 0.40 0.44 0.54 0.64 0.70 0.80]
          ys (map f xs)]
      (is (close? (q/trapezoid-data xs ys) 1.594801 1e-6))
      (is (close? (q/simpson-data xs ys) 1.603641 1e-6))))
  (testing "Example 21.9: the plate's average temperature, 58.66667"
    (let [T (fn [x y] (+ (* 2.0 x y) (* 2.0 x) (- (* x x)) (* -2.0 y y) 72.0))]
      (is (close? (/ (q/double-integral T 0.0 8.0 0.0 6.0 2 2) 48.0) (/ 176.0 3.0) 1e-12)))))

(deftest integration-of-equations
  (testing "Example 22.2: Romberg integration reaches the exact value"
    (is (close? (:value (q/romberg f 0.0 0.8)) exact 1e-12)))
  (testing "adaptive quadrature on the square root's awkward start"
    (is (close? (q/adaptive math/sqrt 0.0 1.0) (/ 2.0 3.0) 1e-9)))
  (testing "Examples 22.3-22.4: two-point Gauss-Legendre gives 1.822578, three points the exact value"
    (is (close? (q/gauss-legendre f 0.0 0.8 2) 1.822578 1e-6))
    (is (close? (q/gauss-legendre f 0.0 0.8 3) exact 1e-12)))
  (testing "Gauss-Legendre of high order on e^x"
    (is (close? (q/gauss-legendre math/exp -1.0 2.0 20) (- (math/exp 2.0) (math/exp -1.0)) 1e-13)))
  (testing "Example 22.6: the cumulative normal N(1) = 0.8413, an improper integral"
    (let [phi #(/ (math/exp (* -0.5 % %)) (math/sqrt (* 2.0 math/PI)))]
      (is (close? (q/improper phi ##-Inf 1.0) 0.8413447460685429 1e-9))))
  (testing "and other improper integrals"
    (is (close? (q/improper #(/ 1.0 (* % %)) 1.0 ##Inf) 1.0 1e-9))
    (is (close? (q/improper #(math/exp (- (* % %))) ##-Inf ##Inf) (math/sqrt math/PI) 1e-9))))

(defn- g [x] (+ (* -0.1 x x x x) (* -0.15 x x x) (* -0.5 x x) (* -0.25 x) 1.2))

(deftest differentiation
  (testing "Example 23.1: at x = 0.5 with h = 0.25, the true -0.9125"
    (is (close? (d/derivative g 0.5 0.25 {:scheme :forward}) -0.859375 1e-7))
    (is (close? (d/derivative g 0.5 0.25 {:scheme :backward}) -0.878125 1e-7))
    (is (close? (d/derivative g 0.5 0.25 {:accuracy 4}) -0.9125 1e-12)))
  (testing "the weights are the familiar formulas"
    (is (close? (nth (d/fd-weights 0.0 [-1 0 1] 2) 1) [-0.5 0.0 0.5] 1e-15))
    (is (close? (nth (d/fd-weights 0.0 [-1 0 1] 2) 2) [1.0 -2.0 1.0] 1e-15))
    (is (close? (nth (d/fd-weights 0.0 [-2 -1 0 1 2] 1) 1) [(/ 1.0 12.0) (/ -2.0 3.0) 0.0 (/ 2.0 3.0) (/ -1.0 12.0)] 1e-15)))
  (testing "second, third and fourth derivatives of a quartic"
    (let [p #(+ (* 3.0 % % % %) (* -2.0 % % %) %)]
      (is (close? (d/derivative p 0.7 1e-3 {:order 2}) (- (* 36.0 0.49) (* 12.0 0.7)) 1e-5))
      (is (close? (d/derivative p 0.7 1e-2 {:order 3 :accuracy 4}) (- (* 72.0 0.7) 12.0) 1e-7))
      (is (close? (d/derivative p 0.7 0.1 {:order 4}) 72.0 1e-9))))
  (testing "Example 23.2: Richardson extrapolation of the centered differences at h = 0.5 and 0.25"
    (is (close? (d/richardson #(d/derivative g 0.5 %) 0.5) -0.9125 1e-12)))
  (testing "Example 23.3: the soil's temperature gradient at the surface, -1.333333 C/cm"
    (is (close? (d/data-derivative [0.0 1.25 3.75] [13.5 12.0 10.0] 0.0) -1.333333 1e-6)))
  (testing "a smoothed derivative of noisy data"
    (let [xs (range 0.0 5.01 0.25)
          ys (map-indexed (fn [i x] (+ (* 2.0 x x) (* 0.05 (if (even? i) 1.0 -1.0)))) xs)]
      (is (close? (reg/smoothed-derivative xs ys 2 2.0) 8.0 0.01))))
  (testing "partial and mixed partial derivatives of x^2 y^3"
    (let [h (fn [[x y]] (* x x y y y))]
      (is (close? (d/partial-derivative h [1.5 2.0] 0) (* 2.0 1.5 8.0) 1e-7))
      (is (close? (d/partial-derivative h [1.5 2.0] 1) (* 2.25 3.0 4.0) 1e-7))
      (is (close? (d/mixed-partial h [1.5 2.0] 0 1) (* 2.0 1.5 3.0 4.0) 1e-5)))))
