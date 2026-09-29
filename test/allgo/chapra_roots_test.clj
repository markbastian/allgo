(ns allgo.chapra-roots-test
  "Roots of equations and of polynomials (Chapra and Canale, *Numerical
  Methods for Engineers*, 7th ed., chapters 5-7): the book's worked
  examples, and exact roots."
  (:require [allgo.numerics.polynomial :as poly]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (< (abs (- a b)) tol))

(def ^:private f (fn [x] (- (math/exp (- x)) x)))
(def ^:private df (fn [x] (- (- (math/exp (- x))) 1.0)))
(def ^:private root 0.56714329040978387)

(defn- parachutist
  "Example 5.1: the drag coefficient c at which a 68.1 kg parachutist
  reaches 40 m/s after 10 s."
  [c]
  (- (* (/ (* 9.81 68.1) c) (- 1.0 (math/exp (- (* (/ c 68.1) 10.0))))) 40.0))

(deftest bracketing
  (testing "Example 5.5: false position from 12 and 16, the parachutist's drag coefficient"
    (is (close? (roots/false-position parachutist 12.0 16.0) 14.8011 1e-4))
    (is (close? (roots/false-position parachutist 12.0 16.0 {:modified? false}) 14.8011 1e-4)))
  (testing "the first chord at 14.9309, as the example computes it (its prose says 14.9113, a slip)"
    (is (close? (- 16.0 (/ (* (parachutist 16.0) (- 12.0 16.0)) (- (parachutist 12.0) (parachutist 16.0)))) 14.9309 1e-4)))
  (testing "incremental search finds every sign change of sin 10x + cos 3x on [3, 6]"
    (let [g (fn [x] (+ (math/sin (* 10.0 x)) (math/cos (* 3.0 x))))
          brackets (roots/incremental-search g 3.0 6.0 300)]
      (is (= 9 (count brackets)))
      (is (every? (fn [[a b]] (neg? (* (g a) (g b)))) brackets))))
  (testing "a slow case for plain false position, x^10 - 1 on [0, 1.3], which the modification fixes"
    (let [g #(- (math/pow % 10.0) 1.0)]
      (is (close? (roots/false-position g 0.0 1.3) 1.0 1e-12)))))

(deftest open-methods
  (testing "Example 6.3: Newton-Raphson for e^-x - x from 0"
    (is (close? (roots/newton f df 0.0) root 1e-15)))
  (testing "Example 6.6: the secant method from 0 and 1"
    (is (close? (roots/secant f 0.0 1.0) root 1e-15)))
  (testing "Example 6.8: the modified secant, delta 0.01, from 1"
    (is (close? (roots/modified-secant f 1.0 0.01) root 1e-15)))
  (testing "Brent's method from the bracket [0, 1], and on the parachutist"
    (is (close? (roots/brent f 0.0 1.0) root 1e-15))
    (is (close? (roots/brent parachutist 12.0 16.0) (roots/newton parachutist
                                                                  #(/ (- (parachutist (+ % 1e-6)) (parachutist (- % 1e-6))) 2e-6)
                                                                  14.0)
                1e-9)))
  (testing "and nil without a sign change"
    (is (nil? (roots/brent f 1.0 2.0)))))

(deftest multiple-roots
  (testing "Section 6.5: the double root of (x - 3)(x - 1)^2, from 0 -- quadratically"
    (let [g (fn [x] (* (- x 3.0) (- x 1.0) (- x 1.0)))
          dg (fn [x] (+ (* (- x 1.0) (- x 1.0)) (* 2.0 (- x 3.0) (- x 1.0))))
          d2g (fn [x] (+ (* 4.0 (- x 1.0)) (* 2.0 (- x 3.0))))]
      (is (close? (roots/newton-multiple g dg d2g 0.0) 1.0 1e-12)))))

(deftest systems
  (let [F (fn [[x y]] [(+ (* x x) (* x y) -10.0) (+ y (* 3.0 x y y) -57.0)])]
    (testing "Example 6.12: Newton-Raphson for x^2 + xy = 10, y + 3xy^2 = 57 from (1.5, 3.5)"
      (is (every? true? (map #(close? %1 %2 1e-12) (roots/newton-system F [1.5 3.5]) [2.0 3.0])))
      (testing "with the Jacobian given"
        (is (every? true? (map #(close? %1 %2 1e-12)
                               (roots/newton-system F (fn [[x y]] [[(+ (* 2.0 x) y) x] [(* 3.0 y y) (+ 1.0 (* 6.0 x y))]]) [1.5 3.5])
                               [2.0 3.0])))))
    (testing "Example 6.11: fixed-point iteration, the convergent rearrangement"
      (let [g (fn [[x y]] [(math/sqrt (- 10.0 (* x y))) (math/sqrt (/ (- 57.0 y) (* 3.0 x)))])]
        (is (every? true? (map #(close? %1 %2 1e-10) (roots/fixed-point-system g [1.5 3.5]) [2.0 3.0])))))
    (testing "and the divergent one gives nil"
      (let [g (fn [[x y]] [(/ (- 10.0 (* x x)) y) (- 57.0 (* 3.0 x y y))])]
        (is (nil? (roots/fixed-point-system g [1.5 3.5])))))))

(deftest polynomials
  (testing "Example 7.1: (x^2 + 2x - 24)/(x - 4) = x + 6"
    (is (= {:quotient [1.0 6.0] :remainder 0.0} (poly/deflate [1.0 2.0 -24.0] 4.0))))
  (testing "long division: x^3 + 4 = (x - 1)(x^2 + x - 1) + 2x + 3"
    (is (= {:quotient [1.0 -1.0] :remainder [2.0 3.0]}
           (update (poly/long-divide [1.0 0.0 0.0 4.0] [1.0 1.0 -1.0]) :remainder #(mapv double %)))))
  (testing "value and derivatives by synthetic division"
    ;; x^3 - 2x^2 + 3x - 4 at 2: 2, 7, 8, 6
    (is (= [2.0 7.0 8.0 6.0] (poly/derivatives [1.0 -2.0 3.0 -4.0] 2.0 3))))
  (testing "Example 7.2: Muller's method on x^3 - 13x - 12 from 4.5, 5.5, 5 reaches 4"
    (let [[re im] (poly/muller [1.0 0.0 -13.0 -12.0] 4.5 5.5 5.0)]
      (is (close? re 4.0 1e-12)) (is (zero? im))))
  (testing "and finds a complex root of x^2 + 1 from real guesses"
    (let [[re im] (poly/muller [1.0 0.0 1.0] 0.5 1.0 1.5)]
      (is (close? re 0.0 1e-12)) (is (close? (abs im) 1.0 1e-12))))
  (testing "Example 7.3: Bairstow's method on x^5 - 3.5x^4 + 2.75x^3 + 2.125x^2 - 3.875x + 1.25"
    (let [rs (sort (poly/bairstow [1.0 -3.5 2.75 2.125 -3.875 1.25]))]
      (is (= 5 (count rs)))
      (is (every? true? (map (fn [[a b] [c d]] (and (close? a c 1e-10) (close? b d 1e-10)))
                             rs [[-1.0 0.0] [0.5 0.0] [1.0 -0.5] [1.0 0.5] [2.0 0.0]])))))
  (testing "Bairstow on a tenth-degree polynomial with known roots"
    (let [want [-3.0 -2.0 -1.0 0.5 1.5 2.5 4.0]
          ;; (x^2 + 1)(x^2 + 2x + 5) times the real factors
          mul (fn [a b] (reduce (fn [acc [i x]] (reduce (fn [acc [j y]] (update acc (+ i j) + (* x y))) acc (map-indexed vector b)))
                                (vec (repeat (+ (count a) (count b) -1) 0.0)) (map-indexed vector a)))
          p (reduce mul [1.0 0.0 1.0] (concat [[1.0 2.0 5.0]] (map (fn [r] [1.0 (- r)]) want)))
          got (poly/bairstow p)]
      (is (= 11 (count got)))
      (is (every? true? (map #(close? %1 %2 1e-8) (poly/real-roots got) want))))))
