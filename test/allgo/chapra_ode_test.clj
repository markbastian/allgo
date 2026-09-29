(ns allgo.chapra-ode-test
  "Ordinary differential equations, boundary-value problems and
  eigenvalues (Chapra and Canale, *Numerical Methods for Engineers*, 7th
  ed., chapters 25-27): the book's running example y' = 4e^(0.8t) - 0.5y,
  y(0) = 2, by every method whose results it prints, the stiff example,
  the heated rod, and the eigenproblems."
  (:require [allgo.numerics.bvp :as bvp]
            [allgo.numerics.eigen :as eigen]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.ode :as ode]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(defn- f [t [y]] [(- (* 4.0 (math/exp (* 0.8 t))) (* 0.5 y))])
(defn- exact [t] (+ (* (/ 4.0 1.3) (- (math/exp (* 0.8 t)) (math/exp (* -0.5 t)))) (* 2.0 (math/exp (* -0.5 t)))))
(defn- values [traj] (mapv (comp first second) (rest traj)))

(deftest one-step
  (testing "Example 25.5: Heun's method, plain and with the corrector iterated"
    (is (close? (values (ode/heun f 0.0 [2.0] 1.0 4)) [6.7010819 16.3197819 37.1992489 83.3377674] 1e-6))
    (is (close? (values (ode/heun f 0.0 [2.0] 1.0 4 {:iterations 15})) [6.3608655 15.3022367 34.7432761 77.7350962] 1e-6)))
  (testing "Example 25.13: one Cash-Karp step of 2 -- the fifth-order 14.83192, the fourth-order 14.83677"
    (let [[[y5] [e]] (rk/tableau-step rk/cash-karp (:c rk/cash-karp) f 0.0 [2.0] 2.0)]
      (is (close? y5 14.83192 1e-5))
      (is (close? (- y5 e) 14.83677 1e-5))))
  (testing "Example 25.12: step halving -- the halves' 14.86249 corrected to 14.84627"
    (let [traj (ode/step-halving f 0.0 [2.0] 2.0 {:h 2.0 :tol 1.0})]
      (is (= 2 (count traj)))
      ;; the book rounds its intermediate values to five decimals
      (is (close? (first (second (peek traj))) 14.84627 2e-5))))
  (testing "and adaptively to a tolerance, against the true solution"
    (let [traj (ode/step-halving f 0.0 [2.0] 4.0 {:tol 1e-9})]
      (is (close? (first (second (peek traj))) (exact 4.0) (* 1e-7 (exact 4.0))))))
  (testing "Butcher's fifth-order method at h = 0.5 against the true solution"
    (let [s (reduce (fn [[t y] _] [(+ t 0.5) (first (rk/tableau-step rk/rk5-butcher (:c rk/rk5-butcher) f t y 0.5))]) [0.0 [2.0]] (range 8))]
      (is (close? (first (second s)) (exact 4.0) 1e-3))))
  (testing "the second-order Taylor method converges at second order"
    (let [g (fn [_ [y]] [y]) gg (fn [_ [y]] [y])
          err (fn [n] (abs (- (first (second (peek (ode/taylor-2 g gg 0.0 [1.0] (/ 1.0 n) n)))) math/E)))]
      (is (< 3.8 (/ (err 100) (err 200)) 4.2)))))

(deftest stiffness
  (let [g (fn [t [y]] [(+ (* -1000.0 y) 3000.0 (* -2000.0 (math/exp (- t))))])
        truth (fn [t] (- 3.0 (* 0.998 (math/exp (* -1000.0 t))) (* 2.002 (math/exp (- t)))))]
    (testing "Example 26.1: implicit Euler at h = 0.05, far past the explicit limit of 0.002, stays on the
              solution -- to its first-order error on the slow part, some 2 percent"
      (let [traj (ode/backward-euler g 0.0 [0.0] 0.05 8)]
        (is (every? (fn [[t [y]]] (< (abs (- y (truth t))) (* 0.02 (truth t)))) (rest traj)))))
    (testing "where explicit Euler, at h = 0.0025, diverges"
      (let [traj (ode/heun g 0.0 [0.0] 0.0025 200)]
        (is (> (abs (first (second (peek traj)))) 1e6))))))

(deftest multistep
  (let [before [[-1.0 [-0.3929953]]]
        three-before [[-3.0 [-4.547302]] [-2.0 [-2.306160]] [-1.0 [-0.3929953]]]]
    (testing "Example 26.2: the non-self-starting Heun method, corrector converged"
      (is (close? (take 2 (values (ode/multistep :heun f 0.0 [2.0] 1.0 4 {:history before}))) [6.360865 15.30224] 1e-5)))
    (testing "Example 26.4: with the modifiers -- 6.210093, then 14.88827"
      (is (close? (take 2 (values (ode/multistep :heun f 0.0 [2.0] 1.0 4 {:history before :modify? true}))) [6.210093 14.88827] 1e-5)))
    (testing "Example 26.5: Milne's method, 6.204855 at 1, then 14.86031, 33.72426, 75.43295"
      (is (close? (values (ode/multistep :milne f 0.0 [2.0] 1.0 4 {:history three-before})) [6.204855 14.86031 33.72426 75.43295] 1e-5)))
    (testing "Example 26.6: the fourth-order Adams method, 6.214424 at 1"
      (is (close? (first (values (ode/multistep :adams f 0.0 [2.0] 1.0 4 {:history three-before}))) 6.214424 1e-6)))
    (testing "started by Runge-Kutta instead, at a finer step, all converge to the truth"
      (doseq [m [:heun :milne :adams]]
        (is (close? (first (second (peek (ode/multistep m f 0.0 [2.0] 0.01 400)))) (exact 4.0) (* 1e-4 (exact 4.0))) (name m))))))

;; Example 27.1's heated rod: T'' = h' (T - Ta), h' = 0.01, Ta = 20, T(0) = 40, T(10) = 200
(def ^:private hp 0.01)
(defn- rod [_ T _] (* hp (- T 20.0)))
(defn- rod-exact [x]
  (let [l (math/sqrt hp)
        ;; T - 20 = A e^(l x) + B e^(-l x) through 20 at 0 and 180 at 10
        A (/ (- 180.0 (* 20.0 (math/exp (* -10.0 l)))) (- (math/exp (* 10.0 l)) (math/exp (* -10.0 l))))
        B (- 20.0 A)]
    (+ 20.0 (* A (math/exp (* l x))) (* B (math/exp (* (- l) x))))))

(deftest boundary-values
  (testing "Example 27.1: linear shooting from slopes 10 and 20, RK4 at h = 2 -- slope 12.6907"
    (let [{:keys [slope solution]} (bvp/shoot-linear rod 0.0 10.0 40.0 200.0 [10.0 20.0] 5)]
      (is (close? slope 12.6907 1e-4))
      (is (close? (map second (subvec solution 1 5)) [65.9520 93.7481 124.5039 159.4538] 1e-3))))
  (testing "Example 27.3: finite differences on four interior nodes"
    (is (close? (map second (subvec (bvp/finite-difference-linear (constantly 0.0) (constantly hp) (constantly (* -20.0 hp)) 0.0 10.0 40.0 200.0 5) 1 5))
                [65.9698 93.7785 124.5382 159.4795] 1e-4)))
  (testing "and both converge to the true solution on fine grids"
    (is (every? (fn [[x T]] (< (abs (- T (rod-exact x))) 1e-8)) (:solution (bvp/shoot rod 0.0 10.0 40.0 200.0 [10.0 20.0] 200))))
    (is (every? (fn [[x T]] (< (abs (- T (rod-exact x))) 1e-3)) (bvp/finite-difference rod 0.0 10.0 40.0 200.0 200))))
  (testing "Example 27.2: the rod losing heat by radiation, T'' = h''(T - Ta)^4 -- shooting and finite differences agree"
    (let [F (fn [_ T _] (* 5e-8 (math/pow (- T 20.0) 4.0)))
          shot (:solution (bvp/shoot F 0.0 10.0 40.0 200.0 [0.0 10.0] 500))
          fd (bvp/finite-difference F 0.0 10.0 40.0 200.0 100)
          at (fn [sol x] (second (first (filter #(< (abs (- (first %) x)) 1e-9) sol))))]
      (doseq [x [2.0 4.0 6.0 8.0]]
        (is (close? (at shot x) (at fd x) 0.05))))))

(deftest eigenproblems
  (let [K [[2.0 -1.0 0.0] [-1.0 2.0 -1.0] [0.0 -1.0 2.0]]
        truth (sort (map #(- 2.0 (* 2.0 (math/cos (/ (* % math/PI) 4.0)))) [1 2 3]))]
    (testing "Example 27.6: the polynomial method on the column's three-node matrix"
      (is (close? (sort (map first (eigen/polynomial-method K))) truth 1e-12)))
    (testing "Examples 27.7-27.8: the power method, highest and (on the inverse) lowest"
      (let [A (mapv (fn [row] (mapv #(* 1.778 %) row)) K)]
        (is (close? (:value (eigen/power A)) (* 1.778 (last truth)) 1e-10))
        (is (close? (:value (eigen/inverse-power A)) (* 1.778 (first truth)) 1e-10))
        (testing "and shifted, the middle one"
          (is (close? (:value (eigen/inverse-power A {:shift 3.0})) (* 1.778 (second truth)) 1e-10)))))
    (testing "Jacobi's method and the QR algorithm on a symmetric 5x5"
      (let [S [[4.0 1.0 -2.0 2.0 0.5] [1.0 2.0 0.0 1.0 -1.0] [-2.0 0.0 3.0 -2.0 1.5]
               [2.0 1.0 -2.0 -1.0 0.0] [0.5 -1.0 1.5 0.0 2.5]]
            {:keys [values vectors]} (eigen/jacobi S)]
        (is (close? values (eigen/qr-algorithm S) 1e-10))
        (is (close? values (sort (map first (eigen/polynomial-method S))) 1e-8))
        (testing "A v = lambda v for each"
          (doseq [i (range 5)]
            (let [v (mapv #(nth % i) vectors)]
              (is (close? (lin/mat-vec S v) (mapv #(* (values i) %) v) 1e-10)))))
        (testing "Householder keeps the eigenvalues and leaves a tridiagonal"
          (let [T (eigen/householder S)]
            (is (every? #(< (abs %) 1e-12) (for [i (range 5) j (range 5) :when (> (abs (- i j)) 1)] (get-in T [i j]))))
            (is (close? (:values (eigen/jacobi T)) values 1e-10))))))))
