(ns allgo.chapra-linear-test
  "Linear algebraic equations (Chapra and Canale, *Numerical Methods for
  Engineers*, 7th ed., chapters 9-11): the book's examples, and the
  identities each decomposition must satisfy."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(def ^:private A [[3.0 -0.1 -0.2] [0.1 7.0 -0.3] [0.3 -0.2 10.0]])
(def ^:private b [7.85 -19.3 71.4])
(def ^:private x [3.0 -2.5 7.0])

(deftest direct
  (testing "Examples 9.5, 9.12 and 10.2: Gauss, Gauss-Jordan, LU and Cramer on the book's system"
    (is (close? (ls/gauss A b) x 1e-12))
    (is (close? (ls/gauss-jordan A b) x 1e-12))
    (is (close? (ls/lu-solve (ls/lu A) b) x 1e-12))
    (is (close? (ls/cramer A b) x 1e-12)))
  (testing "Example 9.9: the tiny first pivot, handled by pivoting"
    (is (close? (ls/gauss [[0.0003 3.0] [1.0 1.0]] [2.0001 1.0]) [(/ 1.0 3.0) (/ 2.0 3.0)] 1e-12)))
  (testing "a zero first pivot, and a singular system"
    (is (close? (ls/gauss [[0.0 2.0] [1.0 1.0]] [4.0 3.0]) [1.0 2.0] 1e-15))
    (is (nil? (ls/gauss [[1.0 2.0] [2.0 4.0]] [1.0 2.0])))))

(deftest decompositions
  (let [M [[2.0 -1.0 4.0 1.0] [4.0 1.0 -2.0 3.0] [-2.0 5.0 1.0 1.0] [1.0 3.0 2.0 -4.0]]]
    (testing "L U reproduces the rows of A in the pivoted order"
      (let [{:keys [lu perm]} (ls/lu M)
            n 4
            L (vec (for [i (range n)] (vec (for [j (range n)] (cond (= i j) 1.0 (> i j) (get-in lu [i j]) :else 0.0)))))
            U (vec (for [i (range n)] (vec (for [j (range n)] (if (<= i j) (get-in lu [i j]) 0.0)))))]
        (is (close? (lin/mat-mul L U) (mapv M perm) 1e-12))))
    (testing "Crout: L U = A, U with a unit diagonal"
      (let [{:keys [L U]} (ls/crout A)]
        (is (close? (lin/mat-mul L U) A 1e-12))
        (is (= [1.0 1.0 1.0] (mapv #(get-in U [% %]) (range 3))))))
    (testing "the inverse, and the determinant against cofactor expansion"
      (is (close? (lin/mat-mul M (ls/inverse M)) (lin/eye 4) 1e-12))
      (is (close? (ls/det A) (lin/det-3 A) 1e-12))
      (is (close? (ls/det [[0.0 1.0] [1.0 0.0]]) -1.0 0.0)))))

(deftest conditioning
  (testing "Example 10.4: the scaled 3x3 Hilbert matrix, row-sum condition number 451.2"
    (let [H [[1.0 0.5 (/ 1.0 3.0)] [1.0 (/ 2.0 3.0) 0.5] [1.0 0.75 0.6]]]
      (is (close? (ls/matrix-norm H :row-sum) 2.35 1e-12))
      (is (close? (ls/condition-number H) 451.2 1e-9))))
  (testing "norms"
    (is (close? (ls/vector-norm [3.0 -4.0]) 5.0 0.0))
    (is (close? (ls/vector-norm [3.0 -4.0] :one) 7.0 0.0))
    (is (close? (ls/vector-norm [3.0 -4.0] :max) 4.0 0.0))
    (is (close? (ls/matrix-norm [[1.0 -2.0] [3.0 4.0]] :column-sum) 6.0 0.0))
    (is (close? (ls/matrix-norm [[1.0 -2.0] [3.0 4.0]]) (Math/sqrt 30.0) 1e-15)))
  (testing "iterative refinement recovers a perturbed solution"
    (is (close? (ls/refine A b [3.01 -2.49 7.02]) x 1e-12))))

(deftest iterative
  (testing "Example 11.3: Gauss-Seidel on the book's system -- its first iterate, then the answer"
    (is (close? (:x (ls/gauss-seidel A b {:max-iter 1 :tol ##Inf})) [2.616667 -2.794524 7.005610] 1e-6))
    (is (close? (:x (ls/gauss-seidel A b)) x 1e-11)))
  (testing "Jacobi converges too, in more sweeps"
    (let [gs (ls/gauss-seidel A b) j (ls/jacobi A b)]
      (is (close? (:x j) x 1e-11))
      (is (> (:iterations j) (:iterations gs)))))
  (testing "over-relaxation speeds up the Laplace stencil (-1, 2, -1), where plain Gauss-Seidel crawls"
    (let [n 20
          T (vec (for [i (range n)] (vec (for [j (range n)] (cond (= i j) 2.0 (= 1 (abs (- i j))) -1.0 :else 0.0)))))
          r (vec (repeat n 1.0))
          plain (ls/gauss-seidel T r {:max-iter 5000})
          sor (ls/gauss-seidel T r {:lambda 1.7})]
      (is (close? (:x sor) (ls/gauss T r) 1e-9))
      (is (< (* 5 (:iterations sor)) (:iterations plain))))))

(deftest complex-systems
  (testing "(1 + i) z1 + 2 z2 = 5 + 3i, z1 - i z2 = 2 - 2i: z = (1, 2 + i)"
    (let [[xr xi] (ls/solve-complex [[1.0 2.0] [1.0 0.0]] [[1.0 0.0] [0.0 -1.0]] [5.0 2.0] [3.0 -2.0])]
      (is (close? xr [1.0 2.0] 1e-12))
      (is (close? xi [0.0 1.0] 1e-12)))))
