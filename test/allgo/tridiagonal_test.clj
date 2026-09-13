(ns allgo.tridiagonal-test
  (:require [allgo.numerics.tridiagonal :as tri]
            [clojure.test :refer [deftest is testing]]))

(defn- apply-tri
  "Multiplies the tridiagonal matrix by `x`, so a solution can be checked
  against the system it claims to solve rather than against an answer
  worked out the same way twice."
  [sub diag super x]
  (let [n (count diag)]
    (mapv (fn [i]
            (+ (* (nth diag i) (nth x i))
               (if (pos? i) (* (nth sub i) (nth x (dec i))) 0.0)
               (if (< i (dec n)) (* (nth super i) (nth x (inc i))) 0.0)))
          (range n))))

(defn- residual [sub diag super rhs x]
  (apply max (map #(abs (- (double %1) (double %2))) (apply-tri sub diag super x) rhs)))

(deftest solves-test
  (testing "a small system worked out by hand"
    (is (= [1.0 1.0 1.0 1.0]
           (tri/solve-vectors [0 1 1 1] [4 4 4 4] [1 1 1 0] [5 6 6 5]))))

  (testing "the identity gives back the right hand side"
    (is (= [3.0 -1.0 7.0]
           (tri/solve-vectors [0 0 0] [1 1 1] [0 0 0] [3 -1 7]))))

  (testing "a single equation"
    (is (= [2.5] (tri/solve-vectors [0] [2] [0] [5]))))

  (testing "the second difference, against the system itself"
    ;; The matrix every implicit one-dimensional scheme produces.
    (doseq [n [4 16 129]]
      (let [sub (vec (repeat n 1.0))
            diag (vec (repeat n -2.0))
            super (vec (repeat n 1.0))
            rhs (mapv #(Math/sin (* 0.1 %)) (range n))
            x (tri/solve-vectors sub diag super rhs)]
        (is (some? x))
        (is (< (residual sub diag super rhs x) 1e-9)))))

  (testing "an asymmetric, strongly dominant system"
    (let [n 50
          sub (vec (repeat n 0.3))
          diag (mapv #(+ 5.0 (* 0.1 %)) (range n))
          super (vec (repeat n -0.7))
          rhs (mapv #(Math/cos (* 0.37 %)) (range n))
          x (tri/solve-vectors sub diag super rhs)]
      (is (< (residual sub diag super rhs x) 1e-9)))))

(deftest singular-test
  (testing "a vanishing pivot is reported, not divided by"
    (is (nil? (tri/solve-vectors [0 1] [0 1] [1 0] [1 1])))
    ;; Rows that sum to zero: the classic singular second difference with
    ;; no boundary condition, which is exactly what the m = 0 mode of a
    ;; spherical Poisson solve looks like before it is pinned.
    (is (nil? (tri/solve-vectors [0 1 1] [-1 -2 -1] [1 1 0] [1 2 3])))))

(deftest scratch-reuse-test
  (testing "solving repeatedly through one scratch buffer gives the same
            answers as solving each on its own"
    (let [n 20
          sub (double-array (repeat n 1.0))
          diag (double-array (repeat n -2.5))
          super (double-array (repeat n 1.0))
          out (double-array n)
          scratch (double-array n)
          runs (for [k (range 5)]
                 (let [rhs (double-array (map #(Math/sin (* 0.1 (+ % k))) (range n)))]
                   (tri/solve sub diag super rhs out scratch)
                   (vec out)))
          separately (for [k (range 5)]
                       (tri/solve-vectors (repeat n 1.0) (repeat n -2.5) (repeat n 1.0)
                                          (map #(Math/sin (* 0.1 (+ % k))) (range n))))]
      (is (= (vec runs) (vec separately))))))
