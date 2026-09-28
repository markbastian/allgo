(ns allgo.numerics.linear
  "Small dense linear algebra: matrices as vectors of row vectors.

  Enough to solve a system and no more. The matrices here are the small
  ones that turn up inside a method -- a six by six damped least squares
  for an arm's joint rates, the normal equations of an orbit fit -- not
  the large sparse ones a field solver works on, which live in flat arrays
  and are solved by iteration rather than factored.

  `cholesky` is the workhorse, because the matrices that arise this way
  are symmetric and positive definite when the problem is well posed, and
  it is the factorization that says so: it returns nil exactly when they
  are not, which is a fact about the problem rather than a numerical
  mishap."
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------- vectors
;;
;; Of any length: a state, a gradient, a row. `allgo.geometry.vec3` is the
;; same vocabulary unrolled for exactly three, and is what to reach for
;; when a vector is a position; these are for when it is not.

(defn add [a b] (mapv + a b))
(defn sub [a b] (mapv - a b))
(defn scale [v s] (mapv #(* % s) v))
(defn dot [a b] (reduce + (map * a b)))
(defn length-squared [v] (dot v v))
(defn length [v] (math/sqrt (dot v v)))
(defn distance-squared [a b] (length-squared (sub a b)))
(defn distance [a b] (length (sub a b)))

(defn normalize
  "A unit vector in the same direction, or `v` itself if it has none."
  [v]
  (let [l (length v)]
    (if (zero? l) v (scale v (/ 1.0 l)))))

;; ---------------------------------------------------------------- matrices

(defn transpose [m] (apply mapv vector m))

(defn mat-mul [a b]
  (let [bt (transpose b)]
    (mapv (fn [row] (mapv #(dot row %) bt)) a)))

(defn mat-vec [m v] (mapv #(dot % v) m))

(defn mat-add [a b] (mapv (fn [r s] (mapv + r s)) a b))
(defn mat-sub [a b] (mapv (fn [r s] (mapv - r s)) a b))
(defn mat-scale [m s] (mapv (fn [r] (mapv #(* % s) r)) m))
(defn eye [n] (mapv (fn [i] (mapv #(if (= i %) 1.0 0.0) (range n))) (range n)))

(defn cholesky
  "Lower-triangular L with L L^T = A, for symmetric positive definite A.

  Returns nil when A is not positive definite. For a normal matrix that is
  not a numerical mishap but a statement about the data: some direction of
  the state is not determined by the observations, and no amount of solving
  will invent it."
  [A]
  (let [n (count A)]
    (loop [i 0 L (vec (repeat n (vec (repeat n 0.0))))]
      (cond
        (nil? L) nil
        (= i n)  L
        :else
        (recur (inc i)
               (loop [j 0 L L]
                 (cond
                   (nil? L) nil
                   (> j i)  L
                   :else
                   (let [s (reduce + (map * (take j (nth L i)) (take j (nth L j))))
                         v (- (nth (nth A i) j) s)]
                     (if (= i j)
                       (if (<= v 0.0)
                         nil
                         (recur (inc j) (assoc-in L [i j] (math/sqrt v))))
                       (recur (inc j) (assoc-in L [i j] (/ v (nth (nth L j) j)))))))))))))

(defn cholesky-solve
  "Solve A x = b for symmetric positive definite A, by forward then back
  substitution through the Cholesky factor."
  [A b]
  (when-let [L (cholesky A)]
    (let [n (count A)
          y (reduce (fn [y i]
                      (conj y (/ (- (nth b i) (reduce + (map * (take i (nth L i)) y)))
                                 (nth (nth L i) i))))
                    [] (range n))
          x (reduce (fn [x i]
                      (let [s (reduce + (map (fn [k] (* (nth (nth L k) i) (get x k 0.0)))
                                             (range (inc i) n)))]
                        (assoc x i (/ (- (nth y i) s) (nth (nth L i) i)))))
                    (vec (repeat n 0.0)) (reverse (range n)))]
      x)))

(defn inverse
  "Inverse of a symmetric positive definite matrix, column by column."
  [A]
  (let [n (count A)]
    (when (cholesky A)
      (transpose (mapv (fn [i] (cholesky-solve A (mapv #(if (= i %) 1.0 0.0) (range n))))
                       (range n))))))

(defn det-3 [[[a b cc] [d e f] [g h i]]]
  (- (+ (* a e i) (* b f g) (* cc d h)) (+ (* cc e g) (* b d i) (* a f h))))

(defn solve-3
  "x of the 3x3 linear system m x = b, by Cramer's rule."
  [m b]
  (let [det (fn [[[a b c] [d e f] [g h i]]]
              (- (+ (* a (- (* e i) (* f h))) (* c (- (* d h) (* e g))))
                 (* b (- (* d i) (* f g)))))
        d (det m)
        col (fn [k] (mapv (fn [row bi] (assoc row k bi)) m b))]
    (mapv #(/ (det (col %)) d) (range 3))))

(defn congruence
  "A B A^T: a matrix `B` carried through the linear map `A`, as a
  covariance is."
  [A B]
  (mat-mul (mat-mul A B) (transpose A)))

(defn inverse-general
  "The inverse of any nonsingular square matrix, by Gauss-Jordan
  elimination with partial pivoting -- for when `inverse`'s Cholesky, which
  needs symmetry and positive definiteness, does not apply. nil if the
  matrix is singular."
  [A]
  (let [n (count A)
        aug (mapv (fn [i row] (into (vec row) (map #(if (= i %) 1.0 0.0) (range n)))) (range n) A)]
    (loop [m aug col 0]
      (if (= col n)
        (mapv #(subvec % n) m)
        (let [pivot (apply max-key #(abs (get-in m [% col])) (range col n))
              pv (get-in m [pivot col])]
          (when-not (zero? pv)
            (let [m (assoc m pivot (m col) col (mapv #(/ % pv) (m pivot)))
                  prow (m col)
                  m (vec (map-indexed (fn [i row]
                                        (if (= i col) row
                                            (let [f (row col)] (mapv - row (map #(* f %) prow)))))
                                      m))]
              (recur m (inc col)))))))))
