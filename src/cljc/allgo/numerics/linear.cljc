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
  (:require [allgo.array :as a]
            [clojure.math :as math]))

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
(defn distance-squared [a b] (reduce + (map (fn [x y] (let [d (- x y)] (* d d))) a b)))
(defn distance [a b] (math/sqrt (distance-squared a b)))
(defn mean
  "The mean of the vectors `vs`."
  [vs]
  (let [n (count vs)] (mapv #(/ % n) (apply mapv + vs))))

(defn normalize
  "A unit vector in the same direction, or `v` itself if it has none."
  [v]
  (let [l (length v)]
    (if (zero? l) v (scale v (/ 1.0 l)))))

;; ---------------------------------------------------------------- matrices

(defn transpose [m] (apply mapv vector m))

(defn- row-dot
  "The dot product of two rows of length `n`, in a loop."
  ^double [a b n]
  (let [n (long n)]
    (loop [i 0 acc 0.0]
      (if (== i n) acc (recur (inc i) (+ acc (* (double (nth a i)) (double (nth b i)))))))))

(defn mat-mul
  "A B, in loops: the product every layer of a network makes, and most of
  the time one takes to train."
  [a b]
  (let [bt (transpose b) n (count (first a))]
    (mapv (fn [row] (mapv #(row-dot row % n) bt)) a)))

(defn mat-vec [m v] (mapv #(dot % v) m))

(defn mat-add [a b] (mapv (fn [r s] (mapv + r s)) a b))
(defn mat-sub [a b] (mapv (fn [r s] (mapv - r s)) a b))
(defn mat-scale [m s] (mapv (fn [r] (mapv #(* % s) r)) m))
(defn eye [n] (mapv (fn [i] (mapv #(if (= i %) 1.0 0.0) (range n))) (range n)))

(defn flat-cholesky!
  "The Cholesky factor L, L L^T = A, of the symmetric positive definite
  `n` by `n` matrix in `a` -- a flat double array, row by row, of which
  only the lower triangle is read -- written over that lower triangle,
  and returned; nil if A is not positive definite.

  The array form is the one to use where speed matters: on vectors of
  rows the joints' and the root's factors of `allgo.physics.articulated`
  were a third of a motorcycle step in a browser. `cholesky` is this on
  vectors."
  [^doubles a n]
  (let [n (long n)]
    (loop [i 0]
      (if (== i n)
        a
        (when (loop [j 0]
                (if (> j i)
                  true
                  (let [s (if (zero? j)
                            0.0
                            (loop [k 1 s (* (aget a (* n i)) (aget a (* n j)))]
                              (if (== k j) s (recur (inc k) (+ s (* (aget a (+ (* n i) k)) (aget a (+ (* n j) k))))))))
                        v (- (aget a (+ (* n i) j)) s)]
                    (if (== i j)
                      (when (> v 0.0) (aset a (+ (* n i) i) (math/sqrt v)) (recur (inc j)))
                      (do (aset a (+ (* n i) j) (/ v (aget a (+ (* n j) j)))) (recur (inc j)))))))
          (recur (inc i)))))))

(defn flat-cholesky-solve
  "x of A x = `b`, through `flat-cholesky!`'s factor `l` of the `n` by `n`
  A: forward substitution through L, then back through L^T. A new array."
  ^doubles [^doubles l n ^doubles b]
  (let [n (long n) y (a/f64 n) x (a/f64 n)]
    (dotimes [i n]
      (let [s (if (zero? i)
                0.0
                (loop [k 1 s (* (aget l (* n i)) (aget y 0))]
                  (if (== k i) s (recur (inc k) (+ s (* (aget l (+ (* n i) k)) (aget y k)))))))]
        (aset y i (/ (- (aget b i) s) (aget l (+ (* n i) i))))))
    (loop [i (dec n)]
      (when (>= i 0)
        (let [s (if (== i (dec n))
                  0.0
                  (loop [k (+ i 2) s (* (aget l (+ (* n (inc i)) i)) (aget x (inc i)))]
                    (if (== k n) s (recur (inc k) (+ s (* (aget l (+ (* n k) i)) (aget x k)))))))]
          (aset x i (/ (- (aget y i) s) (aget l (+ (* n i) i))))
          (recur (dec i)))))
    x))

(defn flat-spd-inverse
  "The inverse of the symmetric positive definite `n` by `n` matrix in the
  flat array `a`, as a new flat array, a column solved for each column of
  the identity; nil if it is not positive definite. `a` is left alone."
  ^doubles [^doubles a n]
  (let [n (long n)]
    (when-let [l (flat-cholesky! (a/f64 (seq a)) n)]
      (let [out (a/f64 (* n n))]
        (dotimes [c n]
          (let [e (a/f64 n)]
            (aset e c 1.0)
            (let [^doubles x (flat-cholesky-solve l n e)]
              (dotimes [r n] (aset out (+ (* n r) c) (aget x r))))))
        out))))

(defn- flat [A] (a/f64 (apply concat A)))

(defn- rows [^doubles m n lower?]
  (let [n (long n)]
    (mapv (fn [i] (mapv (fn [j] (if (and lower? (> j i)) 0.0 (aget m (+ (* n i) j)))) (range n))) (range n))))

(defn cholesky
  "Lower-triangular L with L L^T = A, for symmetric positive definite A,
  on vectors of rows: `flat-cholesky!` on a copy.

  Returns nil when A is not positive definite. For a normal matrix that is
  not a numerical mishap but a statement about the data: some direction of
  the state is not determined by the observations, and no amount of solving
  will invent it."
  [A]
  (let [n (count A)]
    (some-> (flat-cholesky! (flat A) n) (rows n true))))

(defn cholesky-solve
  "Solve A x = b for symmetric positive definite A, by forward then back
  substitution through the Cholesky factor: `flat-cholesky-solve` on
  vectors. nil if A is not positive definite."
  [A b]
  (let [n (count A)]
    (when-let [l (flat-cholesky! (flat A) n)]
      (vec (flat-cholesky-solve l n (a/f64 b))))))

(defn spd-inverse
  "Inverse of a symmetric positive definite matrix, on vectors of rows:
  `flat-spd-inverse`. nil if it is not one. Any invertible matrix is
  `allgo.numerics.linear-systems/inverse`."
  [A]
  (let [n (count A)]
    (some-> (flat-spd-inverse (flat A) n) (rows n false))))

(defn det-3
  "The determinant of a 3x3 matrix, written out -- the closed form for the
  one size that comes up everywhere; any size is
  `allgo.numerics.linear-systems/det`."
  [[[a b cc] [d e f] [g h i]]]
  (- (+ (* a e i) (* b f g) (* cc d h)) (+ (* cc e g) (* b d i) (* a f h))))

(defn solve-3
  "x of the 3x3 linear system m x = b, by Cramer's rule; any size is
  `allgo.numerics.linear-systems/gauss`."
  [m b]
  (let [d (det-3 m)
        col (fn [k] (mapv (fn [row bi] (assoc row k bi)) m b))]
    (mapv #(/ (det-3 (col %)) d) (range 3))))

(defn congruence
  "A B A^T: a matrix `B` carried through the linear map `A`, as a
  covariance is."
  [A B]
  (mat-mul (mat-mul A B) (transpose A)))
