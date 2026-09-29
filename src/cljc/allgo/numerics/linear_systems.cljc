(ns allgo.numerics.linear-systems
  "Linear algebraic equations, A x = b, the direct and iterative ways
  (Chapra and Canale, *Numerical Methods for Engineers*, chapters 9-11):
  Cramer's rule, Gauss elimination with partial pivoting and scaling,
  Gauss-Jordan, LU decomposition (Doolittle with pivoting, and Crout),
  determinants and inverses from it, vector and matrix norms and the
  condition number, iterative refinement, Gauss-Seidel with relaxation,
  Jacobi iteration, and complex systems. The tridiagonal (Thomas) solver
  is `allgo.numerics.tridiagonal`, Cholesky `allgo.numerics.linear`.

  Matrices are vectors of row vectors."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

;; ------------------------------------------------------ small systems

(declare lu det)

(defn cramer
  "The solution of A x = b by Cramer's rule: each x_i the determinant of
  A with its column i replaced by b, over det A. For a handful of
  equations only -- the determinants cost n! the obvious way, n^3 by LU."
  [A b]
  (let [d (det A)]
    (when-not (zero? d)
      (mapv (fn [i] (/ (det (mapv #(assoc %1 i %2) A b)) d)) (range (count b))))))

;; --------------------------------------------------- Gauss elimination

(defn gauss
  "The solution of A x = b by Gauss elimination -- forward elimination to
  an upper triangle, then back substitution -- with partial pivoting: at
  each column the row with the largest pivot, scaled by its own largest
  coefficient when `:scale?` (default true), is swapped up, so that no
  pivot is zero and none tiny where a larger one was to hand. nil if A is
  singular (a pivot below `:tol`, default 1e-14 of the row scale)."
  ([A b] (gauss A b {}))
  ([A b {:keys [scale? tol] :or {scale? true tol 1e-14}}]
   (let [n (count b)
         s (mapv #(if scale? (apply max (map abs %)) 1.0) A)]
     (when (every? pos? s)
       (loop [k 0 A (mapv vec A) b (vec b) s s]
         (if (= k n)
           ;; back substitution
           (reduce (fn [x i]
                     (assoc x i (/ (- (b i) (reduce + (map #(* (get-in A [i %]) (x %)) (range (inc i) n))))
                                   (get-in A [i i]))))
                   (vec (repeat n 0.0))
                   (range (dec n) -1 -1))
           (let [p (apply max-key #(/ (abs (get-in A [% k])) (s %)) (range k n))]
             (when (> (/ (abs (get-in A [p k])) (s p)) tol)
               (let [swap (fn [v] (assoc v k (v p) p (v k)))
                     A (swap A) b (swap b) s (swap s)
                     pivot (get-in A [k k])
                     [A b] (reduce (fn [[A b] i]
                                     (let [f (/ (get-in A [i k]) pivot)]
                                       [(assoc A i (mapv - (A i) (mapv #(* f %) (A k))))
                                        (assoc b i (- (b i) (* f (b k))))]))
                                   [A b]
                                   (range (inc k) n))]
                 (recur (inc k) A b s))))))))))

(defn gauss-jordan
  "The solution of A x = b by Gauss-Jordan: every unknown eliminated from
  every row but its own, the pivot rows normalized, so the right side
  becomes the answer -- no back substitution, at half again the work of
  Gauss elimination. Partial pivoting; nil if singular."
  [A b]
  (let [n (count b)]
    (loop [k 0 M (mapv (fn [row bi] (conj (vec row) bi)) A b)]
      (if (= k n)
        (mapv peek M)
        (let [p (apply max-key #(abs (get-in M [% k])) (range k n))]
          (when-not (zero? (get-in M [p k]))
            (let [M (assoc M k (M p) p (M k))
                  row (mapv #(/ % (get-in M [k k])) (M k))
                  M (assoc M k row)]
              (recur (inc k)
                     (mapv (fn [i r] (if (= i k) r (let [f (r k)] (mapv - r (mapv #(* f %) row)))))
                           (range n) M)))))))))

;; ------------------------------------------------------------------ LU

(defn lu
  "The LU decomposition of A with partial pivoting (Doolittle: L unit
  lower triangular): `{:lu :perm :sign}` -- L below the diagonal and U on
  and above it in one matrix, the row order, and the permutation's sign.
  Decomposed once, A x = b is then two triangular substitutions for any
  b. nil if A is singular."
  [A]
  (let [n (count A)]
    (loop [k 0 M (mapv vec A) perm (vec (range n)) sign 1.0]
      (if (= k n)
        {:lu M :perm perm :sign sign}
        (let [p (apply max-key #(abs (get-in M [% k])) (range k n))]
          (when-not (zero? (get-in M [p k]))
            (let [swapped? (not= p k)
                  M (assoc M k (M p) p (M k))
                  perm (assoc perm k (perm p) p (perm k))
                  pivot (get-in M [k k])
                  M (reduce (fn [M i]
                              (let [f (/ (get-in M [i k]) pivot)
                                    row (M i)]
                                (assoc M i (into (conj (subvec row 0 k) f)
                                                 (map #(- (row %) (* f (get-in M [k %]))) (range (inc k) n))))))
                            M (range (inc k) n))]
              (recur (inc k) M perm (if swapped? (- sign) sign)))))))))

(defn lu-solve
  "The solution of A x = b from `(lu A)`: forward substitution through L
  (L d = P b), then back substitution through U (U x = d)."
  [{:keys [lu perm]} b]
  (let [n (count b)
        pb (mapv b perm)
        d (reduce (fn [d i] (assoc d i (- (pb i) (reduce + (map #(* (get-in lu [i %]) (d %)) (range i))))))
                  (vec (repeat n 0.0)) (range n))]
    (reduce (fn [x i]
              (assoc x i (/ (- (d i) (reduce + (map #(* (get-in lu [i %]) (x %)) (range (inc i) n))))
                            (get-in lu [i i]))))
            (vec (repeat n 0.0)) (range (dec n) -1 -1))))

(defn crout
  "The Crout decomposition of A, `{:L :U}`: L lower triangular, U unit
  upper triangular -- the other way of placing the ones -- computed
  column of L then row of U in turn, without pivoting. nil at a zero
  pivot."
  [A]
  (let [n (count A)
        zero (vec (repeat n (vec (repeat n 0.0))))]
    (loop [j 0 L zero U (mapv #(assoc % %2 1.0) zero (range n))]
      (if (= j n)
        {:L L :U U}
        (let [L (reduce (fn [L i]
                          (assoc-in L [i j] (- (get-in A [i j]) (reduce + (map #(* (get-in L [i %]) (get-in U [% j])) (range j))))))
                        L (range j n))
              ljj (get-in L [j j])]
          (when-not (zero? ljj)
            (recur (inc j) L
                   (reduce (fn [U k]
                             (assoc-in U [j k] (/ (- (get-in A [j k]) (reduce + (map #(* (get-in L [j %]) (get-in U [% k])) (range j))))
                                                  ljj)))
                           U (range (inc j) n)))))))))

(defn det
  "The determinant of A, the product of U's diagonal from the LU
  decomposition, its sign flipped for each row swap."
  [A]
  (if-let [{:keys [lu sign]} (lu A)]
    (* sign (reduce * (map #(get-in lu [% %]) (range (count A)))))
    0.0))

(defn inverse
  "A^-1, column by column from one LU decomposition: each column the
  solution for a column of the identity. nil if A is singular."
  [A]
  (when-let [d (lu A)]
    (let [n (count A)]
      (lin/transpose (mapv (fn [j] (lu-solve d (mapv #(if (= % j) 1.0 0.0) (range n)))) (range n))))))

;; --------------------------------------------------------- conditioning

(defn vector-norm
  "The norm of `x`: `:euclidean` (default), `:one` (the sum of
  magnitudes) or `:max` (the largest)."
  ([x] (vector-norm x :euclidean))
  ([x kind]
   (case kind
     :euclidean (math/sqrt (reduce + (map #(* % %) x)))
     :one (reduce + (map abs x))
     :max (apply max (map abs x)))))

(defn matrix-norm
  "The norm of A: `:frobenius` (default), `:column-sum` (the largest
  column's magnitudes summed, the one-norm) or `:row-sum` (the largest
  row's, the uniform norm)."
  ([A] (matrix-norm A :frobenius))
  ([A kind]
   (case kind
     :frobenius (math/sqrt (reduce + (for [row A x row] (* x x))))
     :column-sum (apply max (map #(reduce + (map abs %)) (lin/transpose A)))
     :row-sum (apply max (map #(reduce + (map abs %)) A)))))

(defn condition-number
  "Cond A = ||A|| ||A^-1||, in the norm `kind` (default `:row-sum`): the
  factor by which A's relative errors may grow in x -- log10 of it the
  digits of precision a solution may lose. ##Inf if A is singular."
  ([A] (condition-number A :row-sum))
  ([A kind]
   (if-let [inv (inverse A)]
     (* (matrix-norm A kind) (matrix-norm inv kind))
     ##Inf)))

(defn refine
  "`x`, a solution of A x = b, improved by iterative refinement: the
  residual b - A x taken, the correction for it solved with the same
  decomposition, and added; `passes` times (default 1)."
  ([A b x] (refine A b x 1))
  ([A b x passes]
   (let [d (lu A)]
     (nth (iterate (fn [x] (mapv + x (lu-solve d (mapv - b (lin/mat-vec A x))))) x) passes))))

;; ----------------------------------------------------- iterative methods

(defn gauss-seidel
  "The solution of A x = b by Gauss-Seidel: each equation solved for its
  diagonal unknown in turn, each new value used at once, the new value
  then blended with the old by the relaxation factor `:lambda` (default
  1; between 1 and 2 over-relaxes, speeding convergence -- successive
  over-relaxation -- and under 1 under-relaxes, damping oscillation).
  Converges when A is diagonally dominant. From `:x0` (default zeros)
  until the largest change is below `:tol` relative (default 1e-12);
  `{:x :iterations}`, or nil after `:max-iter` (1000)."
  ([A b] (gauss-seidel A b {}))
  ([A b {:keys [lambda tol max-iter x0] :or {lambda 1.0 tol 1e-12 max-iter 1000}}]
   (let [n (count b)]
     (loop [x (vec (or x0 (repeat n 0.0))) k 1]
       (when (<= k max-iter)
         (let [[x' change]
               (reduce (fn [[x change] i]
                         (let [row (A i)
                               new (/ (- (b i) (reduce + (map #(if (= % i) 0.0 (* (row %) (x %))) (range n))))
                                      (row i))
                               new (+ (* lambda new) (* (- 1.0 lambda) (x i)))]
                           [(assoc x i new)
                            (max change (if (zero? new) (abs (- new (x i))) (abs (/ (- new (x i)) new))))]))
                       [x 0.0] (range n))]
           (if (< change tol)
             {:x x' :iterations k}
             (recur x' (inc k)))))))))

(defn jacobi
  "The solution of A x = b by Jacobi iteration: like Gauss-Seidel but
  every new value from the previous sweep's, so the sweep could run in
  parallel -- at the cost of converging about half as fast. Options and
  result as `gauss-seidel`'s (without relaxation)."
  ([A b] (jacobi A b {}))
  ([A b {:keys [tol max-iter x0] :or {tol 1e-12 max-iter 2000}}]
   (let [n (count b)]
     (loop [x (vec (or x0 (repeat n 0.0))) k 1]
       (when (<= k max-iter)
         (let [x' (mapv (fn [i]
                          (let [row (A i)]
                            (/ (- (b i) (reduce + (map #(if (= % i) 0.0 (* (row %) (x %))) (range n)))) (row i))))
                        (range n))
               change (apply max (map (fn [a b] (if (zero? b) (abs (- b a)) (abs (/ (- b a) b)))) x x'))]
           (if (< change tol)
             {:x x' :iterations k}
             (recur x' (inc k)))))))))

;; ------------------------------------------------------ complex systems

(defn solve-complex
  "The solution of the complex system (Ar + i Ai)(xr + i xi) = br + i bi,
  as `[xr xi]`: rewritten as the real system of twice the size,
  [[Ar -Ai] [Ai Ar]] [xr xi] = [br bi], and solved by Gauss elimination."
  [Ar Ai br bi]
  (let [n (count br)
        top (mapv (fn [r i] (into (vec r) (map - i))) Ar Ai)
        bottom (mapv (fn [r i] (into (vec i) r)) Ar Ai)]
    (when-let [x (gauss (into top bottom) (into (vec br) bi))]
      [(subvec x 0 n) (subvec x n)])))
