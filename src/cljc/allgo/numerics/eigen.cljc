(ns allgo.numerics.eigen
  "Eigenvalues and eigenvectors (Chapra and Canale, *Numerical Methods for
  Engineers*, 27.2): the polynomial method -- the characteristic
  polynomial det(lambda I - A), by Faddeev and LeVerrier's recurrence,
  and its roots; the power method for the largest eigenvalue, and on the
  inverse for the smallest or, shifted, the one nearest any value; and the
  'other methods' the book describes: Jacobi's rotations for a symmetric
  matrix, Householder's reduction of one to tridiagonal form, and the QR
  algorithm.

  The iterative ones are steps, as in `allgo.search`: a private step, a
  public `-seq` of its states, and the function, which runs the steps by
  `allgo.search/settle` to its tolerance."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.polynomial :as poly]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- trace [M] (reduce + (map-indexed (fn [i row] (row i)) M)))

(defn- leverrier-step
  "One step of Faddeev and LeVerrier's recurrence: M_k = A M_k-1 + c_k-1
  I, and c_k = -tr(A M_k)/k appended to the coefficients."
  [A I {:keys [k M c cs]}]
  (let [M (lin/mat-add (lin/mat-mul A M) (lin/mat-scale I c))
        c' (/ (- (trace (lin/mat-mul A M))) k)]
    {:k (inc k) :M M :c c' :cs (conj cs c')}))

(defn characteristic-polynomial
  "The coefficients of det(lambda I - A), highest power first, by the
  Faddeev-LeVerrier recurrence: M_k = A M_k-1 + c_k-1 I and c_k = -tr(A
  M_k)/k, the determinant never expanded."
  [A]
  (let [n (count A) I (lin/eye n)]
    (:cs (nth (iterate (partial leverrier-step A I) {:k 1 :M (lin/mat-scale I 0.0) :c 1.0 :cs [1.0]}) n))))

(defn polynomial-method
  "The eigenvalues of A as the roots of its characteristic polynomial
  (by Bairstow's method), `[re im]` each -- the direct way for small
  matrices, ill-conditioned for large ones, whose polynomial's roots move
  far with small changes in its coefficients."
  [A]
  (poly/bairstow (characteristic-polynomial A)))

(defn- normalize-max
  "`v` scaled so its largest component is 1, and that component's value
  -- the first of any tied in magnitude, so that an eigenvector like
  (1, 0, -1) is normalized the same way every time and its sign does not
  flip between iterations."
  [v]
  (let [big (apply max (map abs v))
        j (first (filter #(>= (abs (v %)) (* big (- 1.0 1e-9))) (range (count v))))
        m (v j)]
    [m (mapv #(/ % m) v)]))

(defn- power-step
  "x <- `apply` x, normalized by its largest component, which is the new
  estimate `:m`; the one before kept as `:last-m`."
  [apply-A {:keys [x m iterations]}]
  (let [[m' x'] (normalize-max (apply-A x))]
    {:x x' :m m' :last-m m :iterations (inc iterations)}))

(defn- power-settled?
  "Whether the power method's estimate has settled to within `tol`
  relative."
  [tol {:keys [m last-m]}]
  (and last-m (<= (abs (- m last-m)) (* tol (abs m)))))

(defn power-seq
  "The power method's estimates `{:value :vector :iterations}` from `:x0`
  (default ones): those of Chapra and Canale's example 27.7."
  ([A] (power-seq A {}))
  ([A {:keys [x0]}]
   (map (fn [{:keys [x m iterations]}] {:value m :vector x :iterations iterations})
        (rest (iterate (partial power-step (partial lin/mat-vec A)) {:x (vec (or x0 (repeat (count A) 1.0))) :m 0.0 :iterations 0})))))

(defn power
  "The eigenvalue of largest magnitude of A and its eigenvector, by the
  power method: x <- A x, normalized by its largest component, which
  converges to that eigenvalue, at a rate set by the ratio of the next
  largest to it. `{:value :vector :iterations}`, the vector's largest
  component 1."
  ([A] (power A {}))
  ([A {:keys [tol max-iter x0] :or {tol 1e-13 max-iter 10000}}]
   (let [{:keys [x m iterations]} (search/settle (partial power-settled? tol) (inc max-iter)
                                                 (partial power-step (partial lin/mat-vec A))
                                                 {:x (vec (or x0 (repeat (count A) 1.0))) :m 0.0 :iterations 0})]
     {:value m :vector x :iterations iterations})))

(defn- shifted-solver
  "x -> (A - shift I)^-1 x, by one LU decomposition."
  [A shift]
  (let [d (ls/lu (mapv (fn [i row] (update row i - shift)) (range (count A)) A))]
    (partial ls/lu-solve d)))

(defn inverse-power-seq
  "The inverse power method's estimates `{:value :vector :iterations}`,
  `:value` the eigenvalue nearest `:shift` so far: those of Chapra and
  Canale's example 27.8, as 1/value."
  ([A] (inverse-power-seq A {}))
  ([A {:keys [shift] :or {shift 0.0}}]
   (map (fn [{:keys [x m iterations]}] {:value (+ shift (/ 1.0 m)) :vector x :iterations iterations})
        (rest (iterate (partial power-step (shifted-solver A shift)) {:x (vec (repeat (count A) 1.0)) :m 0.0 :iterations 0})))))

(defn inverse-power
  "The eigenvalue of A nearest `shift` (default 0, the smallest in
  magnitude) and its eigenvector, by the power method on (A - shift
  I)^-1 -- each step a solve with one LU decomposition, not an inverse --
  whose largest eigenvalue is 1/(lambda - shift) for that lambda.
  `{:value :vector :iterations}`."
  ([A] (inverse-power A {}))
  ([A {:keys [shift tol max-iter] :or {shift 0.0 tol 1e-13 max-iter 10000}}]
   (let [{:keys [x m iterations]} (search/settle (partial power-settled? tol) (inc max-iter)
                                                 (partial power-step (shifted-solver A shift))
                                                 {:x (vec (repeat (count A) 1.0)) :m 0.0 :iterations 0})]
     {:value (+ shift (/ 1.0 m)) :vector x :iterations iterations})))

(defn- off-diagonal
  "The square root of the sum of squares of M's off-diagonal entries."
  [M]
  (let [n (count M)]
    (math/sqrt (reduce + (for [i (range n) j (range n) :when (not= i j)] (let [x (get-in M [i j])] (* x x)))))))

(defn- jacobi-rotate
  "M and V after the plane rotation that zeroes M's (p, q) entry: M' =
  R^T M R, and V' = V R."
  [[M V] [p q]]
  (let [apq (get-in M [p q])]
    (if (zero? apq)
      [M V]
      (let [theta (/ (- (get-in M [q q]) (get-in M [p p])) (* 2.0 apq))
            t (/ (if (neg? theta) -1.0 1.0) (+ (abs theta) (math/sqrt (+ 1.0 (* theta theta)))))
            c (/ 1.0 (math/sqrt (+ 1.0 (* t t)))) s (* t c)
            rot-cols (fn [X] (mapv (fn [row] (let [xp (row p) xq (row q)]
                                               (assoc row p (- (* c xp) (* s xq)) q (+ (* s xp) (* c xq)))))
                                   X))
            M1 (rot-cols M)
            M2 (let [rp (M1 p) rq (M1 q)]
                 (assoc M1 p (mapv #(- (* c %1) (* s %2)) rp rq) q (mapv #(+ (* s %1) (* c %2)) rp rq)))]
        [M2 (rot-cols V)]))))

(defn- jacobi-sweep
  "One cyclic sweep of Jacobi's rotations over every off-diagonal pair."
  [{:keys [M V sweep]}]
  (let [n (count M)
        [M V] (reduce jacobi-rotate [M V] (for [p (range n) q (range (inc p) n)] [p q]))]
    {:M M :V V :sweep (inc sweep) :off (off-diagonal M)}))

(defn jacobi-seq
  "The sweeps of Jacobi's method on the symmetric A: `{:M :V :sweep
  :off}`, M rotated toward diagonal, V the rotations so far and `:off`
  the off-diagonal part left, falling every sweep."
  [A]
  (let [M (mapv vec A)]
    (iterate jacobi-sweep {:M M :V (lin/eye (count A)) :sweep 0 :off (off-diagonal M)})))

(defn jacobi
  "Every eigenvalue and eigenvector of the symmetric matrix A by Jacobi's
  method: plane rotations, each zeroing an off-diagonal pair, swept
  cyclically until the off-diagonal part is below `:tol` (default 1e-14)
  of the whole -- each rotation undoing some earlier zeros, but the
  off-diagonal sum of squares falling every time. `{:values :vectors}`,
  the eigenvectors the columns of `:vectors`, ascending by value."
  ([A] (jacobi A {}))
  ([A {:keys [tol max-sweeps] :or {tol 1e-14 max-sweeps 100}}]
   (let [n (count A)
         norm (math/sqrt (reduce + (for [row A x row] (* x x))))
         M (mapv vec A)
         {:keys [M V]} (search/settle #(<= (:off %) (* tol (max norm 1e-300))) (inc max-sweeps)
                                      jacobi-sweep {:M M :V (lin/eye n) :sweep 0 :off (off-diagonal M)})
         order (sort-by #(get-in M [% %]) (range n))]
     {:values (mapv #(get-in M [% %]) order)
      :vectors (mapv (fn [row] (mapv row order)) V)})))

(defn- householder-vector
  "The Householder reflection I - 2 v v^T/(v^T v) that maps the vector x
  onto a multiple of the first unit vector: v."
  [x]
  (let [alpha (* (if (neg? (first x)) 1.0 -1.0) (math/sqrt (reduce + (map * x x))))]
    (update (vec x) 0 - alpha)))

(defn- reflect
  "P A P for the reflection of `v` acting on rows and columns from `k` on."
  [A v k]
  (let [n (count A)
        vv (reduce + (map * v v))
        full (vec (concat (repeat k 0.0) v))]
    (if (zero? vv)
      A
      (let [P (vec (for [i (range n)] (vec (for [j (range n)] (- (if (= i j) 1.0 0.0) (/ (* 2.0 (full i) (full j)) vv))))))]
        (lin/mat-mul (lin/mat-mul P A) P)))))

(defn householder
  "The symmetric matrix A reduced to tridiagonal form by n - 2 Householder
  reflections, each zeroing a column below the subdiagonal (and the row,
  by symmetry): the same eigenvalues, and a far cheaper matrix to find
  them in."
  [A]
  (let [n (count A)]
    (reduce (fn [A k]
              (let [x (mapv #(get-in A [% (dec k)]) (range k n))]
                (if (every? #(< (abs %) 1e-300) (rest x))
                  A
                  (reflect A (householder-vector x) k))))
            (mapv vec A) (range 1 (dec n)))))

(defn- qr-decompose
  "Q and R of A by Householder reflections."
  [A]
  (let [n (count A)]
    (loop [k 0 R (mapv vec A) Q (lin/eye n)]
      (if (>= k (dec n))
        [Q R]
        (let [x (mapv #(get-in R [% k]) (range k n))
              v (householder-vector x)
              vv (reduce + (map * v v))]
          (if (zero? vv)
            (recur (inc k) R Q)
            (let [full (vec (concat (repeat k 0.0) v))
                  H (vec (for [i (range n)] (vec (for [j (range n)] (- (if (= i j) 1.0 0.0) (/ (* 2.0 (full i) (full j)) vv))))))]
              (recur (inc k) (lin/mat-mul H R) (lin/mat-mul Q H)))))))))

(defn- deflate
  "T without its last row and column, that eigenvalue added to `out`."
  [{:keys [T out]}]
  (let [n (count T)]
    {:T (mapv #(subvec % 0 (dec n)) (subvec T 0 (dec n))) :out (conj out (get-in T [(dec n) (dec n)])) :k 0 :by :deflate}))

(defn- qr-step
  "One step of the QR algorithm on the tridiagonal T: an eigenvalue
  deflated off the bottom once its row's subdiagonal is below `tol` of
  its neighbors on the diagonal (or after `max-iter` steps on it); else T
  <- R Q + mu I from T - mu I = Q R, mu Wilkinson's shift. nil once
  every eigenvalue is out."
  [tol max-iter {:keys [T out k] :as st}]
  (let [n (count T)]
    (cond
      (= n 0) nil
      (= n 1) (deflate st)
      (or (>= k max-iter)
          (<= (abs (get-in T [(dec n) (- n 2)]))
              (* tol (+ (abs (get-in T [(dec n) (dec n)])) (abs (get-in T [(- n 2) (- n 2)]))))))
      (deflate st)
      :else
      (let [a (get-in T [(- n 2) (- n 2)]) b (get-in T [(dec n) (- n 2)]) c (get-in T [(dec n) (dec n)])
            d (* 0.5 (- a c))
            mu (- c (/ (* b b) (+ d (* (if (neg? d) -1.0 1.0) (math/sqrt (+ (* d d) (* b b)))))))
            shifted (mapv (fn [i row] (update row i - mu)) (range n) T)
            [Q R] (qr-decompose shifted)]
        {:T (mapv (fn [i row] (update row i + mu)) (range n) (lin/mat-mul R Q)) :out out :k (inc k) :by :qr}))))

(defn qr-algorithm-seq
  "The QR algorithm's states on the symmetric A: `{:T :out :by}`, T what
  is left to work on, `:out` the eigenvalues deflated so far, `:by`
  `:qr` or `:deflate`. Options as `qr-algorithm`'s."
  ([A] (qr-algorithm-seq A {}))
  ([A {:keys [tol max-iter] :or {tol 1e-14 max-iter 10000}}]
   (search/step-seq (partial qr-step tol max-iter) {:T (householder A) :out [] :k 0 :by :start})))

(defn qr-algorithm
  "The eigenvalues of the symmetric matrix A by the QR algorithm: reduced
  to tridiagonal form, then A <- R Q + mu I from A - mu I = Q R, the
  shift mu Wilkinson's (the eigenvalue of the trailing 2x2 nearer its
  corner), each eigenvalue deflated off the bottom as its row's
  subdiagonal vanishes. Ascending."
  ([A] (qr-algorithm A {}))
  ([A {:keys [tol max-iter] :or {tol 1e-14 max-iter 10000}}]
   (-> (search/settle #(empty? (:T %)) ##Inf (partial qr-step tol max-iter) {:T (householder A) :out [] :k 0 :by :start})
       :out sort vec)))
