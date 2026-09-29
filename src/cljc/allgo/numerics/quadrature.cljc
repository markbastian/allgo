(ns allgo.numerics.quadrature
  "Definite integrals by the composite Newton-Cotes rules.

  Simpson's rule fits a parabola through each pair of intervals and is
  exact for cubics, its error falling as h^4; the trapezoid rule is exact
  only for lines -- except over a whole period of a smooth periodic
  function, where its error falls faster than any power of h, which makes
  it the rule for orbit averages.

  And the rest (Chapra and Canale, *Numerical Methods for Engineers*,
  chapters 21-22): Simpson's 3/8 rule, the closed and open Newton-Cotes
  formulas to six points, integrals of unequally spaced data, double
  integrals, Romberg integration, adaptive quadrature, Gauss-Legendre
  quadrature of any order, the extended midpoint rule, and improper
  integrals by the change of variable x = 1/t."
  (:require [clojure.math :as math]))

(defn simpson
  "The integral of `f` from `a` to `b` by Simpson's rule on `n` intervals
  (made even if odd)."
  [f a b n]
  (let [n (if (odd? n) (inc n) n)
        h (/ (- b a) n)]
    (* (/ h 3.0)
       (reduce + (for [k (range (inc n))
                       :let [w (cond (or (zero? k) (= k n)) 1.0 (odd? k) 4.0 :else 2.0)]]
                   (* w (f (+ a (* k h)))))))))

(defn trapezoid
  "The integral of `f` from `a` to `b` by the trapezoid rule on `n`
  intervals."
  [f a b n]
  (let [h (/ (- b a) n)]
    (* h (reduce + (for [k (range (inc n))
                         :let [w (if (or (zero? k) (= k n)) 0.5 1.0)]]
                     (* w (f (+ a (* k h)))))))))

;; ------------------------------------------------- Newton-Cotes, continued

(defn simpson-3-8
  "The integral of `f` from `a` to `b` by Simpson's 3/8 rule on `n`
  intervals (rounded up to a multiple of three): a cubic through each
  four points, weights 3h/8 (1 3 3 1) -- the same order as the 1/3 rule,
  and the way to cover an odd number of intervals with it."
  [f a b n]
  (let [n (* 3 (long (math/ceil (/ n 3.0))))
        h (/ (- b a) n)]
    (* (/ (* 3.0 h) 8.0)
       (reduce + (for [k (range (inc n))
                       :let [w (cond (or (zero? k) (= k n)) 1.0 (zero? (mod k 3)) 2.0 :else 3.0)]]
                   (* w (f (+ a (* k h)))))))))

(def ^:private closed-weights
  "Closed Newton-Cotes formulas by number of points: the factor on h and
  the weights (Chapra and Canale table 21.2)."
  {2 [0.5 [1 1]] 3 [(/ 1.0 3.0) [1 4 1]] 4 [0.375 [1 3 3 1]]
   5 [(/ 2.0 45.0) [7 32 12 32 7]] 6 [(/ 5.0 288.0) [19 75 50 50 75 19]]})

(def ^:private open-weights
  "Open Newton-Cotes formulas by number of intervals, the ends not
  evaluated: the factor on h and the weights on the interior points
  (table 21.4)."
  {2 [2.0 [1]] 3 [1.5 [1 1]] 4 [(/ 4.0 3.0) [2 -1 2]]
   5 [(/ 5.0 24.0) [11 1 1 11]] 6 [0.3 [11 -14 26 -14 11]]})

(defn newton-cotes
  "The integral of `f` from `a` to `b` by one application of the
  Newton-Cotes formula on `points` equally spaced points, 2 to 6: closed
  (the ends among them -- trapezoid, Simpson's 1/3 and 3/8, Boole's and
  the six-point rule), or with `:open? true` open, `points` intervals
  whose interior points alone are used, for integrands that cannot be
  evaluated at an end."
  ([f a b points] (newton-cotes f a b points {}))
  ([f a b points {:keys [open?]}]
   (if open?
     (let [[c ws] (open-weights points)
           h (/ (- b a) points)]
       (* c h (reduce + (map-indexed (fn [i w] (* w (f (+ a (* (inc i) h))))) ws))))
     (let [[c ws] (closed-weights points)
           h (/ (- b a) (dec points))]
       (* c h (reduce + (map-indexed (fn [i w] (* w (f (+ a (* i h))))) ws)))))))

;; ----------------------------------------------------------------- data

(defn trapezoid-data
  "The integral of the data `xs` `ys`, any spacing, by the trapezoid
  rule on each interval."
  [xs ys]
  (reduce + (map (fn [x0 x1 y0 y1] (* 0.5 (- x1 x0) (+ y0 y1))) xs (rest xs) ys (rest ys))))

(defn simpson-data
  "The integral of unequally spaced data by the best rule each stretch
  allows, as the book's algorithm takes them from the left (Chapra and
  Canale 21.3): where three intervals in a row are equal, Simpson's 3/8
  rule; where two, the 1/3 rule; elsewhere the trapezoid."
  [xs ys]
  (let [xs (vec xs) ys (vec ys) n (dec (count xs))
        h (fn [i] (- (xs (inc i)) (xs i)))
        eq? (fn [i j] (< (abs (- (h i) (h j))) (* 1e-9 (max (abs (h i)) (abs (h j))))))]
    (loop [i 0 total 0.0]
      (cond
        (>= i n) total
        (and (<= (+ i 3) n) (eq? i (+ i 1)) (eq? i (+ i 2)))
        (recur (+ i 3) (+ total (* 0.375 (h i) (+ (ys i) (* 3.0 (ys (+ i 1))) (* 3.0 (ys (+ i 2))) (ys (+ i 3))))))
        (and (<= (+ i 2) n) (eq? i (+ i 1)))
        (recur (+ i 2) (+ total (* (/ (h i) 3.0) (+ (ys i) (* 4.0 (ys (+ i 1))) (ys (+ i 2))))))
        :else
        (recur (inc i) (+ total (* 0.5 (h i) (+ (ys i) (ys (inc i))))))))))

(defn double-integral
  "The integral of `f` of x and y over the rectangle `ax`..`bx` by
  `ay`..`by`, by Simpson's rule on `nx` and `ny` intervals, iterated: the
  inner integral in y for each x, then the outer."
  [f ax bx ay by nx ny]
  (simpson (fn [x] (simpson #(f x %) ay by ny)) ax bx nx))

;; ------------------------------------------------------------- Romberg

(defn romberg
  "The integral of `f` from `a` to `b` by Romberg integration: trapezoid
  estimates on 1, 2, 4, ... intervals, each new one reusing the old
  points, combined by Richardson's extrapolation (I_jk = (4^k I_j+1,k-1 -
  I_j,k-1)/(4^k - 1)) until two successive diagonal values agree within
  `:tol` relative (default 1e-12). `{:value :table}`."
  ([f a b] (romberg f a b {}))
  ([f a b {:keys [tol max-level] :or {tol 1e-12 max-level 20}}]
   (loop [level 1 n 1 trap (* 0.5 (- b a) (+ (f a) (f b))) rows [[(* 0.5 (- b a) (+ (f a) (f b)))]]]
     (let [n2 (* 2 n) h (/ (- b a) n2)
           trap (+ (* 0.5 trap) (* h (reduce + (for [k (range 1 n2 2)] (f (+ a (* k h)))))))
           prev (peek rows)
           row (reduce (fn [row k]
                         (let [p (math/pow 4.0 k)]
                           (conj row (/ (- (* p (row (dec k))) (prev (dec k))) (- p 1.0)))))
                       [trap] (range 1 (inc level)))
           rows (conj rows row)]
       (if (or (>= level max-level)
               (<= (abs (- (peek row) (peek prev))) (* tol (max 1e-300 (abs (peek row))))))
         {:value (peek row) :table rows}
         (recur (inc level) n2 trap rows))))))

(defn adaptive
  "The integral of `f` from `a` to `b` by adaptive quadrature (Chapra and
  Canale 22.3): Simpson's rule on an interval and on its halves, their
  difference over 15 the error estimate; kept, with that correction added,
  where the error is below the interval's share of `:tol` (default
  1e-10), halved and tried again where not -- the points clustering where
  the integrand is hard."
  ([f a b] (adaptive f a b {}))
  ([f a b {:keys [tol max-depth] :or {tol 1e-10 max-depth 50}}]
   (let [simp (fn [a fa b fb] (let [m (* 0.5 (+ a b)) fm (f m)] [m fm (* (/ (- b a) 6.0) (+ fa (* 4.0 fm) fb))]))
         step (fn step [a fa b fb m fm whole tol depth]
                (let [[lm flm left] (simp a fa m fm)
                      [rm frm right] (simp m fm b fb)
                      delta (- (+ left right) whole)]
                  (if (or (<= depth 0) (<= (abs delta) (* 15.0 tol)))
                    (+ left right (/ delta 15.0))
                    (+ (step a fa m fm lm flm left (* 0.5 tol) (dec depth))
                       (step m fm b fb rm frm right (* 0.5 tol) (dec depth))))))
         fa (f a) fb (f b)
         [m fm whole] (simp a fa b fb)]
     (step a fa b fb m fm whole tol max-depth))))

;; ------------------------------------------------------ Gauss-Legendre

(defn- legendre-nodes
  "The nodes and weights of n-point Gauss-Legendre quadrature on [-1, 1]:
  the roots of P_n by Newton's method from Tricomi's approximation, the
  weights 2 / ((1 - x^2) P_n'(x)^2)."
  [n]
  (let [p (fn [x] (loop [k 1 p0 1.0 p1 x]
                    (if (= k n) [p1 (/ (* n (- (* x p1) p0)) (- (* x x) 1.0))]
                        (recur (inc k) p1 (/ (- (* (+ (* 2 k) 1) x p1) (* k p0)) (inc k))))))]
    (vec (for [i (range 1 (inc n))]
           (let [x0 (math/cos (/ (* math/PI (- i 0.25)) (+ n 0.5)))
                 x (loop [x x0 k 0]
                     (let [[pn dp] (p x) x' (- x (/ pn dp))]
                       (if (or (< (abs (- x' x)) 1e-16) (> k 100)) x' (recur x' (inc k)))))
                 [_ dp] (p x)]
             [x (/ 2.0 (* (- 1.0 (* x x)) dp dp))])))))

(def ^:private gauss-table (memoize legendre-nodes))

(defn gauss-legendre
  "The integral of `f` from `a` to `b` by `n`-point Gauss-Legendre
  quadrature (Chapra and Canale 22.4): the points and weights chosen so
  the rule is exact for polynomials of degree 2n - 1, the interval mapped
  onto [-1, 1]. n = 1 is the midpoint rule, 2 the book's two-point
  formula at +/-1/sqrt(3)."
  [f a b n]
  (let [c (* 0.5 (+ a b)) h (* 0.5 (- b a))]
    (* h (reduce + (map (fn [[x w]] (* w (f (+ c (* h x))))) (gauss-table n))))))

;; ------------------------------------------------------------ improper

(defn extended-midpoint
  "The integral of `f` from `a` to `b` by the extended midpoint rule
  (Chapra and Canale 22.35): f at the middles of `n` equal intervals, the
  ends never evaluated."
  [f a b n]
  (let [h (/ (- b a) n)]
    (* h (reduce + (for [k (range n)] (f (+ a (* (+ k 0.5) h))))))))

(defn- midpoint-romberg
  "The extended midpoint rule on 1, 3, 9, ... intervals (tripling reuses
  every point), extrapolated as Romberg does -- the midpoint rule's error
  also runs in even powers of h."
  [f a b tol]
  (loop [n 1 prev-row [(extended-midpoint f a b 1)] k 1]
    (let [n3 (* 3 n)
          m (extended-midpoint f a b n3)
          row (reduce (fn [row j]
                        (let [p (math/pow 9.0 j)]
                          (conj row (/ (- (* p (row (dec j))) (prev-row (dec j))) (- p 1.0)))))
                      [m] (range 1 (inc (min k 6))))
          row (if (> (count row) 7) (subvec row 0 7) row)]
      (if (or (> k 12) (<= (abs (- (peek row) (peek prev-row))) (* tol (max 1e-300 (abs (peek row))))))
        (peek row)
        (recur n3 row (inc k))))))

(defn improper
  "The integral of `f` from `a` to `b`, either or both possibly infinite
  (##Inf, ##-Inf), by the change of variable x = 1/t on the infinite
  parts (Chapra and Canale 22.5): the integral from A to infinity is that
  of f(1/t)/t^2 from 0 to 1/A, for A > 0, evaluated by the extended
  midpoint rule so t = 0 is never touched, extrapolated. The finite part
  is Romberg's. To within `:tol` relative (default 1e-10)."
  ([f a b] (improper f a b {}))
  ([f a b {:keys [tol] :or {tol 1e-10}}]
   (let [tail-up (fn [A] (midpoint-romberg (fn [t] (/ (f (/ 1.0 t)) (* t t))) 0.0 (/ 1.0 A) tol))
         tail-down (fn [B] (midpoint-romberg (fn [t] (/ (f (/ 1.0 t)) (* t t))) (/ 1.0 B) 0.0 tol))
         finite (fn [lo hi] (if (= lo hi) 0.0 (:value (romberg f lo hi {:tol tol}))))
         inf? #(or (= % ##Inf) (= % ##-Inf))]
     (cond
       (and (not (inf? a)) (not (inf? b))) (finite a b)
       (and (not (inf? a)) (= b ##Inf)) (let [A (max a 1.0)] (+ (finite a A) (tail-up A)))
       (and (= a ##-Inf) (not (inf? b))) (let [B (min b -1.0)] (+ (tail-down B) (finite B b)))
       :else (+ (tail-down -1.0) (finite -1.0 1.0) (tail-up 1.0))))))
