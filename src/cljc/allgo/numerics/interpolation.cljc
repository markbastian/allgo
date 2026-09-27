(ns allgo.numerics.interpolation
  "Interpolation in a table of equally spaced values, and Lagrange's
  formula for unequal ones (Meeus, *Astronomical Algorithms*, chapter 3).

  An almanac tabulates a quantity at noon each day; everything in between
  comes from here. The method is the difference table: the first
  differences of three or five neighboring values, then the differences of
  those, and so on, and the interpolating polynomial is written in them.
  Its advantage over fitting a polynomial directly is that the size of the
  last difference says at a glance whether the table is dense enough --
  if the fourth differences are not small, no amount of arithmetic will
  make the interpolated value good.

  A table here is a map from `table-3` or `table-5`, carrying the first
  and last abscissa and the differences, so the same table can be asked
  for values, an extremum and a zero without recomputing them. The
  interpolating factor n runs from -1 to 1 across a 3-point table (0 at
  the middle value) and from -2 to 2 across a 5-point one."
  (:require [clojure.math :as math]))

(defn horner
  "c0 + c1 x + c2 x^2 + ..., evaluated from the highest power down."
  [x coeffs]
  (reduce (fn [acc c] (+ (* acc x) c)) 0.0 (rseq (vec coeffs))))

(defn- iterate-n
  "Iterate n <- (f n) from zero until it settles. nil if it diverges."
  [f]
  (loop [n0 0.0 i 0]
    (let [n1 (f n0)]
      (cond
        (or (NaN? n1) (infinite? n1) (>= i 50)) nil
        (<= (abs (- n1 n0)) (* 1e-15 (max 1.0 (abs n1)))) n1
        :else (recur n1 (inc i))))))

;; ------------------------------------------------------------ three values

(defn table-3
  "A 3-point table: `ys` at the abscissas `x1`, the midpoint, and `x3`."
  [x1 x3 ys]
  (let [[y1 y2 y3] ys
        a (- y2 y1)
        b (- y3 y2)]
    {:x1 x1 :x3 x3 :y2 y2 :a a :b b :c (- b a)}))

(defn value-n
  "Interpolated value at the factor `n` (Meeus 3.3). For a 3-point table n
  is -1 at the first value and 1 at the last; for a 5-point one, -2 and 2."
  [{:keys [y2 a b c coeffs]} n]
  (if coeffs
    (horner n coeffs)
    (+ y2 (* 0.5 n (+ a b (* n c))))))

(defn- n-of [{:keys [x1 x3 x5]} x]
  (if x5
    (/ (- (* 4.0 x) (* 2.0 (+ x5 x1))) (- x5 x1))
    (/ (- (* 2.0 x) (+ x3 x1)) (- x3 x1))))

(defn- x-of [{:keys [x1 x3 x5]} n]
  (if x5
    (+ (* 0.5 (+ x5 x1)) (* 0.25 (- x5 x1) n))
    (* 0.5 (+ (+ x3 x1) (* (- x3 x1) n)))))

(defn value
  "Interpolated value at abscissa `x`."
  [t x]
  (value-n t (n-of t x)))

(defn- extremum-3 [{:keys [y2 a b c] :as t}]
  (when-not (zero? c)
    (let [n (/ (+ a b) (* -2.0 c))]              ; (3.5)
      (when (<= -1.0 n 1.0)
        [(x-of t n) (- y2 (/ (* (+ a b) (+ a b)) (* 8.0 c)))]))))   ; (3.4)

(defn- zero-3 [{:keys [y2 a b c] :as t} strong?]
  (let [f (if strong?
            ;; Newton's method on the interpolating quadratic (3.7), for
            ;; when the function is too curved for the simpler iteration
            (fn [n] (- n (/ (+ (* 2.0 y2) (* n (+ a b (* c n))))
                            (+ a b (* 2.0 c n)))))
            (fn [n] (/ (* -2.0 y2) (+ a b (* c n)))))    ; (3.6)
        n (iterate-n f)]
    (when (and n (<= -1.0 n 1.0)) (x-of t n))))

;; ------------------------------------------------------------- five values

(defn table-5
  "A 5-point table: `ys` equally spaced from `x1` to `x5`.

  Five points carry the fourth differences, which is what a quantity with
  real curvature over the tabular interval -- the Moon's position, most
  obviously -- needs to interpolate to the precision it is tabulated to."
  [x1 x5 ys]
  (let [[y1 y2 y3 y4 y5] ys
        a (- y2 y1) b (- y3 y2) c (- y4 y3) d (- y5 y4)
        e (- b a) f (- c b) g (- d c)
        h (- f e) j (- g f)
        k (- j h)]
    {:x1 x1 :x5 x5 :y3 y3 :b b :c c :f f :h h :j j :k k
     ;; (3.8) as a polynomial in n
     :coeffs [y3
              (- (/ (+ b c) 2.0) (/ (+ h j) 12.0))
              (- (/ f 2.0) (/ k 24.0))
              (/ (+ h j) 12.0)
              (/ k 24.0)]}))

(defn- extremum-5 [{:keys [b c f h j k coeffs] :as t}]
  (let [den (- k (* 12.0 f))]
    (when-not (zero? den)
      ;; (3.9)
      (when-let [n (iterate-n (fn [n] (/ (horner n [(- (* 6.0 (+ b c)) h j)
                                                    0.0
                                                    (* 3.0 (+ h j))
                                                    (* 2.0 k)])
                                         den)))]
        (when (<= -2.0 n 2.0)
          [(x-of t n) (horner n coeffs)])))))

(defn- zero-5 [{:keys [y3 b c f h j k coeffs] :as t} strong?]
  (let [iter (if strong?
               (let [[_ q p nn m] coeffs
                     d [q (* 2.0 p) (* 3.0 nn) (* 4.0 m)]]
                 (fn [n] (- n (/ (horner n coeffs) (horner n d)))))
               ;; (3.10)
               (let [den (- (* 12.0 (+ b c)) (* 2.0 (+ h j)))
                     num [(* -24.0 y3) 0.0 (- k (* 12.0 f)) (* -2.0 (+ h j)) (- k)]]
                 (fn [n] (/ (horner n num) den))))
        n (iterate-n iter)]
    (when (and n (<= -2.0 n 2.0)) (x-of t n))))

;; ------------------------------------------------------------------ either

(defn extremum
  "`[x y]` of the maximum or minimum of the interpolating polynomial, or
  nil if it has none within the table."
  [t]
  (if (:x5 t) (extremum-5 t) (extremum-3 t)))

(defn zero
  "The abscissa at which the interpolating polynomial crosses zero, or nil
  if it does not do so within the table.

  `strong?` iterates by Newton's method instead of Meeus's simpler
  rearrangement, which fails to converge when the curvature is large
  compared with the slope -- Meeus's example 3.d."
  ([t] (zero t false))
  ([t strong?]
   (if (:x5 t) (zero-5 t strong?) (zero-3 t strong?))))

(defn table-3-around
  "The 3-point table centered on the entry nearest `x`, from a longer
  table of `ys` equally spaced from `x1` to `xn`."
  [x x1 xn ys]
  (let [ys (vec ys)
        cnt (count ys)]
    (if (<= cnt 3)
      (table-3 x1 xn ys)
      (let [step (/ (- xn x1) (dec cnt))
            i    (-> (long (math/floor (+ 0.5 (/ (- x x1) step))))
                     (max 1) (min (- cnt 2)))]
        (table-3 (+ x1 (* (dec i) step)) (+ x1 (* (inc i) step))
                 (subvec ys (dec i) (+ i 2)))))))

(defn half-way
  "The value midway between the two central entries of four equally spaced
  ones -- the one place the general formula simplifies to a
  single line."
  [[y1 y2 y3 y4]]
  (/ (- (* 9.0 (+ y2 y3)) y1 y4) 16.0))

;; ---------------------------------------------------------------- Lagrange

(defn lagrange
  "Lagrange's interpolating polynomial through `points`, `[[x y] ...]` at
  any spacing, evaluated at `x`."
  [points x]
  (reduce + (for [[xi yi] points]
              (* yi (reduce * (for [[xj] points :when (not= xi xj)]
                                (/ (- x xj) (- xi xj))))))))

(defn lagrange-coefficients
  "The same polynomial as coefficients `[c0 c1 c2 ...]` of ascending
  powers of x, rather than evaluated at one point."
  [points]
  (let [n (count points)
        ;; multiply a polynomial by (x - r)
        times (fn [p r] (mapv - (into [0.0] p) (conj (mapv #(* r %) p) 0.0)))]
    (reduce (fn [acc [xi yi]]
              (let [others (for [[xj] points :when (not= xi xj)] xj)
                    basis  (reduce times [1.0] others)
                    den    (reduce * (map #(- xi %) others))]
                (mapv + acc (map #(/ (* yi %) den) basis))))
            (vec (repeat n 0.0))
            points)))
