(ns allgo.numerics.splines
  "Spline interpolation and bilinear interpolation (Chapra and Canale,
  *Numerical Methods for Engineers*, 18.6-18.7).

  A spline is a low-order polynomial on each interval between the data,
  pieced together smoothly: straight lines meeting at the points; parabolas
  whose slopes match there too, the first a straight line to close the
  count of conditions; or cubics whose slopes and curvatures match -- a
  draftsman's flexible strip -- with the ends free of curvature (natural),
  given slopes (clamped), or the first and last pairs of pieces one cubic
  (not-a-knot). Unlike one high-degree polynomial through many points, a
  spline does not oscillate between them.

  Each constructor returns a function of x, extrapolating beyond the ends
  with the end pieces, carrying its pieces as metadata `:pieces`."
  (:require [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.tridiagonal :as tri]))

(defn- piece-index
  "The interval of the sorted knots `xs` that x falls in, 0 to n - 2."
  [xs x]
  (let [n (count xs)]
    (loop [lo 0 hi (dec n)]
      (if (<= (- hi lo) 1)
        (min lo (- n 2))
        (let [mid (quot (+ lo hi) 2)]
          (if (< x (xs mid)) (recur lo mid) (recur mid hi)))))))

(defn- sorted [xs ys]
  (let [pairs (sort-by first (map vector xs ys))]
    [(mapv (comp double first) pairs) (mapv (comp double second) pairs)]))

(defn linear
  "The first-order spline through `xs` `ys`: straight lines between the
  points."
  [xs ys]
  (let [[xs ys] (sorted xs ys)]
    (fn [x]
      (let [i (piece-index xs x)
            t (/ (- x (xs i)) (- (xs (inc i)) (xs i)))]
        (+ (ys i) (* t (- (ys (inc i)) (ys i))))))))

(defn quadratic
  "The second-order spline through `xs` `ys`: a parabola a (x - x_i)^2 +
  b (x - x_i) + c on each interval, meeting its neighbors with the same
  slope, the first one straight (a = 0) -- the book's closing condition.
  The slopes at the knots then follow one from the next: each piece's end
  slope is twice its chord slope less its start slope."
  [xs ys]
  (let [[xs ys] (sorted xs ys)
        n (count xs)
        h (mapv #(- (xs (inc %)) (xs %)) (range (dec n)))
        chord (mapv #(/ (- (ys (inc %)) (ys %)) (h %)) (range (dec n)))
        ;; the slope at each knot, the first piece a line
        slopes (reduce (fn [s i] (conj s (- (* 2.0 (chord i)) (s i)))) [(chord 0)] (range (dec n)))
        pieces (mapv (fn [i] {:a (/ (- (chord i) (slopes i)) (h i)) :b (slopes i) :c (ys i) :x (xs i)}) (range (dec n)))]
    (with-meta
      (fn [x]
        (let [{:keys [a b c] x0 :x} (pieces (piece-index xs x))
              d (- x x0)]
          (+ (* a d d) (* b d) c)))
      {:pieces pieces})))

(defn cubic
  "The cubic spline through `xs` `ys`: cubics meeting with matching slope
  and curvature, the second derivatives at the knots from a tridiagonal
  system (solved by the Thomas algorithm). `:end` is `:natural` (the
  default, zero curvature at both ends), `:not-a-knot`, or `[s0 sn]` for
  given end slopes (clamped)."
  ([xs ys] (cubic xs ys {}))
  ([xs ys {:keys [end] :or {end :natural}}]
   (let [[xs ys] (sorted xs ys)
         n (count xs)
         h (mapv #(- (xs (inc %)) (xs %)) (range (dec n)))
         slope (fn [i] (/ (- (ys (inc i)) (ys i)) (h i)))
         ;; interior rows: h_{i-1} M_{i-1} + 2(h_{i-1} + h_i) M_i + h_i M_{i+1}
         ;;                = 6 (slope_i - slope_{i-1})
         interior (mapv (fn [i] [(h (dec i)) (* 2.0 (+ (h (dec i)) (h i))) (h i)
                                 (* 6.0 (- (slope i) (slope (dec i))))])
                        (range 1 (dec n)))
         M (cond
             (= end :natural)
             (let [sol (if (seq interior)
                         (tri/solve-vectors (map first interior) (map second interior)
                                            (map #(nth % 2) interior) (map #(nth % 3) interior))
                         [])]
               (vec (concat [0.0] sol [0.0])))
             (vector? end)
             (let [[s0 sn] end
                   first-row [0.0 (* 2.0 (h 0)) (h 0) (* 6.0 (- (slope 0) s0))]
                   last-row [(h (- n 2)) (* 2.0 (h (- n 2))) 0.0 (* 6.0 (- sn (slope (- n 2))))]
                   rows (vec (concat [first-row] interior [last-row]))]
               (tri/solve-vectors (map first rows) (map second rows) (map #(nth % 2) rows) (map #(nth % 3) rows)))
             (= end :not-a-knot)
             ;; the third derivative continuous across the second and the
             ;; second-to-last knots; not tridiagonal, so by Gauss elimination
             (let [row (fn [cols] (let [v (vec (repeat n 0.0))] (reduce (fn [v [j c]] (assoc v j c)) v cols)))
                   A (vec (concat [(row [[0 (h 1)] [1 (- (+ (h 0) (h 1)))] [2 (h 0)]])]
                                  (map-indexed (fn [k [a b c]] (row [[k a] [(inc k) b] [(+ k 2) c]])) interior)
                                  [(row [[(- n 3) (h (- n 2))] [(- n 2) (- (+ (h (- n 3)) (h (- n 2))))] [(dec n) (h (- n 3))]])]))
                   b (vec (concat [0.0] (map #(nth % 3) interior) [0.0]))]
               (ls/gauss A b)))
         pieces (mapv (fn [i] {:x0 (xs i) :x1 (xs (inc i)) :y0 (ys i) :y1 (ys (inc i))
                               :m0 (M i) :m1 (M (inc i)) :h (h i)})
                      (range (dec n)))]
     (with-meta
       (fn [x]
         (let [{:keys [x0 x1 y0 y1 m0 m1 h]} (pieces (piece-index xs x))
               a (- x1 x) b (- x x0)]
           (+ (/ (* m0 a a a) (* 6.0 h)) (/ (* m1 b b b) (* 6.0 h))
              (* (- (/ y0 h) (/ (* m0 h) 6.0)) a)
              (* (- (/ y1 h) (/ (* m1 h) 6.0)) b))))
       {:pieces pieces :second-derivatives M}))))

(defn bilinear
  "The value at `x` `y` interpolated bilinearly from the four corners of
  a rectangle, `f11` at (x1, y1), `f21` at (x2, y1), `f12` at (x1, y2) and
  `f22` at (x2, y2): linear in each direction in turn."
  [x1 x2 y1 y2 f11 f21 f12 f22 x y]
  (let [tx (/ (- x x1) (- x2 x1)) ty (/ (- y y1) (- y2 y1))]
    (+ (* (- 1.0 tx) (- 1.0 ty) f11) (* tx (- 1.0 ty) f21)
       (* (- 1.0 tx) ty f12) (* tx ty f22))))
