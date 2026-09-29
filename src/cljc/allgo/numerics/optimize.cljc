(ns allgo.numerics.optimize
  "Minimizing a function of several variables without its derivatives:
  the simplex method of Nelder and Mead (The Computer Journal 7, 1965),
  as Lagarias, Reeds, Wright and Wright state it precisely (SIAM J.
  Optim. 9, 1998) -- reflection 1, expansion 2, contraction 1/2,
  shrinkage 1/2.

  A simplex of n + 1 points walks downhill, each step replacing its worst
  point by one reflected through the others' centroid, stretched further
  if that pays, pulled in if it doesn't, and the whole simplex shrinking
  toward its best point when nothing else works. It needs only values, so
  it copes with functions that are merely continuous -- a launch energy
  over departure and arrival dates, say -- at the price of speed.")

(defn- add [a b] (mapv + a b))
(defn- sub [a b] (mapv - a b))
(defn- scale [a s] (mapv #(* s %) a))

(defn- replace-worst
  "The simplex with its worst point swapped for `p`, sorted best first."
  [simplex p]
  (vec (sort-by second (conj (pop simplex) p))))

(defn- centroid [points]
  (scale (reduce add points) (/ 1.0 (count points))))

(defn nelder-mead
  "The `x` near `x0` where `f` is least: `{:x :f :iterations}`. `f`
  takes a vector; it may return nil or ##Inf where it is undefined, which
  the simplex treats as uphill. Options:

    :step      the initial simplex's edge along each axis, a number or a
               vector (default 0.05 of each coordinate, or 0.00025 where
               the coordinate is zero, as in Matlab's fminsearch)
    :tol-x     stop when every point is this close to the best (1e-8)
    :tol-f     ... and their values this close to its value (1e-10)
    :max-iter  or after this many steps (200 per variable)"
  ([f x0] (nelder-mead f x0 {}))
  ([f x0 {:keys [step tol-x tol-f max-iter]
          :or {tol-x 1e-8 tol-f 1e-10}}]
   (let [n (count x0)
         x0 (mapv double x0)
         max-iter (or max-iter (* 200 n))
         g (fn [x] (let [v (f x)] (if (and v (not (NaN? v))) (double v) ##Inf)))
         steps (cond (vector? step) step
                     (number? step) (vec (repeat n step))
                     :else (mapv #(if (zero? %) 0.00025 (* 0.05 %)) x0))
         start (cons x0 (for [i (range n)] (update x0 i + (steps i))))
         sorted (fn [ps] (vec (sort-by second ps)))]
     (loop [simplex (sorted (map (juxt identity g) start)) k 0]
       (let [[[xb fb]] simplex
             [xw fw] (peek simplex)
             [_ fs] (simplex (dec n))]
         (if (or (>= k max-iter)
                 (and (every? (fn [[x _]] (every? #(<= (abs %) tol-x) (sub x xb))) simplex)
                      (every? (fn [[_ v]] (<= (abs (- v fb)) tol-f)) simplex)))
           {:x xb :f fb :iterations k}
           (let [c (centroid (map first (pop simplex)))
                 at (fn [t] (let [x (add c (scale (sub c xw) t))] [x (g x)]))
                 replace #(replace-worst simplex %)
                 [_ fr :as r] (at 1.0)]
             (recur
              (cond
                (< fr fb) (let [e (at 2.0)] (replace (if (< (second e) fr) e r)))
                (< fr fs) (replace r)
                :else
                ;; contract: outside the simplex if the reflection beat
                ;; the worst point, inside if not
                (let [outside? (< fr fw)
                      [_ fc :as cp] (at (if outside? 0.5 -0.5))]
                  (if (if outside? (<= fc fr) (< fc fw))
                    (replace cp)
                    ;; shrink everything toward the best point
                    (sorted (cons [xb fb] (for [[x _] (rest simplex)]
                                            (let [y (add xb (scale (sub x xb) 0.5))] [y (g y)])))))))
              (inc k)))))))))
