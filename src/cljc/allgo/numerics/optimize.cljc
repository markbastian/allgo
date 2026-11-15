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
  over departure and arrival dates, say -- at the price of speed.

  And the rest of the classical optimizers (Chapra and Canale, *Numerical
  Methods for Engineers*, chapters 13-15): in one variable, parabolic
  golden-section search, parabolic interpolation, Newton's method and
  Brent's method; in several, random search, univariate
  and Powell's pattern searches, steepest descent, Fletcher-Reeves
  conjugate gradients, Newton's method, Marquardt's method and the
  quasi-Newton BFGS and DFP updates; the simplex method for linear
  programs; and penalty functions for nonlinear constraints.

  Everything here minimizes; to maximize f, minimize -f.

  Each iterative method is a step, as in `allgo.search` and
  `allgo.numerics.roots`: a private step function from one state to the
  next, a public `-seq` of the states it makes, and the plain function,
  which runs the same steps by `allgo.search/converge` or `settle` to its
  stopping rule. A step returns nil only when the method cannot go on --
  at a point where the gradient vanishes, at a singular Hessian -- never
  because it has converged, so the seqs show the whole path: steepest
  descent zigzagging, conjugate gradients cutting across, a simplex
  reflecting, expanding and shrinking. The one-variable states are
  records; the many-variable states are maps with `:x`, `:f` and
  `:iterations`, and `:last-x` and `:last-f` for the point before."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.roots :as roots]
            [allgo.random :as random]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- replace-worst
  "The simplex with its worst point swapped for `p`, sorted best first."
  [simplex p]
  (vec (sort-by second (conj (pop simplex) p))))

(defn- centroid [points]
  (lin/scale (reduce lin/add points) (/ 1.0 (count points))))

(defn- sorted-simplex [ps] (vec (sort-by second ps)))

(defn- nelder-mead-step
  "One step of the simplex `[[x f] ...]`, best first, `g` the function: its
  worst point reflected through the others' centroid, the reflection
  stretched if it beats the best, pulled back in if it beats nothing --
  outside the simplex if it at least beat the worst, inside if not -- and
  the whole simplex shrunk toward its best point if even that fails.
  `:by` says which."
  [g n {:keys [simplex iterations]}]
  (let [[[xb fb]] simplex
        [xw fw] (peek simplex)
        [_ fs] (simplex (dec n))
        c (centroid (map first (pop simplex)))
        at (fn [t] (let [x (lin/add c (lin/scale (lin/sub c xw) t))] [x (g x)]))
        [_ fr :as r] (at 1.0)
        [simplex by]
        (cond
          (< fr fb) (let [e (at 2.0)] (if (< (second e) fr) [(replace-worst simplex e) :expand] [(replace-worst simplex r) :reflect]))
          (< fr fs) [(replace-worst simplex r) :reflect]
          :else
          (let [outside? (< fr fw)
                [_ fc :as cp] (at (if outside? 0.5 -0.5))]
            (if (if outside? (<= fc fr) (< fc fw))
              [(replace-worst simplex cp) (if outside? :contract-outside :contract-inside)]
              [(sorted-simplex (cons [xb fb] (for [[x _] (rest simplex)]
                                               (let [y (lin/add xb (lin/scale (lin/sub x xb) 0.5))] [y (g y)]))))
               :shrink])))]
    {:simplex simplex :iterations (inc iterations) :by by}))

(defn- nelder-mead-setup
  "The function Nelder-Mead minimizes -- `f` with nil and NaN made
  uphill -- its first state, and its options filled in."
  [f x0 {:keys [step tol-x tol-f max-iter] :or {tol-x 1e-8 tol-f 1e-10}}]
  (let [n (count x0)
        x0 (mapv double x0)
        g (fn [x] (let [v (f x)] (if (and v (not (NaN? v))) (double v) ##Inf)))
        steps (cond (vector? step) step
                    (number? step) (vec (repeat n step))
                    :else (mapv #(if (zero? %) 0.00025 (* 0.05 %)) x0))
        start (cons x0 (for [i (range n)] (update x0 i + (steps i))))]
    {:g g :n n :tol-x tol-x :tol-f tol-f :max-iter (or max-iter (* 200 n))
     :init {:simplex (sorted-simplex (map (juxt identity g) start)) :iterations 0 :by :start}}))

(defn nelder-mead-seq
  "The simplexes of Nelder-Mead from `x0`: `{:simplex :iterations :by}`,
  the simplex `[[x f] ...]` best first, `:by` how it was reached
  (`:reflect`, `:expand`, `:contract-outside`, `:contract-inside` or
  `:shrink`). Options as `nelder-mead`'s."
  ([f x0] (nelder-mead-seq f x0 {}))
  ([f x0 opts]
   (let [{:keys [g n init]} (nelder-mead-setup f x0 opts)]
     (iterate (partial nelder-mead-step g n) init))))

(defn- simplex-small?
  "Whether every point of the simplex is within `tol-x` of the best, and
  its value within `tol-f` of the best's."
  [tol-x tol-f {:keys [simplex]}]
  (let [[[xb fb]] simplex]
    (and (every? (fn [[x _]] (every? #(<= (abs %) tol-x) (lin/sub x xb))) simplex)
         (every? (fn [[_ v]] (<= (abs (- v fb)) tol-f)) simplex))))

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
  ([f x0 opts]
   (let [{:keys [g n tol-x tol-f max-iter init]} (nelder-mead-setup f x0 opts)
         {:keys [simplex iterations]} (search/settle (partial simplex-small? tol-x tol-f) (inc max-iter)
                                                     (partial nelder-mead-step g n) init)
         [[xb fb]] simplex]
     {:x xb :f fb :iterations iterations})))

;; ---------------------------------------------------------- one variable

(def ^:private golden (/ (- (math/sqrt 5.0) 1.0) 2.0))

;; A bracket `[a b]` about an extremum.
(defrecord Bracket [a b])

(defn- golden-minimize-step
  "Narrow the bracket by the golden section, keeping the side of the
  lower of its interior points."
  [f ^Bracket br]
  (let [a (.-a br) b (.-b br)
        x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
    (if (< (f x1) (f x2)) (->Bracket a x2) (->Bracket x1 b))))

(defn- golden-maximize-step
  "Narrow the bracket by the golden section, keeping the side of the
  higher of its interior points."
  [f ^Bracket br]
  (let [a (.-a br) b (.-b br)
        x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
    (if (< (f x1) (f x2)) (->Bracket x1 b) (->Bracket a x2))))

(defn golden-minimize-seq
  "The brackets `{:a :b}` of golden-section search for the minimum of `f`
  between `a` and `b`."
  [f a b]
  (iterate (partial golden-minimize-step f) (->Bracket a b)))

(defn golden-maximize-seq
  "The brackets `{:a :b}` of golden-section search for the maximum of `f`
  between `a` and `b`: those of Chapra and Canale's example 13.1."
  [f a b]
  (iterate (partial golden-maximize-step f) (->Bracket a b)))

(defn- golden-search [step a b tol max-iter]
  (let [^Bracket br (search/settle (fn [^Bracket br] (< (abs (- (.-b br) (.-a br))) tol)) (+ 2 max-iter)
                                   step (->Bracket a b))]
    (* 0.5 (+ (.-a br) (.-b br)))))

(defn golden-minimize
  "The `x` between `a` and `b` where `f`, which has a single minimum there,
  is least, by golden-section search: until the bracket is narrower than
  `:tol`, or after `:max-iter` narrowings."
  ([f a b] (golden-minimize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (golden-search (partial golden-minimize-step f) a b tol max-iter)))

(defn golden-maximize
  "The `x` between `a` and `b` where `f`, which has a single maximum there,
  is greatest, by golden-section search; options as `golden-minimize`'s."
  ([f a b] (golden-maximize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (golden-search (partial golden-maximize-step f) a b tol max-iter)))

;; Three points bracketing a minimum, `x1` the lowest, and the vertex `x`
;; of the parabola through the three before, the step `dx` to it from
;; the vertex before that.
(defrecord Parabola [x0 x1 x2 f0 f1 f2 x dx])

(defn- parabolic-step
  "The vertex of the parabola through the three points, and the three
  that bracket the lowest of the four. nil where the three are in a
  line."
  [f ^Parabola p]
  (let [x0 (.-x0 p) x1 (.-x1 p) x2 (.-x2 p) f0 (.-f0 p) f1 (.-f1 p) f2 (.-f2 p) x (.-x p)
        num (+ (* f0 (- (* x1 x1) (* x2 x2))) (* f1 (- (* x2 x2) (* x0 x0))) (* f2 (- (* x0 x0) (* x1 x1))))
        den (* 2.0 (+ (* f0 (- x1 x2)) (* f1 (- x2 x0)) (* f2 (- x0 x1))))]
    (when-not (zero? den)
      (let [x3 (/ num den) f3 (f x3) dx (if x (- x3 x) ##Inf)]
        (if (> x3 x1)
          (if (< f3 f1) (->Parabola x1 x3 x2 f1 f3 f2 x3 dx) (->Parabola x0 x1 x3 f0 f1 f3 x3 dx))
          (if (< f3 f1) (->Parabola x0 x3 x1 f0 f3 f1 x3 dx) (->Parabola x3 x1 x2 f3 f1 f2 x3 dx)))))))

(defn parabolic-seq
  "The states of successive parabolic interpolation from `x0` < `x1` <
  `x2`: `{:x0 :x1 :x2 :x}`, `:x` the vertex that chose the three -- those
  of Chapra and Canale's example 13.2."
  [f x0 x1 x2]
  (search/step-seq (partial parabolic-step f) (->Parabola x0 x1 x2 (f x0) (f x1) (f x2) nil ##Inf)))

(defn parabolic
  "The minimum of `f` by successive parabolic interpolation from three
  points bracketing it, `x0` < `x1` < `x2` with f(x1) below both ends: the
  vertex of the parabola through the three is the next point, and the
  worst end on its side is dropped. Until the vertex moves by less than
  `:tol` (default 1e-10); `x` at the minimum."
  ([f x0 x1 x2] (parabolic f x0 x1 x2 {}))
  ([f x0 x1 x2 {:keys [tol max-iter] :or {tol 1e-10 max-iter 200}}]
   (let [done? (fn [^Parabola p] (and (.-x p) (<= (abs (.-dx p)) (* tol (max 1.0 (abs (.-x p)))))))
         p (search/settle done? (inc max-iter) (partial parabolic-step f) (->Parabola x0 x1 x2 (f x0) (f x1) (f x2) nil ##Inf))]
     (if (done? p) (:x p) (:x1 p)))))

(defn newton-1d-seq
  "The estimates `{:x :dx}` of Newton's method on the derivative `df`
  from `x0`: those of Chapra and Canale's example 13.3."
  [df d2f x0]
  (roots/newton-seq df d2f x0))

(defn newton-1d
  "A stationary point of f -- a minimum where f'' > 0 -- by Newton's
  method on its derivative `df`, `d2f` the second derivative: x <- x -
  f'(x)/f''(x). From `x0`; nil if it does not settle in `:max-iter`
  (default 50)."
  ([df d2f x0] (newton-1d df d2f x0 {}))
  ([df d2f x0 opts] (roots/newton df d2f x0 opts)))

(def ^:private cgold (* 0.5 (- 3.0 (math/sqrt 5.0))))

;; Brent's minimizer: the interval `[a b]`, `x` the best point, `w` the
;; second best and `v` the previous w, with their values; `d` the last
;; step and `e` the one before; `by` how the last was taken.
(defrecord BrentMinimum [a b x w v fx fw fv d e by])

(defn- brent-tolerance [tol x] (+ (* tol (abs x)) 1e-12))

(defn- brent-minimum-done?
  "Whether Brent's interval is within its tolerance of the best point."
  [tol ^BrentMinimum s]
  (let [a (.-a s) b (.-b s) x (.-x s)]
    (<= (abs (- x (* 0.5 (+ a b)))) (- (* 2.0 (brent-tolerance tol x)) (* 0.5 (- b a))))))

(defn- brent-minimum-next
  "Brent's state after the step `d` from x, `e` the step before and `by`
  how it was chosen: u = x + d, at least the tolerance from x, and the
  interval and the three best points updated with it."
  [f ^BrentMinimum s t1 d e by]
  (let [a (.-a s) b (.-b s) x (.-x s) w (.-w s) v (.-v s) fx (.-fx s) fw (.-fw s) fv (.-fv s)
        u (+ x (if (>= (abs d) t1) d (if (pos? d) t1 (- t1))))
        fu (f u)]
    (if (<= fu fx)
      (if (< u x) (->BrentMinimum a x u x w fu fx fw d e by) (->BrentMinimum x b u x w fu fx fw d e by))
      (let [a (if (< u x) u a) b (if (< u x) b u)]
        (cond
          (or (<= fu fw) (= w x)) (->BrentMinimum a b x u w fx fu fw d e by)
          (or (<= fu fv) (= v x) (= v w)) (->BrentMinimum a b x w u fx fw fu d e by)
          :else (->BrentMinimum a b x w v fx fw fv d e by))))))

(defn- brent-minimize-step
  "One step of Brent's minimizer: the vertex of the parabola through x,
  w and v where it is safe -- inside the interval and shorter than half
  the step before last -- a golden section of the larger side
  otherwise; `:by` says which."
  [f tol ^BrentMinimum s]
  (let [a (.-a s) b (.-b s) x (.-x s) w (.-w s) v (.-v s) fx (.-fx s) fw (.-fw s) fv (.-fv s) d (.-d s) e (.-e s)
        m (* 0.5 (+ a b))
        t1 (brent-tolerance tol x)
        t2 (* 2.0 t1)
        ge (if (< x m) (- b x) (- a x))]
    (if (> (abs e) t1)
      ;; try the parabola through x, w and v
      (let [r (* (- x w) (- fx fv))
            q (* (- x v) (- fx fw))
            p (- (* (- x v) q) (* (- x w) r))
            q (* 2.0 (- q r))
            p (if (pos? q) (- p) p)
            q (if (pos? q) q (- q))]
        (if (and (< (abs p) (abs (* 0.5 q e))) (> p (* q (- a x))) (< p (* q (- b x))))
          ;; the parabola's step taken; the one before it remembered
          (let [step (/ p q) u (+ x step)]
            (brent-minimum-next f s t1 (if (or (< (- u a) t2) (< (- b u) t2)) (if (< x m) t1 (- t1)) step) d :parabolic))
          (brent-minimum-next f s t1 (* cgold ge) ge :golden)))
      (brent-minimum-next f s t1 (* cgold ge) ge :golden))))

(defn- brent-minimum-start [f a b]
  (let [x0 (+ a (* cgold (- b a))) f0 (f x0)]
    (->BrentMinimum a b x0 x0 x0 f0 f0 f0 0.0 0.0 :start)))

(defn brent-minimize-seq
  "The states of Brent's minimizer between `a` and `b`: `{:a :b :x :fx
  :by}`, `:x` the best point so far, `:by` `:parabolic` or `:golden`.
  `tol` as `brent-minimize`'s."
  [f a b tol]
  (iterate (partial brent-minimize-step f tol) (brent-minimum-start f a b)))

(defn brent-minimize
  "The minimum of `f` between `a` and `b` by Brent's method (Brent,
  Algorithms for Minimization without Derivatives, 1973): a parabolic
  step through the best three points where it is safe -- inside the
  interval and shorter than half the step before last -- golden section
  otherwise. Superlinear on smooth functions, never slower than golden
  section. To within `:tol` relative (default 1e-10); `{:x :f}`. No
  minimizer places x better than about the square root of the machine
  precision, 1e-8: near a minimum f changes only as the step squared."
  ([f a b] (brent-minimize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 500}}]
   (let [^BrentMinimum s (search/settle (partial brent-minimum-done? tol) (inc max-iter)
                                        (partial brent-minimize-step f tol) (brent-minimum-start f a b))]
     {:x (.-x s) :f (.-fx s)})))

;; ---------------------------------------------------------- line search

;; Walking out from 0 by doubling steps: the last three points `a`, `b`
;; and `c`, and the values at the last two.
(defrecord Walk [a b c fb fc])

(defn- walk-step [g ^Walk w]
  (let [c' (* 2.0 (.-c w))] (->Walk (.-b w) (.-c w) c' (.-fc w) (g c'))))

(defn- bracket-forward
  "Walking downhill from 0 by doubling steps from `h`, `g` falling there:
  the first three points whose last rises, `[a c]` its ends."
  [g h]
  (let [^Walk w (search/settle (fn [^Walk w] (or (>= (.-fc w) (.-fb w)) (> (.-c w) 1e12))) ##Inf
                               (partial walk-step g) (->Walk 0.0 h (* 2.0 h) (g h) (g (* 2.0 h))))]
    [(.-a w) (.-c w)]))

;; Walking back from 0 by doubling steps: how many, the point `a` reached
;; and the values there and at the point before.
(defrecord Retreat [k a fa fb])

(defn- retreat-step [g ^Retreat r]
  (let [a' (* 2.0 (.-a r))] (->Retreat (inc (.-k r)) a' (g a') (.-fa r))))

(defn- bracket-backward
  "Walking from 0 the other way, `g` rising forward: doubling steps back
  from -`h` until g stops falling, or 61 doublings. `[a h]`."
  [g h f0]
  (let [^Retreat r (search/settle (fn [^Retreat r] (or (>= (.-fa r) (.-fb r)) (> (.-k r) 60))) ##Inf
                                  (partial retreat-step g) (->Retreat 0 (- h) (g (- h)) f0))]
    [(.-a r) h]))

(defn line-minimize
  "The step t along direction `d` from `x` that minimizes f(x + t d):
  bracketed by stepping out from 0 with growing steps, then found by
  Brent's method. `[t x+td f]`."
  [f x d]
  (let [g (fn [t] (f (lin/add x (lin/scale d t))))
        f0 (g 0.0)
        h (/ 1e-2 (max 1e-12 (math/sqrt (reduce + (map * d d)))))
        [lo hi] (if (< (g h) f0) (bracket-forward g h) (bracket-backward g h f0))
        {t :x fv :f} (brent-minimize g lo hi {:tol 1e-12})]
    [t (lin/add x (lin/scale d t)) fv]))

;; ---------------------------------------------- several variables: direct

(defn random-search
  "The best of `n` points drawn uniformly in the box between the corners
  `lo` and `hi`: `{:x :f}`. Crude and sure -- it finds the global minimum
  of anything, given enough points -- and a way to start the others.
  `:rng` a function of no arguments giving numbers in [0, 1) (default
  `(allgo.random/rng 1)`)."
  ([f lo hi n] (random-search f lo hi n {}))
  ([f lo hi n {:keys [rng]}]
   (let [rng (or rng (random/rng 1))]
     (reduce (fn [best _]
               (let [x (mapv (fn [l h] (+ l (* (- h l) (rng)))) lo hi)
                     fx (f x)]
                 (if (< fx (:f best)) {:x x :f fx} best)))
             {:x nil :f ##Inf}
             (range n)))))

(defn- settled? [x x' fx fx' tol]
  (and (<= (abs (- fx fx')) (* tol (max 1.0 (abs fx'))))
       (every? true? (map #(<= (abs (- %1 %2)) (* (math/sqrt tol) (max 1.0 (abs %2)))) x x'))))

(defn- moved-little?
  "Whether a many-variable state has moved, and changed its value, by
  less than `tol` from the point before -- or has stopped."
  [tol {:keys [x f last-x last-f stopped?]}]
  (or stopped? (and last-x (settled? last-x x last-f f tol))))

(defn- start [f x0] {:x (mapv double x0) :f (f x0) :iterations 0})

(defn- moved
  "The state `st` moved to `x'` with value `fx'`, and the map `more`
  merged in."
  ([st x' fx'] (moved st x' fx' nil))
  ([{:keys [x iterations] fx :f :as st} x' fx' more]
   (into (assoc st :x x' :f fx' :iterations (inc iterations) :last-x x :last-f fx) more)))

(defn- descend
  "The state the many-variable method's `step` reaches from `init`: the
  first that has moved little, or the last within `max-iter` steps --
  `{:x :f :iterations}`."
  [tol max-iter step init]
  (select-keys (search/settle (partial moved-little? tol) (inc max-iter) step init) [:x :f :iterations]))

(defn- axis [n i] (assoc (vec (repeat n 0.0)) i 1.0))

(defn- univariate-step
  "A line minimization along each coordinate in turn."
  [f {:keys [x] fx :f :as st}]
  (let [[x' fx'] (reduce (fn [[x _] i] (let [[_ x' fv] (line-minimize f x (axis (count x) i))] [x' fv]))
                         [x fx] (range (count x)))]
    (moved st x' fx')))

(defn univariate-seq
  "The points `{:x :f :iterations}` of univariate search from `x0`, one
  round of the axes a step."
  [f x0]
  (iterate (partial univariate-step f) (start f x0)))

(defn univariate
  "The minimum of `f` from `x0` by univariate search: a line minimization
  along each coordinate in turn, round and round -- simple, and slow in a
  narrow valley that runs across the axes. `{:x :f :iterations}`."
  ([f x0] (univariate f x0 {}))
  ([f x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 1000}}]
   (descend tol max-iter (partial univariate-step f) (start f x0))))

(defn- powell-step
  "A line minimization along each direction in turn, then along the net
  displacement, which replaces the oldest direction."
  [f {:keys [x dirs] fx :f :as st}]
  (let [[y fy] (reduce (fn [[y _] d] (let [[_ y' fv] (line-minimize f y d)] [y' fv])) [x fx] dirs)
        net (lin/sub y x)]
    (if (every? zero? net)
      (moved st y fy)
      (let [[_ z fz] (line-minimize f y net)]
        (moved st z fz {:dirs (conj (subvec dirs 1) net)})))))

(defn- powell-start [f x0]
  (let [n (count x0)] (assoc (start f x0) :dirs (mapv (partial axis n) (range n)))))

(defn powell-seq
  "The points `{:x :f :iterations :dirs}` of Powell's method from `x0`,
  one round of the directions a step."
  [f x0]
  (iterate (partial powell-step f) (powell-start f x0)))

(defn powell
  "The minimum of `f` from `x0` by Powell's method (Powell 1964): line
  minimizations along a set of directions, starting with the axes; after
  each round the net displacement becomes a new direction and the oldest
  is dropped. The directions it builds are conjugate, and a quadratic in n
  variables is minimized in n rounds. `{:x :f :iterations}`."
  ([f x0] (powell f x0 {}))
  ([f x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 500}}]
   (descend tol max-iter (partial powell-step f) (powell-start f x0))))

;; -------------------------------------------- several variables: gradients

(defn- grad-fn [f grad] (or grad #(d/gradient f %)))

(defn- steepest-descent-step
  "A line minimization down the gradient. nil where it vanishes."
  [f grad {:keys [x] :as st}]
  (let [g (grad x)]
    (when-not (every? zero? g)
      (let [[_ x' fx'] (line-minimize f x (lin/scale g -1.0))]
        (moved st x' fx')))))

(defn steepest-descent-seq
  "The points `{:x :f :iterations}` of steepest descent from `x0`:
  those of Chapra and Canale's example 14.4. `grad` as
  `steepest-descent`'s."
  [f grad x0]
  (search/step-seq (partial steepest-descent-step f (grad-fn f grad)) (start f x0)))

(defn steepest-descent
  "The minimum of `f` from `x0` by steepest descent: a line minimization
  along the negative gradient, again and again -- `grad` its gradient
  function, or central differences if nil. Each step turns a right angle
  from the last, so it zigzags down narrow valleys. `{:x :f
  :iterations}`."
  ([f x0] (steepest-descent f nil x0 {}))
  ([f grad x0] (steepest-descent f grad x0 {}))
  ([f grad x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 5000}}]
   (descend tol max-iter (partial steepest-descent-step f (grad-fn f grad)) (start f x0))))

(defn- conjugate-gradient-step
  "A line minimization along the direction `d`, and the next direction:
  the new negative gradient plus d times |g_new|^2 / |g_old|^2, or the
  negative gradient alone every `n` steps. nil where the gradient
  vanishes."
  [f grad n {:keys [x g d iterations] :as st}]
  (when-not (every? zero? g)
    (let [[_ x' fx'] (line-minimize f x d)
          g' (grad x')
          beta (/ (reduce + (map * g' g')) (reduce + (map * g g)))
          k (inc iterations)
          d' (if (zero? (mod k n)) (lin/scale g' -1.0) (lin/add (lin/scale g' -1.0) (lin/scale d beta)))]
      (moved st x' fx' {:g g' :d d'}))))

(defn- conjugate-gradient-start [f grad x0]
  (assoc (start f x0) :g (grad x0) :d (lin/scale (grad x0) -1.0)))

(defn conjugate-gradient-seq
  "The points `{:x :f :iterations :g :d}` of Fletcher-Reeves conjugate
  gradients from `x0`, `:g` the gradient there and `:d` the next
  direction. Arguments as `conjugate-gradient`'s."
  [f grad x0]
  (let [grad (grad-fn f grad)]
    (search/step-seq (partial conjugate-gradient-step f grad (count x0)) (conjugate-gradient-start f grad x0))))

(defn conjugate-gradient
  "The minimum of `f` from `x0` by the Fletcher-Reeves conjugate gradient
  method: each search direction the negative gradient plus the last
  direction times |g_new|^2 / |g_old|^2, restarted down the gradient
  every n steps -- a quadratic minimized in n line searches, at the cost
  of steepest descent's storage. Arguments as `steepest-descent`'s."
  ([f x0] (conjugate-gradient f nil x0 {}))
  ([f grad x0] (conjugate-gradient f grad x0 {}))
  ([f grad x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 5000}}]
   (let [grad (grad-fn f grad)]
     (descend tol max-iter (partial conjugate-gradient-step f grad (count x0)) (conjugate-gradient-start f grad x0)))))

(defn- solve-sym [H g] (ls/gauss H g))

(defn- hess-fn [f hess] (or hess #(d/hessian f %)))

(defn- newton-step
  "Newton's step -H^-1 g. nil at a singular Hessian."
  [f grad hess {:keys [x] :as st}]
  (when-let [dx (solve-sym (hess x) (grad x))]
    (let [x' (lin/sub x dx)]
      (moved st x' (f x')))))

(defn newton-seq
  "The points `{:x :f :iterations}` of Newton's method from `x0`.
  Arguments as `newton`'s."
  [f grad hess x0]
  (search/step-seq (partial newton-step f (grad-fn f grad) (hess-fn f hess)) (start f x0)))

(defn newton
  "The minimum of `f` from `x0` by Newton's method: each step -H^-1 g, H
  the Hessian and g the gradient (`grad` and `hess` functions of x, or
  central differences if nil) -- quadratic convergence near the minimum,
  none promised away from it, where H may not even be positive definite.
  `{:x :f :iterations}`."
  ([f x0] (newton f nil nil x0 {}))
  ([f grad hess x0] (newton f grad hess x0 {}))
  ([f grad hess x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (descend tol max-iter (partial newton-step f (grad-fn f grad) (hess-fn f hess)) (start f x0))))

(defn- marquardt-step
  "One trial of Marquardt's step -(H + alpha I)^-1 g: taken, and alpha
  shrunk tenfold, if it lowers f; refused, and alpha grown tenfold, if
  not. `:stopped?` once a step moves too little to matter -- the step
  is what says so, since a refused step leaves the point where it was --
  or at `max-iter`; nil after that."
  [f grad hess tol max-iter {:keys [x alpha iterations stopped?] fx :f :as st}]
  (when-not stopped?
    (let [n (count x) k (inc iterations)
          H (hess x)
          dx (solve-sym (mapv (fn [i row] (update row i + alpha)) (range n) H) (grad x))]
      (if (nil? dx)
        (assoc st :alpha (* 10.0 alpha) :iterations k)
        (let [x' (lin/sub x dx) fx' (f x')]
          (cond
            (>= k max-iter) (assoc st :iterations k :stopped? true)
            (< fx' fx) (assoc st :x x' :f fx' :alpha (* 0.1 alpha) :iterations k :stopped? (settled? x x' fx fx' tol))
            (settled? x x' fx fx' tol) (assoc st :iterations k :stopped? true)
            :else (assoc st :alpha (* 10.0 alpha) :iterations k)))))))

(defn marquardt-seq
  "The trials `{:x :f :alpha :iterations}` of Marquardt's method from
  `x0`. Arguments as `marquardt`'s."
  ([f grad hess x0] (marquardt-seq f grad hess x0 {}))
  ([f grad hess x0 {:keys [tol max-iter alpha] :or {tol 1e-14 max-iter 500 alpha 1e3}}]
   (search/step-seq (partial marquardt-step f (grad-fn f grad) (hess-fn f hess) tol max-iter)
                    (assoc (start f x0) :alpha alpha))))

(defn marquardt
  "The minimum of `f` from `x0` by Marquardt's method: each step -(H +
  alpha I)^-1 g, alpha shrunk after a step that lowers f and grown after
  one that doesn't -- steepest descent's reliability far from the
  minimum, Newton's speed near it. `:alpha` starts it (default 1e3).
  Arguments as `newton`'s."
  ([f x0] (marquardt f nil nil x0 {}))
  ([f grad hess x0] (marquardt f grad hess x0 {}))
  ([f grad hess x0 {:keys [tol max-iter alpha] :or {tol 1e-14 max-iter 500 alpha 1e3}}]
   (-> (search/settle :stopped? (inc max-iter)
                      (partial marquardt-step f (grad-fn f grad) (hess-fn f hess) tol max-iter)
                      (assoc (start f x0) :alpha alpha))
       (select-keys [:x :f :iterations]))))

(defn- outer [u v] (mapv (fn [a] (mapv #(* a %) v)) u))
(defn- mat-sum [& ms] (apply mapv (fn [& rows] (apply mapv + rows)) ms))
(defn- mat-times [m c] (mapv (fn [row] (mapv #(* c %) row)) m))
(defn- mat-vec [m v] (mapv #(reduce + (map * % v)) m))

(defn- quasi-newton-step
  "A line minimization along -B g, B the inverse Hessian built so far,
  and B updated from the change in the gradient -- by BFGS, or DFP when
  `update` is `:dfp`; kept as it was where the change gives no
  curvature. nil where the gradient vanishes."
  [f grad update {:keys [x g B] :as st}]
  (when-not (every? zero? g)
    (let [[_ x' fx'] (line-minimize f x (lin/scale (mat-vec B g) -1.0))
          g' (grad x')
          s (lin/sub x' x) y (lin/sub g' g)
          sy (reduce + (map * s y))
          By (mat-vec B y)
          B' (cond
               (<= sy 0.0) B
               (= update :dfp)
               (mat-sum B (mat-times (outer s s) (/ 1.0 sy))
                        (mat-times (outer By By) (/ -1.0 (reduce + (map * y By)))))
               :else
               ;; B' = (I - rho s y^T) B (I - rho y s^T) + rho s s^T
               (let [rho (/ 1.0 sy)
                     yBy (reduce + (map * y By))]
                 (mat-sum B
                          (mat-times (outer s s) (* rho (+ 1.0 (* rho yBy))))
                          (mat-times (mat-sum (outer By s) (outer s By)) (- rho)))))]
      (moved st x' fx' {:g g' :B B'}))))

(defn- quasi-newton-start [f grad x0]
  (let [n (count x0)] (assoc (start f x0) :g (grad x0) :B (mapv (partial axis n) (range n)))))

(defn quasi-newton-seq
  "The points `{:x :f :iterations :g :B}` of a quasi-Newton method from
  `x0`, `:B` the inverse Hessian built so far. Arguments as
  `quasi-newton`'s."
  ([f grad x0] (quasi-newton-seq f grad x0 {}))
  ([f grad x0 {:keys [update] :or {update :bfgs}}]
   (let [grad (grad-fn f grad)]
     (search/step-seq (partial quasi-newton-step f grad update) (quasi-newton-start f grad x0)))))

(defn quasi-newton
  "The minimum of `f` from `x0` by a quasi-Newton method: Newton's step
  with the inverse Hessian built up from the gradients' changes along the
  way, never computed -- by the BFGS update (Broyden, Fletcher, Goldfarb
  and Shanno 1970, the default) or with `:update :dfp` the
  Davidon-Fletcher-Powell -- each step a line minimization along the
  direction it gives. Arguments as `steepest-descent`'s."
  ([f x0] (quasi-newton f nil x0 {}))
  ([f grad x0] (quasi-newton f grad x0 {}))
  ([f grad x0 {:keys [tol max-iter update] :or {tol 1e-14 max-iter 500 update :bfgs}}]
   (let [grad (grad-fn f grad)]
     (descend tol max-iter (partial quasi-newton-step f grad update) (quasi-newton-start f grad x0)))))

;; ----------------------------------------------------------- constrained

(defn- pivot [T r c]
  (let [prow (T r) pv (prow c)
        prow (mapv #(/ % pv) prow)]
    (mapv (fn [i row] (if (= i r) prow (let [f (row c)] (mapv - row (mapv #(* f %) prow))))) (range (count T)) T)))

(defn- simplex-step
  "One pivot of the simplex method on the tableau `T` (last row the
  objective, to be made non-negative), `basis` the basic variable of each
  constraint row, by Bland's rule -- the lowest-indexed entering and
  leaving variable -- so it cannot cycle. nil at the optimum, where no
  column enters; `:unbounded?` where one enters and nothing bounds it."
  [ncols {:keys [T basis unbounded?]}]
  (when-not unbounded?
    (let [obj (peek T)]
      (when-let [c (first (filter #(< (obj %) -1e-12) (range ncols)))]
        (let [rows (keep (fn [i] (let [a (get-in T [i c])]
                                   (when (> a 1e-12) [(/ (peek (T i)) a) (basis i) i])))
                         (range (dec (count T))))]
          (if (empty? rows)
            {:T T :basis basis :unbounded? true}
            ;; the least ratio, ties -- to a tolerance, rounding being
            ;; what it is -- going to the lowest-indexed basic variable
            (let [least (reduce min (map first rows))
                  [_ _ r] (first (sort-by second (filter #(<= (first %) (+ least (* 1e-11 (max 1.0 (abs least))))) rows)))]
              {:T (pivot T r c) :basis (assoc basis r c)})))))))

(defn- simplex-run
  "The simplex method on tableau `T`, `basis` the basic variable of each
  constraint row, to the optimum: `[T basis]`, or nil if unbounded."
  [T basis ncols]
  (let [{:keys [T basis unbounded?]} (search/settle :unbounded? ##Inf (partial simplex-step ncols) {:T T :basis basis})]
    (when-not unbounded? [T basis])))

(defn- drive-out-artificials
  "After phase one, an artificial variable still basic -- at zero, the
  problem being feasible -- pivoted out on any column of its row below
  `limit` with a nonzero entry; a row with none is redundant, a linear
  combination of the others, and is dropped. Otherwise dropping the
  artificial columns would leave a row whose basic variable is gone,
  as degenerate problems do."
  [T basis limit]
  (reduce (fn [[T basis] i]
            (if (< (basis i) limit)
              [T basis]
              (if-let [c (first (filter #(> (abs (get-in T [i %])) 1e-9) (range limit)))]
                [(pivot T i c) (assoc basis i c)]
                [(assoc T i nil) (assoc basis i nil)])))
          [T basis]
          (range (count basis))))

(defn linear-program
  "The linear program: minimize c . x subject to the rows of `A` against
  `b` -- each `<=`, `>=` or `=` by `:kinds` (default all `<=`) -- and x >=
  0, by the two-phase simplex method (Dantzig 1947; Chapra and Canale
  15.1): slack variables make the inequalities equations, artificial ones
  give a first basic feasible solution, the first phase drives them out,
  the second optimizes. `{:x :value}`; `{:status :infeasible}` or
  `{:status :unbounded}` if there is none. To maximize, pass -c and negate
  the value."
  ([c A b] (linear-program c A b {}))
  ([c A b {:keys [kinds]}]
   (let [m (count b) n (count c)
         kinds (or kinds (repeat m :<=))
         ;; rows with a negative right side are flipped
         [A b kinds] (reduce (fn [[A b ks] [row bi k]]
                               (if (neg? bi)
                                 [(conj A (mapv - row)) (conj b (- bi)) (conj ks ({:<= :>= :>= :<= := :=} k))]
                                 [(conj A (vec row)) (conj b bi) (conj ks k)]))
                             [[] [] []] (map vector A b kinds))
         slack (vec (filter #(not= := (kinds %)) (range m)))
         art (vec (filter #(not= :<= (kinds %)) (range m)))
         ns (count slack) na (count art)
         ncols (+ n ns na)
         row (fn [i]
               (vec (concat (A i)
                            (map #(if (= % i) (if (= :<= (kinds i)) 1.0 -1.0) 0.0) slack)
                            (map #(if (= % i) 1.0 0.0) art)
                            [(double (b i))])))
         T0 (mapv row (range m))
         slack-col (zipmap slack (range)) art-col (zipmap art (range))
         basis0 (mapv (fn [i] (if (= :<= (kinds i))
                                (+ n (slack-col i))
                                (+ n ns (art-col i))))
                      (range m))
         ;; phase one: minimize the artificial variables' sum
         phase1-obj (let [base (vec (concat (repeat (+ n ns) 0.0) (repeat na 1.0) [0.0]))]
                      (reduce (fn [o i] (if (>= (basis0 i) (+ n ns)) (mapv - o (T0 i)) o)) base (range m)))
         [T1 basis1] (if (zero? na) [T0 basis0] (simplex-run (conj T0 phase1-obj) basis0 ncols))]
     (cond
       (and (pos? na) (> (abs (peek (peek T1))) 1e-9)) {:status :infeasible}
       :else
       (let [T1 (if (pos? na) (pop T1) T1)
             keep-cols (+ n ns)
             [T1 basis1] (drive-out-artificials T1 basis1 keep-cols)
             T1 (vec (remove nil? T1)) basis1 (vec (remove nil? basis1))
             m (count T1)
             ;; drop the artificial columns
             T1 (mapv #(conj (subvec % 0 keep-cols) (peek %)) T1)
             obj (vec (concat c (repeat ns 0.0) [0.0]))
             ;; the objective row, reduced by the basic columns
             obj (reduce (fn [o i] (let [j (basis1 i)] (if (< j keep-cols) (mapv - o (mapv #(* (o j) %) (T1 i))) o)))
                         obj (range m))
             res (simplex-run (conj T1 obj) basis1 keep-cols)]
         (if (nil? res)
           {:status :unbounded}
           (let [[T basis] res
                 x (reduce (fn [x i] (let [j (basis i)] (if (< j n) (assoc x j (peek (T i))) x)))
                           (vec (repeat n 0.0)) (range m))]
             {:x x :value (reduce + (map * c x))})))))))

(defn- penalized
  "f plus `mu` times the squared violations of the constraints."
  [f inequalities equalities mu]
  (fn [x] (+ (f x) (* mu (+ (reduce + (map #(let [v (% x)] (if (pos? v) (* v v) 0.0)) inequalities))
                            (reduce + (map #(let [v (% x)] (* v v)) equalities)))))))

(defn- violation [inequalities equalities x]
  (max (reduce max 0.0 (map #(% x) inequalities))
       (reduce max 0.0 (map #(abs (% x)) equalities))))

(defn- penalty-step
  "The penalized problem with `:next-mu` minimized from the last point,
  and mu grown tenfold for the next."
  [f inequalities equalities {:keys [x next-mu iterations]}]
  (let [{x' :x} (quasi-newton (penalized f inequalities equalities next-mu) x)]
    {:x x' :mu next-mu :next-mu (* 10.0 next-mu) :violation (violation inequalities equalities x')
     :iterations (inc iterations)}))

(defn penalty-seq
  "The points `{:x :mu :violation :iterations}` of the exterior penalty
  method from `x0`, each the minimum with penalty `:mu`. Arguments as
  `penalty`'s."
  ([f inequalities equalities x0] (penalty-seq f inequalities equalities x0 {}))
  ([f inequalities equalities x0 {:keys [mu] :or {mu 1.0}}]
   (rest (iterate (partial penalty-step f inequalities equalities) {:x (mapv double x0) :next-mu mu :iterations 0}))))

(defn penalty
  "The minimum of `f` subject to inequality constraints g(x) <= 0 (each
  of `inequalities`) and equality constraints h(x) = 0 (`equalities`),
  from `x0`, by an exterior quadratic penalty (Chapra and Canale 15.2;
  Rao 1996): f + mu (sum max(0, g)^2 + sum h^2) minimized without
  constraints (by `quasi-newton`), then again with mu ten times larger
  from there, until the constraints hold to `:tol` (default 1e-8).
  `{:x :f :mu}`."
  ([f inequalities equalities x0] (penalty f inequalities equalities x0 {}))
  ([f inequalities equalities x0 {:keys [tol mu max-iter] :or {tol 1e-8 mu 1.0 max-iter 20}}]
   (let [{:keys [x mu]} (search/settle #(some-> (:violation %) (<= tol)) (+ 2 max-iter)
                                       (partial penalty-step f inequalities equalities)
                                       {:x (mapv double x0) :next-mu mu :iterations 0})]
     {:x x :f (f x) :mu mu})))
