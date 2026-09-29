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
  interpolation, Newton's method and Brent's method (golden section is
  `allgo.numerics.roots/minimize`); in several, random search, univariate
  and Powell's pattern searches, steepest descent, Fletcher-Reeves
  conjugate gradients, Newton's method, Marquardt's method and the
  quasi-Newton BFGS and DFP updates; the simplex method for linear
  programs; and penalty functions for nonlinear constraints.

  Everything here minimizes; to maximize f, minimize -f."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.random :as random]
            [clojure.math :as math]))

(declare line-minimize)

(defn- replace-worst
  "The simplex with its worst point swapped for `p`, sorted best first."
  [simplex p]
  (vec (sort-by second (conj (pop simplex) p))))

(defn- centroid [points]
  (lin/scale (reduce lin/add points) (/ 1.0 (count points))))

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
                 (and (every? (fn [[x _]] (every? #(<= (abs %) tol-x) (lin/sub x xb))) simplex)
                      (every? (fn [[_ v]] (<= (abs (- v fb)) tol-f)) simplex)))
           {:x xb :f fb :iterations k}
           (let [c (centroid (map first (pop simplex)))
                 at (fn [t] (let [x (lin/add c (lin/scale (lin/sub c xw) t))] [x (g x)]))
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
                                            (let [y (lin/add xb (lin/scale (lin/sub x xb) 0.5))] [y (g y)])))))))
              (inc k)))))))))

;; ---------------------------------------------------------- one variable

(defn parabolic
  "The minimum of `f` by successive parabolic interpolation from three
  points bracketing it, `x0` < `x1` < `x2` with f(x1) below both ends: the
  vertex of the parabola through the three is the next point, and the
  worst end on its side is dropped. Until the vertex moves by less than
  `:tol` (default 1e-10); `x` at the minimum."
  ([f x0 x1 x2] (parabolic f x0 x1 x2 {}))
  ([f x0 x1 x2 {:keys [tol max-iter] :or {tol 1e-10 max-iter 200}}]
   (loop [x0 x0 x1 x1 x2 x2 f0 (f x0) f1 (f x1) f2 (f x2) prev nil i 0]
     (let [num (+ (* f0 (- (* x1 x1) (* x2 x2))) (* f1 (- (* x2 x2) (* x0 x0))) (* f2 (- (* x0 x0) (* x1 x1))))
           den (* 2.0 (+ (* f0 (- x1 x2)) (* f1 (- x2 x0)) (* f2 (- x0 x1))))]
       (if (or (zero? den) (>= i max-iter))
         x1
         (let [x3 (/ num den) f3 (f x3)]
           (if (and prev (<= (abs (- x3 prev)) (* tol (max 1.0 (abs x3)))))
             x3
             ;; keep the three that bracket the lowest
             (let [[x0 x1 x2 f0 f1 f2]
                   (if (> x3 x1)
                     (if (< f3 f1) [x1 x3 x2 f1 f3 f2] [x0 x1 x3 f0 f1 f3])
                     (if (< f3 f1) [x0 x3 x1 f0 f3 f1] [x3 x1 x2 f3 f1 f2]))]
               (recur x0 x1 x2 f0 f1 f2 x3 (inc i))))))))))

(defn newton-1d
  "A stationary point of f -- a minimum where f'' > 0 -- by Newton's
  method on its derivative `df`, `d2f` the second derivative: x <- x -
  f'(x)/f''(x). From `x0`; nil if it does not settle in `:max-iter`
  (default 50)."
  ([df d2f x0] (newton-1d df d2f x0 {}))
  ([df d2f x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 50}}]
   (loop [x x0 i 0]
     (when (< i max-iter)
       (let [d2 (d2f x)]
         (when-not (zero? d2)
           (let [x' (- x (/ (df x) d2))]
             (if (<= (abs (- x' x)) (* tol (max 1.0 (abs x')))) x' (recur x' (inc i))))))))))

(def ^:private cgold (* 0.5 (- 3.0 (math/sqrt 5.0))))

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
   (let [x0 (+ a (* cgold (- b a))) f0 (f x0)]
     ;; x the best point, w the second best, v the previous w
     (loop [a a b b x x0 w x0 v x0 fx f0 fw f0 fv f0 d 0.0 e 0.0 i 0]
       (let [m (* 0.5 (+ a b))
             t1 (+ (* tol (abs x)) 1e-12)
             t2 (* 2.0 t1)]
         (if (or (<= (abs (- x m)) (- t2 (* 0.5 (- b a)))) (>= i max-iter))
           {:x x :f fx}
           (let [[d e]
                 (if (> (abs e) t1)
                   ;; try the parabola through x, w and v
                   (let [r (* (- x w) (- fx fv))
                         q (* (- x v) (- fx fw))
                         p (- (* (- x v) q) (* (- x w) r))
                         q (* 2.0 (- q r))
                         [p q] (if (pos? q) [(- p) q] [p (- q)])]
                     (if (and (< (abs p) (abs (* 0.5 q e))) (> p (* q (- a x))) (< p (* q (- b x))))
                       ;; the parabola's step taken; the one before it remembered
                       (let [step (/ p q) u (+ x step)]
                         [(if (or (< (- u a) t2) (< (- b u) t2)) (if (< x m) t1 (- t1)) step) d])
                       (let [e (if (< x m) (- b x) (- a x))] [(* cgold e) e])))
                   (let [e (if (< x m) (- b x) (- a x))] [(* cgold e) e]))
                 u (+ x (if (>= (abs d) t1) d (if (pos? d) t1 (- t1))))
                 fu (f u)]
             (if (<= fu fx)
               (let [[a b] (if (< u x) [a x] [x b])]
                 (recur a b u x w fu fx fw d e (inc i)))
               (let [[a b] (if (< u x) [u b] [a u])]
                 (cond
                   (or (<= fu fw) (= w x)) (recur a b x u w fx fu fw d e (inc i))
                   (or (<= fu fv) (= v x) (= v w)) (recur a b x w u fx fw fu d e (inc i))
                   :else (recur a b x w v fx fw fv d e (inc i))))))))))))

;; ---------------------------------------------------------- line search

(defn line-minimize
  "The step t along direction `d` from `x` that minimizes f(x + t d):
  bracketed by stepping out from 0 with growing steps, then found by
  Brent's method. `[t x+td f]`."
  [f x d]
  (let [g (fn [t] (f (lin/add x (lin/scale d t))))
        f0 (g 0.0)
        h (/ 1e-2 (max 1e-12 (math/sqrt (reduce + (map * d d)))))
        ;; walk downhill from 0, doubling, until f rises
        [lo hi] (if (< (g h) f0)
                  (loop [a 0.0 b h fb (g h)]
                    (let [c (* 2.0 b) fc (g c)]
                      (if (or (>= fc fb) (> c 1e12)) [a c] (recur b c fc))))
                  (loop [a (- h) fb f0 k 0]
                    (let [fa (g a)]
                      (if (or (>= fa fb) (> k 60)) [a h] (recur (* 2.0 a) fa (inc k))))))
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

(defn univariate
  "The minimum of `f` from `x0` by univariate search: a line minimization
  along each coordinate in turn, round and round -- simple, and slow in a
  narrow valley that runs across the axes. `{:x :f :iterations}`."
  ([f x0] (univariate f x0 {}))
  ([f x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 1000}}]
   (let [n (count x0)]
     (loop [x (mapv double x0) fx (f x0) k 1]
       (let [[x' fx'] (reduce (fn [[x _] i]
                                (let [[_ x' fv] (line-minimize f x (assoc (vec (repeat n 0.0)) i 1.0))]
                                  [x' fv]))
                              [x fx] (range n))]
         (if (or (settled? x x' fx fx' tol) (>= k max-iter))
           {:x x' :f fx' :iterations k}
           (recur x' fx' (inc k))))))))

(defn powell
  "The minimum of `f` from `x0` by Powell's method (Powell 1964): line
  minimizations along a set of directions, starting with the axes; after
  each round the net displacement becomes a new direction and the oldest
  is dropped. The directions it builds are conjugate, and a quadratic in n
  variables is minimized in n rounds. `{:x :f :iterations}`."
  ([f x0] (powell f x0 {}))
  ([f x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 500}}]
   (let [n (count x0)]
     (loop [x (mapv double x0) fx (f x0)
            dirs (mapv (fn [i] (assoc (vec (repeat n 0.0)) i 1.0)) (range n)) k 1]
       (let [[y fy] (reduce (fn [[y _] d] (let [[_ y' fv] (line-minimize f y d)] [y' fv])) [x fx] dirs)
             net (lin/sub y x)
             [x' fx' dirs] (if (every? zero? net)
                             [y fy dirs]
                             (let [[_ z fz] (line-minimize f y net)]
                               [z fz (conj (subvec dirs 1) net)]))]
         (if (or (settled? x x' fx fx' tol) (>= k max-iter))
           {:x x' :f fx' :iterations k}
           (recur x' fx' dirs (inc k))))))))

;; -------------------------------------------- several variables: gradients

(defn- grad-fn [f grad] (or grad #(d/gradient f %)))

(defn steepest-descent
  "The minimum of `f` from `x0` by steepest descent: a line minimization
  along the negative gradient, again and again -- `grad` its gradient
  function, or central differences if nil. Each step turns a right angle
  from the last, so it zigzags down narrow valleys. `{:x :f
  :iterations}`."
  ([f x0] (steepest-descent f nil x0 {}))
  ([f grad x0] (steepest-descent f grad x0 {}))
  ([f grad x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 5000}}]
   (let [grad (grad-fn f grad)]
     (loop [x (mapv double x0) fx (f x0) k 1]
       (let [g (grad x)]
         (if (every? zero? g)
           {:x x :f fx :iterations k}
           (let [[_ x' fx'] (line-minimize f x (lin/scale g -1.0))]
             (if (or (settled? x x' fx fx' tol) (>= k max-iter))
               {:x x' :f fx' :iterations k}
               (recur x' fx' (inc k))))))))))

(defn conjugate-gradient
  "The minimum of `f` from `x0` by the Fletcher-Reeves conjugate gradient
  method: each search direction the negative gradient plus the last
  direction times |g_new|^2 / |g_old|^2, restarted down the gradient
  every n steps -- a quadratic minimized in n line searches, at the cost
  of steepest descent's storage. Arguments as `steepest-descent`'s."
  ([f x0] (conjugate-gradient f nil x0 {}))
  ([f grad x0] (conjugate-gradient f grad x0 {}))
  ([f grad x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 5000}}]
   (let [grad (grad-fn f grad) n (count x0)]
     (loop [x (mapv double x0) fx (f x0) g (grad x0) d (lin/scale (grad x0) -1.0) k 1]
       (if (every? zero? g)
         {:x x :f fx :iterations k}
         (let [[_ x' fx'] (line-minimize f x d)
               g' (grad x')
               beta (/ (reduce + (map * g' g')) (reduce + (map * g g)))
               d' (if (zero? (mod k n)) (lin/scale g' -1.0) (lin/add (lin/scale g' -1.0) (lin/scale d beta)))]
           (if (or (settled? x x' fx fx' tol) (>= k max-iter))
             {:x x' :f fx' :iterations k}
             (recur x' fx' g' d' (inc k)))))))))

(defn- solve-sym [H g] (ls/gauss H g))

(defn newton
  "The minimum of `f` from `x0` by Newton's method: each step -H^-1 g, H
  the Hessian and g the gradient (`grad` and `hess` functions of x, or
  central differences if nil) -- quadratic convergence near the minimum,
  none promised away from it, where H may not even be positive definite.
  `{:x :f :iterations}`."
  ([f x0] (newton f nil nil x0 {}))
  ([f grad hess x0] (newton f grad hess x0 {}))
  ([f grad hess x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (let [grad (grad-fn f grad) hess (or hess #(d/hessian f %))]
     (loop [x (mapv double x0) fx (f x0) k 1]
       (if-let [dx (solve-sym (hess x) (grad x))]
         (let [x' (lin/sub x dx) fx' (f x')]
           (if (or (settled? x x' fx fx' tol) (>= k max-iter))
             {:x x' :f fx' :iterations k}
             (recur x' fx' (inc k))))
         {:x x :f fx :iterations k})))))

(defn marquardt
  "The minimum of `f` from `x0` by Marquardt's method: each step -(H +
  alpha I)^-1 g, alpha shrunk after a step that lowers f and grown after
  one that doesn't -- steepest descent's reliability far from the
  minimum, Newton's speed near it. `:alpha` starts it (default 1e3).
  Arguments as `newton`'s."
  ([f x0] (marquardt f nil nil x0 {}))
  ([f grad hess x0] (marquardt f grad hess x0 {}))
  ([f grad hess x0 {:keys [tol max-iter alpha] :or {tol 1e-14 max-iter 500 alpha 1e3}}]
   (let [grad (grad-fn f grad) hess (or hess #(d/hessian f %)) n (count x0)]
     (loop [x (mapv double x0) fx (f x0) alpha alpha k 1]
       (let [g (grad x) H (hess x)
             damped (mapv (fn [i row] (update row i + alpha)) (range n) H)
             dx (solve-sym damped g)]
         (if (nil? dx)
           (recur x fx (* 10.0 alpha) (inc k))
           (let [x' (lin/sub x dx) fx' (f x')]
             (cond
               (>= k max-iter) {:x x :f fx :iterations k}
               (< fx' fx) (if (settled? x x' fx fx' tol)
                            {:x x' :f fx' :iterations k}
                            (recur x' fx' (* 0.1 alpha) (inc k)))
               (settled? x x' fx fx' tol) {:x x :f fx :iterations k}
               :else (recur x fx (* 10.0 alpha) (inc k))))))))))

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
   (let [grad (grad-fn f grad) n (count x0)
         outer (fn [u v] (mapv (fn [a] (mapv #(* a %) v)) u))
         mat-add (fn [& ms] (apply mapv (fn [& rows] (apply mapv + rows)) ms))
         mat-scale (fn [m c] (mapv (fn [row] (mapv #(* c %) row)) m))
         mv (fn [m v] (mapv #(reduce + (map * % v)) m))]
     (loop [x (mapv double x0) fx (f x0) g (grad x0)
            B (mapv (fn [i] (assoc (vec (repeat n 0.0)) i 1.0)) (range n)) k 1]
       (if (every? zero? g)
         {:x x :f fx :iterations k}
         (let [[_ x' fx'] (line-minimize f x (lin/scale (mv B g) -1.0))
               g' (grad x')
               s (lin/sub x' x) y (lin/sub g' g)
               sy (reduce + (map * s y))
               By (mv B y)
               B' (cond
                    (<= sy 0.0) B
                    (= update :dfp)
                    (mat-add B (mat-scale (outer s s) (/ 1.0 sy))
                             (mat-scale (outer By By) (/ -1.0 (reduce + (map * y By)))))
                    :else
                    ;; B' = (I - rho s y^T) B (I - rho y s^T) + rho s s^T
                    (let [rho (/ 1.0 sy)
                          yBy (reduce + (map * y By))]
                      (mat-add B
                               (mat-scale (outer s s) (* rho (+ 1.0 (* rho yBy))))
                               (mat-scale (mat-add (outer By s) (outer s By)) (- rho)))))]
           (if (or (settled? x x' fx fx' tol) (>= k max-iter))
             {:x x' :f fx' :iterations k}
             (recur x' fx' g' B' (inc k)))))))))

;; ----------------------------------------------------------- constrained

(defn- pivot [T r c]
  (let [prow (T r) pv (prow c)
        prow (mapv #(/ % pv) prow)]
    (mapv (fn [i row] (if (= i r) prow (let [f (row c)] (mapv - row (mapv #(* f %) prow))))) (range (count T)) T)))

(defn- simplex-run
  "Runs the simplex method on tableau `T` (last row the objective, to be
  made non-negative), `basis` the basic variable of each constraint row.
  Bland's rule -- the lowest-indexed entering and leaving variable -- so
  it cannot cycle. `[T basis]`, or nil if unbounded."
  [T basis ncols]
  (loop [T T basis basis]
    (let [obj (peek T)
          c (first (filter #(< (obj %) -1e-12) (range ncols)))]
      (if (nil? c)
        [T basis]
        (let [rows (keep (fn [i] (let [a (get-in T [i c])]
                                   (when (> a 1e-12) [(/ (peek (T i)) a) (basis i) i])))
                         (range (dec (count T))))]
          (when (seq rows)
            (let [[_ _ r] (first (sort rows))]
              (recur (pivot T r c) (assoc basis r c)))))))))

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
         basis0 (mapv (fn [i] (if (= :<= (kinds i))
                                (+ n (.indexOf slack i))
                                (+ n ns (.indexOf art i))))
                      (range m))
         ;; phase one: minimize the artificial variables' sum
         phase1-obj (let [base (vec (concat (repeat (+ n ns) 0.0) (repeat na 1.0) [0.0]))]
                      (reduce (fn [o i] (if (>= (basis0 i) (+ n ns)) (mapv - o (T0 i)) o)) base (range m)))
         [T1 basis1] (if (zero? na) [T0 basis0] (simplex-run (conj T0 phase1-obj) basis0 ncols))]
     (cond
       (and (pos? na) (> (abs (peek (peek T1))) 1e-9)) {:status :infeasible}
       :else
       (let [T1 (if (pos? na) (pop T1) T1)
             ;; drop the artificial columns
             keep-cols (+ n ns)
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
   (loop [x (mapv double x0) mu mu k 0]
     (let [P (fn [x] (+ (f x) (* mu (+ (reduce + (map #(let [v (% x)] (if (pos? v) (* v v) 0.0)) inequalities))
                                       (reduce + (map #(let [v (% x)] (* v v)) equalities))))))
           {x' :x} (quasi-newton P x)
           violation (max (reduce max 0.0 (map #(% x') inequalities))
                          (reduce max 0.0 (map #(abs (% x')) equalities)))]
       (if (or (<= violation tol) (>= k max-iter))
         {:x x' :f (f x') :mu mu}
         (recur x' (* 10.0 mu) (inc k)))))))
