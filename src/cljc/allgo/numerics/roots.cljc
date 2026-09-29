(ns allgo.numerics.roots
  "Finding where a function vanishes, where a condition changes, and
  where a function peaks: bisection and fixed-point iteration (Meeus,
  *Astronomical Algorithms*, chapter 5), event location and golden-section
  search -- the slow, sure methods, for when a bracket is known -- and the
  rest of the classical root finders (Chapra and Canale, *Numerical
  Methods for Engineers*, chapters 5 and 6): incremental search, false
  position, Newton-Raphson and its version for multiple roots, the secant
  and modified secant methods, Brent's method, and Newton's method and
  fixed-point iteration for systems."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear-systems :as ls]
            [clojure.math :as math]))

(defn bisect
  "A root of `f` between `lo` and `hi`, where it changes sign, by halving
  the interval to the last bit.

  Slow -- one bit per evaluation -- but it cannot fail where a sign change
  is known, which is why Meeus recommends it for Kepler's equation at
  eccentricities near one, where Newton's method may not converge."
  [f lo hi]
  (let [flo (f lo)]
    (loop [lo lo hi hi flo flo i 0]
      (let [mid (* 0.5 (+ lo hi))
            fm  (f mid)]
        (cond
          (or (zero? fm) (>= i 64) (= mid lo) (= mid hi)) mid
          (= (neg? fm) (neg? flo)) (recur mid hi fm (inc i))
          :else (recur lo mid flo (inc i)))))))

(defn fixed-point
  "Iterate x <- (better x) from `start` until successive values differ by
  less than `tol`, returning the last. nil after `max-iterations` without
  settling.

  This is how most of astronomy's transcendental equations are solved:
  Kepler's equation, the time of a rising, the light time to a planet --
  each rearranged so the unknown appears once on the left, and fed back
  into itself. It converges whenever the rearranged right side changes
  more slowly than its argument."
  ([better start] (fixed-point better start 1e-15 100))
  ([better start tol max-iterations]
   (loop [x start i 0]
     (when (< i max-iterations)
       (let [n (better x)]
         (if (<= (abs (- n x)) (* tol (max 1.0 (abs n))))
           n
           (recur n (inc i))))))))

(defn crossing
  "The instant, to within `tol`, between `a` (where `(pred a)`) and `b`
  (where not) that `pred` changes."
  [pred a b tol]
  (loop [a a b b]
    (if (<= (abs (- b a)) tol)
      (* 0.5 (+ a b))
      (let [m (* 0.5 (+ a b))]
        (if (pred m) (recur m b) (recur a m))))))

(defn transitions
  "Where `pred` changes between `t0` and `t1`, sampling every `step` and
  refining each change to `tol`: `[[t from-value] ...]`. `step` must be
  shorter than the briefest interval worth finding."
  [pred t0 t1 step tol]
  (let [ts (concat (range t0 t1 step) [t1])]
    (->> (map vector ts (rest ts))
         (keep (fn [[a b]]
                 (let [pa (pred a) pb (pred b)]
                   (when (not= pa pb)
                     [(crossing #(= pa (pred %)) a b tol) pa])))))))

(def ^:private golden (/ (- (math/sqrt 5.0) 1.0) 2.0))

(defn minimize
  "The `x` between `a` and `b` where `f`, which has a single minimum there,
  is least, by golden-section search: until the bracket is narrower than
  `:tol`, or after `:max-iter` narrowings."
  ([f a b] (minimize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (loop [a a b b i 0]
     (if (or (> i max-iter) (< (abs (- b a)) tol))
       (* 0.5 (+ a b))
       (let [x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
         (if (< (f x1) (f x2)) (recur a x2 (inc i)) (recur x1 b (inc i))))))))

(defn maximize
  "The `x` between `a` and `b` where `f`, which has a single maximum there,
  is greatest, by golden-section search; options as `minimize`'s."
  ([f a b] (maximize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (loop [a a b b i 0]
     (if (or (> i max-iter) (< (abs (- b a)) tol))
       (* 0.5 (+ a b))
       (let [x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
         (if (< (f x1) (f x2)) (recur x1 b (inc i)) (recur a x2 (inc i))))))))

;; ------------------------------------------------ bracketing, continued

(defn incremental-search
  "The brackets `[[lo hi] ...]` in which `f` changes sign, found by
  stepping from `a` to `b` in `n` equal increments -- the usual way to
  find starting intervals, which misses pairs of roots closer together
  than a step."
  [f a b n]
  (let [h (/ (- b a) n)
        xs (map #(+ a (* % h)) (range (inc n)))]
    (vec (for [[x0 x1] (partition 2 1 xs)
               :let [f0 (f x0) f1 (f x1)]
               :when (or (zero? f0) (neg? (* f0 f1)))]
           [x0 x1]))))

(defn false-position
  "A root of `f` between `lo` and `hi`, where it changes sign, by false
  position (regula falsi): the next estimate where the chord between the
  bracket's ends crosses zero, the bracket kept. Until the estimate moves
  by less than `:tol` relative (default 1e-12), or `:max-iter` (100).
  With `:modified? true` (the default), a bracket end that has stuck for
  two iterations has its function value halved, which stops the one-sided
  crawl plain false position makes on a strongly curved function."
  ([f lo hi] (false-position f lo hi {}))
  ([f lo hi {:keys [tol max-iter modified?] :or {tol 1e-12 max-iter 100 modified? true}}]
   (loop [lo lo hi hi flo (f lo) fhi (f hi) prev nil stuck-lo 0 stuck-hi 0 i 0]
     (let [x (- hi (/ (* fhi (- lo hi)) (- flo fhi)))
           fx (f x)]
       (if (or (zero? fx) (>= i max-iter)
               (and prev (<= (abs (- x prev)) (* tol (max 1e-300 (abs x))))))
         x
         (if (neg? (* flo fx))
           ;; the root is in [lo x]: hi moves, lo stays
           (let [stuck-lo (inc stuck-lo)]
             (recur lo x (if (and modified? (>= stuck-lo 2)) (* 0.5 flo) flo) fx x stuck-lo 0 (inc i)))
           (let [stuck-hi (inc stuck-hi)]
             (recur x hi fx (if (and modified? (>= stuck-hi 2)) (* 0.5 fhi) fhi) x 0 stuck-hi (inc i)))))))))

;; ---------------------------------------------------------- open methods

(defn- converged? [x x' tol] (<= (abs (- x' x)) (* tol (max 1e-300 (abs x')))))

(defn newton
  "A root of `f` from `x0` by Newton-Raphson, x <- x - f(x)/f'(x), `df`
  the derivative: quadratic convergence near a simple root, none promised
  far from one. Until the step is below `:tol` relative (default 1e-14),
  or nil after `:max-iter` (50) or at a zero derivative."
  ([f df x0] (newton f df x0 {}))
  ([f df x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 50}}]
   (loop [x x0 i 0]
     (when (< i max-iter)
       (let [d (df x)]
         (when-not (zero? d)
           (let [x' (- x (/ (f x) d))]
             (if (converged? x x' tol) x' (recur x' (inc i))))))))))

(defn newton-multiple
  "A root of `f` of any multiplicity from `x0`: Newton-Raphson on u =
  f/f', which has only simple roots where f has multiple ones -- x <- x -
  f f' / (f'^2 - f f''), `df` and `d2f` the first two derivatives
  (Ralston and Rabinowitz 1978). Plain Newton-Raphson crawls linearly
  onto a double root; this keeps it quadratic. Options as `newton`'s."
  ([f df d2f x0] (newton-multiple f df d2f x0 {}))
  ([f df d2f x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 50}}]
   (loop [x x0 i 0]
     (when (< i max-iter)
       (let [fx (f x) d (df x)
             den (- (* d d) (* fx (d2f x)))]
         (if (or (zero? fx) (zero? den))
           x
           (let [x' (- x (/ (* fx d) den))]
             (if (converged? x x' tol) x' (recur x' (inc i))))))))))

(defn secant
  "A root of `f` by the secant method from the two estimates `x-1` and
  `x0`: Newton-Raphson with the derivative replaced by the slope through
  the last two points, which need not bracket the root (and so need not
  converge). Options as `newton`'s."
  ([f x-1 x0] (secant f x-1 x0 {}))
  ([f x-1 x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (loop [xp x-1 fp (f x-1) x x0 i 0]
     (when (< i max-iter)
       (let [fx (f x)]
         (if (or (zero? fx) (= fx fp))
           x
           (let [x' (- x (/ (* fx (- xp x)) (- fp fx)))]
             (if (converged? x x' tol) x' (recur x fx x' (inc i))))))))))

(defn modified-secant
  "A root of `f` from `x0` by the modified secant method: Newton-Raphson
  with the derivative the slope over a fractional perturbation `delta` of
  x -- one starting point, no derivative. Too small a delta and round-off
  swamps the slope; too large and it converges slowly. Options as
  `newton`'s."
  ([f x0 delta] (modified-secant f x0 delta {}))
  ([f x0 delta {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (loop [x x0 i 0]
     (when (< i max-iter)
       (let [fx (f x)
             dx (* delta (if (zero? x) 1.0 x))
             slope (/ (- (f (+ x dx)) fx) dx)]
         (when-not (zero? slope)
           (let [x' (- x (/ fx slope))]
             (if (converged? x x' tol) x' (recur x' (inc i))))))))))

(defn brent
  "A root of `f` between `a` and `b`, where it changes sign, by Brent's
  method (Brent, Algorithms for Minimization without Derivatives, 1973,
  after Dekker 1969): each step inverse quadratic interpolation through
  the last three points, or the secant through two, whichever applies --
  kept only when it falls well inside the bracket and shrinks it fast
  enough, bisection otherwise. As sure as bisection, nearly as fast as the
  open methods. To within `:tol` (default 1e-15 relative, 1e-300
  absolute)."
  ([f a b] (brent f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-15 max-iter 200}}]
   (let [fa (f a) fb (f b)]
     (when (<= (* fa fb) 0.0)
       ;; b the best estimate, a the previous one, c the other end of the
       ;; bracket; d the last step and e the one before
       (loop [a a fa fa b b fb fb c a fc fa d (- b a) e (- b a) i 0]
         (let [[a fa b fb c fc] (if (> (* fb fc) 0.0) [a fa b fb a fa] [a fa b fb c fc])
               [d e] (if (> (* fb fc) 0.0) [(- b a) (- b a)] [d e])
               ;; keep b the better end
               [a fa b fb c fc] (if (< (abs fc) (abs fb)) [b fb c fc b fb] [a fa b fb c fc])
               t (+ (* 2.0 tol (abs b)) 1e-300)
               m (* 0.5 (- c b))]
           (if (or (zero? fb) (<= (abs m) t) (>= i max-iter))
             b
             (let [[d e]
                   (if (and (>= (abs e) t) (> (abs fa) (abs fb)))
                     (let [s (/ fb fa)
                           [p q] (if (= a c)
                                   ;; the secant
                                   [(* 2.0 m s) (- 1.0 s)]
                                   ;; inverse quadratic interpolation
                                   (let [q (/ fa fc) r (/ fb fc)]
                                     [(* s (- (* 2.0 m q (- q r)) (* (- b a) (- r 1.0))))
                                      (* (- q 1.0) (- r 1.0) (- s 1.0))]))
                           [p q] (if (pos? p) [p (- q)] [(- p) q])]
                       (if (and (< (* 2.0 p) (- (* 3.0 m q) (abs (* t q))))
                                (< p (abs (* 0.5 e q))))
                         [(/ p q) d]
                         [m m]))
                     [m m])
                   b' (+ b (if (> (abs d) t) d (if (pos? m) t (- t))))]
               (recur b fb b' (f b') c fc d e (inc i))))))))))

;; --------------------------------------------------------------- systems

(defn fixed-point-system
  "The solution of x = (g x), `x0` a vector, by iterating (Jacobi style:
  each new component from the old vector), until the largest change is
  below `:tol` relative (default 1e-12); nil after `:max-iter` (500).
  Converges only when g contracts -- the rearrangement matters, as for one
  equation."
  ([g x0] (fixed-point-system g x0 {}))
  ([g x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 500}}]
   (loop [x (vec x0) i 0]
     (when (< i max-iter)
       (let [x' (vec (g x))]
         (when (every? #(and (number? %) (not (NaN? %)) (< (abs %) 1e300)) x')
           (if (every? true? (map #(converged? %1 %2 tol) x x'))
             x'
             (recur x' (inc i)))))))))

(defn newton-system
  "The solution of F(x) = 0, `x0` a vector, by Newton-Raphson for
  systems: each step solves J dx = -F, J the Jacobian -- `jacobian`, a
  function of x giving the matrix, or central differences if nil. Until
  the step is below `:tol` relative (default 1e-13); nil after
  `:max-iter` (50) or at a singular Jacobian."
  ([F x0] (newton-system F nil x0 {}))
  ([F jacobian x0] (newton-system F jacobian x0 {}))
  ([F jacobian x0 {:keys [tol max-iter] :or {tol 1e-13 max-iter 50}}]
   (let [J (or jacobian #(d/jacobian F % {:steps (mapv (fn [xi] (* 1e-7 (max 1.0 (abs xi)))) %) :richardson? false}))]
     (loop [x (vec x0) i 0]
       (when (< i max-iter)
         (when-let [dx (ls/gauss (J x) (mapv - (F x)))]
           (let [x' (mapv + x dx)]
             (if (every? true? (map #(converged? %1 %2 tol) x x'))
               x'
               (recur x' (inc i))))))))))
