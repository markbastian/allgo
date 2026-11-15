(ns allgo.numerics.roots
  "Finding where a function vanishes and where a condition changes:
  bisection and fixed-point iteration (Meeus, *Astronomical Algorithms*,
  chapter 5) and event location -- the slow, sure methods, for when a
  bracket is known -- and the
  rest of the classical root finders (Chapra and Canale, *Numerical
  Methods for Engineers*, chapters 5 and 6): incremental search, false
  position, Newton-Raphson and its version for multiple roots, the secant
  and modified secant methods, Brent's method, and Newton's method and
  fixed-point iteration for systems.

  Each iterative method is a step, as in `allgo.search`: a step from one
  state to the next, and a `-seq` of the states it makes,
  `allgo.search/step-seq` over it. A state is a record with the estimate
  at `:x` and, for the open methods, the step that reached it at `:dx`;
  the bracketing methods keep their bracket too. A step returns nil only
  when the method cannot go on -- a zero derivative, no doubles left
  between the ends -- and never because it has converged: when to stop
  is the caller's, so that the iteration can be watched, tabulated, cut
  short or handed from one method to another. The plain functions,
  `newton` and the rest, run the same steps by `converge` or `settle`,
  which reduce over them rather than build a seq, and return the root."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear-systems :as ls]
            [allgo.search :as search]))

(defn- steps
  "Room for the steps from the first of `n` states to the last: what the
  drivers below reduce over, `n` infinite for no limit."
  [n]
  (if (infinite? n) (repeat nil) (repeat (dec (long n)) nil)))

(defn converge
  "The first state from `init` by `step`, looking at no more than `n`,
  that is `done?` -- or nil if none is, or a step returns nil first.

  The same states as `(take n (search/step-seq step init))`, and the
  same answer as filtering them, but reduced over rather than built as a
  seq: a root finder's step is often only a few arithmetic operations,
  and the seq would cost more than they do."
  [done? n step init]
  (let [s (reduce (fn [s _] (if (done? s) (reduced s) (or (step s) (reduced nil)))) init (steps n))]
    (when (and s (done? s)) s)))

(defn settle
  "The first state from `init` by `step`, looking at no more than `n`,
  that is `done?` -- or else the last one looked at: for the methods that
  are as good as they will get when they run out. Reduced over, as
  `converge` is."
  [done? n step init]
  (reduce (fn [s _] (if (done? s) (reduced s) (or (step s) (reduced s)))) init (steps n)))

;; An estimate, and the step that reached it -- a number, or for a system
;; a vector of them.
(defrecord Estimate [x dx])

(defn- small-step?
  "Whether the step to the estimate is within `tol` of it, relative, the
  estimate taken as at least `floor`."
  [tol floor ^Estimate e]
  (<= (abs (.-dx e)) (* tol (max floor (abs (.-x e))))))

;; ------------------------------------------------------------ bracketing

;; Bisection's bracket, with its midpoint `x`.
(defrecord Bracket [lo hi flo x fx])

(defn- bisect-start [f lo hi]
  (let [x (* 0.5 (+ lo hi))]
    (->Bracket lo hi (f lo) x (f x))))

(defn- bisect-step
  "Halve the bracket toward the sign change: the half whose ends' signs
  differ, and its midpoint. nil at a zero, or once the midpoint is one of
  the ends."
  [f ^Bracket b]
  (let [lo (.-lo b) hi (.-hi b) flo (.-flo b) x (.-x b) fx (.-fx b)]
    (when-not (or (zero? fx) (= x lo) (= x hi))
      (if (= (neg? fx) (neg? flo))
        (let [x' (* 0.5 (+ x hi))] (->Bracket x hi fx x' (f x')))
        (let [x' (* 0.5 (+ lo x))] (->Bracket lo x flo x' (f x')))))))

(defn bisect-seq
  "Bisection's brackets `{:lo :hi :x :fx}`, `:x` the midpoint, from `lo`
  and `hi`, where `f` changes sign: the estimates of Chapra and Canale's
  examples 5.4 and 5.6."
  [f lo hi]
  (search/step-seq (partial bisect-step f) (bisect-start f lo hi)))

(defn bisect
  "A root of `f` between `lo` and `hi`, where it changes sign, by halving
  the interval to the last bit.

  Slow -- one bit per evaluation -- but it cannot fail where a sign change
  is known, which is why Meeus recommends it for Kepler's equation at
  eccentricities near one, where Newton's method may not converge."
  [f lo hi]
  (:x (settle (constantly false) 65 (partial bisect-step f) (bisect-start f lo hi))))

(defn- fixed-point-step
  "One step of x <- (better x)."
  [better ^Estimate e]
  (let [x (.-x e) x' (better x)]
    (->Estimate x' (- x' x))))

(defn fixed-point-seq
  "The estimates `{:x :dx}` of x <- (better x) from `start`: those of
  Chapra and Canale's example 6.1."
  [better start]
  (iterate (partial fixed-point-step better) (->Estimate start ##Inf)))

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
   (:x (converge (partial small-step? tol 1.0) (inc max-iterations)
                 (partial fixed-point-step better) (->Estimate start ##Inf)))))

;; A bracket `[a b]` about where a condition changes, holding at a.
(defrecord Interval [a b])

(defn- crossing-step
  "Halve the bracket toward where `pred` changes. nil once the midpoint is
  one of the ends."
  [pred ^Interval i]
  (let [a (.-a i) b (.-b i) m (* 0.5 (+ a b))]
    (when-not (or (= m a) (= m b))
      (if (pred m) (->Interval m b) (->Interval a m)))))

(defn crossing-seq
  "The brackets `{:a :b}` about where `pred` changes, from `a` (where
  `(pred a)`) and `b` (where not), halved until doubles run out."
  [pred a b]
  (search/step-seq (partial crossing-step pred) (->Interval a b)))

(defn crossing
  "The instant, to within `tol`, between `a` (where `(pred a)`) and `b`
  (where not) that `pred` changes -- or as near as doubles go, when `tol`
  is finer than they are there. The last place found where `pred` still
  holds, so the answer is always on a's side of the change."
  [pred a b tol]
  (:a (settle (fn [^Interval i] (<= (abs (- (.-b i) (.-a i))) tol)) ##Inf
              (partial crossing-step pred) (->Interval a b))))

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

;; False position's bracket, the estimate `x` where its chord crosses
;; zero and the step `dx` to it, and how many steps each end has stuck.
(defrecord FalsePosition [lo hi flo fhi x fx dx stuck-lo stuck-hi])

(defn- chord
  "Where the chord from (lo, flo) to (hi, fhi) crosses zero."
  [lo hi flo fhi]
  (- hi (/ (* fhi (- lo hi)) (- flo fhi))))

(defn- false-position-start [f lo hi]
  (let [flo (f lo) fhi (f hi) x (chord lo hi flo fhi)]
    (->FalsePosition lo hi flo fhi x (f x) ##Inf 0 0)))

(defn- false-position-step
  "Keep the half of the bracket that holds the sign change, and the next
  estimate where its chord crosses zero -- with `modified?`, an end stuck
  for two steps having its function value halved."
  [f modified? ^FalsePosition s]
  (let [lo (.-lo s) hi (.-hi s) flo (.-flo s) fhi (.-fhi s) x (.-x s) fx (.-fx s)]
    (if (neg? (* flo fx))
      ;; the root is in [lo x]: hi moves, lo stays
      (let [stuck-lo (inc (.-stuck-lo s))
            flo (if (and modified? (>= stuck-lo 2)) (* 0.5 flo) flo)
            x' (chord lo x flo fx)]
        (->FalsePosition lo x flo fx x' (f x') (- x' x) stuck-lo 0))
      (let [stuck-hi (inc (.-stuck-hi s))
            fhi (if (and modified? (>= stuck-hi 2)) (* 0.5 fhi) fhi)
            x' (chord x hi fx fhi)]
        (->FalsePosition x hi fx fhi x' (f x') (- x' x) 0 stuck-hi)))))

(defn false-position-seq
  "False position's estimates `{:lo :hi :x :fx :dx}` from `lo` and `hi`,
  where `f` changes sign, `:modified?` as `false-position`'s (default
  true). Unmodified, the estimates of Chapra and Canale's examples 5.6 and
  6.7."
  ([f lo hi] (false-position-seq f lo hi {}))
  ([f lo hi {:keys [modified?] :or {modified? true}}]
   (iterate (partial false-position-step f modified?) (false-position-start f lo hi))))

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
   (:x (settle (fn [^FalsePosition s]
                 (or (zero? (.-fx s)) (<= (abs (.-dx s)) (* tol (max 1e-300 (abs (.-x s)))))))
               (inc max-iter) (partial false-position-step f modified?) (false-position-start f lo hi)))))

;; ---------------------------------------------------------- open methods

(defn- open
  "The root from the open method's `step` from `init`: the first estimate
  whose step is within `:tol` relative, looking no further than
  `:max-iter` steps; nil if none is."
  [step init {:keys [tol max-iter]}]
  (:x (converge (partial small-step? tol 1e-300) (inc max-iter) step init)))

(defn- newton-step
  "One Newton-Raphson step, x <- x - f(x)/f'(x). nil at a zero
  derivative."
  [f df ^Estimate e]
  (let [x (.-x e) d (df x)]
    (when-not (zero? d)
      (let [x' (- x (/ (f x) d))]
        (->Estimate x' (- x' x))))))

(defn newton-seq
  "Newton-Raphson's estimates `{:x :dx}` from `x0`, `df` the derivative
  of `f`: those of Chapra and Canale's examples 6.3, 6.5 and 6.10."
  [f df x0]
  (search/step-seq (partial newton-step f df) (->Estimate x0 ##Inf)))

(defn newton
  "A root of `f` from `x0` by Newton-Raphson, x <- x - f(x)/f'(x), `df`
  the derivative: quadratic convergence near a simple root, none promised
  far from one. Until the step is below `:tol` relative (default 1e-14),
  or nil after `:max-iter` (50) or at a zero derivative."
  ([f df x0] (newton f df x0 {}))
  ([f df x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 50}}]
   (open (partial newton-step f df) (->Estimate x0 ##Inf) {:tol tol :max-iter max-iter})))

(defn- newton-multiple-step
  "One step of Newton-Raphson on f/f', x <- x - f f' / (f'^2 - f f''):
  none at a zero of f, where it stays. nil where the denominator
  vanishes."
  [f df d2f ^Estimate e]
  (let [x (.-x e) fx (f x)]
    (if (zero? fx)
      (->Estimate x 0.0)
      (let [d (df x)
            den (- (* d d) (* fx (d2f x)))]
        (when-not (zero? den)
          (let [x' (- x (/ (* fx d) den))]
            (->Estimate x' (- x' x))))))))

(defn newton-multiple-seq
  "The estimates `{:x :dx}` of Newton-Raphson for multiple roots from
  `x0`: those of Chapra and Canale's example 6.10."
  [f df d2f x0]
  (search/step-seq (partial newton-multiple-step f df d2f) (->Estimate x0 ##Inf)))

(defn newton-multiple
  "A root of `f` of any multiplicity from `x0`: Newton-Raphson on u =
  f/f', which has only simple roots where f has multiple ones -- x <- x -
  f f' / (f'^2 - f f''), `df` and `d2f` the first two derivatives
  (Ralston and Rabinowitz 1978). Plain Newton-Raphson crawls linearly
  onto a double root; this keeps it quadratic. Options as `newton`'s, nil
  where the denominator vanishes."
  ([f df d2f x0] (newton-multiple f df d2f x0 {}))
  ([f df d2f x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 50}}]
   (open (partial newton-multiple-step f df d2f) (->Estimate x0 ##Inf) {:tol tol :max-iter max-iter})))

;; The secant method's last two estimates and the step between them.
(defrecord Secant [xp fp x dx])

(defn- secant-step
  "One secant step through the last two estimates: none at a zero of f,
  or where the two values agree."
  [f ^Secant s]
  (let [xp (.-xp s) fp (.-fp s) x (.-x s) fx (f x)]
    (if (or (zero? fx) (= fx fp))
      (->Secant xp fp x 0.0)
      (let [x' (- x (/ (* fx (- xp x)) (- fp fx)))]
        (->Secant x fx x' (- x' x))))))

(defn secant-seq
  "The secant method's estimates `{:x :dx}` from `x-1` and `x0`: those of
  Chapra and Canale's example 6.6."
  [f x-1 x0]
  (iterate (partial secant-step f) (->Secant x-1 (f x-1) x0 ##Inf)))

(defn secant
  "A root of `f` by the secant method from the two estimates `x-1` and
  `x0`: Newton-Raphson with the derivative replaced by the slope through
  the last two points, which need not bracket the root (and so need not
  converge). Options as `newton`'s."
  ([f x-1 x0] (secant f x-1 x0 {}))
  ([f x-1 x0 {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (:x (converge (fn [^Secant s] (<= (abs (.-dx s)) (* tol (max 1e-300 (abs (.-x s))))))
                 (inc max-iter) (partial secant-step f) (->Secant x-1 (f x-1) x0 ##Inf)))))

(defn- modified-secant-step
  "One Newton-Raphson step with the slope over the perturbation `delta` x.
  nil at a zero slope."
  [f delta ^Estimate e]
  (let [x (.-x e)
        fx (f x)
        dx (* delta (if (zero? x) 1.0 x))
        slope (/ (- (f (+ x dx)) fx) dx)]
    (when-not (zero? slope)
      (let [x' (- x (/ fx slope))]
        (->Estimate x' (- x' x))))))

(defn modified-secant-seq
  "The modified secant method's estimates `{:x :dx}` from `x0`: those of
  Chapra and Canale's example 6.8."
  [f x0 delta]
  (search/step-seq (partial modified-secant-step f delta) (->Estimate x0 ##Inf)))

(defn modified-secant
  "A root of `f` from `x0` by the modified secant method: Newton-Raphson
  with the derivative the slope over a fractional perturbation `delta` of
  x -- one starting point, no derivative. Too small a delta and round-off
  swamps the slope; too large and it converges slowly. Options as
  `newton`'s."
  ([f x0 delta] (modified-secant f x0 delta {}))
  ([f x0 delta {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (open (partial modified-secant-step f delta) (->Estimate x0 ##Inf) {:tol tol :max-iter max-iter})))

;; ------------------------------------------------------------------ Brent

;; Brent's state: `b` the best estimate (also `x`), `a` the one before,
;; `c` the other end of the bracket; `d` the last step and `e` the one
;; before; `by` how the last step was taken.
(defrecord Brent [a fa b fb c fc d e x by])

(defn- brent-arrange
  "Brent's bookkeeping after a new point `b`, `a` the one before: `c` the
  end of the bracket across the sign change from `b` (a fresh bisection
  if that moved), and `b` the better of the two ends."
  [a fa b fb c fc d e by]
  (let [moved? (> (* fb fc) 0.0)
        c (if moved? a c) fc (if moved? fa fc)
        d (if moved? (- b a) d) e (if moved? (- b a) e)]
    (if (< (abs fc) (abs fb))
      (->Brent b fb c fc b fb d e c by)
      (->Brent a fa b fb c fc d e b by))))

(defn- brent-tolerance [tol b] (+ (* 2.0 tol (abs b)) 1e-300))

(defn- brent-done?
  "Whether Brent's bracket is within `tol` of its best end, or that end
  is a zero."
  [tol ^Brent s]
  (let [b (.-b s)]
    (or (zero? (.-fb s)) (<= (abs (* 0.5 (- (.-c s) b))) (brent-tolerance tol b)))))

(defn- brent-step
  "One step of Brent's method: inverse quadratic interpolation through
  the last three points or the secant through two, kept if it falls well
  inside the bracket and shrinks it fast enough, else bisection; `:by`
  says which."
  [f tol ^Brent st]
  (let [a (.-a st) fa (.-fa st) b (.-b st) fb (.-fb st) c (.-c st) fc (.-fc st) d (.-d st) e (.-e st)
        t (brent-tolerance tol b)
        m (* 0.5 (- c b))
        [d e by]
        (if (and (>= (abs e) t) (> (abs fa) (abs fb)))
          (let [s (/ fb fa)
                [p q by] (if (= a c)
                           [(* 2.0 m s) (- 1.0 s) :secant]
                           (let [q (/ fa fc) r (/ fb fc)]
                             [(* s (- (* 2.0 m q (- q r)) (* (- b a) (- r 1.0))))
                              (* (- q 1.0) (- r 1.0) (- s 1.0))
                              :interpolation]))
                [p q] (if (pos? p) [p (- q)] [(- p) q])]
            (if (and (< (* 2.0 p) (- (* 3.0 m q) (abs (* t q))))
                     (< p (abs (* 0.5 e q))))
              [(/ p q) d by]
              [m m :bisection]))
          [m m :bisection])
        b' (+ b (if (> (abs d) t) d (if (pos? m) t (- t))))]
    (brent-arrange b fb b' (f b') c fc d e by)))

(defn- brent-start
  "Brent's first state from the bracket `a` `b`, or nil if `f` does not
  change sign over it."
  [f a b]
  (let [fa (f a) fb (f b)]
    (when (<= (* fa fb) 0.0)
      (brent-arrange a fa b fb a fa (- b a) (- b a) :start))))

(defn brent-seq
  "Brent's states from the bracket `a` `b`, `f` changing sign between
  them: `:b` (and `:x`) the best estimate, `:c` the other end of the
  bracket, `:a` the estimate before, `:by` how the step was taken -- or
  nil without a sign change. `tol` sets the smallest step, as in `brent`."
  [f a b tol]
  (some->> (brent-start f a b) (iterate (partial brent-step f tol))))

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
   (some->> (brent-start f a b) (settle (partial brent-done? tol) (inc max-iter) (partial brent-step f tol)) :b)))

;; --------------------------------------------------------------- systems

(defn- small-steps?
  "Whether every component of the step to the estimate is within `tol` of
  it, relative."
  [tol ^Estimate e]
  (every? true? (map (fn [xi dxi] (<= (abs dxi) (* tol (max 1e-300 (abs xi))))) (.-x e) (.-dx e))))

(defn- system-start [x0] (->Estimate (vec x0) (mapv (constantly ##Inf) x0)))

(defn- fixed-point-system-step
  "One Jacobi step x <- (g x). nil once a component is not a finite
  number."
  [g ^Estimate e]
  (let [x (.-x e) x' (vec (g x))]
    (when (every? #(and (number? %) (not (NaN? %)) (< (abs %) 1e300)) x')
      (->Estimate x' (mapv - x' x)))))

(defn fixed-point-system-seq
  "The estimates `{:x :dx}` of x <- (g x) from the vector `x0`."
  [g x0]
  (search/step-seq (partial fixed-point-system-step g) (system-start x0)))

(defn fixed-point-system
  "The solution of x = (g x), `x0` a vector, by iterating (Jacobi style:
  each new component from the old vector), until the largest change is
  below `:tol` relative (default 1e-12); nil after `:max-iter` (500).
  Converges only when g contracts -- the rearrangement matters, as for one
  equation."
  ([g x0] (fixed-point-system g x0 {}))
  ([g x0 {:keys [tol max-iter] :or {tol 1e-12 max-iter 500}}]
   (:x (converge (partial small-steps? tol) (inc max-iter) (partial fixed-point-system-step g) (system-start x0)))))

(defn- newton-system-step
  "One Newton-Raphson step for a system: solve J dx = -F. nil at a
  singular Jacobian."
  [F J ^Estimate e]
  (let [x (.-x e)]
    (when-let [dx (ls/gauss (J x) (mapv - (F x)))]
      (let [x' (mapv + x dx)]
        (->Estimate x' (mapv - x' x))))))

(defn- difference-jacobian [F x]
  (d/jacobian F x {:steps (mapv (fn [xi] (* 1e-7 (max 1.0 (abs xi)))) x) :richardson? false}))

(defn- jacobian-or-differences [F jacobian]
  (or jacobian (partial difference-jacobian F)))

(defn newton-system-seq
  "The estimates `{:x :dx}` of Newton-Raphson for F(x) = 0 from the
  vector `x0`, `jacobian` as `newton-system`'s: the first that of Chapra
  and Canale's example 6.12."
  [F jacobian x0]
  (search/step-seq (partial newton-system-step F (jacobian-or-differences F jacobian)) (system-start x0)))

(defn newton-system
  "The solution of F(x) = 0, `x0` a vector, by Newton-Raphson for
  systems: each step solves J dx = -F, J the Jacobian -- `jacobian`, a
  function of x giving the matrix, or central differences if nil. Until
  the step is below `:tol` relative (default 1e-13); nil after
  `:max-iter` (50) or at a singular Jacobian."
  ([F x0] (newton-system F nil x0 {}))
  ([F jacobian x0] (newton-system F jacobian x0 {}))
  ([F jacobian x0 {:keys [tol max-iter] :or {tol 1e-13 max-iter 50}}]
   (:x (converge (partial small-steps? tol) (inc max-iter)
                 (partial newton-system-step F (jacobian-or-differences F jacobian)) (system-start x0)))))
