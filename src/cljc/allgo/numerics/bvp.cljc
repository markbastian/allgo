(ns allgo.numerics.bvp
  "Two-point boundary-value problems, y'' = F(x, y, y') with y(a) and
  y(b) given (Chapra and Canale, *Numerical Methods for Engineers*,
  27.1): by shooting -- guessing the missing slope y'(a), integrating as
  an initial-value problem, and correcting the guess until the far end
  lands on y(b), by linear interpolation between two shots when the
  equation is linear and by the secant method when not -- or by finite
  differences, the derivatives replaced by centered differences on a grid
  and the resulting equations solved together: a tridiagonal system when
  the equation is linear, Newton's method when not.

  Solutions are `[[x y] ...]` at the grid points, both ends included."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.rk :as rk]
            [allgo.numerics.roots :as roots]
            [allgo.numerics.tridiagonal :as tri]))

(defn- fly
  "The ODE y'' = F(x, y, y') from y(a) = ya, y'(a) = z0 to b in `n`
  fourth-order Runge-Kutta steps: `[[x y z] ...]`."
  [F a b ya z0 n]
  (let [h (/ (- b a) n)
        f (fn [x [y z]] [z (F x y z)])]
    (vec (reductions (fn [[x y z] _]
                       (let [[[y' z']] (rk/tableau-step rk/rk4 (:c rk/rk4) f x [y z] h)]
                         [(+ x h) y' z']))
                     [(double a) (double ya) (double z0)] (range n)))))

(defn shoot-linear
  "The linear boundary-value problem y'' = F(x, y, y') between y(`a`) =
  `ya` and y(`b`) = `yb`, by shooting: two trial slopes `z1` `z2` at a,
  each integrated to b in `n` fourth-order Runge-Kutta steps, and -- the
  end value being linear in the slope -- the slope that hits yb found by
  linear interpolation between them. `{:slope :solution}`."
  [F a b ya yb [z1 z2] n]
  (let [end (fn [z] (second (peek (fly F a b ya z n))))
        e1 (end z1) e2 (end z2)
        z (+ z1 (* (- z2 z1) (/ (- yb e1) (- e2 e1))))]
    {:slope z :solution (mapv (fn [[x y]] [x y]) (fly F a b ya z n))}))

(defn shoot
  "The boundary-value problem y'' = F(x, y, y'), linear or not, by
  shooting: the slope at a that lands the far end on yb, found by the
  secant method from the trial slopes `z1` `z2`. `{:slope :solution}`,
  nil if it does not converge."
  [F a b ya yb [z1 z2] n]
  (when-let [z (roots/secant #(- (second (peek (fly F a b ya % n))) yb) z1 z2 {:tol 1e-13})]
    {:slope z :solution (mapv (fn [[x y]] [x y]) (fly F a b ya z n))}))

(defn finite-difference-linear
  "The linear problem y'' = p(x) y' + q(x) y + r(x), y(a) = ya, y(b) =
  yb, on `n` equal intervals: at each interior point (y_i+1 - 2y_i +
  y_i-1)/dx^2 = p (y_i+1 - y_i-1)/(2 dx) + q y_i + r, a tridiagonal system
  solved by the Thomas algorithm."
  [p q r a b ya yb n]
  (let [dx (/ (- b a) n)
        xs (mapv #(+ a (* % dx)) (range (inc n)))
        interior (subvec xs 1 n)
        rows (mapv (fn [x] [(+ (/ 1.0 (* dx dx)) (/ (p x) (* 2.0 dx)))
                            (- (/ -2.0 (* dx dx)) (q x))
                            (- (/ 1.0 (* dx dx)) (/ (p x) (* 2.0 dx)))
                            (r x)])
                   interior)
        ;; the known ends move to the right side
        rhs (-> (mapv #(nth % 3) rows)
                (update 0 - (* (first (first rows)) ya))
                (update (dec (count rows)) - (* (nth (peek rows) 2) yb)))
        ys (tri/solve-vectors (map first rows) (map second rows) (map #(nth % 2) rows) rhs)]
    (mapv vector xs (concat [(double ya)] ys [(double yb)]))))

(defn finite-difference
  "The problem y'' = F(x, y, y'), linear or not, on `n` equal intervals:
  the centered-difference equations at the interior points solved
  together by Newton's method (the Jacobian by central differences), from
  the straight line between the ends. nil if Newton fails."
  [F a b ya yb n]
  (let [dx (/ (- b a) n)
        xs (mapv #(+ a (* % dx)) (range (inc n)))
        residual (fn [ys]
                   (let [full (vec (concat [ya] ys [yb]))]
                     (mapv (fn [i]
                             (let [yl (full (dec i)) y (full i) yr (full (inc i))]
                               (- (/ (+ yl (* -2.0 y) yr) (* dx dx))
                                  (F (xs i) y (/ (- yr yl) (* 2.0 dx))))))
                           (range 1 n))))
        start (mapv #(+ ya (* (- yb ya) (/ % n))) (range 1 n))]
    (loop [ys start k 0]
      (when (< k 50)
        (let [r (residual ys)
              J (d/jacobian residual ys {:steps (mapv #(* 1e-6 (max 1.0 (abs %))) ys) :richardson? false})]
          (when-let [dy (ls/gauss J (mapv - r))]
            (let [ys' (mapv + ys dy)]
              (if (< (apply max (map #(abs %) dy)) (* 1e-12 (max 1.0 (apply max (map abs ys')))))
                (mapv vector xs (concat [(double ya)] ys' [(double yb)]))
                (recur ys' (inc k))))))))))
