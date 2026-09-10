  (ns procedurals.numerics.rk
    "Explicit Runge-Kutta integrators (Montenbruck & Gill 4.1).

  Every method here is a Butcher tableau -- nodes `c`, stage weights `a`,
  and solution weights `b` -- run by one generic stepper, so a method is
  data and adding one is adding a table.

  An embedded pair carries a second weight vector `b-hat` of lower order.
  The two solutions differ by the leading truncation error, which is what
  makes step-size control possible without integrating anything twice."
    (:require [procedurals.numerics.core :as core]))

;; ------------------------------------------------------------- the stepper

(defn tableau-step
  "One step of size `h` from `(t, y)`. Returns `[y-next error-estimate]`,
  the estimate nil for a tableau with no embedded pair."
  [{:keys [a b b-hat]} c f t y h]
  (let [ks (reduce (fn [ks i]
                     (let [row (nth a i)
                           yi  (if (seq row) (core/v+ y (core/v* (core/combine row ks) h)) y)]
                       (conj ks (f (+ t (* h (nth c i))) yi))))
                   []
                   (range (count c)))]
    [(core/v+ y (core/v* (core/combine b ks) h))
     (when b-hat
       (core/v* (core/combine (mapv - b b-hat) ks) h))]))

;; -------------------------------------------------------------- the tables

(defn- tab [m] (assoc m :kind :runge-kutta))

(def euler
  (tab {:name "Euler" :order 1 :stages 1
        :c [0.0] :a [[]] :b [1.0]}))

(def heun
  (tab {:name "Heun" :order 2 :stages 2
        :c [0.0 1.0] :a [[] [1.0]] :b [0.5 0.5]}))

(def rk3
  (tab {:name "RK3" :order 3 :stages 3
        :c [0.0 0.5 1.0]
        :a [[] [0.5] [-1.0 2.0]]
        :b [(/ 1.0 6.0) (/ 2.0 3.0) (/ 1.0 6.0)]}))

(def rk4
  "The classical fourth-order rule (M&G eq. 4.20)."
  (tab {:name "RK4" :order 4 :stages 4
        :c [0.0 0.5 0.5 1.0]
        :a [[] [0.5] [0.0 0.5] [0.0 0.0 1.0]]
        :b [(/ 1.0 6.0) (/ 1.0 3.0) (/ 1.0 3.0) (/ 1.0 6.0)]}))

(def rkf45
  "Runge-Kutta-Fehlberg 4(5) (M&G table 4.3). The fifth-order solution is
  propagated and the fourth-order one used only to size the error, which is
  local extrapolation: more accurate than the order being controlled."
  (tab {:name "RKF4(5)" :order 5 :error-order 4 :stages 6 :adaptive? true
        :c [0.0 (/ 1.0 4.0) (/ 3.0 8.0) (/ 12.0 13.0) 1.0 (/ 1.0 2.0)]
        :a [[]
            [(/ 1.0 4.0)]
            [(/ 3.0 32.0) (/ 9.0 32.0)]
            [(/ 1932.0 2197.0) (/ -7200.0 2197.0) (/ 7296.0 2197.0)]
            [(/ 439.0 216.0) -8.0 (/ 3680.0 513.0) (/ -845.0 4104.0)]
            [(/ -8.0 27.0) 2.0 (/ -3544.0 2565.0) (/ 1859.0 4104.0) (/ -11.0 40.0)]]
        :b     [(/ 16.0 135.0) 0.0 (/ 6656.0 12825.0) (/ 28561.0 56430.0) (/ -9.0 50.0) (/ 2.0 55.0)]
        :b-hat [(/ 25.0 216.0) 0.0 (/ 1408.0 2565.0) (/ 2197.0 4104.0) (/ -1.0 5.0) 0.0]}))

(def dopri54
  "Dormand-Prince 5(4), the workhorse explicit pair. Its last stage repeats
  the solution point, so a step's final evaluation is the next step's first."
  (tab {:name "DOPRI5(4)" :order 5 :error-order 4 :stages 7 :adaptive? true
        :c [0.0 (/ 1.0 5.0) (/ 3.0 10.0) (/ 4.0 5.0) (/ 8.0 9.0) 1.0 1.0]
        :a [[]
            [(/ 1.0 5.0)]
            [(/ 3.0 40.0) (/ 9.0 40.0)]
            [(/ 44.0 45.0) (/ -56.0 15.0) (/ 32.0 9.0)]
            [(/ 19372.0 6561.0) (/ -25360.0 2187.0) (/ 64448.0 6561.0) (/ -212.0 729.0)]
            [(/ 9017.0 3168.0) (/ -355.0 33.0) (/ 46732.0 5247.0) (/ 49.0 176.0) (/ -5103.0 18656.0)]
            [(/ 35.0 384.0) 0.0 (/ 500.0 1113.0) (/ 125.0 192.0) (/ -2187.0 6784.0) (/ 11.0 84.0)]]
        :b     [(/ 35.0 384.0) 0.0 (/ 500.0 1113.0) (/ 125.0 192.0) (/ -2187.0 6784.0) (/ 11.0 84.0) 0.0]
        :b-hat [(/ 5179.0 57600.0) 0.0 (/ 7571.0 16695.0) (/ 393.0 640.0)
                (/ -92097.0 339200.0) (/ 187.0 2100.0) (/ 1.0 40.0)]}))

(def catalog
  "Every tableau here, lowest order first."
  [euler heun rk3 rk4 rkf45 dopri54])

;; ----------------------------------------------------------- the integrator

(defn- advance [{:keys [method f t y h control] :as integ}]
  (let [[y' e] (tableau-step method (:c method) f t y h)]
    (if-not (:adaptive? method)
      (assoc integ :t (+ t h) :y y' :accepted true :error nil)
      (let [err      (core/norm e y' (:tol-abs control) (:tol-rel control))
            [ok? h'] (core/adapt err (:error-order method) control h)]
        (if ok?
          (assoc integ :t (+ t h) :y y' :h h' :accepted true :error err)
          ;; Rejected: time and state stand, only the step shrinks.
          (assoc integ :h h' :accepted false :error err))))))

(defn integrator
  "An integrator for `y' = (f t y)` starting from `y0` at `t0` with step `h`.

  `opts` may carry any of `core/control-defaults`, and `:adaptive? false` to
  hold an embedded method to a fixed step."
  ([method f t0 y0 h] (integrator method f t0 y0 h {}))
  ([method f t0 y0 h opts]
   {:method   (cond-> (assoc method :advance advance)
                (false? (:adaptive? opts)) (assoc :adaptive? false))
    :f        f
    :t        (double t0)
    :y        (mapv double y0)
    :h        (double h)
    :control  (merge core/control-defaults (dissoc opts :adaptive?))
    :accepted true
    :error    nil}))
