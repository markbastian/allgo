(ns allgo.numerics.rkn
  "Runge-Kutta-Nystrom integrators (Montenbruck & Gill 4.1.2).

  For the special second-order form y'' = f(t, y), where the right-hand side
  does not involve y'. Satellite motion is of this form, and so is any
  system driven purely by position-dependent forces.

  Such a system can always be rewritten as twice as many first-order
  equations, but a Nystrom method integrates it as it stands and needs fewer
  force evaluations for the same order -- the classical fourth-order rule
  below takes three, where RK4 on the doubled system takes four.

    k_i     = f(t + c_i h, y + c_i h y' + h^2 sum_j a_ij k_j)
    y_{n+1} = y + h y' + h^2 sum_i b_i k_i
    y'_{n+1}= y'     + h sum_i b'_i k_i"
  (:require [allgo.numerics.core :as core]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn tableau-step
  "One step from `(t, y, y')`. Returns `[y-next y'-next]`."
  [{:keys [a b b-dot c]} f t y dy h]
  (let [h2 (* h h)
        ks (reduce (fn [ks i]
                     (let [ci  (nth c i)
                           row (nth a i)
                           yi  (cond-> (lin/add y (lin/scale dy (* ci h)))
                                 (seq row) (lin/add (lin/scale (core/combine row ks) h2)))]
                       (conj ks (f (+ t (* h ci)) yi))))
                   []
                   (range (count c)))]
    [(lin/add (lin/add y (lin/scale dy h)) (lin/scale (core/combine b ks) h2))
     (lin/add dy (lin/scale (core/combine b-dot ks) h))]))

;; -------------------------------------------------------------- the tables

(def verlet
  "Velocity Verlet as a Nystrom tableau: second order, two force
  evaluations, and the staple of orbit and molecular-dynamics work for its
  even behavior over long integrations."
  {:name "Verlet" :order 2 :stages 2 :kind :nystrom
   :c [0.0 1.0]
   :a [[] [(/ 1.0 2.0)]]
   :b     [(/ 1.0 2.0) 0.0]
   :b-dot [(/ 1.0 2.0) (/ 1.0 2.0)]})

(def rkn4
  "Nystrom's classical fourth-order rule. RK4 applied to the doubled system
  evaluates the same force twice at the midpoint; collapsing that pair is
  what leaves three stages instead of four."
  {:name "RKN4" :order 4 :stages 3 :kind :nystrom
   :c [0.0 (/ 1.0 2.0) 1.0]
   :a [[] [(/ 1.0 8.0)] [0.0 (/ 1.0 2.0)]]
   :b     [(/ 1.0 6.0) (/ 1.0 3.0) 0.0]
   :b-dot [(/ 1.0 6.0) (/ 2.0 3.0) (/ 1.0 6.0)]})

(def catalog [verlet rkn4])

;; ----------------------------------------------------------- the integrator

(defn- raw [{:keys [method f]} t y dy h]
  (let [[y' dy'] (tableau-step method f t y dy h)]
    {:y y' :dy dy'}))

(defn- advance-fixed [{:keys [t y dy h] :as integ}]
  (let [{y' :y dy' :dy} (raw integ t y dy h)]
    (assoc integ :t (+ t h) :y y' :dy dy' :accepted true :error nil)))

(defn- advance-doubled
  "Adaptive control by step doubling (M&G 4.1.3): compare one step of `h`
  against two of `h/2`. Their difference over 2^p - 1 estimates the error,
  and adding it back is Richardson extrapolation -- a solution one order
  better than the one being controlled, for free.

  This is used rather than an embedded Nystrom pair because it needs no
  second weight vector, and so no coefficient table beyond the one already
  verified above."
  [{:keys [t y dy h control method] :as integ}]
  (let [half   (* 0.5 h)
        coarse (raw integ t y dy h)
        mid    (raw integ t y dy half)
        fine   (raw integ (+ t half) (:y mid) (:dy mid) half)
        p      (:order method)
        denom  (- (math/pow 2.0 p) 1.0)
        diff   (lin/sub (:y fine) (:y coarse))
        e      (lin/scale diff (/ 1.0 denom))
        err    (core/norm e (:y fine) (:tol-abs control) (:tol-rel control))]
    (core/settle integ err p
                 {:t  (+ t h)
                  :y  (lin/add (:y fine) e)
                  :dy (lin/add (:dy fine)
                               (lin/scale (lin/sub (:dy fine) (:dy coarse)) (/ 1.0 denom)))})))

(defn integrator
  "An integrator for `y'' = (f t y)` from `y0`, `dy0` at `t0` with step `h`.
  Pass `:adaptive? true` for step-doubling control."
  ([method f t0 y0 dy0 h] (integrator method f t0 y0 dy0 h {}))
  ([method f t0 y0 dy0 h opts]
   (let [adaptive? (boolean (:adaptive? opts))]
     (core/integrator {:method (assoc method
                                      :advance (if adaptive? advance-doubled advance-fixed)
                                      :adaptive? adaptive?)
                       :f f :t t0 :y y0 :h h :opts opts}
                      {:dy (mapv double dy0)}))))
