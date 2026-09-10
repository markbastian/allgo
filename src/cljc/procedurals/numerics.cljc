(ns procedurals.numerics
  "One door onto the integrators of Montenbruck & Gill chapter 4.

  The families differ in what they need to get going -- a tableau method
  needs nothing, a multistep method a history, an extrapolation method a
  sub-step sequence -- so each has its own constructor. This dispatches on
  the method itself, which is what a caller offering a menu of methods
  actually wants."
  (:require [procedurals.numerics.core :as core]
            [procedurals.numerics.extrapolation :as ex]
            [procedurals.numerics.multistep :as ms]
            [procedurals.numerics.rk :as rk]
            [procedurals.numerics.rkn :as rkn]))

(def first-order
  "Everything that integrates y' = f(t, y), in rough order of cost."
  (vec (concat rk/catalog ms/catalog ex/catalog)))

(def second-order
  "Everything that integrates y'' = f(t, y) directly, without first turning
  it into twice as many first-order equations."
  (vec (concat rkn/catalog ms/catalog-2 ex/catalog-2)))

(defn integrator
  "An integrator for `y' = (f t y)`, whichever family `method` belongs to."
  ([method f t0 y0 h] (integrator method f t0 y0 h {}))
  ([method f t0 y0 h opts]
   (case (:kind method)
     :runge-kutta   (rk/integrator method f t0 y0 h opts)
     :multistep     (ms/integrator method f t0 y0 h opts)
     :extrapolation (ex/integrator method f t0 y0 h opts))))

(defn integrator-2
  "An integrator for `y'' = (f t y)`, whichever family `method` belongs to."
  ([method f t0 y0 dy0 h] (integrator-2 method f t0 y0 dy0 h {}))
  ([method f t0 y0 dy0 h opts]
   (case (:kind method)
     :nystrom          (rkn/integrator method f t0 y0 dy0 h opts)
     :multistep-2      (ms/integrator-2 method f t0 y0 dy0 h opts)
     :extrapolation-2  (ex/integrator-2 method f t0 y0 dy0 h opts))))

(def step core/step)
(def step-until core/step-until)
(def trajectory core/trajectory)

;; ------------------------------------------------------------ test systems

(defn lorenz
  "Lorenz's 1963 convection model. Chaotic at the classical parameters, so
  two trajectories a hair apart separate exponentially -- which makes it a
  searching test of an integrator: any error is amplified, not averaged
  away, and different methods visibly part company."
  ([] (lorenz {}))
  ([{:keys [sigma rho beta] :or {sigma 10.0 rho 28.0 beta (/ 8.0 3.0)}}]
   (fn [_ [x y z]]
     [(* sigma (- y x))
      (- (* x (- rho z)) y)
      (- (* x y) (* beta z))])))

(defn harmonic
  "y'' = -w^2 y, as a second-order system for the Nystrom family."
  ([] (harmonic 1.0))
  ([w] (fn [_ y] (core/v* y (- (* w w))))))

(defn kepler
  "Two-body acceleration toward the origin: the problem the whole chapter is
  written for."
  ([] (kepler 1.0))
  ([mu]
   (fn [_ r]
     (let [d2 (reduce + (map * r r))
           d3 (* d2 (Math/sqrt d2))]
       (core/v* r (- (/ mu d3)))))))

;; Kepler motion has conserved quantities, and watching them decay is a
;; sharper test of an integrator than watching the position: a trajectory
;; can look perfectly elliptical while its energy walks steadily away.

(defn specific-energy
  "v^2/2 - mu/r. Constant on an exact orbit, and equal to -mu/(2a)."
  [mu r v]
  (- (* 0.5 (reduce + (map * v v)))
     (/ mu (Math/sqrt (reduce + (map * r r))))))

(defn angular-momentum
  "r x v, conserved in any central field."
  [[rx ry rz] [vx vy vz]]
  [(- (* ry vz) (* rz vy))
   (- (* rz vx) (* rx vz))
   (- (* rx vy) (* ry vx))])

(defn orbital-period
  "2*pi*sqrt(a^3/mu)."
  [mu a]
  (* 2.0 Math/PI (Math/sqrt (/ (* a a a) mu))))

(defn periapsis-state
  "Position and velocity at closest approach of an ellipse of semi-major
  axis `a` and eccentricity `e`, in the xy plane. Starting at periapsis puts
  the fastest, most sharply curving part of the orbit first, which is where
  a fixed step is most exposed."
  [mu a e]
  (let [rp (* a (- 1.0 e))
        vp (Math/sqrt (/ (* mu (+ 1.0 e)) rp))]
    [[rp 0.0 0.0] [0.0 vp 0.0]]))
