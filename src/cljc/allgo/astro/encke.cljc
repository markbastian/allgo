(ns allgo.astro.encke
  "Encke's method of special perturbations (Vallado, *Fundamentals of
  Astrodynamics and Applications*, chapter 8; Battin, *An Introduction to
  the Mathematics and Methods of Astrodynamics*, 1987, section 9.3).

  Cowell's method integrates the whole motion, the central attraction
  and its small perturbations together, so the step must follow the
  orbit's curvature. Encke's integrates only the deviation from a
  reference conic that the two-body solution carries exactly. The
  deviation is small and changes slowly, so the steps can be long -- until
  it grows, when the reference is rectified: restarted as the osculating
  orbit of the moment, and the deviation set back to zero.

  With r = rho + delta, rho the reference conic,

    delta'' = -(mu/rho^3) (f(q) r + delta) + a_p

  where q = delta.(delta - 2r)/r^2 and f(q) = (1+q)^3/2 - 1 is taken as
  q (3 + 3q + q^2) / (1 + (1+q)^3/2), Battin's form, free of the
  cancellation in subtracting two nearly equal numbers.

  km, km/s, seconds."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.universal :as universal]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]))

(defn f-of-q
  "(1 + q)^3/2 - 1, in Battin's form."
  [q]
  (let [s (math/pow (+ 1.0 q) 1.5)]
    (/ (* q (+ 3.0 (* 3.0 q) (* q q))) (+ 1.0 s))))

(defn- deviation-rate
  "d/dt of [delta delta'] for the reference state `[r0 v0]` at `t0`."
  [mu accel [r0 v0] t0]
  (fn [t y]
    (let [delta (subvec y 0 3) ddelta (subvec y 3 6)
          [rho vrho] (universal/propagate mu r0 v0 (- t t0))
          r (v3/add rho delta)
          v (v3/add vrho ddelta)
          q (/ (v3/dot delta (v3/sub delta (v3/scale r 2.0))) (v3/dot r r))
          rho3 (math/pow (v3/length rho) 3)
          a (v3/add (v3/scale (v3/add (v3/scale r (f-of-q q)) delta) (- (/ mu rho3)))
                    (accel t r v))]
      (into ddelta a))))

(defn propagate
  "The state `[r v]` `t` seconds after `[r0 v0]`, under the central body of
  parameter `mu` and the perturbing acceleration `(accel t r v)`, by
  Encke's method: `{:state :rectifications :steps}`. Options:

    :rectify  rectify when |delta| exceeds this fraction of |rho| (1e-3)
    :h        the first step, seconds (60)
    :tol-abs  :tol-rel  the integrator's error control on the deviation

  The deviation is integrated by Dormand and Prince's 5(4) pair."
  ([r0 v0 accel t] (propagate c/GM-earth r0 v0 accel t {}))
  ([mu r0 v0 accel t {:keys [rectify h tol-abs tol-rel] :or {rectify 1e-3 h 60.0 tol-abs 1e-9 tol-rel 1e-12}}]
   (let [zero [0.0 0.0 0.0 0.0 0.0 0.0]
         start (fn [ref t0 h]
                 (rk/integrator rk/dopri54 (deviation-rate mu accel ref t0) t0 zero h
                                {:tol-abs tol-abs :tol-rel tol-rel}))
         state (fn [[r0 v0] t0 {:keys [t y]}]
                 (let [[rho vrho] (universal/propagate mu r0 v0 (- t t0))]
                   [(v3/add rho (subvec y 0 3)) (v3/add vrho (subvec y 3 6))]))]
     (loop [ref [r0 v0] t0 0.0 integ (start [r0 v0] 0.0 h) rects 0 steps 0]
       (let [remaining (- t (:t integ))]
         (if (<= remaining (* 1e-12 (max 1.0 (abs t))))
           {:state (state ref t0 integ) :rectifications rects :steps steps}
           (let [integ' (core/step (assoc integ :h (min (:h integ) remaining)))]
             (if-not (:accepted integ')
               (recur ref t0 integ' rects steps)
               (let [[r :as now] (state ref t0 integ')
                     delta (v3/length (subvec (:y integ') 0 3))]
                 (if (> delta (* rectify (v3/length r)))
                   (recur now (:t integ') (start now (:t integ') (:h integ')) (inc rects) (inc steps))
                   (recur ref t0 integ' rects (inc steps))))))))))))
