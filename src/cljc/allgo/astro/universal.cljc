(ns allgo.astro.universal
  "Two-body motion in universal variables, for every kind of orbit at once
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapter 2).

  Kepler's equation comes in three versions, one per conic, each with its
  own anomaly and each singular where the next one begins. The universal
  variable chi sidesteps that: it measures progress along any conic --
  sqrt(a) times the change in eccentric anomaly on an ellipse, and the
  analogous quantities on a parabola or hyperbola -- and a single equation
  in it, written with the Stumpff functions c2 and c3, covers them all.
  The same chi then gives the Lagrange coefficients f and g, which carry
  the initial position and velocity forward without ever forming the
  orbital elements. That is what makes it the propagator of choice for
  near-parabolic orbits and for Lambert's problem.

  Units are the caller's: kilometers and seconds with `allgo.astro.
  constants/GM-earth`, as everywhere in the package."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.perturbations :as perturbations]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

(def ^:private small 1e-6)

(defn stumpff
  "`[c2 c3]`, the Stumpff functions of psi: (1 - cos sqrt psi)/psi and
  (sqrt psi - sin sqrt psi)/sqrt psi^3 for positive psi, their hyperbolic
  counterparts for negative, and their limits 1/2 and 1/6 near zero --
  where the closed forms lose every digit to cancellation."
  [psi]
  (cond
    (> psi small)
    (let [s (math/sqrt psi)]
      [(/ (- 1.0 (math/cos s)) psi) (/ (- s (math/sin s)) (* s s s))])
    (< psi (- small))
    (let [s (math/sqrt (- psi))]
      [(/ (- 1.0 (math/cosh s)) psi) (/ (- (math/sinh s) s) (* s s s))])
    :else
    [(- 0.5 (/ psi 24.0)) (- (/ 1.0 6.0) (/ psi 120.0))]))

(defn- initial-chi
  "Vallado's starting guesses: proportional to time on an ellipse, Barker's
  solution on a parabola, and the logarithmic asymptote on a hyperbola."
  [mu r0 v0 dt alpha]
  (let [rm (v3/length r0) rdv (v3/dot r0 v0) smu (math/sqrt mu)]
    (cond
      (> alpha small) (* smu dt alpha)
      (< (abs alpha) small)
      (let [p (/ (v3/length-squared (v3/cross r0 v0)) mu)
            s (* 0.5 (- (/ math/PI 2) (math/atan (* 3.0 (math/sqrt (/ mu (* p p p))) dt))))
            w (math/atan (math/cbrt (math/tan s)))]
        (* (math/sqrt p) (/ 2.0 (math/tan (* 2.0 w)))))
      :else
      (let [a (/ 1.0 alpha)
            sg (math/signum dt)]
        (* sg (math/sqrt (- a))
           (math/log (/ (* -2.0 mu alpha dt)
                        (+ rdv (* sg (math/sqrt (* (- mu) a)) (- 1.0 (* rm alpha)))))))))))

(defn chi
  "`{:chi :psi :c2 :c3 :r}`: the universal variable reached after `dt`
  seconds from `r0` `v0`, the quantities it came with, and the new radius
  -- found by Newton's method on the universal Kepler equation. nil if it
  does not converge."
  [mu r0 v0 dt]
  (let [rm (v3/length r0)
        rdv (/ (v3/dot r0 v0) (math/sqrt mu))
        alpha (- (/ 2.0 rm) (/ (v3/length-squared v0) mu))
        smu-dt (* (math/sqrt mu) dt)]
    (loop [x (initial-chi mu r0 v0 dt alpha) i 0]
      (let [psi (* x x alpha)
            [c2 c3] (stumpff psi)
            r (+ (* x x c2) (* rdv x (- 1.0 (* psi c3))) (* rm (- 1.0 (* psi c2))))
            dx (/ (- smu-dt (* x x x c3) (* rdv x x c2) (* rm x (- 1.0 (* psi c3)))) r)
            x' (+ x dx)]
        (cond
          (< (abs dx) 1e-10)
          (let [psi (* x' x' alpha) [c2 c3] (stumpff psi)]
            {:chi x' :psi psi :c2 c2 :c3 c3
             :r (+ (* x' x' c2) (* rdv x' (- 1.0 (* psi c3))) (* rm (- 1.0 (* psi c2))))})
          (> i 50) nil
          :else (recur x' (inc i)))))))

(defn f-and-g
  "`[f g fdot gdot]`, the Lagrange coefficients that carry a state forward:
  r = f r0 + g v0 and v = fdot r0 + gdot v0. From the universal variable
  and its radius, as `chi` returns them."
  [mu r0-mag dt {:keys [chi psi c2 c3 r]}]
  (let [smu (math/sqrt mu)]
    [(- 1.0 (/ (* chi chi c2) r0-mag))
     (- dt (/ (* chi chi chi c3) smu))
     (/ (* smu chi (- (* psi c3) 1.0)) (* r r0-mag))
     (- 1.0 (/ (* chi chi c2) r))]))

(defn propagate
  "`[r v]` `dt` seconds after `[r0 v0]` on any conic, by universal
  variables (Vallado's KEPLER) -- ellipse, parabola or hyperbola alike,
  and backward for negative `dt`. nil if the iteration fails, as it can
  for a near-parabola far out."
  [mu r0 v0 dt]
  (if (zero? dt)
    [r0 v0]
    (when-let [u (chi mu r0 v0 dt)]
      (let [[f g fd gd] (f-and-g mu (v3/length r0) dt u)]
        [(v3/add (v3/scale r0 f) (v3/scale v0 g))
         (v3/add (v3/scale r0 fd) (v3/scale v0 gd))]))))

(defn f-and-g-series
  "`[f g fdot gdot]` for a short `dt` from the Taylor series in the
  invariants u = mu/r^3, p = (r.v)/r^2 and q = (v.v)/r^2 - u, to the
  sixth power of dt (as Bate, Mueller and White give it): no iteration,
  and good while dt is a small fraction of the period."
  [mu r0 v0 dt]
  (let [r (v3/length r0)
        u (/ mu (* r r r))
        p (/ (v3/dot r0 v0) (* r r))
        q (- (/ (v3/dot v0 v0) (* r r)) u)
        p2 (* p p) p3 (* p2 p) p4 (* p3 p) u2 (* u u) q2 (* q q)
        t2 (* dt dt) t3 (* t2 dt) t4 (* t3 dt) t5 (* t4 dt) t6 (* t5 dt)
        f (+ 1.0 (* -0.5 u t2) (* 0.5 u p t3)
             (* (/ u 24.0) (+ (* -15.0 p2) (* 3.0 q) u) t4)
             (* (/ (* p u) 8.0) (- (* 7.0 p2) (* 3.0 q) u) t5)
             (* (/ u 720.0) (+ (* -945.0 p4) (* 630.0 p2 q) (* -45.0 q2) (* 210.0 p2 u) (* -24.0 q u) (- u2)) t6))
        g (+ dt (* (/ -1.0 6.0) u t3) (* 0.25 u p t4)
             (* (/ u 120.0) (+ (* -45.0 p2) (* 9.0 q) u) t5)
             (* (/ (* p u) 360.0) (- (* 210.0 p2) (* 90.0 q) (* 15.0 u)) t6))
        fd (+ (* (- u) dt) (* 1.5 u p t2)
              (* (/ u 6.0) (+ (* -15.0 p2) (* 3.0 q) u) t3)
              (* (/ (* 5.0 p u) 8.0) (- (* 7.0 p2) (* 3.0 q) u) t4)
              (* (/ u 120.0) (+ (* -945.0 p4) (* 630.0 p2 q) (* -45.0 q2) (* 210.0 p2 u) (* -24.0 q u) (- u2)) t5))
        gd (+ 1.0 (* -0.5 u t2) (* u p t3)
              (* (/ u 24.0) (+ (* -45.0 p2) (* 9.0 q) u) t4)
              (* (/ (* p u) 60.0) (- (* 210.0 p2) (* 90.0 q) (* 15.0 u)) t5))]
    [f g fd gd]))

;; ------------------------------------------------------ time of flight

(defn time-of-flight
  "Seconds to travel from `r0` to `r` on the conic of semilatus rectum `p`
  (Vallado's FINDTOF), the short way round. The semi-major axis follows
  from the two radii, the angle between them and p; a parabola is taken
  where it is infinite, and a vanishing one -- no such orbit -- gives NaN."
  [mu r0 r p]
  (let [r0m (v3/length r0) rm (v3/length r)
        cdn (/ (v3/dot r0 r) (* r0m rm))
        dn (math/acos (am/clamp cdn -1.0 1.0))
        k (* r0m rm (- 1.0 cdn))
        l (+ r0m rm)
        m (* r0m rm (+ 1.0 cdn))
        t1 (* (- (* 2.0 m) (* l l)) p p) t2 (* 2.0 k l p) t3 (* k k)
        den (- (+ t1 t2) t3)]
    ;; a parabola where den vanishes -- next to the size of its terms,
    ;; which run to 1e16 and more, so rounding alone leaves a residue
    (if (< (abs den) (* 1e-12 (max (abs t1) (abs t2) t3)))
      (let [chord (math/sqrt (- (+ (* r0m r0m) (* rm rm)) (* 2.0 r0m rm cdn)))
            s (* 0.5 (+ r0m rm chord))]
        (* (/ 2.0 3.0) (math/sqrt (/ (* s s s) (* 2.0 mu)))
           (- 1.0 (math/pow (/ (- s chord) s) 1.5))))
      (let [a (/ (* m k p) den)
            f (- 1.0 (* (/ rm p) (- 1.0 cdn)))
            g (/ (* r0m rm (math/sin dn)) (math/sqrt (* mu p)))]
        (cond
          (zero? a) ##NaN
          (pos? a)
          (let [fd (* (math/sqrt (/ mu p)) (math/tan (* 0.5 dn))
                      (- (/ (- 1.0 cdn) p) (/ 1.0 r0m) (/ 1.0 rm)))
                cde (- 1.0 (* (/ r0m a) (- 1.0 f)))
                sde (/ (* (- r0m) rm fd) (math/sqrt (* mu a)))
                de (math/atan2 sde cde)]
            (+ g (* (math/sqrt (/ (* a a a) mu)) (- (am/wrap-2pi de) (math/sin de)))))
          :else
          (let [dh (am/acosh (+ 1.0 (* (- f 1.0) (/ r0m a))))]
            (+ g (* (math/sqrt (/ (- (* a a a)) mu)) (- (math/sinh dh) dh)))))))))

;; ------------------------------------------------------ hit the Earth?

(defn hits-earth?
  "Whether the transfer from `r1` `v1` to `r2` `v2`, with `revs` complete
  revolutions between, passes within `pad` km of the Earth's surface
  (Vallado's CHECKHITEARTH): either end already below it, or the perigee
  crossed on the way. Returns the reason as a keyword, or nil."
  ([r1 v1 r2 v2 revs pad] (hits-earth? c/GM-earth r1 v1 r2 v2 revs pad))
  ([mu r1 v1 r2 v2 revs pad]
   (let [rpad (+ c/R-earth pad)]
     (cond
       (or (< (v3/length r1) rpad) (< (v3/length r2) rpad)) :at-an-end
       :else
       (let [h (v3/cross r1 v1)
             p (/ (v3/length-squared h) mu)
             energy (kepler/specific-energy mu r1 v1)
             a (- (/ mu (* 2.0 energy)))
             e (math/sqrt (max 0.0 (- 1.0 (/ p a))))
             rp (/ p (+ 1.0 e))
             moving-in (neg? (v3/dot r1 v1))
             moving-out (pos? (v3/dot r2 v2))]
         (cond
           (and (pos? revs) (< rp rpad) (neg? energy)) :during-revolutions
           ;; perigee between the two points: inbound at the first,
           ;; outbound at the second
           (and moving-in moving-out (< rp rpad)) :through-perigee
           :else nil))))))

;; ---------------------------------------------------- J2 secular drift

(defn propagate-j2
  "`[r v]` `dt` seconds on, with the secular effect of the Earth's
  oblateness (Vallado's PKEPLER): the node regresses and the perigee
  advances at the averaged rates J2 gives,

    dOmega/dt = -k cos i,   domega/dt = k (2 - 5/2 sin^2 i),
    k = 3/2 n J2 R^2 / p^2,

  and the mean anomaly runs at the mean motion. `ndot` and `nddot` -- the
  first and second derivatives of the mean motion, rad/s^2 and rad/s^3 --
  decay the orbit as drag does: the mean anomaly takes their Taylor
  terms, a follows from n = sqrt(mu / a^3), so da = -2a dn / 3n, and e
  shrinks with it at a fixed perigee radius a (1 - e), so de = (1 - e)
  da / a.

  Only the drift: the short-periodic wobble J2 also causes is left out,
  so this is the right tool for where an orbit will be in a week and the
  wrong one for where it will be in ten minutes."
  ([r v dt] (propagate-j2 c/GM-earth r v dt 0.0 0.0))
  ([mu r v dt ndot nddot]
   (let [{:keys [a e i raan argp M]} (kepler/state->elements mu r v)
         n (kepler/mean-motion mu a)
         rates (perturbations/j2-secular mu c/R-earth geo/J2 a e i)
         raan' (am/wrap-2pi (+ raan (* (:raan rates) dt)))
         argp' (am/wrap-2pi (+ argp (* (:argp rates) dt)))
         M' (am/wrap-2pi (+ M (* n dt) (* 0.5 ndot dt dt) (/ (* nddot dt dt dt) 6.0)))
         a' (- a (* (/ (* 2.0 a) (* 3.0 n)) ndot dt))
         e' (- e (* (/ (* 2.0 (- 1.0 e)) (* 3.0 n)) ndot dt))]
     (kepler/elements->state mu {:a a' :e e' :i i :raan raan' :argp argp' :M M'}))))
