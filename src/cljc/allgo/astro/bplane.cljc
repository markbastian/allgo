(ns allgo.astro.bplane
  "The B-plane: how a hyperbolic approach misses its target (Kizner, \"A
  method of describing miss distances for lunar and interplanetary
  trajectories\", JPL External Publication 674, 1959, and Planetary and
  Space Science 7, 1961; Vallado, chapter 12).

  An approach hyperbola has an incoming asymptote along S, the direction
  of the arriving v-infinity. The B-plane is the plane through the
  target's center normal to S, and B the vector in it to where the
  asymptote pierces it -- its length the impact parameter, the distance
  by which the trajectory would miss were there no gravity. B is resolved
  along T = S x k/|S x k|, parallel to the reference plane whose pole is k
  (the ecliptic, say, or the target's equator), and R = S x T. B.T and
  B.R are nearly linear in small changes to the approach, which makes
  them the coordinates to target in: `target` finds the velocity change
  that brings them to what is wanted.

  States are `[r v]` relative to the target, km and km/s; `mu` the
  target's."
  (:require [allgo.astro.kepler :as kepler]
            [allgo.astro.universal :as universal]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn b-plane
  "The B-plane of the hyperbolic state `[r v]` about a body of
  gravitational parameter `mu`, the reference pole `k` (default z):
  `{:b :bt :br :s :t :r :v-inf :theta}` -- B, its components B.T and B.R,
  the unit vectors S, T, R, the hyperbolic excess speed and the angle of
  B from T. With e the eccentricity vector and h the angular momentum,
  S = e/e^2 + (sqrt(e^2 - 1)/e) (h x e)/(h e), and B = b S x h/h with b =
  |a| sqrt(e^2 - 1)."
  ([mu s] (b-plane mu s [0.0 0.0 1.0]))
  ([mu [r v] k]
   (let [h (v3/cross r v)
         hhat (v3/normalize h)
         rm (v3/length r)
         evec (v3/scale (v3/sub (v3/scale r (- (v3/dot v v) (/ mu rm))) (v3/scale v (v3/dot r v))) (/ 1.0 mu))
         e (v3/length evec)
         _ (when (<= e 1.0) (throw (ex-info "not a hyperbola" {:e e})))
         ehat (v3/scale evec (/ 1.0 e))
         a (/ 1.0 (- (/ 2.0 rm) (/ (v3/dot v v) mu)))
         ;; the incoming asymptote
         s (v3/add (v3/scale ehat (/ 1.0 e)) (v3/scale (v3/cross hhat ehat) (/ (math/sqrt (- (* e e) 1.0)) e)))
         t (v3/normalize (v3/cross s k))
         rhat (v3/cross s t)
         b (* (abs a) (math/sqrt (- (* e e) 1.0)))
         bvec (v3/scale (v3/cross s hhat) b)
         bt (v3/dot bvec t) br (v3/dot bvec rhat)]
     {:b bvec :bt bt :br br :s s :t t :r rhat
      :v-inf (math/sqrt (/ mu (abs a))) :theta (math/atan2 br bt)})))

(defn periapsis-state
  "The state at periapsis of the hyperbola arriving along `v-inf` (the
  vector, km/s) with B-plane components `bt` and `br` about the pole `k`:
  the inverse of `b-plane`. The periapsis radius is |a| (e - 1), e =
  sqrt(1 + b^2/a^2); periapsis lies at the angle arccos(1/e) from S
  toward B."
  ([mu v-inf bt br] (periapsis-state mu v-inf bt br [0.0 0.0 1.0]))
  ([mu v-inf bt br k]
   (let [vinf (v3/length v-inf)
         s (v3/scale v-inf (/ 1.0 vinf))
         t (v3/normalize (v3/cross s k))
         rhat (v3/cross s t)
         bvec (v3/add (v3/scale t bt) (v3/scale rhat br))
         b (v3/length bvec)
         bhat (v3/scale bvec (/ 1.0 b))
         a (/ mu (* vinf vinf))                               ; |a|
         e (math/sqrt (+ 1.0 (/ (* b b) (* a a))))
         hhat (v3/cross bhat s)
         p (v3/add (v3/scale s (/ 1.0 e)) (v3/scale bhat (/ (math/sqrt (- (* e e) 1.0)) e)))
         rp (* a (- e 1.0))
         vp (math/sqrt (* mu (+ (/ 2.0 rp) (/ 1.0 a))))]
     [(v3/scale p rp) (v3/scale (v3/cross hhat p) vp)])))

(defn periapsis-radius
  "The periapsis radius of a hyperbola of excess speed `v-inf` and impact
  parameter `b`: b^2 = rp^2 (1 + 2 mu/(rp v-inf^2)) solved for rp."
  [mu v-inf b]
  (let [q (/ mu (* v-inf v-inf))]
    (- (math/sqrt (+ (* q q) (* b b))) q)))

(defn target
  "The velocity change at state `s0` that makes the approach, once flown
  by `propagate` (a function of a state and a time, default two-body
  about `mu`) for `dt` seconds, cross the B-plane at `bt` and `br`: the
  least-squares, least-norm solution of the two targets in the three
  components, found by Newton's method with the partials by central
  differences. `{:dv :bt :br :iterations}`, nil if it does not converge."
  ([mu s0 dt bt br] (target mu s0 dt bt br {}))
  ([mu [r0 v0] dt bt br {:keys [k propagate tol] :or {k [0.0 0.0 1.0] tol 1e-8}}]
   (let [fly (or propagate (fn [[r v] dt] (universal/propagate mu r v dt)))
         miss (fn [dv] (let [{:keys [bt br] :as bp} (b-plane mu (fly [r0 (v3/add v0 dv)] dt) k)]
                         [[bt br] bp]))
         wanted [bt br]]
     (loop [dv [0.0 0.0 0.0] i 0]
       (let [[got bp] (miss dv)
             res (mapv - got wanted)]
         (cond
           (< (apply max (map abs res)) tol) {:dv dv :bt (:bt bp) :br (:br bp) :iterations i}
           (> i 30) nil
           :else
           (let [;; the 2 x 3 Jacobian of (B.T, B.R) in the velocity change
                 jac (d/jacobian #(first (miss %)) dv {:steps [1e-7 1e-7 1e-7] :richardson? false})
                 ;; the least-norm step: J^T (J J^T)^-1 res
                 jjt (lin/mat-mul jac (lin/transpose jac))
                 [[p q] [r s]] jjt
                 det (- (* p s) (* q r))
                 y [(/ (- (* s (res 0)) (* q (res 1))) det) (/ (- (* p (res 1)) (* r (res 0))) det)]
                 step (lin/mat-vec (lin/transpose jac) y)]
             (recur (v3/sub dv step) (inc i)))))))))

(defn hyperbolic-time-to-periapsis
  "Seconds until periapsis on the hyperbola `[r v]`, negative once past:
  from the hyperbolic anomaly, t = sqrt(-a^3/mu) (e sinh F - F)."
  [mu [r v]]
  (let [{:keys [a e nu]} (kepler/state->elements mu r v)
        x (* (math/sqrt (/ (- e 1.0) (+ e 1.0))) (math/tan (* 0.5 nu)))
        ;; F = 2 atanh x
        f (math/log (/ (+ 1.0 x) (- 1.0 x)))]
    (- (* (math/sqrt (/ (- (* a a a)) mu)) (- (* e (math/sinh f)) f)))))
