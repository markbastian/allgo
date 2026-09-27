(ns allgo.astro.interplanetary
  "Getting from one planet to another by patched conics (Vallado,
  *Fundamentals of Astrodynamics and Applications*, chapter 12; the method
  is the standard one of Bate, Mueller and White and of every mission
  design text since).

  The trip is cut into conics that each feel one body: a hyperbola
  leaving the first planet, an ellipse about the Sun, a hyperbola arriving
  at the second. They are patched at each planet's sphere of influence,
  which is small enough against the heliocentric orbit that the
  heliocentric leg can treat the planets as points and large enough that
  the planetocentric legs can treat the Sun as far away. What passes
  between them is the hyperbolic excess velocity v-infinity: the
  planetocentric speed at the sphere's edge, which is the heliocentric
  velocity less the planet's.

  km, km/s, seconds, radians; heliocentric vectors in equatorial J2000."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.iod :as iod]
            [allgo.astro.planets :as planets]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(def ^:private mu-sun c/GM-sun)

;; ------------------------------------------------------ spheres of influence

(defn sphere-of-influence
  "Laplace's sphere of influence of a body of parameter `mu` orbiting one
  of `mu-primary` at distance `a`: a (mu/mu-primary)^(2/5), where the
  primary's perturbation of motion about the body equals the body's of
  motion about the primary."
  ([planet] (sphere-of-influence (* c/AU (first (get-in planets/elements [planet :a])))
                                 (c/GM-planet planet) mu-sun))
  ([a mu mu-primary] (* a (math/pow (/ mu mu-primary) 0.4))))

(defn hill-radius
  "The Hill radius, a (mu/(3 mu-primary))^(1/3): the distance to the
  first Lagrange point, inside which a moon can keep an orbit."
  [a mu mu-primary]
  (* a (math/cbrt (/ mu (* 3.0 mu-primary)))))

;; ---------------------------------------------------------- hyperbolas

(defn departure
  "Leaving a circular parking orbit of radius `rp` about a body of
  parameter `mu` with hyperbolic excess speed `v-inf`: `{:dv :e
  :asymptote}`, the burn at periapsis (vis-viva on the hyperbola,
  sqrt(v-inf^2 + 2mu/rp), less the circular speed), the hyperbola's
  eccentricity 1 + rp v-inf^2/mu, and the true anomaly of its outgoing
  asymptote, acos(-1/e)."
  [mu rp v-inf]
  (let [e (+ 1.0 (/ (* rp v-inf v-inf) mu))]
    {:dv (- (math/sqrt (+ (* v-inf v-inf) (/ (* 2.0 mu) rp))) (math/sqrt (/ mu rp)))
     :e e
     :asymptote (math/acos (/ -1.0 e))}))

(defn capture
  "Arriving with excess speed `v-inf` and braking at periapsis `rp` into an
  orbit of apoapsis `ra` (circular by default): the burn."
  ([mu rp v-inf] (capture mu rp v-inf rp))
  ([mu rp v-inf ra]
   (let [a (* 0.5 (+ rp ra))]
     (- (math/sqrt (+ (* v-inf v-inf) (/ (* 2.0 mu) rp)))
        (math/sqrt (- (/ (* 2.0 mu) rp) (/ mu a)))))))

(defn turn-angle
  "How far a flyby of periapsis `rp` bends the excess velocity of speed
  `v-inf` about a body of parameter `mu`: 2 asin(1/e), e = 1 + rp v-inf^2/mu."
  [mu rp v-inf]
  (* 2.0 (math/asin (/ 1.0 (+ 1.0 (/ (* rp v-inf v-inf) mu))))))

(defn flyby
  "A gravity assist: the planet moving at `v-planet`, the spacecraft
  arriving at heliocentric velocity `v-in`, passing at periapsis `rp`
  about the planet (parameter `mu`) on the side that turns its excess
  velocity about the axis `n` (its part perpendicular to the excess
  velocity: the normal of the flyby's plane). Returns `{:v-out :v-inf :turn :dv}`,
  the heliocentric velocity after -- the excess velocity keeps its size
  and only turns, so the gain is 2 v-inf sin(turn/2), taken from the
  planet's motion."
  [mu v-planet v-in rp n]
  (let [vi (v3/sub v-in v-planet)
        v-inf (v3/length vi)
        d (turn-angle mu rp v-inf)
        ;; the plane of the flyby holds v-inf: only n's part across it counts
        u (v3/normalize vi)
        k (v3/normalize (v3/add-scaled n u (- (v3/dot n u))))
        ;; Rodrigues's rotation of vi by d about k
        vo (v3/add (v3/add (v3/scale vi (math/cos d)) (v3/scale (v3/cross k vi) (math/sin d)))
                   (v3/scale k (* (v3/dot k vi) (- 1.0 (math/cos d)))))
        v-out (v3/add v-planet vo)]
    {:v-out v-out :v-inf v-inf :turn d :dv (v3/distance v-out v-in)}))

;; ------------------------------------------------------ between planets

(defn hohmann
  "The Hohmann transfer between two planets' orbits taken as circles of
  their mean distances: `{:v-inf-depart :v-inf-arrive :tof}`, the excess
  speeds at each end and the half-orbit's time -- the least energetic
  way there, and the reference every real trajectory is compared with."
  [from to]
  (let [au #(* c/AU (first (get-in planets/elements [% :a])))
        r1 (au from) r2 (au to)
        at (* 0.5 (+ r1 r2))
        speed (fn [r a] (math/sqrt (* mu-sun (- (/ 2.0 r) (/ 1.0 a)))))]
    {:v-inf-depart (abs (- (speed r1 at) (speed r1 r1)))
     :v-inf-arrive (abs (- (speed r2 r2) (speed r2 at)))
     :tof (* math/PI (math/sqrt (/ (* at at at) mu-sun)))}))

(defn transfer
  "The heliocentric transfer from `from` at TT MJD `depart` to `to` at
  `arrive`, the short way unless `:long?`: Lambert's problem between the
  planets' places. Returns `{:v1 :v2 :v-inf-depart :v-inf-arrive :c3
  :dla}` -- the transfer's velocities, the excess velocities at each
  planet, the launch energy C3 = v-inf^2 and the declination of the
  departure asymptote, which a launch site's latitude must reach."
  ([from depart to arrive] (transfer from depart to arrive {}))
  ([from depart to arrive opts]
   (let [[r1 vp1] (planets/heliocentric-state from depart)
         [r2 vp2] (planets/heliocentric-state to arrive)
         [v1 v2] (iod/lambert mu-sun r1 r2 (* 86400.0 (- arrive depart)) opts)
         vinf1 (v3/sub v1 vp1)
         vinf2 (v3/sub v2 vp2)]
     {:v1 v1 :v2 v2 :v-inf-depart vinf1 :v-inf-arrive vinf2
      :c3 (v3/dot vinf1 vinf1)
      :dla (math/asin (/ (nth vinf1 2) (v3/length vinf1)))})))
