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
  (:require [allgo.astro.bplane :as bplane]
            [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.iod :as iod]
            [allgo.astro.planets :as planets]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
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

;; ------------------------------------------------- hyperbolas in space

(defn- asymptote-angle
  "The true anomaly of the asymptotes of the hyperbola of periapsis `rp`
  and excess speed `v-inf`, acos(-1/e), and its speed at periapsis."
  [mu rp v-inf]
  (let [e (+ 1.0 (/ (* rp v-inf v-inf) mu))]
    [(math/acos (/ -1.0 e)) (math/sqrt (+ (* v-inf v-inf) (/ (* 2.0 mu) rp)))]))

(defn escape-hyperbola
  "The state `[r v]` at periapsis, radius `rp`, of the hyperbola about a
  body of parameter `mu` whose outgoing excess velocity is the vector
  `v-inf`, in the plane of unit normal `h` (which must be normal to
  `v-inf`; the motion is counterclockwise about it). The outgoing
  asymptote S lies at the true anomaly nu = acos(-1/e) from periapsis, so
  periapsis is at cos nu S - sin nu (h x S)."
  [mu rp v-inf h]
  (let [s (v3/normalize v-inf)
        [nu vp] (asymptote-angle mu rp (v3/length v-inf))
        p (v3/sub (v3/scale s (math/cos nu)) (v3/scale (v3/cross h s) (math/sin nu)))]
    [(v3/scale p rp) (v3/scale (v3/cross h p) vp)]))

(defn approach-hyperbola
  "The state `[r v]` at periapsis, radius `rp`, of the hyperbola about a
  body of parameter `mu` arriving with the excess velocity `v-inf`, in the
  plane of unit normal `h` (normal to `v-inf`). The incoming asymptote
  comes from the true anomaly -nu, travelling along -cos nu P + sin nu Q,
  so periapsis is at -cos nu S - sin nu (h x S)."
  [mu rp v-inf h]
  (let [s (v3/normalize v-inf)
        [nu vp] (asymptote-angle mu rp (v3/length v-inf))
        p (v3/sub (v3/scale s (- (math/cos nu))) (v3/scale (v3/cross h s) (math/sin nu)))]
    [(v3/scale p rp) (v3/scale (v3/cross h p) vp)]))

(defn planes
  "The orbit normals of inclination `i` to the equator of pole `k` whose
  planes contain the direction `v-inf`: two, one when `i` just reaches the
  direction's declination, none when it falls short -- no orbit less
  inclined than an asymptote's declination can hold it. With S's
  declination delta, the normal's node lies where cos(Omega - alpha) =
  -cot i tan delta, alpha S's right ascension."
  [v-inf i k]
  (let [s (v3/normalize v-inf)
        sk (v3/dot s k)
        perp (v3/sub s (v3/scale k sk))
        pm (v3/length perp)]
    (if (< pm 1e-12)
      []
      (let [e1 (v3/scale perp (/ 1.0 pm))
            e2 (v3/cross k e1)
            si (math/sin i)
            c (if (< (abs si) 1e-15) ##Inf (/ (* -1.0 (math/cos i) sk) (* si pm)))]
        (cond
          (> (abs c) (+ 1.0 1e-12)) []
          :else
          (let [w (math/acos (max -1.0 (min 1.0 c)))
                normal #(v3/add (v3/scale k (math/cos i))
                                (v3/scale (v3/add (v3/scale e1 (math/cos %)) (v3/scale e2 (math/sin %))) si))]
            (if (< w 1e-9) [(normal 0.0)] [(normal w) (normal (- w))])))))))

(defn departure-orbits
  "Leaving a circular parking orbit of radius `rp` and inclination `i`
  (to the equator of pole `k`, default z -- the Earth's in EME2000) with
  the excess velocity vector `v-inf`: for each parking-orbit plane that
  holds the asymptote (see `planes`), `{:h :r :v :dv}` -- the plane's
  normal, the injection state at periapsis, and the burn there, which is
  along the circular velocity since both are normal to the radius."
  ([mu rp v-inf i] (departure-orbits mu rp v-inf i [0.0 0.0 1.0]))
  ([mu rp v-inf i k]
   (let [vc (math/sqrt (/ mu rp))]
     (vec (for [h (planes v-inf i k)
                :let [[r v] (escape-hyperbola mu rp v-inf h)]]
            {:h h :r r :v v :dv (- (v3/length v) vc)})))))

(defn arrival-orbits
  "Arriving with the excess velocity vector `v-inf` and braking at
  periapsis `rp` into an orbit of apoapsis `ra` and inclination `i` to the
  equator of pole `k` (for a planet, `allgo.astro.rotation/pole`): for
  each plane that holds the approach asymptote, `{:h :r :v :dv :bt :br}`
  -- the plane's normal, the state at periapsis on the approach
  hyperbola, the braking burn there, and the B-plane aim point (about `k`)
  that the approach must be steered to."
  [mu rp ra v-inf i k]
  (let [dv (capture mu rp (v3/length v-inf) ra)]
    (vec (for [h (planes v-inf i k)
               :let [[r v] (approach-hyperbola mu rp v-inf h)
                     {:keys [bt br]} (bplane/b-plane mu [r v] k)]]
           {:h h :r r :v v :dv dv :bt bt :br br}))))

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

(defn asymptote
  "Right ascension and declination `[ra dec]`, radians, of an excess
  velocity `v-inf` given in EME2000 -- or, with `mjd-tt`, referred to the
  mean equator and equinox of that date, as launch asymptotes are usually
  quoted (NASA's mission design handbooks among them)."
  ([v-inf]
   (let [[x y z] v-inf]
     [(am/wrap-2pi (math/atan2 y x)) (math/asin (/ z (v3/length v-inf)))]))
  ([v-inf mjd-tt] (asymptote (lin/mat-vec (frames/precession mjd-tt) v-inf))))

(defn transfer
  "The heliocentric transfer from `from` at TT MJD `depart` to `to` at
  `arrive`, the short way (type I) unless `:long?` (type II): Lambert's
  problem between the planets' places. Returns `{:v1 :v2 :v-inf-depart
  :v-inf-arrive :c3 :rla :dla}` -- the transfer's velocities, the excess
  velocities at each planet, the launch energy C3 = v-inf^2 and the right
  ascension and declination of the departure asymptote in EME2000 (the
  declination is what a launch site's latitude must reach).

  The planets' places come from `:ephemeris`, a function of planet and TT
  MJD giving heliocentric `[r v]` in EME2000; by default the mean elements
  of `allgo.astro.planets`, whose Earth is the Earth-Moon barycenter --
  which is what NASA's mission design handbooks use too, and against which
  this reproduces them to their printed precision.

  nil when no transfer fits: arrival not after departure, or too little
  time for the revolutions asked."
  ([from depart to arrive] (transfer from depart to arrive {}))
  ([from depart to arrive {:keys [ephemeris] :or {ephemeris planets/heliocentric-state} :as opts}]
   (when (> arrive depart)
     (let [[r1 vp1] (ephemeris from depart)
           [r2 vp2] (ephemeris to arrive)]
       (when-let [[v1 v2] (iod/lambert mu-sun r1 r2 (* 86400.0 (- arrive depart)) opts)]
         (let [vinf1 (v3/sub v1 vp1)
               vinf2 (v3/sub v2 vp2)
               [rla dla] (asymptote vinf1)]
           {:v1 v1 :v2 v2 :v-inf-depart vinf1 :v-inf-arrive vinf2
            :c3 (v3/dot vinf1 vinf1) :rla rla :dla dla}))))))
