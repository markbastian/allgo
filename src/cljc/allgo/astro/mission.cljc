(ns allgo.astro.mission
  "Mission geometry (Vallado, *Fundamentals of Astrodynamics and
  Applications*, chapter 11): how far apart two places are and in which
  direction, where a satellite is over the ground, the orbits that retrace
  their tracks or hold their shape, and what a sensor can see.

  The Earth here is a sphere of its equatorial radius unless said
  otherwise -- `allgo.astro.geodesy` has the ellipsoid -- and angles are
  radians, longitudes positive east, azimuths from north through east."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as geodesy]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.perturbations :as perturbations]
            [allgo.astro.reduction :as reduction]
            [allgo.geometry.sphere :as sphere]
            [allgo.math :as am]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

;; ----------------------------------------------------- between two places

(defn range-azimuth
  "`{:range :angle :azimuth}` from `[lat1 lon1]` to `[lat2 lon2]` on a
  sphere of radius `R` (default the Earth's): the great-circle distance,
  the angle it subtends at the center, and the initial azimuth."
  ([p1 p2] (range-azimuth p1 p2 c/R-earth))
  ([p1 p2 R]
   (let [{:keys [angle azimuth]} (sphere/great-circle p1 p2)]
     {:range (* R angle) :angle angle :azimuth azimuth})))

(defn destination
  "`[lat lon]` reached from `[lat lon]` by going `range` along a great
  circle that sets out at `azimuth`, on a sphere of radius `R`."
  ([p range azimuth] (destination p range azimuth c/R-earth))
  ([p range azimuth R] (sphere/destination p (/ range R) azimuth)))

;; ---------------------------------------------------------- ground track

(defn ground-point
  "`[lat lon height]`, geodetic, of the point beneath the inertial state
  `s` at TT `mjd-tt` and UT1 `mjd-ut1`: the state taken Earth-fixed by the
  FK5 reduction, then onto the ellipsoid."
  ([s mjd-tt mjd-ut1] (ground-point s mjd-tt mjd-ut1 {}))
  ([s mjd-tt mjd-ut1 eop]
   (geodesy/cartesian->geodetic (first (reduction/eci->ecef s mjd-tt mjd-ut1 eop)))))

(defn nodal-period
  "Seconds from one ascending node to the next, with J2 turning the perigee
  and speeding the mean anomaly: 2 pi / (domega/dt + dM/dt)."
  [a e i]
  (let [{:keys [argp M]} (perturbations/j2-secular a e i)]
    (/ (* 2.0 math/PI) (+ argp M))))

(defn repeat-period
  "Seconds for `revs` nodal periods, and the days of Earth rotation they
  span, relative to a node J2 is turning: `{:revs-time :days}`. An orbit
  repeats its ground track when :days is a whole number."
  [a e i revs]
  (let [{:keys [raan]} (perturbations/j2-secular a e i)
        t (* revs (nodal-period a e i))]
    {:revs-time t :days (/ (* t (- c/omega-earth raan)) (* 2.0 math/PI))}))

(defn orbit-for-period
  "The semi-major axis of a two-body orbit of period `T` seconds."
  [T]
  (math/cbrt (* c/GM-earth (am/sq (/ T (* 2.0 math/PI))))))

(defn repeat-ground-track
  "The semi-major axis at which an orbit of eccentricity `e` and
  inclination `i` makes `revs` revolutions in `days` turns of the Earth
  under its node, retracing its ground track: found by bisection on the
  node-to-node timing, with J2 moving the node and speeding the orbit."
  [revs days e i]
  (let [miss (fn [a] (- (:days (repeat-period a e i revs)) days))
        a-guess (orbit-for-period (/ (* days 86164.0905) revs))]
    (roots/bisect miss (* 0.8 a-guess) (* 1.2 a-guess))))

(defn ground-track-drift
  "How fast the ground track of an orbit at `a` + da slides against that
  of the orbit at `a` (eccentricity `e`, inclination `i`), per km of da:
  rad/s of longitude at the node, east positive. Each revolution the node
  moves (dOmega/dt - omega_earth) P in longitude; a higher orbit takes
  longer over it, the Earth turning further beneath, and its node
  regresses more slowly -- the derivative of that per-revolution shift,
  divided by the period, taken by central difference."
  [a e i]
  (let [shift (fn [a] (* (- (:raan (perturbations/j2-secular a e i)) c/omega-earth) (nodal-period a e i)))
        h 1e-3]
    (/ (/ (- (shift (+ a h)) (shift (- a h))) (* 2.0 h)) (nodal-period a e i))))

(defn ground-track-maintenance
  "Keeping a repeat ground track within `tolerance` km east or west at the
  equator while drag lowers the orbit at `a-rate` km/s (negative; see
  `perturbations/circular-decay-rate`), for the reference orbit `a` `e`
  `i` (Vallado, algorithm 71; Wertz, Mission Geometry, 2001).

  An orbit da above the reference slides west at K da (see
  `ground-track-drift`); decaying, da falls through zero and the slide
  turns round, so the track's error traces a parabola in time. Started at
  the east edge da0 above the reference, it just grazes the west edge
  when da0 = sqrt(4 tolerance a-rate / (K R)) and comes back to the east
  edge after 2 da0 / |a-rate|, da0 below -- when a burn lifts the orbit
  2 da0 again. Returns `{:drift :da :cycle :dv :dv-per-year}`: K, da0 in
  km, the cycle in s, each burn's speed km/s (v da/a, the two-impulse
  raise of a near-circular orbit taken as one) and a year's worth."
  [a e i a-rate tolerance]
  (let [k (ground-track-drift a e i)
        band (/ tolerance c/R-earth)
        da (math/sqrt (/ (* 4.0 band a-rate) k))
        cycle (/ (* 2.0 da) (abs a-rate))
        dv (* (math/sqrt (/ c/GM-earth a)) (/ da a))]
    {:drift k :da da :cycle cycle :dv dv :dv-per-year (* dv (/ (* 365.25 86400.0) cycle))}))

;; --------------------------------------------------------------- frozen

(def J3
  "The Earth's J3, unnormalized: the pear shape, the southern hemisphere a
  little fuller."
  -2.5324e-6)

(defn frozen-eccentricity
  "The eccentricity at which, with the perigee at 90 degrees, J3's pull on
  the perigee cancels J2's and the orbit's shape holds: -(J3/J2)(R/a)
  sin i / 2, some 0.001 for a low orbit. (At 270 degrees the sign would
  flip, and no eccentricity is negative.)"
  [a i]
  (* -0.5 (/ J3 geo/J2) (/ c/R-earth a) (math/sin i)))

;; ---------------------------------------------------------- field of view

(defn field-of-view
  "What a sensor at height `h` looking `eta` off nadir sees on a spherical
  Earth: `{:elevation :central-angle :slant-range :swath}` -- the
  elevation of the satellite seen from the point viewed, the Earth
  central angle from the sub-satellite point to it, the distance to it,
  and the width of ground a cone of half-angle eta covers. nil beyond the
  horizon, where sin eta > R/(R + h)."
  [h eta]
  (let [R c/R-earth
        rs (+ R h)
        ce (/ (* rs (math/sin eta)) R)]
    (when (<= ce 1.0)
      (let [elevation (math/acos ce)
            lam (- (* 0.5 math/PI) eta elevation)]
        {:elevation elevation
         :central-angle lam
         :slant-range (/ (* R (math/sin lam)) (math/sin eta))
         :swath (* 2.0 R lam)}))))

(defn horizon
  "The nadir angle and Earth central angle of the horizon from height `h`."
  [h]
  (let [s (/ c/R-earth (+ c/R-earth h))]
    {:nadir (math/asin s) :central-angle (math/acos s)}))
