(ns allgo.astro.forces
  "The complete force model of Montenbruck & Gill chapter 3, assembled into
  an acceleration the integrators of chapter 4 can take.

  A caveat that decides which integrator is usable. Nystrom and
  Stoermer-Cowell methods integrate y'' = f(t, y) -- acceleration as a
  function of position alone. Gravity, third-body and radiation pressure are
  of that form; drag and the relativistic correction are not, since both
  depend on velocity. Switch either on and the system leaves the special
  form, and only the general first-order integrators still apply. That is
  not a limitation of this code but of the method, and `second-order`
  refuses rather than silently dropping the terms it cannot represent."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.drag :as drag]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.frames :as frames]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.relativity :as rel]
            [allgo.astro.srp :as srp]
            [allgo.astro.tides :as tides]
            [allgo.astro.time :as time]
            [clojure.math :as math]))

(defn earth-fixed
  "The rotation from the celestial frame of J2000 to the Earth-fixed frame
  at `mjd-tt`, as a matrix.

  The full chain of chapter 5, not sidereal time alone. Rotating by GMST
  only would ignore that the pole and equinox have moved since J2000, which
  is nil at the epoch and tens of kilometres at the surface two decades on.

  UT1 is taken as UTC unless `dut1` is given: the difference is under a
  second by construction, and cannot be computed, only looked up."
  ([mjd-tt] (earth-fixed mjd-tt 0.0))
  ([mjd-tt dut1]
   (frames/celestial->terrestrial-cached
    mjd-tt (time/utc->ut1 (time/tt->utc mjd-tt) dut1))))

(defn eci->ecef [r m] (frames/apply-m m r))
(defn ecef->eci [r m] (frames/apply-m (frames/transpose m) r))

(def defaults
  {:degree       4
   :sun?         true
   :moon?        true
   :tides?       false
   :relativity?  false
   :srp          nil        ; {:area-to-mass m^2/kg :cr reflectivity}
   :drag         nil        ; {:area-to-mass m^2/kg :cd coefficient}
   :field        nil})      ; defaults to geopotential/earth

(defn velocity-dependent?
  "Whether the configured model needs velocity, and so cannot be integrated
  by a Nystrom or Stoermer-Cowell method."
  [config]
  (let [{:keys [drag relativity?]} (merge defaults config)]
    (boolean (or drag relativity?))))

(defn acceleration
  "Total inertial acceleration on a spacecraft, km/s^2.

  `mjd` is the epoch as a Modified Julian Date, `r` and `v` the inertial
  position and velocity. Gravity is evaluated in the Earth-fixed frame and
  rotated back, since that is the frame the harmonic coefficients are tied
  to -- the field turns with the planet."
  [config mjd r v]
  (let [{:keys [degree sun? moon? tides? relativity? field]
         srp-cfg :srp drag-cfg :drag} (merge defaults config)
        u     (earth-fixed mjd)
        r-sun (eph/sun mjd)
        r-moon (eph/moon mjd)
        base  (or field geo/earth)
        ;; The tidal bulge is fixed to the Earth, so the bodies raising it
        ;; must be given in Earth-fixed coordinates like the field itself.
        fld   (if tides?
                (tides/perturb base (cond-> []
                                      moon? (conj [c/GM-moon (eci->ecef r-moon u)])
                                      sun?  (conj [c/GM-sun (eci->ecef r-sun u)])))
                base)
        grav  (ecef->eci (geo/acceleration fld (eci->ecef r u) degree) u)]
    (cond-> grav
      sun?     (as-> a (mapv + a (eph/third-body c/GM-sun r r-sun)))
      moon?    (as-> a (mapv + a (eph/third-body c/GM-moon r r-moon)))
      srp-cfg  (as-> a (mapv + a (srp/acceleration r r-sun
                                                   (:area-to-mass srp-cfg)
                                                   (:cr srp-cfg))))
      drag-cfg (as-> a (mapv + a (drag/acceleration r v r-sun
                                                    (:area-to-mass drag-cfg)
                                                    (:cd drag-cfg))))
      relativity? (as-> a (mapv + a (rel/acceleration r v))))))

(defn breakdown
  "Each contribution separately, for inspection: a map from force to its
  acceleration vector. Useful for seeing which terms actually matter at a
  given altitude, which varies enormously between low orbit and
  geostationary."
  [config mjd r v]
  (let [{:keys [degree sun? moon? relativity? field] srp-cfg :srp drag-cfg :drag}
        (merge defaults config)
        u      (earth-fixed mjd)
        r-sun  (eph/sun mjd)
        r-moon (eph/moon mjd)
        base   (or field geo/earth)
        ecef   (eci->ecef r u)
        pm     (geo/point-mass c/GM-earth c/R-earth)]
    (cond-> {:two-body   (ecef->eci (geo/acceleration pm ecef 0) u)
             :harmonics  (ecef->eci (mapv - (geo/acceleration base ecef degree)
                                          (geo/acceleration pm ecef 0)) u)}
      sun?        (assoc :sun (eph/third-body c/GM-sun r r-sun))
      moon?       (assoc :moon (eph/third-body c/GM-moon r r-moon))
      srp-cfg     (assoc :srp (srp/acceleration r r-sun (:area-to-mass srp-cfg) (:cr srp-cfg)))
      drag-cfg    (assoc :drag (drag/acceleration r v r-sun (:area-to-mass drag-cfg) (:cd drag-cfg)))
      relativity? (assoc :relativity (rel/acceleration r v)))))

;; ------------------------------------------- adapters for the integrators

(defn first-order
  "`(fn [t y])` for the general integrators, where `y` is the six-vector of
  position and velocity and `t` is seconds from `epoch-mjd`. Every force is
  available in this form."
  [config epoch-mjd]
  (fn [t y]
    (let [r [(nth y 0) (nth y 1) (nth y 2)]
          v [(nth y 3) (nth y 4) (nth y 5)]
          a (acceleration config (+ epoch-mjd (/ t 86400.0)) r v)]
      (into v a))))

(defn second-order
  "`(fn [t r])` for the Nystrom and Stoermer-Cowell integrators.

  Throws if the model includes a velocity-dependent force, because there is
  no honest way to supply one: the method's whole efficiency comes from
  assuming the acceleration does not need a velocity, and quietly dropping
  drag to satisfy that would give a plausible-looking wrong answer."
  [config epoch-mjd]
  (when (velocity-dependent? config)
    (throw (ex-info "Drag and the relativistic correction depend on velocity, so this model is not of the form y'' = f(t, y) that Nystrom and Stoermer-Cowell methods require. Integrate it with a general first-order method instead."
                    {:config config})))
  (fn [t r]
    (acceleration config (+ epoch-mjd (/ t 86400.0)) r nil)))

(defn circular-state
  "Position and velocity for a circular orbit of the given radius and
  inclination, starting at the ascending node."
  [radius inclination]
  (let [v (math/sqrt (/ c/GM-earth radius))]
    [[radius 0.0 0.0]
     [0.0 (* v (math/cos inclination)) (* v (math/sin inclination))]]))
