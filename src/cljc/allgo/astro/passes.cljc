(ns allgo.astro.passes
  "Predicting a satellite's passes over a site (Vallado, *Fundamentals of
  Astrodynamics and Applications*, chapter 11): when it rises above the
  horizon (or a mask), when it culminates and how high, and when it sets.

  A pass is found where the elevation crosses the mask, the elevation
  sampled finely enough not to step over a pass and each crossing refined
  by bisection; the culmination is the maximum between, by golden-section
  search. Any elevation function of time will do -- `sgp4-elevation`
  makes one for a two-line element set."
  (:require [allgo.astro.geodesy :as geodesy]
            [allgo.astro.reduction :as reduction]
            [allgo.astro.sgp4 :as sgp4]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

(defn passes
  "The passes between `t0` and `t1` of a satellite whose elevation at time
  `t` is `(elevation t)`: `[{:rise :culmination :set :max-elevation}
  ...]`, :rise nil for a pass already under way at t0 and :set nil for one
  still going at t1. Options: `:mask` the elevation a pass must clear
  (default 0), `:step` the sampling (default 30, in t's units -- shorter
  than the briefest pass worth finding) and `:tol` how finely each time is
  found (default 1e-3)."
  ([elevation t0 t1] (passes elevation t0 t1 {}))
  ([elevation t0 t1 {:keys [mask step tol] :or {mask 0.0 step 30.0 tol 1e-3}}]
   (let [up? #(> (elevation %) mask)
         edges (roots/transitions up? t0 t1 step tol)
         bounds (loop [edges edges current (when (up? t0) {:rise nil}) out []]
                  (if-let [[[t was-up?] & more] (seq edges)]
                    (if was-up?
                      (recur more nil (conj out (assoc (or current {:rise nil}) :set t)))
                      (recur more {:rise t} out))
                    (if current (conj out (assoc current :set nil)) out)))]
     (mapv (fn [{:keys [rise set] :as p}]
             (let [c (roots/maximize elevation (or rise t0) (or set t1) {:tol tol})]
               (assoc p :culmination c :max-elevation (elevation c))))
           bounds))))

(defn sgp4-elevation
  "The elevation, radians, of the satellite of the element record `satrec`
  seen from geodetic `lat` `lon` `alt` (km), as a function of minutes from
  the elements' epoch: SGP4's TEME state taken Earth-fixed by the mean
  sidereal time (UT1 taken as UTC), then the look angles. -pi/2 where SGP4
  fails."
  [satrec lat lon alt]
  (let [site (geodesy/geodetic->cartesian lat lon alt)
        epoch (- (+ (:jdsatepoch satrec) (:jdsatepochF satrec)) 2400000.5)]
    (fn [minutes]
      (let [{:keys [r v error]} (sgp4/sgp4 satrec minutes)]
        (if (or (nil? r) (not (zero? error)))
          (- (/ math/PI 2))
          (let [[recef] (reduction/teme->ecef [r v] (+ epoch (/ minutes 1440.0)))]
            (:elevation (geodesy/look-angles site recef))))))))
