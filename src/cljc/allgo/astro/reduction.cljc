(ns allgo.astro.reduction
  "The FK5 reduction between inertial and Earth-fixed coordinates, for
  positions and velocities, through every intermediate frame (Vallado,
  *Fundamentals of Astrodynamics and Applications*, chapter 3).

    GCRF/J2000  --precession-->  MOD  --nutation-->  TOD
      --sidereal time-->  PEF  --polar motion-->  ITRF (Earth-fixed)

  and TEME, the frame SGP4 works in, hung off TOD: the true equator with
  the mean equinox. `allgo.astro.frames` has each rotation for the orbit
  integrator; this chains them with what the IERS publishes to make the
  reduction exact to its millimeters -- corrections to the nutation, the
  1997 terms of the equation of the equinoxes, and the length of day,
  which sets the Earth's spin rate and so the velocity it adds.

  Earth orientation parameters are a map: `:xp` `:yp` polar motion and
  `:ddpsi` `:ddeps` nutation corrections, radians; `:lod`, the excess
  length of day, seconds. Any left out are zero. Times are MJD, the Earth's
  rotation in UT1 and everything else in TT. States are `[r v]`, km and
  km/s."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.time :as time]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn- mv [m v] (lin/mat-vec m v))
(defn- tr [m] (lin/transpose m))

(def ^:private omega-nominal
  "The Earth's mean rotation rate, rad/s, that the length of day adjusts."
  7.292115146706979e-5)

(defn nutation-angles
  "`[dpsi deps eps]`: the IAU 1980 nutation with the IERS corrections, and
  the true obliquity."
  [mjd-tt {:keys [ddpsi ddeps] :or {ddpsi 0.0 ddeps 0.0}}]
  (let [[dpsi deps] (frames/nutation-angles mjd-tt)
        dpsi (+ dpsi ddpsi) deps (+ deps ddeps)]
    [dpsi deps (+ (frames/mean-obliquity mjd-tt) deps)]))

(defn- eqe-1982
  "The equation of the equinoxes as TEME is defined with it: dpsi times
  the cosine of the mean obliquity, as the IAU's 1994 form has it."
  [mjd-tt eop]
  (* (first (nutation-angles mjd-tt eop)) (math/cos (frames/mean-obliquity mjd-tt))))

(defn equation-of-equinoxes
  "Apparent less mean sidereal time, the IAU's 1994 form with the two terms
  added in 1997 for the kinematic effect of the node's regression: dpsi cos
  eps_A + 0.00264'' sin Omega + 0.000063'' sin 2 Omega, eps_A the mean
  obliquity."
  [mjd-tt eop]
  (let [om (nth (frames/delaunay mjd-tt) 4)]
    (+ (eqe-1982 mjd-tt eop)
       (* c/arcsec (+ (* 0.00264 (math/sin om)) (* 0.000063 (math/sin (* 2.0 om))))))))

(defn- nutation-matrix [mjd-tt eop]
  (let [[dpsi _ eps] (nutation-angles mjd-tt eop)
        eps0 (frames/mean-obliquity mjd-tt)]
    (frames/chain (frames/rx (- eps)) (frames/rz (- dpsi)) (frames/rx eps0))))

(defn- polar-matrix
  "PEF to ITRF, Vallado's full rotation rather than its small-angle form."
  [{:keys [xp yp] :or {xp 0.0 yp 0.0}}]
  (let [cx (math/cos xp) sx (math/sin xp) cy (math/cos yp) sy (math/sin yp)]
    (tr [[cx 0.0 (- sx)]
         [(* sx sy) cy (* cx sy)]
         [(* sx cy) (- sy) (* cx cy)]])))

(defn- omega [{:keys [lod] :or {lod 0.0}}]
  [0.0 0.0 (* omega-nominal (- 1.0 (/ lod 86400.0)))])

;; ------------------------------------------------ inertial to rotating

(defn eci->mod
  "GCRF (J2000) to mean-of-date."
  [[r v] mjd-tt]
  (let [P (frames/precession mjd-tt)] [(mv P r) (mv P v)]))

(defn eci->tod
  "GCRF to true-of-date."
  ([s mjd-tt] (eci->tod s mjd-tt {}))
  ([s mjd-tt eop]
   (let [N (nutation-matrix mjd-tt eop) [r v] (eci->mod s mjd-tt)] [(mv N r) (mv N v)])))

(defn eci->teme
  "GCRF to TEME, the true equator and mean equinox SGP4's elements are
  referred to: true of date turned along the equator by the equation of
  the equinoxes in its 1982 form -- without the 1997 terms, which TEME's
  definition predates."
  ([s mjd-tt] (eci->teme s mjd-tt {}))
  ([s mjd-tt eop]
   (let [[r v] (eci->tod s mjd-tt eop)
         R (frames/rz (eqe-1982 mjd-tt eop))]
     [(mv R r) (mv R v)])))

(defn eci->pef
  "GCRF to the pseudo-Earth-fixed frame: turned with the Earth by the
  apparent sidereal time, but not yet for polar motion. The velocity
  loses the Earth's spin."
  ([s mjd-tt mjd-ut1] (eci->pef s mjd-tt mjd-ut1 {}))
  ([s mjd-tt mjd-ut1 eop]
   (let [[r v] (eci->tod s mjd-tt eop)
         R (frames/rz (+ (time/gmst mjd-ut1) (equation-of-equinoxes mjd-tt eop)))
         rp (mv R r)]
     [rp (lin/sub (mv R v) (v3/cross (omega eop) rp))])))

(defn eci->ecef
  "GCRF to the Earth-fixed ITRF."
  ([s mjd-tt mjd-ut1] (eci->ecef s mjd-tt mjd-ut1 {}))
  ([s mjd-tt mjd-ut1 eop]
   (let [[r v] (eci->pef s mjd-tt mjd-ut1 eop)
         W (polar-matrix eop)]
     [(mv W r) (mv W v)])))

;; ------------------------------------------------ rotating to inertial

(defn mod->eci [[r v] mjd-tt]
  (let [P (tr (frames/precession mjd-tt))] [(mv P r) (mv P v)]))

(defn tod->eci
  ([s mjd-tt] (tod->eci s mjd-tt {}))
  ([[r v] mjd-tt eop]
   (let [N (tr (nutation-matrix mjd-tt eop))] (mod->eci [(mv N r) (mv N v)] mjd-tt))))

(defn teme->eci
  "TEME to GCRF: how an SGP4 state is brought into the inertial frame."
  ([s mjd-tt] (teme->eci s mjd-tt {}))
  ([[r v] mjd-tt eop]
   (let [R (frames/rz (- (eqe-1982 mjd-tt eop)))]
     (tod->eci [(mv R r) (mv R v)] mjd-tt eop))))

(defn pef->eci
  ([s mjd-tt mjd-ut1] (pef->eci s mjd-tt mjd-ut1 {}))
  ([[r v] mjd-tt mjd-ut1 eop]
   (let [R (tr (frames/rz (+ (time/gmst mjd-ut1) (equation-of-equinoxes mjd-tt eop))))
         v' (lin/add v (v3/cross (omega eop) r))]
     (tod->eci [(mv R r) (mv R v')] mjd-tt eop))))

(defn ecef->eci
  "ITRF to GCRF: the inverse of `eci->ecef`."
  ([s mjd-tt mjd-ut1] (ecef->eci s mjd-tt mjd-ut1 {}))
  ([[r v] mjd-tt mjd-ut1 eop]
   (let [W (tr (polar-matrix eop))]
     (pef->eci [(mv W r) (mv W v)] mjd-tt mjd-ut1 eop))))

(defn teme->ecef
  "TEME to ITRF directly, as a ground track from SGP4 needs: turned by the
  mean sidereal time, the frame being on the mean equinox."
  ([s mjd-ut1] (teme->ecef s mjd-ut1 {}))
  ([[r v] mjd-ut1 eop]
   (let [R (frames/rz (time/gmst mjd-ut1))
         rp (mv R r)
         vp (lin/sub (mv R v) (v3/cross (omega eop) rp))
         W (polar-matrix eop)]
     [(mv W rp) (mv W vp)])))
