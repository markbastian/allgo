(ns allgo.astro.states
  "The other ways to write down a state vector, and the local frames it is
  seen in (Vallado, *Fundamentals of Astrodynamics and Applications*,
  chapters 2 to 4).

  Position and velocity are six numbers, and there are many useful sixes:
  right ascension, declination and their rates, as a telescope reports
  them; range, azimuth, elevation and their rates, as a radar does; the
  flight elements a launch vehicle's guidance speaks in; and equinoctial
  elements, which have no singularity at a circular or equatorial orbit.
  Each is here with its inverse. So are the frames a satellite's own
  motion defines -- RSW, radial-along-cross, for orbit errors; NTW, along
  the velocity, for burns; PQW, the perifocal frame of the orbit -- and
  SEZ, the south-east-zenith frame of a site on the ground.

  Kilometers, seconds, radians; azimuth from north through east."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as geodesy]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.reduction :as reduction]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn- mv [m v] (lin/mat-vec m v))

;; ------------------------------------------------------- right ascension

(defn state->radec
  "`[range ra dec range-rate ra-rate dec-rate]` of a state, about the
  origin of its frame."
  [[r v]]
  (let [[x y z] r [vx vy vz] v
        rr (v3/length r)
        h (math/hypot x y)
        ra (if (< h 1e-12) (math/atan2 vy vx) (math/atan2 y x))
        dec (math/asin (/ z rr))
        drr (/ (v3/dot r v) rr)]
    [rr ra dec drr
     (/ (- (* vx y) (* vy x)) (- (* h h)))
     (/ (- vz (* drr (/ z rr))) h)]))

(defn radec->state
  "The state from `[range ra dec range-rate ra-rate dec-rate]`."
  [[rr ra dec drr dra ddec]]
  (let [cd (math/cos dec) sd (math/sin dec) ca (math/cos ra) sa (math/sin ra)]
    [[(* rr cd ca) (* rr cd sa) (* rr sd)]
     [(- (* drr cd ca) (* rr sd ca ddec) (* rr cd sa dra))
      (+ (- (* drr cd sa) (* rr sd sa ddec)) (* rr cd ca dra))
      (+ (* drr sd) (* rr cd ddec))]]))

(def ^:private obliquity
  "Vallado's constant obliquity, 23.439291 degrees, for the ecliptic forms."
  (* 23.439291 c/degrees))

(defn- ecl-matrix [sign]
  (let [ce (math/cos obliquity) se (math/sin obliquity)]
    [[1.0 0.0 0.0] [0.0 ce (* sign se)] [0.0 (- (* sign se)) ce]]))

(defn state->ecliptic
  "`[range lon lat range-rate lon-rate lat-rate]` in ecliptic coordinates."
  [[r v]]
  (let [m (ecl-matrix 1.0)] (state->radec [(mv m r) (mv m v)])))

(defn ecliptic->state [e]
  (let [[r v] (radec->state e) m (ecl-matrix -1.0)] [(mv m r) (mv m v)]))

(defn state->topocentric-radec
  "`[range ra dec ...rates]` of a state seen from a site whose state is
  `site` (same inertial frame): the direction an optical tracker reports."
  [[r v] [rs vs]]
  (state->radec [(v3/sub r rs) (v3/sub v vs)]))

(defn topocentric-radec->state [tr [rs vs]]
  (let [[r v] (radec->state tr)] [(v3/add r rs) (v3/add v vs)]))

(defn state->adbar
  "`[r v ra dec fpa-vertical azimuth]`: the spherical elements a launch is
  planned in -- position by its right ascension and declination, velocity
  by its angle from the local vertical and its azimuth."
  [[r v]]
  (let [rm (v3/length r) vm (v3/length v)
        [_ ra dec] (state->radec [r v])
        fpav (math/atan2 (v3/length (v3/cross r v)) (v3/dot r v))]
    [rm vm ra dec fpav
     (let [;; the east and north unit vectors at the position
           e (v3/normalize [(- (math/sin ra)) (math/cos ra) 0.0])
           n (v3/normalize (v3/cross (v3/normalize r) e))]
       (am/wrap-2pi (math/atan2 (v3/dot v e) (v3/dot v n))))]))

(defn adbar->state [[rm vm ra dec fpav az]]
  (let [sa (math/sin ra) ca (math/cos ra) sd (math/sin dec) cd (math/cos dec)
        sf (math/sin fpav) cf (math/cos fpav) sz (math/sin az) cz (math/cos az)]
    [[(* rm cd ca) (* rm cd sa) (* rm sd)]
     [(* vm (- (* ca (+ (* (- cz) sf sd) (* cf cd))) (* sz sf sa)))
      (* vm (+ (* sa (+ (* (- cz) sf sd) (* cf cd))) (* sz sf ca)))
      (* vm (+ (* cz cd sf) (* cf sd)))]]))

;; ---------------------------------------------------- classical elements

(defn state->classical
  "Vallado's RV2COE: the classical elements with the extra ones that stand
  in where an orbit makes some undefined -- `:p` the semilatus rectum,
  `:arglat` the argument of latitude (for circular inclined orbits),
  `:truelon` the true longitude (circular equatorial) and `:lonper` the
  longitude of periapsis (elliptical equatorial), NaN where not needed --
  and `:type`."
  ([s] (state->classical c/GM-earth s))
  ([mu [r v]]
   (let [{:keys [e i argp nu] :as el} (kepler/state->elements mu r v)
         p (/ (v3/length-squared (v3/cross r v)) mu)
         circular? (< e 1e-8)
         equatorial? (or (< i 1e-8) (< (abs (- i math/PI)) 1e-8))
         type (cond (and circular? equatorial?) :circular-equatorial
                    circular? :circular-inclined
                    equatorial? :elliptical-equatorial
                    :else :elliptical-inclined)]
     (assoc el :p p :type type
            :arglat (if equatorial? ##NaN (am/wrap-2pi (+ argp nu)))
            :truelon (if (= type :circular-equatorial) nu ##NaN)
            :lonper (if (= type :elliptical-equatorial) argp ##NaN)
            :M (let [[_ M] (kepler/true->anomaly-and-mean e nu)] M)))))

(defn classical->state
  "The state from `{:p :e :i :raan :argp :nu}` (Vallado's COE2RV)."
  ([el] (classical->state c/GM-earth el))
  ([mu {:keys [p e] :as el}]
   (kepler/elements->state mu (assoc el :a (/ p (- 1.0 (* e e)))))))

(defn perifocal-axes
  "`[P Q W]`, the perifocal frame's axes in the inertial frame: toward
  periapsis, 90 degrees on in the orbit, and along the orbit normal."
  [i raan argp]
  (let [co (math/cos raan) so (math/sin raan) ci (math/cos i) si (math/sin i)
        cw (math/cos argp) sw (math/sin argp)]
    [[(- (* co cw) (* so sw ci)) (+ (* so cw) (* co sw ci)) (* sw si)]
     [(- (- (* co sw)) (* so cw ci)) (+ (- (* so sw)) (* co cw ci)) (* cw si)]
     [(* so si) (- (* co si)) ci]]))

;; --------------------------------------------------- equinoctial elements

(defn state->equinoctial
  "Vallado's equinoctial elements `{:a :n :af :ag :chi :psi :meanlon
  :truelon :fr}`: the eccentricity vector's components af, ag and the
  node vector's chi, psi, measured from the equinoctial frame, with `fr`
  +1 for prograde orbits and -1 for retrograde, which keeps the one
  remaining singularity at i = 180 degrees (or 0) out of the way."
  ([s] (state->equinoctial c/GM-earth s))
  ([mu [r v]]
   (let [{:keys [a e i raan argp nu M]} (state->classical mu [r v])
         fr (if (> i (/ math/PI 2)) -1.0 1.0)
         t (math/pow (math/tan (* 0.5 i)) fr)]
     {:a a :n (math/sqrt (/ mu (* a a a)))
      :af (* e (math/cos (+ argp (* fr raan)))) :ag (* e (math/sin (+ argp (* fr raan))))
      :chi (* t (math/sin raan)) :psi (* t (math/cos raan))
      :meanlon (am/wrap-2pi (+ M argp (* fr raan)))
      :truelon (am/wrap-2pi (+ nu argp (* fr raan))) :fr fr})))

(defn equinoctial->state
  "The state from `{:a :af :ag :chi :psi :meanlon :fr}`."
  ([el] (equinoctial->state c/GM-earth el))
  ([mu {:keys [a af ag chi psi meanlon fr] :or {fr 1.0}}]
   (let [e (math/hypot af ag)
         raan (math/atan2 chi psi)
         i (* 2.0 (math/atan (math/pow (math/hypot chi psi) fr)))
         lonper (math/atan2 ag af)
         argp (- lonper (* fr raan))
         M (- meanlon lonper)]
     (kepler/elements->state mu {:a a :e e :i i :raan raan :argp argp :M (am/wrap-2pi M)}))))

;; ------------------------------------------------------- satellite frames

(defn rsw
  "The radial, along-track, cross-track frame of a state: rows R, S, W."
  [[r v]]
  (let [R (v3/normalize r) W (v3/normalize (v3/cross r v)) S (v3/cross W R)] [R S W]))

(defn ntw
  "The frame along the velocity: N in-plane normal to it, T along it, W
  the orbit normal."
  [[r v]]
  (let [T (v3/normalize v) W (v3/normalize (v3/cross r v)) N (v3/cross T W)] [N T W]))

(defn into-frame
  "A state's components in the frame whose rows are `axes`."
  [axes [r v]]
  [(mv axes r) (mv axes v)])

(defn state->pqw
  "A state in its own perifocal frame."
  ([s] (state->pqw c/GM-earth s))
  ([mu [r v]]
   (let [{:keys [i raan argp]} (state->classical mu [r v])]
     (into-frame (perifocal-axes i raan argp) [r v]))))

(defn eci->hill
  "The state of `interceptor` relative to `target`, both inertial `[r v]`,
  in the target's rotating Hill frame -- x radial, y along-track, z
  cross-track, the RSW axes turning with the target (Vallado's ECI2HILL):
  rho = R (r_i - r_t), and rho' = R (v_i - v_t) - omega x rho, omega =
  h/r^2 about the orbit normal the frame's rate. The linear
  Clohessy-Wiltshire solution (`maneuvers/hill`) lives in this frame."
  [target interceptor]
  (let [[rt vt] target [ri vi] interceptor
        axes (rsw target)
        rho (mv axes (v3/sub ri rt))
        w (/ (v3/length (v3/cross rt vt)) (v3/dot rt rt))]
    [rho (v3/sub (mv axes (v3/sub vi vt)) (v3/cross [0.0 0.0 w] rho))]))

(defn hill->eci
  "The inertial state of an interceptor at `rel`, `[rho rho']` in the
  Hill frame of `target` (Vallado's HILL2ECI): the inverse of `eci->hill`."
  [target [rho drho]]
  (let [[rt vt] target
        axes (rsw target)
        back (lin/transpose axes)
        w (/ (v3/length (v3/cross rt vt)) (v3/dot rt rt))]
    [(v3/add rt (mv back rho))
     (v3/add vt (mv back (v3/add drho (v3/cross [0.0 0.0 w] rho))))]))

;; -------------------------------------------------- site, SEZ and radar

(defn sez-axes
  "Rows: south, east and zenith at geodetic latitude `lat`, longitude `lon`."
  [lat lon]
  (let [sl (math/sin lat) cl (math/cos lat) so (math/sin lon) co (math/cos lon)]
    [[(* sl co) (* sl so) (- cl)] [(- so) co 0.0] [(* cl co) (* cl so) sl]]))

(defn razel->sez
  "Range, azimuth, elevation and their rates as a relative state in SEZ."
  [[rho az el drho daz del]]
  (let [ce (math/cos el) se (math/sin el) ca (math/cos az) sa (math/sin az)]
    [[(* (- rho) ce ca) (* rho ce sa) (* rho se)]
     [(+ (* (- drho) ce ca) (* rho se ca del) (* rho ce sa daz))
      (+ (* drho ce sa) (* (- rho) se sa del) (* rho ce ca daz))
      (+ (* drho se) (* rho ce del))]]))

(defn sez->razel
  "`[range az el range-rate az-rate el-rate]` of a relative state in SEZ."
  [[r v]]
  (let [[s e z] r [ds de dz] v
        rho (v3/length r)
        h (math/hypot s e)
        az (am/wrap-angle (math/atan2 e (- s)))
        el (math/asin (/ z rho))
        drho (/ (v3/dot r v) rho)]
    [rho az el drho
     (if (< h 1e-12) 0.0 (/ (- (* ds e) (* de s)) (* h h)))
     (/ (- dz (* drho (math/sin el))) h)]))

(defn razel->ecef
  "The Earth-fixed state of a target seen at range, azimuth, elevation and
  their rates from a site at geodetic `lat` `lon` and height `alt` km."
  [razel lat lon alt]
  (let [[rs vs] (razel->sez razel)
        m (lin/transpose (sez-axes lat lon))]
    [(v3/add (geodesy/geodetic->cartesian lat lon alt) (mv m rs)) (mv m vs)]))

(defn site-track
  "The inertial state of a target a site sees at `razel` -- range,
  azimuth, elevation and their rates -- from geodetic `lat` `lon` and
  height `alt` km, at `mjd-tt` `mjd-ut1` with Earth orientation `eop`
  (Vallado's SITE-TRACK, algorithm 50): the site-relative state in SEZ,
  turned to the Earth-fixed frame and added to the site's position, then
  reduced to the GCRF, the Earth's turning adding its omega x r."
  ([razel lat lon alt mjd-tt mjd-ut1] (site-track razel lat lon alt mjd-tt mjd-ut1 {}))
  ([razel lat lon alt mjd-tt mjd-ut1 eop]
   (reduction/ecef->eci (razel->ecef razel lat lon alt) mjd-tt mjd-ut1 eop)))

(defn ecef->razel
  "Range, azimuth, elevation and rates of an Earth-fixed state from a site."
  [[r v] lat lon alt]
  (let [m (sez-axes lat lon)
        rho (v3/sub r (geodesy/geodetic->cartesian lat lon alt))]
    (sez->razel [(mv m rho) (mv m v)])))

;; ------------------------------------------------------ flight elements

(defn flight->state
  "The inertial state from Vallado's flight elements -- magnitudes `rm`
  `vm`, geocentric latitude and longitude, flight-path angle from the
  horizontal and azimuth, all Earth-fixed -- at `mjd-tt` `mjd-ut1` with
  Earth orientation `eop`."
  [[rm vm latgc lon fpa az] mjd-tt mjd-ut1 eop]
  (let [cl (math/cos latgc) sl (math/sin latgc) co (math/cos lon) so (math/sin lon)
        cf (math/cos fpa) sf (math/sin fpa) ca (math/cos az) sa (math/sin az)
        r [(* rm cl co) (* rm cl so) (* rm sl)]
        ;; velocity in SEZ-like local axes: up, east, north
        v [(* vm (+ (* co (- (* cl sf) (* sl cf ca))) (* (- so) cf sa)))
           (* vm (+ (* so (- (* cl sf) (* sl cf ca))) (* co cf sa)))
           (* vm (+ (* sl sf) (* cl cf ca)))]]
    (reduction/ecef->eci [r v] mjd-tt mjd-ut1 eop)))

(defn state->flight
  "`[rm vm latgc lon fpa az ra dec]`: flight elements of an inertial state,
  with its right ascension and declination."
  [s mjd-tt mjd-ut1 eop]
  (let [[r v] (reduction/eci->ecef s mjd-tt mjd-ut1 eop)
        rm (v3/length r) vm (v3/length v)
        [x y z] r
        latgc (math/asin (/ z rm))
        lon (math/atan2 y x)
        [_ ra dec] (state->radec s)
        up (v3/normalize r)
        east (v3/normalize [(- y) x 0.0])
        north (v3/cross up east)
        fpa (math/asin (/ (v3/dot v up) vm))
        az (math/atan2 (v3/dot v east) (v3/dot v north))]
    [rm vm latgc lon fpa az ra dec]))
