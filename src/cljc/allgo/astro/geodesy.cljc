(ns allgo.astro.geodesy
  "Geodetic coordinates and local frames (Montenbruck & Gill 5.3).

  The Earth is an ellipsoid, not a sphere, and the distinction is not a
  refinement: geodetic latitude is measured from the local vertical -- the
  normal to the ellipsoid -- while geocentric latitude is measured from the
  center, and at 45 degrees the two differ by 11.5 arcminutes, some 21 km
  along the surface. A position quoted in one and read as the other is
  wrong by that much.

  Local frames follow: east-north-up for a ground station, and the
  radial-transverse-normal frame that orbit work uses to describe where a
  satellite is relative to where it was expected."
  (:require [allgo.astro.constants :as c]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def a-earth c/R-earth)
(def b-earth (* c/R-earth (- 1.0 c/flattening)))

(def e2
  "First eccentricity squared of the reference ellipsoid."
  (- (* 2.0 c/flattening) (* c/flattening c/flattening)))

(def ep2
  "Second eccentricity squared, which appears in the inverse problem."
  (/ (- (* a-earth a-earth) (* b-earth b-earth)) (* b-earth b-earth)))

(defn- prime-vertical
  "Radius of curvature in the prime vertical -- the distance from the
  surface to where the local normal cuts the spin axis. Not the distance to
  the center, which is why the ellipsoid needs its own radius at all."
  [sin-lat]
  (/ a-earth (math/sqrt (- 1.0 (* e2 sin-lat sin-lat)))))

(defn geodetic->cartesian
  "Earth-fixed Cartesian position from geodetic latitude and longitude
  (radians) and height above the ellipsoid (km)."
  [lat lon h]
  (let [sl (math/sin lat) cl (math/cos lat)
        N  (prime-vertical sl)]
    [(* (+ N h) cl (math/cos lon))
     (* (+ N h) cl (math/sin lon))
     (* (+ (* N (- 1.0 e2)) h) sl)]))

(defn cartesian->geodetic
  "Geodetic latitude, longitude and height from an Earth-fixed position.

  Solved by iteration rather than in closed form: latitude and height are
  coupled, since the height is measured along a normal whose direction
  depends on the latitude. It converges in a handful of passes anywhere
  from the surface to geostationary."
  [[x y z]]
  (let [p (math/sqrt (+ (* x x) (* y y)))]
    (if (< p 1e-9)
      ;; on the axis: longitude is undefined and latitude is a pole
      [(if (neg? z) (- (/ math/PI 2.0)) (/ math/PI 2.0)) 0.0 (- (abs z) b-earth)]
      (loop [lat (math/atan2 z (* p (- 1.0 e2))) n 0]
        (let [N (prime-vertical (math/sin lat))
              h (- (/ p (math/cos lat)) N)
              lat' (math/atan2 z (* p (- 1.0 (* e2 (/ N (+ N h))))))]
          (if (or (>= n 8) (< (abs (- lat' lat)) 1e-14))
            [lat' (math/atan2 y x) (- (/ p (math/cos lat')) (prime-vertical (math/sin lat')))]
            (recur lat' (inc n))))))))

(defn geocentric-latitude
  "Latitude as seen from the center, which is what a spherical model would
  give. It differs from the geodetic value by up to 11.5 arcminutes."
  [[x y z]]
  (math/atan2 z (math/sqrt (+ (* x x) (* y y)))))

;; ------------------------------------------------ the Earth's globe (Meeus 11)

(defn parallax-constants
  "`[rho-sin rho-cos]` -- rho sin phi' and rho cos phi', the observer's
  distance from the Earth's axis and from the equatorial plane in
  equatorial radii -- for geodetic latitude `lat` and height `h` km. They
  are what diurnal parallax and eclipse work need, and are simply the
  observer's Earth-fixed position scaled and projected."
  [lat h]
  (let [[x y z] (geodetic->cartesian lat 0.0 h)]
    [(/ z a-earth) (/ (math/hypot x y) a-earth)]))

(defn radius-of-parallel
  "Radius of the circle of latitude `lat` on the ellipsoid, km. One degree
  of longitude there is this times pi/180."
  [lat]
  (* (math/cos lat) (prime-vertical (math/sin lat))))

(defn radius-of-curvature
  "Radius of curvature of the meridian at latitude `lat`, km. One degree
  of latitude is this times pi/180 -- longer at the poles than at the
  equator, which is how the Earth's flattening was first measured."
  [lat]
  (let [s (math/sin lat)]
    (/ (* a-earth (- 1.0 e2)) (math/pow (- 1.0 (* e2 s s)) 1.5))))

(defn geodetic->geocentric-latitude
  "Geocentric latitude of a point on the surface at geodetic latitude
  `lat`."
  [lat]
  (math/atan (* (- 1.0 e2) (math/tan lat))))

(defn distance
  "Distance along the ellipsoid between two places given as geodetic
  `[lat lon]`, km, by Andoyer's formula as Meeus gives it. Good to about
  50 meters on intercontinental distances -- the flattening enters only to
  first order -- and unusable for points nearly antipodal."
  [[lat1 lon1] [lat2 lon2]]
  (let [sq (fn [x] (* x x))
        F (* 0.5 (+ lat1 lat2)) G (* 0.5 (- lat1 lat2)) L (* 0.5 (- lon1 lon2))
        s2F (sq (math/sin F)) c2F (sq (math/cos F))
        s2G (sq (math/sin G)) c2G (sq (math/cos G))
        s2L (sq (math/sin L)) c2L (sq (math/cos L))
        S (+ (* s2G c2L) (* c2F s2L))
        C (+ (* c2G c2L) (* s2F s2L))
        w (math/atan (math/sqrt (/ S C)))
        R (/ (math/sqrt (* S C)) w)
        H1 (/ (- (* 3.0 R) 1.0) (* 2.0 C))
        H2 (/ (+ (* 3.0 R) 1.0) (* 2.0 S))]
    (* 2.0 w a-earth (+ 1.0 (* c/flattening (- (* H1 s2F c2G) (* H2 c2F s2G)))))))

;; ---------------------------------------------------------------- local frames

(defn east-north-up
  "The local horizon frame at a geodetic latitude and longitude: rows are
  the east, north and up directions in Earth-fixed coordinates.

  Up is the ellipsoid normal, not the direction to the center. On an
  ellipsoid a plumb line does not point at the middle of the Earth."
  [lat lon]
  (let [sl (math/sin lat) cl (math/cos lat)
        so (math/sin lon) co (math/cos lon)]
    [[(- so) co 0.0]
     [(* (- sl) co) (* (- sl) so) cl]
     [(* cl co) (* cl so) sl]]))

(defn look-angles
  "Azimuth, elevation and range from a ground station to a target, both
  given in Earth-fixed coordinates.

  Azimuth is measured from north through east, the surveying convention.
  Elevation is negative when the target is below the horizon, which for a
  satellite is most of the time."
  [station target]
  (let [[lat lon _] (cartesian->geodetic station)
        d     (v3/sub target station)
        [e n u] (lin/mat-vec (east-north-up lat lon) d)
        rng   (math/sqrt (+ (* e e) (* n n) (* u u)))]
    {:azimuth   (let [az (math/atan2 e n)] (if (neg? az) (+ az c/two-pi) az))
     :elevation (math/asin (/ u rng))
     :range     rng}))

(defn visible?
  "Whether a target clears a station's horizon by at least `mask` radians.
  Real stations use five or ten degrees, since the atmosphere near the
  horizon is a long slant path and refraction is worst there."
  ([station target] (visible? station target 0.0))
  ([station target mask]
   (> (:elevation (look-angles station target)) mask)))

(defn rtn-frame
  "The radial-transverse-normal frame of an orbit: rows are the outward
  radial direction, the along-track direction and the orbit normal.

  This is the frame orbit errors are quoted in, because they are not
  isotropic. A small error in speed grows along-track and hardly at all
  radially, so after a day a prediction is typically wrong by kilometers in
  one direction and meters in the others."
  [r v]
  (let [R (v3/normalize r)
        N (v3/normalize (v3/cross r v))
        T (v3/cross N R)]
    [R T N]))

(defn to-rtn
  "Express a difference vector in the RTN frame of an orbit."
  [r v d]
  (lin/mat-vec (rtn-frame r v) d))
