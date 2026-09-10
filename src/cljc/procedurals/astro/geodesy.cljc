(ns procedurals.astro.geodesy
  "Geodetic coordinates and local frames (Montenbruck & Gill 5.3).

  The Earth is an ellipsoid, not a sphere, and the distinction is not a
  refinement: geodetic latitude is measured from the local vertical -- the
  normal to the ellipsoid -- while geocentric latitude is measured from the
  centre, and at 45 degrees the two differ by 11.5 arcminutes, some 21 km
  along the surface. A position quoted in one and read as the other is
  wrong by that much.

  Local frames follow: east-north-up for a ground station, and the
  radial-transverse-normal frame that orbit work uses to describe where a
  satellite is relative to where it was expected."
  (:require [procedurals.astro.constants :as c]
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
  the centre, which is why the ellipsoid needs its own radius at all."
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
  "Latitude as seen from the centre, which is what a spherical model would
  give. It differs from the geodetic value by up to 11.5 arcminutes."
  [[x y z]]
  (math/atan2 z (math/sqrt (+ (* x x) (* y y)))))

;; ---------------------------------------------------------------- local frames

(defn- unit [v]
  (let [m (math/sqrt (reduce + (map * v v)))]
    (if (zero? m) v (mapv #(/ % m) v))))

(defn- cross [[a b cc] [d e f]]
  [(- (* b f) (* cc e)) (- (* cc d) (* a f)) (- (* a e) (* b d))])

(defn east-north-up
  "The local horizon frame at a geodetic latitude and longitude: rows are
  the east, north and up directions in Earth-fixed coordinates.

  Up is the ellipsoid normal, not the direction to the centre. On an
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
        d     (mapv - target station)
        [e n u] (mapv (fn [row] (reduce + (map * row d))) (east-north-up lat lon))
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
  radially, so after a day a prediction is typically wrong by kilometres in
  one direction and metres in the others."
  [r v]
  (let [R (unit r)
        N (unit (cross r v))
        T (cross N R)]
    [R T N]))

(defn to-rtn
  "Express a difference vector in the RTN frame of an orbit."
  [r v d]
  (mapv (fn [row] (reduce + (map * row d))) (rtn-frame r v)))
