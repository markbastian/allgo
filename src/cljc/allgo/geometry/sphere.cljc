(ns allgo.geometry.sphere
  "Geometry on and about a sphere: great circles between two places given
  by latitude and longitude, and where lines meet or miss a sphere
  centered at the origin. Angles radians, longitudes positive east,
  azimuths from north through east."
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

(defn great-circle
  "`{:angle :azimuth}` from `[lat1 lon1]` to `[lat2 lon2]`: the angle the
  great-circle arc subtends at the center, by the haversine, and the
  initial azimuth, by the four-part formula."
  [[lat1 lon1] [lat2 lon2]]
  (let [dlon (- lon2 lon1)
        hav (+ (am/sq (math/sin (* 0.5 (- lat2 lat1))))
               (* (math/cos lat1) (math/cos lat2) (am/sq (math/sin (* 0.5 dlon)))))
        angle (* 2.0 (math/asin (math/sqrt (min 1.0 hav))))
        az (math/atan2 (* (math/sin dlon) (math/cos lat2))
                       (- (* (math/cos lat1) (math/sin lat2))
                          (* (math/sin lat1) (math/cos lat2) (math/cos dlon))))]
    {:angle angle :azimuth (am/wrap-2pi az)}))

(defn destination
  "`[lat lon]` reached from `[lat lon]` by going through `angle` along the
  great circle that sets out at `azimuth`."
  [[lat lon] angle azimuth]
  (let [d angle
        lat2 (math/asin (+ (* (math/sin lat) (math/cos d))
                           (* (math/cos lat) (math/sin d) (math/cos azimuth))))
        lon2 (+ lon (math/atan2 (* (math/sin azimuth) (math/sin d) (math/cos lat))
                                (- (math/cos d) (* (math/sin lat) (math/sin lat2)))))]
    [lat2 (am/wrap-angle lon2)]))

(defn segment-clears?
  "Whether the segment from `a` to `b` stays outside the sphere of radius
  `radius` about the origin: both ends outside, and the point of the
  segment nearest the center too."
  [a b radius]
  (let [d (v3/sub b a)
        t (let [dd (v3/dot d d)] (if (zero? dd) 0.0 (- (/ (v3/dot a d) dd))))
        nearest (if (<= 0.0 t 1.0) (v3/add-scaled a d t) (if (< t 0.0) a b))]
    (and (>= (v3/length a) radius) (>= (v3/length b) radius)
         (>= (v3/length nearest) radius))))

(defn far-intersection
  "The distance along the unit direction `u` from `origin` at which the
  line reaches radius `radius` from the center -- the far root of
  |origin + t u| = radius -- nil where it never gets that far out."
  [origin u radius]
  (let [c (* 2.0 (v3/dot u origin))
        disc (- (* c c) (* 4.0 (- (v3/dot origin origin) (* radius radius))))]
    (when (>= disc 0.0) (* 0.5 (+ (- c) (math/sqrt disc))))))
