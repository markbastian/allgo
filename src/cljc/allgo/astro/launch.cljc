(ns allgo.astro.launch
  "From a launch site onto an interplanetary asymptote: the direction to
  fly and the times of day it can be done (Burke, Falck and McGuire,
  NASA/TM-2010-216764, equations 3-6; the geometry is that of any mission
  design text).

  A departure must end on the hyperbola whose outgoing asymptote points
  along the excess velocity, at right ascension RLA and declination DLA.
  The parking orbit holds that asymptote in its plane, and at launch it
  holds the site too, so its plane is the one through the Earth's center,
  the site and the asymptote's direction on the sky. The azimuth the
  rocket flies is where that great circle leaves the site, and the
  parking orbit's inclination follows from it: cos i = cos(latitude)
  sin(azimuth). As the Earth turns the site through the day the azimuth
  swings, and a range's safety limits -- 40 to 115 degrees at Cape
  Canaveral -- leave the hours when it is allowed: the daily launch
  window.

  The plane can be flown either way round. Heading toward the asymptote
  or away from it, the parking orbit coasts to the injection point ahead
  of the asymptote by the hyperbola's asymptote angle; each is a
  separate opportunity.

  Angles in radians; the asymptote and the site referred to the equator
  and equinox of date (`interplanetary/asymptote` with a date); times UT1
  as MJD."
  (:require [allgo.astro.time :as time]
            [allgo.math :as am]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

(defn azimuth
  "The azimuth, east of north, at which the great circle from a site at
  latitude `lat` and right ascension `ra-site` heads toward the asymptote
  at `rla`, `dla` (handbook equation 3): cot az = (cos lat tan dla - sin
  lat cos(rla - ra-site)) / sin(rla - ra-site). nil when the site sits on
  the asymptote's line, where every plane holds both."
  [lat ra-site rla dla]
  (let [d (- rla ra-site)
        y (math/sin d)
        x (- (* (math/cos lat) (math/tan dla)) (* (math/sin lat) (math/cos d)))]
    (when (or (> (abs x) 1e-15) (> (abs y) 1e-15))
      (am/wrap-2pi (math/atan2 y x)))))

(defn inclination
  "The inclination of the orbit launched from latitude `lat` at `azimuth`:
  cos i = cos lat sin azimuth. Due east gives the least, the latitude
  itself."
  [lat azimuth]
  (math/acos (* (math/cos lat) (math/sin azimuth))))

(defn site-right-ascension
  "The right ascension of a site at east longitude `lon` at `mjd-ut1`: its
  longitude plus Greenwich sidereal time (equation 5, with the sidereal
  time from `allgo.astro.time/gmst` rather than the handbook's linear
  approximation)."
  [lon mjd-ut1]
  (am/wrap-2pi (+ lon (time/gmst mjd-ut1))))

(defn- in-range? [az [lo hi]]
  (and az (<= lo az hi)))

(defn headings
  "The two azimuths along the parking-orbit plane at `mjd-ut1` -- toward
  the asymptote and away from it -- that fall within `limits`
  `[least most]`: `{:toward az :away az}` with those outside left out."
  [lat lon rla dla mjd-ut1 limits]
  (when-let [az (azimuth lat (site-right-ascension lon mjd-ut1) rla dla)]
    (let [away (am/wrap-2pi (+ az math/PI))]
      (cond-> {}
        (in-range? az limits) (assoc :toward az)
        (in-range? away limits) (assoc :away away)))))

(defn windows
  "The launch windows between `start` and `end` (UT1 MJD) from a site at
  latitude `lat` and east longitude `lon` onto the asymptote `rla`, `dla`,
  flying within the azimuth `limits` `[least most]`: `[{:heading :open
  :close} ...]`, heading `:toward` or `:away` from the asymptote, each
  open and close to about a second. The asymptote is taken as fixed
  through the span, as it nearly is over a day (handbook, launch/injection
  geometry)."
  ([lat lon rla dla start end limits] (windows lat lon rla dla start end limits {}))
  ([lat lon rla dla start end limits {:keys [step] :or {step (/ 1.0 1440.0)}}]
   (let [tol (/ 1.0 86400.0)]
     (vec (for [heading [:toward :away]
                :let [open? #(contains? (headings lat lon rla dla % limits) heading)
                      edges (roots/transitions open? start end step tol)
                      ;; a window already open at the start opens there
                      opens (cond->> (keep (fn [[t was]] (when-not was t)) edges)
                              (open? start) (cons start))
                      closes (concat (keep (fn [[t was]] (when was t)) edges)
                                     (when (open? end) [end]))]
                [o c] (map vector opens closes)]
            {:heading heading :open o :close c})))))
