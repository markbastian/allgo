(ns allgo.astro.planets
  "Planetary positions from low-precision Keplerian elements (Vallado; the
  table originates with Standish and appears in the Explanatory Supplement
  and the Astronomical Almanac).

  Each planet is described by six elements and six linear rates per century.
  Evaluate the elements at the epoch, solve Kepler's equation, and the
  result is a heliocentric position -- so this leans directly on chapter 2
  rather than carrying its own machinery. Good to a few arcminutes across
  1800-2050, which is the accuracy the table is fitted for and ample for
  drawing the solar system.

  Everything returned here is referred to the mean equator and equinox of
  J2000, the same frame as `allgo.astro.ephemeris` and everything
  else in this package. The elements themselves are ecliptic; the rotation
  to equatorial uses the J2000 obliquity, not the obliquity of date, since
  the longitudes are measured from the fixed J2000 equinox. Mixing those is
  worth 8000 km at the Sun by 2025."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.kepler :as kep]
            [allgo.astro.time :as time]
            [clojure.math :as math]))

(def elements
  "Semi-major axis (AU), eccentricity, inclination, mean longitude,
  longitude of perihelion and longitude of the ascending node -- each as a
  value at J2000 and a rate per Julian century. Angles in degrees.

  Note the Earth row is the Earth-Moon barycentre, which is what the fit
  actually tracks; the Earth itself wanders about it by some 4700 km."
  {:mercury {:a [0.38709927  0.00000037] :e [0.20563593  0.00001906]
             :i [7.00497902 -0.00594749] :L [252.25032350 149472.67411175]
             :peri [77.45779628 0.16047689] :node [48.33076593 -0.12534081]}
   :venus   {:a [0.72333566  0.00000390] :e [0.00677672 -0.00004107]
             :i [3.39467605 -0.00078890] :L [181.97909950 58517.81538729]
             :peri [131.60246718 0.00268329] :node [76.67984255 -0.27769418]}
   :earth   {:a [1.00000261  0.00000562] :e [0.01671123 -0.00004392]
             :i [-0.00001531 -0.01294668] :L [100.46457166 35999.37244981]
             :peri [102.93768193 0.32327364] :node [0.0 0.0]}
   :mars    {:a [1.52371034  0.00001847] :e [0.09339410  0.00007882]
             :i [1.84969142 -0.00813131] :L [-4.55343205 19140.30268499]
             :peri [-23.94362959 0.44441088] :node [49.55953891 -0.29257343]}
   :jupiter {:a [5.20288700 -0.00011607] :e [0.04838624 -0.00013253]
             :i [1.30439695 -0.00183714] :L [34.39644051 3034.74612775]
             :peri [14.72847983 0.21252668] :node [100.47390909 0.20469106]}
   :saturn  {:a [9.53667594 -0.00125060] :e [0.05386179 -0.00050991]
             :i [2.48599187  0.00193609] :L [49.95424423 1222.49362201]
             :peri [92.59887831 -0.41897216] :node [113.66242448 -0.28867794]}
   :uranus  {:a [19.18916464 -0.00196176] :e [0.04725744 -0.00004397]
             :i [0.77263783 -0.00242939] :L [313.23810451 428.48202785]
             :peri [170.95427630 0.40805281] :node [74.01692503 0.04240589]}
   :neptune {:a [30.06992276  0.00026291] :e [0.00859048  0.00005105]
             :i [1.77004347  0.00035372] :L [-55.12002969 218.45945325]
             :peri [44.96476227 -0.32241464] :node [131.78422574 -0.00508664]}})

(def order [:mercury :venus :earth :mars :jupiter :saturn :uranus :neptune])

(defn- wrap-180 [deg] (- (mod (+ deg 180.0) 360.0) 180.0))

(defn elements-at
  "The classical elements of `planet` at `mjd-tt`, in the units the rest of
  this package uses: kilometres and radians."
  [planet mjd-tt]
  (let [T   (time/centuries-J2000 mjd-tt)
        at  (fn [k] (let [[v0 rate] (get-in elements [planet k])] (+ v0 (* rate T))))
        a   (* (at :a) c/AU)
        e   (at :e)
        i   (* (at :i) c/degrees)
        L   (at :L)
        per (at :peri)
        nod (at :node)]
    {:a a :e e :i i
     :raan (* (mod nod 360.0) c/degrees)
     ;; argument of periapsis is the longitude of perihelion less the node
     :argp (* (mod (- per nod) 360.0) c/degrees)
     :M    (* (wrap-180 (- L per)) c/degrees)}))

(defn heliocentric
  "Position of `planet` relative to the Sun, km, in equatorial J2000."
  [planet mjd-tt]
  (let [{:keys [a e i raan argp M]} (elements-at planet mjd-tt)
        nu (kep/mean->true M e)
        ;; Kepler's own routine, in the ecliptic frame the elements live in.
        ;; Position does not depend on mu -- only the discarded velocity does
        ;; -- but passing the real one keeps the call honest.
        [r _] (kep/elements->state c/GM-sun {:a a :e e :i i :raan raan :argp argp :nu nu})
        ce (math/cos eph/obliquity-J2000)
        se (math/sin eph/obliquity-J2000)
        [x y z] r]
    [x (- (* ce y) (* se z)) (+ (* se y) (* ce z))]))

(defn geocentric
  "Position of `planet` as seen from the Earth, km, in equatorial J2000.

  The Earth's own place is taken from the same table, so the two share a
  frame and their difference is meaningful even though neither is precise."
  [planet mjd-tt]
  (mapv - (heliocentric planet mjd-tt) (heliocentric :earth mjd-tt)))

(defn sun-from-earth
  "The Sun's geocentric position implied by this table -- simply the
  negative of the Earth's heliocentric position.

  Worth having because it duplicates `ephemeris/sun`, which comes from an
  entirely separate series. Two independent routes agreeing is a check that
  neither alone can give."
  [mjd-tt]
  (mapv - (heliocentric :earth mjd-tt)))

(defn period
  "Orbital period in days, from the semi-major axis and Kepler's third law."
  [planet mjd-tt]
  (/ (kep/period c/GM-sun (:a (elements-at planet mjd-tt))) 86400.0))
