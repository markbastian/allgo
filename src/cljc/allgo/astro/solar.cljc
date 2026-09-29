(ns allgo.astro.solar
  "The Sun: its position at two levels of accuracy, its rectangular
  coordinates, the equinoxes and solstices, the equation of time, and the
  orientation of its disk (Meeus, *Astronomical Algorithms*, chapters 25 to
  29).

  The low-accuracy position treats the Earth's orbit as a Kepler ellipse
  with slowly changing elements and is good to 0.01 degree. The
  high-accuracy one is the Earth's VSOP87 position turned around, with
  the corrections that make it apparent -- nutation, which moves the
  equinox it is measured from, and aberration, which displaces the Sun by
  20.5 arcseconds against the Earth's orbital motion. Good to about an
  arcsecond.

  `allgo.astro.ephemeris/sun` is a third, cruder series from Montenbruck
  and Gill, kept for the force models that need nothing better.

  Times are MJD in TT; longitudes and latitudes are radians, distances
  AU."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.frames :as frames]
            [allgo.astro.time :as time]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

(defn- deg [x] (* x c/degrees))

;; --------------------------------------------------- low accuracy (25.2-25.8)

(defn mean-anomaly [mjd-tt]
  (deg (poly/horner-ascending [357.52911 35999.05029 -0.0001537] (time/centuries-J2000 mjd-tt))))

(defn eccentricity
  "Eccentricity of the Earth's orbit, slowly decreasing."
  [mjd-tt]
  (poly/horner-ascending [0.016708634 -0.000042037 -0.0000001267] (time/centuries-J2000 mjd-tt)))

(defn- low-accuracy [mjd-tt]
  (let [T  (time/centuries-J2000 mjd-tt)
        L0 (deg (poly/horner-ascending [280.46646 36000.76983 0.0003032] T))
        M  (mean-anomaly mjd-tt)
        C  (deg (+ (* (poly/horner-ascending [1.914602 -0.004817 -0.000014] T) (math/sin M))
                   (* (- 0.019993 (* 0.000101 T)) (math/sin (* 2.0 M)))
                   (* 0.000289 (math/sin (* 3.0 M)))))]
    {:lon (am/wrap-2pi (+ L0 C)) :anomaly (am/wrap-2pi (+ M C))
     :node (deg (- 125.04 (* 1934.136 T)))}))

(defn true-longitude
  "The Sun's geometric longitude referred to the mean equinox of date, by
  the equation of the center -- 0.01 degree."
  [mjd-tt]
  (:lon (low-accuracy mjd-tt)))

(defn true-anomaly [mjd-tt] (:anomaly (low-accuracy mjd-tt)))

(defn radius
  "The Earth-Sun distance in AU, from the low-accuracy orbit."
  [mjd-tt]
  (let [e (eccentricity mjd-tt)]
    (/ (* 1.000001018 (- 1.0 (* e e))) (+ 1.0 (* e (math/cos (true-anomaly mjd-tt)))))))

(defn apparent-longitude-low
  "Apparent longitude from the low-accuracy theory: corrected for
  nutation and aberration together, by the dominant nutation term."
  [mjd-tt]
  (let [{:keys [lon node]} (low-accuracy mjd-tt)]
    (- lon (deg 0.00569) (* (deg 0.00478) (math/sin node)))))

(defn true-longitude-J2000
  "`true-longitude` referred to the equinox of J2000 instead of the date."
  [mjd-tt]
  (- (true-longitude mjd-tt) (* (deg 0.01397) (time/centuries-J2000 mjd-tt) 100.0)))

(defn equatorial-low
  "`[ra dec]` from the low-accuracy theory: geometric and mean, or when
  `apparent?`, apparent -- the obliquity then also corrected by the main
  nutation term."
  ([mjd-tt] (equatorial-low mjd-tt false))
  ([mjd-tt apparent?]
   (let [{:keys [node]} (low-accuracy mjd-tt)
         eps (cond-> (frames/mean-obliquity mjd-tt)
               apparent? (+ (* (deg 0.00256) (math/cos node))))
         lon (if apparent? (apparent-longitude-low mjd-tt) (true-longitude mjd-tt))]
     (coord/ecliptic->equatorial [lon 0.0] eps))))

;; ---------------------------------------------------- high accuracy (25.9)

(defn geometric
  "The Sun's geometric `[lon lat r]`, referred to the mean equinox of date
  and the FK5 frame: the Earth's heliocentric VSOP87 position, turned to
  point the other way."
  [mjd-tt]
  (let [[l b r] (vsop87/heliocentric :earth mjd-tt)
        [l b] (vsop87/->fk5 [(+ l math/PI) (- b)] mjd-tt)]
    [(am/wrap-2pi l) b r]))

(defn aberration
  "The Sun's aberration in longitude at distance `r` AU: the Earth moves
  some 30 km/s across the line of sight, so the Sun is seen where it was
  8.3 minutes ago, 20.5 arcseconds behind."
  [r]
  (/ (* -20.4898 c/arcsec) r))

(defn apparent
  "The Sun's apparent `[lon lat r]`: geometric, plus nutation and
  aberration."
  [mjd-tt]
  (let [[l b r] (geometric mjd-tt)
        [dpsi] (frames/nutation-angles mjd-tt)]
    [(am/wrap-2pi (+ l dpsi (aberration r))) b r]))

(defn apparent-equatorial
  "The Sun's apparent `[ra dec r]`, on the true equator of date."
  [mjd-tt]
  (let [[l b r] (apparent mjd-tt)
        [ra dec] (coord/ecliptic->equatorial [l b] (frames/true-obliquity mjd-tt))]
    [ra dec r]))

;; -------------------------------------------------------- rectangular (26)

(defn rectangular
  "Geocentric equatorial rectangular coordinates of the Sun, AU, referred
  to the mean equator and equinox of date."
  [mjd-tt]
  (let [[l b r] (geometric mjd-tt)
        eps (frames/mean-obliquity mjd-tt)
        se (math/sin eps) ce (math/cos eps)
        sl (math/sin l) sb (math/sin b) cb (math/cos b)]
    [(* r cb (math/cos l))
     (* r (- (* cb sl ce) (* sb se)))
     (* r (+ (* cb sl se) (* sb ce)))]))

(defn rectangular-J2000
  "The same, referred to the mean equator and equinox of J2000 -- the
  frame of the star catalogs and of everything else in this package."
  [mjd-tt]
  (lin/mat-vec (lin/transpose (frames/precession mjd-tt)) (rectangular mjd-tt)))

(defn rectangular-at
  "The same, referred to the mean equator and equinox of `epoch` (MJD):
  B1950, say, to compare with an old ephemeris."
  [mjd-tt epoch]
  (lin/mat-vec (frames/precession epoch) (rectangular-J2000 mjd-tt)))

;; ------------------------------------------------ equinoxes and solstices (27)

(def ^:private seasons-before-1000
  {:march     [1721139.29189 365242.13740 0.06134 0.00111 -0.00071]
   :june      [1721233.25401 365241.72562 -0.05232 0.00907 0.00025]
   :september [1721325.70455 365242.49558 -0.11677 -0.00297 0.00074]
   :december  [1721414.39987 365242.88257 -0.00769 -0.00933 -0.00006]})

(def ^:private seasons-after-1000
  {:march     [2451623.80984 365242.37404 0.05169 -0.00411 -0.00057]
   :june      [2451716.56767 365241.62603 0.00325 0.00888 -0.00030]
   :september [2451810.21715 365242.01767 -0.11575 0.00337 0.00078]
   :december  [2451900.05952 365242.74049 -0.06223 -0.00823 0.00032]})

(def ^:private seasons-periodic
  "Meeus's table 27.C: amplitude in 1e-5 day, phase and rate in degrees."
  [[485 324.96 1934.136] [203 337.23 32964.467] [199 342.08 20.186]
   [182 27.85 445267.112] [156 73.14 45036.886] [136 171.52 22518.443]
   [77 222.54 65928.934] [74 296.72 3034.906] [70 243.58 9037.513]
   [58 119.81 33718.147] [52 297.17 150.678] [50 21.02 2281.226]
   [45 247.54 29929.562] [44 325.15 31555.956] [29 60.93 4443.417]
   [18 155.12 67555.328] [17 288.79 4562.452] [16 198.04 62894.029]
   [14 199.76 31436.921] [12 95.39 14577.848] [12 287.11 31931.756]
   [12 320.81 34777.259] [9 227.73 1222.114] [8 15.45 16859.074]])

(defn- mean-season [year season]
  (let [[table y] (if (< year 1000) [seasons-before-1000 year] [seasons-after-1000 (- year 2000)])]
    (- (poly/horner-ascending (table season) (* 0.001 y)) c/jd-mjd-offset)))

(defn season
  "MJD (TT) of the March or September equinox or the June or December
  solstice of `year`: `season` is :march, :june, :september or :december.

  From Meeus's mean instants and 24 periodic terms, good to a minute
  between 1951 and 2050. `season-exact` iterates on the VSOP87 Sun
  instead, for a few seconds anywhere in its range."
  [year season]
  (let [J0 (mean-season year season)
        T  (time/centuries-J2000 J0)
        W  (deg (- (* 35999.373 T) 2.47))
        dl (+ 1.0 (* 0.0334 (math/cos W)) (* 0.0007 (math/cos (* 2.0 W))))
        S  (reduce (fn [s [a b cc]] (+ s (* a (math/cos (deg (+ b (* cc T))))))) 0.0
                   seasons-periodic)]
    (+ J0 (/ (* 0.00001 S) dl))))

(defn season-exact
  "`season`, found by correcting the mean instant until the Sun's apparent
  longitude is exactly 0, 90, 180 or 270 degrees (Meeus 27.1)."
  [year season]
  (let [q (deg ({:march 0 :june 90 :september 180 :december 270} season))]
    (loop [J (mean-season year season) i 0]
      (let [dJ (* 58.0 (math/sin (- q (first (apparent J)))))]
        (if (or (< (abs dJ) 5e-6) (>= i 20))
          (+ J dJ)
          (recur (+ J dJ) (inc i)))))))

;; ---------------------------------------------------- equation of time (28)

(defn- sun-mean-longitude [mjd-tt]
  (deg (poly/horner-ascending [280.4664567 360007.6982779 0.03032028 (/ 1.0 49931) (/ -1.0 15300)
                               (/ -1.0 2000000)] (vsop87/millennia mjd-tt))))

(defn equation-of-time
  "Apparent minus mean solar time, radians of hour angle (15 degrees an
  hour): how far a sundial runs ahead of a clock. Its two humps -- 16
  minutes fast in early November, 14 slow in February -- are the
  eccentricity of the orbit and the tilt of the axis, beating."
  [mjd-tt]
  (let [[ra] (apparent-equatorial mjd-tt)
        [dpsi deps] (frames/nutation-angles mjd-tt)
        eps (+ (frames/mean-obliquity mjd-tt) deps)
        E (- (sun-mean-longitude mjd-tt) (deg 0.0057183) ra (- (* dpsi (math/cos eps))))]
    (am/wrap-angle E)))

(defn equation-of-time-smart
  "The equation of time from W. M. Smart's series in the orbital elements,
  good to a few seconds: no position of the Sun needed at all."
  [mjd-tt]
  (let [eps (frames/mean-obliquity mjd-tt)
        y   (let [t (math/tan (* 0.5 eps))] (* t t))
        L0  (sun-mean-longitude mjd-tt)
        e   (eccentricity mjd-tt)
        M   (mean-anomaly mjd-tt)
        s2L (math/sin (* 2.0 L0)) c2L (math/cos (* 2.0 L0))
        sM  (math/sin M)]
    (- (+ (* y s2L) (* -2.0 e sM) (* 4.0 e y sM c2L))
       (* y y s2L c2L) (* 1.25 e e (math/sin (* 2.0 M))))))

;; --------------------------------------- physical observations of the Sun (29)

(defn disk
  "Orientation of the Sun's disk at `mjd` (TT): `[P B0 L0]` -- the position
  angle of its rotation axis, measured east from the north point of the
  disk; the heliographic latitude of the disk's center, which swings by
  7.25 degrees a year because the solar equator is tilted to the
  ecliptic; and the heliographic longitude of the center, which runs
  backward through a full turn every 27.3 days as the Sun rotates."
  [mjd]
  (let [jd  (+ mjd c/jd-mjd-offset)
        th  (am/wrap-2pi (* (/ (- jd 2398220.0) 25.38) c/two-pi))
        I   (deg 7.25)
        K   (+ (deg 73.6667) (* (deg 1.3958333) (/ (- jd 2396758.0) 36525.0)))
        [L _ R] (geometric mjd)
        [dpsi deps] (frames/nutation-angles mjd)
        eps (+ (frames/mean-obliquity mjd) deps)
        lam (+ L (aberration R))
        lp  (+ lam dpsi)
        slk (math/sin (- lam K)) clk (math/cos (- lam K))
        P   (+ (math/atan (- (* (math/cos lp) (math/tan eps))))
               (math/atan (- (* clk (math/tan I)))))
        B0  (math/asin (* slk (math/sin I)))
        eta (math/atan2 (- (* slk (math/cos I))) (- clk))]
    [P B0 (am/wrap-2pi (- eta th))]))

(defn carrington-rotation
  "MJD (TT) at which Carrington rotation `n` begins -- the numbering of
  solar rotations from 1853 November 9, when Carrington started counting,
  and still the way sunspot records are indexed."
  [n]
  (let [m (deg (+ 281.96 (* 26.882476 n)))]
    (- (+ 2398140.227 (* 27.2752316 n)
          (* 0.1454 (math/sin m)) (* -0.0085 (math/sin (* 2.0 m))) (* -0.0141 (math/cos (* 2.0 m))))
       c/jd-mjd-offset)))
