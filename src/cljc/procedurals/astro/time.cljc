(ns procedurals.astro.time
  "Dates and sidereal time, enough to place a satellite relative to a
  rotating Earth (Montenbruck & Gill chapter 5, in service of chapter 3).

  The geopotential is defined in an Earth-fixed frame while the equation of
  motion is integrated in an inertial one, so every gravity evaluation needs
  to know how far the Earth has turned. That is what Greenwich sidereal time
  supplies."
  (:require [procedurals.astro.constants :as c]
            [clojure.math :as math]))

(defn calendar->mjd
  "Modified Julian Date from a calendar date and fractional hour, valid across
  the Gregorian reform. Uses the standard integer-arithmetic form rather than
  a table of month lengths."
  ([year month day] (calendar->mjd year month day 0.0))
  ([year month day hour]
   (let [[y m] (if (<= month 2) [(dec year) (+ month 12)] [year month])
         ;; Gregorian calendar from 1582-10-15; Julian before it
         gregorian? (>= (+ (* 10000.0 year) (* 100.0 month) day) 15821015.0)
         b (if gregorian?
             (let [a (long (math/floor (/ y 100.0)))]
               (+ (- 2 a) (long (math/floor (/ a 4.0)))))
             -2)
         jd (+ (long (math/floor (* 365.25 (+ y 4716))))
               (long (math/floor (* 30.6001 (inc m))))
               day b -1524.5
               (/ hour 24.0))]
     (- jd c/jd-mjd-offset))))

(defn mjd->jd [mjd] (+ mjd c/jd-mjd-offset))

(defn centuries-J2000
  "Julian centuries since J2000.0 -- the argument every polynomial in the
  reduction formulae is expanded in."
  [mjd]
  (/ (- mjd c/mjd-J2000) 36525.0))

(defn gmst
  "Greenwich Mean Sidereal Time in radians, from UT1 as an MJD.

  The linear term is 360.98564736629 degrees a day rather than 360: the
  Earth turns once relative to the stars in slightly less than a solar day,
  and that excess is the whole reason sidereal time exists."
  [mjd-ut1]
  (let [t   (centuries-J2000 mjd-ut1)
        deg (+ 280.46061837
               (* 360.98564736629 (- (mjd->jd mjd-ut1) 2451545.0))
               (* 0.000387933 t t)
               (- (/ (* t t t) 38710000.0)))
        rad (* (rem deg 360.0) c/degrees)]
    (if (neg? rad) (+ rad c/two-pi) rad)))

;; --------------------------------------------------------------- time scales
;;
;; Montenbruck & Gill 5.1. Four different questions get four different
;; answers, and orbit work needs all of them:
;;
;;   TAI  atomic time, the one that just counts seconds
;;   UTC  civil time, atomic but stepped by leap seconds to track the Earth
;;   UT1  the Earth's actual rotation angle, which wanders unpredictably
;;   TT   the argument of the ephemerides, TAI + 32.184 s
;;
;; Gravity and the ephemerides run on TT; the Earth's orientation runs on
;; UT1. Confusing them is a a fraction of a second, which at 465 m/s of
;; equatorial rotation is a few hundred metres of position.

(def leap-seconds
  "TAI - UTC in whole seconds, from each MJD at which it changed.

  Leap seconds are inserted by decree, not by formula: they keep UTC within
  0.9 s of the Earth's rotation, which is slowing at a rate nobody can
  predict far ahead. So this is a table that has to be maintained, and any
  epoch past its end is an extrapolation."
  [[41317.0 10] [41499.0 11] [41683.0 12] [42048.0 13] [42413.0 14]
   [42778.0 15] [43144.0 16] [43509.0 17] [43874.0 18] [44239.0 19]
   [44786.0 20] [45151.0 21] [45516.0 22] [46247.0 23] [47161.0 24]
   [47892.0 25] [48257.0 26] [48804.0 27] [49169.0 28] [49534.0 29]
   [50083.0 30] [50630.0 31] [51179.0 32] [53736.0 33] [54832.0 34]
   [56109.0 35] [57204.0 36] [57754.0 37]])

(defn tai-utc
  "TAI - UTC in seconds at a UTC epoch. Zero before 1972, when UTC began
  stepping in whole seconds."
  [mjd-utc]
  (or (some (fn [[at n]] (when (>= mjd-utc at) n)) (reverse leap-seconds))
      0))

(def tt-tai
  "TT - TAI, seconds. Fixed by definition, not measured: the offset was
  chosen in 1977 so terrestrial time would continue the ephemeris time
  scale that preceded it without a jump."
  32.184)

(def gps-tai
  "GPS time - TAI, seconds. GPS started from UTC in January 1980 and has
  ignored leap seconds ever since, so it runs 19 s behind TAI forever."
  -19.0)

(defn- shift [mjd seconds] (+ mjd (/ seconds 86400.0)))

(defn utc->tai [mjd-utc] (shift mjd-utc (tai-utc mjd-utc)))
(defn tai->tt  [mjd-tai] (shift mjd-tai tt-tai))
(defn utc->tt  [mjd-utc] (tai->tt (utc->tai mjd-utc)))
(defn utc->gps [mjd-utc] (shift (utc->tai mjd-utc) gps-tai))

(defn tt->utc
  "TT back to UTC. The leap-second count is looked up from an approximate
  UTC first, since the table is indexed by UTC and the offset is what we are
  trying to find -- one pass is ample, the count being constant except
  within a second of a leap."
  [mjd-tt]
  (let [approx (shift mjd-tt (- tt-tai))]
    (shift approx (- (tai-utc approx)))))

(defn utc->ut1
  "UT1 from UTC given the observed difference in seconds.

  `dut1` cannot be computed, only measured and published: it is the Earth's
  rotation running fast or slow against atomic time, and it drifts by a
  millisecond a day in ways that depend on the weather and the core. Zero is
  a reasonable default at the metre level, and wrong at the centimetre one."
  ([mjd-utc] mjd-utc)
  ([mjd-utc dut1] (shift mjd-utc dut1)))

(defn tt->tdb
  "Barycentric Dynamical Time from Terrestrial Time, MJD.

  The two differ by at most 1.7 ms, a relativistic effect: a clock on Earth
  runs at a varying rate relative to one at the solar system barycentre,
  because the Earth's distance from the Sun and its speed both vary over the
  year. The dominant term is annual and follows Earth's mean anomaly."
  [mjd-tt]
  (let [g (* c/degrees (+ 357.53 (* 0.9856003 (- (mjd->jd mjd-tt) 2451545.0))))]
    (shift mjd-tt (+ (* 0.001658 (math/sin g))
                     (* 0.000014 (math/sin (* 2.0 g)))))))
