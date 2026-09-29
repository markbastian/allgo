(ns allgo.astro.time
  "Dates and sidereal time, enough to place a satellite relative to a
  rotating Earth (Montenbruck & Gill chapter 5, in service of chapter 3).

  The geopotential is defined in an Earth-fixed frame while the equation of
  motion is integrated in an inertial one, so every gravity evaluation needs
  to know how far the Earth has turned. That is what Greenwich sidereal time
  supplies."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.interpolation :as interp]
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
             0)
         jd (+ (long (math/floor (* 365.25 (+ y 4716))))
               (long (math/floor (* 30.6001 (inc m))))
               day b -1524.5
               (/ hour 24.0))]
     (- jd c/jd-mjd-offset))))

(defn mjd->jd [mjd] (+ mjd c/jd-mjd-offset))

(defn mjd->calendar
  "Calendar date and fractional hour from a Modified Julian Date -- the
  inverse of `calendar->mjd`, by the standard integer reduction. Returns
  `[year month day hour]`."
  [mjd]
  (let [jd (+ mjd c/jd-mjd-offset 0.5)
        a  (long (math/floor jd))
        f  (- jd a)
        cc (if (< a 2299161)
             (+ a 1524)
             (let [b (long (math/floor (/ (- a 1867216.25) 36524.25)))]
               (+ a b (- (long (math/floor (/ b 4.0)))) 1525)))
        d  (long (math/floor (/ (- cc 122.1) 365.25)))
        e  (long (math/floor (* 365.25 d)))
        g  (long (math/floor (/ (- cc e) 30.6001)))
        day   (- cc e (long (math/floor (* 30.6001 g))))
        month (- g 1 (* 12 (long (math/floor (/ g 14.0)))))
        year  (- d 4715 (long (math/floor (/ (+ 7 month) 10.0))))]
    [year month day (* 24.0 f)]))

;; ------------------------------------------------ the calendar (Meeus 7)

(defn julian-leap-year?
  "Every fourth year, which is all the Julian calendar knew."
  [year]
  (zero? (mod year 4)))

(defn gregorian-leap-year?
  "Every fourth year except centuries, except every fourth century -- the
  correction of 1582, which drops three days in four hundred years to keep
  the equinox near March 21."
  [year]
  (or (and (zero? (mod year 4)) (not (zero? (mod year 100))))
      (zero? (mod year 400))))

(defn leap-year?
  "Leap year in whichever calendar was in force: Julian through 1582,
  Gregorian from 1583, the switch `calendar->mjd` also makes."
  [year]
  (if (< year 1583) (julian-leap-year? year) (gregorian-leap-year? year)))

(defn day-of-week
  "0 for Sunday through 6 for Saturday. Unaffected by any calendar reform:
  the week has run unbroken through all of them."
  [mjd]
  ;; MJD 0, 1858 November 17, was a Wednesday
  (long (mod (+ (math/floor mjd) 3) 7)))

(defn- whole-months [month k]
  (- (quot (* 275 month) 9) (* k (quot (+ month 9) 12)) 30))

(defn day-number
  "Day of the year, 1 for January 1, given whether the year is a leap year."
  [month day leap?]
  (+ (whole-months month (if leap? 1 2)) day))

(defn day-of-year
  "Day of the year, 1 for January 1."
  [year month day]
  (day-number month day (leap-year? year)))

(defn day-of-year->date
  "`[month day]` of day `n` of the year, the inverse of `day-number`."
  [n leap?]
  (let [k     (if leap? 1 2)
        month (if (< n 32) 1 (quot (+ (* 900 (+ k n)) (* 98 275)) 27500))]
    [month (- n (whole-months month k))]))

(defn decimal-year
  "The calendar year and the fraction of it elapsed at `mjd`, as ΔT
  tables and precession formulae are indexed."
  [mjd]
  (let [[y m d] (mjd->calendar mjd)]
    (+ y (/ (dec (+ (day-number m (long d) (leap-year? y))
                    (- d (math/floor d))))
            (if (leap-year? y) 366.0 365.0)))))

(defn centuries-J2000
  "Julian centuries since J2000.0 -- the argument every polynomial in the
  reduction formulae is expanded in."
  [mjd]
  (/ (- mjd c/mjd-J2000) 36525.0))

(defn julian-epoch->mjd
  "MJD of a Julian epoch such as 2000.0 or 2050.5: years of exactly
  365.25 days from J2000.0, the unit precession is expanded in."
  [year]
  (+ c/mjd-J2000 (* 365.25 (- year 2000.0))))

(defn mjd->julian-epoch [mjd] (+ 2000.0 (/ (- mjd c/mjd-J2000) 365.25)))

(defn besselian-epoch->mjd
  "MJD of a Besselian epoch such as 1950.0, the older unit of tropical
  years that the FK4 catalog and the galactic system are referred to."
  [year]
  (+ 33281.923459 (* 365.242198781 (- year 1950.0))))

(defn gmst
  "Greenwich Mean Sidereal Time in radians, from UT1 as an MJD.

  The linear term is 360.98564736629 degrees a day rather than 360: the
  Earth turns once relative to the stars in slightly less than a solar day,
  and that excess is the whole reason sidereal time exists."
  [mjd-ut1]
  (let [t   (centuries-J2000 mjd-ut1)
        deg (+ 280.46061837
               ;; days from J2000 taken in MJD: going through the JD first rounds
               ;; them to 5e-10 days, a few 1e-9 radians of sidereal time
               (* 360.98564736629 (- mjd-ut1 51544.5))
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
;; equatorial rotation is a few hundred meters of position.

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
  a reasonable default at the meter level, and wrong at the centimeter one."
  ([mjd-utc] mjd-utc)
  ([mjd-utc dut1] (shift mjd-utc dut1)))

;; ------------------------------------------------------ coordinate times

;; TCG and TCB are the coordinate times of the geocentric and barycentric
;; reference systems, the times their metrics are written in; TT and TDB
;; are them rescaled to tick with clocks on the geoid. The scalings are
;; defined constants: TT = TCG - L_G (TCG - T0) (IAU 2000 Resolution
;; B1.9) and TDB = TCB - L_B (TCB - T0) + TDB0 (IAU 2006 Resolution B3),
;; T0 1977 January 1, 0h TAI, when all four agreed but for TDB0. The rates
;; are why TCG gains a second on TT every 45 years and TCB one on TDB
;; every two.

(def ^:private L-G 6.969290134e-10)
(def ^:private L-B 1.550519768e-8)
(def ^:private TDB0 (/ -6.55e-5 86400.0))
(def ^:private T0 43144.0003725)

(defn tt->tcg
  "Geocentric Coordinate Time from TT, MJD."
  [mjd-tt]
  (+ mjd-tt (* (/ L-G (- 1.0 L-G)) (- mjd-tt T0))))

(defn tcg->tt
  "TT from Geocentric Coordinate Time, MJD."
  [mjd-tcg]
  (- mjd-tcg (* L-G (- mjd-tcg T0))))

(defn tdb->tcb
  "Barycentric Coordinate Time from TDB, MJD."
  [mjd-tdb]
  (+ T0 (/ (- mjd-tdb T0 TDB0) (- 1.0 L-B))))

(defn tcb->tdb
  "TDB from Barycentric Coordinate Time, MJD."
  [mjd-tcb]
  (+ (- mjd-tcb (* L-B (- mjd-tcb T0))) TDB0))

(defn tt->tdb
  "Barycentric Dynamical Time from Terrestrial Time, MJD.

  The two differ by at most 1.7 ms, a relativistic effect: a clock on Earth
  runs at a varying rate relative to one at the solar system barycenter,
  because the Earth's distance from the Sun and its speed both vary over the
  year. The dominant term is annual and follows Earth's mean anomaly."
  [mjd-tt]
  (let [g (* c/degrees (+ 357.53 (* 0.9856003 (- (mjd->jd mjd-tt) 2451545.0))))]
    (shift mjd-tt (+ (* 0.001658 (math/sin g))
                     (* 0.000014 (math/sin (* 2.0 g)))))))

;; ------------------------------------------------------------ ΔT (Meeus 10)
;;
;; TT - UT1, the gap between uniform time and the Earth's rotation. Since
;; 1972 it is what the leap-second table above plus dUT1 gives exactly;
;; before that, and after the last published value, it can only be
;; estimated -- from eclipse and occultation timings going back, and by
;; extrapolation going forward.

(def delta-t-table
  "ΔT in seconds every second year from 1620 to 2010 -- Meeus's table 10.A,
  continued past 1998 with the observed values."
  {:first 1620.0 :last 2010.0
   :values [121.0 112.0 103.0 95.0 88.0 82.0 77.0 72.0 68.0 63.0
            60.0 56.0 53.0 51.0 48.0 46.0 44.0 42.0 40.0 38.0
            35.0 33.0 31.0 29.0 26.0 24.0 22.0 20.0 18.0 16.0
            14.0 12.0 11.0 10.0 9.0 8.0 7.0 7.0 7.0 7.0
            7.0 7.0 8.0 8.0 9.0 9.0 9.0 9.0 9.0 10.0
            10.0 10.0 10.0 10.0 10.0 10.0 10.0 11.0 11.0 11.0
            11.0 11.0 12.0 12.0 12.0 12.0 13.0 13.0 13.0 14.0
            14.0 14.0 14.0 15.0 15.0 15.0 15.0 15.0 16.0 16.0
            16.0 16.0 16.0 16.0 16.0 16.0 15.0 15.0 14.0 13.0
            13.1 12.5 12.2 12.0 12.0 12.0 12.0 12.0 12.0 11.9
            11.6 11.0 10.2 9.2 8.2 7.1 6.2 5.6 5.4 5.3
            5.4 5.6 5.9 6.2 6.5 6.8 7.1 7.3 7.5 7.6
            7.7 7.3 6.2 5.2 2.7 1.4 -1.2 -2.8 -3.8 -4.8
            -5.5 -5.3 -5.6 -5.7 -5.9 -6.0 -6.3 -6.5 -6.2 -4.7
            -2.8 -0.1 2.6 5.3 7.7 10.4 13.3 16.0 18.2 20.2
            21.1 22.4 23.5 23.8 24.3 24.0 23.9 23.9 23.7 24.0
            24.3 25.3 26.2 27.3 28.2 29.1 30.0 30.7 31.4 32.2
            33.1 34.0 35.0 36.5 38.3 40.2 42.2 44.5 46.5 48.5
            50.5 52.2 53.8 54.9 55.8 56.9 58.3 60.0 61.6 63.0
            63.8 64.3 64.6 64.8 65.5 66.1]})

(defn delta-t
  "ΔT = TT - UT in seconds, for a decimal `year` (see `decimal-year`).

  Table 10.A where it reaches, interpolated to the second difference;
  outside it, the parabolas Meeus quotes -- Morrison and Stephenson's
  before 948, Stephenson and Houlden's after, the latter with a linear
  correction through 2100 that Meeus added to join it to the present. All
  of them are fits to a quantity that varies irregularly, so the
  uncertainty grows fast away from the table: some minutes at the time of
  Christ, and of the order of the correction itself for the late 21st
  century. The parabola does not meet the table's 2010 value, which marks
  where the observations end, and a caller wanting the present day should
  use the leap-second table and published dUT1 instead."
  [year]
  (let [t (/ (- year 2000.0) 100.0)
        {y0 :first yn :last values :values} delta-t-table]
    (cond
      (< year 948.0)   (+ 2177.0 (* 497.0 t) (* 44.1 t t))
      (<= y0 year yn)
      (interp/value (interp/table-3-around year y0 yn values) year)
      :else
      (cond-> (+ 102.0 (* 102.0 t) (* 25.3 t t))
        (and (> year 2000.0) (< year 2100.0)) (+ (* 0.37 (- year 2100.0)))))))

(def delta-t-polynomials
  "Polynomials in Julian centuries since 1900.0 that Meeus gives as
  alternatives to table 10.A over the years each covers, ascending powers,
  seconds."
  {:1800-1997 [-1.02 91.02 265.90 -839.16 -1545.20 3603.62 4385.98
               -6993.23 -6090.04 6298.12 4102.86 -2137.64 -1081.51]
   :1800-1899 [-2.50 228.95 5218.61 56282.84 324011.78 1061660.75
               2087298.89 2513807.78 1818961.41 727058.63 123563.95]
   :1900-1997 [-2.44 87.24 815.20 -2637.80 -18756.33 124906.15
               -303191.19 372919.88 -232424.66 58353.42]})

(defn delta-t-polynomial
  "ΔT in seconds at `mjd` from one of `delta-t-polynomials`."
  [which mjd]
  (interp/horner (/ (- (mjd->jd mjd) 2415020.0) 36525.0)
                 (delta-t-polynomials which)))
