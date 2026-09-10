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
