(ns allgo.meeus-support
  "Shared helpers for the tests that check Meeus's worked examples: his
  numbers are printed in sexagesimal, and these read them."
  (:require [allgo.astro.constants :as c]))

(defn close? [a b tol] (< (abs (double (- a b))) tol))

(defn deg [d] (* d c/degrees))

(defn dms
  "Degrees, minutes and seconds to radians; the sign of the first nonzero
  part is the sign of the whole, as Meeus prints -0 17 56.9."
  ([d m s] (dms 1 d m s))
  ([sign d m s] (* sign c/degrees (+ (abs d) (/ m 60.0) (/ s 3600.0)))))

(defn hms
  "Hours, minutes and seconds of right ascension to radians."
  [h m s]
  (* 15.0 c/degrees (+ h (/ m 60.0) (/ s 3600.0))))

(defn ->deg [r] (/ r c/degrees))
(defn ->arcsec [r] (/ r c/arcsec))
(defn ->hours [r] (/ r c/degrees 15.0))

(def arcsec c/arcsec)

(defn jd->mjd [jd] (- jd c/jd-mjd-offset))
(defn mjd->jd [mjd] (+ mjd c/jd-mjd-offset))
