(ns procedurals.astro-test
  (:require [procedurals.astro.constants :as c]
            [procedurals.astro.time :as t]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (< (abs (double (- a b))) tol))

(deftest julian-dates-match-known-epochs
  (doseq [[y m d h want] [[2000 1 1 12.0 2451545.0]      ; J2000.0 by definition
                          [1999 1 1  0.0 2451179.5]
                          [1858 11 17 0.0 2400000.5]     ; MJD zero
                          [2024 2 29  6.0 2460369.75]    ; a leap day
                          [1900 1 1  0.0 2415020.5]
                          [1582 10 15 0.0 2299160.5]]]   ; the Gregorian reform
    (is (close? want (t/mjd->jd (t/calendar->mjd y m d h)) 1e-6)
        (str y "-" m "-" d " " h "h"))))

(deftest mjd-and-jd-are-inverse
  (doseq [mjd [0.0 51544.5 60000.25]]
    (is (close? mjd (- (t/mjd->jd mjd) c/jd-mjd-offset) 1e-12))))

(deftest gmst-at-j2000-is-the-defining-value
  ;; 280.46061837 degrees, or 18h 41m 50.55s, is where the series is anchored.
  (is (close? 280.46061837 (/ (t/gmst c/mjd-J2000) c/degrees) 1e-8)))

(deftest gmst-closes-over-a-sidereal-day
  ;; A sidereal day is 86164.0905 s, not 86400: the Earth turns once relative
  ;; to the stars in less than a solar day, and that gap is what separates
  ;; sidereal time from civil time.
  (let [before (t/gmst c/mjd-J2000)
        after  (t/gmst (+ c/mjd-J2000 (/ 86164.0905 86400.0)))]
    (is (< (abs (- after before)) (* 0.001 c/arcsec))
        "returns to the same angle to within a milliarcsecond")))

(deftest gmst-stays-in-range-and-advances
  (doseq [mjd [40000.0 51544.5 60000.0 70000.123]]
    (let [g (t/gmst mjd)]
      (is (<= 0.0 g) (str "mjd " mjd))
      (is (< g c/two-pi) (str "mjd " mjd))))
  (testing "it advances by roughly one turn per day"
    (let [a (t/gmst 60000.0) b (t/gmst 60001.0)]
      (is (< (abs (- b a)) (* 2.0 c/degrees))
          "a day later, back within a couple of degrees of the same angle"))))

(deftest earth-rotation-rate-agrees-with-the-sidereal-day
  ;; omega must be 2*pi over a sidereal day, not a solar one. Getting this
  ;; wrong tilts every drag calculation by the atmosphere's co-rotation.
  (is (close? c/omega-earth (/ c/two-pi 86164.0905) 1e-11)))

(deftest constants-are-self-consistent
  (testing "solar radiation pressure is the solar constant over c"
    ;; M&G's 4.560e-6 implies a solar constant of 1367 W/m^2, the value in
    ;; use when the book was written. Modern TSI is about 1361, which is a
    ;; real 0.4% difference in the resulting acceleration, not a rounding.
    (is (close? c/solar-pressure (/ 1367.0 (* c/c-light 1000.0)) 2e-9))
    (is (close? 1367.0 (* c/solar-pressure c/c-light 1000.0) 0.5)))
  (testing "the Sun dominates the Earth by the expected factor"
    (is (close? 332946.0 (/ c/GM-sun c/GM-earth) 20.0)))
  (testing "the Moon is about 1/81 of the Earth"
    (is (close? 81.30 (/ c/GM-earth c/GM-moon) 0.02))))
