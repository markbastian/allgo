(ns procedurals.astro-test
  (:require [procedurals.astro.constants :as c]
            [procedurals.astro.geopotential :as geo]
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

;; ------------------------------------------------------------ geopotential

(def ^:private test-field
  "A field with terms at every shape the recursion has to handle: zonal
  (m=0), sectorial (m=n) and tesseral (0<m<n)."
  {:GM c/GM-earth :R c/R-earth :normalised? true
   :C {[0 0] 1.0 [2 0] -4.841654e-4 [2 2] 2.43926e-6
       [3 0] 9.5717e-7 [3 1] 2.02929e-6 [3 3] 7.2114e-7 [4 0] 5.3997e-7}
   :S {[2 2] -1.40027e-6 [3 1] 2.4892e-7 [3 3] 1.41437e-6}})

(defn- mag [v] (Math/sqrt (reduce + (map * v v))))
(defn- sub [a b] (mapv - a b))

(def ^:private sample-points
  [[7000.0 0.0 0.0] [5000.0 3000.0 4000.0] [100.0 200.0 7000.0]
   [0.0 0.0 7500.0] [20000.0 5000.0 9000.0] [-6000.0 -3000.0 1500.0]])

(deftest degree-zero-is-a-point-mass
  (let [pm (geo/point-mass c/GM-earth c/R-earth)]
    (doseq [r sample-points]
      (let [d (mag r)
            want (mapv #(* (- (/ c/GM-earth (* d d d))) %) r)]
        (is (< (/ (mag (sub (geo/acceleration pm r 0) want)) (mag want)) 1e-14)
            (str "at " r))))))

(deftest j2-matches-its-closed-form
  ;; The one term with a standard closed form, so it pins the recursion,
  ;; the normalisation and the acceleration formula together.
  (let [only {:GM c/GM-earth :R c/R-earth :normalised? true
              :C {[0 0] 1.0 [2 0] (- (/ geo/J2 (Math/sqrt 5.0)))} :S {}}
        pm   (geo/point-mass c/GM-earth c/R-earth)]
    (doseq [r sample-points]
      (let [d (mag r) [x y z] r
            k (/ (* -1.5 geo/J2 c/GM-earth c/R-earth c/R-earth) (Math/pow d 5))
            s (/ (* 5.0 z z) (* d d))
            want [(* k x (- 1.0 s)) (* k y (- 1.0 s)) (* k z (- 3.0 s))]
            got  (sub (geo/acceleration only r 2) (geo/acceleration pm r 0))]
        (is (< (/ (mag (sub got want)) (mag want)) 1e-12) (str "at " r))))))

(deftest acceleration-is-the-gradient-of-the-potential
  ;; The two are derived independently -- the acceleration is not
  ;; differentiated from the potential -- so agreement checks both. The
  ;; residual must fall as h^2; a mismatch that plateaus instead is a real
  ;; error hiding under finite-difference noise, which is exactly how the
  ;; sign slip in the m>0 branch was found.
  (doseq [r sample-points]
    (let [a  (geo/acceleration test-field r 4)
          at (fn [h] (mapv (fn [i]
                             (/ (- (geo/potential test-field (update r i + h) 4)
                                   (geo/potential test-field (update r i - h) 4))
                                (* 2.0 h)))
                           (range 3)))
          e1 (/ (mag (sub a (at 2.0))) (mag a))
          e2 (/ (mag (sub a (at 1.0))) (mag a))]
      (is (< e2 1e-7) (str "at " r " residual " e2))
      (is (< 3.0 (/ e1 e2) 5.0)
          (str "at " r " the residual must fall as h^2, got ratio " (/ e1 e2))))))

(deftest every-harmonic-shape-is-exercised
  ;; Regression: the sectorial and tesseral terms went through a branch that
  ;; J2 never touches, and were wrong there while every zonal test passed.
  (doseq [[label C S] [["sectorial C22" {[2 2] 2.43926e-6} {}]
                       ["sectorial S22" {} {[2 2] -1.40027e-6}]
                       ["tesseral C31"  {[3 1] 2.0e-6} {}]
                       ["sectorial C33" {[3 3] 7.0e-7} {}]]]
    (testing label
      (let [m  {:GM c/GM-earth :R c/R-earth :normalised? true :C (assoc C [0 0] 1.0) :S S}
            pm (geo/point-mass c/GM-earth c/R-earth)
            r  [5000.0 3000.0 4000.0]
            a  (sub (geo/acceleration m r 4) (geo/acceleration pm r 0))
            n  (mapv (fn [i]
                       (/ (- (- (geo/potential m (update r i + 0.5) 4)
                                (geo/potential pm (update r i + 0.5) 0))
                             (- (geo/potential m (update r i - 0.5) 4)
                                (geo/potential pm (update r i - 0.5) 0)))
                          1.0))
                     (range 3))]
        (is (pos? (mag a)) "the term must actually contribute")
        (is (< (/ (mag (sub a n)) (mag a)) 1e-5) label)))))

(deftest higher-degrees-contribute-and-diminish
  (let [r [5000.0 3000.0 4000.0]
        step (fn [d] (mag (sub (geo/acceleration test-field r d)
                               (geo/acceleration test-field r (dec d)))))]
    (is (pos? (step 2)))
    (is (pos? (step 3)))
    (is (pos? (step 4)))
    (is (> (step 2) (* 100 (step 3))) "J2 dominates everything above it")))

(deftest normalisation-factors-are-right
  (is (close? 1.0 (geo/normalisation-factor 0 0) 1e-14))
  (is (close? (Math/sqrt 5.0) (geo/normalisation-factor 2 0) 1e-14))
  (is (close? (Math/sqrt 7.0) (geo/normalisation-factor 3 0) 1e-14))
  (is (close? 3.0 (geo/normalisation-factor 4 0) 1e-14))
  ;; N_22 = sqrt(2 * 5 * (0!/4!)) = sqrt(10/24)
  (is (close? (Math/sqrt (/ 10.0 24.0)) (geo/normalisation-factor 2 2) 1e-14)))
