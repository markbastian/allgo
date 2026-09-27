(ns allgo.meeus-basics-test
  "Meeus chapters 3 to 12 -- interpolation, fitting, iteration, the
  calendar, ΔT, the Earth's globe and sidereal time -- against the book's
  worked examples."
  (:require [allgo.astro.calendar :as cal]
            [allgo.astro.frames :as fr]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.time :as t]
            [allgo.meeus-support :refer [->deg ->hours close? deg dms jd->mjd mjd->jd]]
            [allgo.numerics.fit :as fit]
            [allgo.numerics.interpolation :as interp]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- d [dd m s] (+ dd (/ m 60.0) (/ s 3600.0)))

;; ---------------------------------------------------------------- chapter 3

(deftest interpolation-examples
  (testing "3.a: the Moon's distance at 4h21m on the 8th"
    (is (close? 0.876125 (interp/value (interp/table-3 7 9 [0.884226 0.877366 0.870531])
                                       (+ 8 (/ (d 4 21 0) 24.0)))
                1e-6)))
  (testing "3.b: Mars's least distance"
    (let [[x y] (interp/extremum (interp/table-3 12 20 [1.3814294 1.3812213 1.3812453]))]
      (is (close? 17.5864 x 1e-4))
      (is (close? 1.3812030 y 1e-7))))
  (testing "3.c: Mercury's latitude crosses zero"
    (is (close? 26.79873 (interp/zero (interp/table-3 26 28 [(- (d 0 28 13.4)) (d 0 6 46.3)
                                                             (d 0 38 23.2)]))
                1e-5)))
  (testing "3.d: the curved case needs Newton"
    (is (close? -0.720759220056 (interp/zero (interp/table-3 -1 1 [-2 3 2]) true) 1e-12)))
  (testing "3.e: five points"
    (is (close? (d 0 54 13.369)
                (interp/value (interp/table-5 27 29 (map #(apply d 0 %) [[54 36.125] [54 24.606]
                                                                         [54 15.486] [54 8.694]
                                                                         [54 4.133]]))
                              (+ 28 (/ (d 3 20 0) 24.0)))
                (/ 0.001 3600))))
  (testing "exercise p. 30: a zero from five points"
    (is (close? 26.638587
                (interp/zero (interp/table-5 25 29 [(- (d 1 11 21.23)) (- (d 0 28 12.31))
                                                    (d 0 16 7.02) (d 1 1 0.13) (d 1 45 46.33)]))
                1e-6)))
  (testing "3.f: the half-way value"
    (is (close? (d 10 25 40.001)
                (interp/half-way (map #(apply d %) [[10 18 48.732] [10 23 22.835]
                                                    [10 27 57.247] [10 32 31.983]]))
                (/ 0.001 3600))))
  (testing "3.g: Lagrange's polynomial"
    (is (every? true? (map #(close? %1 (* 5 %2) 1e-9) [-87 69 -13 1]
                           (interp/lagrange-coefficients [[1 -6] [3 6] [4 9] [6 15]]))))
    (is (close? 0.5 (interp/lagrange [[29.43 0.4913598528] [30.97 0.5145891926]
                                      [27.69 0.4646875083] [28.11 0.4711658342]
                                      [31.58 0.5236885653] [33.05 0.5453707057]] 30)
                1e-9))))

;; ---------------------------------------------------------------- chapter 4

(deftest fitting-examples
  (testing "4.a: a comet's magnitudes"
    (let [[a b] (fit/linear [[0.2982 10.92] [0.2969 11.01] [0.2918 10.99] [0.2905 10.78]
                             [0.2707 10.87] [0.2574 10.80] [0.2485 10.75] [0.2287 10.14]
                             [0.2238 10.21] [0.2156 9.97] [0.1992 9.69] [0.1948 9.57]
                             [0.1931 9.66] [0.1889 9.63] [0.1781 9.65] [0.1772 9.44]
                             [0.1770 9.44] [0.1755 9.32] [0.1746 9.20]])]
      (is (close? 13.67 a 0.005))
      (is (close? 7.03 b 0.005))))
  (testing "4.b: correlation"
    (let [pts [[73 90.4] [38 125.3] [35 161.8] [42 143.4] [78 52.5] [68 50.8] [74 71.5]
               [42 152.8] [52 131.3] [54 98.5] [39 144.8] [61 78.1] [42 89.5] [49 63.9]
               [50 112.1] [62 82.0] [44 119.8] [39 161.2] [43 208.4] [54 111.6] [44 167.1]
               [37 162.1]]
          [a b] (fit/linear pts)]
      (is (close? -2.49 a 0.005))
      (is (close? 244.18 b 0.005))
      (is (close? -0.767 (fit/correlation pts) 0.0005))))
  (testing "p. 40: a parabola through exact points, and as three functions"
    (let [pts [[-4 -6] [-3 -1] [-2 2] [-1 3] [0 2] [1 -1] [2 -6]]]
      (is (every? true? (map #(close? %1 %2 1e-12) [-1 -2 2] (fit/quadratic pts))))
      (is (every? true? (map #(close? %1 %2 1e-12) [-1 -2 2]
                             (fit/functions-3 pts #(* % %) identity (constantly 1.0)))))))
  (testing "4.c: three sines"
    (let [pts (map (fn [[x y]] [(deg x) y])
                   [[3 0.0433] [20 0.2532] [34 0.3386] [50 0.3560] [75 0.4983] [88 0.7577]
                    [111 1.4585] [129 1.8628] [143 1.8264] [160 1.2431] [183 -0.2043]
                    [200 -1.2431] [218 -1.8422] [230 -1.8726] [248 -1.4889] [269 -0.8372]
                    [290 -0.4377] [303 -0.3640] [320 -0.3508] [344 -0.2126]])]
      (is (every? true? (map #(close? %1 %2 5e-5) [1.2 -0.77 0.39]
                             (fit/functions-3 pts math/sin #(math/sin (* 2 %)) #(math/sin (* 3 %))))))))
  (testing "one function"
    (is (close? 1.016 (fit/function-1 [[0 0] [1 1.2] [2 1.4] [3 1.7] [4 2.1] [5 2.2]] math/sqrt)
                0.0005))))

;; ---------------------------------------------------------------- chapter 5

(deftest fixed-point-and-bisection
  (is (close? (math/sqrt 2.0) (fit/bisect #(- (* % %) 2.0) 0.0 2.0) 1e-15))
  (is (close? 0.7390851332151607 (fit/fixed-point math/cos 1.0) 1e-14))
  (is (nil? (fit/fixed-point #(* 2.0 %) 1.0 1e-15 20)) "diverging gives nil"))

;; ---------------------------------------------------------------- chapter 7

(deftest julian-days
  (testing "7.a: Sputnik"
    (is (close? 2436116.31 (mjd->jd (t/calendar->mjd 1957 10 4.81)) 1e-6)))
  (testing "7.b: a Julian date"
    (is (close? 1842713.0 (mjd->jd (t/calendar->mjd 333 1 27.5)) 1e-6)))
  (testing "p. 62: more dates, both calendars"
    (doseq [[y m dd jd] [[2000 1 1.5 2451545.0] [1987 6 19.5 2446966.0] [1988 1 27 2447187.5]
                         [1600 12 31 2305812.5] [837 4 10.3 2026871.8] [-123 12 31 1676496.5]
                         [-1000 2 29 1355866.5] [-1001 8 17.9 1355671.4] [-4712 1 1.5 0.0]]]
      (is (close? jd (mjd->jd (t/calendar->mjd y m dd)) 1e-6) (str y "-" m "-" dd))))
  (testing "7.c: back again"
    (let [[y m dd] (t/mjd->calendar (jd->mjd 2436116.31))]
      (is (= [1957 10 4] [y m dd])))
    (let [[y m dd h] (cal/mjd->julian-date (jd->mjd 1507900.13))]
      (is (= [-584 5] [y m]))
      (is (close? 28.63 (+ dd (/ (or h 0.0) 24.0)) 0.01))))
  (testing "7.c: Halley's periods"
    (is (= 27689.0 (- (t/calendar->mjd 1986 2 9) (t/calendar->mjd 1910 4 20)))))
  (testing "leap years"
    (is (every? t/julian-leap-year? [900 1236]))
    (is (not-any? t/julian-leap-year? [750 1429]))
    (is (every? t/gregorian-leap-year? [1600 2000 2400]))
    (is (not-any? t/gregorian-leap-year? [1700 1800 1900 2100])))
  (testing "7.e: 1954 June 30 was a Wednesday"
    (is (= 3 (t/day-of-week (jd->mjd 2434923.5)))))
  (testing "7.f, 7.g: day of the year, and back"
    (is (= 318 (t/day-of-year 1978 11 14)))
    (is (= 113 (t/day-of-year 1988 4 22)))
    (is (= [11 14] (t/day-of-year->date 318 false)))
    (is (= [4 22] (t/day-of-year->date 113 true)))))

;; ---------------------------------------------------------------- chapter 8

(deftest easter
  (testing "Gregorian, p. 68"
    (doseq [[y md] [[1991 [3 31]] [1992 [4 19]] [1993 [4 11]] [1954 [4 18]] [2000 [4 23]]
                    [1818 [3 22]] [2285 [3 22]] [1700 [4 11]] [2038 [4 25]]]]
      (is (= md (cal/easter y)) (str y))))
  (testing "Julian, p. 69"
    (doseq [y [179 711 1243]]
      (is (= [4 12] (cal/julian-easter y)) (str y)))))

;; ---------------------------------------------------------------- chapter 9

(deftest jewish-and-moslem-calendars
  (testing "9.a"
    (is (= {:year 5751 :pesach [4 10] :new-year [9 20] :months 12 :days 354}
           (cal/jewish-year 1990))))
  (testing "9.b: 1 Muharram 1421"
    (is (= [2000 3 24] (cal/moslem->julian 1421 1 1)))
    (is (= [2000 4 6] (cal/moslem->gregorian 1421 1 1)))
    (is (not (cal/moslem-leap-year? 1421))))
  (testing "9.c: 1991 August 13"
    (is (= [1991 7 31] (cal/gregorian->julian 1991 8 13)))
    (is (= [1412 2 2] (cal/julian->moslem 1991 7 31)))
    (is (= [1412 2 2] (cal/gregorian->moslem 1991 8 13))))
  (testing "the two calendars round trip"
    (doseq [[y m dd] [[1582 10 15] [1900 2 28] [2024 2 29] [-500 3 1]]]
      (is (= [y m dd] (apply cal/julian->gregorian (cal/gregorian->julian y m dd)))))))

;; --------------------------------------------------------------- chapter 10

(deftest delta-t
  (testing "10.a: 1977 February 18, from the table"
    (is (close? 47.6 (t/delta-t (t/decimal-year (t/calendar->mjd 1977 2 18))) 0.05)))
  (testing "10.a: and from the polynomial"
    (is (close? 47.1 (t/delta-t-polynomial :1900-1997 (t/calendar->mjd 1977 2 18 3.628)) 0.05)))
  (testing "10.b: AD 333"
    (is (close? 6146.0 (t/delta-t 333.1) 0.5))))

;; --------------------------------------------------------------- chapter 11

(deftest the-earths-globe
  ;; Meeus's ellipsoid is IAU 1976 (a = 6378.14 km, f = 1/298.257); this
  ;; package uses WGS-84, which differs by 0.8 m in radius -- a few parts in
  ;; ten million: below the printed precision of most answers, and a couple
  ;; of meters on the radii of 11.b.
  (testing "11.a: Palomar"
    (let [[s c] (gd/parallax-constants (dms 33 21 22) 1.706)]
      (is (close? 0.546861 s 1e-6))
      (is (close? 0.836339 c 1e-6))))
  (testing "p. 83: geodetic minus geocentric latitude"
    (let [lat (dms 45 5 46.36)]
      (is (close? (dms 0 11 32.73) (- lat (gd/geodetic->geocentric-latitude lat)) (dms 0 0 0.05)))))
  (testing "11.b: at latitude 42"
    (is (close? 4747.001 (gd/radius-of-parallel (deg 42)) 0.005))
    (is (close? 6364.033 (gd/radius-of-curvature (deg 42)) 0.005)))
  (testing "11.c: Paris to Washington"
    (is (close? 6181.63 (gd/distance [(dms 48 50 11) (dms 2 20 14)]
                                     [(dms 38 55 17) (dms -1 77 3 56)])
                0.01))))

;; --------------------------------------------------------------- chapter 12

(deftest sidereal-time
  ;; allgo.astro.time/gmst and frames/gast are this chapter already
  (testing "12.a: 1987 April 10, 0h UT"
    (let [mjd (jd->mjd 2446895.5)]
      (is (close? (d 13 10 46.3668) (->hours (t/gmst mjd)) (/ 0.0001 3600)))
      (is (close? (d 13 10 46.1351) (->hours (fr/gast mjd mjd)) (/ 0.0001 3600)))))
  (testing "12.b: and at 19h21m"
    (is (close? (d 8 34 57.0896) (->hours (t/gmst (t/calendar->mjd 1987 4 10 (d 19 21 0))))
                (/ 0.0001 3600)))))

(deftest obliquity-and-nutation
  (testing "22.a: 1987 April 10"
    (let [mjd (jd->mjd 2446895.5)
          [dpsi deps] (fr/nutation-angles mjd)]
      (is (close? -3.788 (/ dpsi (dms 0 0 1)) 0.001))
      (is (close? 9.443 (/ deps (dms 0 0 1)) 0.001))
      (is (close? (d 23 26 27.407) (->deg (fr/mean-obliquity mjd)) (/ 0.001 3600)))
      (is (close? (d 23 26 36.850) (->deg (fr/true-obliquity mjd)) (/ 0.001 3600)))
      (is (close? (fr/mean-obliquity mjd) (fr/mean-obliquity-laskar mjd) (dms 0 0 0.01)))
      (testing "and the four-term approximation is good to half an arcsecond"
        (let [[lp le] (fr/nutation-angles-low-precision mjd)]
          (is (close? dpsi lp (dms 0 0 0.5)))
          (is (close? deps le (dms 0 0 0.1))))))))
