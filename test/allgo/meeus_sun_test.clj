(ns allgo.meeus-sun-test
  "Meeus chapters 25 to 29 and 32 -- VSOP87, the Sun, the seasons and the
  equation of time -- against the book's worked examples.

  Where Meeus prints both his truncated-series result and the full VSOP87
  one, the full one is the target: the truncation here is not his term
  for term, so it agrees with the complete theory to about an arcsecond
  rather than with his rounding of it."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.moon :as moon]
            [allgo.astro.planets :as pl]
            [allgo.astro.solar :as sun]
            [allgo.astro.time :as t]
            [allgo.astro.vsop87 :as vsop]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg dms hms jd->mjd]]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private arcsec (dms 0 0 1))

(def ^:private oct-13 (t/calendar->mjd 1992 10 13))

(deftest vsop87
  (testing "32.a: Venus, 1992 December 20"
    (let [[l b r] (vsop/heliocentric :venus (t/calendar->mjd 1992 12 20))]
      (is (close? 26.11412 (->deg l) 1e-4))
      (is (close? -2.62060 (->deg b) 1e-4))
      (is (close? 0.724602 r 2e-6))))
  (testing "33.a: the FK5 correction"
    (let [l (deg 313.07689) b (deg -2.08489)
          [l5 b5] (vsop/->fk5 [l b] (jd->mjd 2448976.5))]
      (is (close? -0.09027 (->arcsec (- l5 l)) 1e-5))
      (is (close? 0.05535 (->arcsec (- b5 b)) 1e-5))))
  (testing "every planet stays near its own distance"
    (doseq [[p lo hi] [[:mercury 0.30 0.47] [:venus 0.71 0.73] [:earth 0.98 1.02]
                       [:mars 1.38 1.67] [:jupiter 4.9 5.5] [:saturn 9.0 10.1]
                       [:uranus 18.2 20.1] [:neptune 29.7 30.4]]
            mjd [40000.0 51544.5 62000.0]]
      (let [[_ _ r] (vsop/heliocentric p mjd)]
        (is (< lo r hi) (str p " at " mjd))))))

(deftest low-accuracy-sun
  (testing "25.a: 1992 October 13"
    (is (close? -2241.00603 (->deg (sun/mean-anomaly oct-13)) 1e-5))
    (is (close? 0.016711668 (sun/eccentricity oct-13) 1e-9))
    (is (close? 199.90987 (->deg (sun/true-longitude oct-13)) 1e-5))
    (is (close? 0.99766 (sun/radius oct-13) 1e-5))
    (is (close? (dms 199 54 32) (sun/apparent-longitude-low oct-13) (* 0.5 arcsec)))
    (let [[ra dec] (sun/equatorial-low oct-13 true)]
      (is (close? (hms 13 13 31.4) ra (* 15 0.05 arcsec)))
      (is (close? (dms -1 7 47 6) dec (* 0.5 arcsec))))))

(deftest high-accuracy-sun
  (testing "25.b, full VSOP87 values on p. 169"
    (let [[ra dec] (sun/apparent-equatorial oct-13)]
      (is (close? (hms 13 13 30.749) ra (* 15 0.07 arcsec)))
      (is (close? (dms -1 7 47 1.74) dec arcsec))))
  (testing "26.a: rectangular, equinox of date"
    (let [[x y z] (sun/rectangular oct-13)]
      (is (close? -0.9379963 x 5e-6))
      (is (close? -0.3116537 y 5e-6))
      (is (close? -0.1351207 z 5e-6))))
  (testing "26.b: J2000, B1950 and 2044"
    (doseq [[pos want] [[(sun/rectangular-J2000 oct-13) [-0.93739707 -0.31316724 -0.13577841]]
                        [(sun/rectangular-at oct-13 (t/besselian-epoch->mjd 1950.0))
                         [-0.94148805 -0.30266488 -0.13121349]]
                        [(sun/rectangular-at oct-13 (t/julian-epoch->mjd 2044.0))
                         [-0.93368100 -0.32237347 -0.13977803]]]]
      (is (every? true? (map #(close? %1 %2 5e-6) want pos)) (str want)))))

(deftest seasons
  (testing "27.a: the June solstice of 1962"
    (is (close? 2437837.39245 (+ 2400000.5 (sun/season 1962 :june)) 1e-5)))
  (testing "27.b: and by iterating on the Sun, 21h24m42s"
    (is (close? (+ (t/calendar->mjd 1962 6 21) (/ (+ 21 (/ 24 60.0) (/ 42 3600.0)) 24))
                (sun/season-exact 1962 :june)
                (/ 5 86400.0))))
  (testing "table 27.E: the mean instants agree with VSOP87 to a minute"
    (doseq [[year season month day h m s]
            [[1996 :march 3 20 8 4 7] [2000 :march 3 20 7 36 19] [2005 :march 3 20 12 34 29]
             [1996 :june 6 21 2 24 46] [2003 :june 6 21 19 11 32]
             [1998 :september 9 23 5 38 15] [2004 :september 9 22 16 30 54]
             [1996 :december 12 21 14 6 56] [2000 :december 12 21 13 38 30]]]
      (let [want (t/calendar->mjd year month day (+ h (/ m 60.0) (/ s 3600.0)))]
        (is (close? want (sun/season year season) (/ 1 1440.0)) (str year season))
        (is (close? want (sun/season-exact year season) (/ 10 86400.0)) (str year season " exact"))))))

(deftest equation-of-time
  (testing "28.a: 13m42.6s"
    (is (close? (+ 13 (/ 42.6 60)) (* 4 (->deg (sun/equation-of-time oct-13))) (/ 0.1 60))))
  (testing "28.b: Smart's formula"
    (is (close? 0.0598256 (sun/equation-of-time-smart oct-13) 1e-7))))

(deftest the-suns-disk
  (testing "29.a"
    (let [[P B0 L0] (sun/disk (jd->mjd 2448908.50068))]
      (is (close? 26.27 (->deg P) 0.005))
      (is (close? 5.99 (->deg B0) 0.005))
      (is (close? 238.63 (->deg L0) 0.005))))
  (testing "Carrington rotation 1699"
    (is (close? 2444480.7230 (+ 2400000.5 (sun/carrington-rotation 1699)) 1e-4))))

(deftest into-EME2000
  (testing "the Earth from VSOP87 is the Sun of chapter 26, turned around"
    (doseq [mjd [40000.0 51544.5 60000.0 70000.0]]
      (let [earth (vsop/equatorial-J2000 :earth mjd)
            sun   (sun/rectangular-J2000 mjd)]
        ;; the Sun's is corrected to FK5, a tenth of an arcsecond
        (is (every? true? (map #(close? (- %1) %2 1e-6) earth sun)) (str mjd)))))
  (testing "and agrees with Standish's elements to their arcminutes"
    (doseq [p [:mercury :venus :mars :jupiter :saturn :uranus :neptune]]
      (let [a (vsop/equatorial-J2000 p 60000.0)
            b (mapv #(/ % c/AU) (pl/heliocentric p 60000.0))
            cosang (/ (reduce + (map * a b))
                      (math/sqrt (* (reduce + (map * a a)) (reduce + (map * b b)))))]
        (is (< (math/acos (min 1.0 cosang)) (deg (/ 5 60.0))) (str p)))))
  (testing "the ELP Moon agrees with Montenbruck and Gill's to their few arcminutes"
    (doseq [mjd [51544.5 60000.0 60010.3]]
      (let [a (moon/geocentric-J2000 mjd)
            b (eph/moon mjd)
            cosang (/ (reduce + (map * a b))
                      (math/sqrt (* (reduce + (map * a a)) (reduce + (map * b b)))))]
        (is (< (math/acos (min 1.0 cosang)) (deg 0.2)) (str mjd))))))
