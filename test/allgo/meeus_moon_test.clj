(ns allgo.meeus-moon-test
  "Meeus chapters 47 to 54 -- the Moon's position, phase, physical
  ephemeris and events, and eclipses -- against the book's worked
  examples."
  (:require [allgo.astro.eclipse :as ecl]
            [allgo.astro.lunar-events :as ev]
            [allgo.astro.moon :as moon]
            [allgo.astro.solar :as sun]
            [allgo.astro.time :as t]
            [allgo.math :as am]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg mjd->jd]]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private apr-12 (t/calendar->mjd 1992 4 12))

(deftest position
  (testing "47.a: 1992 April 12"
    (let [[l b r] (moon/position apr-12)]
      (is (close? 133.162655 (->deg l) 1e-6))
      (is (close? -3.229126 (->deg b) 1e-6))
      (is (close? 368409.7 r 0.05))
      (is (close? 0.991990 (->deg (moon/horizontal-parallax r)) 1e-6))))
  (testing "p. 344: the mean node at the dates it passes 0 and 180 degrees"
    (doseq [[y m d] [[1913 5 27] [1950 8 17] [2006 6 19] [2099 7 13]]]
      (let [n (moon/mean-node (t/calendar->mjd y m d))]
        (is (< (min n (- (* 2 math/PI) n)) 1e-3) (str y))))
    (doseq [[y m d] [[1922 9 16] [1978 7 19] [2034 5 21]]]
      (is (close? math/PI (moon/mean-node (t/calendar->mjd y m d)) 1e-3) (str y)))))

(deftest illuminated-fraction
  (testing "48.a: from the positions"
    (let [psi (moon/elongation [(deg 134.6885) (deg 13.7684)] [(deg 20.6579) (deg 8.6964)])
          i (moon/phase-angle psi 368410.0 149971520.0)]
      (is (close? 69.0756 (->deg i) 5e-5))
      (is (close? 0.6786 (moon/illuminated-fraction i) 5e-5))))
  (testing "48.a: from the arguments alone"
    (let [i (moon/phase-angle-approximate apr-12)]
      (is (close? 68.88 (->deg i) 0.005))
      (is (close? 0.6801 (moon/illuminated-fraction i) 5e-5))))
  (testing "and from the full theories"
    (let [{:keys [i k]} (moon/phase apr-12)]
      (is (close? 69.0756 (->deg i) 5e-3))
      (is (close? 0.6786 k 1e-4)))))

(deftest physical-ephemeris
  (testing "53.a: 1992 April 12"
    (let [{:keys [l b p]} (moon/libration apr-12)
          [l0 b0] (moon/selenographic-sun apr-12)]
      (is (close? -1.23 (->deg l) 0.005))
      (is (close? 4.20 (->deg b) 0.005))
      (is (close? 15.08 (->deg p) 0.005))
      (is (close? 67.90 (->deg l0) 0.005))
      (is (close? 1.46 (->deg b0) 0.005))
      (is (close? 2.318 (->deg (moon/sun-altitude (deg -20) (deg 9.7) [l0 b0])) 5e-4))))
  (testing "53.b: sunrise at a crater"
    (is (close? (t/calendar->mjd 1992 4 11.8069)
                (moon/sunrise (deg -20) (deg 9.7) (t/calendar->mjd 1992 4 15)) 1e-4))))

(deftest phases
  (testing "49.a: the new moon of 1977 February"
    (is (close? 2443192.94102 (mjd->jd (ev/mean-phase 1977.13 0)) 1e-5))
    (is (close? 2443192.65118 (mjd->jd (ev/moon-phase 1977.13 0)) 1e-5)))
  (testing "49.b: the last quarter of 2044 January"
    (is (close? 2467636.88597 (mjd->jd (ev/mean-phase 2044.04 0.75)) 1e-5))
    (is (close? 2467636.49186 (mjd->jd (ev/moon-phase 2044.04 0.75)) 1e-5)))
  (testing "at each phase the Moon is 0, 90, 180 or 270 degrees east of the Sun"
    ;; to 20 arcseconds, which the Moon covers in 40 seconds of time
    (doseq [q [0 0.25 0.5 0.75]]
      (let [mjd (ev/moon-phase 2024.3 q)
            [lm] (moon/apparent mjd)
            [ls] (sun/apparent mjd)]
        (is (close? 0.0 (am/wrap-angle (- lm ls (* q 2 math/PI))) (* 20 (deg (/ 1 3600.0))))
            (str q))))))

(deftest perigee-and-apogee
  (testing "50.a: the apogee of 1988 October 7"
    (is (close? 2447442.8191 (mjd->jd (ev/mean-perigee-time 1988.75 true)) 1e-4))
    (let [[mjd parallax] (ev/apogee 1988.75)]
      (is (close? 2447442.3543 (mjd->jd mjd) 1e-4))
      (is (close? 3240.679 (->arcsec parallax) 5e-4))))
  (testing "p. 361: perigees against the full theory"
    (doseq [[y m d year] [[1997 12 (+ 9 (/ 16.9 24)) 1997.93] [1998 1 (+ 3 (/ 8.5 24)) 1998.01]
                          [1990 12 (+ 2 (/ 10.8 24)) 1990.92] [1990 12 (+ 30 (/ 23.8 24)) 1991.0]]]
      (is (close? (t/calendar->mjd y m d) (first (ev/perigee year)) 0.1) (str year))))
  (testing "the perigee parallax matches the Moon's distance then"
    (let [[mjd parallax] (ev/perigee 1997.93)
          [_ _ r] (moon/position mjd)]
      (is (close? (->arcsec (moon/horizontal-parallax r)) (->arcsec parallax) 0.1)))))

(deftest nodes-and-declinations
  (testing "51.a: the ascending node of 1987 May 23"
    (is (close? 2446938.76803 (mjd->jd (ev/node-passage 1987.37)) 1e-5)))
  (testing "52.a, 52.b, 52.c: greatest declinations"
    (doseq [[year south? jd dec] [[1988.95 false 2447518.3346 28.1562]
                                  [2049.3 true 2469553.0834 -22.1384]
                                  [-3.8 false 1719672.1412 28.9739]]]
      (let [[mjd d] (ev/greatest-declination year south?)]
        (is (close? jd (mjd->jd mjd) 1e-4) (str year))
        (is (close? dec (->deg d) 1e-4) (str year))))))

(deftest eclipses
  (testing "54.a: the partial solar eclipse of 1993 May 21"
    (let [{:keys [type central? mjd gamma u penumbra magnitude]} (ecl/solar 1993.38)]
      (is (= :partial type))
      (is (not central?))
      (is (close? 2449129.0978 (mjd->jd mjd) 1e-4))
      (is (close? 1.1348 gamma 1e-4))
      (is (close? 0.0097 u 1e-4))
      (is (close? 0.5558 penumbra 1e-4))
      (is (close? 0.740 magnitude 5e-4))))
  (testing "54.b: the total solar eclipse of 2009 July 22"
    (let [{:keys [type central? mjd gamma u]} (ecl/solar 2009.56)]
      (is (= :total type))
      (is central?)
      (is (close? 2455034.6088 (mjd->jd mjd) 1e-4))
      (is (close? 0.0695 gamma 1e-4))
      (is (close? -0.0157 u 1e-4))))
  (testing "54.c: the penumbral lunar eclipse of 1973 June 15"
    (let [{:keys [type mjd gamma rho magnitude semiduration]} (ecl/lunar 1973.46)]
      (is (= :penumbral type))
      (is (close? 0.4625 magnitude 1e-4))
      (is (close? 2441849.3687 (mjd->jd mjd) 1e-4))
      (is (close? -1.3249 gamma 1e-4))
      (is (close? 1.3045 rho 1e-4))
      (is (close? 101 (* 1440 (:penumbral semiduration)) 0.5))))
  (testing "54.d: the total lunar eclipse of 1997 September 16"
    (let [{:keys [type mjd gamma rho sigma magnitude semiduration]} (ecl/lunar 1997.7)]
      (is (= :total type))
      (is (close? 1.1868 magnitude 1e-4))
      (is (close? 2450708.2835 (mjd->jd mjd) 1e-4))
      (is (close? -0.3791 gamma 1e-4))
      (is (close? 0.7534 sigma 1e-4))
      (is (close? 1.2717 rho 1e-4))
      (is (close? 30 (* 1440 (:total semiduration)) 0.5))
      (is (close? 98 (* 1440 (:partial semiduration)) 0.5))
      (is (close? 153 (* 1440 (:penumbral semiduration)) 0.5))))
  (testing "a year has two to five solar eclipses and a few lunar ones"
    (let [n (count (ecl/solar-eclipses 2024.0 2024.99))]
      (is (<= 2 n 5)))
    (is (= [:total :total] (mapv :type (ecl/lunar-eclipses 2025.0 2025.99))))))

(def ^:private erfa-moon98
  "ERFA's eraMoon98 (BSD; a separate implementation of Meeus chapter 47,
  compiled locally) at dates from 1900 to 2100: MJD (TT), and the ecliptic
  longitude and latitude of date (rad) and distance (km) before its
  rotation into the GCRS."
  [[15020.0 -1.528690662834109 0.019344726755660 368391.602286300]
   [21000.25 -2.246540590652771 -0.082881678850777 364762.774523616]
   [28123.5 -3.867155862732842 -0.050744422807410 400113.436639407]
   [33282.0 -5.211349976951466 0.065996586303254 399601.725356845]
   [40587.75 -2.788867082698899 -0.052519018690097 387811.387772077]
   [45000.1 -5.906128698529827 -0.091834499635346 377311.920170246]
   [51544.5 3.897650395562911 0.090255863898317 402444.812387244]
   [55197.0 1.801660086238916 0.012648109200478 359367.306516674]
   [58849.33 6.109711764070158 -0.087595348743545 404229.512299345]
   [62502.0 4.160953377008764 -0.042066535789128 364522.473705170]
   [70000.9 0.926089860145161 -0.006238239667531 382613.268358961]
   [88069.0 2.747152738199672 0.019061726155735 371715.470667818]])

(deftest against-erfa
  (testing "chapter 47's series against ERFA's implementation of it, 1900-2100: identical but for ERFA's
            mean-longitude constant, Simon et al.'s 218.31665436 degrees where Meeus's 218.3164477 includes
            the -0.70 arcsecond light-time term"
    (let [offset (deg (- 218.31665436 218.3164477))]
      (doseq [[mjd l b r] erfa-moon98]
        (let [[l' b' r'] (moon/position mjd)]
          (is (< (abs (->arcsec (am/wrap-angle (- (+ l' offset) l)))) 1e-4) (str mjd))
          (is (< (abs (->arcsec (- b' b))) 1e-3) (str mjd))
          (is (< (abs (- r' r)) 1e-5) (str mjd)))))))
