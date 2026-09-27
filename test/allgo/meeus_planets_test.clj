(ns allgo.meeus-planets-test
  "Meeus chapters 31 and 36 to 46 and 55 -- the planets' orbits, their
  phenomena, Pluto, phases, magnitudes, physical ephemerides and
  satellites -- against the book's worked examples."
  (:require [allgo.astro.calendar :as cal]
            [allgo.astro.constants :as c]
            [allgo.astro.illumination :as il]
            [allgo.astro.jupiter-moons :as jm]
            [allgo.astro.phenomena :as ph]
            [allgo.astro.physical :as phys]
            [allgo.astro.planet-orbits :as po]
            [allgo.astro.precession :as pr]
            [allgo.astro.rotation :as rot]
            [allgo.astro.saturn-moons :as sm]
            [allgo.astro.time :as t]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg dms hms jd->mjd mjd->jd]]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private arcsec (dms 0 0 1))

(deftest mean-elements
  (testing "31.a: Mercury, 2065 June 24"
    (let [mjd (t/calendar->mjd 2065 6 24)
          {:keys [L a e i raan peri]} (po/mean-elements :mercury mjd)]
      (is (close? 203.494701 (->deg L) 1e-6))
      (is (close? 0.387098310 a 1e-9))
      (is (close? 0.20564510 e 1e-8))
      (is (close? 7.006171 (->deg i) 1e-6))
      (is (close? 49.107650 (->deg raan) 1e-6))
      (is (close? 78.475382 (->deg peri) 1e-6))
      (testing "referred to J2000, the orbit agrees with chapter 24's reduction"
        (let [j2000 (po/mean-elements-J2000 :mercury mjd)
              reduced (pr/elements (po/mean-elements :mercury mjd) mjd c/mjd-J2000)]
          (doseq [k [:i :raan :argp]]
            (is (close? (k reduced) (k j2000) 1e-12) (str k)))
          (is (close? (:M (po/mean-elements :mercury mjd))
                      (- (:L j2000) (:peri j2000)) 1e-12)
              "and the mean anomaly does not depend on the frame")))))
  (testing "the Earth's orbit tips out of the J2000 ecliptic by 47'' a century"
    (let [{:keys [i]} (po/mean-elements-J2000 :earth (t/julian-epoch->mjd 2100.0))]
      (is (close? 47.0 (->arcsec i) 0.1)))))

(deftest phenomena
  (testing "36.a: Mercury's inferior conjunction, 1993 November 6 at 3h"
    (is (close? 2449297.645 (mjd->jd (ph/phenomenon :mercury-inferior-conjunction 1993.75)) 5e-4)))
  (testing "36.b: Saturn's conjunction, 2125 August 26"
    (is (close? 2497437.904 (mjd->jd (ph/phenomenon :saturn-conjunction 2125.5)) 5e-4)))
  (testing "36.c: Mercury's western elongation, 1993 November 22"
    (let [[mjd el] (ph/phenomenon :mercury-greatest-west-elongation 1993.9)]
      (is (close? 2449314.14 (mjd->jd mjd) 0.005))
      (is (close? 19.7506 (->deg el) 5e-5))))
  (testing "36.d: Mars stationary, 1997 April 27"
    (is (close? 2450566.255 (mjd->jd (ph/phenomenon :mars-second-station 1997.3)) 5e-4)))
  (testing "p. 255: far from the present, to the hour"
    (doseq [[event mjd hour] [[:mercury-inferior-conjunction (t/calendar->mjd 1631 11 7) 7]
                              [:venus-inferior-conjunction (t/calendar->mjd 1882 12 6) 17]
                              [:mars-opposition (t/calendar->mjd 2729 9 9) 3]
                              [:jupiter-opposition (cal/julian-date->mjd -6 9 15) 7]
                              [:saturn-opposition (cal/julian-date->mjd -6 9 14) 9]
                              [:uranus-opposition (t/calendar->mjd 1780 12 17) 14]
                              [:neptune-opposition (t/calendar->mjd 1846 8 20) 4]]]
      (let [got (ph/phenomenon event (t/mjd->julian-epoch mjd))]
        (is (= hour (long (math/floor (+ 0.5 (* 24 (- got mjd)))))) (str event))))))

(deftest pluto
  (testing "37.a: 1992 October 13"
    (let [[l b r] (po/pluto (jd->mjd 2448908.5))]
      (is (close? 232.74071 (->deg l) 5e-6))
      (is (close? 14.58782 (->deg b) 5e-6))
      (is (close? 29.711111 r 5e-7)))
    (let [[ra dec] (po/pluto-astrometric (jd->mjd 2448908.5))]
      (is (close? (hms 15 31 43.8) ra (* 15 0.1 arcsec)))
      (is (close? (dms -1 4 27 29) dec arcsec))))
  (testing "outside 1885-2099 there is no answer"
    (is (nil? (po/pluto (t/calendar->mjd 1850 1 1))))))

(deftest perihelion-and-aphelion
  (testing "38.a: Venus, 1978 December 31"
    (is (close? 2443873.704 (mjd->jd (po/perihelion :venus 1978.79)) 5e-4)))
  (testing "38.b: Mars, 2032 October 24"
    (is (close? 2463530.456 (mjd->jd (po/aphelion :mars 2032.5)) 5e-4)))
  (testing "p. 270: Jupiter and Saturn, and how far the mean orbit misses"
    (let [[y m dd] (t/mjd->calendar (po/aphelion :jupiter 1981.5))]
      (is (= [1981 7 19] [y m dd])))
    (let [[y m dd] (t/mjd->calendar (first (po/apsis-exact :jupiter 1981.5 true)))]
      (is (= [1981 7 28] [y m dd])))
    (let [[y m dd] (t/mjd->calendar (po/perihelion :saturn 1944.5))]
      (is (= [1944 7 30] [y m dd])))
    ;; The distance is flat at an apsis -- it changes by 1e-7 AU in a day
    ;; -- so the truncation of VSOP87 moves the instant by hours: this
    ;; comes out at 22h on the 7th.
    (let [[mjd r] (po/apsis-exact :saturn 1944.5 false)]
      (is (close? (t/calendar->mjd 1944 9 8.5) mjd 0.6))
      (is (close? 9.0288 r 1e-4))))
  (testing "p. 273: the Earth, with and without the Moon"
    (let [[y m dd] (t/mjd->calendar (po/perihelion :earth-moon 1990.0))]
      (is (= [1990 1 3] [y m dd])))
    (let [[y m dd h] (t/mjd->calendar (po/perihelion :earth 1990.0))]
      (is (= [1990 1 4 16] [y m dd (long (math/floor (+ 0.5 h)))])))))

(deftest node-passages
  (let [perihelion (t/calendar->mjd 1986 2 9.45891)]
    (testing "39.a: Halley"
      (let [[mjd r] (po/node-passage 17.9400782 0.96727426 (deg 111.84644) perihelion)]
        (is (close? (t/calendar->mjd 1985 11 9.16) mjd 0.005))
        (is (close? 1.8045 r 5e-5)))
      (let [[mjd r] (po/node-passage 17.9400782 0.96727426 (deg 111.84644) perihelion true)]
        (is (close? (t/calendar->mjd 1986 3 10.37) mjd 0.005))
        (is (close? 0.8493 r 5e-5)))))
  (testing "39.b: a parabolic comet"
    (let [perihelion (t/calendar->mjd 1989 8 20.291)
          [mjd r] (po/node-passage-parabolic 1.324502 (deg 154.9103) perihelion)
          [y m dd] (t/mjd->calendar mjd)]
      (is (= [1977 9 17] [y m dd]))
      (is (close? 28.07 r 0.005))
      (let [[mjd r] (po/node-passage-parabolic 1.324502 (deg 154.9103) perihelion true)]
        (is (close? (t/calendar->mjd 1989 9 17.636) mjd 5e-4))
        (is (close? 1.3901 r 5e-5)))))
  (testing "39.c: Venus's ascending node, 1978 November 27"
    (let [{:keys [a e argp]} (po/mean-elements :venus (t/calendar->mjd 1979 1 1))
          [mjd] (po/node-passage a e argp (po/perihelion :venus 1979.0))]
      (is (close? (t/calendar->mjd 1978 11 27.409) mjd 5e-4)))))

(deftest illumination
  (testing "41.a: Venus"
    (is (close? 0.29312 (math/cos (il/phase-angle 0.724604 0.910947 0.983824)) 5e-6))
    (is (close? 0.647 (il/illuminated-fraction 0.724604 0.910947 0.983824) 5e-4))
    (is (close? 0.29312 (math/cos (il/phase-angle-from-positions
                                   [(deg 26.10588) (deg -2.62102) 0.724604]
                                   [(deg 88.35704) 0.0 0.983824] 0.910947))
                5e-6)))
  (testing "41.b: without an ephemeris"
    (is (close? 0.640 (il/venus-illuminated-fraction (jd->mjd 2448976.5)) 5e-4)))
  (testing "41.c, 41.d: magnitudes"
    (is (close? -3.8 (il/magnitude :venus 0.724604 0.910947 (deg 72.96)) 0.05))
    (is (close? 0.9 (il/magnitude :saturn 9.867882 10.464606 0.0 {:b (deg 16.442) :du (deg 4.198)})
                0.05)))
  (testing "55: semidiameters scale inversely with distance"
    (is (close? (* 2 (il/semidiameter :jupiter 4.0)) (il/semidiameter :jupiter 2.0) 1e-15))
    (is (close? (il/semidiameters :saturn) (il/saturn-polar-semidiameter 1.0 (deg 90)) 1e-15))))

(deftest physical-ephemerides
  (testing "42.a: Mars"
    (let [{:keys [de ds omega p q d k defect]} (phys/mars (jd->mjd 2448935.500683))]
      (is (close? 12.44 (->deg de) 0.005))
      (is (close? -2.76 (->deg ds) 0.005))
      (is (close? 111.55 (->deg omega) 0.005))
      (is (close? 347.64 (->deg p) 0.005))
      (is (close? 279.91 (->deg q) 0.005))
      (is (close? 10.75 (->arcsec d) 0.005))
      (is (close? 0.9012 k 5e-5))
      (is (close? 1.06 (->arcsec defect) 0.005))))
  (testing "43.a: Jupiter"
    (let [{:keys [ds de omega1 omega2 p]} (phys/jupiter (jd->mjd 2448972.50068))]
      (is (close? -2.20 (->deg ds) 0.005))
      (is (close? -2.48 (->deg de) 0.005))
      (is (close? 268.06 (->deg omega1) 0.005))
      (is (close? 72.74 (->deg omega2) 0.005))
      (is (close? 24.80 (->deg p) 0.005))))
  (testing "43.b: Jupiter from mean orbits"
    (let [{:keys [ds de omega1 omega2]} (phys/jupiter-low-precision (jd->mjd 2448972.50068))]
      (is (close? -2.194 (->deg ds) 5e-4))
      (is (close? -2.50 (->deg de) 0.005))
      (is (close? 268.12 (->deg omega1) 0.005))
      (is (close? 72.79 (->deg omega2) 0.005))))
  (testing "45.a: Saturn's ring"
    (let [{:keys [b b' du p a b-axis]} (phys/saturn-ring (jd->mjd 2448972.50068))]
      (is (close? 16.442 (->deg b) 5e-4))
      (is (close? 14.679 (->deg b') 5e-4))
      (is (close? 4.198 (->deg du) 5e-4))
      (is (close? 6.741 (->deg p) 5e-4))
      (is (close? 35.87 (->arcsec a) 0.005))
      (is (close? 10.15 (->arcsec b-axis) 0.005)))))

(deftest satellites
  (testing "44.a: the Galilean moons, the quick way"
    (doseq [[[x y] [x' y']] (map vector (jm/positions-low-precision (jd->mjd 2448972.50068))
                                 [[-3.44 0.21] [7.44 0.25] [1.24 0.65] [7.08 1.10]])]
      (is (close? x' x 0.005))
      (is (close? y' y 0.005))))
  (testing "44.b: and by E5"
    (doseq [[[x y] [x' y']] (map vector (jm/positions (jd->mjd 2448972.50068))
                                 [[-3.4503 0.2137] [7.4418 0.2752] [1.2010 0.5900] [7.0720 1.0290]])]
      (is (close? x' x 2e-4))
      (is (close? y' y 2e-4))))
  (testing "exercise p. 314: Ganymede and Callisto cross the center"
    (let [mjd (t/calendar->mjd 1988 11 23)
          tt  (+ mjd (/ (t/delta-t (t/decimal-year mjd)) 86400.0))
          [_ _ [x3 y3] _] (jm/positions (+ tt (/ (+ 7 (/ 28 60.0)) 24.0)))
          [_ _ _ [x4 y4]] (jm/positions (+ tt (/ (+ 5 (/ 15 60.0)) 24.0)))]
      (is (close? 0.0032 x3 2e-4)) (is (close? -0.8042 y3 2e-4))
      (is (close? 0.0002 x4 2e-4)) (is (close? 1.3990 y4 2e-4))))
  (testing "46.a: Saturn's eight"
    (let [pos (sm/positions (jd->mjd 2451439.50074))]
      (doseq [[moon [x' y']] (map vector sm/names
                                  [[3.102 -0.204] [3.823 0.318] [4.027 -1.061] [-5.365 -1.148]
                                   [-0.972 -3.136] [14.568 4.738] [-18.001 -5.328] [-48.760 4.137]])]
        (let [[x y] (pos moon)]
          (is (close? x' x 0.002) (str moon))
          (is (close? y' y 0.002) (str moon)))))))

(defn- angle-to-plane
  "Angle of vector `v` out of the plane whose pole is the unit vector `n`."
  [v n]
  (math/asin (/ (reduce + (map * v n)) (math/sqrt (reduce + (map * v v))))))

(deftest satellites-in-three-dimensions
  (testing "the Galilean moons circle in Jupiter's equator, as the IAU has its pole"
    (doseq [mjd [51544.5 55000.0 60000.0]]
      (let [pole (rot/pole :jupiter mjd)]
        (doseq [[moon r] (map vector jm/names (jm/positions-3d mjd))]
          (is (< (abs (angle-to-plane r pole)) (deg 0.6)) (str moon " at " mjd))))))
  (testing "at their own distances"
    (doseq [[r lo hi] (map vector (jm/positions-3d 60000.0) [5.8 9.3 14.9 26.2] [6.0 9.5 15.1 26.6])]
      (is (< lo (math/sqrt (reduce + (map * r r))) hi))))
  (testing "Saturn's inner seven circle in its equator; Iapetus is tilted to it"
    (doseq [mjd [51544.5 55000.0 60000.0]]
      (let [pole (rot/pole :saturn mjd)
            pos (sm/positions-3d mjd)]
        (doseq [moon [:mimas :enceladus :tethys :dione :rhea :titan :hyperion]]
          (is (< (abs (angle-to-plane (pos moon) pole)) (deg 2.0)) (str moon " at " mjd)))
        (is (< (abs (angle-to-plane (pos :iapetus) pole)) (deg 20.0)) (str "iapetus at " mjd))))))
