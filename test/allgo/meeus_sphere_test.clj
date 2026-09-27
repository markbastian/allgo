(ns allgo.meeus-sphere-test
  "Meeus chapters 13 to 24 and 40 -- coordinate transformations, angles on
  the sky, refraction, precession and parallax -- against the book's
  worked examples."
  (:require [allgo.astro.calendar :as cal]
            [allgo.astro.coordinates :as co]
            [allgo.astro.frames :as fr]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.precession :as pr]
            [allgo.astro.time :as t]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg dms hms]]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private eps-J2000 (deg 23.4392911))

(def ^:private arcsec (dms 0 0 1))

(deftest ecliptic-and-equatorial
  (testing "13.a: Pollux"
    (let [[l b] (co/equatorial->ecliptic [(hms 7 45 18.946) (dms 28 1 34.26)] eps-J2000)]
      (is (close? 113.21563 (->deg l) 1e-5))
      (is (close? 6.684170 (->deg b) 1e-6))
      (let [[a d] (co/ecliptic->equatorial [l b] eps-J2000)]
        (is (close? (hms 7 45 18.946) a 1e-12))
        (is (close? (dms 28 1 34.26) d 1e-12))))))

(deftest horizontal
  (testing "13.b: Venus from Washington, 1987 April 10, 19h21m UT"
    ;; Meeus takes the apparent sidereal time 8h34m56.853s. Washington is
    ;; at 77 03 56 west, so -77 03 56 east here, and Meeus's azimuth of
    ;; 68.034 from the south is 248.034 from the north.
    (let [mjd   (t/calendar->mjd 1987 4 10 (+ 19 (/ 21 60.0)))
          theta (fr/gast mjd (+ mjd (/ 56.0 86400)))
          lat   (dms 38 55 17)
          lon   (dms -1 77 3 56)
          [az alt] (co/equatorial->horizontal [(hms 23 9 16.641) (dms -1 6 43 11.61)]
                                              lat lon theta)]
      (is (close? (hms 8 34 56.853) theta (* 15 0.001 arcsec)))
      (is (close? 248.034 (->deg az) 0.0005))
      (is (close? 15.125 (->deg alt) 0.0005))
      (let [[a d] (co/horizontal->equatorial [az alt] lat lon theta)]
        (is (close? (hms 23 9 16.641) a 1e-10))
        (is (close? (dms -1 6 43 11.61) d 1e-10))))))

(deftest galactic
  (testing "exercise p. 96: Nova Serpentis 1978"
    (let [[l b] (co/equatorial->galactic [(hms 17 48 59.74) (dms -1 14 43 8.2)])]
      (is (close? 12.9593 (->deg l) 5e-5))
      (is (close? 6.0463 (->deg b) 5e-5))
      (let [[a d] (co/galactic->equatorial [l b])]
        (is (close? (hms 17 48 59.74) a 1e-10))
        (is (close? (dms -1 14 43 8.2) d 1e-10))))))

(deftest chapter-14
  (testing "the ecliptic on the horizon, p. 100"
    (let [[l1 l2 i] (co/ecliptic-at-horizon (deg 23.44) (deg 51) (hms 5 0 0))]
      (is (close? (dms 169 21 30) l1 (* 0.5 arcsec)))
      (is (close? (dms 349 21 30) l2 (* 0.5 arcsec)))
      (is (close? (dms 61 53 14) i (* 0.5 arcsec)))))
  (testing "the diurnal path at the horizon"
    (is (close? (- (/ math/PI 2) (deg 40)) (co/diurnal-path-at-horizon 0.0 (deg 40)) 1e-15))
    (is (close? (dms 45 31 0) (co/diurnal-path-at-horizon (deg 23.44) (deg 40)) (deg 0.05)))))

(deftest refraction
  (testing "16.a: the setting Sun"
    (let [h0 (deg 0.5)
          R  (co/refraction-from-apparent h0)
          lower (- h0 R)
          upper (+ lower (dms 0 32 0))]
      (is (close? 28.754 (* 60 (->deg R)) 0.0005))
      (is (close? 24.618 (* 60 (->deg (co/refraction upper))) 0.0005))))
  (testing "p. 106: at the zenith Bennett is 0.08'' out, and his correction 0.89''"
    (is (close? -0.08 (->arcsec (co/refraction-from-apparent (/ math/PI 2))) 0.01))
    (is (close? -0.89 (->arcsec (co/refraction-from-apparent (/ math/PI 2) true)) 0.01)))
  (testing "above 15 degrees the series agrees with Bennett, to Bennett's 0.07'"
    (is (close? (co/refraction-from-apparent (deg 30)) (co/refraction-high (deg 30) true)
                (* 4.2 arcsec)))))

(deftest separation
  (let [arcturus [(hms 14 15 39.7) (dms 19 10 57)]
        spica    [(hms 13 25 11.6) (dms -1 11 9 41)]]
    (testing "17.a and 17.b"
      (is (close? (dms 32 47 35) (co/separation arcturus spica) (* 0.5 arcsec)))
      (is (close? (dms 32 47 35) (co/separation-haversine arcturus spica) (* 0.5 arcsec))))
    (testing "exercise p. 110"
      (is (close? (dms 169 58 0) (co/separation [(hms 4 35 55.2) (dms 16 30 33)]
                                                [(hms 16 29 24) (dms -1 26 25 55)])
                  1e-4))))
  (testing "p. 111 and 113: Mercury and Saturn, 1978 September"
    (let [mercury [[(hms 10 29 44.27) (dms 11 2 5.9)] [(hms 10 36 19.63) (dms 10 29 51.7)]
                   [(hms 10 43 1.75) (dms 9 55 16.7)]]
          saturn  [[(hms 10 33 29.64) (dms 10 40 13.2)] [(hms 10 33 57.97) (dms 10 37 33.4)]
                   [(hms 10 34 26.22) (dms 10 34 53.9)]]]
      (is (close? 0.5017 (->deg (co/minimum-separation 13 15 mercury saturn)) 5e-4))
      (let [[_ sep] (co/minimum-separation-rectangular 13 15 mercury saturn)]
        (is (close? 224 (->arcsec sep) 2)))))
  (testing "position angle: due north is 0, due east 90"
    (is (close? 0.0 (co/position-angle [1.0 0.2] [1.0 0.3]) 1e-12))
    (is (close? (deg 90) (co/position-angle [1.0 0.0] [1.01 0.0]) 1e-9))))

(deftest conjunctions
  (testing "18.a: Mercury and Venus, 1991 August"
    (let [venus   [[(hms 10 27 27.175) (dms 4 4 41.83)] [(hms 10 26 32.410) (dms 3 55 54.66)]
                   [(hms 10 25 29.042) (dms 3 48 3.51)] [(hms 10 24 17.191) (dms 3 41 10.25)]
                   [(hms 10 22 57.024) (dms 3 35 16.61)]]
          mercury [[(hms 10 24 30.125) (dms 6 26 32.05)] [(hms 10 25 0.342) (dms 6 10 57.72)]
                   [(hms 10 25 12.515) (dms 5 57 33.08)] [(hms 10 25 6.235) (dms 5 46 27.07)]
                   [(hms 10 24 41.185) (dms 5 37 48.45)]]
          [day dd] (co/conjunction 5 9 venus mercury)]
      (is (close? 7.23797 day 1e-5))
      (is (close? (dms 2 8 22) dd (* 0.5 arcsec))))))

(deftest straight-lines
  (testing "19.a: Mars between Castor and Pollux"
    (let [castor [(deg 113.56833) (deg 31.89756)]
          pollux [(deg 116.25042) (deg 28.03681)]
          mars   (map (fn [a d] [(deg a) (deg d)])
                      [118.98067 119.59396 120.20413 120.81108 121.41475]
                      [21.68417 21.58983 21.49394 21.39653 21.29761])
          t (co/collinear-time castor pollux (t/calendar->mjd 1994 9 29)
                               (t/calendar->mjd 1994 10 3) mars)]
      (is (close? (t/calendar->mjd 1994 10 1.2233) t 1e-4))))
  (let [delta [(hms 5 32 0.40) (dms -1 0 17 56.9)]
        eps   [(hms 5 36 12.81) (dms -1 1 12 7.0)]
        zeta  [(hms 5 40 45.52) (dms -1 1 56 33.3)]]
    (testing "p. 123: Orion's belt"
      (is (close? 172.4830 (->deg (co/collinear-angle delta eps zeta)) 5e-5)))
    (testing "p. 124"
      (is (close? 0.089876 (abs (->deg (co/great-circle-distance zeta delta eps))) 5e-7)))
    (testing "p. 125"
      (let [[psi omega] (co/collinearity delta eps zeta)]
        (is (close? (dms 7 31 0) psi (dms 0 0 30)))
        (is (close? (dms -1 0 5 24) omega (* 0.5 arcsec)))))))

(deftest smallest-circles
  (testing "20.a"
    (let [[dia pair?] (co/smallest-circle [(hms 12 41 8.64) (dms -1 5 37 54.2)]
                                          [(hms 12 52 5.21) (dms -1 4 22 26.2)]
                                          [(hms 12 39 28.11) (dms -1 1 50 3.7)])]
      (is (close? 4.26363 (->deg dia) 5e-6))
      (is (false? pair?))))
  (testing "exercise p. 128"
    (let [[dia pair?] (co/smallest-circle [(hms 9 5 41.44) (dms 18 30 30)]
                                          [(hms 9 9 29) (dms 17 43 56.7)]
                                          [(hms 8 59 47.14) (dms 17 49 36.8)])]
      (is (close? (dms 2 19 0) dia (dms 0 0 30)))
      (is (true? pair?)))))

(deftest precession
  (testing "21.a: Regulus, the quick way"
    (let [pos  [(hms 10 8 22.3) (dms 11 58 2)]
          from (t/julian-epoch->mjd 2000.0)
          to   (t/julian-epoch->mjd 1978.0)
          [ra dec] (pr/equatorial-approximate pos from to [(* 15 -0.0169 arcsec) (* 0.006 arcsec)])]
      (let [[a d] (pr/annual-precession pos from to)]
        (is (close? 3.207 (/ (->arcsec a) 15) 0.001))
        (is (close? -17.71 (->arcsec d) 0.01)))
      (is (close? (hms 10 7 12.1) ra (* 15 0.05 arcsec)))
      (is (close? (dms 12 4 32) dec (* 0.5 arcsec)))))
  (testing "21.b: theta Persei to 2028 November 13.19"
    (let [[ra dec] (pr/equatorial [(hms 2 44 11.986) (dms 49 13 42.48)]
                                  (t/julian-epoch->mjd 2000.0)
                                  (t/calendar->mjd 2028 11 13.19)
                                  [(* 15 0.03425 arcsec) (* -0.0895 arcsec)])]
      (is (close? (hms 2 46 11.331) ra (* 15 0.001 arcsec)))
      (is (close? (dms 49 20 54.54) dec (* 0.01 arcsec)))))
  (testing "exercise p. 136: Polaris, near the pole"
    (doseq [[to ra dec] [[(t/besselian-epoch->mjd 1900.0) (hms 1 22 33.90) (dms 88 46 26.18)]
                         [(t/julian-epoch->mjd 2050.0) (hms 3 48 16.43) (dms 89 27 15.38)]
                         [(t/julian-epoch->mjd 2100.0) (hms 5 53 29.17) (dms 89 32 22.18)]]]
      (let [[a d] (pr/equatorial [(hms 2 31 48.704) (dms 89 15 50.72)]
                                 (t/julian-epoch->mjd 2000.0) to
                                 [(* 15 0.19877 arcsec) (* -0.0152 arcsec)])]
        (is (close? ra a (* 15 0.01 arcsec)))
        (is (close? dec d (* 0.01 arcsec))))))
  (testing "21.c: Venus's ecliptic position back to -214"
    (let [[l b] (pr/ecliptic [(deg 149.48194) (deg 1.76549)]
                             (t/julian-epoch->mjd 2000.0)
                             (cal/julian-date->mjd -214 6 30))]
      (is (close? 118.704 (->deg l) 0.0005))
      (is (close? 1.615 (->deg b) 0.0005))))
  (testing "21.d: Sirius's proper motion in three dimensions"
    (let [pos [(hms 6 45 8.871) (dms -1 16 42 57.99)]
          pm  [(* 15 -0.03847 arcsec) (* -1.2053 arcsec)]
          from (t/julian-epoch->mjd 2000.0)]
      (doseq [[year ra dec] [[1000 (hms 6 45 47.16) (dms -1 16 22 56.0)]
                             [-2000 (hms 6 47 39.91) (dms -1 15 23 30.6)]
                             [-10000 (hms 6 52 25.72) (dms -1 12 50 6.7)]]]
        (let [[a d] (pr/proper-motion-3d pos from (t/julian-epoch->mjd year)
                                         2.64 (/ -7.6 977792) pm)]
          (is (close? ra a (* 15 0.01 arcsec)) (str year))
          (is (close? dec d (* 0.1 arcsec)) (str year)))))))

(deftest reduction-of-elements
  (let [el {:i (deg 47.122) :argp (deg 151.4486) :raan (deg 45.7481)}]
    (testing "24.a: comet Klinkenberg, 1744 to 1950"
      (let [{:keys [i argp raan]} (pr/elements el (t/besselian-epoch->mjd 1744.0)
                                               (t/besselian-epoch->mjd 1950.0))]
        (is (close? 47.1380 (->deg i) 5e-5))
        (is (close? 48.6037 (->deg raan) 5e-5))
        (is (close? 151.4782 (->deg argp) 5e-5)))))
  (let [el {:i (deg 11.93911) :raan (deg 334.04096) :argp (deg 186.24444)}]
    (testing "24.b: B1950 to J2000"
      (let [{:keys [i argp raan]} (pr/elements-B1950->J2000 el)]
        (is (close? 11.94524 (->deg i) 5e-6))
        (is (close? 334.75006 (->deg raan) 5e-6))
        (is (close? 186.23352 (->deg argp) 5e-6))))
    (testing "24.c: FK4 to FK5"
      (let [{:keys [i argp raan]} (pr/elements-FK4->FK5 el)]
        (is (close? 11.94521 (->deg i) 5e-6))
        (is (close? 334.75043 (->deg raan) 5e-6))
        (is (close? 186.23327 (->deg argp) 5e-6))))))

(deftest parallax
  (let [pos   [(deg 339.530208) (deg -15.771083)]
        rho   [0.546861 0.836339]
        lon   (- (hms 7 47 27))                 ; Palomar, west
        mjd   (t/calendar->mjd 2003 8 28 (+ 3 (/ 17 60.0)))
        theta (fr/gast mjd mjd)]
    (testing "40.a: Mars from Palomar"
      (is (close? 23.592 (->arcsec (co/horizontal-parallax 0.37276)) 0.0005))
      (let [[a d] (co/topocentric pos 0.37276 rho lon theta)]
        (is (close? (hms 22 38 8.54) a (* 15 0.01 arcsec)))
        (is (close? (dms -1 15 46 30.0) d (* 0.1 arcsec)))))
    (testing "40.a, by the hour angle"
      (let [H (co/hour-angle (first pos) lon theta)
            [H' d] (co/topocentric-hour-angle H (second pos) 0.37276 rho)]
        (is (close? (- (* 2 math/PI) (hms 4 44 50.28)) H' (* 15 0.01 arcsec)))
        (is (close? (dms -1 15 46 30.0) d (* 0.1 arcsec))))))
  (testing "exercise p. 282: the Moon in ecliptic coordinates"
    (let [rho (gd/parallax-constants (dms 50 5 7.8) 0.0)
          [l b s] (co/topocentric-ecliptic [(dms 181 46 22.5) (dms 2 17 26.2)] (dms 0 16 15.5)
                                           rho (dms 23 28 0.8) (dms 209 46 7.9)
                                           (dms 0 59 27.7))]
      (is (close? (dms 181 48 5.0) l (* 0.1 arcsec)))
      (is (close? (dms 1 29 7.1) b (* 0.1 arcsec)))
      (is (close? (dms 0 16 25.5) s (* 0.1 arcsec))))))
