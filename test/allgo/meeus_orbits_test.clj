(ns allgo.meeus-orbits-test
  "Meeus chapters 23, 30 and 33 to 35 -- the apparent place of a star,
  Kepler's equation, elliptic, parabolic and near-parabolic motion --
  against the book's worked examples."
  (:require [allgo.astro.apparent :as ap]
            [allgo.astro.constants :as c]
            [allgo.astro.elliptic :as ell]
            [allgo.astro.kepler :as kep]
            [allgo.astro.precession :as pr]
            [allgo.astro.time :as t]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg dms hms]]
            [clojure.test :refer [deftest is testing]]))

(def ^:private arcsec (dms 0 0 1))

(def ^:private nov-2028 (t/calendar->mjd 2028 11 13.19))

(deftest apparent-place-of-a-star
  (let [theta-persei [(hms 2 44 11.986) (dms 49 13 42.48)]
        pm [(* 15 0.03425 arcsec) (* -0.0895 arcsec)]
        mean [(hms 2 46 11.331) (dms 49 20 54.54)]]
    (testing "23.a: nutation and aberration"
      (let [[da dd] (ap/nutation mean nov-2028)]
        (is (close? 15.843 (->arcsec da) 0.001))
        (is (close? 6.217 (->arcsec dd) 0.001)))
      (let [[da dd] (ap/aberration mean nov-2028)]
        (is (close? 30.045 (->arcsec da) 0.001))
        (is (close? 6.697 (->arcsec dd) 0.001))))
    (testing "23.a: the whole reduction"
      (let [[ra dec] (ap/apparent-place theta-persei c/mjd-J2000 nov-2028 pm)]
        (is (close? (hms 2 46 14.390) ra (* 15 0.002 arcsec)))
        (is (close? (dms 49 21 7.45) dec (* 0.02 arcsec)))))
    (testing "23.b: Ron and Vondrák"
      (let [[da dd] (ap/aberration-ron-vondrak [(hms 2 44 12.9747) (dms 49 13 39.896)] nov-2028)]
        (is (close? 0.000145252 da 1e-9))
        (is (close? 0.000032723 dd 1e-9)))
      (let [[ra dec] (ap/apparent-place-ron-vondrak theta-persei nov-2028 pm)]
        (is (close? (hms 2 46 14.392) ra (* 15 0.002 arcsec)))
        (is (close? (dms 49 21 7.45) dec (* 0.02 arcsec)))))))

(deftest keplers-equation
  (testing "30.a, 30.b: e = 0.1, M = 5 degrees"
    (is (close? 5.554589254 (->deg (kep/kepler-equation (deg 5) 0.1)) 1e-9))
    (is (close? 5.554589254 (->deg (kep/kepler-bisection (deg 5) 0.1)) 1e-9)))
  (testing "p. 205: e = 0.99, M = 0.2, where Newton needs care"
    (is (close? 1.066997365282 (kep/kepler-equation 0.2 0.99) 1e-12))
    (is (close? 1.066997365282 (kep/kepler-bisection 0.2 0.99) 1e-12)))
  (testing "the bisection handles the second half of the orbit"
    (let [E (kep/kepler-bisection 4.0 0.5)]
      (is (close? 4.0 (- E (* 0.5 (Math/sin E))) 1e-13))))
  (testing "p. 207: the first-order formula"
    (is (close? 5.554599 (->deg (kep/kepler-approximate (deg 5) 0.1)) 1e-6))))

(deftest elliptic-motion
  (testing "33.a: Venus, 1992 December 20, full VSOP87 values on p. 227"
    (let [[ra dec] (ell/apparent :venus (t/calendar->mjd 1992 12 20))]
      (is (close? (hms 21 4 41.454) ra (* 15 0.05 arcsec)))
      (is (close? (dms -1 18 53 16.84) dec arcsec))))
  (testing "33.b: comet Encke, 1990 October 6"
    (let [el (pr/elements-B1950->J2000 {:i (deg 11.93911) :raan (deg 334.04096)
                                        :argp (deg 186.24444)})
          orbit (ell/orbit (assoc el :a 2.2091404 :e 0.8502196
                                  :perihelion (t/calendar->mjd 1990 10 28.54502)))
          [ra dec _ psi] (ell/astrometric orbit (t/calendar->mjd 1990 10 6))]
      (is (close? (hms 10 34 14.2) ra (* 15 0.1 arcsec)))
      (is (close? (dms 19 9 31) dec arcsec))
      (is (close? 40.51 (->deg psi) 0.005))))
  (testing "33.c, 33.d: Halley's speeds and the length of its orbit"
    (let [a 17.9400782 e 0.96727426
          v (fn [r a] (kep/vis-viva (* c/GM-sun 1.0) (* r c/AU) (* a c/AU)))]
      (is (close? 41.53 (v 1.0 a) 0.01))
      (is (close? 54.52 (v (* a (- 1 e)) a) 0.01))
      (is (close? 0.91 (v (* a (+ 1 e)) a) 0.01))
      (is (close? 77.06 (kep/ellipse-circumference a e) 0.005))
      (is (close? 77.07 (kep/ellipse-circumference a e true) 0.005)))))

(deftest parabolic-motion
  (testing "34.a"
    (let [[nu r] (kep/parabolic 1.487469 (- (t/calendar->mjd 1998 8 5)
                                            (t/calendar->mjd 1998 4 14.4358)))]
      (is (close? 66.78862 (->deg nu) 1e-5))
      (is (close? 2.133911 r 1e-6)))))

(deftest near-parabolic-motion
  (testing "p. 247"
    (doseq [[q e dt nu r] [[0.921326 1.0 138.4783 102.74426 2.364192]
                           [0.1 0.987 254.9 164.50029 4.063777]
                           [0.123456 0.99997 -30.47 221.91190 0.965053]
                           [3.363943 1.05731 1237.1 109.40598 10.668551]
                           [0.5871018 0.9672746 20.0 52.85331 0.729116]
                           [0.5871018 0.9672746 0.0 0.0 0.5871018]]]
      (let [[nu' r'] (kep/near-parabolic q e dt)]
        (is (close? nu (->deg nu') 1e-5) (str [q e dt]))
        (is (close? r r' 1e-6) (str [q e dt])))))
  (testing "p. 248: where it converges, and where it does not"
    (doseq [[q e dt nu tol] [[0.1 0.9 10 126 1] [0.1 0.987 400 167 1] [0.1 0.999 5000 174 1]
                             [1 0.99999 100000 172.5 0.1] [1 0.99999 17000000 178.68 0.01]]]
      (is (close? nu (->deg (first (kep/near-parabolic q e dt))) tol) (str [q e dt])))
    (doseq [[q e dt] [[0.1 0.9 30] [0.1 0.987 500] [1 0.99999 18000000]]]
      (is (nil? (kep/near-parabolic q e dt)) (str [q e dt])))))
