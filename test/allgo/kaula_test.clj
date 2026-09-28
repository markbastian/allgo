(ns allgo.kaula-test
  "Kaula's functions against their closed forms, and his expansion against
  the geopotential evaluated where the satellite is: each harmonic's
  potential, summed over p and q in the elements, the same as the
  spherical harmonic at the satellite's position."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kaula :as kaula]
            [allgo.astro.kepler :as kep]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(def ^:private R c/R-earth)

(deftest inclination-functions
  (testing "Kaula's Table 1 for degree 2"
    (doseq [i [0.3 1.0 1.7 2.9]]
      (let [s (math/sin i) co (math/cos i)
            F (fn [m p] (kaula/inclination-function 2 m p i))]
        (is (< (abs (- (F 0 0) (* -0.375 s s))) 1e-14))
        (is (< (abs (- (F 0 1) (- (* 0.75 s s) 0.5))) 1e-14))
        (is (< (abs (- (F 0 2) (* -0.375 s s))) 1e-14))
        (is (< (abs (- (F 1 0) (* 0.75 s (+ 1.0 co)))) 1e-14))
        (is (< (abs (- (F 1 1) (* -1.5 s co))) 1e-14))
        (is (< (abs (- (F 1 2) (* -0.75 s (- 1.0 co)))) 1e-14))
        (is (< (abs (- (F 2 0) (* 0.75 (math/pow (+ 1.0 co) 2)))) 1e-14))
        (is (< (abs (- (F 2 1) (* 1.5 s s))) 1e-14))
        (is (< (abs (- (F 2 2) (* 0.75 (math/pow (- 1.0 co) 2)))) 1e-14))))))

(deftest eccentricity-functions
  (testing "closed forms: G_lp0 with l = 2p is the mean of (a/r)^(l+1) -- (1 - e^2)^-3/2 for l = 2,
            (1 + 3e^2/2)(1 - e^2)^-7/2 for l = 4 -- and every G is 1 or 0 at e = 0"
    (doseq [e [0.01 0.1 0.4]]
      (is (< (abs (- (kaula/eccentricity-function 2 1 0 e) (math/pow (- 1.0 (* e e)) -1.5))) 1e-12))
      (is (< (abs (- (kaula/eccentricity-function 4 2 0 e) (* (+ 1.0 (* 1.5 e e)) (math/pow (- 1.0 (* e e)) -3.5)))) 1e-12)))
    (is (< (abs (- (kaula/eccentricity-function 3 1 0 1e-9) 1.0)) 1e-8))
    (is (< (abs (kaula/eccentricity-function 3 1 2 1e-9)) 1e-8)))
  (testing "Kaula's series for G_200 = 1 - 5/2 e^2 + 13/16 e^4 and G_20,-1 ~ -e/2"
    (let [e 0.02]
      (is (< (abs (- (kaula/eccentricity-function 2 0 0 e) (+ 1.0 (* -2.5 e e) (* (/ 13.0 16.0) (math/pow e 4))))) 1e-10))
      (is (< (abs (- (kaula/eccentricity-function 2 0 -1 e) (* -0.5 e))) 1e-5)))))

(defn- legendre
  "The unnormalized associated Legendre function P_lm(x), without the
  Condon-Shortley phase, by the standard recurrence."
  [l m x]
  (let [pmm (* (reduce * 1.0 (range 1 (* 2 m) 2)) (math/pow (math/sqrt (- 1.0 (* x x))) m))]
    (if (= l m)
      pmm
      (loop [p0 pmm p1 (* x (inc (* 2 m)) pmm) n (inc m)]
        (if (= n l)
          p1
          (recur p1 (/ (- (* x (inc (* 2 n)) p1) (* (+ n m) p0)) (- (inc n) m)) (inc n)))))))

(defn- direct
  "V_lm at inertial position `r` with the Earth turned by `theta`."
  [l m C S r theta]
  (let [rm (v3/length r)
        sinphi (/ (nth r 2) rm)
        lam (- (math/atan2 (nth r 1) (nth r 0)) theta)]
    (* (/ mu rm) (math/pow (/ R rm) l) (legendre l m sinphi)
       (+ (* C (math/cos (* m lam))) (* S (math/sin (* m lam)))))))

(deftest expansion
  (testing "Kaula's expansion in the elements is the harmonic at the satellite"
    (doseq [[l m C S] [[2 0 -1.08263e-3 0.0] [2 2 1.57e-6 -9.0e-7] [3 1 2.19e-6 2.7e-7] [4 3 5.9e-8 -1.2e-8] [5 5 1.7e-10 -9.7e-10]]
            el [{:a 7000.0 :e 0.0 :i 1.0 :raan 0.4 :argp 0.0 :M 1.3}
                {:a 8000.0 :e 0.05 :i 0.6 :raan 2.1 :argp 1.1 :M 4.0}
                {:a 26560.0 :e 0.2 :i 2.4 :raan 5.0 :argp 3.0 :M 0.7}]]
      (let [theta 1.9
            [r] (kep/elements->state mu el)
            expected (direct l m C S r theta)
            got (kaula/term mu R l m C S el theta)]
        (is (< (abs (- got expected)) (* 1e-9 (+ (abs expected) (* (/ mu (:a el)) (math/pow (/ R (:a el)) l) (math/hypot C S)))))
            (str l " " m " " el))))))

(deftest resonance
  (testing "a geostationary orbit makes the (2,2,0,0) term's argument stand still"
    (let [n (/ (* 2 math/PI) 86164.0905)]
      (is (< (abs (kaula/argument-rate 2 2 0 0 {:argp 0.0 :M n :raan 0.0} n)) 1e-18)))))
