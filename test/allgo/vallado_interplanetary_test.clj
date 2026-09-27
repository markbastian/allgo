(ns allgo.vallado-interplanetary-test
  "Vallado chapter 12's patched conics, each flown: departure and flyby
  hyperbolas propagated out to where their asymptotes lie, the Hohmann
  ellipse and the Lambert transfer between planets flown to arrival --
  and the sizes checked against the widely published ones: Earth's sphere
  of influence of some 925,000 km, a Hohmann trip to Mars of 259 days,
  and Mars 2020's launch energy."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.interplanetary :as ip]
            [allgo.astro.planets :as planets]
            [allgo.astro.time :as time]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))

(deftest spheres-of-influence
  (testing "the Earth's, Mars's and Jupiter's, as they are usually quoted"
    (is (< 9.0e5 (ip/sphere-of-influence :earth) 9.5e5))
    (is (< 5.6e5 (ip/sphere-of-influence :mars) 5.9e5))
    (is (< 4.7e7 (ip/sphere-of-influence :jupiter) 4.9e7)))
  (testing "the Hill sphere is the larger, reaching to the Lagrange point"
    (let [a c/AU mu (:earth c/GM-planet)]
      (is (< (ip/sphere-of-influence a mu c/GM-sun) (ip/hill-radius a mu c/GM-sun)))
      (is (< 1.4e6 (ip/hill-radius a mu c/GM-sun) 1.6e6)))))

(deftest hyperbolas
  (let [mu c/GM-earth rp 6678.0 v-inf 3.0
        {:keys [dv e asymptote]} (ip/departure mu rp v-inf)
        vp (+ (math/sqrt (/ mu rp)) dv)
        ;; fly out from periapsis for a year
        [r v] (u/propagate mu [rp 0.0 0.0] [0.0 vp 0.0] (* 365 86400.0))]
    (testing "the burn leaves the excess speed asked for: v^2 - 2mu/r = v-inf^2 all the way out"
      (is (close? (- (v3/dot v v) (/ (* 2 mu) (v3/length r))) (* v-inf v-inf) 1e-10)))
    ;; the velocity, which runs parallel to the asymptote far out -- the
    ;; position stays off the line through the focus by the impact parameter
    (testing "and the body heads out along the asymptote"
      (is (close? e (+ 1.0 (/ (* rp v-inf v-inf) mu)) 1e-12))
      (is (< (abs (- (math/atan2 (second v) (first v)) asymptote)) 1e-4))))
  (testing "capture into a circle takes what departure gave"
    (is (close? (ip/capture c/GM-earth 6678.0 3.0) (:dv (ip/departure c/GM-earth 6678.0 3.0)) 1e-12))))

(deftest flybys
  (let [mu (:jupiter c/GM-planet) rp 200000.0 v-inf 6.0
        d (ip/turn-angle mu rp v-inf)
        ;; a hyperbola with that periapsis and excess speed, flown both ways
        vp (math/sqrt (+ (* v-inf v-inf) (/ (* 2 mu) rp)))
        t (* 3650 86400.0)
        [_ vout] (u/propagate mu [rp 0.0 0.0] [0.0 vp 0.0] t)
        [_ vin] (u/propagate mu [rp 0.0 0.0] [0.0 vp 0.0] (- t))
        bend (math/acos (/ (v3/dot vin vout) (* (v3/length vin) (v3/length vout))))]
    (testing "the turn angle is the bend of the flown hyperbola, far out"
      (is (< (abs (- bend d)) 1e-4)))
    (testing "a gravity assist keeps the excess speed and turns it, gaining 2 v-inf sin(turn/2)"
      (let [v-planet [0.0 13.0 0.0] v-in [3.0 8.0 1.0]
            {:keys [v-out turn dv]} (ip/flyby mu v-planet v-in rp [0.0 0.0 1.0])]
        (is (close? (v3/distance v-out v-planet) (v3/distance v-in v-planet) 1e-12))
        (is (close? dv (* 2 (v3/distance v-in v-planet) (math/sin (* 0.5 turn))) 1e-12))))))

(deftest between-planets
  (testing "Hohmann to Mars: 259 days, and flown it arrives at Mars's distance"
    (let [{:keys [v-inf-depart v-inf-arrive tof]} (ip/hohmann :earth :mars)
          r1 (* c/AU 1.00000261) r2 (* c/AU 1.52371034)
          v1 (+ (math/sqrt (/ c/GM-sun r1)) v-inf-depart)
          [r v] (u/propagate c/GM-sun [r1 0.0 0.0] [0.0 v1 0.0] tof)]
      (is (< 258.0 (/ tof 86400.0) 260.0))
      (is (< 2.9 v-inf-depart 3.0))
      (is (< 2.6 v-inf-arrive 2.7))
      (is (close? (v3/length r) r2 1e-9))
      (is (close? (- (math/sqrt (/ c/GM-sun r2)) (v3/length v)) v-inf-arrive 1e-9))))
  (testing "Mars 2020: launched 2020-07-30, arrived 2021-02-18"
    (let [depart (time/calendar->mjd 2020 7 30) arrive (time/calendar->mjd 2021 2 18)
          {:keys [v1 c3 v-inf-arrive]} (ip/transfer :earth depart :mars arrive)
          [r1] (planets/heliocentric-state :earth depart)
          [r2] (planets/heliocentric-state :mars arrive)
          [r] (u/propagate c/GM-sun r1 v1 (* 86400.0 (- arrive depart)))]
      ;; the mission's launch energy was about 14 km^2/s^2
      (is (< 13.0 c3 16.0))
      (is (< 2.0 (v3/length v-inf-arrive) 4.0))
      (is (< (v3/distance r r2) (* 1e-6 (v3/length r2))) "the transfer arrives at Mars"))))
