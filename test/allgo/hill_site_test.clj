(ns allgo.hill-site-test
  "Vallado's HILL2ECI and ECI2HILL (algorithms 48 and 49) and SITE-TRACK
  (algorithm 50), checked by the physics they must agree with: the
  Clohessy-Wiltshire solution for nearby motion, and the Earth's rotation
  for a target that hangs still over a site."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.maneuvers :as maneuvers]
            [allgo.astro.reduction :as reduction]
            [allgo.astro.states :as states]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(< (abs (- %1 %2)) tol) (flatten a) (flatten b))))

(def ^:private a 6778.137)
(def ^:private target
  (let [v (math/sqrt (/ c/GM-earth a)) i (math/to-radians 51.6)]
    [[a 0.0 0.0] [0.0 (* v (math/cos i)) (* v (math/sin i))]]))

(deftest round-trip
  (let [rel [[0.3 -2.5 0.8] [0.001 -0.0004 0.0007]]]
    (is (close? (states/eci->hill target (states/hill->eci target rel)) rel 1e-12))))

(deftest co-orbiting
  (testing "a chaser on the target's own circular orbit, 2 degrees behind,
            is still in the rotating frame"
    (let [th (math/to-radians -2.0)
          [r v] target
          [[_ rel-v]] [(states/eci->hill target (u/propagate c/GM-earth r v (/ th (math/sqrt (/ c/GM-earth (* a a a))))))]]
      (is (every? #(< (abs %) 1e-12) rel-v)))))

(deftest clohessy-wiltshire
  (testing "a nearby chaser flown by two-body motion follows the CW solution"
    (let [r0 [0.1 -0.5 0.2] v0 [0.001 -0.0002 0.0005]
          t 3000.0
          [rt vt] target
          [ri vi] (states/hill->eci target [r0 v0])
          target' (u/propagate c/GM-earth rt vt t)
          chaser' (u/propagate c/GM-earth ri vi t)
          [rho drho] (states/eci->hill target' chaser')
          [cw-r cw-v] (maneuvers/hill a r0 v0 t)
          ;; CW is linear: its error grows as the separation squared over
          ;; a, 2.6 m at the 4.2 km the chaser has drifted to
          sep (v3/length rho)
          tol (/ (* sep sep) a)]
      (is (close? rho cw-r tol))
      (is (close? drho cw-v (* tol (math/sqrt (/ c/GM-earth (* a a a)))))))))

(deftest site-track
  (let [lat (math/to-radians 39.007) lon (math/to-radians -104.883) alt 2.187
        mjd-tt 53157.5 mjd-ut1 53157.4996
        razel [604.68 (math/to-radians 205.6) (math/to-radians 30.7) 2.08 (math/to-radians 0.15) (math/to-radians 0.17)]]
    (testing "back through the Earth-fixed frame to what the site saw"
      (let [s (states/site-track razel lat lon alt mjd-tt mjd-ut1)
            back (states/ecef->razel (reduction/eci->ecef s mjd-tt mjd-ut1) lat lon alt)
            ;; azimuth comes back wrapped to within half a turn
            back (update back 1 #(mod % (* 2.0 math/PI)))]
        (is (close? back razel 1e-9))))
    (testing "a target hanging still over the site moves with the Earth's
              turning alone, at the IERS's nominal rate the reduction uses"
      (let [still (assoc razel 3 0.0 4 0.0 5 0.0)
            [r v] (states/site-track still lat lon alt mjd-tt mjd-ut1)
            [re] (reduction/eci->ecef [r v] mjd-tt mjd-ut1)]
        (is (< (abs (- (v3/length v) (* 7.292115146706979e-5 (math/hypot (re 0) (re 1))))) 1e-12))))))
