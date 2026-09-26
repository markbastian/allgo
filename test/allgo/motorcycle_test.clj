(ns allgo.motorcycle-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.physics.rigid :as rigid]
            [allgo.simulation.motorcycle :as m]
            [clojure.test :refer [deftest is testing]]))

(def ^:private cfg m/defaults)
(def ^:private dt (/ 1.0 240.0))
(def ^:private floor (m/ground 400.0))

(defn- ride [scene controls seconds]
  (m/run cfg scene controls dt (long (/ (double seconds) dt))))

(defn- tel [scene] (m/telemetry cfg (m/bike-pose scene)))

(defn- leaning [speed lean]
  (-> (m/start-pose cfg speed)
      (assoc-in [:base :rot] (q/from-axis-angle [1.0 0.0 0.0] (- (double lean))))))

(deftest sits-at-its-sag-test
  (testing "the springs are preloaded so the bike rolls along with the suspension near zero travel"
    (let [t (tel (ride (m/scene cfg [floor] (m/start-pose cfg 6.0)) {:lean 0.0} 2.0))]
      (is (< (abs (double (:rear-travel t))) 0.01))
      (is (< (abs (double (:front-travel t))) 0.01)))))

(deftest rolls-without-slipping-test
  (testing "at a steady cruise both tyres turn at road speed"
    (let [t (tel (ride (m/scene cfg [floor] (m/start-pose cfg 8.0)) {:throttle 0.05 :lean 0.0} 2.0))]
      (is (< (abs (- (double (:rear-wheel-speed t)) (double (:speed t)))) 0.05))
      (is (< (abs (- (double (:front-wheel-speed t)) (double (:speed t)))) 0.2)))))

(deftest rider-rights-the-bike-test
  (testing "started leaning, a rider holding it upright stands it back up"
    ;; Steering into the fall. A tenth of a radian is a bike already well
    ;; on its way over; the rider has it straight inside two seconds.
    (let [t (tel (ride (m/scene cfg [floor] (leaning 8.0 0.1)) {:throttle 0.05 :lean 0.0} 3.0))]
      (is (:upright? t))
      (is (< (abs (double (:lean t))) 0.01) (str "still leaning " (:lean t))))))

(deftest nobody-steering-test
  (testing "with nobody on the bars the same bike goes over"
    ;; The control for the test above: it is the rider holding the bike
    ;; up, not something in the model that would have held it anyway.
    (let [t (tel (ride (m/scene cfg [floor] (leaning 8.0 0.1)) {:throttle 0.05 :hands? false} 4.0))]
      (is (not (:upright? t))))))

(deftest leans-into-a-turn-test
  (testing "asked for a lean, the rider countersteers into it and holds a steady turn"
    (let [s0 (m/scene cfg [floor] (m/start-pose cfg 8.0))
          early (ride s0 {:throttle 0.08 :lean 0.35} 0.08)
          later (ride early {:throttle 0.08 :lean 0.35} 3.5)
          t (tel later)]
      ;; To lean left the bars first go *right*.
      (is (neg? (double (:steer (tel early)))) "countersteered")
      (is (< (abs (- (double (:lean t)) 0.35)) 0.03) (str "lean " (:lean t)))
      (is (pos? (double (:steer t))) "and then steered into the turn")
      (is (> (double (:heading t)) 0.8) "and has come round"))))

(deftest brakes-test
  (testing "hard on both brakes from 15 m/s, it stops upright in a plausible distance"
    (let [s (ride (m/scene cfg [floor] (m/start-pose cfg 15.0)) {:brake 1.0 :lean 0.0} 3.0)
          t (tel s)
          x (double (first (:pos (:base (m/bike-pose s)))))]
      (is (< (abs (double (:speed t))) 0.05))
      (is (:upright? t))
      ;; v^2 / 2a with a around a g, less than a car's length either way.
      (is (< 9.0 x 15.0) (str "stopped at " x)))))

(deftest suspension-test
  (testing "a kerb compresses the fork, and the bike rides over it and settles"
    (let [kerb (rigid/box {:pos [6.0 0.04 0.0] :size [0.5 0.08 3.0]})
          s0 (m/scene cfg [floor kerb] (m/start-pose cfg 8.0))
          travel (loop [s s0 k 0 lowest 0.0]
                   (if (= k 240)
                     [s lowest]
                     (recur (m/step cfg s {:throttle 0.05 :lean 0.0} dt)
                            (inc k)
                            (min lowest (double (:front-travel (tel s)))))))
          [s lowest] travel
          settled (tel (ride s {:throttle 0.05 :lean 0.0} 1.0))]
      (is (< lowest -0.02) "the fork took the hit")
      (is (> lowest -0.12) "without going through its stop")
      (is (:upright? settled))
      (is (< (abs (double (:front-travel settled))) 0.01)))))

(deftest top-speed-test
  (testing "a steady throttle levels off rather than accelerating for ever"
    ;; A minute near 30 m/s is well over a kilometre, which is further
    ;; than the usual floor goes.
    (let [s1 (ride (m/scene cfg [(m/ground 4000.0)] (m/start-pose cfg 20.0))
                   {:throttle 0.3 :lean 0.0} 50.0)
          s2 (ride s1 {:throttle 0.3 :lean 0.0} 5.0)
          v1 (double (:speed (tel s1)))
          v2 (double (:speed (tel s2)))
          ;; Where the motor's power at this throttle would exactly meet
          ;; the air and the rolling loss, with nothing else lost.
          ideal (loop [v 10.0]
                  (let [push (/ (* 0.3 (double (:peak-power cfg))) v)
                        v' (+ v (* 0.001 (- push (m/resistance cfg v))))]
                    (if (< (abs (- v' v)) 1e-9) v' (recur v'))))]
      (is (< (abs (- v2 v1)) 0.2) (str "still changing: " v1 " then " v2))
      ;; Short of the ideal, and not by much: the rest is lost in the
      ;; tyres, which the contact solver does not model as lossless.
      (is (< (* 0.75 ideal) v2 ideal) (str v2 " against an ideal of " ideal)))))
