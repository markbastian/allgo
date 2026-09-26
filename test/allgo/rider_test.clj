(ns allgo.rider-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.world :as pw]
            [allgo.simulation.motorcycle :as m]
            [allgo.simulation.rider :as r]
            [clojure.test :refer [deftest is testing]]))

(def ^:private cfg (m/with-rider m/defaults))
(def ^:private dt (/ 1.0 240.0))
(def ^:private floor (m/ground 1000.0))

(defn- ride [scene controls seconds]
  (m/run cfg scene controls dt (long (/ (double seconds) dt))))

(defn- on-bike [speed]
  (m/scene cfg [floor] (m/start-pose cfg speed) {:rider? true}))

(defn- pin-gaps [scene]
  (let [configs (mapv (fn [{:keys [model pose]}] (ab/configuration model (:q pose) (:base pose)))
                      (:models scene))]
    (mapv #(apply v/distance (pw/pin-points scene configs %)) (:pins scene))))

(deftest seated-test
  (testing "the figure is built already in the saddle"
    ;; Its rest pose is the riding pose, so at the start every pin is
    ;; closed and every muscle is at rest.
    (let [s (on-bike 8.0)
          rider (second (:models s))]
      (is (= 5 (count (:pins s))))
      (is (every? #(< (double %) 1e-9) (pin-gaps s)) (str (pin-gaps s)))
      (is (every? #(< (abs (double %)) 1e-12) (r/muscles (:model rider) (:pose rider) 1.0))))))

(deftest reaches-test
  (testing "the arms and legs are the lengths they are said to be, and reach"
    (let [{:keys [shoulder elbow hand hip knee foot]} (r/layout cfg)]
      (doseq [s [:left :right]]
        (is (< (abs (- 0.30 (v/distance (shoulder s) (elbow s)))) 1e-9))
        (is (< (abs (- 0.30 (v/distance (elbow s) (hand s)))) 1e-9))
        (is (< (abs (- 0.42 (v/distance (hip s) (knee s)))) 1e-9))
        (is (< (abs (- 0.42 (v/distance (knee s) (foot s)))) 1e-9))))))

(deftest stays-on-test
  (testing "riding and turning, the rider stays on and holds a posture"
    (let [s (ride (ride (on-bike 8.0) {:throttle 0.08 :lean 0.0} 2.0)
                  {:throttle 0.08 :lean 0.3} 3.0)
          rider (:pose (second (:models s)))
          spine (nth (:q rider) r/torso)
          ;; How far the spine has turned from the riding pose.
          bent (* 2.0 (Math/acos (min 1.0 (abs (double (nth spine 3))))))]
      (is (get-in s [:rider :attached?]))
      (is (= 5 (count (:pins s))))
      (is (every? #(< (double %) 0.01) (pin-gaps s)) (str "pins apart by " (pin-gaps s)))
      (is (< bent 0.3) (str "spine bent " bent " from the riding pose"))
      (is (:upright? (m/telemetry cfg (m/bike-pose s)))))))

(deftest turns-carrying-a-rider-test
  (testing "with a rider aboard the bike still holds the lean it is asked for"
    ;; A rider is most of the weight and sits high, which is exactly what
    ;; a steer worked out for a rigid bike gets wrong; the rider's feel for
    ;; the error -- the integral -- takes it out.
    (let [s (ride (on-bike 8.0) {:throttle 0.08 :lean 0.3} 7.0)
          lean (double (:lean (m/telemetry cfg (m/bike-pose s))))]
      (is (< (abs (- lean 0.3)) 0.03) (str "leaning " lean)))))

(deftest thrown-test
  (testing "when the bike goes over the rider comes off, goes slack and lands"
    (let [tilted (assoc-in (m/start-pose cfg 10.0) [:base :rot]
                           (q/from-axis-angle [1.0 0.0 0.0] -0.15))
          s (ride (m/scene cfg [floor] tilted {:rider? true}) {:throttle 0.05 :hands? false} 5.0)
          {:keys [model pose]} (second (:models s))
          parts (ab/collision-bodies model (:q pose) (:base pose))
          lowest (apply min (map #(double (second (:pos (:body %)))) parts))]
      (is (not (get-in s [:rider :attached?])))
      (is (empty? (:pins s)))
      (is (< (double (get-in s [:rider :tone])) 0.5))
      ;; Centres of the parts: nothing has gone through the floor.
      (is (pos? lowest) (str "lowest part at " lowest))
      (is (< (double (second (:pos (:base pose)))) 0.6) "and is lying on it"))))
