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
      ;; Centers of the parts: nothing has gone through the floor.
      (is (pos? lowest) (str "lowest part at " lowest))
      (is (< (double (second (:pos (:base pose)))) 0.6) "and is lying on it"))))

(deftest jolt-test
  (testing "a jolt through the bike does not throw a rider off it"
    ;; A wheel striking something changes the bike's velocity in a step,
    ;; and a rigid pin changes the rider's in the same step: a meter and a
    ;; half a second, upward, is thirty kilonewtons through the seat for
    ;; 1/240 s. The seat used to let go on that and throw the rider off a
    ;; bike that was still upright. Held over a twentieth of a second, as
    ;; a person's body takes it, it is nothing like enough.
    (let [s (ride (on-bike 8.0) {:throttle 0.05 :lean 0.0} 1.0)
          {:keys [model pose]} (first (:models s))
          jolted (assoc-in s [:models 0 :pose]
                           (ab/apply-impulse model pose -1 (:pos (:base pose)) [0.0 1.0 0.0]
                                             (* 1.5 (+ (double (#'m/total-mass cfg)) (r/mass (:model (second (:models s))))))))
          after (ride jolted {:throttle 0.05 :lean 0.0} 1.0)]
      (is (:upright? (m/telemetry cfg (m/bike-pose after))) "the bike rides it out")
      (is (get-in after [:rider :attached?]) "and so does the rider")))
  (testing "but a pin held past its limit for long enough still lets go"
    (let [s (on-bike 8.0)
          weak (update s :pins (fn [ps] (mapv #(if (= :seat (:part %)) (assoc % :break-force 100.0) %) ps)))
          after (ride weak {:throttle 0.05 :lean 0.0} 0.5)]
      (is (not (get-in after [:rider :attached?]))))))

(deftest pins-that-throw-test
  (testing "a foot knocked off its peg leaves the rider on the bike"
    (let [s (update (on-bike 8.0) :pins
                    (fn [ps] (mapv #(if (= [:foot :left] (:part %)) (assoc % :break-force 1.0) %) ps)))
          after (ride s {:throttle 0.05 :lean 0.0} 0.5)]
      (is (not-any? #(= [:foot :left] (:part %)) (:pins after)) "the foot is off")
      (is (get-in after [:rider :attached?]) "and the rider rides on")))
  (testing "but losing a grip throws them"
    (let [s (update (on-bike 8.0) :pins
                    (fn [ps] (mapv #(if (= [:hand :right] (:part %)) (assoc % :break-force 1.0) %) ps)))
          after (ride s {:throttle 0.05 :lean 0.0} 0.5)]
      (is (not (get-in after [:rider :attached?]))))))

(deftest cannonball-test
  (testing "a cannonball at the rider takes them off, and nothing comes apart"
    (let [s (m/fire-at cfg (ride (on-bike 10.0) {:throttle 0.05 :lean 0.0} 1.0))
          ball (fn [s] (some #(when (= :cannonball (:kind %)) %) (:bodies s)))
          states (vec (take 720 (iterate #(m/step cfg % {:throttle 0.05 :lean 0.0} dt) s)))
          closest (apply min (map #(v/distance (:pos (ball %)) (:pos (:base (m/rider-pose %)))) (take 120 states)))
          end (peek states)]
      (is (< closest 0.6) (str "it reached the rider: " closest " m"))
      (is (not (get-in end [:rider :attached?])) "and knocked them off")
      (is (every? #(Double/isFinite (double %))
                  (concat (:qd (m/rider-pose end)) (:qd (m/bike-pose end)) (:pos (ball end))))))))

(deftest leans-into-it-test
  (testing "asked for a turn, the rider tips their upper body into it"
    ;; A lean toward -z, the positive way, is a turn of the spine about
    ;; the figure's forward axis the negative way.
    (let [s (ride (on-bike 8.0) {:throttle 0.08 :lean 0.3} 3.0)
          [x _ _ w] (nth (:q (m/rider-pose s)) r/torso)
          about-forward (* 2.0 (Math/atan2 (double x) (double w)))]
      (is (< about-forward -0.1) (str "spine turned " about-forward " about forward")))))

(defn- deepest-into-bike
  "How far the rider's parts have sunk into the bike's, at worst."
  [scene]
  (let [configs (mapv (fn [{:keys [model pose]}] (ab/configuration model (:q pose) (:base pose)))
                      (:models scene))
        between (filter #(and (= :link (first (:a %))) (= :link (first (:b %)))
                              (not= (second (:a %)) (second (:b %))))
                        (pw/contacts (assoc scene :collide-models? true) configs))]
    (reduce max 0.0 (map #(double (:depth %)) between))))

(deftest thrown-off-test
  (testing "released on a moving bike, the rider is a ragdoll at once"
    (let [s0 (m/throw-rider (ride (on-bike 6.0) {:throttle 0.05 :lean 0.0} 1.0))
          rider (second (:models s0))]
      (is (empty? (:pins s0)) "every pin let go at once")
      (is (zero? (double (get-in s0 [:rider :tone]))) "and every muscle")
      (is (:collide-models? s0) "and the rider now hits the bike")
      (is (every? zero? (r/muscles (:model rider) (:pose rider) 0.0)))))
  (testing "bailing out, the rider goes over the side and ends up on the road, clear of the bike"
    (let [s (ride (m/throw-rider (ride (on-bike 8.0) {:throttle 0.05 :lean 0.0} 1.0)
                                 m/bail-impulse)
                  {:throttle 1.0 :lean 0.0} 4.0)
          rp (:base (m/rider-pose s))
          bp (:base (m/bike-pose s))]
      ;; Whatever is asked of it now, nobody is riding: the bike, left to
      ;; itself, has gone over.
      (is (not (:upright? (m/telemetry cfg (m/bike-pose s)))))
      (is (> (v/distance (:pos rp) (:pos bp)) 1.0) "the rider is clear of the bike")
      (is (< (double (second (:pos rp))) 0.4) "and lying on the road")
      ;; Nothing passing through the bike on the way down.
      (is (< (deepest-into-bike s) 0.05) (str "sunk " (deepest-into-bike s) " into the bike")))))
