(ns allgo.courses-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]
            [allgo.simulation.courses :as c]
            [allgo.simulation.motorcycle :as m]
            [clojure.test :refer [deftest is testing]]))

(def ^:private cfg (m/with-rider m/defaults))
(def ^:private dt (/ 1.0 240.0))

(defn- ride
  "Ride `course` on its own autopilot for up to `seconds`, stopping at the
  finish or when the rider comes off. Returns the scene and the time.

  The autopilot is asked once a frame and held for the frame's four
  substeps, as the demo asks it -- the pitch control in the air depends
  on how often it is asked, so this is the rate that has to work."
  [id seconds]
  (let [crs (c/course id)]
    (loop [s (c/scene cfg crs {:rider? true}) k 0 controls nil]
      (let [pose (m/bike-pose s)
            t (* k dt)
            controls (if (zero? (mod k 4))
                       ((:autopilot crs) pose (m/telemetry cfg pose) {:cfg cfg :max-lean 0.7})
                       controls)]
        (if (or (c/finished? crs pose)
                (not (get-in s [:rider :attached?]))
                (> t seconds))
          [s t]
          (recur (m/step cfg s controls dt) (inc k) controls))))))

(deftest slope-test
  (testing "a ramp's riding surface runs between the two points it was given"
    (let [r (c/slope [10.0 0.0] [16.0 1.5] 4.0 2.0)
          [len thick _] (:size r)]
      (is (< (v/distance (rigid/local->world r [(* -0.5 len) (* 0.5 thick) 0.0]) [10.0 0.0 0.0]) 1e-9))
      (is (< (v/distance (rigid/local->world r [(* 0.5 len) (* 0.5 thick) 0.0]) [16.0 1.5 0.0]) 1e-9)))))

(deftest flight-test
  (testing "a faster approach needs a longer gap"
    (is (< (:gap (c/flight 6.0 1.2 12.0 0.4)) (:gap (c/flight 6.0 1.2 16.0 0.4)))))
  (testing "the landing slopes a little shallower than the bike is falling as it arrives"
    ;; Fly the arc by hand and compare.
    (let [run 6.0 h 1.5 speed 15.0 touch 0.4
          {:keys [landing]} (c/flight run h speed touch)
          angle (Math/atan2 h run)
          v (Math/sqrt (- (* speed speed) (* 2.0 9.81 h)))
          vx (* v (Math/cos angle)) vy (* v (Math/sin angle))
          t (/ (+ vy (Math/sqrt (+ (* vy vy) (* 2.0 9.81 touch h)))) 9.81)
          falling (/ (- (* 9.81 t) vy) vx)]
      (is (< (abs (- (/ h landing) (* 0.7 falling))) 1e-9)))))

(deftest excitebike-test
  (testing "the Excitebike lane can be ridden flat out, every jump landed, rider aboard"
    (let [[s t] (ride :excitebike 40.0)]
      (is (c/finished? (c/course :excitebike) (m/bike-pose s)) (str "stopped at " t "s"))
      (is (get-in s [:rider :attached?]) "the rider is still on")
      (is (< t 26.0) (str "in " t "s")))))

(deftest trials-test
  (testing "the trials yard can be ridden at a walk: logs, steps and rocks"
    (let [[s t] (ride :trials 40.0)]
      (is (c/finished? (c/course :trials) (m/bike-pose s)) (str "stopped at " t "s"))
      (is (get-in s [:rider :attached?])))))

(deftest loop-test
  (testing "the loop is lapped over its curbs without drifting off it"
    (let [[s _] (ride :loop 12.0)
          [x _ z] (:pos (:base (m/bike-pose s)))]
      (is (get-in s [:rider :attached?]))
      (is (< (abs (- (v/length [x 0.0 (+ (double z) 35.0)]) 35.0)) 4.0)))))

(deftest hard-crash-test
  (testing "leaned hard over at speed, the bike and rider come down and stay in one piece"
    ;; The loop, cranked over further than the autopilot would go. This
    ;; is the ride that took the rider apart while still pinned to an
    ;; upright bike -- a forearm's damper overshooting every step -- and
    ;; whose crash then took the solver apart too. It has to fall over,
    ;; land, and stop.
    (let [crs (c/course :loop)
          free {:throttle 0.0 :brake 0.0 :lean 0.0 :hands? false}
          finite? (fn [s] (every? #(Double/isFinite (double %))
                                  (concat (:qd (m/bike-pose s)) (:qd (m/rider-pose s))
                                          (:vel (:base (m/bike-pose s))) (:vel (:base (m/rider-pose s))))))
          fastest (fn [s] (max (v/length (subvec (vec (:vel (:base (m/bike-pose s)))) 3 6))
                               (v/length (subvec (vec (:vel (:base (m/rider-pose s)))) 3 6))))
          crashed (loop [s (c/scene cfg crs {:rider? true}) k 0 controls nil]
                    (let [pose (m/bike-pose s)
                          controls (if (zero? (mod k 4))
                                     (assoc ((:autopilot crs) pose (m/telemetry cfg pose) {:cfg cfg :max-lean 0.7})
                                            :lean 0.9 :throttle 0.8)
                                     controls)]
                      (if (or (not (get-in s [:rider :attached?])) (> k 960))
                        s
                        (recur (m/step cfg s controls dt) (inc k) controls))))
          after (nth (iterate #(m/step cfg % free dt) crashed) 720)]
      (is (not (get-in crashed [:rider :attached?])) "it did crash")
      (is (finite? after))
      (is (< (fastest after) 1.0) (str "and came to rest, fastest " (fastest after))))))
