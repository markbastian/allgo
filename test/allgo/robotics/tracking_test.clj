(ns allgo.robotics.tracking-test
  (:require [allgo.robotics.tracking :as track]
            [clojure.test :refer [deftest is testing]]))

;; Trajectory tracking control, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 26.5.3: a planar
;; two-link arm driven along a retimed path.

(def ^:private arm {:m1 1.0 :m2 1.0 :l1 1.0 :l2 1.0 :g 9.81})

(deftest dynamics-test
  (testing "with no torque the arm swings and its energy stays what it was: the dynamics are consistent"
    (let [s0 [[0.3 0.5] [0.0 0.0]]
          s (nth (iterate #(track/step arm % [0.0 0.0] 0.001) s0) 3000)]
      (is (< (abs (- (apply track/energy arm s0) (apply track/energy arm s))) 1e-6))))
  (testing "inverse dynamics undoes forward dynamics"
    (let [q [0.4 -0.7] qd [0.3 1.1] u [2.0 -1.5]]
      (is (every? #(< (abs %) 1e-9) (map - u (track/inverse-dynamics arm q qd (track/forward-dynamics arm q qd u))))))))

(deftest retime-test
  (testing "long enough to reach the top speed: accelerate, cruise, brake -- D / v + v / a"
    (let [{:keys [duration at]} (track/retime [[0.0 0.0] [3.0 4.0]] 1.0 0.5)]
      (is (= 7.0 duration))
      (is (= [[0.0 0.0] [0.0 0.0]] (subvec (at 0.0) 0 2)))
      (is (every? #(< (abs %) 1e-12) (map - [3.0 4.0] (first (at duration)))))
      (is (every? #(< (abs %) 1e-12) (map - [0.6 0.8] (second (at 3.0)))))))
  (testing "too short: a triangle, 2 sqrt(D / a)"
    (is (= 2.0 (:duration (track/retime [[0.0] [1.0]] 2.0 1.0))))))

(deftest tracking-test
  (let [traj (track/retime [[-0.5 0.8] [0.9 -0.4]] 1.0 2.0)
        ;; a model whose masses are a fifth too light
        wrong (assoc arm :m1 0.8 :m2 0.8)
        n (long (/ (:duration traj) 0.002))
        error (fn [controller] (track/tracking-error (take n (track/track-seq arm controller traj 0.002))))]
    (testing "open loop, the inverse dynamics of the trajectory: off by a tenth of a radian even with the exact model, as errors grow unchecked; by radians with the wrong one"
      (is (< 0.05 (error (track/feedforward arm)) 0.3))
      (is (> (error (track/feedforward wrong)) 2.0)))
    (testing "PD feedback alone, no model, lags and sags under gravity: a third of a radian at these gains"
      (is (> (error (track/pd 100 20)) 0.3)))
    (testing "computed torque, feedforward from the wrong model plus feedback scaled by its inertia: a tenth of PD's error"
      (is (< (error (track/computed-torque wrong 100 20)) 0.05)))
    (testing "and with the exact model, almost none"
      (is (< (error (track/computed-torque arm 100 20)) 1e-3)))))
