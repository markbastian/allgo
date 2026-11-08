(ns allgo.robotics.control-test
  (:require [allgo.robotics.control :as control]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Control, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, chapter 26: what each PID term does to a mass
;; on a line, and the LQR gain against cases worked out by hand.

(defn- positions [gains disturbance]
  (mapv (comp first :x)
        (take 3000 (control/simulate {:plant (control/mass disturbance) :measure first
                                      :target (constantly 1.0) :gains gains :dt 0.01}
                                     [0.0 0.0]))))

(deftest pid-test
  (testing "proportional alone oscillates without end, from 0 to twice the target"
    (let [xs (drop 2500 (positions {:kp 4.0} 0.0))]
      (is (< (apply min xs) 0.01))
      (is (> (apply max xs) 1.99))))
  (testing "adding the derivative damps it onto the target"
    (is (< (abs (- 1.0 (peek (positions {:kp 4.0 :kd 3.0} 0.0)))) 1e-9)))
  (testing "a constant push leaves PD off by push / kp"
    (is (< (abs (- 1.5 (peek (positions {:kp 4.0 :kd 3.0} 2.0)))) 1e-9)))
  (testing "and the integral term takes that error away"
    (is (< (abs (- 1.0 (peek (positions {:kp 4.0 :kd 3.0 :ki 2.0} 2.0)))) 1e-9))))

(deftest lqr-test
  (testing "x' = x + u with unit costs: P solves P^2 = P + 1, the golden ratio, and K = 1/P"
    (let [{:keys [P K]} (control/lqr [[1.0]] [[1.0]] [[1.0]] [[1.0]])
          phi (/ (+ 1.0 (math/sqrt 5.0)) 2.0)]
      (is (< (abs (- phi (get-in P [0 0]))) 1e-10))
      (is (< (abs (- (/ 1.0 phi) (get-in K [0 0]))) 1e-10))))
  (testing "a double integrator brought to rest, and the gain a fixed point of the recursion"
    (let [dt 0.1
          A [[1.0 dt] [0.0 1.0]] B [[(* 0.5 dt dt)] [dt]]
          Q [[1.0 0.0] [0.0 0.0]] R [[0.1]]
          {:keys [P K steps]} (control/lqr A B Q R)
          {P' :P} (nth (control/riccati-seq A B Q R) (inc steps))
          x (nth (control/closed-loop A B K [1.0 0.0]) 200)]
      (is (every? #(< (abs %) 1e-6) x))
      (is (every? pos? (first K)) "pushing back against displacement and velocity")
      (is (every? #(< (abs %) 1e-9) (map - (flatten P) (flatten P'))) "one more step leaves P where it is"))))

(deftest ilqr-test
  (testing "linear dynamics and a quadratic cost: the expansion is exact, so one iteration is the optimum"
    (let [dt 0.1
          problem {:dynamics (fn [[p v] [u]] [(+ p (* dt v)) (+ v (* dt u))]) :x0 [1.0 0.0] :goal [0.0 0.0]
                   :Q [[1.0 0.0] [0.0 0.1]] :R [[0.01]] :Qf [[10.0 0.0] [0.0 1.0]]}
          [c0 c1 c2] (map :cost (control/ilqr-seq problem (repeat 40 [0.0])))]
      (is (< c1 c0))
      (is (< (abs (- c1 c2)) (* 1e-9 c1)))))
  (let [pendulum (fn [m] (fn [[th w] [u]] (let [a (+ (- (* 9.81 (math/sin th))) (* -0.1 w) (/ u m))]
                                            [(+ th (* 0.05 w)) (+ w (* 0.05 a))])))
        swing {:dynamics (pendulum 1.0) :x0 [0.0 0.0] :goal [math/PI 0.0]
               :Q [[0.1 0.0] [0.0 0.01]] :R [[0.001]] :Qf [[100.0 0.0] [0.0 10.0]]}
        states (vec (take 30 (control/ilqr-seq swing (repeat 60 [0.0]))))
        plan (:us (peek states))]
    (testing "a pendulum swung up from hanging: the cost never rises, and the plan ends upright and still"
      (is (apply >= (map :cost states)))
      (let [[th w] (peek (control/rollout (:dynamics swing) [0.0 0.0] plan))]
        (is (< (abs (- th math/PI)) 1e-3))
        (is (< (abs w) 1e-3))))
    (testing "the plan run open loop on a pendulum 30% heavier falls short"
      (let [[th] (peek (control/rollout (pendulum 1.3) [0.0 0.0] plan))]
        (is (> (abs (- th math/PI)) 1.0))))
    (testing "model predictive control with the same wrong model, replanning every step, gets it up and holds it there"
      (let [states (vec (take 101 (control/mpc-seq (pendulum 1.3) swing 40 3)))]
        (is (every? #(< (abs (- (first (:x %)) math/PI)) 1e-3) (take-last 20 states)))))))
