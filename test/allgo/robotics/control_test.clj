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
