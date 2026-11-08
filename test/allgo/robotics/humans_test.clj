(ns allgo.robotics.humans-test
  (:require [allgo.random :as random]
            [allgo.robotics.humans :as humans]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Robots among people, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 26.8.

(def ^:private actions (vec (for [a (range 8)] [(* 0.5 (math/cos (* a (/ math/PI 4)))) (* 0.5 (math/sin (* a (/ math/PI 4))))])))
(def ^:private goals {:window [10.0 5.0] :kitchen [0.0 5.0] :hallway [5.0 10.0]})

(deftest goal-prediction-test
  (let [steps (for [i (range 6)] [[(+ 5.0 (* 0.5 i)) 5.0] [0.5 0.0]])
        beliefs (vec (humans/belief-seq goals steps actions 2.0))]
    (testing "a person heading right from the middle of the room: each step makes the window likelier and the kitchen, behind them, less"
      (is (apply < (map :window beliefs)))
      (is (apply > (map :kitchen beliefs)))
      (is (> (:window (peek beliefs)) 0.99)))
    (testing "the belief predicts the next step: onward, toward the window, most likely"
      (let [p (humans/predict-step goals (peek beliefs) [8.0 5.0] actions 2.0)]
        (is (= [0.5 0.0] (key (apply max-key val p))))
        (is (< (abs (- 1.0 (reduce + (vals p)))) 1e-12))))))

(def ^:private dt 0.1)
(defn- car [[y psi] [u]] (let [u (max -1.0 (min 1.0 u))] [(+ y (* dt 2.0 (math/sin psi))) (+ psi (* dt u))]))
(defn- expert [[y psi]] [(max -1.0 (min 1.0 (- (- y) (* 2.5 psi))))])
(defn- noisy-car [rng] (fn [x u] (let [[y psi] (car x u)] [(+ y (random/gaussian rng 0.0 0.05)) (+ psi (random/gaussian rng 0.0 0.03))])))

(deftest dagger-test
  ;; keeping a car in its lane, buffeted by noise, with two demonstrations
  ;; and a nearest-neighbor learner
  (doseq [seed (range 1 5)]
    (let [rng (random/rng seed)
          draw (fn [] [(* 2.0 (- (* 2 (rng)) 1)) (* 0.4 (- (* 2 (rng)) 1))])
          starts (vec (repeatedly 2 draw))
          tests (vec (repeatedly 20 draw))
          offset (fn [policy] (let [dynamics (noisy-car (random/rng 99))]
                                (/ (reduce + (map (fn [x0] (reduce + (map #(abs (first %)) (humans/rollout dynamics policy x0 80)))) tests))
                                   (* 81.0 (count tests)))))
          train (noisy-car (random/rng (+ 100 seed)))
          cloned (humans/behavioral-cloning train expert humans/nearest-neighbor starts 80)
          dagger (:policy (nth (humans/dagger-seq train expert humans/nearest-neighbor starts 80) 4))]
      (testing "behavioral cloning drifts where it has no examples; four rounds of DAgger stay within 30% of the expert, and do no worse than cloning"
        (is (< (offset dagger) (* 1.3 (offset expert))))
        (is (<= (offset dagger) (offset cloned)))))))
