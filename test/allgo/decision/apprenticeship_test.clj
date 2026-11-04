(ns allgo.decision.apprenticeship-test
  (:require [allgo.decision.apprenticeship :as app]
            [allgo.decision.mdp :as mdp]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Apprenticeship and inverse reinforcement learning, after Russell and
;; Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
;; section 22.6.

;; A 5 x 5 world: +1 at the top right, -1 at [4 1], both terminal, and a
;; puddle costing 0.5 between them and the start at the bottom left. Its
;; reward is theta . f(s), the features being in the goal, in the pit, in
;; the puddle, and 1.

(def ^:private puddle #{[1 1] [2 1] [1 2] [2 2] [3 2] [2 3]})

(defn- features [s]
  [(if (= s [4 4]) 1.0 0.0) (if (= s [4 1]) 1.0 0.0) (if (puddle s) 1.0 0.0) 1.0])

(def ^:private theta* [1.0 -1.0 -0.5 -0.04])
(def ^:private base (mdp/grid-mdp (vec (repeat 5 (vec (repeat 5 0.0)))) [[4 4] [4 1]] 0.9))
(def ^:private world (app/with-reward base features theta*))
(def ^:private expert (app/optimal-policy world))

(defn- value [pi] (mdp/policy-evaluation world pi (zipmap (:states world) (repeat 0.0)) 1000))

(defn- mean-value [pi]
  (let [U (value pi) ss (remove (:terminals world) (:states world))] (/ (reduce + (map U ss)) (count ss))))

(def ^:private first-heading (fn [_] [1 0]))
(def ^:private always-left (into {} (for [s (:states base)] [s (when-not ((:terminals base) s) [-1 0])])))

(deftest feature-expectations-test
  (testing "U(s) = theta . mu(pi): the utility of a policy is the reward weights times its feature expectations"
    (is (< (abs (- ((value expert) [0 0]) (reduce + (map * theta* (app/feature-expectations base expert features [0 0]))))) 1e-9)))
  (testing "and the expectations estimated from 3000 trajectories come near the exact ones"
    (let [demos (app/demonstrations world expert [0 0] 3000 (random/rng 1))]
      (is (every? #(< (abs %) 0.02) (map - (app/empirical-feature-expectations world demos features)
                                         (app/feature-expectations base expert features [0 0])))))))

(deftest feature-matching-test
  (let [demos (app/demonstrations world expert [0 0] 20 (random/rng 20))
        mu-e (app/empirical-feature-expectations world demos features)
        steps (vec (take 15 (app/feature-matching-seq base features mu-e always-left [0 0])))
        learned (app/closest-policy (peek steps) mu-e)]
    (testing "the margin |theta| by which the expert still beats the policies found never grows"
      (is (apply >= (map :t steps))))
    (testing "from 20 of the expert's trajectories, and never its reward, a policy as good as the expert's by that reward -- from the start, and from every state"
      (is (< (abs (- ((value learned) [0 0]) ((value expert) [0 0]))) 0.01))
      (is (> (mean-value learned) (- (mean-value expert) 0.05))))
    (testing "imitation learning from the same trajectories does worse, and worst off the expert's path: brittle"
      (let [imitated (app/imitation-policy world demos first-heading)]
        (is (< ((value imitated) [0 0]) ((value learned) [0 0])))
        (is (< (mean-value imitated) (- (mean-value learned) 0.1)))))))

(deftest boltzmann-test
  (let [demos (app/demonstrations world expert [0 0] 20 (random/rng 20))
        choices (count (filter second (apply concat demos)))]
    (testing "R = 0 makes every behavior optimal and so explains none: each choice a uniform one of four"
      (is (< (abs (- (app/boltzmann-log-likelihood base (constantly 0.0) demos 10.0) (* choices (Math/log 0.25)))) 1e-9)))
    (testing "the true reward explains the expert's choices far better"
      (is (> (app/boltzmann-log-likelihood base #(reduce + (map * theta* (features %))) demos 10.0)
             (+ 50 (app/boltzmann-log-likelihood base (constantly 0.0) demos 10.0)))))))
