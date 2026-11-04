(ns allgo.decision.hierarchical-test
  (:require [allgo.decision.hierarchical :as hrl]
            [allgo.decision.mdp :as mdp]
            [allgo.decision.reinforcement :as rl]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Hierarchical reinforcement learning, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, section
;; 22.4.5.

(defn- learn
  "The agent's memory after `n` trials from `start`, and the moves each
  took."
  [w agent start n rng]
  (loop [i 0 [p m] agent moves []]
    (if (= i n)
      {:memory m :moves moves}
      (let [t (rl/trial-seq w [p m] start rng 20000)]
        (recur (inc i) [p (:memory (last t))] (conj moves (count t)))))))

(deftest flat-program-test
  (testing "HRL includes ordinary RL: with the program `while true do choose(A(s))`, the 4 x 3 world's greedy policy nearly the best"
    (let [w mdp/four-by-three
          optimal (mdp/value-iteration w 0.0001)]
      (doseq [seed [1 2 3]]
        (let [{:keys [Q]} (:memory (learn w (hrl/hrl-agent w (hrl/flat-program w) [:top] {:rng (random/rng seed)}) [0 0] 600 (random/rng (+ 10 seed))))
              pi (into {} (for [s (:states w) :let [as ((:actions w) s)]]
                            [s (when (seq as) (apply max-key #(get Q [[s [:top]] %] 0.0) as))]))]
          (is (> ((mdp/policy-evaluation w pi (zipmap (:states w) (repeat 0.0)) 500) [0 0]) (- (optimal [0 0]) 0.1)) seed))))))

;; Two rooms, the wall at x = 5 with its door at the bottom, +1 at the top
;; right and no other reward; from the top left.

(def ^:private rooms
  (mdp/grid-mdp (vec (for [y (range 6 -1 -1)]
                       (vec (for [x (range 11)] (cond (and (= x 5) (pos? y)) nil (= [x y] [10 6]) 1.0 :else 0.0)))))
                [[10 6]] 0.95))

(def ^:private door [5 0])
(def ^:private goal [10 6])

(defn- distance [[x y] [u v]] (+ (abs (- x u)) (abs (- y v))))

(defn- rooms-program
  "At the top level, choose a target, the door or the goal; NAVIGATE(t)
  chooses moves until it is at t, then returns."
  [{:keys [actions]}]
  (fn [s [routine t a]]
    (case routine
      :top {:choose [door goal] :then (fn [t] [:nav t])}
      :nav (if (= s t) {:then [:top]} {:choose (actions s) :then (fn [a] [:step t a])})
      :step {:act a :then [:nav t]})))

(defn- toward-target
  "Shaping in the joint state space: inside NAVIGATE(t), gamma^(distance
  to t). The target is no part of the physical state."
  [s [routine t]]
  (if (= routine :nav) (Math/pow 0.95 (distance s t)) 0.0))

(deftest rooms-test
  (doseq [seed [1 2 3]]
    (let [flat (learn rooms (hrl/hrl-agent rooms (hrl/flat-program rooms) [:top] {:rng (random/rng seed)}) [0 6] 60 (random/rng (+ 10 seed)))
          hier (learn rooms (hrl/hrl-agent rooms (rooms-program rooms) [:top] {:rng (random/rng seed) :potential toward-target})
                      [0 6] 60 (random/rng (+ 10 seed)))
          Q (:Q (:memory hier))]
      (testing "the partial program, with shaping by the distance to the routine's target, learns in under half the moves flat Q-learning takes"
        (is (< (* 2 (reduce + (:moves hier))) (reduce + (:moves flat))) seed))
      (testing "the same square in two machine states: in the right room, NAVIGATE(door) values moving down or left most, NAVIGATE(goal) up or right"
        (let [best (fn [t] (apply max-key #(get Q [[[7 2] [:nav t]] %] -1.0) ((:actions rooms) [7 2])))]
          (is (#{[-1 0] [0 -1]} (best door)) seed)
          (is (#{[1 0] [0 1]} (best goal)) seed))))))
