(ns allgo.decision.reinforcement-test
  (:require [allgo.decision.mdp :as mdp]
            [allgo.decision.reinforcement :as rl]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The 4 x 3 world of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 22, and the fixed policy of
;; aima-python's (MIT) tests.

(def ^:private w mdp/four-by-three)

(def ^:private pi
  {[0 2] [1 0] [1 2] [1 0] [2 2] [1 0] [3 2] nil [0 1] [0 1] [2 1] [0 1] [3 1] nil
   [0 0] [0 1] [1 0] [-1 0] [2 0] [-1 0] [3 0] [-1 0]})

(def ^:private truth (mdp/policy-evaluation w pi (zipmap (:states w) (repeat 0.0)) 500))

(defn- learned [agent n seed] (last (rl/trials w agent [0 0] n (random/rng seed))))

(deftest passive-test
  (testing "each passive learner comes near the policy's true utilities; ADP, which learns a model, nearest"
    (doseq [[agent n tolerance] [[(rl/passive-adp-agent w pi) 100 0.1]
                                 [(rl/passive-td-agent w pi) 300 0.25]
                                 [(rl/direct-utility-agent w pi) 300 0.25]]
            :let [U (:U (learned agent n 3))]]
      (is (= 9 (count U)) "every state the policy reaches")
      (doseq [[s u] U] (is (< (abs (- u (truth s))) tolerance) s)))))

(deftest active-test
  (let [optimal (mdp/value-iteration w 0.0001)]
    ;; Q-learning with more tries before it trusts its Q values than the
    ;; default 5, which on some seeds settles too soon; SARSA with
    ;; exploration that never stops, which it needs.
    (doseq [[make opts] [[rl/q-learning-agent {:ne 30}] [rl/sarsa-agent {:epsilon 0.1}]] seed [1 2 3]]
      (let [rng (random/rng seed)
            agent (make w (assoc opts :rng (random/rng (+ 10 seed))))
            lens (loop [i 0 [p m] agent lens []]
                   (if (= i 600)
                     {:lens lens :memory m}
                     (let [t (vec (rl/trial-seq w [p m] [0 0] rng 5000))]
                       (recur (inc i) [p (:memory (peek t))] (conj lens (count t))))))
            policy (rl/greedy-policy w (:Q (:memory lens)))
            value (mdp/policy-evaluation w policy (zipmap (:states w) (repeat 0.0)) 500)]
        (testing "every trial ends: ties broken at random, no agent circles for ever"
          (is (every? #(< % 5000) (:lens lens))))
        (testing "the greedy policy it learns is worth nearly the best from the start"
          (is (> (value [0 0]) (- (optimal [0 0]) 0.1))))))))

(deftest policy-search-test
  (testing "REINFORCE's softmax policy, made greedy, is worth nearly the best from the start"
    (let [optimal (mdp/value-iteration w 0.0001)]
      (doseq [seed [1 2]]
        (let [theta (nth (rl/reinforce-seq w [0 0] 0.1 (random/rng seed)) 2000)
              greedy (into {} (for [s (:states w)]
                                [s (let [d (rl/softmax-policy w theta s)] (when (seq d) (key (apply max-key val d))))]))
              value (mdp/policy-evaluation w greedy (zipmap (:states w) (repeat 0.0)) 500)]
          (is (> (value [0 0]) (- (optimal [0 0]) 0.02))))))))
