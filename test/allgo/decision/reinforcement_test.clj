(ns allgo.decision.reinforcement-test
  (:require [allgo.decision.mdp :as mdp]
            [allgo.decision.reinforcement :as rl]
            [allgo.learning.models :as models]
            [allgo.learning.neural :as nn]
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

;; Reward shaping, section 22.4.4.

(defn- q-iteration
  "Q(s, a) = R(s) + sum over s' of P(s' | s, a) [F(s, s') + gamma max Q(s',
  .)], a terminal state worth its reward: the Q the learners converge to,
  with shaping F."
  [{:keys [states actions transitions reward gamma]} F sweeps]
  (nth (iterate (fn [Q]
                  (into {} (for [s states a (actions s)]
                             [[s a] (+ (reward s)
                                       (reduce + (for [[p s'] (transitions s a)]
                                                   (* p (+ (F s s') (* gamma (if (seq (actions s')) (apply max (map #(Q [s' %]) (actions s'))) (reward s'))))))))])))
                (into {} (for [s states a (actions s)] [[s a] 0.0])))
       sweeps))

(def ^:private sparse
  "A 10 x 10 world, its only reward +1 at the top right."
  (mdp/grid-mdp (vec (for [y (range 9 -1 -1)] (vec (for [x (range 10)] (if (= [x y] [9 9]) 1.0 0.0))))) [[9 9]] 0.95))

(defn- to-goal [[x y]] (+ (- 9 x) (- 9 y)))

(deftest shaping-test
  (testing "shaping by a potential, R + gamma Phi(s') - Phi(s), changes Q by exactly -Phi(s), and so not the optimal policy (Phi of a terminal zero)"
    (let [{:keys [states actions gamma terminals]} w
          phi (fn [[x y]] (+ (* 0.3 x) (* -0.2 y) 0.5))
          shaped (q-iteration w (fn [s s'] (- (* gamma (if (terminals s') 0.0 (phi s'))) (phi s))) 300)
          plain (q-iteration w (constantly 0.0) 300)]
      (doseq [s states a (actions s)]
        (is (< (abs (- (shaped [s a]) (- (plain [s a]) (phi s)))) 1e-9)))
      (doseq [s states :when (seq (actions s))]
        (is (= (apply max-key #(plain [s %]) (actions s)) (apply max-key #(shaped [s %]) (actions s)))))))
  (testing "with a sparse reward, shaping by Phi(s) = gamma^(distance to the goal) cuts the moves epsilon-greedy Q-learning makes in its first 50 trials by more than half"
    (doseq [seed [1 2 3]]
      (let [moves (fn [opts]
                    (let [rng (random/rng seed)]
                      (loop [i 0 [p m] (rl/q-learning-agent sparse (merge {:ne 0 :epsilon 0.2 :rng (random/rng (+ 10 seed))} opts)) total 0]
                        (if (= i 50)
                          total
                          (let [t (vec (rl/trial-seq sparse [p m] [0 0] rng 20000))]
                            (recur (inc i) [p (:memory (peek t))] (+ total (count t))))))))]
        (is (< (* 2 (moves {:potential #(Math/pow 0.95 (to-goal %))})) (moves {})) seed))))
  (testing "a pseudoreward that is not of that form changes what is optimal: a bonus beside the goal, and the best policy stays by it rather than finish"
    (let [bonus (assoc-in sparse [:reward [9 8]] 0.3)
          policy (mdp/best-policy bonus (mdp/value-iteration bonus 1e-6))
          true-value (fn [pi] ((mdp/policy-evaluation sparse pi (zipmap (:states sparse) (repeat 0.0)) 1000) [0 0]))]
      (is (not= [0 1] (policy [9 8])))
      (is (< (true-value policy) (* 0.5 (true-value (mdp/best-policy sparse (mdp/value-iteration sparse 1e-6)))))))))

;; Function approximation, section 22.4.

(defn- world
  "A 10 x 10 world, -0.04 a step and +1 at `goal`, gamma 0.99."
  [goal]
  (mdp/grid-mdp (vec (for [y (range 9 -1 -1)] (vec (for [x (range 10)] (if (= [x y] goal) 1.0 -0.04))))) [goal] 0.99))

(defn- xy
  "The chapter's features, the coordinates from (1, 1)."
  [[x y]]
  [1.0 (inc x) (inc y)])

(defn- rms [f g ss] (Math/sqrt (/ (reduce + (map #(let [d (- (f %) (g %))] (* d d)) ss)) (count ss))))

(defn- run-trials
  "The agent's memory after `n` trials, each from a random state and
  cut off after 500 moves."
  [w agent n rng]
  (let [starts (vec (remove (:terminals w) (:states w)))]
    (loop [i 0 [p m] agent]
      (if (= i n) m (recur (inc i) [p (:memory (last (rl/trial-seq w [p m] (random/pick rng starts) rng 500)))])))))

(deftest widrow-hoff-test
  (testing "the chapter's example: theta = (0.5, 0.2, 0.1) gives U(1,1) = 0.8, and an observed 0.4 lowers each parameter by 0.4 alpha"
    (is (< (abs (- 0.8 (rl/linear-utility [0.5 0.2 0.1] xy [0 0]))) 1e-12))
    (is (every? #(< (abs %) 1e-12) (map - (rl/widrow-hoff [0.5 0.2 0.1] (xy [0 0]) 0.4 0.1) [0.46 0.16 0.06])))))

(deftest function-approximation-test
  (let [examples (fn [w pi n rng]
                   ;; direct utility estimation's training examples: each state of a trial, and its reward-to-go
                   (let [starts (vec (remove (:terminals w) (:states w)))]
                     (vec (mapcat (fn [_]
                                    (let [t (vec (rl/trial-seq w [(fn [m [s _]] [(pi s) m]) nil] (random/pick rng starts) rng))]
                                      (map vector (map :state t) (reverse (reductions (fn [g r] (+ r (* (:gamma w) g))) 0.0 (map :reward (rseq t)))))))
                                  (range n)))))
        fit (fn [ex features]
              (let [theta (models/linear-regression (mapv (comp vec rest features first) ex) (mapv second ex))]
                #(rl/linear-utility theta features %)))]
    (testing "least squares on the trials' examples: with +1 at (10, 10) the true utility is nearly linear in x and y; with +1 at (5, 5), a pyramid, it is not -- until a feature for the distance to the goal is added"
      (doseq [[goal features good?] [[[9 9] xy true]
                                     [[4 4] xy false]
                                     [[4 4] (fn [[x y]] [1.0 (inc x) (inc y) (Math/sqrt (+ (* (- x 4) (- x 4)) (* (- y 4) (- y 4))))]) true]]]
        (let [w (world goal)
              U (mdp/value-iteration w 1e-6)
              err (rms (fit (examples w (mdp/best-policy w U) 300 (random/rng 1)) features) U (remove (:terminals w) (:states w)))]
          (if good? (is (< err 0.05) goal) (is (> err 0.1) goal)))))
    (let [w (world [9 9])
          U (mdp/value-iteration w 1e-6)
          pi (mdp/best-policy w U)
          states (remove (:terminals w) (:states w))]
      (testing "online, by the Widrow-Hoff rule after each trial and by TD after each move, near that fit"
        (is (< (rms #(rl/linear-utility (:theta (run-trials w (rl/approximate-direct-utility-agent w pi xy [0.0 0.0 0.0] 0.0005) 2000 (random/rng 1))) xy %) U states) 0.06))
        (is (< (rms #(rl/linear-utility (:theta (run-trials w (rl/approximate-td-agent w pi xy [0.0 0.0 0.0] 0.0005) 2000 (random/rng 1))) xy %) U states) 0.08)))
      (testing "generalization: five trials from (1, 1) and the linear estimate is good even where the table has nothing"
        (let [table (:U (last (rl/trials w (rl/direct-utility-agent w pi) [0 0] 5 (random/rng 2))))
              theta (:theta (last (rl/trials w (rl/approximate-direct-utility-agent w pi xy [0.0 0.0 0.0] 0.002) [0 0] 5 (random/rng 2))))
              unseen (remove #(contains? table %) states)]
          (is (> (count unseen) 40))
          (is (< (rms #(rl/linear-utility theta xy %) U unseen) 0.1)))))))

(def ^:private headings [[1 0] [0 1] [-1 0] [0 -1]])

(defn- mean-value [w pi]
  (let [U (mdp/policy-evaluation w pi (zipmap (:states w) (repeat 0.0)) 1000)
        ss (remove (:terminals w) (:states w))]
    (/ (reduce + (map U ss)) (count ss))))

(deftest approximate-q-test
  (testing "Q-learning with Q linear in x and y for each action, equation 22.12: with +1 at (10, 10), a greedy policy nearly the best"
    (let [w (world [9 9])
          lin (rl/linear-q (fn [[x y] a] (vec (mapcat #(if (= % a) [1.0 (/ (inc x) 10.0) (/ (inc y) 10.0)] [0.0 0.0 0.0]) headings))))
          m (run-trials w (rl/approximate-q-agent w lin (vec (repeat 12 0.0)) {:alpha 0.05 :rng (random/rng 5)}) 1000 (random/rng 6))]
      (is (> (mean-value w (rl/approximate-greedy-policy w lin (:params m)))
             (- (mean-value w (mdp/best-policy w (mdp/value-iteration w 1e-6))) 0.02)))))
  (testing "deep Q-learning, section 22.4.3, on a 5 x 5 world with +1 in the middle: a network from (x, y), with experience replay, learns what no Q linear in x and y can"
    (let [w (mdp/grid-mdp (vec (for [y (range 4 -1 -1)] (vec (for [x (range 5)] (if (= [x y] [2 2]) 1.0 -0.04))))) [[2 2]] 0.9)
          best (mean-value w (mdp/best-policy w (mdp/value-iteration w 1e-6)))
          centered (fn [[x y]] [(- (/ x 2.0) 1) (- (/ y 2.0) 1)])
          net (rl/network-q centered headings)
          lin (rl/linear-q (fn [s a] (vec (mapcat #(if (= % a) (into [1.0] (centered s)) [0.0 0.0 0.0]) headings))))]
      (doseq [seed [1 2]]
        (let [deep (run-trials w (rl/approximate-q-agent w net [(nn/dense 2 32 :relu (random/rng seed)) (nn/dense 32 4 :identity (random/rng (+ 10 seed)))]
                                                         {:alpha 0.01 :replay 8 :rng (random/rng seed)})
                               300 (random/rng 6))
              linear (run-trials w (rl/approximate-q-agent w lin (vec (repeat 12 0.0)) {:alpha 0.01 :replay 8 :rng (random/rng seed)}) 300 (random/rng 6))]
          (is (> (mean-value w (rl/approximate-greedy-policy w net (:params deep))) (- best 0.02)) seed)
          (is (< (mean-value w (rl/approximate-greedy-policy w lin (:params linear))) (- best 0.3)) seed))))))

;; Safe exploration, section 22.3.2.

(deftest safe-exploration-test
  (let [h1 w
        h2 (assoc-in w [:reward [3 1]] -20.0)
        optimal #(mdp/best-policy % (mdp/value-iteration % 1e-6))
        risky (optimal h1)
        careful (optimal h2)]
    (testing "two hypotheses about the 4 x 3 world, the pit costing 1 or 20: their optimal policies differ by the pit, the careful one going the long way round"
      (is (= [0 1] (risky [2 1])))
      (is (= [-1 0] (careful [2 1]))))
    (testing "the Bayesian policy takes the risk while the costly pit is unlikely, and not once it is likely enough"
      (is (= risky (rl/bayesian-policy [h1 h2] [0.95 0.05] [risky careful])))
      (is (= careful (rl/bayesian-policy [h1 h2] [0.7 0.3] [risky careful]))))
    (testing "the robust policy, the best in the worst case, is careful whatever the probabilities"
      (is (= careful (rl/robust-policy [h1 h2] [risky careful]))))
    (testing "U_h^pi itself: in the 4 x 3 world, the risky policy is worth more by its own model, and much less by the other"
      (is (> (rl/policy-utility h1 risky) (rl/policy-utility h1 careful)))
      (is (< (rl/policy-utility h2 risky) (- (rl/policy-utility h2 careful) 0.5))))))
