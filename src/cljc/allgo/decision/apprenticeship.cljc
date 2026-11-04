(ns allgo.decision.apprenticeship
  "Apprenticeship learning: Russell and Norvig, *Artificial Intelligence:
  A Modern Approach*, 4th edition, section 22.6 -- learning to behave
  well from an expert's behavior.

    imitation learning   supervised learning of pi(s) from the expert's
                         state-action pairs: here the action the expert
                         took most often in each state, and in a state it
                         never visited, a default -- brittle, since one
                         slip into such a state can end in failure
    feature matching     inverse reinforcement learning: the reward taken
                         to be theta . f(s), and a policy sought whose
                         feature expectations mu(pi) = E[sum over t of
                         gamma^t f(S_t)] match the expert's, estimated
                         from its trajectories -- for then it is as good
                         as the expert by the expert's own reward. Each
                         step finds theta from the policies so far, then
                         the optimal policy for theta . f. The chapter
                         finds theta as the one by which the expert most
                         outperforms them, a quadratic program; this is
                         the projection method from the same paper
                         (Abbeel and Ng, \"Apprenticeship Learning via
                         Inverse Reinforcement Learning\", ICML 2004),
                         which needs none and has the same guarantee
    Boltzmann            the likelihood of the expert's choices under a
    rationality          reward R, each made with probability softmax(beta
                         Q_R(s, .)): why R = 0, which makes every
                         behavior optimal, explains no behavior well

  An MDP is `allgo.decision.mdp`'s, its rewards a function of the state."
  (:require [allgo.decision.mdp :as mdp]
            [allgo.decision.reinforcement :as rl]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn demonstrations
  "`n` trajectories of the policy `pi` -- the expert -- from `start`,
  each `[[s a] ...]`, the action nil in the terminal state."
  [mdp pi start n rng]
  (vec (repeatedly n (fn [] (mapv (juxt :state :action) (rl/trial-seq mdp [(fn [m [s _]] [(pi s) m]) nil] start rng 1000))))))

(defn imitation-policy
  "The action the demonstrations took most often in each state they
  visited, and in every other state of `mdp` with actions, `(default
  s)`."
  [{:keys [states actions]} demos default]
  (let [seen (reduce (fn [acc [s a]] (if a (update-in acc [s a] (fnil inc 0)) acc)) {} (apply concat demos))]
    (into {} (for [s states]
               [s (cond (empty? (actions s)) nil
                        (seen s) (key (apply max-key val (seen s)))
                        :else (default s))]))))

(defn with-reward
  "`mdp` with the reward theta . f(s)."
  [mdp features theta]
  (assoc mdp :reward (fn [s] (lin/dot theta (features s)))))

(defn optimal-policy [mdp] (mdp/best-policy mdp (mdp/value-iteration mdp 1e-6)))

(defn feature-expectations
  "mu(pi) from `start`, exactly: for each feature, the policy's utility
  with that feature as the reward."
  [{:keys [states] :as mdp} pi features start]
  (mapv (fn [i] ((mdp/policy-evaluation (assoc mdp :reward #(nth (features %) i)) pi (zipmap states (repeat 0.0)) 1000) start))
        (range (count (features start)))))

(defn empirical-feature-expectations
  "mu estimated from trajectories: the mean over them of sum over t of
  gamma^t f(s_t)."
  [{:keys [gamma]} demos features]
  (let [per (map (fn [traj] (reduce lin/add (map-indexed (fn [t [s]] (lin/scale (features s) (math/pow gamma t))) traj))) demos)]
    (lin/scale (reduce lin/add per) (/ 1.0 (count demos)))))

(defn feature-matching-seq
  "Inverse reinforcement learning by feature matching, a step at a time,
  toward the expert's feature expectations `mu-expert`, from the policy
  `pi0`: each step

    mu-bar  moved toward mu-expert along the line to the newest policy's
            mu -- the nearest point to mu-expert among the policies'
            mixtures that the method tracks
    theta   mu-expert - mu-bar, the reward weights
    pi      the optimal policy for theta . f, and its mu

  `{:theta :t :policy :mu :policies :mus}`, `t` = |theta|, the margin by
  which the expert's reward still beats the policies found: when it is
  small, one of them is about as good as the expert."
  [mdp features mu-expert pi0 start]
  (let [mu-of #(feature-expectations mdp % features start)
        step (fn [mu-bar policies mus]
               (let [theta (lin/sub mu-expert mu-bar)
                     pi (optimal-policy (with-reward mdp features theta))
                     mu (mu-of pi)]
                 {:mu-bar mu-bar :theta theta :t (lin/length theta) :policy pi :mu mu
                  :policies (conj policies pi) :mus (conj mus mu)}))
        mu0 (mu-of pi0)]
    (iterate (fn [{:keys [mu-bar mu policies mus]}]
               (let [d (lin/sub mu mu-bar)
                     dd (lin/dot d d)
                     mu-bar (if (pos? dd) (lin/add mu-bar (lin/scale d (/ (lin/dot d (lin/sub mu-expert mu-bar)) dd))) mu-bar)]
                 (step mu-bar policies mus)))
             (step mu0 [pi0] [mu0]))))

(defn closest-policy
  "Of the policies feature matching has found, the one whose feature
  expectations are nearest the expert's."
  [{:keys [policies mus]} mu-expert]
  (policies (apply min-key #(lin/distance (mus %) mu-expert) (range (count mus)))))

(defn boltzmann-log-likelihood
  "log P(d | h_R): the log probability of the demonstrations' choices if
  the expert's Q-function is R's and it chooses by softmax(beta Q(s, .))."
  [{:keys [actions gamma] :as mdp} reward demos beta]
  (let [m (assoc mdp :reward reward)
        U (mdp/value-iteration m 1e-6)
        q (fn [s a] (+ (reward s) (* gamma (mdp/q-value m s a U))))]
    (reduce + (for [traj demos [s a] traj :when a
                    :let [qs (map #(* beta (q s %)) (actions s))
                          top (apply max qs)
                          z (reduce + (map #(math/exp (- % top)) qs))]]
                (- (* beta (q s a)) top (math/log z))))))
