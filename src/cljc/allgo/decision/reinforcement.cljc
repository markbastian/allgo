(ns allgo.decision.reinforcement
  "Reinforcement learning: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 22, as steps.

    direct utility    the chapter's first passive learner, which the book
    estimation        gives no figure: each state's utility the average of
                      the rewards-to-go observed from it
    passive ADP       PASSIVE-ADP-AGENT (figure 22.2): learns the
                      transition model by counting, and evaluates the
                      policy on it
    passive TD        PASSIVE-TD-AGENT (figure 22.4): nudges U(s) toward
                      r + gamma U(s') after each move
    Q-learning        Q-LEARNING-AGENT (figure 22.8): action values
                      learned without a model, actions chosen through an
                      exploration function
    SARSA             the chapter's on-policy variant, without a figure
    reward shaping    section 22.4.4: Q-learning on R + gamma Phi(s') -
                      Phi(s) for a potential Phi, which leaves the optimal
                      policy as it was
    function          section 22.4: U or Q a function of parameters theta
    approximation     rather than a table -- linear, theta . f, or a
                      network -- learned by direct utility estimation and
                      the Widrow-Hoff rule (equation 22.10), by TD
                      (equation 22.11), or by Q-learning (equation 22.12),
                      with experience replay if wanted
    policy search     REINFORCE, a softmax policy's preferences moved up
                      the gradient of the return

  An agent is `(fn [memory percept]) -> [action memory']`, a percept
  `[state reward]`, as for every agent in this library; `trial-seq` runs
  one through an MDP of `allgo.decision.mdp` until a terminal state, a
  step a move, and `trials` runs it again and again, its memory carried
  over. The learning rate and the exploration function are aima-python's
  (MIT): alpha(n) = 60 / (59 + n), and f(u, n) = R+ until a state-action
  pair has been tried Ne times."
  (:require [allgo.decision.mdp :as mdp]
            [allgo.learning.neural :as nn]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn alpha [n] (/ 60.0 (+ 59.0 n)))

;; ---------------------------------------------------------------------------
;; Running trials

(defn- sample
  "An outcome of `outcomes`, `[p s]` pairs as an MDP's transitions are."
  [rng outcomes]
  (random/weighted rng (map (fn [[p s]] [s p]) outcomes)))

(defn- trial-step
  "A move: the action's outcome sampled and perceived; nil at a
  terminal state or after `max-moves`."
  [{:keys [actions transitions reward]} program rng max-moves {:keys [state action] :as st}]
  (when (and (seq (actions state)) action (< (long (:moves st)) (long max-moves)))
    (let [s' (sample rng (transitions state action))
          [a m] (program (:memory st) [s' (reward s')])]
      {:state s' :reward (reward s') :action a :memory m :moves (inc (long (:moves st)))})))

(defn trial-seq
  "One trial of `agent`, a `[program memory]` pair, in `mdp` from `start`:
  each state `{:state :reward :action :memory}`, ending after the agent
  perceives a terminal state, or after `max-moves` (default 10000) -- an
  agent can learn to like a wall."
  ([mdp agent start rng] (trial-seq mdp agent start rng 10000))
  ([{:keys [reward] :as mdp} [program memory] start rng max-moves]
   (let [[a m] (program memory [start (reward start)])]
     (search/step-seq (partial trial-step mdp program rng max-moves)
                      {:state start :reward (reward start) :action a :memory m :moves 0}))))

(defn- trial-memory
  "The agent's memory after a trial that starts with memory `m`."
  [mdp program start rng m]
  (:memory (last (trial-seq mdp [program m] start rng))))

(defn trials
  "The agent's memory after each of `n` trials from `start`, memory carried
  from one to the next."
  [mdp [program memory] start n rng]
  (take (inc n)
        (iterate (partial trial-memory mdp program start rng) memory)))

;; ---------------------------------------------------------------------------
;; Passive learners

(defn direct-utility-agent
  "Direct utility estimation for the fixed policy `pi`: the rewards of a
  trial are remembered, and when it ends each state's estimate becomes
  the average of every reward-to-go seen from it."
  [{:keys [actions gamma]} pi]
  [(fn [{:keys [trial totals counts] :as m} [s r]]
     (let [trial (conj trial [s r])]
       (if (seq (actions s))
         [(pi s) (assoc m :trial trial)]
         ;; The trial is over: the reward-to-go from each step, backward.
         (let [[totals counts] (loop [[[s r] & more] (rseq trial) to-go 0.0 totals totals counts counts]
                                 (if s
                                   (let [g (+ r (* gamma to-go))]
                                     (recur more g (update totals s (fnil + 0.0) g) (update counts s (fnil inc 0))))
                                   [totals counts]))]
           [nil {:trial [] :totals totals :counts counts
                 :U (into {} (map (fn [[s t]] [s (/ t (counts s))]) totals))}]))))
   {:trial [] :totals {} :counts {} :U {}}])

(defn passive-adp-agent
  "PASSIVE-ADP-AGENT for the fixed policy `pi` in an MDP whose transitions
  it does not know: it counts the outcomes of each state-action pair to
  estimate P(s' | s, a), and after each move evaluates the policy on its
  model (`k` sweeps, default 20, as aima-python does)."
  ([world pi] (passive-adp-agent world pi 20))
  ([{:keys [states actions gamma]} pi k]
   [(fn [{:keys [U R n-sa n-s'sa s a]} [s' r']]
      (let [U (if (contains? U s') U (assoc U s' r'))
            R (assoc R s' r')
            [n-sa n-s'sa] (if s
                            [(update n-sa [s a] (fnil inc 0)) (update n-s'sa [s a s'] (fnil inc 0))]
                            [n-sa n-s'sa])
            model {:states (filter #(contains? R %) states)
                   :actions actions
                   :gamma gamma
                   :reward R
                   :transitions (fn [s a]
                                  (let [n (n-sa [s a])]
                                    (if n
                                      (for [[[s0 a0 t] c] n-s'sa :when (and (= s0 s) (= a0 a))] [(/ (double c) n) t])
                                      [])))}
            U (mdp/policy-evaluation model #(when (seq (actions %)) (pi %)) (merge (zipmap (keys R) (repeat 0.0)) U) k)
            terminal? (empty? (actions s'))]
        [(when-not terminal? (pi s'))
         {:U U :R R :n-sa n-sa :n-s'sa n-s'sa
          :s (when-not terminal? s') :a (when-not terminal? (pi s'))}]))
    {:U {} :R {} :n-sa {} :n-s'sa {} :s nil :a nil}]))

(defn passive-td-agent
  "PASSIVE-TD-AGENT for the fixed policy `pi`: after each move,
  U(s) += alpha(N(s)) (r + gamma U(s') - U(s))."
  [{:keys [actions gamma]} pi]
  [(fn [{:keys [U n-s s r]} [s' r']]
     (let [U (if (contains? U s') U (assoc U s' r'))
           [U n-s] (if s
                     (let [n-s (update n-s s (fnil inc 0))]
                       [(update U s + (* (alpha (n-s s)) (+ r (* gamma (U s')) (- (U s))))) n-s])
                     [U n-s])
           terminal? (empty? (actions s'))]
       [(when-not terminal? (pi s'))
        {:U U :n-s n-s :s (when-not terminal? s') :r (when-not terminal? r')}]))
   {:U {} :n-s {} :s nil :r nil}])

;; ---------------------------------------------------------------------------
;; Active learners

(defn- explore [ne r+ u n] (if (< (long n) (long ne)) r+ u))

(defn- q-agent
  [{:keys [actions gamma]} {:keys [ne r+ rng epsilon potential] :or {ne 5 r+ 2.0 epsilon 0.0}} on-policy?]
  [(fn [{:keys [Q n-sa s a r]} [s' r']]
     (let [terminal? (empty? (actions s'))
           Q (if terminal? (assoc Q [s' nil] r') Q)
           q (fn [s a] (get Q [s a] 0.0))
           as' (if terminal? [nil] (actions s'))
           ;; The next action: the exploration function's favorite, ties
           ;; broken at random. Broken the same way every time, an agent
           ;; whose actions all sit at the value of wandering for ever --
           ;; -0.04 / (1 - 0.9) = -0.4 in the 4 x 3 world -- circles there
           ;; and never finds its way out: SARSA did, every trial after its
           ;; seventeenth.
           a' (when-not terminal?
                (let [scored (map (fn [a] [a (explore ne r+ (q s' a) (get n-sa [s' a] 0))]) as')
                      top (apply max (map second scored))
                      best (vec (keep (fn [[a v]] (when (>= v (- top 1e-12)) a)) scored))]
                  ;; And with probability epsilon, any action at all.
                  (if (< (double (rng)) epsilon) (random/pick rng (vec as')) (random/pick rng best))))
           [Q n-sa] (if s
                      (let [n-sa (update n-sa [s a] (fnil inc 0))
                            ;; potential-based shaping, the potential of a
                            ;; terminal state taken as zero
                            shaping (if potential (- (* gamma (if terminal? 0.0 (potential s'))) (potential s)) 0.0)
                            target (if on-policy?
                                     (q s' a')
                                     (apply max (map #(q s' %) as')))]
                        [(assoc Q [s a] (+ (q s a) (* (alpha (n-sa [s a])) (+ r shaping (* gamma target) (- (q s a))))))
                         n-sa])
                      [Q n-sa])]
       [a' {:Q Q :n-sa n-sa :s (when-not terminal? s') :a a' :r (when-not terminal? r')}]))
   {:Q {} :n-sa {} :s nil :a nil :r nil}])

(defn q-learning-agent
  "Q-LEARNING-AGENT: Q(s, a) += alpha(N(s, a)) (r + gamma max over a' of
  Q(s', a') - Q(s, a)), and the next action the one the exploration
  function f(Q, N) favors -- `:ne` (default 5) tries of each pair before
  trusting its Q, valued `:r+` (default 2) until then -- or, with
  probability `:epsilon` (default 0), any action. With `:potential`, a
  function Phi of states, its rewards are shaped (section 22.4.4): R'(s,
  a, s') = R + gamma Phi(s') - Phi(s), which changes Q to Q - Phi(s) but
  not the optimal policy; Phi of a terminal state is taken as zero, as
  that needs."
  ([mdp] (q-learning-agent mdp {}))
  ([mdp opts] (q-agent mdp (update opts :rng #(or % (random/rng 1))) false)))

(defn sarsa-agent
  "SARSA: as Q-learning, but the target is the Q of the action actually
  chosen next, not the best one -- learning the policy it follows. That
  makes it need exploration that never stops: with the book's
  exploration function alone, once every pair has had its `:ne` tries the
  policy it follows is greedy, and in the 4 x 3 world it settles, more
  often than not, on wandering at -0.4. Give it `:epsilon`."
  ([mdp] (sarsa-agent mdp {}))
  ([mdp opts] (q-agent mdp (update opts :rng #(or % (random/rng 1))) true)))

(defn greedy-policy
  "The action Q values most in each state of `mdp`."
  [{:keys [states actions]} Q]
  (into {} (for [s states :let [as (actions s)]]
             [s (when (seq as) (apply max-key #(get Q [s %] 0.0) as))])))

;; ---------------------------------------------------------------------------
;; Safe exploration

(defn policy-utility
  "U_h^pi: the expected utility of `pi` in the model `mdp`, averaged over
  every start state (but the terminal ones)."
  [{:keys [states actions] :as mdp} pi]
  (let [U (mdp/policy-evaluation mdp pi (zipmap states (repeat 0.0)) 1000)
        starts (filter #(seq (actions %)) states)]
    (/ (reduce + (map U starts)) (count starts))))

(defn bayesian-policy
  "Section 22.3.2's Bayesian choice for an agent that has stopped
  learning: of `policies`, the one with the highest expected utility over
  the `models` it might be in, weighted by their posterior probabilities
  -- argmax over pi of sum over h of P(h | e) U_h^pi. (Over every policy,
  the chapter says, this can be found only in special cases; here the
  candidates are given.)"
  [models posterior policies]
  (apply max-key (fn [pi] (reduce + (map #(* %2 (policy-utility %1 pi)) models posterior))) policies))

(defn robust-policy
  "The robust control choice: of `policies`, the one whose worst utility
  over the `models` is best -- argmax over pi of min over h of U_h^pi,
  a game against an adversary choosing the model."
  [models policies]
  (apply max-key (fn [pi] (apply min (map #(policy-utility % pi) models))) policies))

;; ---------------------------------------------------------------------------
;; Function approximation

(defn linear-utility
  "U^theta(s) = theta . f(s), equation 22.9's form, `features` giving the
  vector f(s)."
  [theta features s]
  (lin/dot theta (features s)))

(defn widrow-hoff
  "The Widrow-Hoff rule -- the delta rule, for online least squares --
  equation 22.10 for a linear approximator: each theta_i moved by alpha
  (u - U^theta) f_i, `f` the features of the state and `u` the value
  observed for it."
  [theta f u alpha]
  (let [err (- u (lin/dot theta f))]
    (mapv (fn [t fi] (+ t (* alpha err fi))) theta f)))

(defn approximate-direct-utility-agent
  "Direct utility estimation with a linear function approximator, section
  22.4.1, for the fixed policy `pi`: a trial's rewards remembered, and
  when it ends each state's reward-to-go a training example, theta moved
  toward each by the Widrow-Hoff rule. Changing theta for one state
  changes U^theta for every other: the agent generalizes."
  [{:keys [actions gamma]} pi features theta alpha]
  [(fn [{:keys [trial theta] :as m} [s r]]
     (let [trial (conj trial [s r])]
       (if (seq (actions s))
         [(pi s) (assoc m :trial trial)]
         [nil {:trial []
               :theta (loop [[[s r] & more] (rseq trial) to-go 0.0 theta theta]
                        (if s
                          (let [g (+ r (* gamma to-go))] (recur more g (widrow-hoff theta (features s) g alpha)))
                          theta))}])))
   {:trial [] :theta theta}])

(defn approximate-td-agent
  "Passive TD learning with a linear function approximator, equation
  22.11, for the fixed policy `pi`: after each move theta_i += alpha [r +
  gamma U^theta(s') - U^theta(s)] f_i(s). A terminal state's utility is
  its reward, and theta is moved toward that too."
  [{:keys [actions gamma]} pi features theta alpha]
  [(fn [{:keys [theta s r]} [s' r']]
     (let [terminal? (empty? (actions s'))
           u' (if terminal? r' (lin/dot theta (features s')))
           theta (if s (widrow-hoff theta (features s) (+ r (* gamma u')) alpha) theta)
           theta (if terminal? (widrow-hoff theta (features s') r' alpha) theta)]
       [(when-not terminal? (pi s'))
        {:theta theta :s (when-not terminal? s') :r (when-not terminal? r')}]))
   {:theta theta :s nil :r nil}])

(defn linear-q
  "A linear Q-function, Q^theta(s, a) = theta . f(s, a), for
  `approximate-q-agent`: `{:q :learn}`, `learn` the Widrow-Hoff rule
  toward a target."
  [features]
  {:q (fn [theta s a] (lin/dot theta (features s a)))
   :learn (fn [theta s a target alpha] (widrow-hoff theta (features s a) target alpha))})

(defn network-q
  "A deep Q-function, section 22.4.3: a network from the `features` of a
  state to a Q value for each of `actions`, moved toward a target by
  back-propagation -- the target the output for the action taken, the
  other outputs their own values, so that only that output's error is
  propagated. (The network's loss is the squared error, not half of it,
  so its step is twice alpha's.)"
  [features actions]
  (let [index (zipmap actions (range))]
    {:q (fn [net s a] (nth (nn/predict net (features s)) (index a)))
     :learn (fn [net s a target alpha]
              (let [x (features s)]
                (nn/sgd-step net x (assoc (nn/predict net x) (index a) target) alpha :squared)))}))

(defn approximate-q-agent
  "Q-learning with a function approximator, equation 22.12: after each
  move the parameters moved by alpha [R + gamma max over a' of Q(s', a')
  - Q(s, a)] times the gradient of Q(s, a) -- `approximator` from
  `linear-q` or `network-q`, `params` its parameters to start. Actions
  are epsilon-greedy (`:epsilon`, default 0.1): the exploration
  function's counts are per state-action pair, which a function
  approximator is meant to do without. With `:replay k`, experience
  replay: the last `:buffer` (default 10000) transitions are kept, and
  after each move `k` of them, drawn at random, learned from again --
  the remedy for catastrophic forgetting."
  [{:keys [actions gamma]} {:keys [q learn]} params {:keys [alpha epsilon rng replay buffer]
                                                     :or {alpha 0.05 epsilon 0.1 replay 0 buffer 10000}}]
  (let [rng (or rng (random/rng 1))
        learn-from (fn [params [s a r s' terminal? r']]
                     (learn params s a (+ r (* gamma (if terminal? r' (apply max (map #(q params s' %) (actions s')))))) alpha))]
    [(fn [{:keys [params s a r transitions]} [s' r']]
       (let [terminal? (empty? (actions s'))
             t (when s [s a r s' terminal? r'])
             params (if t (learn-from params t) params)
             transitions (cond-> transitions
                           (and t (>= (count transitions) buffer)) (subvec 1)
                           t (conj t))
             params (if (and (pos? replay) (seq transitions))
                      (reduce learn-from params (repeatedly replay #(random/pick rng transitions)))
                      params)
             a' (when-not terminal?
                  (let [as (vec (actions s'))]
                    (if (< (double (rng)) epsilon)
                      (random/pick rng as)
                      (let [scored (map (fn [a] [a (q params s' a)]) as)
                            top (apply max (map second scored))]
                        (random/pick rng (vec (keep (fn [[a v]] (when (>= v (- top 1e-12)) a)) scored)))))))]
         [a' {:params params :s (when-not terminal? s') :a a' :r (when-not terminal? r') :transitions transitions}]))
     {:params params :s nil :a nil :r nil :transitions []}]))

(defn approximate-greedy-policy
  "The action the approximate Q values most in each state of `mdp`."
  [{:keys [states actions]} {:keys [q]} params]
  (into {} (for [s states :let [as (actions s)]]
             [s (when (seq as) (apply max-key #(q params s %) as))])))

;; ---------------------------------------------------------------------------
;; Policy search

(defn softmax-policy
  "The policy that `theta` -- a preference for each `[s a]`, default 0 --
  stands for: in s, action a with probability proportional to
  e^theta(s, a)."
  [{:keys [actions]} theta s]
  (let [as (actions s)
        prefs (map #(get theta [s %] 0.0) as)
        top (if (seq prefs) (apply max prefs) 0.0)
        ws (map #(math/exp (- % top)) prefs)
        z (reduce + ws)]
    (into {} (map (fn [a w] [a (/ w z)]) as ws))))

(defn- reinforce-step
  "One trial under the softmax policy of `theta`, and each preference
  moved along the return times the log probability's gradient."
  [{:keys [actions transitions reward gamma] :as mdp} start alpha rng theta]
  (let [episode (loop [s start out [] k 0]
                  (if (or (empty? (actions s)) (> k 200))
                    (conj out [s nil (reward s)])
                    (let [a (random/weighted rng (softmax-policy mdp theta s))
                          s' (sample rng (transitions s a))]
                      (recur s' (conj out [s a (reward s)]) (inc k)))))
        returns (reverse (reductions (fn [g [_ _ r]] (+ r (* gamma g))) 0.0 (reverse episode)))]
    (reduce (fn [theta [[s a] g]]
              (if (nil? a)
                theta
                (let [pi (softmax-policy mdp theta s)]
                  (reduce (fn [theta b]
                            (update theta [s b] (fnil + 0.0) (* alpha g (- (if (= a b) 1.0 0.0) (pi b)))))
                          theta (actions s)))))
            theta (map vector episode returns))))

(defn reinforce-seq
  "Policy search by REINFORCE, which the chapter describes without a
  figure: each step one trial under the softmax policy, then each
  preference moved by `alpha` times the return that followed each choice
  times the gradient of the log probability of the choice --
  d log pi(a | s) / d theta(s, b) = [a = b] - pi(b | s). Each state is the
  preferences `theta`."
  [mdp start alpha rng]
  (iterate (partial reinforce-step mdp start alpha rng) {}))
