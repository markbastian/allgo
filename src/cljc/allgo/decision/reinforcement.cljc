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

  An agent is `(fn [memory percept]) -> [action memory']`, a percept
  `[state reward]`, as for every agent in this library; `trial-seq` runs
  one through an MDP of `allgo.decision.mdp` until a terminal state, a
  step a move, and `trials` runs it again and again, its memory carried
  over. The learning rate and the exploration function are aima-python's
  (MIT): alpha(n) = 60 / (59 + n), and f(u, n) = R+ until a state-action
  pair has been tried Ne times."
  (:require [allgo.decision.mdp :as mdp]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn alpha [n] (/ 60.0 (+ 59.0 n)))

;; ---------------------------------------------------------------------------
;; Running trials

(defn- sample [rng outcomes]
  (let [r (double (rng))]
    (loop [[[p s] & more] outcomes acc 0.0]
      (if (or (empty? more) (< r (+ acc p))) s (recur more (+ acc p))))))

(defn trial-seq
  "One trial of `agent`, a `[program memory]` pair, in `mdp` from `start`:
  each state `{:state :reward :action :memory}`, ending after the agent
  perceives a terminal state, or after `max-moves` (default 10000) -- an
  agent can learn to like a wall."
  ([mdp agent start rng] (trial-seq mdp agent start rng 10000))
  ([{:keys [actions transitions reward]} [program memory] start rng max-moves]
   (let [[a m] (program memory [start (reward start)])]
     (take-while some?
                 (iterate (fn [{:keys [state action] :as st}]
                            (when (and st (seq (actions state)) action (< (long (:moves st)) (long max-moves)))
                              (let [s' (sample rng (transitions state action))
                                    [a m] (program (:memory st) [s' (reward s')])]
                                {:state s' :reward (reward s') :action a :memory m :moves (inc (long (:moves st)))})))
                          {:state start :reward (reward start) :action a :memory m :moves 0})))))

(defn trials
  "The agent's memory after each of `n` trials from `start`, memory carried
  from one to the next."
  [mdp [program memory] start n rng]
  (take (inc n)
        (iterate (fn [m] (:memory (last (trial-seq mdp [program m] start rng))))
                 memory)))

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
  [{:keys [actions gamma]} {:keys [ne r+ rng epsilon] :or {ne 5 r+ 2.0 epsilon 0.0}} on-policy?]
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
                            target (if on-policy?
                                     (q s' a')
                                     (apply max (map #(q s' %) as')))]
                        [(assoc Q [s a] (+ (q s a) (* (alpha (n-sa [s a])) (+ r (* gamma target) (- (q s a))))))
                         n-sa])
                      [Q n-sa])]
       [a' {:Q Q :n-sa n-sa :s (when-not terminal? s') :a a' :r (when-not terminal? r')}]))
   {:Q {} :n-sa {} :s nil :a nil :r nil}])

(defn q-learning-agent
  "Q-LEARNING-AGENT: Q(s, a) += alpha(N(s, a)) (r + gamma max over a' of
  Q(s', a') - Q(s, a)), and the next action the one the exploration
  function f(Q, N) favors -- `:ne` (default 5) tries of each pair before
  trusting its Q, valued `:r+` (default 2) until then -- or, with
  probability `:epsilon` (default 0), any action."
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

(defn reinforce-seq
  "Policy search by REINFORCE, which the chapter describes without a
  figure: each step one trial under the softmax policy, then each
  preference moved by `alpha` times the return that followed each choice
  times the gradient of the log probability of the choice --
  d log pi(a | s) / d theta(s, b) = [a = b] - pi(b | s). Each state is the
  preferences `theta`."
  [{:keys [actions transitions reward gamma] :as mdp} start alpha rng]
  (let [choose (fn [dist] (let [r (double (rng))]
                            (loop [[[a p] & more] (seq dist) acc 0.0]
                              (if (or (empty? more) (< r (+ acc p))) a (recur more (+ acc p))))))]
    (iterate (fn [theta]
               (let [episode (loop [s start out [] k 0]
                               (if (or (empty? (actions s)) (> k 200))
                                 (conj out [s nil (reward s)])
                                 (let [a (choose (softmax-policy mdp theta s))
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
             {})))
