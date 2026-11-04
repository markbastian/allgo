(ns allgo.decision.hierarchical
  "Hierarchical reinforcement learning: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 22.4.5.

  A partial program outlines the agent's behavior and leaves choices
  open; learning fills them in. Here a partial program is a function of
  the physical state s and the machine state m -- whatever the program
  needs to know where it is: a routine and its arguments, a call stack --
  returning what the program does next:

    {:choose cs :then f}   a choice point: one of `cs` to be chosen, and
                           (f c) the machine state after choosing c
    {:act a :then m'}      a primitive action, done in the world
    {:then m'}             a computational step, the machine state
                           changing and nothing else

  The choice states sigma = (s, m) of the joint state space and the
  choices there make a decision problem -- semi-Markov, since any number
  of physical actions may come between one choice and the next -- which
  `hrl-agent` learns by Q-learning: Q(sigma, c) moved toward rho + gamma^k
  max over c' of Q(sigma', c'), rho the discounted rewards of the k moves
  between the two choice states (the extra discount bookkeeping the
  chapter's footnote mentions). The policy it learns is the best
  consistent with the program. `flat-program` is the trivial program --
  while true do choose(A(s)) -- with which this is ordinary Q-learning.

  The agent is `(fn [memory percept]) -> [action memory']` as in
  `allgo.decision.reinforcement`, and runs with its `trial-seq`."
  (:require [allgo.decision.reinforcement :as rl]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn flat-program
  "while true do choose(A(s)): HRL with this program is ordinary
  Q-learning. Its machine states are `[:top]` and `[:do a]`."
  [{:keys [actions]}]
  (fn [s [routine a]]
    (case routine
      :top {:choose (actions s) :then (fn [a] [:do a])}
      :do {:act a :then [:top]})))

(defn- next-node
  "The program run from `s` and `m` through its computational steps to
  the next choice or action: `[node m]`."
  [program s m]
  (loop [m m steps 0]
    (let [node (program s m)]
      (if (or (:choose node) (contains? node :act))
        [node m]
        (if (< steps 10000)
          (recur (:then node) (inc steps))
          (throw (ex-info "the partial program neither chooses nor acts" {:s s :m m})))))))

(defn hrl-agent
  "Hierarchical Q-learning with the partial `program` from machine state
  `m0`: at each choice state the last choice's Q moved toward what it led
  to, by alpha(N) as `allgo.decision.reinforcement` has it, and the next
  choice the best by Q, ties broken at random, or with probability
  `:epsilon` (default 0.1) any. With `:potential`, a function Phi of the
  joint state (s, m), the choices' rewards are shaped, gamma^k
  Phi(sigma') - Phi(sigma) added -- advice that may depend on the
  machine state, which no flat learner's shaping can (Phi of a terminal
  state taken as zero). With `:learn? false` it only acts."
  [{:keys [actions gamma]} program m0 {:keys [epsilon rng potential learn? Q]
                                       :or {epsilon 0.1 learn? true}}]
  (let [rng (or rng (random/rng 1))
        phi (fn [sigma] (if (and potential sigma) (potential (first sigma) (second sigma)) 0.0))
        settle (fn [{:keys [Q n pending] :as mem} sigma' v']
                 ;; the pending choice's update, toward rho + gamma^k v'
                 (if (and pending learn?)
                   (let [{:keys [sigma c rho k]} pending
                         key [sigma c]
                         n (update n key (fnil inc 0))
                         g (math/pow gamma k)
                         target (+ rho (- (* g (phi sigma')) (phi sigma)) (* g v'))
                         q (get Q key 0.0)]
                     (assoc mem :Q (assoc Q key (+ q (* (rl/alpha (n key)) (- target q)))) :n n))
                   mem))]
    [(fn [{:keys [m] :as mem} [s r]]
       (if (empty? (actions s))
         ;; a terminal state, worth its reward, k moves after the pending choice
         [nil (-> (settle mem nil r) (assoc :m m0 :pending nil))]
         (loop [mem mem m m]
           (let [[node m] (next-node program s m)]
             (if-let [cs (seq (:choose node))]
               (let [sigma [s m]
                     qs (map #(get (:Q mem) [sigma %] 0.0) cs)
                     mem (settle mem sigma (apply max qs))
                     c (if (< (double (rng)) epsilon)
                         (random/pick rng (vec cs))
                         (let [top (apply max qs)]
                           (random/pick rng (vec (keep (fn [[c q]] (when (>= q (- top 1e-12)) c)) (map vector cs qs))))))]
                 (recur (assoc mem :pending {:sigma sigma :c c :rho 0.0 :k 0}) ((:then node) c)))
               [(:act node)
                (-> mem
                    (update :pending #(some-> % (update :rho + (* (math/pow gamma (:k %)) r)) (update :k inc)))
                    (assoc :m (:then node)))])))))
     {:Q (or Q {}) :n {} :m m0 :pending nil}]))
