(ns allgo.robotics.humans
  "Robots among people: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 26.8.

    predicting     a person modeled as noisily rational: the better an
                   action for their goal, the likelier -- P(u | x, goal)
                   proportional to e^-Q(x, u; goal), equation 26.8, Q
                   here the step's length plus the straight-line distance
                   left -- so each step they take is evidence about the
                   goal, b'(goal) proportional to b(goal) P(u | x, goal),
                   equation 26.9, and the belief predicts their next steps
    imitating      a policy learned from demonstrations, by any supervised
                   learner: behavioral cloning learns from the expert's
                   own states only, so a small error takes the robot
                   where it has no examples and the errors compound;
                   DAgger (data aggregation) rolls its own policy out,
                   asks the expert what it should have done at the states
                   it actually reached, adds those to all the data so
                   far, and learns again

  Positions are `[x y]`."
  (:require [clojure.math :as math]))

(defn- dist [[ax ay] [bx by]] (math/hypot (- ax bx) (- ay by)))

(defn step-likelihoods
  "P(u | x, goal) for each of `actions`, displacements `[dx dy]`: a
  softmax of -`beta` Q(x, u; goal), Q the step's length plus the distance
  from where it lands to the goal. Larger beta, more nearly optimal."
  [x goal actions beta]
  (let [q (fn [[dx dy :as u]] (+ (dist [0 0] u) (dist [(+ (x 0) dx) (+ (x 1) dy)] goal)))
        ws (map #(math/exp (- (* beta (q %)))) actions)
        z (reduce + ws)]
    (zipmap actions (map #(/ % z) ws))))

(defn- nearest-action [actions u]
  (apply min-key #(dist % u) actions))

(defn belief-seq
  "The belief over `goals` (a map of goal to position) as the person's
  `steps` -- `[x u]`, where they were and how they moved -- come in, from
  `prior` (default uniform): each step weighs each goal by how likely it
  made that step, the step taken as the nearest of `actions`."
  ([goals steps actions beta]
   (belief-seq goals (zipmap (keys goals) (repeat (/ 1.0 (count goals)))) steps actions beta))
  ([goals prior steps actions beta]
   (reductions (fn [b [x u]]
                 (let [a (nearest-action actions u)
                       w (into {} (map (fn [[g p]] [g (* p ((step-likelihoods x (goals g) actions beta) a))]) b))
                       z (reduce + (vals w))]
                   (into {} (map (fn [[g p]] [g (/ p z)]) w))))
               prior steps)))

(defn predict-step
  "The distribution over the person's next step from `x` under `belief`:
  each goal's step likelihoods, weighted by the belief in it."
  [goals belief x actions beta]
  (apply merge-with + (map (fn [[g p]] (update-vals (step-likelihoods x (goals g) actions beta) #(* p %))) belief)))

;; ---------------------------------------------------------------------------
;; Imitation

(defn rollout
  "The states from `x0` under `policy` (fn [x]) -> u, by `dynamics` (fn [x
  u]) -> x', for `n` steps."
  [dynamics policy x0 n]
  (vec (take (inc n) (iterate #(dynamics % (policy %)) x0))))

(defn behavioral-cloning
  "A policy fitted by `fit` (fn [states actions]) -> policy to the
  expert's own demonstrations: rollouts of `expert` from `starts`."
  [dynamics expert fit starts n]
  (let [xs (vec (mapcat #(rollout dynamics expert % n) starts))]
    (fit xs (mapv expert xs))))

(defn dagger-seq
  "DAgger, an iteration a step: from the expert's demonstrations from
  `starts`, each step rolls the current policy out from the same starts,
  labels every state it reached with the `expert`'s action there, adds
  them to all the data so far, and fits the next policy with `fit`. Each
  state `{:policy :states :actions}`."
  [dynamics expert fit starts n]
  (let [xs0 (vec (mapcat #(rollout dynamics expert % n) starts))
        init {:states xs0 :actions (mapv expert xs0) :policy (fit xs0 (mapv expert xs0))}]
    (iterate (fn [{:keys [policy states actions]}]
               (let [visited (vec (mapcat #(rollout dynamics policy % n) starts))
                     states (into states visited)
                     actions (into actions (map expert visited))]
                 {:states states :actions actions :policy (fit states actions)}))
             init)))

(defn nearest-neighbor
  "The 1-nearest-neighbor policy of `states` labeled `actions`: the action
  of the nearest state seen -- exact where it has examples, and no better
  than its nearest example anywhere else."
  [states actions]
  (let [pairs (mapv vector states actions)]
    (fn [x] (second (apply min-key (fn [[s]] (reduce + (map #(let [d (- %1 %2)] (* d d)) s x))) pairs)))))
