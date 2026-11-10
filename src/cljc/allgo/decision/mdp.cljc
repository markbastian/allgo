(ns allgo.decision.mdp
  "Sequential decisions: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 17, as steps.

    value iteration    VALUE-ITERATION (figure 17.6): a sweep of Bellman
                       updates a step, until no utility changes by more
                       than epsilon (1 - gamma) / gamma
    policy iteration   POLICY-ITERATION (figure 17.9): a step an
                       evaluation of the policy -- `k` sweeps of it, as in
                       modified policy iteration -- and an improvement,
                       until no state's action changes
    POMDP value        POMDP-VALUE-ITERATION (figure 17.16): utility as a
    iteration          set of conditional plans, each an alpha vector; a
                       step the plans one action deeper, those no belief
                       prefers removed -- checked over a fine lattice of
                       beliefs, see `prune`

  An MDP is a map:

    :states       the states
    :actions      (fn [s]) -> the actions possible there; a terminal
                  state has none
    :transitions  (fn [s a]) -> `[[p s'] ...]`
    :reward       (fn [s]) -> R(s), the reward for being in s
    :gamma        the discount

  `grid-mdp` makes one of the chapter's grid worlds; `four-by-three` is
  its 4 x 3 world as aima-python (MIT) states it."
  (:require [allgo.numerics.linear :as lin]
            [allgo.search :as search]))

;; ---------------------------------------------------------------------------
;; Grid worlds

(def ^:private headings [[1 0] [0 1] [-1 0] [0 -1]])
(def ^:private heading-index (zipmap headings (range)))
(defn- turn [a by] (headings (mod (+ (long (heading-index a)) (long by)) 4)))

(defn grid-mdp
  "A grid world from `rows`, top row first, each cell a reward or nil for
  a wall; `terminals` are `[x y]` cells, `[0 0]` the bottom left. An
  action is a heading `[dx dy]`; it goes that way with probability 0.8 and
  veers to either side with 0.1 each, and a move into a wall or off the
  grid stays put."
  [rows terminals gamma]
  (let [h (count rows)
        cells (into {} (for [[r row] (map-indexed vector rows) [x v] (map-indexed vector row) :when (some? v)]
                         [[x (- h 1 r)] v]))
        terminal? (set terminals)
        go (fn [[x y] [dx dy]] (let [c [(+ x dx) (+ y dy)]] (if (contains? cells c) c [x y])))]
    {:states (vec (sort (keys cells)))
     :terminals terminal?
     :reward cells
     :gamma gamma
     :actions (fn [s] (if (terminal? s) [] headings))
     :transitions (fn [s a]
                    (if (nil? a)
                      [[1.0 s]]
                      [[0.8 (go s a)] [0.1 (go s (turn a -1))] [0.1 (go s (turn a 1))]]))}))

(def four-by-three
  "The chapter's 4 x 3 world as aima-python's `sequential_decision_environment`
  has it: -0.04 a step, +1 and -1 at the terminals, a wall at [1 1], and
  gamma 0.9."
  (grid-mdp [[-0.04 -0.04 -0.04 1.0]
             [-0.04 nil -0.04 -1.0]
             [-0.04 -0.04 -0.04 -0.04]]
            [[3 2] [3 1]]
            0.9))

;; ---------------------------------------------------------------------------
;; Value iteration

(defn q-value
  "The expected utility of doing `a` in `s` and acting on `U` from then
  on, rewards counted as R(s): sum over s' of P(s' | s, a) U(s')."
  [{:keys [transitions]} s a U]
  (reduce + (map (fn [[p s']] (* p (U s'))) (transitions s a))))

(defn- bellman [{:keys [actions reward gamma] :as mdp} s U]
  (let [as (actions s)]
    (+ (reward s) (if (seq as) (* gamma (apply max (map #(q-value mdp s % U) as))) 0.0))))

(defn- value-iteration-step
  "A sweep of Bellman updates, until the largest change is under `bound`."
  [{:keys [states] :as mdp} bound {:keys [U delta sweeps]}]
  (when (or (nil? delta) (> delta bound))
    (let [U' (into {} (map (fn [s] [s (bellman mdp s U)]) states))]
      {:U U' :previous U :sweeps (inc (long sweeps))
       :delta (apply max (map #(abs (- (U' %) (U %))) states))})))

(defn value-iteration-seq
  "Every state of value iteration on `mdp`: `{:U :delta :sweeps}`, a sweep
  of Bellman updates a step from utilities of zero, ending when the
  largest change is under epsilon (1 - gamma) / gamma -- then the
  utilities are within epsilon of the true ones."
  [{:keys [states gamma] :as mdp} epsilon]
  (search/step-seq (partial value-iteration-step mdp (/ (* epsilon (- 1.0 gamma)) gamma))
                   {:U (zipmap states (repeat 0.0)) :delta nil :sweeps 0}))

(defn value-iteration
  "The utilities of `mdp`'s states, within `epsilon` of the true ones:
  as the figure has it, U, the utilities the last sweep started from, not
  the U' it made."
  [mdp epsilon]
  (:previous (last (value-iteration-seq mdp epsilon))))

(defn best-policy
  "The action in each state that maximizes expected utility under `U`; nil
  in a terminal state."
  [{:keys [states actions] :as mdp} U]
  (into {} (for [s states]
             [s (let [as (actions s)] (when (seq as) (apply max-key #(q-value mdp s % U) as)))])))

;; ---------------------------------------------------------------------------
;; Policy iteration

(defn policy-evaluation
  "`U` updated `k` times by the Bellman equation for the fixed policy `pi`:
  modified policy iteration's approximate evaluation."
  [{:keys [states reward gamma] :as mdp} pi U k]
  (loop [U U i 0]
    (if (= i k)
      U
      (recur (into {} (map (fn [s] [s (+ (reward s) (if (pi s) (* gamma (q-value mdp s (pi s) U)) 0.0))]) states))
             (inc i)))))

(defn- policy-iteration-step
  "An evaluation of `k` sweeps and an improvement, until no action
  changes."
  [{:keys [states actions] :as mdp} k {:keys [pi U changed] :as st}]
  (when (pos? (long changed))
    (let [U (policy-evaluation mdp pi U k)
          [pi changed] (reduce (fn [[pi n] s]
                                 (let [as (actions s)]
                                   (if (empty? as)
                                     [pi n]
                                     (let [best (apply max-key #(q-value mdp s % U) as)]
                                       (if (> (q-value mdp s best U) (q-value mdp s (pi s) U))
                                         [(assoc pi s best) (inc n)]
                                         [pi n])))))
                               [pi 0] states)]
      {:pi pi :U U :changed changed :iterations (inc (long (:iterations st)))})))

(defn policy-iteration-seq
  "Every state of policy iteration on `mdp`: `{:pi :U :changed}`, a step an
  evaluation (`k` sweeps, default 20) and an improvement, from the first
  action in every state, until no action changes."
  ([mdp] (policy-iteration-seq mdp 20))
  ([{:keys [states actions] :as mdp} k]
   (search/step-seq (partial policy-iteration-step mdp k)
                    {:pi (into {} (for [s states] [s (first (actions s))]))
                     :U (zipmap states (repeat 0.0)) :changed 1 :iterations 0})))

(defn policy-iteration
  "The optimal policy of `mdp`."
  [mdp]
  (:pi (last (policy-iteration-seq mdp))))

;; ---------------------------------------------------------------------------
;; POMDPs

;; A POMDP here, as aima-python states its examples: states 0..n-1,
;; `:transition` a [a][s][s'] array of P(s' | s, a), `:sensor` a [a][s'][e]
;; array of P(e | s', a), `:reward` a [a][s] array of R(s, a), and
;; `:gamma`. An alpha vector is a plan's expected utility from each state.

(defn update-belief
  "POMDP filtering: the belief after doing `a` in `belief` and perceiving
  `e`, b'(s') = alpha P(e | s') sum over s of P(s' | s, a) b(s)."
  [{:keys [transition sensor]} belief a e]
  (let [n (count belief)
        predicted (mapv (fn [s'] (reduce + (map (fn [s] (* (belief s) (get-in transition [a s s']))) (range n)))) (range n))
        updated (mapv (fn [s'] (* (predicted s') (get-in sensor [a s' e]))) (range n))
        total (reduce + updated)]
    (if (pos? total) (mapv #(/ % total) updated) predicted)))

(defn- lattice
  "Every belief over `n` states whose probabilities are multiples of
  1/`r`: the corners of the simplex, its midpoints, and a fine net
  between."
  [n r]
  (letfn [(parts [k left]
            (if (= k 1)
              [[left]]
              (for [i (range (inc left)) rest (parts (dec k) (- left i))] (cons i rest))))]
    (mapv (fn [p] (mapv #(/ (double %) r) p)) (parts n r))))

(defn prune
  "REMOVE-DOMINATED-PLANS, approximately: the plans, `{:alpha :action
  ...}`, that are the best at some belief of a lattice of spacing 1/`r`
  (default 200). Exact pruning asks a linear program of each plan whether
  any belief prefers it; over hundreds of plans that was slow here, and
  failed, so this checks a net of beliefs instead -- exact up to the
  spacing for two states, an approximation beyond."
  ([plans] (prune plans 200))
  ([plans r]
   (let [plans (mapv first (vals (group-by :alpha plans)))
         n (count (:alpha (first plans)))
         r (if (> n 3) (max 4 (quot 60 n)) r)
         best (fn [b] (apply max-key #(lin/dot (:alpha %) b) plans))]
     (vec (distinct (map best (lattice n r)))))))

(defn- deeper
  "Every plan one step deeper: an action, then for each percept one of
  `plans` -- alpha(s) = R(s, a) + gamma sum over s' and e of P(s' | s, a)
  P(e | s', a) alpha_e(s')."
  [{:keys [transition sensor reward gamma]} plans]
  (let [n (count (first reward))
        percepts (count (get-in sensor [0 0]))]
    (for [a (range (count reward))
          choice (reduce (fn [acc _] (for [c acc p plans] (conj c p))) [[]] (range percepts))]
      {:action a
       :then (mapv :id choice)
       :alpha (mapv (fn [s] (+ (get-in reward [a s])
                               (* gamma (reduce + (for [s' (range n) e (range percepts)]
                                                    (* (get-in transition [a s s']) (get-in sensor [a s' e])
                                                       ((:alpha (choice e)) s')))))))
                    (range n))})))

(defn pomdp-value
  "The utility of `belief` under a set of alpha-vector plans."
  [plans belief]
  (apply max (map #(lin/dot (:alpha %) belief) plans)))

(defn- pomdp-step
  "The plans a step deeper, pruned, until the utility at every probe
  belief changes by less than `bound`."
  [pomdp probes bound {:keys [plans depth diff]}]
  (when (or (nil? diff) (>= diff bound))
    (let [ids (map-indexed (fn [i p] (assoc p :id i)) plans)
          plans' (prune (deeper pomdp (vec ids)))]
      {:plans (mapv #(dissoc % :id) plans') :depth (inc (long depth))
       :diff (apply max (map #(abs (- (pomdp-value plans' %) (pomdp-value plans %))) probes))})))

(defn pomdp-value-iteration-seq
  "Every state of POMDP value iteration: `{:plans :depth}`, the
  undominated conditional plans of each depth, from the empty plan of
  utility zero, until the utility at every belief changes by less than
  epsilon (1 - gamma) / gamma -- measured as the largest change of the
  best alpha at each corner and midpoint of the belief simplex."
  [{:keys [reward gamma] :as pomdp} epsilon]
  (let [n (count (first reward))
        probes (concat (for [i (range n)] (assoc (vec (repeat n 0.0)) i 1.0)) [(vec (repeat n (/ 1.0 n)))])]
    (search/step-seq (partial pomdp-step pomdp probes (/ (* epsilon (- 1.0 gamma)) gamma))
                     {:plans [{:alpha (vec (repeat n 0.0)) :action nil}] :depth 0 :diff nil})))
