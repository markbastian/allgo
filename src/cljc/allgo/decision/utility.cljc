(ns allgo.decision.utility
  "Decisions under uncertainty: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapters 12 and 16.

    DT-AGENT                 figure 12.1: keep a belief state, and act to
                             maximize expected utility under it
    value of information     the expected gain in utility from learning
                             an evidence variable before deciding
    INFORMATION-GATHERING    figure 16.9: ask for the evidence whose value
    -AGENT                   per unit cost is greatest, if it is worth
                             its cost; otherwise act

  A decision problem here is over a discrete hidden state: `:prior`, a
  map from states to probabilities; `:actions`; and `:utility`, (fn [a
  s]) -> U. An evidence variable is `{:name :cost :likelihood}`, its
  likelihood (fn [e s]) -> P(e | s) over values `:values`.")

(defn expected-utility [utility belief a]
  (reduce + (map (fn [[s p]] (* p (utility a s))) belief)))

(defn best-action
  "The action of greatest expected utility under `belief`, and that
  utility: `[a eu]`."
  [{:keys [actions utility]} belief]
  (let [a (apply max-key #(expected-utility utility belief %) actions)]
    [a (expected-utility utility belief a)]))

(defn- posterior [belief likelihood e]
  (let [w (into {} (map (fn [[s p]] [s (* p (likelihood e s))]) belief))
        total (reduce + (vals w))]
    [total (if (pos? total) (into {} (map (fn [[s v]] [s (/ v total)]) w)) belief)]))

(defn vpi
  "The value of perfect information about `evidence`: the expected utility
  of deciding after seeing it, averaged over what it might show, minus
  that of deciding now. Never negative."
  [problem belief {:keys [values likelihood]}]
  (let [[_ now] (best-action problem belief)
        later (reduce + (for [e values
                              :let [[pe post] (posterior belief likelihood e)]
                              :when (pos? pe)]
                          (* pe (second (best-action problem post)))))]
    (- later now)))

(defn information-gathering-agent
  "INFORMATION-GATHERING-AGENT for `problem` and the `evidence` it may ask
  for, as `(fn [memory percept]) -> [action memory']`: a percept is nil or
  `[name value]`, an observation it asked for. Its action is `[:request
  name]` while some evidence is worth more than it costs, else `[:act a]`."
  [problem evidence]
  [(fn [{:keys [belief asked]} percept]
     (let [belief (if-let [[n v] percept]
                    (second (posterior belief (:likelihood (first (filter #(= n (:name %)) evidence))) v))
                    belief)
           asked (if percept (conj asked (first percept)) asked)
           open (remove #(asked (:name %)) evidence)
           scored (map (fn [ev] [ev (vpi problem belief ev)]) open)
           [best value] (when (seq scored) (apply max-key (fn [[ev v]] (/ v (max 1e-12 (:cost ev)))) scored))]
       (if (and best (> value (:cost best)))
         [[:request (:name best)] {:belief belief :asked asked}]
         [[:act (first (best-action problem belief))] {:belief belief :asked asked}])))
   {:belief (:prior problem) :asked #{}}])

(defn dt-agent
  "DT-AGENT: `update` folds the last action and a percept into the belief
  state, and the action is the one of highest expected utility under it."
  [problem update]
  [(fn [{:keys [belief action]} percept]
     (let [belief (update belief action percept)
           [a] (best-action problem belief)]
       [a {:belief belief :action a}]))
   {:belief (:prior problem) :action nil}])
