(ns allgo.planning.satplan
  "SATPlan: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 7 (SATPLAN, figure 7.22). A planning
  problem translated to propositional logic for plans of 0, 1, 2, ...
  steps in turn, each translation handed to DPLL, until one is
  satisfiable; the plan is read off the action symbols true in the
  model.

  A problem is a map of propositional STRIPS-like parts:

    :fluents   the facts that may change, as keywords
    :init      the set of fluents true at the start; the rest are false
    :goal      literals that must hold at the end: a fluent, or
               `[:not fluent]`
    :actions   `{:name :pre :effect}`, `:pre` and `:effect` sets of
               literals

  The translation for a plan of `t` steps, with a symbol for each fluent
  and each action at each step, `Have@0`, `Eat@0`, `Have@1` ...:

    the initial state at step 0, and the goal at step `t`;
    each action at step `i` implies its preconditions at step `i`;
    successor-state axioms, one per fluent and step: true at `i + 1` if
      an action at `i` made it true, or it was true at `i` and no action
      made it false;
    and no two actions at the same step -- the book's simplest exclusion
      axiom, which makes plans totally ordered."
  (:require [allgo.logic.propositional :as pl]))

(defn- at [k i] (keyword (str (name k) "@" i)))

(defn- lit-at [lit i]
  (if (keyword? lit) (at lit i) [:not (at (second lit) i)]))

(defn- adds? [{:keys [effect]} f] (contains? effect f))
(defn- deletes? [{:keys [effect]} f] (contains? effect [:not f]))

(defn translate
  "TRANSLATE-TO-SAT: the clauses whose models are the plans of exactly
  `t` steps."
  [{:keys [fluents init goal actions]} t]
  (let [steps (range t)
        sentences
        (concat
         (for [f fluents] (if (contains? init f) (at f 0) [:not (at f 0)]))
         (for [g goal] (lit-at g t))
         (for [i steps {:keys [name pre]} actions p pre]
           [:implies (at name i) (lit-at p i)])
         (for [i steps f fluents
               :let [adders (filter #(adds? % f) actions)
                     deleters (filter #(deletes? % f) actions)]]
           [:iff (at f (inc i))
            (into [:or]
                  (concat (map #(at (:name %) i) adders)
                          [(into [:and (at f i)] (map #(vector :not (at (:name %) i)) deleters))]))])
         (for [i steps
               a (range (count actions)) b (range (inc a) (count actions))]
           [:or [:not (at (:name (actions a)) i)] [:not (at (:name (actions b)) i)]]))]
    (apply pl/cnf sentences)))

(defn- extract
  "EXTRACT-SOLUTION: the actions true in `model`, in step order."
  [{:keys [actions]} model t]
  (vec (for [i (range t) {:keys [name]} actions
             :when (true? (get model (at name i)))]
         name)))

(defn satplan-seq
  "Every state of SATPlan on `problem`: DPLL's machine states for the plan
  of 0 steps, then 1, ... up to `t-max`, each tagged with its `:t`, until
  one is satisfiable. Its last state has the `:plan`, or none if no plan
  of `t-max` steps or fewer exists."
  [problem t-max]
  (let [problem (update problem :actions vec)]
    (letfn [(tries [t]
              (lazy-seq
               (when (<= t t-max)
                 (let [states (map #(assoc % :t t) (pl/dpll-seq (translate problem t)))
                       end (last states)]
                   (if (true? (:value end))
                     (concat (butlast states)
                             [(assoc end :plan (extract problem (:model (:context end)) t))])
                     (concat states (tries (inc t))))))))]
      (tries 0))))

(defn satplan
  "The shortest plan for `problem` of `t-max` steps or fewer, or nil."
  [problem t-max]
  (:plan (last (satplan-seq problem t-max))))
