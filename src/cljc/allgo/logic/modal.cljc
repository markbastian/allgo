(ns allgo.logic.modal
  "Knowledge as modal logic over possible worlds: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 10
  (chapter 12 of the 3rd edition, whose text this follows).

  A model is a set of possible worlds, each with the atoms true in it, and
  for each agent an accessibility relation: the worlds consistent with
  what that agent knows in each world. `[:knows a p]` -- K_a P -- is true
  in a world when p is true in every world accessible to a from it, so
  knowledge nests: what Lois knows about what Superman knows. Besides the
  modal operator, sentences are as `allgo.logic.first-order` writes them:
  atoms keywords or `[:Spy :Bond]`, then `:not :and :or :implies :iff`,
  and `[:forall [x] s]`, `[:exists [x] s]` over the model's `:domain`, the
  same objects in every world -- so `[:exists [x] [:knows :Bond [:Spy
  x]]]` (someone Bond knows to be a spy) and `[:knows :Bond [:exists [x]
  [:Spy x]]]` (Bond knows there is a spy) come apart.

  Model checking is by labeling: quantifiers expanded over the domain,
  then every subsentence, innermost first, a step, labeled with the set
  of worlds where it holds -- a knowledge sentence's worlds those whose
  accessible worlds all carry its sentence's label.

  The frame conditions the chapter's axioms need are here too: K_a P =>
  P holds in every model whose accessibility is reflexive, K_a P => K_a
  K_a P in every transitive one."
  (:require [allgo.logic.first-order :as fo]))

(defn model
  "A model: `worlds` a map from each world to the atoms true in it,
  `access` a map from agent to a map from world to the worlds it can't
  tell from it, `domain` the objects quantifiers range over. With
  `:reflexive? true` (the default) every world is also accessible from
  itself, as the chapter's diagrams take for granted."
  [worlds access & {:keys [domain reflexive?] :or {reflexive? true}}]
  {:worlds (into {} (map (fn [[w atoms]] [w (set atoms)])) worlds)
   :domain (vec domain)
   :access (into {} (for [[a rel] access]
                      [a (into {} (for [w (keys worlds)]
                                    [w (cond-> (set (get rel w)) reflexive? (conj w))]))]))})

(defn- expand
  "`s` with its quantifiers replaced by the conjunction or disjunction of
  their instances over `domain`."
  [domain s]
  (if (vector? s)
    (case (first s)
      (:forall :exists)
      (let [[q vs body] s
            instances (for [objs (reduce (fn [acc _] (for [a acc o domain] (conj a o))) [[]] vs)]
                        (expand domain (fo/subst (zipmap vs objs) body)))]
        (into [(if (= q :forall) :and :or)] instances))
      :knows (let [[_ a p] s] [:knows a (expand domain p)])
      (:not :and :or :implies :iff) (into [(first s)] (map #(expand domain %) (rest s)))
      s)
    s))

(defn- atom? [s] (or (keyword? s) (and (vector? s) (not (#{:not :and :or :implies :iff :knows} (first s))))))

(defn- subsentences
  "Every subsentence of `s`, innermost first, each once."
  [s]
  (distinct (if (atom? s)
              [s]
              (concat (mapcat subsentences (if (= :knows (first s)) [(nth s 2)] (rest s))) [s]))))

(defn- label [{:keys [worlds access]} labels s]
  (let [ws (set (keys worlds))
        of #(labels %)]
    (if (atom? s)
      (into #{} (filter #(contains? (worlds %) s)) ws)
      (let [[op & args] s]
        (case op
          :not (reduce disj ws (of (first args)))
          :and (reduce (fn [acc a] (into #{} (filter (of a)) acc)) ws args)
          :or (into #{} (mapcat of) args)
          :implies (into #{} (filter #(or (not ((of (first args)) %)) ((of (second args)) %))) ws)
          :iff (into #{} (filter #(= (contains? (of (first args)) %) (contains? (of (second args)) %))) ws)
          :knows (let [[a p] args, rel (access a)]
                   (into #{} (filter #(every? (of p) (rel %))) ws)))))))

(defn check-seq
  "Model checking `s` in `m` by labeling, a step a subsentence, innermost
  first: each state `{:sentence :worlds :labels}`, `:worlds` where that
  step's subsentence holds and `:labels` every label so far. The last
  state's `:worlds` are where `s` holds."
  [m s]
  (let [subs (subsentences (expand (:domain m) s))]
    (rest (reductions (fn [{:keys [labels]} sub]
                        (let [w (label m labels sub)]
                          {:sentence sub :worlds w :labels (assoc labels sub w)}))
                      {:labels {}} subs))))

(defn worlds-where
  "The worlds of `m` where `s` holds."
  [m s]
  (:worlds (last (check-seq m s))))

(defn holds?
  "Whether `s` holds in world `w` of `m`."
  [m w s]
  (contains? (worlds-where m s) w))

(defn valid?
  "Whether `s` holds in every world of `m`."
  [m s]
  (= (set (keys (:worlds m))) (worlds-where m s)))

;; ---------------------------------------------------------------------------
;; Frames

(defn reflexive? [m a] (every? (fn [[w ws]] (contains? ws w)) ((:access m) a)))

(defn transitive? [m a]
  (let [rel ((:access m) a)]
    (every? (fn [[_ ws]] (every? #(every? ws (rel %)) ws)) rel)))

(defn symmetric? [m a]
  (let [rel ((:access m) a)]
    (every? (fn [[w ws]] (every? #(contains? (rel %) w) ws)) rel)))
