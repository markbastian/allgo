(ns allgo.logic.tms
  "Truth maintenance: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 10 (chapter 12 of the 3rd
  edition, whose text this follows). The chapter describes both systems
  and leaves out their algorithms; these are the usual ones.

    JTMS   justification-based: every sentence keeps the justifications
           it was inferred from -- Q from {P, P => Q} -- and is in while
           some justification has all its antecedents in. Retracting P
           takes out exactly what has no justification left without it:
           first everything that depended on P is marked out, then each
           is brought back in if another justification still holds,
           which is what keeps two sentences that justify each other from
           holding each other up. Nothing is forgotten, so re-asserting P
           brings its consequences back without inferring them again
    ATMS   assumption-based: every sentence is labeled with the minimal
           sets of assumptions under which it holds, all at once, so
           switching between hypothetical worlds needs no retraction. A
           justification for a contradiction makes its assumption sets
           nogoods, struck from every label. A label is a sentence's
           explanations: the car won't start because the battery is
           dead, or because the tank is empty

  Both run a step at a time: a node, or a justification, off a queue.
  Sentences are any values -- keywords here, `:P` and `:P=>Q` -- and a
  justification `{:antecedents #{nodes} :consequent node}`; `rules`
  makes them from definite clauses, each rule a node of its own, and
  `horn` with the rules left as background."
  (:require [allgo.logic.propositional :as pl]
            [allgo.search :as search]
            [clojure.set :as set]))

(defn rules
  "Justifications from definite clauses in text, `A & B => C`, each
  justified by its premises and the rule itself, named by its text."
  [& texts]
  (vec (for [t texts
             :let [{:keys [premise conclusion]} (pl/definite t)]
             :when (seq premise)]
         {:antecedents (conj premise t) :consequent conclusion})))

(defn horn
  "Justifications from definite clauses in text, `A & B => C`, justified
  by their premises alone: background knowledge, as an ATMS takes it,
  rather than sentences that may be retracted."
  [& texts]
  (vec (for [t texts
             :let [{:keys [premise conclusion]} (pl/definite t)]
             :when (seq premise)]
         {:antecedents premise :consequent conclusion})))

;; ---------------------------------------------------------------------------
;; JTMS

(defn jtms
  "A JTMS with `justifications` and nothing yet enabled."
  [justifications]
  {:justifications (vec justifications) :enabled #{} :in #{}})

(defn- supported? [{:keys [justifications in]} n]
  (some #(and (= n (:consequent %)) (every? in (:antecedents %))) justifications))

(defn- consequents-of [{:keys [justifications]} n]
  (distinct (map :consequent (filter #(contains? (:antecedents %) n) justifications))))

(defn- propagate-in
  "One node off the queue: each consequent it newly supports comes in."
  [{:keys [queue] :as st}]
  (let [n (first queue)
        ins (remove (:in st) (filter #(supported? st %) (consequents-of st n)))]
    (assoc st :node n :queue (into (subvec queue 1) ins) :in (into (:in st) ins) :changed ins)))

(defn- jtms-step
  [{:keys [phase queue] :as st}]
  (cond
    (seq queue)
    (case phase
      :in (propagate-in st)
      ;; Out: every node that depends on this one is marked.
      :out (let [n (first queue)
                 outs (filter (:in st) (consequents-of st n))]
             (assoc st :node n :queue (into (subvec queue 1) outs)
                    :in (reduce disj (:in st) outs) :marked (into (:marked st) outs) :changed outs))
      ;; Back in: a marked node with a justification that holds.
      :support (let [n (first queue)]
                 (if (and (not ((:in st) n)) (supported? st n))
                   (assoc st :node n :phase :in :queue [n] :in (conj (:in st) n) :changed [n]
                          :pending (subvec queue 1))
                   (assoc st :node n :queue (subvec queue 1) :changed []))))
    (= phase :out) (assoc st :phase :support :queue (vec (:marked st)) :node nil :changed [])
    (seq (:pending st)) (assoc st :phase :support :queue (:pending st) :pending [] :node nil :changed [])
    :else nil))

(defn enable-seq
  "Asserting `node` -- a premise or an assumption -- a step a node: what
  it and the rest now support comes in. Each state's `:in` the nodes in
  and `:changed` those this step brought in."
  [tms node]
  (search/step-seq jtms-step (-> tms
                                 (update :enabled conj node) (update :in conj node)
                                 (assoc :phase :in :queue [node] :marked #{} :pending [] :node node :changed [node]))))

(defn retract-seq
  "Retracting `node`, a step a node: first it and everything that
  depended on it are marked out (`:phase :out`), then each marked node
  with a justification still holding comes back in (`:support`, then
  `:in` for what that brings with it)."
  [tms node]
  (let [tms (update tms :enabled disj node)]
    (if (supported? tms node)
      [tms]
      (search/step-seq jtms-step (-> tms
                                     (update :in disj node)
                                     (assoc :phase :out :queue [node] :marked #{} :pending [] :node node :changed [node]))))))

(defn enable [tms & nodes] (reduce #(-> (last (enable-seq %1 %2)) (dissoc :queue :phase :marked :pending)) tms nodes))
(defn retract [tms & nodes] (reduce #(-> (last (retract-seq %1 %2)) (dissoc :queue :phase :marked :pending)) tms nodes))
(defn in? [tms node] (contains? (:in tms) node))

;; ---------------------------------------------------------------------------
;; ATMS

(defn- minimize
  "The sets in `envs` that contain no other."
  [envs]
  (let [envs (set envs)]
    (set (remove (fn [e] (some #(and (not= % e) (set/subset? % e)) envs)) envs))))

(defn- consistent [envs nogoods] (remove (fn [e] (some #(set/subset? % e) nogoods)) envs))

(defn atms
  "An ATMS: `assumptions`, each labeled with itself; `premises`, labeled
  with the empty set, true under any assumptions; `justifications`, all
  queued to propagate. A justification whose consequent is
  `:contradiction` makes nogoods."
  [assumptions premises justifications]
  {:justifications (vec justifications)
   :labels (merge (zipmap premises (repeat #{#{}})) (into {} (map (fn [a] [a #{#{a}}]) assumptions)))
   :nogoods #{} :queue (vec justifications)})

(defn- atms-step
  [{:keys [queue labels nogoods justifications] :as st}]
  (when-let [{:keys [antecedents consequent] :as j} (first queue)]
    (let [envs (reduce (fn [acc a] (for [e acc l (get labels a)] (set/union e l))) [#{}] antecedents)
          envs (minimize (consistent envs nogoods))
          queue (subvec queue 1)]
      (if (= :contradiction consequent)
        (let [nogoods (minimize (into nogoods envs))
              labels (into {} (map (fn [[n ls]] [n (set (consistent ls nogoods))])) labels)]
          (assoc st :justification j :queue queue :nogoods nogoods :labels labels :new envs))
        (let [old (get labels consequent #{})
              label (minimize (into old envs))]
          (assoc st :justification j :labels (assoc labels consequent label) :new (set/difference label old)
                 :queue (if (= label old)
                          queue
                          (into queue (filter #(contains? (:antecedents %) consequent)) justifications))))))))

(defn atms-seq
  "Label propagation in `atms`, a justification a step: each state's
  `:justification` the one just used, `:new` the environments it added
  to its consequent's label (or to the nogoods), `:labels` and
  `:nogoods` so far. A consequent whose label grows queues every
  justification it is an antecedent of."
  [atms]
  (search/step-seq atms-step atms))

(defn labels
  "Each node's label once propagation is done."
  [atms]
  (:labels (last (atms-seq atms))))

(defn holds-under?
  "Whether `node` holds under the assumptions `env`: some environment of
  its label is within it, and no nogood is."
  [{:keys [labels nogoods]} node env]
  (and (not-any? #(set/subset? % env) nogoods)
       (boolean (some #(set/subset? % env) (get labels node)))))
