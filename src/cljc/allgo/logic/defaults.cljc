(ns allgo.logic.defaults
  "Reasoning with default information: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 10 (chapter 12 of
  the 3rd edition, whose text this follows). Both logics here are
  nonmonotonic: something believed by default can stop being believed
  when more is learned.

    circumscription  some predicates -- Abnormal1, say, in Bird(x) &
                     ~Abnormal1(x) => Flies(x) -- assumed as false as they
                     can be: what holds in every preferred model, a model
                     preferred when no other makes a proper subset of
                     those atoms true. Prioritized, the predicates
                     minimized one group after another, so religion can
                     outrank politics in the Nixon diamond
    default logic    rules P : J1, ..., Jn / C -- given P, and with each J
                     consistent with what is believed, conclude C. An
                     extension is a set of beliefs closed under the facts
                     and the rules applied, with no rule's justification
                     contradicted. Extensions are found by trying the
                     rules in every order that stays consistent, a
                     process a step; the Nixon diamond has two

  Sentences are first-order, written as `allgo.logic.first-order` reads
  them, and function-free: they are made propositional by replacing every
  variable with each constant they mention, each ground atom becoming a
  symbol, and then decided by `allgo.logic.propositional`'s DPLL."
  (:require [allgo.logic.first-order :as fo]
            [allgo.logic.propositional :as pl]
            [allgo.search :as search]
            [clojure.string :as str]))

;; ---------------------------------------------------------------------------
;; Grounding

(defn constants
  "The constants among the arguments of the atoms in `sentences`."
  [sentences]
  (set (for [s sentences
             x (tree-seq vector? rest s)
             :when (fo/compound? x)
             a (rest x)
             :when (keyword? a)]
         a)))

(defn- ground [objs s]
  (cond
    (fo/compound? s) s
    (keyword? s) s
    :else
    (let [[op & args] s]
      (case op
        (:forall :exists)
        (let [[vs body] args
              instances (for [combo (reduce (fn [acc _] (for [a acc o objs] (conj a o))) [[]] vs)]
                          (ground objs (fo/subst (zipmap vs combo) body)))]
          (into [(if (= op :forall) :and :or)] instances))
        (into [op] (map #(ground objs %) args))))))

(defn symbol-of
  "The propositional symbol for a ground atom: `[:Flies :Tweety]` is
  `:Flies(Tweety)`."
  [atom]
  (if (keyword? atom) atom (keyword (fo/show atom))))

(defn- propositional [s]
  (cond
    (fo/compound? s) (symbol-of s)
    (keyword? s) s
    :else (into [(first s)] (map propositional (rest s)))))

(defn propositionalize
  "`sentences`, read from text if they are strings, with free variables
  taken as universal, made propositional over `objs` (default: their
  constants)."
  ([sentences] (propositionalize sentences nil))
  ([sentences objs]
   (let [ss (map #(if (string? %) (fo/parse %) %) sentences)
         objs (or objs (constants ss))]
     (vec (for [s ss
                :let [vs (fo/variables s)
                      s (if (seq vs) [:forall (vec vs) s] s)]]
            (propositional (ground objs s)))))))

;; ---------------------------------------------------------------------------
;; Models

(defn- completions [model symbols]
  (reduce (fn [ms p] (if (contains? model p) ms (for [m ms v [true false]] (assoc m p v))))
          [model] symbols))

(defn models-seq
  "Every model of propositional `sentences` over their symbols, a model a
  step: DPLL finds one, and a clause ruling out that assignment is added
  before the next search. Each state `{:model :count}`."
  [sentences]
  (let [symbols (sort (apply pl/symbols sentences))]
    (letfn [(go [clauses pending n]
                (lazy-seq
                 (if (seq pending)
                   (cons {:model (first pending) :count (inc (long n))} (go clauses (rest pending) (inc (long n))))
                   (when-let [m (pl/dpll clauses)]
                   ;; DPLL may leave symbols unset: every way of setting
                   ;; them is a model too.
                     (go (conj clauses (set (map (fn [[p v]] (if v [:not p] p)) m)))
                         (completions m symbols) n)))))]
      (go (apply pl/cnf sentences) nil 0))))

(defn models [sentences] (mapv :model (models-seq sentences)))

;; ---------------------------------------------------------------------------
;; Circumscription

(defn- true-of [model preds]
  (set (filter (fn [p] (and (model p) (some #(str/starts-with? (name p) (str (name %) "(")) preds)))
               (keys model))))

(defn- minimal [models preds]
  (let [ab (map #(true-of % preds) models)]
    (vec (for [[m a] (map vector models ab)
               :when (not-any? #(and (not= % a) (every? a %)) ab)]
           m))))

(defn preferred-models
  "The preferred models of `kb` when the predicates `circumscribed` are
  made as false as they can be: those no other model betters by making
  a proper subset of their atoms true. `circumscribed` is a set of
  predicate names, or for prioritized circumscription a vector of sets,
  minimized first to last."
  [kb circumscribed]
  (let [ms (models (propositionalize kb))]
    (reduce minimal ms (if (set? circumscribed) [circumscribed] circumscribed))))

(defn circumscription-entails?
  "Whether `query` holds in every preferred model of `kb`."
  [kb circumscribed query]
  (let [q (first (propositionalize [query] (constants (map #(if (string? %) (fo/parse %) %) (cons query kb)))))]
    (every? #(pl/pl-true? q %) (preferred-models kb circumscribed))))

;; ---------------------------------------------------------------------------
;; Default logic

(defn- split-top
  "`text` split at each `sep` outside parentheses."
  [text sep]
  (loop [cs (seq text) depth 0 cur [] out []]
    (if-let [c (first cs)]
      (cond
        (and (= c sep) (zero? depth)) (recur (rest cs) depth [] (conj out (apply str cur)))
        :else (recur (rest cs) (case c \( (inc depth) \) (dec depth) depth) (conj cur c) out))
      (conj out (apply str cur)))))

(defn default-rule
  "A default rule from text, `P : J1, ..., Jn / C`: `{:prerequisite
  :justifications :conclusion}`, its sentences parsed."
  [text]
  (let [[p rest-text] (str/split text #":" 2)
        [js c] (let [i (str/last-index-of rest-text "/")] [(subs rest-text 0 i) (subs rest-text (inc i))])]
    {:prerequisite (fo/parse (str/trim p))
     :justifications (mapv (comp fo/parse str/trim) (split-top js \,))
     :conclusion (fo/parse (str/trim c))}))

(defn- ground-defaults [defaults objs]
  (distinct
   (for [d defaults
         :let [vs (fo/variables [(:prerequisite d) (:justifications d) (:conclusion d)])]
         combo (reduce (fn [acc _] (for [a acc o objs] (conj a o))) [[]] vs)
         :let [th (zipmap vs combo)
               g #(propositional (ground objs (fo/subst th %)))]]
     {:prerequisite (g (:prerequisite d)) :justifications (mapv g (:justifications d))
      :conclusion (g (:conclusion d)) :text (str (fo/show (fo/subst th (:prerequisite d))) " : "
                                                 (str/join ", " (map #(fo/show (fo/subst th %)) (:justifications d)))
                                                 " / " (fo/show (fo/subst th (:conclusion d))))})))

(defn- consistent? [beliefs s] (pl/satisfiable? (into [:and] (conj beliefs s))))
(defn- entails? [beliefs s] (pl/entails? (into [:and] beliefs) s))

(defn- process-step
  [{:keys [facts defaults stack extensions] :as st}]
  (when-let [applied (peek stack)]
    (let [stack (pop stack)
          beliefs (into (vec facts) (map :conclusion applied))
          failed? (or (not (pl/satisfiable? (into [:and] beliefs)))
                      (some (fn [j] (not (consistent? beliefs j))) (mapcat :justifications applied)))
          applicable (when-not failed?
                       (filter (fn [d] (and (not (some #{d} applied))
                                            (entails? beliefs (:prerequisite d))
                                            (every? #(consistent? beliefs %) (:justifications d))))
                               defaults))
          closed? (and (not failed?) (empty? applicable))
          conclusions (set (map :conclusion applied))
          new? (and closed? (not-any? #(= conclusions (:conclusions %)) extensions))]
      (assoc st :process applied :failed? (boolean failed?) :closed? closed?
             :stack (into stack (map #(conj applied %) (reverse applicable)))
             :extensions (if new? (conj extensions {:conclusions conclusions :applied applied :beliefs beliefs})
                             extensions)))))

(defn extensions-seq
  "The extensions of the default theory of `facts` and `defaults`
  (`default-rule`s, or text), found by trying the rules in every order, a
  process a step, depth first: each state's `:process` the rules applied
  in order, `:failed?` if a later conclusion contradicts an earlier
  rule's justification, `:closed?` if no rule is left to apply -- a
  successful closed process makes an extension, gathered in
  `:extensions` with its `:conclusions` and `:beliefs`."
  [facts defaults]
  (let [defaults (map #(if (string? %) (default-rule %) %) defaults)
        facts (map #(if (string? %) (fo/parse %) %) facts)
        objs (constants (concat facts (mapcat (juxt :prerequisite :conclusion) defaults)))]
    (search/step-seq process-step (process-step {:facts (propositionalize facts objs)
                                                 :defaults (vec (ground-defaults defaults objs))
                                                 :stack [[]] :extensions []}))))

(defn extensions
  "Each extension: `{:conclusions :applied :beliefs}`."
  [facts defaults]
  (:extensions (last (extensions-seq facts defaults))))

(defn believes?
  "Whether `extension`'s beliefs entail `query`."
  [extension query]
  (entails? (:beliefs extension) (first (propositionalize [query] #{}))))
