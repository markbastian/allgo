(ns allgo.logic.first-order
  "First-order logic and inference in it: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapters 8 and 9, as
  steps.

    sentences      variables are lowercase symbols (`x`), constants
                   capitalized keywords (`:West`), and a compound term or
                   atomic sentence a vector `[:Owns :Nono :M1]`; then
                   `[:not s]`, `[:and ...]`, `[:or ...]`, `[:implies a b]`,
                   `[:iff a b]`, `[:forall [x y] s]`, `[:exists [x] s]`.
                   `parse` reads the book's notation, `show` writes it
    unification    UNIFY (figure 9.1), with the occurs check
    forward        FOL-FC-ASK (figure 9.3) on definite clauses: a round at
    chaining       a time, each step one rule applied to every way its
                   premises match the facts so far
    backward       FOL-BC-ASK (figure 9.6), the book's generators of
    chaining       substitutions run as one explicit stack of goal lists,
                   each step the first goal of the top entry resolved
                   against every rule -- depth first, rules in order,
                   which is what the generators do
    resolution     clausal form (implications out, negation in, variables
                   standardized, existentials Skolemized, universals
                   dropped, or distributed over and), then resolution
                   refutation with factoring, a resolvent a step, the set
                   of support the negated query; each clause remembers its
                   parents, so a refutation is a proof tree

  A substitution is a map from variables to terms."
  (:require [clojure.string :as str])
  #?(:clj (:import (clojure.lang PersistentQueue))))

(def ^:private empty-queue #?(:clj PersistentQueue/EMPTY :cljs #queue []))

;; ---------------------------------------------------------------------------
;; Terms

(defn variable? [x] (symbol? x))
(defn compound? [x] (and (vector? x) (keyword? (first x))
                         (not (#{:not :and :or :implies :iff :forall :exists} (first x)))))

(defn variables
  "The variables in `x`, in order of first appearance."
  [x]
  (distinct (filter variable? (tree-seq coll? seq x))))

(defn subst
  "SUBST: `x` with each variable `theta` binds replaced, through chains of
  bindings."
  [theta x]
  (cond
    (variable? x) (if (contains? theta x) (subst theta (theta x)) x)
    (vector? x) (if (#{:forall :exists} (first x))
                  [(first x) (second x) (subst (apply dissoc theta (second x)) (nth x 2))]
                  (mapv #(subst theta %) x))
    ;; A rule, `{:premises :conclusion}`.
    (map? x) (into {} (map (fn [[k v]] [k (subst theta v)]) x))
    :else x))

(defn- occurs? [v x theta]
  (cond
    (= v x) true
    (variable? x) (and (contains? theta x) (occurs? v (theta x) theta))
    (vector? x) (some #(occurs? v % theta) x)
    :else false))

(declare unify)

(defn- unify-var [v x theta]
  (cond
    (contains? theta v) (unify (theta v) x theta)
    (and (variable? x) (contains? theta x)) (unify v (theta x) theta)
    (occurs? v x theta) nil
    :else (assoc theta v x)))

(defn unify
  "UNIFY: the most general substitution extending `theta` (default empty)
  that makes `x` and `y` identical, or nil."
  ([x y] (unify x y {}))
  ([x y theta]
   (cond
     (nil? theta) nil
     (= x y) theta
     (variable? x) (unify-var x y theta)
     (variable? y) (unify-var y x theta)
     (and (vector? x) (vector? y) (= (count x) (count y)))
     (reduce (fn [th [a b]] (or (unify a b th) (reduced nil))) theta (map vector x y))
     :else nil)))

(def ^:private counter (atom 0))

(defn standardize-apart
  "`x` with every variable renamed to one used nowhere else."
  [x]
  (let [n (swap! counter inc)
        renaming (into {} (map (fn [v] [v (symbol (str (name v) "_" n))]) (variables x)))]
    (subst renaming x)))

(defn- canonical
  "`x` with its variables renamed in order of appearance: equal for two
  sentences exactly when one is a renaming of the other."
  [x]
  (subst (into {} (map-indexed (fn [i v] [v (symbol (str "v" i))]) (variables x))) x))

;; ---------------------------------------------------------------------------
;; Reading and writing

(defn- tokens [text]
  (re-seq #"<=>|=>|[~&|(),:\[\]]|[A-Za-z_][A-Za-z0-9_]*|-?[0-9]+(?:\.[0-9]+)?" text))

(defn- term-of [t]
  (cond
    (re-matches #"-?[0-9]+(\.[0-9]+)?" t) #?(:clj (read-string t) :cljs (js/parseFloat t))
    (re-matches #"[A-Z].*" t) (keyword t)
    :else (symbol t)))

(defn parse
  "A sentence in the book's notation: `Owns(Nono, M1)`, `~`, `&`, `|`,
  `=>`, `<=>`, and `forall x, y: S` or `exists x: S`, which reach as far
  right as they can. Capitalized names are constants, predicates and
  functions; lowercase ones are variables."
  [text]
  (let [toks (vec (tokens text))
        pos (volatile! 0)
        peek-tok #(get toks @pos)
        next-tok #(let [t (peek-tok)] (vswap! pos inc) t)
        expect (fn [t] (assert (= t (next-tok)) (str "expected " t)))]
    (letfn [(term []
              (let [t (next-tok)
                    head (term-of t)]
                (if (= "(" (peek-tok))
                  (do (next-tok)
                      (loop [args [(term)]]
                        (if (= "," (peek-tok))
                          (do (next-tok) (recur (conj args (term))))
                          (do (expect ")") (into [head] args)))))
                  head)))
            (unary []
              (let [t (peek-tok)]
                (cond
                  (= t "~") (do (next-tok) [:not (unary)])
                  (= t "(") (do (next-tok) (let [s (sentence)] (expect ")") s))
                  (#{"forall" "exists"} t)
                  (do (next-tok)
                      (let [vs (loop [vs [(symbol (next-tok))]]
                                 (if (= "," (peek-tok)) (do (next-tok) (recur (conj vs (symbol (next-tok))))) vs))]
                        (expect ":")
                        [(keyword t) vs (sentence)]))
                  :else (term))))
            (chain [sub op tok]
              (loop [args [(sub)]]
                (if (= tok (peek-tok))
                  (do (next-tok) (recur (conj args (sub))))
                  (if (= 1 (count args)) (first args) (into [op] args)))))
            (conj* [] (chain unary :and "&"))
            (disj* [] (chain conj* :or "|"))
            (implies [] (let [a (disj*)] (if (= "=>" (peek-tok)) (do (next-tok) [:implies a (implies)]) a)))
            (sentence [] (let [a (implies)] (if (= "<=>" (peek-tok)) (do (next-tok) [:iff a (sentence)]) a)))]
      (let [s (sentence)]
        (assert (= @pos (count toks)) (str "unparsed: " (str/join " " (subvec toks @pos))))
        s))))

(defn show
  "A sentence or term in the book's notation."
  [x]
  (cond
    (variable? x) (name x)
    (keyword? x) (name x)
    (compound? x) (str (name (first x)) "(" (str/join ", " (map show (rest x))) ")")
    (vector? x)
    (let [[op & args] x
          wrap #(if (or (not (vector? %)) (compound? %) (= :not (first %))) (show %) (str "(" (show %) ")"))]
      (case op
        :not (str "~" (wrap (first args)))
        :and (str/join " & " (map wrap args))
        :or (str/join " | " (map wrap args))
        :implies (str (wrap (first args)) " => " (wrap (second args)))
        :iff (str (wrap (first args)) " <=> " (wrap (second args)))
        (:forall :exists) (str (name op) " " (str/join ", " (map name (first args))) ": " (show (second args)))))
    :else (str x)))

(defn show-theta [theta vars]
  (str "{" (str/join ", " (for [v vars] (str (name v) "/" (show (subst theta v))))) "}"))

;; ---------------------------------------------------------------------------
;; Definite clauses

(defn definite
  "A definite clause from text, `p1 & ... & pn => q`, or a fact: `{:premises
  [atoms] :conclusion atom}`."
  [text]
  (let [s (parse text)]
    (if (and (vector? s) (= :implies (first s)))
      (let [[_ p q] s]
        {:premises (if (= :and (first p)) (vec (rest p)) [p]) :conclusion q})
      {:premises [] :conclusion s})))

(defn- matches
  "Every substitution extending `theta` under which each of `premises`
  equals some fact."
  [premises facts theta]
  (if-let [p (first premises)]
    (mapcat (fn [f] (when-let [th (unify (subst theta p) f theta)] (matches (rest premises) facts th)))
            facts)
    [theta]))

;; ---------------------------------------------------------------------------
;; Forward chaining

(defn- fc-step
  [{:keys [rules facts seen new todo query] :as st}]
  (if-let [rule (first todo)]
    (let [{:keys [premises conclusion]} (standardize-apart rule)
          inferred (for [th (matches premises facts {})
                         :let [q (subst th conclusion)]]
                     q)
          [new seen fresh] (reduce (fn [[new seen fresh] q]
                                     (let [c (canonical q)]
                                       (if (contains? seen c) [new seen fresh] [(conj new q) (conj seen c) (conj fresh q)])))
                                   [new seen []] inferred)]
      (assoc st :todo (rest todo) :new new :seen seen :rule rule :inferred fresh
             :answers (into (:answers st) (keep #(unify query % {}) fresh))))
    (when (seq new)
      ;; The round is over: what it found joins the facts.
      (assoc st :facts (into facts new) :new [] :todo (filter #(seq (:premises %)) rules)
             :round (inc (long (:round st))) :rule nil :inferred []))))

(defn fol-fc-ask-seq
  "Every state of FOL-FC-ASK on definite clauses `kb` (as `definite`
  makes them) and an atomic `query`: each step a rule applied to all the
  facts so far, `:inferred` what it added; a round ends when every rule
  has had its turn, and the next begins with what the round inferred. It
  ends when a round infers nothing. `:answers` gathers every substitution
  under which an inferred fact -- or a fact given -- matches the query;
  the book's stops at the first."
  [kb query]
  (let [facts (vec (map :conclusion (filter #(empty? (:premises %)) kb)))
        rules (filter #(seq (:premises %)) kb)]
    (take-while some?
                (iterate (fn [s] (when s (fc-step s)))
                         {:rules rules :facts facts :seen (set (map canonical facts)) :new [] :todo rules
                          :round 1 :query query :answers (vec (keep #(unify query % {}) facts))}))))

(defn fol-fc-ask
  "Every substitution for `query`'s variables forward chaining proves."
  [kb query]
  (let [vs (variables query)]
    (vec (distinct (map #(select-keys (into {} (map (fn [v] [v (subst % v)]) vs)) vs)
                        (:answers (last (fol-fc-ask-seq kb query))))))))

;; ---------------------------------------------------------------------------
;; Backward chaining

(defn- bc-step
  [{:keys [stack rules max-depth] :as st}]
  (when-let [[goals theta depth] (peek stack)]
    (let [stack (pop stack)
          st (assoc st :stack stack :steps (inc (long (:steps st))))]
      (if (empty? goals)
        (assoc st :answer theta :goal nil)
        (let [goal (subst theta (first goals))
              alternatives (when (< (long depth) (long max-depth))
                             (for [rule rules
                                   :let [{:keys [premises conclusion]} (standardize-apart rule)
                                         th (unify conclusion goal theta)]
                                   :when th]
                               [(into (vec premises) (rest goals)) th (inc (long depth))]))]
          ;; Reversed, so that the first rule is tried first.
          (assoc st :stack (into stack (reverse alternatives)) :answer nil :goal goal
                 :cut-off? (>= (long depth) (long max-depth))))))))

(defn fol-bc-ask-seq
  "Every state of FOL-BC-ASK on definite clauses `kb` and an atomic
  `query`: the stack of goal lists waiting, each with its substitution and
  depth; each step takes the top one and replaces it with one entry for
  each rule whose conclusion unifies with its first goal, that rule's
  premises put in front of the goals left. A state with an empty goal
  list popped has an `:answer`. Proofs deeper than `:max-depth` (default
  50) are cut off, since a rule like `Human(h) & Mother(m, h) =>
  Human(m)` in the wrong order sends depth-first search down forever."
  ([kb query] (fol-bc-ask-seq kb query {}))
  ([kb query {:keys [max-depth] :or {max-depth 50}}]
   (take-while some?
               (iterate (fn [s] (when s (bc-step s)))
                        {:rules (vec kb) :stack [[[query] {} 0]] :steps 0 :max-depth max-depth}))))

(defn fol-bc-ask
  "The substitutions for `query`'s variables backward chaining proves, in
  the order it proves them, as a lazy seq."
  ([kb query] (fol-bc-ask kb query {}))
  ([kb query opts]
   (let [vs (variables query)]
     (distinct (keep (fn [{:keys [answer]}]
                       (when answer (into {} (map (fn [v] [v (subst answer v)]) vs))))
                     (fol-bc-ask-seq kb query opts))))))

;; ---------------------------------------------------------------------------
;; Clausal form

(defn- eliminate [s]
  (cond
    (not (vector? s)) s
    (compound? s) s
    :else
    (let [[op & args] s]
      (case op
        :iff (let [[a b] (map eliminate args)] [:and [:or [:not a] b] [:or [:not b] a]])
        :implies (let [[a b] (map eliminate args)] [:or [:not a] b])
        (:forall :exists) [op (first args) (eliminate (second args))]
        (into [op] (map eliminate args))))))

(defn- push-not [s]
  (cond
    (or (not (vector? s)) (compound? s)) s
    :else
    (let [[op & args] s]
      (case op
        :not (let [a (first args)]
               (if (or (not (vector? a)) (compound? a))
                 s
                 (let [[op2 & args2] a]
                   (case op2
                     :not (push-not (first args2))
                     :and (into [:or] (map #(push-not [:not %]) args2))
                     :or (into [:and] (map #(push-not [:not %]) args2))
                     :forall [:exists (first args2) (push-not [:not (second args2)])]
                     :exists [:forall (first args2) (push-not [:not (second args2)])]))))
        (:forall :exists) [op (first args) (push-not (second args))]
        (into [op] (map push-not args))))))

(defn- skolemize
  "Each quantified variable renamed apart, each existential replaced by a
  new function of the universals around it, and the universals dropped."
  [s universals]
  (cond
    (or (not (vector? s)) (compound? s)) s
    :else
    (let [[op & args] s]
      (case op
        :forall (let [renaming (into {} (map (fn [v] [v (symbol (str (name v) "_" (swap! counter inc)))]) (first args)))]
                  (skolemize (subst renaming (second args)) (into universals (vals renaming))))
        :exists (let [skolems (into {} (map (fn [v]
                                              (let [f (keyword (str "Sk" (swap! counter inc)))]
                                                [v (if (seq universals) (into [f] universals) f)]))
                                            (first args)))]
                  (skolemize (subst skolems (second args)) universals))
        :not s
        (into [op] (map #(skolemize % universals) args))))))

(defn- literal? [s] (or (compound? s) (keyword? s) (and (vector? s) (= :not (first s)))))
(defn- negate [lit] (if (and (vector? lit) (= :not (first lit))) (second lit) [:not lit]))

(defn- clauses-of [s]
  (cond
    (literal? s) [#{s}]
    (= :and (first s)) (mapcat clauses-of (rest s))
    (= :or (first s)) (reduce (fn [acc cs] (for [a acc c cs] (into a c))) [#{}] (map clauses-of (rest s)))))

(defn clausal-form
  "The clauses of `sentences`, each a set of literals, variables implicitly
  universal and existentials Skolemized."
  [& sentences]
  (vec (distinct (mapcat #(clauses-of (skolemize (push-not (eliminate %)) [])) sentences))))

;; ---------------------------------------------------------------------------
;; Resolution

(defn- rename-clause [c]
  (let [n (swap! counter inc)]
    (set (subst (into {} (map (fn [v] [v (symbol (str (name v) "_" n))]) (variables (vec c)))) (vec c)))))

(defn- factors
  "The clause with pairs of its literals that unify merged: factoring."
  [clause]
  (let [lits (vec clause)]
    (distinct
     (for [i (range (count lits)) j (range (inc i) (count lits))
           :let [th (unify (lits i) (lits j))]
           :when th]
       (set (subst th lits))))))

(defn- resolvents
  "Every clause from resolving `a` and `b` on one complementary pair."
  [a b]
  (let [b (rename-clause b)]
    (for [la a lb b
          :let [th (unify la (negate lb))]
          :when th
          :let [c (set (subst th (vec (concat (disj a la) (disj b lb)))))]]
      c)))

(defn- tautology? [c] (some #(contains? c (negate %)) c))

(defn- resolution-step
  [{:keys [usable queue] :as st}]
  (when (nil? (:proved? st))
    (if-let [given (peek queue)]
      ;; The given clause: resolve it against everything kept so far.
      (let [queue (pop queue)
            made (for [other (conj usable given)
                       r (resolvents (:clause given) (:clause other))
                       c (cons r (factors r))
                       :when (not (tautology? c))]
                   {:clause c :parents [given other]})
            seen (:seen st)
            [fresh seen] (reduce (fn [[fresh seen] {:keys [clause] :as m}]
                                   (let [k (canonical (vec (sort-by str clause)))]
                                     (if (contains? seen k) [fresh seen] [(conj fresh m) (conj seen k)])))
                                 [[] seen] made)
            empty-clause (first (filter #(empty? (:clause %)) fresh))]
        (cond-> (assoc st :usable (conj usable given) :queue (into queue (sort-by #(count (:clause %)) fresh))
                       :seen seen :given given :made fresh :steps (inc (long (:steps st))))
          empty-clause (assoc :proved? true :proof empty-clause)))
      (assoc st :proved? false))))

(defn resolution-seq
  "Every state of resolution refutation of `query` from `kb`, sentences
  both: the clauses of the KB kept as usable, the clauses of the negated
  query queued as the set of support; each step takes the next queued
  clause -- the given clause -- and resolves it with everything usable,
  queueing what is new, shortest first. It ends with `:proved?` true and
  the empty clause as `:proof`, or false when the queue runs dry; and may
  not end at all, first-order logic being semidecidable -- take as many
  steps as there is patience for."
  [kb query]
  (let [kb-clauses (apply clausal-form kb)
        goal (clausal-form [:not query])]
    (take-while some?
                (iterate (fn [s] (when s (resolution-step s)))
                         {:usable (mapv (fn [c] {:clause c :parents nil}) kb-clauses)
                          :queue (into empty-queue (map (fn [c] {:clause c :parents nil :goal? true}) goal))
                          :seen (set (map #(canonical (vec (sort-by str %))) (concat kb-clauses goal)))
                          :steps 0}))))

(defn resolution
  "Whether resolution refutes ~`query` from `kb` within `max-steps`
  (default 500): true, false, or nil for undecided."
  ([kb query] (resolution kb query 500))
  ([kb query max-steps]
   (:proved? (last (take max-steps (resolution-seq kb query))))))

(defn proof-lines
  "A refutation as lines, each a clause and the two it came from, leaves
  first."
  [node]
  (let [lines (volatile! [])
        ids (volatile! {})
        id-of (fn walk [n]
                (or (@ids n)
                    (let [ps (mapv walk (:parents n))
                          i (count @lines)]
                      (vswap! lines conj {:clause (:clause n) :from ps :goal? (:goal? n)})
                      (vswap! ids assoc n i)
                      i)))]
    (id-of node)
    @lines))
