(ns allgo.logic.propositional
  "Propositional logic: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 7, as steps.

    sentences        keywords for symbols, and `[:not s]`, `[:and s ...]`,
                     `[:or s ...]`, `[:implies a b]`, `[:iff a b]`; `parse`
                     reads them from text, `~ & | => <=>` by precedence
    model checking   TT-ENTAILS? (figure 7.10): every model of the symbols,
                     one a step
    resolution       PL-RESOLUTION (figure 7.13) on the clauses of
                     KB & ~query, one pair resolved a step
    forward chaining PL-FC-ENTAILS? (figure 7.15) on definite clauses, one
                     symbol off the agenda a step; and backward chaining,
                     which the chapter describes without a figure, its
                     recursion on `allgo.search/call-seq`
    DPLL             DPLL-SATISFIABLE? (figure 7.17), with pure symbols and
                     unit clauses, its recursion on `call-seq`, a step a
                     call or return
    WalkSAT          figure 7.18: a flip a step

  A clause is a set of literals, a literal a symbol or `[:not symbol]`; a
  model a map from symbols to booleans. Every algorithm is a seq of
  states, run by `iterate` or `call-seq`, and has a plain form that runs
  it to its answer."
  (:require [allgo.random :as random]
            [allgo.search :as search]
            [clojure.string :as str])
  #?(:clj (:import (clojure.lang PersistentQueue))))

(def ^:private empty-queue #?(:clj PersistentQueue/EMPTY :cljs #queue []))

;; ---------------------------------------------------------------------------
;; Sentences

(defn literal? [s] (or (keyword? s) (and (vector? s) (= :not (first s)) (keyword? (second s)))))

(defn negate
  "The negation of a literal, or of any sentence."
  [s]
  (if (and (vector? s) (= :not (first s))) (second s) [:not s]))

(defn symbol-of [lit] (if (keyword? lit) lit (second lit)))

(defn symbols
  "The proposition symbols of `sentences`, in order of first appearance."
  [& sentences]
  (distinct (filter keyword? (remove #{:not :and :or :implies :iff} (mapcat #(tree-seq vector? rest %) sentences)))))

(defn pl-true?
  "PL-TRUE?: whether `s` holds in `model`; nil if `model` leaves it
  undecided."
  [s model]
  (if (keyword? s)
    (get model s)
    (let [[op & args] s]
      (case op
        :not (let [v (pl-true? (first args) model)] (when (some? v) (not v)))
        :and (let [vs (map #(pl-true? % model) args)]
               (cond (some false? vs) false (every? true? vs) true :else nil))
        :or (let [vs (map #(pl-true? % model) args)]
              (cond (some true? vs) true (every? false? vs) false :else nil))
        :implies (pl-true? [:or [:not (first args)] (second args)] model)
        :iff (let [[a b] (map #(pl-true? % model) args)]
               (when (and (some? a) (some? b)) (= a b)))))))

;; ---------------------------------------------------------------------------
;; Reading sentences

(defn- tokens [text]
  (re-seq #"<=>|=>|[~&|()]|[A-Za-z_][A-Za-z0-9_,]*" text))

(defn parse
  "A sentence from text: symbols are words (letters, digits, `_` and `,`,
  so `P1,2` is one), `~` is not, then `&`, `|`, `=>` and `<=>` by
  precedence, `=>` and `<=>` grouping to the right; parentheses as usual."
  [text]
  (let [toks (vec (tokens text))
        pos (volatile! 0)
        peek-tok #(get toks @pos)
        next-tok #(let [t (peek-tok)] (vswap! pos inc) t)]
    (letfn [(atom* []
              (let [t (next-tok)]
                (cond
                  (= t "~") [:not (atom*)]
                  (= t "(") (let [e (iff)] (next-tok) e)
                  :else (keyword t))))
            (chain [sub op tok]
              (loop [args [(sub)]]
                (if (= tok (peek-tok))
                  (do (next-tok) (recur (conj args (sub))))
                  (if (= 1 (count args)) (first args) (into [op] args)))))
            (conj* [] (chain atom* :and "&"))
            (disj* [] (chain conj* :or "|"))
            (implies [] (let [a (disj*)] (if (= "=>" (peek-tok)) (do (next-tok) [:implies a (implies)]) a)))
            (iff [] (let [a (implies)] (if (= "<=>" (peek-tok)) (do (next-tok) [:iff a (iff)]) a)))]
      (let [s (iff)]
        (when-not (= @pos (count toks))
          (throw (ex-info (str "unparsed: " (str/join " " (subvec toks @pos))) {:text text})))
        s))))

(defn show
  "A sentence written out as `parse` reads it."
  [s]
  (if (keyword? s)
    (name s)
    (let [[op & args] s
          wrap #(if (keyword? %) (show %) (str "(" (show %) ")"))]
      (case op
        :not (str "~" (wrap (first args)))
        :and (str/join " & " (map wrap args))
        :or (str/join " | " (map wrap args))
        :implies (str (wrap (first args)) " => " (wrap (second args)))
        :iff (str (wrap (first args)) " <=> " (wrap (second args)))))))

;; ---------------------------------------------------------------------------
;; Conjunctive normal form

(defn- eliminate
  "Biconditionals and implications rewritten with and, or and not."
  [s]
  (if (keyword? s)
    s
    (let [[op & args] s
          args (map eliminate args)]
      (case op
        :iff (let [[a b] args] [:and [:or [:not a] b] [:or [:not b] a]])
        :implies (let [[a b] args] [:or [:not a] b])
        (into [op] args)))))

(defn- push-not
  "Negations moved inward to the symbols, by De Morgan's laws."
  [s]
  (if (keyword? s)
    s
    (let [[op & args] s]
      (case op
        :not (let [a (first args)]
               (if (keyword? a)
                 s
                 (let [[op2 & args2] a]
                   (case op2
                     :not (push-not (first args2))
                     :and (into [:or] (map #(push-not [:not %]) args2))
                     :or (into [:and] (map #(push-not [:not %]) args2))))))
        (into [op] (map push-not args))))))

(defn- clauses-of
  "The clauses of a sentence in negation normal form: a conjunction's are
  its parts', a disjunction's every way of picking one clause from each
  part, joined -- distributing or over and. Tautologies are dropped."
  [s]
  (let [tautology? (fn [c] (some #(contains? c (negate %)) c))]
    (cond
      (literal? s) [#{s}]
      (= :and (first s)) (vec (distinct (mapcat clauses-of (rest s))))
      (= :or (first s)) (->> (rest s)
                             (map clauses-of)
                             (reduce (fn [acc cs] (for [a acc c cs] (into a c))) [#{}])
                             (remove tautology?)
                             distinct
                             vec))))

(defn cnf
  "The clauses of `sentences`, each a set of literals: their conjunction
  in conjunctive normal form."
  [& sentences]
  (vec (distinct (mapcat (comp clauses-of push-not eliminate) sentences))))

(defn clause-true?
  "Whether `model` satisfies `clause`: true, false, or nil if undecided."
  [clause model]
  (let [vs (map (fn [lit] (let [v (get model (symbol-of lit))]
                            (when (some? v) (if (keyword? lit) v (not v)))))
                clause)]
    (cond (some true? vs) true (every? false? vs) false :else nil)))

;; ---------------------------------------------------------------------------
;; Model checking

(defn tt-entails-seq
  "Every state of TT-ENTAILS?: one model of the symbols of `kb` and
  `query` a step, until one makes the KB true and the query false. The
  last state's `:entails?` is the answer; `:model` the model checked, and
  `:checked` how many."
  [kb query]
  (let [syms (vec (symbols kb query))
        n (count syms)
        total (long (reduce * 1 (repeat n 2)))
        model-of (fn [k] (into {} (map-indexed (fn [i s] [s (odd? (quot k (long (reduce * 1 (repeat i 2)))))]) syms)))]
    (take-while some?
                (iterate (fn [{:keys [k entails?] :as st}]
                           (when (and (< (long k) total) (not (false? entails?)))
                             (let [m (model-of k)
                                   counter? (and (pl-true? kb m) (not (pl-true? query m)))]
                               (assoc st :k (inc (long k)) :model m :checked (inc (long k))
                                      :entails? (if counter? false (if (= (inc (long k)) total) true nil))))))
                         {:k 0 :checked 0 :entails? (when (zero? total) (boolean (pl-true? query {})))}))))

(defn tt-entails?
  "Whether `kb` entails `query`, by checking every model."
  [kb query]
  (true? (:entails? (last (tt-entails-seq kb query)))))

;; ---------------------------------------------------------------------------
;; Resolution

(defn resolve-clauses
  "PL-RESOLVE: every clause from resolving `a` and `b` on a complementary
  pair of literals, tautologies dropped."
  [a b]
  (for [lit a
        :when (contains? b (negate lit))
        :let [c (into (disj a lit) (disj b (negate lit)))]
        :when (not (some #(contains? c (negate %)) c))]
    c))

(defn- resolution-step
  [{:keys [clauses pairs new] :as st}]
  (when-not (some? (:entails? st))
    (if-let [[a b] (first pairs)]
      (let [rs (resolve-clauses a b)
            st (assoc st :pairs (rest pairs) :resolved [a b] :steps (inc (long (:steps st))))]
        (if (some empty? rs)
          (assoc st :entails? true :new (into new rs))
          (assoc st :new (into new rs))))
      ;; The round is over: nothing new means no proof.
      (if (every? (set clauses) new)
        (assoc st :entails? false :pairs nil)
        (let [old (set clauses)
              fresh (vec (remove old new))
              all (into (vec clauses) fresh)]
          ;; Only pairs with a new clause in them: the old pairs were
          ;; all resolved last round.
          (assoc st :clauses all :new #{} :round (inc (long (:round st)))
                 :pairs (for [i (range (count all)) j (range (max (inc i) (count clauses)) (count all))
                              :let [a (all i) b (all j)]
                              :when (some #(contains? b (negate %)) a)]
                          [a b])))))))

(defn pl-resolution-seq
  "Every state of PL-RESOLUTION: the clauses of `kb` and the negated
  `query`, resolved a pair a step, a round at a time, adding what a round
  makes, until the empty clause appears -- the query is entailed -- or a
  round makes nothing new -- it is not. `:entails?` is set at the end."
  [kb query]
  (let [clauses (cnf kb [:not query])]
    (take-while some?
                (iterate (fn [s] (when s (resolution-step s)))
                         {:clauses clauses :new #{} :round 1 :steps 0
                          :pairs (for [i (range (count clauses)) j (range (inc i) (count clauses))
                                       :let [a (clauses i) b (clauses j)]
                                       :when (some #(contains? b (negate %)) a)]
                                   [a b])}))))

(defn pl-resolution
  "Whether `kb` entails `query`, by resolution refutation."
  [kb query]
  (true? (:entails? (last (pl-resolution-seq kb query)))))

;; ---------------------------------------------------------------------------
;; Forward and backward chaining

(defn definite
  "A definite clause from text: `A & B => C`, or a fact `C`."
  [text]
  (let [s (parse text)]
    (if (and (vector? s) (= :implies (first s)))
      (let [[_ p c] s]
        {:premise (if (keyword? p) #{p} (set (rest p))) :conclusion c})
      {:premise #{} :conclusion s})))

(defn pl-fc-entails-seq
  "Every state of PL-FC-ENTAILS? on definite clauses `kb` (as `definite`
  makes them) and the symbol `q`: one symbol off the agenda a step,
  `:inferred` growing, until `q` comes off -- `:entails?` true -- or the
  agenda is empty -- false."
  [kb q]
  (let [kb (vec kb)]
    (take-while some?
                (iterate (fn [{:keys [agenda count inferred] :as st}]
                           (when (nil? (:entails? st))
                             (if-let [p (peek agenda)]
                               (let [st (assoc st :agenda (pop agenda) :popped p)]
                                 (cond
                                   (= p q) (assoc st :entails? true)
                                   (inferred p) st
                                   :else
                                   (let [hit (filter #(contains? (:premise (kb %)) p) (range (clojure.core/count kb)))
                                         count (reduce (fn [c i] (update c i dec)) count hit)
                                         ready (filter #(zero? (count %)) hit)]
                                     (assoc st :inferred (conj inferred p) :count count
                                            :agenda (into (:agenda st) (map #(:conclusion (kb %)) ready))))))
                               (assoc st :entails? false))))
                         {:agenda (into empty-queue (map :conclusion (filter #(empty? (:premise %)) kb)))
                          :count (mapv #(clojure.core/count (:premise %)) kb)
                          :inferred #{}}))))

(defn pl-fc-entails?
  "Whether definite clauses `kb` entail the symbol `q`, by forward
  chaining."
  [kb q]
  (true? (:entails? (last (pl-fc-entails-seq kb q)))))

(defn pl-bc-entails-seq
  "Every machine state of backward chaining on definite clauses `kb` for
  the symbol `q`: to prove a goal, find a clause concluding it and prove
  each of its premises in turn, failing on a goal already being proved
  further up. The recursion runs on `allgo.search/call-seq`; each frame is
  a goal, and the last state's `:value` is whether `q` is proved."
  [kb q]
  (let [kb (vec kb)
        frame (fn [{:keys [goal path rules premises returned] :as f}]
                (cond
                  ;; First visit: the clauses that conclude the goal.
                  (nil? rules)
                  (if (contains? path goal)
                    [:return false]
                    (recur (assoc f :rules (filter #(= goal (:conclusion %)) kb))))

                  ;; A premise just proved: on to the next; one failed:
                  ;; on to the next clause.
                  (some? returned)
                  (if returned
                    (if (seq premises)
                      [:call {:goal (first premises) :path (conj path goal)}
                       (-> f (assoc :premises (rest premises)) (dissoc :returned))]
                      [:return true])
                    (recur (-> f (assoc :rules (rest rules) :premises nil) (dissoc :returned))))

                  :else
                  (if-let [r (first rules)]
                    (if-let [p (first (:premise r))]
                      [:call {:goal p :path (conj path goal)}
                       (assoc f :premises (rest (:premise r)) :rules (cons r (rest rules)))]
                      [:return true])
                    [:return false])))]
    (search/call-seq frame {:goal q :path #{}})))

(defn pl-bc-entails?
  "Whether definite clauses `kb` entail the symbol `q`, by backward
  chaining."
  [kb q]
  (true? (:value (last (pl-bc-entails-seq kb q)))))

;; ---------------------------------------------------------------------------
;; DPLL

(defn- assign
  "The clauses with `p` = `v`: those it makes true dropped, the literal it
  makes false removed from the rest -- the clauses under a partial model,
  which is what DPLL's tests of every clause amount to, kept as it goes
  rather than recomputed. nil if a clause is left empty: false."
  [clauses p v]
  (let [yes (if v p [:not p])
        no (negate yes)]
    (loop [cs clauses out (transient [])]
      (if-let [c (first cs)]
        (cond
          (contains? c yes) (recur (rest cs) out)
          (contains? c no) (let [c (disj c no)]
                             (when (seq c) (recur (rest cs) (conj! out c))))
          :else (recur (rest cs) (conj! out c)))
        (persistent! out)))))

(defn- find-pure-symbol
  "A symbol that appears with one sign only in the clauses not yet true,
  and that sign."
  [clauses]
  (let [lits (persistent! (reduce (fn [s c] (reduce conj! s c)) (transient #{}) clauses))]
    (some (fn [lit] (when-not (contains? lits (negate lit)) [(symbol-of lit) (keyword? lit)])) lits)))

(defn- find-unit-clause
  "A clause with one literal left: its symbol, and the value that makes it
  true."
  [clauses]
  (some (fn [c] (when (= 1 (count c)) (let [lit (first c)] [(symbol-of lit) (keyword? lit)]))) clauses))

(defn- dpll-frame
  [{:keys [clauses model returned tried split] :as f} ctx]
  (cond
    (true? returned) [:return true nil ctx]
    (false? returned) (if (= tried :split-true)
                        ;; The first branch failed: the other one.
                        (if-let [cs (assign clauses split false)]
                          [:call {:clauses cs :model (assoc model split false)}
                           (-> f (assoc :tried :split-false) (dissoc :returned)) ctx]
                          [:return false nil ctx])
                        [:return false nil ctx])
    :else
    (let [ctx (update ctx :calls inc)
          go (fn [p v kind]
               (if-let [cs (assign clauses p v)]
                 [:call {:clauses cs :model (assoc model p v)}
                  (assoc f :tried kind :split p) (update ctx kind inc)]
                 ;; The value empties a clause: for a pure symbol or a unit
                 ;; clause that cannot happen; for a split, the other value.
                 (if (= kind :split-true)
                   (recur p false :split-false)
                   [:return false nil ctx])))]
      (cond
        (empty? clauses) [:return true nil (assoc ctx :model model)]
        :else
        (if-let [[p v] (find-pure-symbol clauses)]
          (go p v :pure)
          (if-let [[p v] (find-unit-clause clauses)]
            (go p v :unit)
            (go (symbol-of (ffirst clauses)) true :split-true)))))))

(defn dpll-seq
  "Every machine state of DPLL on `clauses`, a step a call or a return:
  each frame a partial `:model` and the `:clauses` left under it, which
  a pure symbol, a unit clause or a split on a symbol -- true, then false
  -- extends. The context counts `:calls`, `:pure`, `:unit` and
  `:split-true`, and holds the satisfying `:model` when one is found; the
  last state's `:value` says whether.

  The figure tests every clause against the partial model on each call;
  this keeps the clauses simplified by the model instead, which is the
  same test made once, and splits on a symbol of the first clause left
  rather than the next in a fixed list -- a symbol that still matters."
  [clauses]
  (search/call-seq dpll-frame
                   {:clauses (vec (distinct clauses)) :model {}}
                   {:calls 0 :pure 0 :unit 0 :split-true 0 :split-false 0 :model nil}))

(defn dpll
  "A model satisfying `clauses`, or nil. Symbols the search never needed
  to set are left out of it."
  [clauses]
  (let [end (last (dpll-seq clauses))]
    (when (true? (:value end)) (:model (:context end)))))

(defn satisfiable? [& sentences] (some? (dpll (apply cnf sentences))))

(defn entails?
  "Whether `kb` entails `query`: whether KB & ~query is unsatisfiable."
  [kb query]
  (nil? (dpll (cnf kb [:not query]))))

;; ---------------------------------------------------------------------------
;; WalkSAT

(defn- walksat-step
  [{:keys [clauses model rng p flips max-flips] :as st}]
  (when (< (long flips) (long max-flips))
    (let [false-clauses (filterv #(false? (clause-true? % model)) clauses)]
      (when (seq false-clauses)
        (let [clause (random/pick rng false-clauses)
              syms (vec (map symbol-of clause))
              flip (fn [s] (update model s not))
              s (if (< (double (rng)) p)
                  (random/pick rng syms)
                  (apply max-key (fn [s] (let [m (flip s)] (count (filter #(true? (clause-true? % m)) clauses))))
                         syms))]
          (assoc st :model (flip s) :flipped s :flips (inc (long flips)) :unsatisfied (count false-clauses)))))))

(defn walksat-seq
  "Every state of WalkSAT on `clauses`: from a random model, a flip a
  step -- a symbol of a false clause, at random with probability `:p`
  (default 0.5), else the one whose flip satisfies most clauses -- until
  every clause is true or `:max-flips` (default 10000) are spent."
  ([clauses] (walksat-seq clauses {}))
  ([clauses {:keys [p max-flips rng]}]
   (let [rng (or rng (random/rng 1))
         syms (distinct (map symbol-of (mapcat identity clauses)))]
     (take-while some?
                 (iterate (fn [s] (when s (walksat-step s)))
                          {:clauses (vec clauses) :rng rng :p (or p 0.5) :max-flips (or max-flips 10000)
                           :flips 0 :model (into {} (map (fn [s] [s (< (double (rng)) 0.5)]) syms))})))))

(defn walksat
  "A model satisfying `clauses`, or nil if the flips run out -- which
  does not prove there is none."
  ([clauses] (walksat clauses {}))
  ([clauses opts]
   (let [m (:model (last (walksat-seq clauses opts)))]
     (when (every? #(true? (clause-true? % m)) clauses) m))))

(defn random-3-sat
  "`m` clauses of three distinct symbols each, out of `n`, each negated
  with probability one half."
  [rng n m]
  (vec (repeatedly m (fn []
                       (set (map (fn [k] (let [s (keyword (str "X" k))] (if (< (double (rng)) 0.5) s [:not s])))
                                 (take 3 (random/shuffle rng (vec (range n))))))))))
