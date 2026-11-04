(ns allgo.language.augmented
  "Augmented grammars: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 23.4 -- nonterminals that are
  not atomic symbols but terms with arguments, NP(Sbj, 1S, Speaker) where
  a context-free grammar has NP.

    agreement       arguments shared between a rule's categories must
                    unify: S(v) -> NP(Sbj, pn, n) VP(pn, v) makes the
                    subject's case Sbj and its person and number the
                    verb's (figure 23.9)
    lexicalized     a rule's probability a function of what its arguments
    PCFGs           are bound to -- VP(v) -> Verb(v) NP(n) [P1(v, n)], the
                    heads of the phrases -- so that \"ate a banana\" can be
                    likelier than \"ate a bandanna\"
    semantics       an argument for each phrase's meaning, built from its
                    parts' (compositional semantics): the value of an
                    arithmetic expression (figure 23.10), or a sentence of
                    first-order logic made by lambda-expressions applied
                    and beta-reduced (figure 23.12)

  A rule is `{:lhs term :rhs [item ...] :p p}`, each item a term or a word
  (a string), `:p` a probability or `(fn [bindings]) -> probability`
  (default 1). Terms are `allgo.logic.first-order`'s: a category
  `[:NP :Sbj pn n]`, variables lowercase symbols, constants keywords,
  strings or numbers, and any vector a compound term -- `[op x1 x2]`,
  `op` a variable, is the operator applied, as the chapter writes
  op(x1, x2). A grammar is `{:rules [...] :simplify f}`, `f` applied to
  every term built: an evaluator, a beta-reducer.

  Parsing is a chart parser, CYK's dynamic programming without its normal
  form: for each span, shortest first, every rule's right-hand side laid
  over every way of splitting the span among its items, each piece
  unified with an entry of the chart already found for it -- and, for a
  rule of one item, with entries of the same span, until none is new. Each
  entry is a term, the tree that made it, and its probability."
  (:require [allgo.logic.first-order :as fo]
            [clojure.walk :as walk]))

(defn rule
  "A rule `lhs -> rhs...` with probability 1."
  [lhs & rhs]
  {:lhs lhs :rhs (vec rhs) :p 1.0})

(defn- probability [p theta] (if (fn? p) (p theta) p))

(defn- splits
  "Every way to split [i, j) into `k` nonempty contiguous pieces, as the
  pieces' boundaries."
  [i j k]
  (if (= k 1)
    [[[i j]]]
    (for [m (range (inc i) (inc (- j (dec k)))) rest-splits (splits m j (dec k))]
      (into [[i m]] rest-splits))))

(def ^:private counter (atom 0))

(defn- standardized
  "Rule `r` with its variables renamed to fresh ones, and the renaming."
  [r]
  (let [k (swap! counter inc)
        renaming (into {} (map (fn [v] [v (symbol (str (name v) "_" k))]) (fo/variables [(:lhs r) (:rhs r)])))]
    [(assoc r :lhs (fo/subst renaming (:lhs r)) :rhs (fo/subst renaming (:rhs r))) renaming]))

(defn- matches
  "The ways rule `r`, standardized apart, covers [i, j): each `[r theta
  children p renaming]`."
  [r words chart i j]
  (let [[{:keys [rhs] :as r} renaming] (standardized r)]
    (for [pieces (splits i j (count rhs))
          [theta kids p] (reduce (fn [partials [item [a b]]]
                                   (if (string? item)
                                     (for [[theta kids p] partials :when (and (= (inc a) b) (= item (words a)))]
                                       [theta (conj kids item) p])
                                     (for [[theta kids p] partials
                                           {:keys [term tree] q :p} (chart [a b])
                                           :let [theta' (fo/unify item term theta)]
                                           :when theta']
                                       [theta' (conj kids tree) (* p q)])))
                                 [[{} [] 1.0]]
                                 (map vector rhs pieces))]
      [r theta kids p renaming])))

(defn- entry
  "The chart entry rule `r` makes with `theta`: its term, simplified, its
  tree, and its probability, the rule's -- given the bindings of the
  rule's own variables -- times its children's."
  [{:keys [simplify] :or {simplify identity}} r theta kids p renaming]
  (let [term (simplify (fo/subst theta (:lhs r)))
        bindings (into {} (map (fn [[v fresh]] [v (simplify (fo/subst theta fresh))]) renaming))]
    {:term term :tree (into [term] kids) :p (* p (double (probability (:p r) bindings)))}))

(defn chart-seq
  "Every state of chart parsing `words` with `grammar`: `{:chart :span}`,
  the chart mapping `[i j]` to the entries found for that span, a span a
  step, shortest first."
  [words {:keys [rules] :as grammar}]
  (let [words (vec words)
        n (count words)
        unary (filter #(= 1 (count (:rhs %))) rules)
        ;; Unary rules over a span's own entries, until nothing new: a
        ;; cycle of them would never end, and is cut off after ten rounds.
        close (fn [chart span]
                (loop [chart chart round 0]
                  (let [found (set (map :tree (chart span)))
                        new (for [r unary :when (not (string? (first (:rhs r))))
                                  [r theta kids p renaming] (matches r words chart (first span) (second span))
                                  :let [e (entry grammar r theta kids p renaming)]
                                  :when (not (found (:tree e)))]
                              e)]
                    (if (or (empty? new) (= round 10))
                      chart
                      (recur (update chart span (fnil into []) (distinct new)) (inc round))))))]
    (reductions (fn [{:keys [chart]} [i j]]
                  (let [direct (for [r rules
                                     :when (or (> (count (:rhs r)) 1) (string? (first (:rhs r))))
                                     [r theta kids p renaming] (matches r words chart i j)]
                                 (entry grammar r theta kids p renaming))]
                    {:chart (close (assoc chart [i j] (vec direct)) [i j]) :span [i j]}))
                {:chart {} :span nil}
                (for [len (range 1 (inc n)) i (range 0 (inc (- n len)))] [i (+ i len)]))))

(defn parses
  "Every parse of the whole of `words` whose term unifies with `root`
  (default `[:S x]`), most probable first: `[{:term :tree :p} ...]`."
  ([words grammar] (parses words grammar '[:S x]))
  ([words grammar root]
   (let [n (count words)
         chart (:chart (last (chart-seq words grammar)))]
     (sort-by (comp - :p) (filter #(fo/unify root (:term %)) (get chart [0 n]))))))

(defn parse
  "The most probable parse of `words`, or nil."
  ([words grammar] (first (parses words grammar)))
  ([words grammar root] (first (parses words grammar root))))

;; ---------------------------------------------------------------------------
;; Semantics

(defn evaluate-arithmetic
  "Evaluates every ground arithmetic term, `[op a b]` with `op` one of
  :+ :- :× :÷ and `a` and `b` numbers, to its number."
  [term]
  (walk/postwalk (fn [t]
                   (if (and (vector? t) (= 3 (count t)) (#{:+ :- :× :÷} (first t)) (number? (second t)) (number? (nth t 2)))
                     (let [[op a b] t] (case op :+ (+ a b) :- (- a b) :× (* a b) :÷ (/ a b)))
                     t))
                 term))

(defn- substitute
  "`body` with the variable `v` replaced by `x`, but not inside a lambda
  that binds `v` again."
  [body v x]
  (cond
    (= body v) x
    (and (vector? body) (= :lambda (first body)) (= v (second body))) body
    (vector? body) (mapv #(substitute % v x) body)
    :else body))

(defn beta-reduce
  "Every application of a lambda-expression in `term` -- `[[:lambda x
  body] arg]` -- replaced by `body` with `arg` for `x`, until none is
  left: (lambda x Loves(x, Bo))(Ali) becomes Loves(Ali, Bo)."
  [term]
  (if (vector? term)
    (let [t (mapv beta-reduce term)]
      (if (and (= 2 (count t)) (vector? (first t)) (= :lambda (ffirst t)))
        (let [[[_ v body] arg] t] (beta-reduce (substitute body v arg)))
        t))
    term))

;; ---------------------------------------------------------------------------
;; The chapter's grammars

(def arithmetic
  "Figure 23.10's grammar for arithmetic expressions, each phrase's
  semantics its value: Exp(op(x1, x2)) -> Exp(x1) Operator(op) Exp(x2),
  parenthesized expressions, numbers of digits."
  {:simplify evaluate-arithmetic
   :rules (vec (concat [(rule '[:Exp [op x1 x2]] '[:Exp x1] '[:Operator op] '[:Exp x2])
                        (rule '[:Exp x] "(" '[:Exp x] ")")
                        (rule '[:Exp x] '[:Number x])
                        (rule '[:Number x] '[:Digit x])
                        (rule '[:Number [:+ [:× 10 x1] x2]] '[:Number x1] '[:Digit x2])
                        (rule [:Operator :+] "+") (rule [:Operator :-] "-") (rule [:Operator :×] "×") (rule [:Operator :÷] "÷")]
                       (for [d (range 10)] (rule [:Digit d] (str d)))))})

(def ali-loves-bo
  "Figure 23.12's grammar, each phrase's semantics a term of first-order
  logic or a lambda-expression: S(pred(n)) -> NP(n) VP(pred), VP(pred(n))
  -> Verb(pred) NP(n), NP(n) -> Name(n), and \"loves\" lambda y lambda x
  Loves(x, y)."
  {:simplify beta-reduce
   :rules [(rule '[:S [pred n]] '[:NP n] '[:VP pred])
           (rule '[:VP [pred n]] '[:Verb pred] '[:NP n])
           (rule '[:NP n] '[:Name n])
           (rule [:Name :Ali] "Ali")
           (rule [:Name :Bo] "Bo")
           (rule '[:Verb [:lambda y [:lambda x [:Loves x y]]]] "loves")]})
