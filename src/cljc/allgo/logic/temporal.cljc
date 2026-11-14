(ns allgo.logic.temporal
  "Linear temporal logic: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 10.4.1.

    [:next p]        X P, P is true in the next time step
    [:eventually p]  F P, P is eventually (finally) true
    [:always p]      G P, P is always (globally) true
    [:until p q]     P U Q, P remains true until Q occurs
    [:release p q]   P R Q, Q holds until and including when P first does
                     (or forever) -- one of the operators the section
                     says can be derived: not (not P U not Q)

  and the propositional connectives `:not :and :or :implies`, atoms
  keywords.

  A behavior is infinite, so a trace is a lasso: a finite `:prefix` of
  states followed by a `:loop` repeated forever, each state the set of
  atoms true in it. On a lasso every formula is decided exactly, by
  labeling the finitely many positions -- each subformula innermost
  first, Until as the least set closed under q or (p and the next
  position in it), Release as the greatest, both reached by iterating to
  a fixed point.

  A system is a transition system: initial states, a successor function,
  and the atoms of each state. A formula holds of it when it holds of
  every behavior. `counterexample` looks for a behavior that falsifies it
  among every lasso up to a length -- bounded model checking: a found
  counterexample is genuine, none found is proof only up to that bound."
  (:require [allgo.search :as search]))

(defn- connective? [s] (and (vector? s) (#{:not :and :or :implies :next :eventually :always :until :release} (first s))))

(defn- subformulas
  "Every subformula of `f`, innermost first, each once."
  [f]
  (distinct (if (connective? f) (concat (mapcat subformulas (rest f)) [f]) [f])))

(defn positions
  "A lasso's positions, prefix then loop, and each one's successor: the
  last loop position's is the loop's first."
  [{:keys [prefix loop]}]
  (let [n (+ (count prefix) (count loop))]
    {:states (vec (concat prefix loop))
     :next (mapv #(if (= % (dec n)) (count prefix) (inc %)) (range n))}))

(defn- fixpoint-step
  "One step toward the fixed point of `f` from `xs`: nil once it stops
  changing."
  [f xs]
  (let [xs' (f xs)] (when (not= xs xs') xs')))

(defn- fixpoint [f start]
  (last (search/step-seq (partial fixpoint-step f) start)))

(defn- label
  "The positions where `f` holds, its subformulas' already in `labels`."
  [{:keys [states next]} labels f]
  (let [all (set (range (count states)))
        of labels]
    (if-not (connective? f)
      (into #{} (filter #(contains? (states %) f)) all)
      (let [[op p q] f]
        (case op
          :not (reduce disj all (of p))
          :and (reduce (fn [acc a] (into #{} (filter (of a)) acc)) all (rest f))
          :or (into #{} (mapcat of) (rest f))
          :implies (into #{} (filter #(or (not ((of p) %)) ((of q) %))) all)
          :next (into #{} (filter #((of p) (next %))) all)
          :until (fixpoint (fn [xs] (into (of q) (filter #(and ((of p) %) (xs (next %)))) all)) (of q))
          :release (fixpoint (fn [xs] (into #{} (filter #(and ((of q) %) (or ((of p) %) (xs (next %))))) all)) all)
          :eventually (fixpoint (fn [xs] (into (of p) (filter #(xs (next %))) all)) (of p))
          :always (fixpoint (fn [xs] (into #{} (filter #(and ((of p) %) (xs (next %)))) all)) all))))))

(defn- label-step [lasso labels f] (assoc labels f (label lasso labels f)))

(defn holds?
  "Whether `formula` holds of the behavior `lasso`, `{:prefix [states]
  :loop [states]}`, at its start."
  [lasso formula]
  (let [ps (positions lasso)
        labels (reduce (partial label-step ps) {} (subformulas formula))]
    (contains? (labels formula) 0)))

;; ---------------------------------------------------------- systems

(defn- lasso-of
  "The lasso a path of states makes when its last step returns to the
  state at index `k`."
  [system path k]
  (let [atoms (:atoms system)]
    {:prefix (mapv atoms (subvec path 0 k)) :loop (mapv atoms (subvec path k))}))

(defn- extend-step
  "One step of a depth-first search over paths: the next path off the
  stack extended by each successor -- closing a lasso if the successor is
  already on it -- or dropped at the bound."
  [{:keys [successors] :as system} formula bound {:keys [stack] :as st}]
  (when-let [path (peek stack)]
    (let [st (update st :stack pop)
          found (some (fn [s] (when-let [k (first (keep-indexed #(when (= %2 s) %1) path))]
                                (let [l (lasso-of system path k)]
                                  (when-not (holds? l formula) l))))
                      (successors (peek path)))]
      (cond
        found (assoc st :found found :stack [])
        (>= (count path) bound) st
        :else (update st :stack into (for [s (successors (peek path)) :when (not (some #{s} path))] (conj path s)))))))

(defn counterexample
  "A behavior of `system` -- `{:initial [states] :successors (fn [state])
  :atoms (fn [state])}` -- that falsifies `formula`, as a lasso of the
  atoms of its states, among every lasso of at most `bound` distinct
  states; nil if none is."
  [system formula bound]
  (:found (last (search/step-seq (partial extend-step system formula bound)
                                 {:stack (mapv vector (:initial system))}))))
