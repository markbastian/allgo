(ns allgo.search.csp
  "Constraint satisfaction: Russell and Norvig's chapter 6, as steps.

    AC-3                 figure 6.3: arcs revised one at a time from a
                         queue until every one is consistent, or a
                         domain is empty
    backtracking         figure 6.5, its recursion run by
                         `allgo.search/call-seq`, with the section's
                         choices: the first variable or the one with the
                         fewest values left (MRV, ties to the one in most
                         constraints); values in order or least
                         constraining first; and after each assignment no
                         inference, forward checking, or arc consistency
                         maintained by AC-3 (MAC)
    min-conflicts        figure 6.9: local search over complete
                         assignments, each step a conflicted variable set
                         to its least conflicted value
    tree solver          figure 6.11: a tree-structured CSP made arc
                         consistent from the leaves up, then assigned from
                         the root down, without backtracking
    cutset conditioning  the chapter's cycle cutsets: every consistent assignment to a
                         cycle cutset, each leaving a tree for the tree
                         solver -- a composition of its seqs

  Each is a step function run by `iterate`, its states a seq to watch,
  stop or count; each state has `:checks`, the constraint checks so far.

  A CSP is a map, made by `csp`:

    :variables    the variables, in order
    :domains      variable -> the values it may take, in order
    :neighbors    variable -> the variables it shares a binary constraint
                  with
    :constraint   (fn [a x b y]) -> whether a = x and b = y may hold
                  together, for neighbors a and b
    :constraints  higher-order constraints, each `{:scope [vars] :ok?
                  (fn [values])}`, checked once their scope is assigned;
                  AC-3 and the tree solver see only the binary ones"
  (:require [allgo.graph :as graph]
            [allgo.random :as random]
            [allgo.search :as search])
  #?(:clj (:import (clojure.lang PersistentQueue))))

(defn csp
  "A CSP from `:variables`, `:domains` (a map, or one seq for every
  variable), `:edges`, the pairs of variables bound by `:constraint`, and
  any higher-order `:constraints`."
  [{:keys [variables domains edges constraint constraints]}]
  {:variables (vec variables)
   :domains (if (map? domains) domains (zipmap variables (repeat (vec domains))))
   :neighbors (graph/adjacency variables edges)
   :constraint constraint
   :constraints (vec constraints)})

(def ^:private empty-queue #?(:clj PersistentQueue/EMPTY :cljs #queue []))

;; ---------------------------------------------------------------------------
;; Consistency

(defn- nary-ok?
  "Whether every higher-order constraint on `var` whose scope `assignment`
  fills holds."
  [{:keys [constraints]} assignment var]
  (every? (fn [{:keys [scope ok?]}]
            (or (not (some #{var} scope))
                (not (every? #(contains? assignment %) scope))
                (ok? (mapv assignment scope))))
          constraints))

(defn conflicts
  "How many constraints `var` = `value` breaks against `assignment`."
  [{:keys [neighbors constraint constraints]} assignment var value]
  (let [a (assoc assignment var value)]
    (+ (count (filter (fn [n] (and (contains? assignment n) (not= n var)
                                   (not (constraint var value n (assignment n)))))
                      (neighbors var)))
       (if (seq constraints)
         (count (remove (fn [{:keys [scope ok?]}]
                          (or (not (some #{var} scope))
                              (not (every? #(contains? a %) scope))
                              (ok? (mapv a scope))))
                        constraints))
         0))))

(defn consistent?
  "Whether `var` = `value` breaks no constraint against `assignment`."
  [{:keys [neighbors constraint] :as p} assignment var value]
  (and (every? (fn [n] (or (not (contains? assignment n)) (constraint var value n (assignment n))))
               (neighbors var))
       (nary-ok? p (assoc assignment var value) var)))

(defn solution?
  "Whether `assignment` gives every variable a value and breaks nothing."
  [{:keys [variables] :as p} assignment]
  (and (every? #(contains? assignment %) variables)
       (every? #(consistent? p (dissoc assignment %) % (assignment %)) variables)))

;; ---------------------------------------------------------------------------
;; AC-3

(defn- revise
  "REVISE: `domains` with every value of `xi` that no value of `xj`
  supports removed, whether any was, and the checks it took."
  [{:keys [constraint]} domains xi xj]
  (let [dj (domains xj)
        checks (volatile! 0)
        kept (filterv (fn [x] (some (fn [y] (vswap! checks inc) (constraint xi x xj y)) dj))
                      (domains xi))]
    [(assoc domains xi kept) (< (count kept) (count (domains xi))) @checks]))

(defn- ac-3-step
  [{:keys [queue queued neighbors] :as state}]
  (when-let [[xi xj :as arc] (peek queue)]
    (when-not (:failed? state)
      (let [[domains revised? n] (revise state (:domains state) xi xj)
            state (assoc state :domains domains :queue (pop queue) :queued (disj queued arc)
                         :arc arc :revised? revised? :checks (+ (:checks state) n))]
        (cond
          (not revised?) state
          (empty? (domains xi)) (assoc state :failed? true)
          :else (let [new (for [xk (neighbors xi) :when (not= xk xj)
                                :let [a [xk xi]] :when (not (contains? (:queued state) a))]
                            a)]
                  (assoc state :queue (into (:queue state) new) :queued (into (:queued state) new))))))))

(defn ac-3-seq
  "Every state of AC-3 on `p`, from `domains` (default its own) and the
  queue `arcs` (default every arc): one revision a step, `:arc` the one
  revised and `:revised?` whether it removed anything. The last has the
  arc-consistent `:domains`, or `:failed?`."
  ([p] (ac-3-seq p (:domains p) nil))
  ([{:keys [variables neighbors] :as p} domains arcs]
   (let [arcs (or arcs (for [x variables y (neighbors x)] [x y]))]
     (take-while some?
                 (iterate (fn [s] (when s (ac-3-step s)))
                          (assoc p :domains domains :queue (into empty-queue arcs) :queued (set arcs)
                                 :checks 0 :failed? false))))))

(defn ac-3
  "`p`'s domains made arc consistent, or nil if one comes up empty."
  ([p] (ac-3 p (:domains p) nil))
  ([p domains arcs]
   (let [end (last (ac-3-seq p domains arcs))]
     (when-not (:failed? end) (:domains end)))))

;; ---------------------------------------------------------------------------
;; Backtracking

(defn- select-variable
  "SELECT-UNASSIGNED-VARIABLE: the first in order, or with `:mrv` the one
  with fewest values left, ties broken by the most constraints on
  unassigned variables -- the degree heuristic."
  [{:keys [variables neighbors]} select assignment domains]
  (let [open (remove #(contains? assignment %) variables)]
    (if (= select :mrv)
      (apply min-key (fn [v] (+ (* 1000 (count (domains v)))
                                (- (count (remove #(contains? assignment %) (neighbors v))))))
             (reverse open))
      (first open))))

(defn- order-values
  "ORDER-DOMAIN-VALUES: as they are, or with `:lcv` the value that rules
  out fewest values of the unassigned neighbors first."
  [{:keys [neighbors constraint]} order var assignment domains]
  (if (= order :lcv)
    (sort-by (fn [x] (reduce + (for [n (neighbors var) :when (not (contains? assignment n))]
                                 (count (remove #(constraint var x n %) (domains n))))))
             (domains var))
    (domains var)))

(defn inference
  "INFERENCE after `var` = `value` (`kind` `:forward-checking`,
  `:mac` or anything else for none): the domains, narrowed, and the checks
  it took, or nil for a domain emptied. Forward checking removes from each
  unassigned neighbor the values that contradict the new one; MAC runs
  AC-3 from the arcs into `var`."
  [{:keys [neighbors constraint] :as p} kind assignment domains var value]
  (let [domains (assoc domains var [value])]
    (case kind
      :forward-checking
      (let [checks (volatile! 0)
            domains (reduce (fn [d n]
                              (if (contains? assignment n)
                                d
                                (let [kept (filterv (fn [y] (vswap! checks inc) (constraint var value n y)) (d n))]
                                  (if (empty? kept) (reduced nil) (assoc d n kept)))))
                            domains (neighbors var))]
        (when domains [domains @checks]))
      :mac
      (let [end (last (ac-3-seq p domains (for [n (neighbors var) :when (not (contains? assignment n))] [n var])))]
        (when-not (:failed? end) [(:domains end) (:checks end)]))
      [domains 0])))

(defn- backtrack-frame
  [{:keys [variables] :as p} {:keys [select order all?] infer :inference}]
  (fn [{:keys [assignment domains values var returned] :as frame} ctx]
    (cond
      ;; A child's solution is passed up.
      (and (some? returned) (not= ::failure returned)) [:return returned nil ctx]

      (and (nil? values) (= (count assignment) (count variables)))
      (if all?
        [:return ::failure nil (update ctx :solutions conj assignment)]
        [:return assignment nil ctx])

      :else
      (let [[var values] (if (nil? values)
                           (let [v (select-variable p select assignment domains)]
                             [v (order-values p order v assignment domains)])
                           [var values])]
        (loop [values values ctx ctx]
          (if-let [x (first values)]
            (let [ctx (update ctx :assignments inc)
                  ok? (consistent? p assignment var x)
                  ctx (update ctx :checks + (count (filter #(contains? assignment %) ((:neighbors p) var))))
                  inferred (when ok? (inference p infer assignment domains var x))]
              (if inferred
                (let [[d n] inferred]
                  [:call {:assignment (assoc assignment var x) :domains d}
                   (-> frame (assoc :var var :values (rest values)) (dissoc :returned))
                   (update ctx :checks + n)])
                (recur (rest values) ctx)))
            [:return ::failure nil (update ctx :backtracks inc)]))))))

(defn backtracking-seq
  "Every machine state of a backtracking search on `p`: one assignment
  tried a step, the stack of frames each holding its partial
  `:assignment` and pruned `:domains`. Options: `:select` `:first` or
  `:mrv`; `:order` `:given` or `:lcv`; `:inference` `:none`,
  `:forward-checking` or `:mac`; `:all?` to go on past the first solution
  and gather every one in the context's `:solutions`. The context counts
  `:assignments`, `:backtracks` and `:checks`."
  ([p] (backtracking-seq p {}))
  ([p opts]
   (map (fn [m] (assoc m :checks (:checks (:context m))))
        (search/call-seq (backtrack-frame p opts)
                         {:assignment {} :domains (:domains p)}
                         {:assignments 0 :backtracks 0 :checks 0 :solutions []}))))

(defn backtracking
  "A solution of `p`, or nil."
  ([p] (backtracking p {}))
  ([p opts]
   (let [v (:value (last (backtracking-seq p opts)))]
     (when (not= ::failure v) v))))

(defn all-solutions
  "Every solution of `p`, by backtracking with `opts` (default MRV and
  MAC: plain backtracking ran past ten minutes on a sudoku)."
  ([p] (all-solutions p {:select :mrv :inference :mac}))
  ([p opts] (:solutions (:context (last (backtracking-seq p (assoc opts :all? true)))))))

;; ---------------------------------------------------------------------------
;; Min-conflicts

(defn- conflicted [p assignment]
  (filterv #(pos? (conflicts p assignment % (assignment %))) (:variables p)))

(defn- min-conflicts-step
  [{:keys [current rng domains] :as state}]
  (let [bad (conflicted state current)]
    (if (empty? bad)
      nil
      (let [var (random/pick rng bad)
            scored (map (fn [x] [x (conflicts state current var x)]) (domains var))
            least (apply min (map second scored))
            x (random/pick rng (vec (keep (fn [[x c]] (when (= c least) x)) scored)))]
        (assoc state :current (assoc current var x) :var var :steps (inc (long (:steps state)))
               :conflicts (count bad))))))

(defn min-conflicts-seq
  "Every state of min-conflicts on `p`, one variable changed a step, from
  `:initial` -- default, each variable in turn given its least conflicted
  value, the book's greedy start -- until nothing conflicts or
  `:max-steps` (default 100000) pass. `:current` is the assignment, `:var`
  the variable changed and `:conflicts` how many were conflicted before."
  ([p] (min-conflicts-seq p {}))
  ([{:keys [variables domains] :as p} {:keys [rng initial max-steps]}]
   (let [rng (or rng (random/rng 1))
         initial (or initial
                     (reduce (fn [a v]
                               (let [scored (map (fn [x] [x (conflicts p a v x)]) (domains v))
                                     least (apply min (map second scored))]
                                 (assoc a v (random/pick rng (vec (keep (fn [[x c]] (when (= c least) x)) scored))))))
                             {} variables))
         max-steps (or max-steps 100000)]
     (take-while some?
                 (iterate (fn [s] (when (and s (< (long (:steps s)) (long max-steps))) (min-conflicts-step s)))
                          (assoc p :current initial :rng rng :steps 0 :var nil))))))

(defn min-conflicts
  "A solution of `p` by min-conflicts, or nil if `:max-steps` run out."
  ([p] (min-conflicts p {}))
  ([p opts]
   (let [end (:current (last (min-conflicts-seq p opts)))]
     (when (solution? p end) end))))

;; ---------------------------------------------------------------------------
;; Tree-structured CSPs

(defn tree?
  "Whether `p`'s constraint graph has no cycles -- a forest, which the tree
  solver takes one tree at a time."
  [{:keys [variables neighbors]}]
  (let [edges (set (for [x variables y (neighbors x)] #{x y}))]
    (= (count edges) (- (count variables) (count (graph/components variables edges))))))

(defn- topological
  "The variables root first, each tree breadth first, and each one's
  parent."
  [{:keys [variables neighbors]}]
  (loop [order [] parent {} todo variables]
    (if-let [root (first (remove #(contains? parent %) todo))]
      (let [[order parent]
            (loop [queue (conj empty-queue root) order order parent (assoc parent root nil)]
              (if-let [x (peek queue)]
                (let [kids (remove #(contains? parent %) (neighbors x))]
                  (recur (into (pop queue) kids) (conj order x) (into parent (map (fn [k] [k x]) kids))))
                [order parent]))]
        (recur order parent (rest (drop-while #(not= % root) todo))))
      [order parent])))

(defn- tree-step
  [{:keys [order parent domains phase j assignment constraint] :as state}]
  (when-not (or (:failed? state) (:solved? state))
    (case phase
      ;; From the last variable up, each parent made consistent with its
      ;; child.
      :arcs
      (if (zero? (long j))
        (assoc state :phase :assign :j 0)
        (let [x (order j)
              par (parent x)]
          (if (nil? par)
            (assoc state :j (dec (long j)) :arc nil)
            (let [[d _ n] (revise state domains par x)]
              (cond-> (assoc state :domains d :j (dec (long j)) :arc [par x] :checks (+ (:checks state) n))
                (empty? (d par)) (assoc :failed? true))))))
      ;; From the root down, any value consistent with the parent's.
      :assign
      (if (= (long j) (count order))
        (assoc state :solved? true)
        (let [x (order j)
              par (parent x)
              v (first (filter (fn [y] (or (nil? par) (constraint par (assignment par) x y))) (domains x)))]
          (if (nil? v)
            (assoc state :failed? true)
            (assoc state :assignment (assoc assignment x v) :j (inc (long j)) :var x
                   :checks (+ (:checks state) (if par 1 0)))))))))

(defn tree-solver-seq
  "Every state of the tree solver on a tree-structured `p`: the arcs made
  consistent from the leaves up, one a step, then the variables assigned
  from the roots down, one a step, never backtracking. The last state is
  `:solved?`, with the `:assignment`, or `:failed?`."
  ([p] (tree-solver-seq p (:domains p)))
  ([p domains]
   (let [[order parent] (topological p)]
     (take-while some?
                 (iterate (fn [s] (when s (tree-step s)))
                          (assoc p :order order :parent parent :domains domains :phase :arcs
                                 :j (dec (count order)) :assignment {} :checks 0
                                 :solved? false :failed? false))))))

(defn tree-solver
  "A solution of a tree-structured `p`, or nil."
  [p]
  (let [end (last (tree-solver-seq p))]
    (when (:solved? end) (:assignment end))))

;; ---------------------------------------------------------------------------
;; Cutset conditioning

(defn cycle-cutset
  "Variables whose removal leaves `p`'s constraint graph a forest: the most
  constrained taken greedily until it is."
  [{:keys [variables neighbors] :as p}]
  (loop [cut [] left (set variables)]
    (let [sub (assoc p :variables (filterv left variables)
                     :neighbors (into {} (for [x left] [x (set (filter left (neighbors x)))])))]
      (if (tree? sub)
        cut
        (let [x (apply max-key #(count (filter left (neighbors %))) (filter left variables))]
          (recur (conj cut x) (disj left x)))))))

(defn cutset-conditioning-seq
  "Every state of cutset conditioning on `p`: for each consistent
  assignment to `cutset` (default `cycle-cutset`), in turn, the rest's
  domains pruned to agree with it and the tree solver run on what is
  left, its states tagged with the `:cutset-assignment` -- until one
  succeeds. A state that is `:solved?` has the whole `:assignment`."
  ([p] (cutset-conditioning-seq p (cycle-cutset p)))
  ([{:keys [variables neighbors domains constraint] :as p} cutset]
   (let [cut (set cutset)
         rest-vars (filterv (complement cut) variables)
         sub (assoc p :variables rest-vars
                    :neighbors (into {} (for [x rest-vars] [x (set (remove cut (neighbors x)))])))
         ;; Whether the values of `x` outside the cutset leave any agreeing
         ;; with every cutset neighbor assigned so far.
         supported? (fn [ca x]
                      (some (fn [y] (every? (fn [c] (or (not (contains? ca c)) (constraint x y c (ca c))))
                                            (neighbors x)))
                            (domains x)))
         ;; The consistent assignments to the cutset, depth first from an
         ;; explicit stack, a branch dropped as soon as it leaves a variable
         ;; outside the cutset with no value -- the space is exponential in
         ;; the cutset, and this is what keeps a small one cheap.
         assignments (fn assignments [stack]
                       (lazy-seq
                        (when-let [[vs a] (peek stack)]
                          (let [stack (pop stack)]
                            (if-let [v (first vs)]
                              (assignments
                               (into stack
                                     (for [x (reverse (domains v))
                                           :let [a (assoc a v x)]
                                           :when (and (consistent? p (dissoc a v) v x)
                                                      (every? #(or (cut %) (supported? a %)) (neighbors v)))]
                                       [(rest vs) a])))
                              (cons a (assignments stack)))))))]
     (letfn [(runs [cas]
               (lazy-seq
                (when-let [ca (first cas)]
                  (let [pruned (into domains
                                     (for [x rest-vars]
                                       [x (filterv (fn [y] (every? (fn [c] (or (not (contains? ca c))
                                                                               (constraint x y c (ca c))))
                                                                   (neighbors x)))
                                                   (domains x))]))
                        states (map #(assoc % :cutset-assignment ca
                                            :assignment (merge ca (:assignment %)))
                                    (tree-solver-seq sub pruned))]
                    (if (:solved? (last states))
                      states
                      (concat states (runs (rest cas))))))))]
       (runs (assignments [[cutset {}]]))))))

(defn cutset-conditioning
  "A solution of `p` by cutset conditioning, or nil."
  [p]
  (let [end (last (cutset-conditioning-seq p))]
    (when (:solved? end) (:assignment end))))
