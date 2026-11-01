(ns allgo.planning
  "Automated planning: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 11.

    classical planning    PDDL-style problems -- states sets of ground
                          atoms, action schemas with variables, preconditions
                          and effects -- grounded by the unifier of
                          `allgo.logic.first-order` and searched forward
                          (progression) by `allgo.search`'s A*, the
                          heuristic the number of goals unmet
    hierarchical search   HIERARCHICAL-SEARCH (figure 11.8): plans of
                          high-level actions refined breadth first until
                          one of primitive actions reaches the goal
    angelic search        ANGELIC-SEARCH (figure 11.11): high-level plans
                          judged by the states they can reach -- some of
                          them for certain, others perhaps -- committed to
                          once they surely reach the goal, then decomposed
    scheduling            the critical path method, and scheduling under
                          resource constraints by minimum slack (the
                          chapter's job-shop problem, figure 11.13)

  Atoms are vectors `[:At :C1 :SFO]`, as in `allgo.logic.first-order`;
  variables lowercase symbols. A literal is an atom or `[:not atom]`, and
  `[:neq a b]` in a precondition requires two terms to differ."
  (:require [allgo.logic.first-order :as fo]
            [allgo.search :as search]))

;; ---------------------------------------------------------------------------
;; Classical planning

(defn- negative? [lit] (= :not (first lit)))

(defn- ground-actions
  "Every ground instance of `schema` whose positive preconditions hold in
  `state`, its variables bound by matching them against the state's
  atoms."
  [{:keys [name pre effect]} state]
  (let [positive (remove #(or (negative? %) (= :neq (first %))) pre)
        others (filter #(or (negative? %) (= :neq (first %))) pre)
        thetas (reduce (fn [thetas p]
                         (for [th thetas a state :let [th' (fo/unify (fo/subst th p) a th)] :when th'] th'))
                       [{}] positive)]
    (distinct
     (for [th thetas
           :let [ok? (every? (fn [lit]
                               (let [l (fo/subst th lit)]
                                 (cond (= :neq (first l)) (not= (nth l 1) (nth l 2))
                                       :else (not (contains? state (second l))))))
                             others)]
           :when (and ok? (empty? (fo/variables (fo/subst th name))))]
       {:name (fo/subst th name) :effect (mapv #(fo/subst th %) effect)}))))

(defn result
  "The state after a ground action: delete list first, then add list."
  [state {:keys [effect]}]
  (let [dels (map second (filter negative? effect))
        adds (remove negative? effect)]
    (into (apply disj state dels) adds)))

(defn holds? [state goal]
  (every? (fn [lit] (if (negative? lit) (not (contains? state (second lit))) (contains? state lit))) goal))

(defn applicable
  "Every ground action applicable in `state`."
  [{:keys [actions]} state]
  (mapcat #(ground-actions % state) actions))

(defn plan
  "A plan for `problem` by forward state-space search with A*, each action
  costing one and the heuristic the number of goal literals unmet -- the
  'ignore preconditions' heuristic when an action meets one goal at a
  time. A vector of ground action names, or nil."
  [{:keys [init goal] :as problem}]
  (let [succ (memoize (fn [state] (map (fn [a] [(result state a) (:name a)]) (applicable problem state))))
        path (search/a-star {:start (set init)
                             :goal? #(holds? % goal)
                             :neighbors (fn [s] (map first (succ s)))
                             :heuristic (fn [s _] (count (remove (fn [lit] (holds? s [lit])) goal)))})]
    (when path
      (vec (for [[a b] (partition 2 1 path)]
             (some (fn [[s name]] (when (= s b) name)) (succ a)))))))

(defn execute
  "The state after the plan of ground action `names` from `init`, or nil if
  some action's preconditions do not hold."
  [{:keys [init] :as problem} names]
  (reduce (fn [state n]
            (if-let [a (first (filter #(= n (:name %)) (applicable problem state)))]
              (result state a)
              (reduced nil)))
          (set init) names))

;; ---------------------------------------------------------------------------
;; The chapter's problems, as aima-pseudocode (CC BY 4.0) states them

(defn- parse-all [texts] (mapv fo/parse texts))

(def air-cargo
  {:init (set (parse-all ["At(C1, SFO)" "At(C2, JFK)" "At(P1, SFO)" "At(P2, JFK)" "Cargo(C1)" "Cargo(C2)"
                          "Plane(P1)" "Plane(P2)" "Airport(JFK)" "Airport(SFO)"]))
   :goal (parse-all ["At(C1, JFK)" "At(C2, SFO)"])
   :actions [{:name (fo/parse "Load(c, p, a)")
              :pre (parse-all ["At(c, a)" "At(p, a)" "Cargo(c)" "Plane(p)" "Airport(a)"])
              :effect (parse-all ["~At(c, a)" "In(c, p)"])}
             {:name (fo/parse "Unload(c, p, a)")
              :pre (parse-all ["In(c, p)" "At(p, a)" "Cargo(c)" "Plane(p)" "Airport(a)"])
              :effect (parse-all ["At(c, a)" "~In(c, p)"])}
             {:name (fo/parse "Fly(p, from, to)")
              :pre (into (parse-all ["At(p, from)" "Plane(p)" "Airport(from)" "Airport(to)"]) [[:neq 'from 'to]])
              :effect (parse-all ["~At(p, from)" "At(p, to)"])}]})

(def spare-tire
  {:init (set (parse-all ["Tire(Flat)" "Tire(Spare)" "At(Flat, Axle)" "At(Spare, Trunk)"]))
   :goal (parse-all ["At(Spare, Axle)"])
   :actions [{:name (fo/parse "Remove(obj, loc)")
              :pre (into (parse-all ["At(obj, loc)"]) [[:neq 'loc :Ground]])
              :effect (parse-all ["~At(obj, loc)" "At(obj, Ground)"])}
             {:name (fo/parse "PutOn(t, Axle)")
              :pre (parse-all ["Tire(t)" "At(t, Ground)" "~At(Flat, Axle)"])
              :effect (parse-all ["~At(t, Ground)" "At(t, Axle)"])}
             {:name [:LeaveOvernight]
              :pre []
              :effect (parse-all ["~At(Spare, Ground)" "~At(Spare, Axle)" "~At(Spare, Trunk)"
                                  "~At(Flat, Ground)" "~At(Flat, Axle)" "~At(Flat, Trunk)"])}]})

(def blocks-world
  {:init (set (parse-all ["On(A, Table)" "On(B, Table)" "On(C, A)" "Block(A)" "Block(B)" "Block(C)" "Clear(B)" "Clear(C)"]))
   :goal (parse-all ["On(A, B)" "On(B, C)"])
   :actions [{:name (fo/parse "Move(b, x, y)")
              :pre (into (parse-all ["On(b, x)" "Clear(b)" "Clear(y)" "Block(b)" "Block(y)"])
                         [[:neq 'b 'x] [:neq 'b 'y] [:neq 'x 'y]])
              :effect (parse-all ["On(b, y)" "Clear(x)" "~On(b, x)" "~Clear(y)"])}
             {:name (fo/parse "MoveToTable(b, x)")
              :pre (into (parse-all ["On(b, x)" "Clear(b)" "Block(b)"]) [[:neq 'b 'x]])
              :effect (parse-all ["On(b, Table)" "Clear(x)" "~On(b, x)"])}]})

;; ---------------------------------------------------------------------------
;; Hierarchical search

;; A hierarchy maps the name of each high-level action to a function of
;; the HLA's arguments and the state, returning its refinements whose
;; preconditions hold there -- each a sequence of actions, primitive or
;; not. A primitive action is any name not in the hierarchy; `step` (fn
;; [state name]) -> the next state, or nil if it cannot be done.

(defn- hla? [hierarchy action] (contains? hierarchy (first action)))

(defn- run-primitive [step state actions]
  (reduce (fn [s a] (or (step s a) (reduced nil))) state actions))

(defn hierarchical-search-seq
  "HIERARCHICAL-SEARCH, a step a plan taken from the queue: the first HLA in
  it refined every way its hierarchy allows in the state the primitive
  prefix reaches, until a primitive plan reaches a state satisfying
  `goal?`. Each state has the plan popped, `:plan`, and at the end
  `:solution`."
  [{:keys [init goal? step]} hierarchy top]
  (take-while some?
              (iterate (fn [{:keys [frontier] :as st}]
                         (when (and (seq frontier) (not (:solution st)))
                           (let [plan (peek frontier)
                                 frontier (pop frontier)
                                 i (first (keep-indexed (fn [i a] (when (hla? hierarchy a) i)) plan))
                                 prefix (if i (subvec plan 0 i) plan)
                                 outcome (run-primitive step init prefix)]
                             (cond
                               (nil? outcome) (assoc st :frontier frontier :plan plan)
                               (nil? i) (cond-> (assoc st :frontier frontier :plan plan)
                                          (goal? outcome) (assoc :solution plan))
                               :else
                               (let [hla (plan i)
                                     suffix (subvec plan (inc i))]
                                 (assoc st :plan plan
                                        :frontier (into frontier (for [r ((hierarchy (first hla)) (vec (rest hla)) outcome)]
                                                                   (vec (concat prefix r suffix))))))))))
                       {:frontier (conj #?(:clj clojure.lang.PersistentQueue/EMPTY :cljs #queue []) [top])})))

(defn hierarchical-search
  "A primitive plan refined from `top`, or nil."
  [problem hierarchy top]
  (:solution (last (hierarchical-search-seq problem hierarchy top))))

;; ---------------------------------------------------------------------------
;; Angelic search

(defn reachable
  "The states a plan can reach from `state`, each HLA standing for every
  primitive refinement of it -- found by refining, down to `depth`
  levels -- so that REACH+ and REACH- are this exact set. Supplying
  approximate descriptions instead (an optimistic superset, a pessimistic
  subset) is what makes angelic search cheaper than this."
  [{:keys [step]} hierarchy state plan depth]
  (if-let [a (first plan)]
    (if (hla? hierarchy a)
      (if (zero? depth)
        #{}
        (set (mapcat #(reachable {:step step} hierarchy state (concat % (rest plan)) (dec depth))
                     ((hierarchy (first a)) (vec (rest a)) state))))
      (if-let [s (step state a)] (reachable {:step step} hierarchy s (rest plan) depth) #{}))
    #{state}))

(defn angelic-search
  "ANGELIC-SEARCH: high-level plans refined breadth first, but judged by
  their reachable sets: a plan that cannot reach the goal is dropped, one
  that surely reaches it is committed to and decomposed -- each HLA in it
  solved, back to front, as a problem of its own, from a state its prefix
  surely reaches to the state the rest needs. Here the reachable sets are
  exact, `reachable` to `depth` levels, so REACH+ = REACH-. A primitive
  plan, or nil."
  [{:keys [init goal?] :as problem} hierarchy initial-plan depth]
  (loop [frontier (conj #?(:clj clojure.lang.PersistentQueue/EMPTY :cljs #queue []) (vec initial-plan))]
    (when-let [plan (peek frontier)]
      (let [frontier (pop frontier)
            reach (reachable problem hierarchy init plan depth)
            goals (filter goal? reach)]
        (cond
          (empty? goals) (recur frontier)
          (not-any? #(hla? hierarchy %) plan) plan
          ;; MAKING-PROGRESS: a plan refined since the call began is
          ;; committed to and decomposed -- the last action first, from a
          ;; state the rest of the plan reaches and from which it reaches
          ;; the final state. Decomposing the initial plan itself would
          ;; only call this again on the same HLA, for ever.
          (not= plan (vec initial-plan))
          (let [final (first goals)]
            (or (loop [plan plan sf final solution ()]
                  (if (empty? plan)
                    (vec (apply concat solution))
                    (let [action (peek plan)
                          rest-plan (pop plan)
                          si (first (filter #(contains? (reachable problem hierarchy % [action] depth) sf)
                                            (reachable problem hierarchy init rest-plan depth)))
                          sub (when si (if (hla? hierarchy action)
                                         (angelic-search (assoc problem :init si :goal? #(= % sf)) hierarchy [action] depth)
                                         [action]))]
                      (when sub (recur rest-plan si (cons sub solution))))))
                (recur frontier)))

          ;; Otherwise refine the first HLA, its primitive prefix run to
          ;; find the state it is refined in.
          :else
          (let [i (first (keep-indexed (fn [i a] (when (hla? hierarchy a) i)) plan))
                prefix (subvec plan 0 i)
                outcome (run-primitive (:step problem) init prefix)
                hla (plan i)]
            (recur (if outcome
                     (into frontier (for [r ((hierarchy (first hla)) (vec (rest hla)) outcome)]
                                      (vec (concat prefix r (subvec plan (inc i))))))
                     frontier))))))))

;; ---------------------------------------------------------------------------
;; Scheduling

(defn critical-path
  "The critical path method: for actions `{name {:duration :after
  #{names}}}`, each one's earliest and latest start `[ES LS]`, and the
  length of the shortest schedule. The actions whose ES = LS are the
  critical path; LS - ES is the slack."
  [actions]
  (let [order (loop [done [] left (set (keys actions))]
                (if (empty? left) done
                    (let [ready (filter #(every? (set done) (:after (actions %))) (sort left))]
                      (recur (into done ready) (apply disj left ready)))))
        es (reduce (fn [es a] (assoc es a (reduce max 0 (map #(+ (es %) (:duration (actions %))) (:after (actions a))))))
                   {} order)
        length (reduce max (map #(+ (es %) (:duration (actions %))) order))
        successors (fn [a] (filter #(contains? (:after (actions %)) a) order))
        ls (reduce (fn [ls a] (assoc ls a (- (reduce min length (map ls (successors a))) (:duration (actions a)))))
                   {} (reverse order))]
    {:start (into {} (for [a order] [a [(es a) (ls a)]])) :length length}))

(defn schedule
  "Scheduling with reusable resources by the minimum-slack heuristic: on
  each step, of the actions whose predecessors are all scheduled, the one
  of least slack (by the critical path method on what remains) is placed
  at the earliest time its predecessors are done and its resources --
  `{:uses {resource n}}` against `capacity` -- are free. `{name start}`
  and the schedule's length."
  [actions capacity]
  (loop [placed {} busy []]
    (if (= (count placed) (count actions))
      {:start placed :length (reduce max (map (fn [[a s]] (+ s (:duration (actions a)))) placed))}
      (let [ready (filter (fn [a] (and (not (contains? placed a)) (every? #(contains? placed %) (:after (actions a)))))
                          (keys actions))
            slack (:start (critical-path actions))
            a (apply min-key (fn [a] (let [[es ls] (slack a)] (- ls es))) (sort ready))
            {:keys [duration after uses]} (actions a)
            earliest (reduce max 0 (map #(+ (placed %) (:duration (actions %))) after))
            free? (fn [t] (every? (fn [[r n]]
                                    (<= (+ n (reduce + (for [[s e rs] busy :when (and (< s (+ t duration)) (< t e))] (get rs r 0))))
                                        (capacity r)))
                                  uses))
            t (first (filter free? (sort (distinct (cons earliest (filter #(>= % earliest) (map second busy)))))))]
        (recur (assoc placed a t) (conj busy [t (+ t duration) uses]))))))

(defn- place
  "`a` placed at the earliest time its predecessors are done and its
  resources free, given what is `placed` and `busy`."
  [actions capacity placed busy a]
  (let [{:keys [duration after uses]} (actions a)
        earliest (reduce max 0 (map #(+ (placed %) (:duration (actions %))) after))
        free? (fn [t] (every? (fn [[r n]]
                                (<= (+ n (reduce + (for [[s e rs] busy :when (and (< s (+ t duration)) (< t e))] (get rs r 0))))
                                    (capacity r)))
                              uses))]
    (first (filter free? (sort (distinct (cons earliest (filter #(>= % earliest) (map second busy)))))))))

(defn optimal-schedule
  "The shortest schedule, by trying every order of the actions that their
  precedences allow, each action placed as early as it can go in turn,
  and pruning an order once it is already no shorter than the best found:
  exponential, for small problems, and the yardstick for `schedule`."
  [actions capacity]
  (let [best (volatile! {:length ##Inf})]
    (letfn [(go [placed busy]
                (let [length (reduce max 0 (map (fn [[a s]] (+ s (:duration (actions a)))) placed))]
                  (cond
                    (>= length (:length @best)) nil
                    (= (count placed) (count actions)) (vreset! best {:start placed :length length})
                    :else (doseq [a (sort (filter (fn [a] (and (not (contains? placed a))
                                                               (every? #(contains? placed %) (:after (actions a)))))
                                                  (keys actions)))]
                            (let [t (place actions capacity placed busy a)
                                  {:keys [duration uses]} (actions a)]
                              (go (assoc placed a t) (conj busy [t (+ t duration) uses])))))))]
      (go {} [])
      @best)))

(def car-assembly
  "The chapter's job-shop problem: two cars, each an engine, then wheels,
  then an inspection, as aima-pseudocode states it. The lug nuts it
  consumes are plentiful and left out."
  {:actions {:AddEngine1 {:duration 30 :after #{} :uses {:EngineHoists 1}}
             :AddEngine2 {:duration 60 :after #{} :uses {:EngineHoists 1}}
             :AddWheels1 {:duration 30 :after #{:AddEngine1} :uses {:WheelStations 1}}
             :AddWheels2 {:duration 15 :after #{:AddEngine2} :uses {:WheelStations 1}}
             :Inspect1 {:duration 10 :after #{:AddWheels1} :uses {:Inspectors 1}}
             :Inspect2 {:duration 10 :after #{:AddWheels2} :uses {:Inspectors 1}}}
   :capacity {:EngineHoists 1 :WheelStations 1 :Inspectors 2}})
