(ns allgo.search.nondeterministic
  "Search where an action's outcome is not known, or the state is not:
  Russell and Norvig's chapter 4, on nondeterministic actions and partial
  observability.

    and-or-seq    AND-OR search (figure 4.10, AND-OR-GRAPH-SEARCH): a
                  plan that reaches a goal whatever the actions turn out
                  to do. Its OR-SEARCH and AND-SEARCH call each other; the
                  recursion is run a step at a time by
                  `allgo.search/call-seq`, so it is a seq of machine
                  states like every other search here
    cyclic plans  the same search allowed to loop back to a state on its
                  path, for actions that may do nothing at all
    partial       percepts split belief states, and AND-OR search plans
    observability on each
    sensorless    a problem with no percepts at all, made a search over
                  belief states -- the sets of states the agent might be
                  in -- which `allgo.search` then solves like any other

  A problem is a map:

    :initial    the state to begin from
    :actions    (fn [state]) -> the actions possible in it
    :results    (fn [state action]) -> the set of states it may lead to
    :goal?      (fn [state]) -> whether it is a goal

  A conditional plan is `[]`, nothing left to do, or `[action plans]`,
  where `plans` maps each state the action may lead to onto the plan to
  follow from there -- the book's [action, if s1 then plan1 else ...]."
  (:require [allgo.search :as search]))

;; ---------------------------------------------------------------------------
;; AND-OR search

(defn- or-frame [{:keys [actions]} state path]
  {:kind :or :state state :path path :actions (seq (actions state))})

(defn- reaches-goal?
  "Whether a goal can be reached from every point of a cyclic plan: some
  outcome of every action leads on without only looping -- which, if
  every outcome eventually happens, makes it a solution."
  [plan]
  (cond
    (empty? plan) true
    (= :loop (first plan)) false
    :else (let [[_ plans] plan
                subs (vals plans)]
            (and (some reaches-goal? subs)
                 (every? #(or (= :loop (first %)) (reaches-goal? %)) subs)))))

(defn- and-or-frame
  "OR-SEARCH and AND-SEARCH, as frames. An OR frame tries its state's
  actions in turn, each an AND frame over the states the action may lead
  to, and returns the first action every one of those can be planned
  from. An AND frame plans from each of its states in turn and fails if
  any fails."
  [{:keys [results goal? cyclic?] :as problem}]
  (fn [{:keys [kind state path actions action returned states i plans] :as frame}]
    (case kind
      :or
      (cond
        (and (nil? action) (goal? state)) [:return []]
        ;; A state already on the path: a cycle. Plain AND-OR search fails
        ;; on it; allowing cyclic plans, it is a loop back to that state.
        (and (nil? action) (contains? path state)) [:return (if cyclic? [:loop state] ::failure)]
        ;; With cycles allowed, a plan whose outcomes all loop back is
        ;; no plan -- sucking a clean square for ever -- so the next
        ;; action is tried instead.
        (and returned (not= ::failure returned)
             (or (not cyclic?) (reaches-goal? [action returned])))
        [:return [action returned]]
        (seq actions)
        (let [a (first actions)]
          [:call {:kind :and :states (vec (results state a)) :path (conj path state) :i 0 :plans {}}
           (-> frame (assoc :actions (rest actions) :action a) (dissoc :returned))])
        :else [:return ::failure])

      :and
      (cond
        (= ::failure returned) [:return ::failure]
        :else
        (let [plans (if (some? returned) (assoc plans (states (dec i)) returned) plans)]
          (if (= i (count states))
            [:return plans]
            [:call (or-frame problem (states i) path)
             (-> frame (assoc :i (inc i) :plans plans) (dissoc :returned))]))))))

(defn and-or-seq
  "Every machine state of an AND-OR search, one per call or return; the
  last has the plan as its `:value`, or `:allgo.search.nondeterministic/failure`.
  It fails on a cycle rather than plan one, as the book's does: a plan
  that says 'try again until it works' is the chapter's cyclic plan,
  not this."
  [problem]
  (search/call-seq (and-or-frame problem) (or-frame problem (:initial problem) #{})))

(defn and-or
  "A conditional plan that reaches a goal however the actions turn out,
  or nil if there is none without cycles."
  [problem]
  (let [v (:value (last (and-or-seq problem)))]
    (when (not= ::failure v) v)))

(defn- labels
  "Every state a plan acts in, mapped to the plan from there: the labels
  a loop jumps back to."
  [plan state]
  (if (or (empty? plan) (= :loop (first plan)))
    {}
    (let [[_ plans] plan]
      (into {state plan} (mapcat (fn [[s p]] (labels p s)) plans)))))

(defn execute
  "The states a plan passes through when the world picks each outcome
  with `choose` -- (fn [states]) -> one of them -- from `state`, at most
  `max-steps` (default 1000) of them: a cyclic plan's loops jump back to
  the plan for the state they name."
  ([problem plan state choose] (execute problem plan state choose 1000))
  ([{:keys [results]} plan state choose max-steps]
   (let [label (labels plan state)]
     (loop [plan plan state state seen [state]]
       (cond
         (empty? plan) seen
         (>= (count seen) max-steps) seen
         (= :loop (first plan)) (recur (label (second plan)) state seen)
         :else (let [[action plans] plan
                     nxt (choose (vec (results state action)))]
                 (recur (plans nxt) nxt (conj seen nxt))))))))

(defn cyclic-and-or
  "A cyclic plan: AND-OR search allowed to loop back to a state already on
  its path -- `[:loop state]` -- for worlds where an action may do
  nothing and the only plan is to try again. nil if there is none, or the
  one found could loop for ever without reaching the goal."
  [problem]
  (let [v (:value (last (and-or-seq (assoc problem :cyclic? true))))]
    (when (and (not= ::failure v) (reaches-goal? v)) v)))

;; ---------------------------------------------------------------------------
;; Partial observability

(defn partially-observable
  "A problem whose agent senses something of each state -- `:percept`
  (fn [state]) -> what it perceives there -- as an AND-OR problem over
  belief states: from `:initial`, a set of states, an action leads to the
  states it may (PREDICT), which the percepts it might then get split into
  beliefs (POSSIBLE-PERCEPTS, UPDATE); a belief is a goal when every state
  in it is. `and-or` then plans what to do on each percept."
  [{:keys [initial actions results percept goal?]}]
  {:initial (set initial)
   :actions (fn [belief] (distinct (mapcat actions belief)))
   :results (fn [belief a]
              (let [predicted (set (mapcat #(results % a) belief))]
                (set (map set (vals (group-by percept predicted))))))
   :goal? (fn [belief] (every? goal? belief))})

;; ---------------------------------------------------------------------------
;; Sensorless problems

(defn sensorless
  "A deterministic problem whose agent cannot see -- `:actions` the same
  in every state, `:result` (fn [state action]) -> the state it leads to
  -- as a search over belief states for `allgo.search`: from the set of
  every state it might start in, `:initial`, an action leads to the set
  of where each might go, and a belief is a goal when every state in it
  is. Breadth-first search on it gives the shortest sequence of actions
  that works from any of them; `actions` reads them off its path."
  [{:keys [initial actions result goal?]}]
  (let [step (fn [belief action] (set (map #(result % action) belief)))]
    {:start (set initial)
     :goal? (fn [belief] (every? goal? belief))
     :neighbors (fn [belief] (distinct (map #(step belief %) actions)))
     ::step step
     ::actions actions}))

(defn actions
  "The actions that take a sensorless problem along `path`, a sequence of
  belief states."
  [problem path]
  (let [step (::step problem)]
    (mapv (fn [[b b']] (first (filter #(= b' (step b %)) (::actions problem))))
          (partition 2 1 path))))
