(ns allgo.search.nondeterministic
  "Search where an action's outcome is not known, or the state is not:
  Russell and Norvig's sections 4.3 and 4.4.

    and-or-seq    AND-OR search (figure 4.11, AND-OR-GRAPH-SEARCH): a
                  plan that reaches a goal whatever the actions turn out
                  to do. Its OR-SEARCH and AND-SEARCH call each other; the
                  recursion is run a step at a time by
                  `allgo.search/call-seq`, so it is a seq of machine
                  states like every other search here
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

(defn- and-or-frame
  "OR-SEARCH and AND-SEARCH, as frames. An OR frame tries its state's
  actions in turn, each an AND frame over the states the action may lead
  to, and returns the first action every one of those can be planned
  from. An AND frame plans from each of its states in turn and fails if
  any fails."
  [{:keys [results goal?] :as problem}]
  (fn [{:keys [kind state path actions action returned states i plans] :as frame}]
    (case kind
      :or
      (cond
        (and (nil? action) (goal? state)) [:return []]
        (and (nil? action) (contains? path state)) [:return ::failure]
        (and returned (not= ::failure returned)) [:return [action returned]]
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
  that says 'try again until it works' is section 4.3.3's cyclic plan,
  not this."
  [problem]
  (search/call-seq (and-or-frame problem) (or-frame problem (:initial problem) #{})))

(defn and-or
  "A conditional plan that reaches a goal however the actions turn out,
  or nil if there is none without cycles."
  [problem]
  (let [v (:value (last (and-or-seq problem)))]
    (when (not= ::failure v) v)))

(defn execute
  "The states a plan passes through when the world picks each outcome
  with `choose` -- (fn [states]) -> one of them -- from `state`."
  [{:keys [results]} plan state choose]
  (loop [plan plan state state seen [state]]
    (if (empty? plan)
      seen
      (let [[action plans] plan
            nxt (choose (vec (results state action)))]
        (recur (plans nxt) nxt (conj seen nxt))))))

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
