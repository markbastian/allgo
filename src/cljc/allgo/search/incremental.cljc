(ns allgo.search.incremental
  "Incremental search: replanning as the world turns out different from
  the map -- LaValle, *Planning Algorithms*, Cambridge University Press,
  2006, section 12.3.2, Stentz's D*; here in the form Koenig and Likhachev
  gave it, D* Lite (\"D* Lite\", AAAI 2002), which is simpler and does the
  same.

  A robot plans a path to the goal on the map it has, assuming what it has
  not seen is clear, and sets off. When it finds an edge costs more than
  it thought -- a wall where the map had floor -- it must replan; but most
  of the old plan is still right, so rather than search again from
  nothing, D* Lite repairs it. It searches backward, from the goal, so
  that the robot's moving does not invalidate the costs it has; it keeps
  for each state g, its cost-to-goal as last settled, and rhs, the cost
  its neighbors now imply; and only the states where the two disagree go
  on the queue, ordered by a key that adds a growing `km` for how far the
  robot has moved, rather than re-sorting everything.

  A problem is `{:start :goal :neighbors :cost :heuristic}`: `(neighbors
  s)` the states joined to s (the graph undirected), `(cost a b)` the
  cost of the edge as the map has it, `(heuristic a b)` an admissible
  estimate between two states. States must be comparable, as vectors are."
  (:require [allgo.search :as search]))

(def ^:private inf ##Inf)

(defn- g [st s] (get-in st [:g s] inf))
(defn- rhs [st s] (get-in st [:rhs s] inf))
(defn- edge [st a b] (get-in st [:changed [a b]] ((:cost st) a b)))

(defn- key-of [{:keys [start heuristic km] :as st} s]
  (let [m (min (g st s) (rhs st s))]
    [(+ m (heuristic start s) km) m]))

(defn- dequeue [st s]
  (if-let [k (get-in st [:keys s])]
    (-> st (update :open disj (conj k s)) (update :keys dissoc s))
    st))

(defn- enqueue [st s]
  (let [k (key-of st s)]
    (-> st (update :open conj (conj k s)) (assoc-in [:keys s] k))))

(defn- update-vertex [{:keys [goal neighbors] :as st} u]
  (let [st (if (= u goal)
             st
             (assoc-in st [:rhs u] (reduce min inf (map #(+ (edge st u %) (g st %)) (neighbors u)))))
        st (dequeue st u)]
    (if (not= (g st u) (rhs st u)) (enqueue st u) st)))

(defn- key< [[a1 a2] [b1 b2]] (or (< a1 b1) (and (= a1 b1) (< a2 b2))))

(defn- settle
  "Expands states until the start's cost is settled: the heart of D*
  Lite, a state taken from the queue when its key is stale re-queued,
  when its g is too high lowered to its rhs, and when too low raised to
  infinity, its neighbors updated either way. Counts the expansions in
  `:expanded`."
  [{:keys [neighbors start] :as st}]
  (loop [st st]
    (let [[k1 k2 u :as top] (first (:open st))]
      (if (and top (or (key< [k1 k2] (key-of st start)) (not= (rhs st start) (g st start))))
        (let [st (-> st (update :open disj top) (update :keys dissoc u) (update :expanded inc))
              k-new (key-of st u)]
          (recur (cond (key< [k1 k2] k-new) (enqueue st u)
                       (> (g st u) (rhs st u)) (reduce update-vertex (assoc-in st [:g u] (rhs st u)) (neighbors u))
                       :else (reduce update-vertex (assoc-in st [:g u] inf) (cons u (neighbors u))))))
        st))))

(defn plan
  "D* Lite's state for `problem`, its first plan made: `{:start :g :rhs
  :expanded ...}`."
  [{:keys [goal heuristic] :as problem}]
  (settle (merge problem {:g {} :rhs {goal 0.0} :km 0.0 :changed {} :expanded 0 :last (:start problem)
                          :open (sorted-set (conj [(heuristic (:start problem) goal) 0.0] goal))
                          :keys {goal [(heuristic (:start problem) goal) 0.0]}})))

(defn path-cost
  "The cost of the plan from the robot's position, as the map now has it."
  [st]
  (g st (:start st)))

(defn next-state
  "Where the plan goes next from the robot's position: the neighbor
  minimizing the edge's cost plus the neighbor's cost-to-goal; nil if
  every way is blocked."
  [{:keys [start neighbors] :as st}]
  (let [best (apply min-key #(+ (edge st start %) (g st %)) (neighbors start))]
    (when (< (+ (edge st start best) (g st best)) inf) best)))

(defn replan
  "The plan repaired after the robot, now at `position`, learns the edge
  costs `changes`, each `[a b cost]` (both directions): `km` grows by the
  heuristic distance it has moved since the last repair, each changed
  edge's ends updated, and the queue settled again."
  [{:keys [heuristic last] :as st} position changes]
  (if (empty? changes)
    (assoc st :start position)
    (let [st (-> st
                 (assoc :start position)
                 (update :km + (heuristic last position))
                 (assoc :last position))
          st (reduce (fn [st [a b c]] (-> st (assoc-in [:changed [a b]] c) (assoc-in [:changed [b a]] c))) st changes)
          st (reduce update-vertex st (distinct (mapcat (fn [[a b]] [a b]) changes)))]
      (settle st))))

(defn- news
  "The edges seen from `position` whose costs differ from the plan's."
  [st sense position]
  (filter (fn [[a b c]] (not= c (edge st a b))) (sense position)))

(defn- navigate-step
  "A move to the plan's next state and a repair for what is seen there;
  nil at the goal or when stuck."
  [sense st]
  (when (not= (:start st) (:goal st))
    (when-let [nxt (next-state st)]
      (replan st nxt (news st sense nxt)))))

(defn navigate-seq
  "The robot driving to the goal, a move a step: from each position it
  calls `(sense position)` for the edges whose true costs it can now see,
  `[[a b cost] ...]`, repairs the plan with those that differ from what it
  believed, and moves to the plan's next state. Each state D* Lite's,
  `:start` the robot's position, until it reaches the goal or is stuck."
  [problem sense]
  (let [first-state (plan problem)]
    (search/step-seq (partial navigate-step sense)
                     (replan first-state (:start problem) (news first-state sense (:start problem))))))
