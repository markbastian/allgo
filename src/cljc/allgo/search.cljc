(ns allgo.search
  "State-space search: breadth-first, depth-first, Dijkstra, greedy
  best-first and A*.

  Adapted from Mark Bastian's https://github.com/markbastian/planning,
  whose organising idea is that these five are one algorithm. Each keeps a
  frontier of states to expand and a record of how it reached them; they
  differ only in which state comes off the frontier next. So there is one
  step function shape, one notion of a search state, and the algorithms
  are picked apart by the frontier they carry and the priority they assign.

  A search is a map:

    :start       the state to begin from
    :goal        the state to reach
    :neighbours  (fn [state]) -> the states reachable in one move
    :cost        (fn [from to]) -> the cost of that move (default 1)
    :heuristic   (fn [state goal]) -> estimated remaining cost (default 0)

  `:neighbours` is a function rather than a graph, so the state space can
  be implicit and unbounded -- board positions, puzzle configurations,
  hexes on a grid that is never materialised.

  Every algorithm is exposed twice. The `-seq` forms return the lazy seq
  of search states, one per expansion, which is what you want to watch a
  search work or to stop it early; the plain forms run to the goal and
  return the path. Because the seq is lazy and each element is a complete
  search state, animating a search is just rendering its elements.

  A note on the two priorities: Dijkstra orders by cost so far, greedy
  best-first by estimated cost remaining, and A* by their sum -- which is
  why A* is written here as one line of arithmetic on top of Dijkstra."
  #?(:clj (:import (clojure.lang PersistentQueue PersistentVector))))

;; ---------------------------------------------------------------------------
;; Frontiers
;;
;; The one thing the algorithms vary. A queue makes the search breadth
;; first, a stack depth first, and a priority queue makes it one of the
;; three cost-directed searches.
;;
;; Every frontier holds `[state priority]` entries even where the priority
;; is ignored, so the step functions never have to know which kind they
;; hold. The priority queue is implemented here rather than pulled in:
;; org.clojure/data.priority-map is Clojure-only, and a search that cannot
;; run in the browser would be half a library.

(defprotocol Frontier
  (put [frontier entries] "Adds `[state priority]` entries.")
  (take-next [frontier]
    "`[[state priority] frontier']`, or nil when the frontier is empty."))

(defrecord PriorityFrontier [by-priority priority-of]
  Frontier
  (put [_ entries]
    (reduce (fn [{:keys [by-priority priority-of] :as frontier} [state priority]]
              (if-let [old (priority-of state)]
                (if (<= old priority)
                  frontier
                  ;; A cheaper route to a state already queued: move it
                  ;; rather than queue it twice.
                  (let [pruned (disj (by-priority old) state)]
                    (PriorityFrontier.
                     (-> (if (seq pruned)
                           (assoc by-priority old pruned)
                           (dissoc by-priority old))
                         (update priority (fnil conj #{}) state))
                     (assoc priority-of state priority))))
                (PriorityFrontier.
                 (update by-priority priority (fnil conj #{}) state)
                 (assoc priority-of state priority))))
            (PriorityFrontier. by-priority priority-of)
            entries))
  (take-next [_]
    (when-let [[priority states] (first by-priority)]
      (let [state  (first states)
            others (disj states state)]
        [[state priority]
         (PriorityFrontier.
          (if (seq others) (assoc by-priority priority others) (dissoc by-priority priority))
          (dissoc priority-of state))]))))

(defn priority-frontier [] (->PriorityFrontier (sorted-map) {}))

(extend-protocol Frontier
  #?(:clj PersistentQueue :cljs cljs.core/PersistentQueue)
  (put [frontier entries] (into frontier entries))
  (take-next [frontier]
    (when (seq frontier) [(peek frontier) (pop frontier)]))

  #?(:clj PersistentVector :cljs cljs.core/PersistentVector)
  (put [frontier entries] (into frontier entries))
  (take-next [frontier]
    (when (seq frontier) [(peek frontier) (pop frontier)])))

(def empty-queue #?(:clj PersistentQueue/EMPTY :cljs #queue []))

;; ---------------------------------------------------------------------------
;; Steps

(defn- expand
  "Records that `state` was reached from `from`, for each state not seen."
  [search from states]
  (update search :came-from into (zipmap states (repeat from))))

(defn- unvisited [{:keys [neighbours came-from]} state]
  (remove #(contains? came-from %) (neighbours state)))

(defn- cheaper-neighbours
  "The neighbours of `state` this move improves on, with their new costs.

  Unlike the uninformed searches this does not skip states already seen: a
  state can be reached again more cheaply, and refusing to revisit it is
  what turns Dijkstra back into breadth-first search."
  [{:keys [neighbours cost costs]} state]
  (for [neighbour (neighbours state)
        :let  [new-cost (+ (costs state) (cost state neighbour))]
        :when (< new-cost (get costs neighbour ##Inf))]
    [neighbour new-cost]))

(defn- uninformed-step
  "Breadth-first or depth-first, depending only on the frontier passed in."
  [{:keys [frontier] :as search}]
  (let [[[state] frontier] (take-next frontier)
        discovered         (unvisited search state)]
    (-> search
        (assoc :frontier (put frontier (map (fn [s] [s 0]) discovered)))
        (expand state discovered))))

(defn- greedy-step
  "Best-first on the heuristic alone: always chase whatever looks closest.
  Fast and frequently wrong, since nothing accounts for distance covered."
  [{:keys [frontier heuristic goal] :as search}]
  (let [[[state] frontier] (take-next frontier)
        discovered         (unvisited search state)]
    (-> search
        (assoc :frontier (put frontier (map (fn [s] [s (heuristic s goal)]) discovered)))
        (expand state discovered))))

(defn- cost-directed-step
  "Dijkstra when `priority` is the cost so far, A* when it also carries the
  heuristic."
  [priority]
  (fn [{:keys [frontier] :as search}]
    (let [[[state] frontier] (take-next frontier)
          improved           (cheaper-neighbours search state)]
      (-> search
          (assoc :frontier (put frontier (map (fn [[s c]] [s (priority search s c)]) improved)))
          (update :costs into improved)
          (expand state (map first improved))))))

(def ^:private dijkstra-step
  (cost-directed-step (fn [_ _ cost] cost)))

(def ^:private a-star-step
  (cost-directed-step (fn [{:keys [heuristic goal]} state cost]
                        (+ cost (heuristic state goal)))))

;; ---------------------------------------------------------------------------
;; Running a search

(defn- initialize [{:keys [start frontier] :as search}]
  (-> search
      (update :cost (fn [c] (or c (constantly 1))))
      (update :heuristic (fn [h] (or h (constantly 0))))
      (assoc :came-from {start nil}
             :costs     {start 0}
             :frontier  (put frontier [[start 0]]))))

(defn- searcher [step frontier]
  (fn [search]
    (->> (initialize (assoc search :frontier frontier))
         (iterate step)
         (take-while (comp some? :frontier))
         ;; The seq ends when the frontier runs dry. A state cannot be
         ;; reached after that, so stopping here loses nothing.
         (take-while (fn [{:keys [frontier]}] (some? (take-next frontier)))))))

(def breadth-first-seq
  "Every search state of a breadth-first search, one per expansion."
  (searcher uninformed-step empty-queue))

(def depth-first-seq
  "Every search state of a depth-first search."
  (searcher uninformed-step []))

(def greedy-seq
  "Every search state of a greedy best-first search."
  (searcher greedy-step (priority-frontier)))

(def dijkstra-seq
  "Every search state of a uniform-cost (Dijkstra) search."
  (searcher dijkstra-step (priority-frontier)))

(def a-star-seq
  "Every search state of an A* search."
  (searcher a-star-step (priority-frontier)))

(defn goal-state
  "The first search state that has reached the goal, or nil."
  [states]
  (first (filter (fn [{:keys [goal came-from]}] (contains? came-from goal)) states)))

(defn path
  "The route from start to goal recorded in a search state, as a vector of
  states, or nil when the goal was never reached.

  Walks the came-from chain back from the goal; the start maps to nil,
  which is what ends the walk."
  [{:keys [goal came-from] :as search}]
  (when (and search (contains? came-from goal))
    (vec (reverse (take-while some? (iterate came-from goal))))))

(defn- solver [search-seq] (comp path goal-state search-seq))

(def breadth-first
  "Fewest moves from start to goal, ignoring cost."
  (solver breadth-first-seq))

(def depth-first
  "*A* route from start to goal. Rarely the shortest; cheap to find."
  (solver depth-first-seq))

(def greedy
  "A route found by always heading toward the goal. Fast, not optimal."
  (solver greedy-seq))

(def dijkstra
  "The cheapest route, exploring outward by cost in every direction."
  (solver dijkstra-seq))

(def a-star
  "The cheapest route, exploring toward the goal first. Optimal whenever
  the heuristic never overestimates what remains."
  (solver a-star-seq))
