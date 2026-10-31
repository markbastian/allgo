(ns allgo.search
  "State-space search: Russell and Norvig's chapter 3 (*Artificial
  Intelligence: A Modern Approach*, 4th ed.), as steps.

  Adapted from Mark Bastian's https://github.com/markbastian/planning,
  whose organizing idea is that breadth-first, depth-first, Dijkstra,
  greedy and A* are one algorithm. Each keeps a frontier of states to
  expand and a record of how it reached them; they differ only in which
  state comes off the frontier next. So there is one step function shape,
  one notion of a search state, and the algorithms are picked apart by the
  frontier they carry and the priority they assign.

  The same idea carries the searches that are usually written as loops
  within loops or as recursion. Every algorithm here is a step function
  from one search state to the next, run by `iterate`:

    best-first      breadth-first, depth-first, greedy, uniform-cost
                    (Dijkstra), A*, weighted A*, and beam search, which is
                    A* on a frontier that forgets all but its best few
    deepening       depth-limited and iterative-deepening search, and
                    IDA*: one step over an explicit stack of nodes, the
                    two deepening searches differing only in what they
                    bound -- depth, or g + h
    bidirectional   best-first from both ends at once (figure 3.14)
    recursive       recursive best-first search (figure 3.22), its
                    recursion reified by `call-seq` as a stack of frames,
                    each step a call or a return

  A search is a map:

    :start        the state to begin from
    :goal         the state to reach, or
    :goal?        (fn [state]) -> whether it is a goal, for problems with
                  many; the searches that need a goal state -- the
                  heuristics' second argument, bidirectional -- need
                  `:goal` too
    :neighbors    (fn [state]) -> the states reachable in one move
    :cost         (fn [from to]) -> the cost of that move (default 1)
    :heuristic    (fn [state goal]) -> estimated remaining cost (default 0)

  and, for the searches that use them, `:weight` (weighted A*, default 2),
  `:beam-width` (beam, default 16), `:limit` (depth-limited) and
  `:predecessors` (bidirectional, default `:neighbors`, right for any
  graph whose moves can be undone at the same cost).

  `:neighbors` is a function rather than a graph, so the state space can
  be implicit and unbounded -- board positions, puzzle configurations,
  hexes on a grid that is never materialized.

  Every algorithm is exposed twice. The `-seq` forms return the lazy seq
  of search states, one per step, which is what you want to watch a
  search work or to stop it early; the plain forms run to the goal and
  return the path. Every search state has `:current`, the state the step
  looked at, and `:solved?`; `goal-state` finds the first solved one and
  `path` reads the route out of it, whichever algorithm made it.

  A note on the priorities: Dijkstra orders by cost so far, greedy
  best-first by estimated cost remaining, A* by their sum and weighted A*
  by the sum with the estimate counted `:weight` times -- which is why
  they are all one step function, `best-first-seq`, given different
  priorities."
  #?(:clj (:import (clojure.lang PersistentQueue PersistentVector))))

;; ---------------------------------------------------------------------------
;; Frontiers
;;
;; The one thing the best-first algorithms vary. A queue makes the search
;; breadth first, a stack depth first, and a priority queue makes it one
;; of the cost-directed searches; a priority queue with a limit makes it a
;; beam.
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

(defn- trim
  "Drops the worst states until no more than `limit` are left."
  [by-priority priority-of limit]
  (if (and limit (> (count priority-of) (long limit)))
    (let [[worst states] (first (rseq by-priority))
          state (first states)
          others (disj states state)]
      (recur (if (seq others) (assoc by-priority worst others) (dissoc by-priority worst))
             (dissoc priority-of state)
             limit))
    [by-priority priority-of]))

(defrecord PriorityFrontier [by-priority priority-of limit]
  Frontier
  (put [_ entries]
    (let [{:keys [by-priority priority-of]}
          (reduce (fn [{:keys [by-priority priority-of] :as frontier} [state priority]]
                    (if-let [old (priority-of state)]
                      (if (<= old priority)
                        frontier
                        ;; A cheaper route to a state already queued: move
                        ;; it rather than queue it twice.
                        (let [pruned (disj (by-priority old) state)]
                          {:by-priority (-> (if (seq pruned)
                                              (assoc by-priority old pruned)
                                              (dissoc by-priority old))
                                            (update priority (fnil conj #{}) state))
                           :priority-of (assoc priority-of state priority)}))
                      {:by-priority (update by-priority priority (fnil conj #{}) state)
                       :priority-of (assoc priority-of state priority)}))
                  {:by-priority by-priority :priority-of priority-of}
                  entries)
          [by-priority priority-of] (trim by-priority priority-of limit)]
      (PriorityFrontier. by-priority priority-of limit)))
  (take-next [_]
    (when-let [[priority states] (first by-priority)]
      (let [state  (first states)
            others (disj states state)]
        [[state priority]
         (PriorityFrontier.
          (if (seq others) (assoc by-priority priority others) (dissoc by-priority priority))
          (dissoc priority-of state)
          limit)]))))

(defn priority-frontier
  "An empty priority queue; with `limit`, one that keeps only that many of
  its lowest-priority states."
  ([] (priority-frontier nil))
  ([limit] (->PriorityFrontier (sorted-map) {} limit)))

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

(defn- goal?
  "Whether `state` is a goal of `search`: its `:goal?`, or equal to its
  `:goal`."
  [{:keys [goal] :as search} state]
  (if-let [test (:goal? search)] (test state) (= state goal)))

(defn- found
  "The first of `states` that is a goal, or nil."
  [search states]
  (first (filter #(goal? search %) states)))

(defn- defaults [search]
  (-> search
      (update :cost (fn [c] (or c (constantly 1))))
      (update :heuristic (fn [h] (or h (constantly 0))))))

;; ---------------------------------------------------------------------------
;; Best-first steps

(defn- expand
  "Records that `state` was reached from `from`, for each state not seen."
  [search from states]
  (update search :came-from into (zipmap states (repeat from))))

(defn- unvisited [{:keys [neighbors came-from]} state]
  (remove #(contains? came-from %) (neighbors state)))

(defn- cheaper-neighbors
  "The neighbors of `state` this move improves on, with their new costs.

  Unlike the uninformed searches this does not skip states already seen: a
  state can be reached again more cheaply, and refusing to revisit it is
  what turns Dijkstra back into breadth-first search."
  [{:keys [neighbors cost costs]} state]
  (for [neighbor (neighbors state)
        :let  [new-cost (+ (costs state) (cost state neighbor))]
        :when (< new-cost (get costs neighbor ##Inf))]
    [neighbor new-cost]))

(defn- uninformed-step
  "Breadth-first or depth-first, depending only on the frontier passed in.
  Both test for the goal as it is generated -- the book's early goal
  test, which costs breadth-first nothing, since every move costs the
  same."
  [{:keys [frontier] :as search}]
  (when-let [[[state] frontier] (take-next frontier)]
    (let [discovered (unvisited search state)
          hit        (found search (cons state discovered))]
      (-> search
          (assoc :frontier (put frontier (map (fn [s] [s 0]) discovered))
                 :current state
                 :solved? (some? hit)
                 :found hit)
          (expand state discovered)))))

(defn- greedy-step
  "Best-first on the heuristic alone: always chase whatever looks closest.
  Fast and frequently wrong, since nothing accounts for distance covered."
  [{:keys [frontier heuristic goal] :as search}]
  (when-let [[[state] frontier] (take-next frontier)]
    (let [discovered (unvisited search state)
          hit        (found search (cons state discovered))]
      (-> search
          (assoc :frontier (put frontier (map (fn [s] [s (heuristic s goal)]) discovered))
                 :current state
                 :solved? (some? hit)
                 :found hit)
          (expand state discovered)))))

(defn- best-first-step
  "The book's BEST-FIRST-SEARCH (figure 3.7): expand the state of lowest
  `priority`, a function of the search, a state and the cost of reaching
  it, and queue every neighbor this reaches more cheaply than before.

  It is solved when the goal comes off the frontier, not when it is first
  reached: a costlier route can reach the goal first, and only once it is
  the cheapest thing left is nothing cheaper still to come."
  [{:keys [frontier priority] :as search}]
  (when-let [[[state] frontier] (take-next frontier)]
    (let [improved (cheaper-neighbors search state)]
      (-> search
          (assoc :frontier (put frontier (map (fn [[s c]] [s (priority search s c)]) improved))
                 :current state
                 :solved? (goal? search state)
                 :found (when (goal? search state) state))
          (update :costs into improved)
          (expand state (map first improved))))))

;; ---------------------------------------------------------------------------
;; Running a best-first search

(defn- initialize [{:keys [start frontier] :as search}]
  (-> (defaults search)
      (assoc :came-from {start nil}
             :costs     {start 0}
             :current   nil
             :solved?   false
             :frontier  (put frontier [[start 0]]))))

(defn- run
  "The seq of search states from `init`: each the step of the one before,
  until a step finds nothing left to do and returns nil. A solved state
  that a step is not asked past ends nothing; `goal-state` stops there."
  [step init]
  (take-while some? (iterate (fn [s] (when s (step s))) init)))

(defn- searcher [step frontier]
  (fn [search]
    (run step (initialize (assoc search :frontier (frontier search))))))

(def breadth-first-seq
  "Every search state of a breadth-first search, one per expansion."
  (searcher uninformed-step (constantly empty-queue)))

(def depth-first-seq
  "Every search state of a depth-first search.

  Requires a finite state space. See `depth-first`."
  (searcher uninformed-step (constantly [])))

(def greedy-seq
  "Every search state of a greedy best-first search."
  (searcher greedy-step (fn [_] (priority-frontier))))

(defn best-first-seq
  "Every search state of a best-first search ordered by `:priority`, a
  function of the search, a state and the cost of reaching it. The
  searches below are this, with a priority each."
  [search]
  ((searcher best-first-step (fn [_] (priority-frontier))) search))

(defn dijkstra-seq
  "Every search state of a uniform-cost (Dijkstra) search."
  [search]
  (best-first-seq (assoc search :priority (fn [_ _ cost] cost))))

(defn- a-star-priority [{:keys [heuristic goal weight]} state cost]
  (+ cost (* (or weight 1) (heuristic state goal))))

(defn a-star-seq
  "Every search state of an A* search."
  [search]
  (best-first-seq (-> search (dissoc :weight) (assoc :priority a-star-priority))))

(defn weighted-a-star-seq
  "Every search state of a weighted A* search: A* with the heuristic
  counted `:weight` times (default 2). It goes more directly for the goal
  and expands fewer states, for a path that may cost up to `:weight` times
  the best."
  [search]
  (best-first-seq (-> search (update :weight #(or % 2)) (assoc :priority a-star-priority))))

(defn beam-seq
  "Every search state of a beam search: A* whose frontier keeps only its
  `:beam-width` best states (default 16) and forgets the rest. It runs in
  bounded memory, and can lose the only way to the goal."
  [search]
  ((searcher best-first-step (fn [s] (priority-frontier (or (:beam-width s) 16))))
   (-> search (dissoc :weight) (assoc :priority a-star-priority))))

;; ---------------------------------------------------------------------------
;; Nodes, for the searches that keep paths rather than a table of states

(defn- node [state parent g]
  {:state state :parent parent :g g :depth (if parent (inc (long (:depth parent))) 0)})

(defn- node-path [n]
  (vec (reverse (map :state (take-while some? (iterate :parent n))))))

(defn- on-path? [n state]
  (some #(= state (:state %)) (take-while some? (iterate :parent n))))

(defn- children
  "The nodes one move from `n`, leaving out any that would go back to a
  state already on its own path."
  [{:keys [neighbors cost]} n]
  (for [s (neighbors (:state n))
        :when (not (on-path? n s))]
    (node s n (+ (:g n) (cost (:state n) s)))))

;; ---------------------------------------------------------------------------
;; Deepening: depth-limited, iterative deepening and IDA*
;;
;; One step over an explicit stack of nodes. Each pops a node: a node past
;; the bound is cut off, and the smallest amount by which any was past is
;; remembered; the goal within it ends the search; any other node is
;; replaced by its children. When the stack runs dry with something cut off, a
;; deepening search starts again from the root with the bound raised to
;; that smallest amount. Bound by depth, that is iterative deepening
;; (figures 3.12 and 3.14) -- raised by one each time; bound by g + h, it
;; is IDA*, raised to the next contour.

(defn- deepening-step
  [{:keys [stack bound limit next-limit deepen? root] :as search}]
  (cond
    (:solved? search) nil
    (:failed search) nil

    (empty? stack)
    (if (and deepen? (< next-limit ##Inf))
      (assoc search :stack [root] :limit next-limit :next-limit ##Inf
             :iteration (inc (long (:iteration search))) :current nil)
      ;; The book's two failures: cutoff, a solution may lie deeper, and
      ;; failure, there is none.
      (assoc search :failed (if (< next-limit ##Inf) :cutoff :failure)))

    :else
    (let [n (peek stack)
          stack (pop stack)
          s (assoc search :stack stack :current (:state n) :node n)
          b (bound search n)]
      ;; The bound before the goal: a node past it is not in this pass's
      ;; search at all, and a goal past it may not be the cheapest.
      (cond
        (> b limit) (assoc s :next-limit (min next-limit b))
        (goal? search (:state n)) (assoc s :solved? true :solution n)
        ;; Reversed, so that the first neighbor is the first popped.
        :else (assoc s :stack (into stack (reverse (children search n))))))))

(defn- deepening-seq [search bound limit deepen?]
  (let [{:keys [start] :as search} (defaults search)
        root (node start nil 0)]
    (run deepening-step
         (assoc search :root root :stack [root] :bound bound
                :limit limit :next-limit ##Inf :deepen? deepen?
                :iteration 0 :solved? false :current nil))))

(defn- depth-bound [_ n] (:depth n))

(defn- f-bound [{:keys [heuristic goal]} n] (+ (:g n) (heuristic (:state n) goal)))

(defn depth-limited-seq
  "Every search state of a depth-limited search to depth `:limit`. If it
  fails, its last state's `:failed` says how: `:cutoff` if something lay
  deeper, `:failure` if nothing did."
  [search]
  (deepening-seq search depth-bound (:limit search) false))

(defn iterative-deepening-seq
  "Every search state of an iterative-deepening search: depth-limited
  searches to depths 0, 1, 2, ... until one finds the goal, each state's
  `:limit` being the depth it is searching to. Breadth-first's shallowest
  solution in depth-first's memory, for a little repeated work."
  [search]
  (deepening-seq search depth-bound 0 true))

(defn ida-star-seq
  "Every search state of IDA*: iterative deepening bounded by g + h
  instead of depth, each pass's bound the smallest that cut something off
  the pass before."
  [search]
  (let [{:keys [start goal heuristic]} (defaults search)]
    (deepening-seq search f-bound (heuristic start goal) true)))

;; ---------------------------------------------------------------------------
;; Recursion as a seq

(defn call-seq
  "A recursive computation run a step at a time, as a seq of machine
  states `{:stack frames}`. `run-frame` takes the frame on top of the
  stack and returns `[:call child frame']`, to suspend the frame as
  `frame'` and run `child`, or `[:return value]`, to finish it. A frame
  that called another gets that call's value as its `:returned` when it
  next runs. The last state has an empty stack and the root's `:value`.

  With `context`, the computation has state beyond its frames -- a count,
  a table of results -- in each machine state's `:context`: `run-frame`
  then takes the frame and the context, and may return a new context as
  the fourth element of either answer, `[:call child frame' context']` or
  `[:return value nil context']`; nil leaves it as it was."
  ([run-frame root] (call-seq (fn [frame _] (run-frame frame)) root nil))
  ([run-frame root context]
   (run (fn [{:keys [stack context]}]
          (when (seq stack)
            (let [[op x frame' ctx] (run-frame (peek stack) context)
                  context (if (some? ctx) ctx context)
                  stack (pop stack)]
              (case op
                :call {:stack (conj stack frame' x) :event :call :context context}
                :return (if (seq stack)
                          {:stack (conj (pop stack) (assoc (peek stack) :returned x))
                           :event :return :context context}
                          {:stack [] :value x :event :return :context context})))))
        {:stack [root] :event :call :context context})))

;; ---------------------------------------------------------------------------
;; Recursive best-first search

(defn- rbfs-frame
  "One frame of RBFS (figure 3.22): a node, the f-cost it may not exceed,
  and once expanded its successors, each with the f it was last known to
  cost -- raised to what a search below it found when it gave up."
  [{:keys [goal heuristic] :as search}]
  (letfn [(decide [{:keys [succ f-limit] :as frame}]
            (let [order (sort-by (fn [i] (:f (succ i))) (range (count succ)))
                  best (first order)
                  best-f (:f (succ best))]
              (if (> best-f f-limit)
                [:return {:f best-f}]
                (let [alternative (if-let [i (second order)] (:f (succ i)) ##Inf)]
                  [:call {:node (succ best) :f-limit (min f-limit alternative)}
                   (assoc frame :chosen best)]))))]
    (fn [{:keys [node succ returned chosen] :as frame}]
      (cond
        (and (nil? succ) (goal? search (:state node))) [:return {:solution node}]

        (nil? succ)
        (let [kids (mapv (fn [c] (assoc c :f (max (+ (:g c) (heuristic (:state c) goal)) (:f node))))
                         (children search node))]
          (if (empty? kids)
            [:return {:f ##Inf}]
            (decide (assoc frame :succ kids))))

        (:solution returned) [:return returned]

        :else (decide (-> frame
                          (assoc-in [:succ chosen :f] (:f returned))
                          (dissoc :returned)))))))

(defn recursive-best-first-seq
  "Every machine state of a recursive best-first search: best-first in
  linear memory, keeping only the path it is on and, beside it, the best
  alternative, to which it backs up when the path costs more. Each state
  has the `:stack` of frames, and `:current`, the state on top."
  [search]
  (let [{:keys [start goal heuristic] :as search} (defaults search)
        root (assoc (node start nil 0) :f (heuristic start goal))]
    (map (fn [{:keys [stack value] :as m}]
           (assoc m
                  :current (:state (:node (peek stack)))
                  :solved? (some? (:solution value))
                  :solution (:solution value)))
         (call-seq (rbfs-frame search) {:node root :f-limit ##Inf}))))

;; ---------------------------------------------------------------------------
;; Bidirectional best-first search

(defn- top-priority [frontier]
  (if-let [[[_ p]] (take-next frontier)] p ##Inf))

(defn- bidirectional-step
  "Expands one state from whichever end has the cheaper state on top,
  noting each state reached from both ends as a possible meeting place,
  and is solved when no unexpanded state could lead to a cheaper meeting
  than the best found -- figure 3.14's TERMINATED."
  [{:keys [fwd bwd best start goal heuristic plain?] :as search}]
  (let [ff (top-priority (:frontier fwd))
        fb (top-priority (:frontier bwd))
        ;; Lower bounds on any route not yet found: with a consistent
        ;; heuristic, either end's smallest f; with none, their g's sum.
        ;; An end with nothing left has found all it can, and bounds
        ;; nothing: any route still to find runs through the other.
        bound (cond
                (and (= ff ##Inf) (= fb ##Inf)) ##Inf
                (= ff ##Inf) fb
                (= fb ##Inf) ff
                plain? (+ ff fb)
                :else (max ff fb))]
    (cond
      (:solved? search) nil
      (>= bound best) (when (< best ##Inf) (assoc search :solved? true :current nil))
      :else
      (let [dir (if (<= ff fb) :fwd :bwd)
            other (if (= dir :fwd) :bwd :fwd)
            {:keys [frontier costs came-from step-to cost-of]} (search dir)
            [[state] frontier] (take-next frontier)
            target (if (= dir :fwd) goal start)
            improved (for [n (step-to state)
                           :let [c (+ (costs state) (cost-of state n))]
                           :when (< c (get costs n ##Inf))]
                       [n c])
            side (assoc (search dir)
                        :frontier (put frontier (map (fn [[n c]] [n (+ c (heuristic n target))]) improved))
                        :costs (into costs improved)
                        :came-from (into came-from (map (fn [[n]] [n state]) improved)))
            other-costs (:costs (search other))
            meeting (reduce (fn [[b m] [n c]]
                              (let [total (+ c (get other-costs n ##Inf))]
                                (if (< total b) [total n] [b m])))
                            [best (:meet search)]
                            improved)]
        (assoc search dir side :best (first meeting) :meet (second meeting)
               :current state :direction dir)))))

(defn bidirectional-seq
  "Every search state of a bidirectional best-first search: uniform-cost
  from both ends, or with a heuristic, A* from both -- forward toward the
  goal along `:neighbors`, backward toward the start along
  `:predecessors`, the heuristic measuring to whichever end it heads for.
  Optimal for a consistent heuristic. Each state has `:direction`, the
  end it last expanded from."
  [search]
  (let [plain? (nil? (:heuristic search))
        {:keys [start goal cost] :as search} (defaults search)
        pred (or (:predecessors search) (:neighbors search))
        side (fn [s step-to cost-of]
               {:frontier (put (priority-frontier) [[s 0]]) :costs {s 0} :came-from {s nil}
                :step-to step-to :cost-of cost-of})]
    (run bidirectional-step
         (assoc search
                :plain? plain?
                :fwd (side start (:neighbors search) cost)
                ;; Backward, a move from n to s costs what s to n would.
                :bwd (side goal pred (fn [s n] (cost n s)))
                :best (if (= start goal) 0 ##Inf)
                :meet (when (= start goal) start)
                :solved? false
                :current nil))))

;; ---------------------------------------------------------------------------
;; Answers

(defn goal-state
  "The first search state that has reached the goal, or nil."
  [states]
  (first (filter :solved? states)))

(defn path
  "The route from start to goal recorded in a search state, as a vector of
  states, or nil when the goal was never reached -- read from whatever
  the algorithm keeps: the node it stopped on, the two halves of a
  bidirectional search, or the came-from table, walked back from the goal
  to the start, which maps to nil."
  [{:keys [goal came-from solution meet solved? fwd bwd] :as search}]
  (let [goal (or (:found search) goal)]
    (cond
      (nil? search) nil
      solution (node-path solution)
      (and fwd solved?) (let [walk (fn [table s] (take-while some? (iterate table s)))]
                          (into (vec (reverse (walk (:came-from fwd) meet)))
                                (rest (walk (:came-from bwd) meet))))
      fwd nil
      (contains? came-from goal) (vec (reverse (take-while some? (iterate came-from goal)))))))

(defn path-cost
  "What `path` costs under `search`'s `:cost`."
  [search path]
  (let [cost (or (:cost search) (constantly 1))]
    (reduce + 0 (map (fn [[a b]] (cost a b)) (partition 2 1 path)))))

(defn- solver [search-seq] (comp path goal-state search-seq))

(def breadth-first
  "Fewest moves from start to goal, ignoring cost."
  (solver breadth-first-seq))

(def depth-first
  "*A* route from start to goal. Rarely the shortest; cheap to find.

  It needs a finite state space. It commits to a branch and follows it to
  the end before considering any alternative, so where the others fan out
  and meet the goal, this one walks away from it forever: on an unbounded
  grid it runs to an OutOfMemoryError rather than returning a bad answer.
  Bound `:neighbors` -- to a board, a depth limit, anything -- or use
  `breadth-first` or `iterative-deepening`."
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

(def weighted-a-star
  "A route at most `:weight` times the cheapest, found with less work."
  (solver weighted-a-star-seq))

(def beam
  "A route found in bounded memory, or nil if the beam lost it."
  (solver beam-seq))

(def depth-limited
  "A route of at most `:limit` moves, or nil."
  (solver depth-limited-seq))

(def iterative-deepening
  "Fewest moves from start to goal, in memory proportional to the depth."
  (solver iterative-deepening-seq))

(def ida-star
  "The cheapest route, in memory proportional to its length, for an
  admissible heuristic."
  (solver ida-star-seq))

(def recursive-best-first
  "The cheapest route, in memory proportional to its length, for an
  admissible heuristic."
  (solver recursive-best-first-seq))

(def bidirectional
  "The cheapest route, searched for from both ends."
  (solver bidirectional-seq))
