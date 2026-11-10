(ns allgo.robotics.coordination
  "Planning with time and with other robots: LaValle, *Planning
  Algorithms*, Cambridge University Press, 2006, sections 7.1 and 7.2.

    velocity tuning        the path is fixed, the obstacles move: choose
                           when to go and when to wait. A search over
                           (s, t), how far along the path and when, which
                           may only go forward in time
    coordination diagram   two robots on fixed paths: the square of their
                           two progresses, a cell blocked where they would
                           collide, searched for a way from corner to
                           corner -- each robot may wait, neither go back
    prioritized planning   many robots on a grid, planned one at a time in
                           order of priority, each in space and time
                           around the reservations of those before it.
                           Fast, and incomplete: a robot planned early
                           can block one planned later for good

  All three are A* (`allgo.search`'s) over a space with time in it, a
  step a tick."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn samples
  "Points along `path` at arc length 0, `ds`, 2 `ds`, ... and its end."
  [path ds]
  (let [legs (mapv (fn [[a b]] [a b (polygon/distance a b)]) (partition 2 1 path))
        total (reduce + (map peek legs))
        at (fn [s] (loop [[[a b l] & more] legs s s]
                     (if (or (empty? more) (<= s l)) (polygon/lerp a b (if (pos? l) (min 1.0 (/ s l)) 0.0)) (recur more (- s l)))))]
    (conj (mapv at (range 0.0 total ds)) (peek path))))

(defn velocity-tuning
  "When to move along `path` so as to miss the moving `obstacles`, each
  a function from time to its center, by at least `clearance`: the
  quickest schedule, `[[i t] ...]`, the robot at sample i (of
  `samples`, `ds` apart) at time t. Each tick of `dt` the robot moves
  forward up to `max-step` samples or waits. Nil if no schedule within
  `horizon` ticks."
  [path obstacles {:keys [ds dt max-step clearance horizon] :or {max-step 1 horizon 1000}}]
  (let [pts (samples path ds)
        n (dec (count pts))
        clear? (fn [i k] (let [p (pts i) t (* k dt)] (every? #(> (polygon/distance p (% t)) clearance) obstacles)))]
    (when (clear? 0 0)
      (some->> (search/a-star {:start [0 0]
                               :goal? (fn [[i]] (= i n))
                               :neighbors (fn [[i k]]
                                            (when (< k horizon)
                                              (for [d (range (inc max-step)) :let [j (min n (+ i d))] :when (clear? j (inc k))] [j (inc k)])))
                               :heuristic (fn [[i] _] (math/ceil (/ (- n i) (double max-step))))})
               (mapv (fn [[i k]] [i (* k dt)]))))))

(defn coordination-diagram
  "The coordination diagram of robots on `path-a` and `path-b`, sampled
  `ds` apart: the set of cells `[i j]` where robot a at its sample i and
  b at its j are within `clearance` -- the obstacle region of the
  square."
  [path-a path-b ds clearance]
  (let [a (samples path-a ds) b (samples path-b ds)]
    {:size [(count a) (count b)]
     :blocked (set (for [i (range (count a)) j (range (count b))
                         :when (<= (polygon/distance (a i) (b j)) clearance)]
                     [i j]))}))

(defn coordinate
  "The quickest schedule through the diagram from both at their starts
  to both at their ends: each tick one robot, or both, moves a sample
  on. `[[i j] ...]`, or nil if each blocks the other."
  [{[na nb] :size blocked :blocked}]
  (let [goal [(dec na) (dec nb)]]
    (search/a-star {:start [0 0] :goal goal
                    :neighbors (fn [[i j]]
                                 (for [[di dj] [[1 0] [0 1] [1 1]]
                                       :let [c [(+ i di) (+ j dj)]]
                                       :when (and (< (c 0) na) (< (c 1) nb) (not (blocked c)))]
                                   c))
                    :heuristic (fn [[i j] [gi gj]] (max (- gi i) (- gj j)))})))

(def ^:private moves [[0 0] [1 0] [-1 0] [0 1] [0 -1]])

(defn- plan-one
  "One robot's way from `start` to `goal` through `free?` cells and
  around `reserved`: `{:cells #{[cell t]} :swaps #{[from to t]} :parked
  {cell t}}`, those before it, the parked ones staying put forever.
  Each tick it moves a cell or waits; the cells it is in, by time."
  [free? reserved start goal horizon]
  (let [{:keys [cells swaps parked]} reserved
        taken? (fn [c t] (or (cells [c t]) (when-let [t0 (parked c)] (>= t t0))))
        later? (fn [c t] (some (fn [[c' t']] (and (= c c') (> t' t))) cells))
        manhattan (fn [[a b] [c d]] (+ (abs (- a c)) (abs (- b d))))]
    (some->> (search/a-star {:start [start 0]
                             :goal? (fn [[c t]] (and (= c goal) (not (later? c t))))
                             :neighbors (fn [[[i j :as c] t]]
                                          (when (< t horizon)
                                            (for [[di dj] moves
                                                  :let [n [(+ i di) (+ j dj)] t' (inc t)]
                                                  :when (and (free? n) (not (taken? n t')) (not (swaps [n c t'])))]
                                              [n t'])))
                             :heuristic (fn [[c] _] (manhattan c goal))})
             (mapv first))))

(defn- reserve [reserved path]
  (-> reserved
      (update :cells into (map-indexed (fn [t c] [c t]) path))
      (update :swaps into (map-indexed (fn [t [a b]] [a b (inc t)]) (partition 2 1 path)))
      (assoc-in [:parked (peek path)] (dec (count path)))))

(defn- prioritized-step
  "The next robot planned around the reservations so far, and its path
  reserved."
  [free? horizon {:keys [paths reserved]} [start goal]]
  (let [p (plan-one free? reserved start goal horizon)]
    {:paths (conj paths p) :reserved (if p (reserve reserved p) reserved)}))

(defn prioritized-seq
  "Prioritized planning, a robot a step: each robot of `robots`, `[start
  goal]` cells in order of priority, planned on the grid of `free?`
  cells around the robots before it. Each state `{:paths :reserved}`,
  the paths the cells by tick (nil for a robot that found none)."
  ([free? robots] (prioritized-seq free? robots 200))
  ([free? robots horizon]
   (reductions (partial prioritized-step free? horizon)
               {:paths [] :reserved {:cells #{} :swaps #{} :parked {}}}
               robots)))

(defn prioritized
  "Every robot's path by prioritized planning."
  ([free? robots] (:paths (last (prioritized-seq free? robots))))
  ([free? robots horizon] (:paths (last (prioritized-seq free? robots horizon)))))

(defn conflicts
  "Where `paths` collide: `[a b t]` for robots a and b in one cell at
  tick t, or swapping cells between t - 1 and t; a robot stays at its
  goal once there. Robots without a path are left out."
  [paths]
  (let [at (fn [p t] (p (min t (dec (count p)))))
        T (reduce max 0 (map count paths))]
    (vec (for [a (range (count paths)) b (range (inc a) (count paths)) t (range T)
               :let [pa (paths a) pb (paths b)]
               :when (and pa pb)
               :when (or (= (at pa t) (at pb t))
                         (and (pos? t) (= (at pa t) (at pb (dec t))) (= (at pb t) (at pa (dec t)))))]
           [a b t]))))
