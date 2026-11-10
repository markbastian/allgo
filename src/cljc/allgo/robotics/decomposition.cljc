(ns allgo.robotics.decomposition
  "Exact cell decomposition and coverage: LaValle, *Planning Algorithms*,
  Cambridge University Press, 2006, sections 6.2.2 and 7.6.

  The vertical decomposition cuts the free space of a polygonal world
  into trapezoids by a vertical line up and down from every corner, as
  far as the nearest edge. Each cell is convex, so a point robot crosses
  any cell in a straight line; two cells are neighbors where they share a
  piece of a vertical line, and the middle of that piece is a safe place
  to cross. A path is then a search over a handful of cells, complete and
  exact -- not sampled, not down to a grid's resolution.

  Coverage asks for a path that passes over every point, as a mower or a
  cleaning robot needs. The boustrophedon (\"ox-turning\") plan visits the
  cells depth first and sweeps each up and down in strips no wider than
  the tool -- the cells' convexity is what makes a strip of back-and-forth
  passes inside each one safe.

  Worlds are `allgo.robotics.motion`'s, the robot a point: grow the
  obstacles by its radius first (`motion/grown`). The obstacles are
  disjoint, inside the bounds. A cell is `{:x0 :x1 :lower :upper}`, the
  lower and upper sides each an edge `[p q]` (or the room's floor and
  ceiling) evaluated between x0 and x1."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- y-at
  "Where the line through edge `[p q]` crosses the vertical at `x`."
  ^double [[[px py] [qx qy]] x]
  (+ py (* (- qy py) (/ (- x px) (- qx px)))))

(defn- slab-cells
  "The free intervals of the slab between `a` and `b`, no corner strictly
  inside: each `{:x0 :x1 :lower :upper}` with the edges that bound it."
  [{[[x0 y0] [x1 y1]] :bounds obstacles :obstacles} a b]
  (let [xm (* 0.5 (+ a b))
        floor [[x0 y0] [x1 y0]]
        ceiling [[x0 y1] [x1 y1]]
        crossing (->> (mapcat polygon/edges obstacles)
                      (filter (fn [[[px] [qx]]] (< (min px qx) xm (max px qx))))
                      (sort-by #(y-at % xm)))
        sides (concat [floor] crossing [ceiling])]
    (vec (for [[lo hi] (partition 2 1 sides)
               :let [m [xm (* 0.5 (+ (y-at lo xm) (y-at hi xm)))]]
               ;; an obstacle standing on the floor leaves no room below it
               :when (> (- (y-at hi xm) (y-at lo xm)) 1e-9)
               :when (not-any? #(polygon/contains-point? % m) obstacles)]
           {:x0 a :x1 b :lower lo :upper hi}))))

(defn- merge-slab
  "The cells so far, with a new slab's joined to the cell left of each
  that has the same two sides -- no corner between them -- and the rest
  begun."
  [cells slab]
  (reduce (fn [cells {:keys [x0 x1 lower upper] :as c}]
            (if-let [i (first (keep-indexed (fn [i d] (when (and (= x0 (:x1 d)) (= lower (:lower d)) (= upper (:upper d))) i)) cells))]
              (assoc-in cells [i :x1] x1)
              (conj cells c)))
          cells slab))

(defn vertical-decomposition
  "The trapezoids of the free space of `world`."
  [{[[x0] [x1]] :bounds obstacles :obstacles :as world}]
  (let [xs (->> (mapcat identity obstacles) (map first) (concat [x0 x1]) (filter #(<= x0 % x1)) distinct sort)]
    (reduce merge-slab [] (map (fn [[a b]] (slab-cells world a b)) (partition 2 1 xs)))))

(defn interval
  "The cell's free interval `[lo hi]` on the vertical at `x`."
  [{:keys [lower upper]} x]
  [(y-at lower x) (y-at upper x)])

(defn corners
  "The cell's trapezoid, counterclockwise."
  [{:keys [x0 x1] :as c}]
  (let [[a0 b0] (interval c x0) [a1 b1] (interval c x1)]
    [[x0 a0] [x1 a1] [x1 b1] [x0 b0]]))

(defn- center [c] (polygon/centroid (corners c)))

(defn- door
  "The middle of the piece of vertical cells `c` and `d` share, c left of
  d, or nil if they share none."
  [c d]
  (when (= (:x1 c) (:x0 d))
    (let [x (:x1 c) [a b] (interval c x) [e f] (interval d x)
          lo (max a e) hi (min b f)]
      (when (> hi lo) [x (* 0.5 (+ lo hi))]))))

(defn adjacency
  "Each cell's neighbors by index, and the door between them: a map
  from `[i j]` to the point."
  [cells]
  (into {} (for [i (range (count cells)) j (range (count cells))
                 :let [p (or (door (cells i) (cells j)) (door (cells j) (cells i)))]
                 :when p]
             [[i j] p])))

(defn cell-of
  "The index of the cell that holds `p`, or nil."
  [cells [x y]]
  (first (keep-indexed (fn [i {:keys [x0 x1] :as c}]
                         (when (<= x0 x x1) (let [[lo hi] (interval c x)] (when (<= lo y hi) i))))
                       cells)))

(defn- cell-route
  "The cells from `i` to `j`, by A* over the cells' centers."
  [cells doors i j]
  (let [centers (mapv center cells)
        neighbors (reduce (fn [m [a b]] (update m a (fnil conj []) b)) {} (keys doors))
        d (fn [a b] (polygon/distance (centers a) (centers b)))]
    (search/a-star {:start i :goal j :neighbors #(neighbors % []) :cost d :heuristic d})))

(defn path
  "A path from `start` to `goal` through the decomposition: through the
  door of each cell boundary crossed and, between two doors, the cell's
  center -- two doors can lie on one side of a cell, and the straight
  way between them along its edge could graze an obstacle's corner. So
  the path is straight inside every cell and free. Nil if either is in
  no cell or no cells join them."
  ([world start goal] (let [cells (vertical-decomposition world)] (path cells (adjacency cells) start goal)))
  ([cells doors start goal]
   (when-let [route (some->> [(cell-of cells start) (cell-of cells goal)]
                             (#(when (every? some? %) %))
                             (apply cell-route cells doors))]
     (let [ds (map doors (partition 2 1 route))
           centers (map (comp center cells) (rest (butlast route)))]
       (into [start] (if (seq ds) (interleave ds (concat centers [goal])) [goal]))))))

(defn- sweep
  "Back-and-forth passes over cell `c` in strips no wider than `width`,
  kept `margin` off its floor and ceiling: the points, from the left."
  [c width margin]
  (let [{:keys [x0 x1]} c
        n (max 1 (long (math/ceil (/ (- x1 x0) width))))]
    (vec (mapcat (fn [k]
                   (let [x (+ x0 (* (+ k 0.5) (/ (- x1 x0) n)))
                         [lo hi] (interval c x)
                         pass [[x (+ lo margin)] [x (- hi margin)]]]
                     (if (odd? k) (reverse pass) pass)))
                 (range n)))))

(defn depth-first-order
  "The cells in the order a depth-first walk from cell `i` first
  reaches them."
  [cells doors i]
  (let [neighbors (reduce (fn [m [a b]] (update m a (fnil conj []) b)) {} (keys doors))]
    (loop [stack [i] seen #{} order []]
      (if-let [c (peek stack)]
        (if (seen c)
          (recur (pop stack) seen order)
          (recur (into (pop stack) (reverse (sort-by #(:x0 (cells %)) (neighbors c [])))) (conj seen c) (conj order c)))
        order))))

(defn boustrophedon
  "A coverage path from `start`: the cells in depth-first order, each
  swept in passes no farther apart than `width`, the way from one cell's
  last pass to the next's first through the doors. `{:order :path}`,
  the cells' indices and the points."
  ([world start width] (boustrophedon world start width 1e-6))
  ([world start width margin]
   (let [cells (vertical-decomposition world)
         doors (adjacency cells)
         order (depth-first-order cells doors (cell-of cells start))]
     {:order order
      :path (reduce (fn [p i]
                      (let [s (sweep (cells i) width margin)]
                        (into p (concat (rest (path cells doors (peek p) (first s))) (rest s)))))
                    [start] order)})))
