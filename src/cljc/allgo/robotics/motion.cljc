(ns allgo.robotics.motion
  "Motion planning: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 26, which teaches it in the text rather
  than in figures of pseudocode.

  A world is a map:

    :bounds     `[[x0 y0] [x1 y1]]`, the room
    :obstacles  polygons, each a vector of `[x y]` corners, counterclockwise
    :radius     the robot's: a disk, so that its configuration space is
                the plane less the obstacles grown by the radius

  The planners, from exact to sampled:

    visibility graph     the corners of the grown obstacles, joined where
                         they see each other, searched by A*: the shortest
                         path, which grazes the obstacles
    Voronoi graph        the points equidistant from two obstacles,
                         searched by A*: a path that keeps as far from
                         them as it can, and may be long
    cell decomposition   a grid whose cells are wholly free, searched by A*:
                         complete only down to the cell size
    PRM                  a probabilistic roadmap: free samples joined to
                         their nearest neighbors, a step a sample, then
                         searched for each query
    RRT, RRT*            a tree grown toward random samples, a step a
                         sample; RRT* also chooses each new node's parent
                         and rewires its neighbors, so the path it holds
                         shortens toward the optimum as it grows
    potential field      descent on an attraction to the goal and a
                         repulsion from obstacles: fast, and stuck in any
                         local minimum

  and `shortcut`, which straightens a sampled path. The two-link arm
  shows what configuration space is for: `arm-cspace` marks which pairs
  of joint angles collide, and planning is a grid search there.

  Paths are vectors of points. The A* searches are `allgo.search`'s."
  (:require [allgo.geometry.delaunay :as delaunay]
            [allgo.geometry.polygon :as polygon]
            [allgo.geometry.quaternion :as q]
            [allgo.random :as random]
            [allgo.robotics.kinematics.chain :as chain]
            [allgo.search :as search]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Collision checking

(defn- edges-of [obstacles] (mapcat polygon/edges obstacles))

(defn clearance
  "How far `p` is from the nearest obstacle edge: `##Inf` with none."
  ^double [{:keys [obstacles]} p]
  (reduce (fn [^double m [a b]] (min m (polygon/segment-distance p a b))) ##Inf (edges-of obstacles)))

(defn- in-bounds? [{:keys [bounds radius]} [x y]]
  (let [[[x0 y0] [x1 y1]] bounds r (double (or radius 0.0))]
    (and (<= (+ x0 r) x (- x1 r)) (<= (+ y0 r) y (- y1 r)))))

(defn free?
  "Whether the robot fits at `p`: inside the room and every obstacle at
  least its radius away."
  [{:keys [obstacles radius] :as world} p]
  (and (in-bounds? world p)
       (not-any? #(polygon/contains-point? % p) obstacles)
       (>= (clearance world p) (double (or radius 0.0)))))

(defn segment-free?
  "Whether the robot can slide straight from `a` to `b`: both ends free
  and the segment at least its radius from every obstacle edge. With
  both ends outside every obstacle and no edge crossed, the segment is
  outside them all."
  [{:keys [obstacles radius] :as world} a b]
  (let [r (double (or radius 0.0))]
    (and (free? world a) (free? world b)
         (every? (fn [[c d]] (> (polygon/segment-segment-distance a b c d) (if (pos? r) r 0.0)))
                 (edges-of obstacles)))))

(defn path-free?
  "Whether every leg of `path` is."
  [world path]
  (every? (fn [[a b]] (segment-free? world a b)) (partition 2 1 path)))

(defn path-length
  "The length of `path`."
  ^double [path]
  (reduce + 0.0 (map (fn [[a b]] (polygon/distance a b)) (partition 2 1 path))))

(defn- graph-path
  "The shortest path through `nodes` (points) joined as `neighbors` says,
  by A* with straight-line distance, from node index `s` to `g`."
  [nodes neighbors s g]
  (when-let [p (search/a-star {:start s :goal g :neighbors neighbors
                               :cost (fn [i j] (polygon/distance (nodes i) (nodes j)))
                               :heuristic (fn [i j] (polygon/distance (nodes i) (nodes j)))})]
    (mapv nodes p)))

;; ---------------------------------------------------------------------------
;; Visibility graph

(defn grown
  "The obstacles grown by the robot's radius, corners mitered -- a little
  more than the configuration-space obstacle, whose corners are round --
  so that a point may go wherever the disk may."
  [{:keys [obstacles radius]}]
  (let [r (double (or radius 0.0))]
    (mapv #(if (pos? r) (polygon/inset % (- r)) %) obstacles)))

(defn- visible?
  "Whether the straight line `a`-`b` stays out of every polygon's
  interior: it crosses no edge properly, and no point along it lies inside
  a polygon shrunk by a hair, which lets it run along an edge or touch a
  corner."
  [polys shrunk a b]
  (and (not-any? (fn [[c d]]
                   (when-let [[t u] (polygon/segment-intersection a b c d)]
                     (and (< 1e-9 t (- 1.0 1e-9)) (< 1e-9 u (- 1.0 1e-9)))))
                 (edges-of polys))
       (not-any? (fn [k]
                   (let [p (polygon/lerp a b (/ k 16.0))]
                     (some #(polygon/contains-point? % p) shrunk)))
                 (range 1 16))))

(defn visibility-graph
  "The visibility graph of `start`, `goal` and the grown obstacles'
  corners that are inside the room: `{:nodes :edges}`, the edges index
  pairs."
  [world start goal]
  (let [polys (grown world)
        shrunk (keep #(polygon/inset % 1e-6) polys)
        corners (filter #(in-bounds? (assoc world :radius 0.0) %) (apply concat polys))
        nodes (vec (concat [start goal] corners))
        n (count nodes)]
    {:nodes nodes
     :edges (vec (for [i (range n) j (range (inc i) n)
                       :when (visible? polys shrunk (nodes i) (nodes j))]
                   [i j]))}))

(defn- adjacency [n edges]
  (reduce (fn [m [i j]] (-> m (update i conj j) (update j conj i))) (vec (repeat n [])) edges))

(defn visibility-path
  "The shortest path from `start` to `goal` for a robot of the world's
  radius, through the visibility graph -- shortest up to the mitered
  corners -- or nil."
  [world start goal]
  (let [{:keys [nodes edges]} (visibility-graph world start goal)
        adj (adjacency (count nodes) edges)]
    (graph-path nodes adj 0 1)))

;; ---------------------------------------------------------------------------
;; Cell decomposition

(defn cell-grid
  "The room cut into square cells of side `h`: `{:h :origin :cols :rows
  :free}`, `:free` a set of `[i j]` cells wholly free -- their centers
  clear of the obstacles by the radius plus half a diagonal. Cells that
  are partly free are left out, which makes the planner sound but
  incomplete: a passage narrower than a cell is not seen."
  [{:keys [bounds radius] :as world} h]
  (let [[[x0 y0] [x1 y1]] bounds
        cols (long (math/floor (/ (- x1 x0) h)))
        rows (long (math/floor (/ (- y1 y0) h)))
        margin (+ (double (or radius 0.0)) (* h (math/sqrt 0.5)))
        center (fn [[i j]] [(+ x0 (* h (+ i 0.5))) (+ y0 (* h (+ j 0.5)))])]
    {:h h :origin [x0 y0] :cols cols :rows rows
     :free (set (for [i (range cols) j (range rows)
                      :let [c (center [i j])]
                      :when (free? (assoc world :radius margin) c)]
                  [i j]))}))

(defn- cell-center [{:keys [h origin]} [i j]]
  (let [[x0 y0] origin] [(+ x0 (* h (+ i 0.5))) (+ y0 (* h (+ j 0.5)))]))

(defn- cell-of [{:keys [h origin]} [x y]]
  (let [[x0 y0] origin] [(long (math/floor (/ (- x x0) h))) (long (math/floor (/ (- y y0) h)))]))

(defn cell-seq
  "Every search state of A* through the free cells of `grid`, eight ways,
  from the cell holding `start` to the one holding `goal`."
  [grid start goal]
  (let [{:keys [free]} grid
        dist (fn [[i j] [k l]] (math/sqrt (+ (* (- i k) (- i k)) (* (- j l) (- j l)))))]
    (search/a-star-seq
     {:start (cell-of grid start) :goal (cell-of grid goal)
      :neighbors (fn [[i j]]
                   (for [di [-1 0 1] dj [-1 0 1]
                         :let [c [(+ i di) (+ j dj)]]
                         :when (and (not= 0 di dj) (free c)
                                    ;; no corner cutting past a blocked cell
                                    (free [(+ i di) j]) (free [i (+ j dj)]))]
                     c))
      :cost dist :heuristic dist})))

(defn cell-path
  "A path from `start` to `goal` through the free cells of a grid of side
  `h`, or nil: start, the cell centers, goal."
  [world start goal h]
  (let [grid (cell-grid world h)]
    (when-let [cells (search/path (search/goal-state (cell-seq grid start goal)))]
      (vec (concat [start] (map #(cell-center grid %) (rest (butlast cells))) [goal])))))

;; ---------------------------------------------------------------------------
;; Voronoi graph

(defn- boundaries
  "Every obstacle, and each wall of the room, as a separate owner of its
  edges: `[[owner [a b]] ...]`."
  [{:keys [bounds obstacles]}]
  (let [[[x0 y0] [x1 y1]] bounds
        walls [[[x0 y0] [x1 y0]] [[x1 y0] [x1 y1]] [[x1 y1] [x0 y1]] [[x0 y1] [x0 y0]]]]
    (concat (map-indexed (fn [i w] [[:wall i] w]) walls)
            (for [[i poly] (map-indexed vector obstacles) e (polygon/edges poly)] [i e]))))

(defn- boundary-samples
  "Points every `h` or less along each boundary edge: a map from each to
  its owner, a corner two owners share kept once."
  [world h]
  (into {} (for [[owner [a b]] (reverse (boundaries world))
                 :let [n (max 1 (long (math/ceil (/ (polygon/distance a b) h))))]
                 k (range n)]
             [(polygon/lerp a b (/ (double k) n)) owner])))

(defn voronoi-graph
  "The Voronoi graph of the free space: the points equidistant from two or
  more obstacles -- each wall of the room counting as one -- and at least
  the robot's radius from all of them. `{:nodes :edges}`, the edges index
  pairs.

  Made from the Voronoi diagram of points spaced `h` (default 0.2) along
  every boundary, by `allgo.geometry.delaunay`: of its edges, those that
  part two points on different obstacles. Points on one obstacle part
  only from each other, so those edges are dropped; what is left
  approximates the graph to within about `h`."
  ([world] (voronoi-graph world 0.2))
  ([world h]
   (let [owner (boundary-samples world h)
         tris (vec (delaunay/triangulate (keys owner)))
         ;; Each Delaunay edge, with the triangles either side of it.
         sides (reduce (fn [m [t {[a b c] :points}]]
                         (reduce (fn [m e] (update m (set e) (fnil conj []) t)) m [[a b] [b c] [c a]]))
                       {} (map-indexed vector tris))
         center #(get-in tris [% :circle :center])
         free-center? (memoize #(free? world (center %)))
         links (for [[e ts] sides
                     :when (= 2 (count ts))
                     :let [[a b] (vec e) [t u] ts]
                     :when (and (not= (owner a) (owner b)) (free-center? t) (free-center? u)
                                (segment-free? world (center t) (center u)))]
                 [t u])
         used (vec (distinct (mapcat identity links)))
         index (zipmap used (range))]
     {:nodes (mapv center used)
      :edges (vec (distinct (map (fn [[t u]] (vec (sort [(index t) (index u)]))) links)))})))

(defn voronoi-path
  "A path from `start` to `goal` along `graph` (default the world's
  Voronoi graph), as the chapter plans one: a straight move onto the
  graph, at the nearest point reachable that way, along it, and a
  straight move off it to the goal. Nil if either end cannot reach it or
  the graph does not join them. It keeps away from the obstacles rather
  than being short."
  ([world start goal] (voronoi-path world (voronoi-graph world) start goal))
  ([world {:keys [nodes edges]} start goal]
   (let [n (count nodes)
         onto (fn [p] (first (filter #(segment-free? world p (nodes %))
                                     (sort-by #(polygon/distance p (nodes %)) (range n)))))
         s (onto start) g (onto goal)]
     (when (and s g)
       (when-let [p (graph-path nodes (adjacency n edges) s g)]
         (vec (concat [start] p [goal])))))))

;; ---------------------------------------------------------------------------
;; Sampling

(defn- sample-free
  "A uniformly random free point."
  [{:keys [bounds] :as world} rng]
  (let [[[x0 y0] [x1 y1]] bounds]
    (loop []
      (let [p [(random/uniform rng x0 x1) (random/uniform rng y0 y1)]]
        (if (free? world p) p (recur))))))

(defn- nearest-k [nodes p k]
  (take k (sort-by #(polygon/distance (nodes %) p) (range (count nodes)))))

(defn prm-seq
  "A probabilistic roadmap growing a sample at a time: each state
  `{:nodes :edges}`, every new free sample joined to whichever of its `k`
  nearest earlier nodes it can reach straight (default 10)."
  ([world rng] (prm-seq world rng 10))
  ([world rng k]
   (iterate (fn [{:keys [nodes edges]}]
              (let [p (sample-free world rng)
                    i (count nodes)
                    joined (filter #(segment-free? world (nodes %) p) (nearest-k nodes p k))]
                {:nodes (conj nodes p) :edges (into edges (map (fn [j] [j i]) joined))}))
            {:nodes [] :edges []})))

(defn prm-path
  "A path from `start` to `goal` through `roadmap`: each joined to its
  `k` nearest nodes it can reach (default 10), then A*. Nil when the
  roadmap does not connect them."
  ([world roadmap start goal] (prm-path world roadmap start goal 10))
  ([world {:keys [nodes edges]} start goal k]
   (let [n (count nodes)
         nodes' (conj nodes start goal)
         links (fn [p i] (map (fn [j] [j i]) (filter #(segment-free? world (nodes %) p) (nearest-k nodes p k))))
         direct (when (segment-free? world start goal) [[n (inc n)]])
         edges' (concat edges (links start n) (links goal (inc n)) direct)]
     (graph-path nodes' (adjacency (count nodes') edges') n (inc n)))))

(defn- steer [from to step]
  (let [d (polygon/distance from to)]
    (if (<= d step) to (polygon/lerp from to (/ step d)))))

(defn- subtree [parent i]
  (let [children (reduce (fn [m [c p]] (if p (update m p conj c) m)) {} (map-indexed vector parent))]
    (loop [todo [i] out []]
      (if-let [[x & more] (seq todo)]
        (recur (concat more (children x)) (conj out x))
        out))))

(defn rrt-seq
  "A rapidly-exploring random tree, a sample a step: each state `{:nodes
  :parent :cost :goal}`, `:goal` the index of the node at the goal once
  the tree has reached it. Each step draws a point -- the goal itself with
  probability `:goal-bias` (default 0.05) -- finds the nearest node,
  steps toward the point by at most `:step` (default 0.5), and adds the
  new point if the way there is free; reaching within `:step` of the goal
  with a free line to it adds the goal.

  With `:star? true` it is RRT*: a new node takes as parent whichever node
  within `:radius` (default 4 x step x (log n / n)^(1/2), shrinking as the
  tree fills in, never below the step) gives it the shortest path, and
  those neighbors are rewired through it when that is shorter for them --
  so the path to the goal keeps shortening, toward the optimum."
  [world start goal rng & {:keys [step goal-bias star? radius]
                           :or {step 0.5 goal-bias 0.05}}]
  (let [dist polygon/distance]
    (iterate
     (fn [{:keys [nodes cost] :as tree}]
       (let [target (if (< (double (rng)) goal-bias) goal (sample-free world rng))
             near-i (apply min-key #(dist (nodes %) target) (range (count nodes)))
             p (steer (nodes near-i) target step)
             new? (and (not= p (nodes near-i)) (segment-free? world (nodes near-i) p))]
         (if-not new?
           tree
           (let [n (count nodes)
                 r (or radius (max step (* 4.0 step (math/sqrt (/ (math/log (inc n)) (inc n))))))
                 near (if star?
                        (filter #(and (<= (dist (nodes %) p) r) (segment-free? world (nodes %) p)) (range n))
                        [near-i])
                 best (apply min-key #(+ (cost %) (dist (nodes %) p)) (cons near-i near))
                 tree (-> tree
                          (update :nodes conj p)
                          (update :parent conj best)
                          (update :cost conj (+ (cost best) (dist (nodes best) p))))
                 ;; RRT*: rewire the neighbors through the new node when it
                 ;; shortens their paths, and their subtrees with them.
                 tree (if-not star?
                        tree
                        (reduce (fn [{:keys [cost parent] :as t} j]
                                  (let [c (+ (cost n) (dist p (nodes j)))]
                                    (if (< c (- (cost j) 1e-12))
                                      (let [delta (- c (cost j))]
                                        (-> t
                                            (assoc-in [:parent j] n)
                                            (update :cost (fn [cs] (reduce #(update %1 %2 + delta) cs (subtree (assoc parent j n) j))))))
                                      t)))
                                tree (remove #{best} near)))
                 reached? (and (nil? (:goal tree)) (<= (dist p goal) step)
                               (or (= p goal) (segment-free? world p goal)))]
             (cond
               (= p goal) (assoc tree :goal n)
               reached? (-> tree
                            (update :nodes conj goal)
                            (update :parent conj n)
                            (update :cost conj (+ ((:cost tree) n) (dist p goal)))
                            (assoc :goal (inc n)))
               :else tree)))))
     {:nodes [start] :parent [nil] :cost [0.0] :goal nil})))

(defn tree-path
  "The path from the root to the tree's goal node, or nil."
  [{:keys [nodes parent goal]}]
  (when goal
    (vec (reverse (map nodes (take-while some? (iterate parent goal)))))))

(defn shortcut-seq
  "`path` straightened a try at a time: each step picks two points on it
  at random and, if the straight line between them is free, cuts out what
  lies between."
  [world path rng]
  (iterate (fn [path]
             (let [n (count path)]
               (if (< n 3)
                 path
                 (let [i (random/below rng (- n 2))
                       j (+ i 2 (random/below rng (- n i 2)))]
                   (if (segment-free? world (path i) (path j))
                     (into (subvec path 0 (inc i)) (subvec path j))
                     path)))))
           path))

;; ---------------------------------------------------------------------------
;; Potential fields

(defn- nearest-obstacle-point [{:keys [obstacles]} p]
  (reduce (fn [[bd bp] [a b]]
            (let [q (polygon/lerp a b (polygon/closest-on-segment p a b))
                  d (polygon/distance p q)]
              (if (< d bd) [d q] [bd bp])))
          [##Inf nil] (edges-of obstacles)))

(defn potential-gradient
  "The gradient at `p` of the potential: an attraction of strength
  `:attract` to `goal`, U = attract |p - goal|^2 / 2, and a repulsion of
  strength `:repulse` from the nearest obstacle within `:influence`, U =
  repulse (1/d - 1/influence)^2 / 2, d its distance less the radius."
  [{:keys [radius] :as world} p goal {:keys [attract repulse influence]
                                      :or {attract 1.0 repulse 1.0 influence 1.0}}]
  (let [att (polygon/scale (polygon/sub p goal) attract)
        [d q] (nearest-obstacle-point world p)
        d (- (double d) (double (or radius 0.0)))]
    (if (and q (< 1e-9 d influence))
      (let [away (polygon/unit (polygon/sub p q))
            mag (* repulse (- (/ 1.0 d) (/ 1.0 influence)) (/ 1.0 (* d d)))]
        (polygon/sub att (polygon/scale away mag)))
      att)))

(defn potential-seq
  "Gradient descent on the potential from `start`, a step of length
  `:step` (default 0.05) at a time down the gradient: the positions. It
  arrives, or settles in a local minimum short of the goal."
  [world start goal & {:keys [step] :or {step 0.05} :as opts}]
  (iterate (fn [p]
             (if (<= (polygon/distance p goal) step)
               goal
               (let [g (potential-gradient world p goal opts)]
                 (if (< (polygon/norm g) 1e-12) p (polygon/sub p (polygon/scale (polygon/unit g) step))))))
           start))

;; ---------------------------------------------------------------------------
;; A two-link arm in configuration space

(defn arm
  "A planar two-link arm: its shoulder at `base`, links `lengths` long,
  as an `allgo.robotics.kinematics.chain` of two revolute joints about z -- in
  Craig's parameters, the second link a distance l1 along the first, the
  hand a tool l2 along the second."
  [base [l1 l2]]
  {:base base
   :chain (chain/chain [{:a 0.0} {:a l1}] {:rot q/identity-q :pos [l2 0.0 0.0]})})

(defn arm-links
  "The arm's two links as segments, `[[shoulder elbow] [elbow hand]]`, at
  joint angles `[t1 t2]`, the elbow angle relative to the upper arm."
  [{:keys [base chain]} angles]
  (let [at (fn [{[x y] :pos}] (polygon/add base [x y]))
        frames (chain/frames chain angles)
        shoulder (at (frames 1)) elbow (at (frames 2)) hand (at (chain/pose chain angles))]
    [[shoulder elbow] [elbow hand]]))

(defn arm-collides?
  "Whether the arm at `angles` touches an obstacle."
  [arm {:keys [obstacles]} angles]
  (boolean (some (fn [[a b]]
                   (or (some #(polygon/contains-point? % b) obstacles)
                       (some (fn [[c d]] (zero? (polygon/segment-segment-distance a b c d))) (edges-of obstacles))))
                 (arm-links arm angles))))

(defn- angle-of-cell ^double [n i] (- (* 2.0 math/PI (/ (+ i 0.5) n)) math/PI))

(defn arm-cspace
  "The arm's configuration space as an `n` by `n` grid over both joint
  angles, each in [-pi, pi): `(grid i j)` true where the arm at those
  angles collides."
  [arm world n]
  (mapv (fn [i] (mapv (fn [j] (arm-collides? arm world [(angle-of-cell n i) (angle-of-cell n j)])) (range n)))
        (range n)))

(defn arm-angles
  "The joint angles at the center of cell `[i j]` of an `n` grid."
  [n [i j]]
  [(angle-of-cell n i) (angle-of-cell n j)])

(defn arm-cell
  "The cell of an `n` grid holding `angles`."
  [n [t1 t2]]
  (let [c (fn [t] (mod (long (math/floor (* n (/ (+ t math/PI) (* 2.0 math/PI))))) n))]
    [(c t1) (c t2)]))

(defn arm-plan
  "A path of joint angles from `from` to `to` through the free cells of
  configuration space `grid`, by A* -- the angles wrap round, so the grid
  is a torus -- or nil."
  [grid from to]
  (let [n (count grid)
        free? (fn [[i j]] (not (get-in grid [i j])))
        wrap (fn [d] (let [d (abs d)] (min d (- n d))))
        dist (fn [[i j] [k l]] (math/sqrt (+ (math/pow (wrap (- i k)) 2) (math/pow (wrap (- j l)) 2))))
        s (arm-cell n from) g (arm-cell n to)]
    (when (and (free? s) (free? g))
      (when-let [cells (search/a-star {:start s :goal g
                                       :neighbors (fn [[i j]]
                                                    (for [di [-1 0 1] dj [-1 0 1]
                                                          :let [c [(mod (+ i di) n) (mod (+ j dj) n)]]
                                                          :when (and (not= 0 di dj) (free? c))]
                                                      c))
                                       :cost dist :heuristic dist})]
        (mapv #(arm-angles n %) cells)))))
