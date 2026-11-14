(ns allgo.robotics.clearance
  "The maximum-clearance roadmap of a polygonal world, built the way
  LaValle, *Planning Algorithms*, Cambridge University Press, 2006,
  section 6.2.3, describes it (and his Exercise 14 asks for): the
  generalized Voronoi diagram, every point of it equidistant from two
  points of the obstacles' boundary and nearer to them than to any other.

  The features are the corners and edges of the obstacles and the room.
  Every pair of them has a curve of points equidistant from the two,
  Figure 6.9: two edges, a line bisecting the angle between them (or the
  midline of two parallel ones); two corners, the perpendicular bisector;
  a corner and an edge, a parabola, the corner its focus and the edge's
  line its directrix. The roadmap is the parts of those curves where the
  two features are the nearest of all -- found by walking each curve and
  narrowing each end, where some third feature comes as near, by
  bisection: the curves' intersections. This is the naive construction
  the section gives, quadratic in the pairs and linear in the features
  for each test, not the n log n sweep it cites; the pairs on one convex
  obstacle are skipped, their curves lying inside it.

  Unlike `allgo.robotics.motion`'s Voronoi graph, approximated from
  points along the boundaries, this one is exact up to rounding, its
  curves straight or parabolic. A robot that is a disk of radius r keeps
  only the parts with clearance at least r -- no obstacle need be grown.

  A query point joins the roadmap by retraction, the section's other
  name for the method: it moves straight away from its nearest boundary
  point, its clearance growing, until another feature is as near -- a
  point of the roadmap.

  Worlds are `allgo.robotics.motion`'s: `:bounds`, the room, and
  `:obstacles`, disjoint counterclockwise polygons inside it."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- room [{[[x0 y0] [x1 y1]] :bounds}]
  [[x0 y0] [x1 y0] [x1 y1] [x0 y1]])

(defn- diagonal [{[a b] :bounds}] (polygon/distance a b))

(defn features
  "The corners and edges of the room and the obstacles, each
  `{:kind :id :owner}` with `:p` for a corner, `:a` `:b` and the ids of
  its `:ends` for an edge."
  [{:keys [obstacles] :as world}]
  (vec (for [[owner poly] (cons [:room (room world)] (map-indexed vector obstacles))
             :let [n (count poly)]
             [k [a b]] (map-indexed vector (polygon/edges poly))
             f [{:kind :corner :owner owner :id [owner :corner k] :p a}
                {:kind :edge :owner owner :id [owner :edge k] :a a :b b
                 :ends #{[owner :corner k] [owner :corner (mod (inc k) n)]}}]]
         f)))

(defn distance-to
  "How far `p` is from feature `f`: from an edge, along a perpendicular
  meeting it strictly between its ends, and infinitely far otherwise --
  then a corner is nearer."
  ^double [f p]
  (if (= :corner (:kind f))
    (polygon/distance p (:p f))
    (let [{:keys [a b]} f
          ab (polygon/sub b a)
          ap (polygon/sub p a)
          t (/ (polygon/dot ap ab) (polygon/dot ab ab))]
      (if (< 0.0 t 1.0) (/ (abs (polygon/cross ab ap)) (polygon/norm ab)) ##Inf))))

(defn free-point?
  "Whether `p` is in the room and in no obstacle."
  [{:keys [obstacles] :as world} p]
  (let [[[x0 y0] [x1 y1]] (:bounds world) [x y] p]
    (and (<= x0 x x1) (<= y0 y y1)
         (not-any? #(polygon/contains-point? % p) obstacles))))

;; ------------------------------------------------------------- bisectors

(defn- perp [[x y]] [(- y) x])

(defn- bisectors
  "The curves of points equidistant from features `f` and `g`, Figure
  6.9: lines `{:o :d}` and parabolas `{:a :u :s :h}` -- the point at
  parameter s on the directrix a + s u, raised along its normal by the
  height at which it is as far from the focus, (s - s_v)^2 + d^2 over
  2d."
  [f g]
  (let [[f g] (sort-by #(if (= :corner (:kind %)) 0 1) [f g])]
    (case [(:kind f) (:kind g)]
      [:corner :corner]
      [{:o (polygon/lerp (:p f) (:p g) 0.5) :d (polygon/unit (perp (polygon/sub (:p g) (:p f))))}]

      [:edge :edge]
      (let [u1 (polygon/unit (polygon/sub (:b f) (:a f)))
            u2 (polygon/unit (polygon/sub (:b g) (:a g)))
            c (polygon/cross u1 u2)
            gap (polygon/cross u1 (polygon/sub (:a g) (:a f)))]
        (cond
          (> (abs c) 1e-9)
          (let [x (polygon/add (:a f) (polygon/scale u1 (/ (polygon/cross (polygon/sub (:a g) (:a f)) u2) c)))]
            [{:o x :d (polygon/unit (polygon/add u1 u2))}
             {:o x :d (polygon/unit (polygon/sub u1 u2))}])
          (> (abs gap) 1e-9)
          [{:o (polygon/add (:a f) (polygon/scale (perp u1) (* 0.5 gap))) :d u1}]
          :else []))

      [:corner :edge]
      (let [u (polygon/unit (polygon/sub (:b g) (:a g)))
            v (polygon/sub (:p f) (:a g))
            d (polygon/cross u v)]
        (if (or ((:ends g) (:id f)) (< (abs d) 1e-9))
          []
          [{:a (:a g) :u u :sv (polygon/dot u v) :dv d}])))))

(defn- point-on [{:keys [o d a u sv dv]} t]
  (if o
    (polygon/add o (polygon/scale d t))
    (let [h (/ (+ (* (- t sv) (- t sv)) (* dv dv)) (* 2.0 dv))]
      (polygon/add (polygon/add a (polygon/scale u t)) (polygon/scale (perp u) h)))))

(defn- span
  "The parameters for which `curve` might be in the room: a line clipped
  to it, a parabola as far either side of its vertex as rises a room's
  diagonal."
  [{:keys [o d sv dv]} world]
  (let [[[x0 y0] [x1 y1]] (:bounds world)]
    (if o
      (let [slab (fn [o d lo hi]
                   (if (< (abs d) 1e-12)
                     (when (<= lo o hi) [##-Inf ##Inf])
                     (sort [(/ (- lo o) d) (/ (- hi o) d)])))
            [ax bx] (slab (o 0) (d 0) x0 x1)
            [ay by] (slab (o 1) (d 1) y0 y1)]
        (when (and ax ay)
          (let [t0 (max ax ay) t1 (min bx by)] (when (< t0 t1) [t0 t1]))))
      (let [r (math/sqrt (* 2.0 (abs dv) (diagonal world)))]
        [(- sv r) (+ sv r)]))))

;; ------------------------------------------------------------ the roadmap

(defn- nearest-pair?
  "Whether `p` is as far from `f` as from `g`, at least `radius` from
  both, in the free space, and no nearer any other feature. The slack
  for that is tight: where a curve ends a third feature's distance
  departs from the two's only to second order, and a slack of e lets the
  curve run on by sqrt(2 r e)."
  [world fs radius f g p]
  (let [df (distance-to f p) dg (distance-to g p) r (max df dg)]
    (and (< r ##Inf)
         (<= (abs (- df dg)) (* 1e-9 (+ 1.0 r)))
         (>= r radius)
         (free-point? world p)
         (every? (fn [h] (or (identical? h f) (identical? h g)
                             (>= (distance-to h p) (- r (* 1e-14 (+ 1.0 r))))))
                 fs))))

(defn- narrow-step
  "Halve the bracket [in out] between a parameter that passes and one
  that fails."
  [ok? [in out]]
  (when (> (abs (- out in)) (* 1e-12 (max 1.0 (abs in))))
    (let [m (* 0.5 (+ in out))]
      (if (ok? m) [m out] [in m]))))

(defn- narrow [ok? in out]
  (first (last (search/step-seq (partial narrow-step ok?) [in out]))))

(defn- runs
  "The intervals of [t0, t1] where `ok?` holds, from `n` samples, each end
  narrowed by bisection."
  [ok? t0 t1 n]
  (let [ts (mapv #(+ t0 (* (- t1 t0) (/ % n))) (range (inc n)))
        flags (mapv ok? ts)]
    (->> (range (inc n))
         (partition-by flags)
         (filter #(flags (first %)))
         (map (fn [ks]
                (let [i (first ks) j (last ks)]
                  [(if (pos? i) (narrow ok? (ts i) (ts (dec i))) (ts i))
                   (if (< j n) (narrow ok? (ts j) (ts (inc j))) (ts j))]))))))

(defn- polyline-length [pts]
  (reduce + (map (fn [[a b]] (polygon/distance a b)) (partition 2 1 pts))))

(defn- skip-pair?
  "Pairs whose curve cannot be on the roadmap: a feature with itself, an
  edge with its own corner, two features of one convex obstacle."
  [convex f g]
  (or (= (:id f) (:id g))
      (and (= (:owner f) (:owner g)) (convex (:owner f)))))

(defn pieces
  "The roadmap's curves, each `{:points :between}`: the polyline of a
  straight piece's two ends or a parabolic piece's samples, and the two
  features it keeps equidistant from. `:radius` (default 0) is the least
  clearance kept, `:samples` (default 160) how finely each curve is
  walked before its ends are narrowed."
  [world & {:keys [radius samples] :or {radius 0.0 samples 160}}]
  (let [fs (features world)
        convex (set (keep-indexed (fn [i p] (when (polygon/convex? p) i)) (:obstacles world)))]
    (vec (for [[i f] (map-indexed vector fs)
               g (subvec fs (inc i))
               :when (not (skip-pair? convex f g))
               curve (bisectors f g)
               :let [[t0 t1] (span curve world)]
               :when t0
               :let [ok? #(nearest-pair? world fs radius f g (point-on curve %))]
               [ta tb] (runs ok? t0 t1 samples)
               :let [pts (if (:o curve)
                           [(point-on curve ta) (point-on curve tb)]
                           (mapv #(point-on curve (+ ta (* (- tb ta) (/ % 24)))) (range 25)))]
               :when (> (polyline-length pts) (* 1e-4 (diagonal world)))]
           {:points pts :between [(:id f) (:id g)]}))))

(defn- node-of
  "The index of the node at `p` among `nodes`, adding it if none is
  within `tol`."
  [nodes p tol]
  (if-let [i (first (filter #(<= (polygon/distance (nodes %) p) tol) (range (count nodes))))]
    [nodes i]
    [(conj nodes p) (count nodes)]))

(defn roadmap
  "The maximum-clearance roadmap: `:pieces`, `:nodes` -- the pieces' ends,
  where three or more features are equally near or a piece stops -- and
  `:graph`, from each node to its neighbors with the length between them
  and the piece that goes there."
  [world & opts]
  (let [ps (apply pieces world opts)
        tol (* 1e-5 (diagonal world))
        [nodes ends] (reduce (fn [[nodes ends] {:keys [points]}]
                               (let [[nodes a] (node-of nodes (first points) tol)
                                     [nodes b] (node-of nodes (peek points) tol)]
                                 [nodes (conj ends [a b])]))
                             [[] []] ps)
        ps (mapv (fn [p [a b]] (assoc p :from a :to b :length (polyline-length (:points p)))) ps ends)
        graph (reduce (fn [g [k {:keys [from to length]}]]
                        (if (and (not= from to) (< length (get-in g [from to :length] ##Inf)))
                          (-> g (assoc-in [from to] {:length length :piece k})
                              (assoc-in [to from] {:length length :piece k}))
                          g))
                      (into {} (map (fn [i] [i {}])) (range (count nodes)))
                      (map-indexed vector ps))]
    {:world world :features (features world) :pieces ps :nodes nodes :graph graph}))

;; ------------------------------------------------------------- retraction

(defn- nearest-feature [fs p]
  (apply min-key #(distance-to % p) fs))

(defn- foot [f p]
  (if (= :corner (:kind f))
    (:p f)
    (polygon/lerp (:a f) (:b f) (polygon/closest-on-segment p (:a f) (:b f)))))

(defn retract
  "Where `p` meets the roadmap moving straight away from its nearest
  boundary point, or nil from a point on the boundary or not free."
  [{:keys [world features]} p]
  (when (free-point? world p)
    (let [f (nearest-feature features p)
          c (foot f p)]
      (when (> (polygon/distance p c) 1e-12)
        (let [w (polygon/unit (polygon/sub p c))
              at #(polygon/add p (polygon/scale w %))
              ok? (fn [l] (let [q (at l) r (distance-to f q)]
                            (and (free-point? world q)
                                 (every? #(or (identical? % f) (> (distance-to % q) r)) features))))
              step (/ (diagonal world) 400.0)
              out (first (filter #(not (ok? %)) (map #(* step %) (range 1 401))))]
          (if (or (nil? out) (not (ok? 0.0)))
            p
            (at (narrow ok? 0.0 out))))))))

;; ----------------------------------------------------------------- routes

(defn- attach
  "Where on the roadmap's pieces `r` lies: the piece, the segment of it
  and the point there."
  [{:keys [pieces]} r]
  (apply min-key :gap
         (for [[k {:keys [points]}] (map-indexed vector pieces)
               [i [a b]] (map-indexed vector (partition 2 1 points))
               :let [q (polygon/lerp a b (polygon/closest-on-segment r a b))]]
           {:piece k :segment i :point q :gap (polygon/distance q r)})))

(defn- toward
  "The part of piece `k`'s polyline from the attachment `at` to its end
  node `n`."
  [{:keys [pieces]} {:keys [piece segment point]} n]
  (let [{:keys [points from]} (pieces piece)]
    (vec (cons point (if (= n from)
                       (rseq (subvec points 0 (inc segment)))
                       (subvec points (inc segment)))))))

(defn- between
  "Piece `k`'s polyline from node `a` to node `b`."
  [{:keys [pieces graph]} a b]
  (let [{:keys [points from]} (pieces (get-in graph [a b :piece]))]
    (if (= a from) points (vec (rseq points)))))

(defn- within
  "Along one piece from attachment `ta` to attachment `tb`."
  [{:keys [pieces]} ta tb]
  (let [points (:points (pieces (:piece ta))) i (:segment ta) j (:segment tb)]
    (concat [(:point ta)]
            (cond (< i j) (subvec points (inc i) (inc j))
                  (> i j) (rseq (subvec points (inc j) (inc i)))
                  :else [])
            [(:point tb)])))

(defn route
  "The shortest way along the roadmap from `a` to `b`, retracting each
  onto it, as points; nil if either cannot reach it or it does not join
  them."
  [{:keys [graph pieces] :as rm} a b]
  (let [ra (retract rm a) rb (retract rm b)]
    (when (and ra rb (seq pieces))
      (let [ta (attach rm ra) tb (attach rm rb)
            ends (fn [t] (let [{:keys [from to]} (pieces (:piece t))] [from to]))
            leg (fn [t n] (polyline-length (toward rm t n)))
            g (as-> (update-vals graph #(update-vals % :length)) g
                (reduce (fn [g n] (assoc-in g [:a n] (leg ta n))) g (ends ta))
                (reduce (fn [g n] (assoc-in g [n :b] (leg tb n))) g (ends tb))
                (if (= (:piece ta) (:piece tb))
                  (let [{:keys [from]} (pieces (:piece ta))]
                    (assoc-in g [:a :b] (abs (- (leg ta from) (leg tb from)))))
                  g))]
        (when-let [nodes (some->> (search/dijkstra-seq {:start :a :goal :b
                                                        :neighbors #(keys (g %))
                                                        :cost (fn [x y] (get-in g [x y]))})
                                  search/goal-state search/path)]
          (vec (dedupe
                (concat [a ra]
                        (if (= nodes [:a :b])
                          (within rm ta tb)
                          (concat (toward rm ta (second nodes))
                                  (mapcat (fn [[x y]] (between rm x y)) (partition 2 1 (butlast (rest nodes))))
                                  (reverse (toward rm tb (last (butlast nodes))))))
                        [rb b]))))))))

(defn clearance
  "How far `p` is from the nearest feature."
  [{:keys [features]} p]
  (reduce min ##Inf (map #(distance-to % p) features)))
