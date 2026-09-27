(ns allgo.geometry.dual-mesh
  "The Delaunay triangulation and its Voronoi dual, as one graph you can
  walk.

  `allgo.geometry.delaunay` gives triangles, and `bounded-cells` gives
  polygons. Neither gives you the thing a polygon map is actually built
  on, which is *adjacency*: which cell borders which, which corners a
  cell has, which corners are joined to which, and which cells lie on
  either side of a given corner-to-corner edge. A heightmap does not need
  any of that because a grid's neighbors are arithmetic. A polygon map
  needs all of it, because every interesting pass over one -- flood-fill
  the ocean, raise the land away from the coast, run a river downhill,
  spread moisture inland -- is a graph traversal.

  ## Three kinds of thing, after Amit Patel

  The naming follows the polygon-map article this exists to support, so
  that the passes in `allgo.procedural.island` read like the description
  of them:

      centers   one per input point. The Voronoi cells, the polygons a
                map is drawn as, and where biomes live.
      corners   the Voronoi vertices, which are Delaunay circumcenters.
                Elevation and water flow live here, not on the centers --
                water runs along cell boundaries, not through middles.
      edges     each one joins two centers *and* two corners, because it
                is a Delaunay edge and a Voronoi edge at the same time.
                That is the whole point of the dual: a river crosses the
                Voronoi edge, and the two centers it separates are the
                two banks.

  Everything is a vector indexed by id, and every reference between them
  is an id. No cycles to print, no records, and the whole mesh is an
  ordinary value you can hold.

  ## The outside

  A Voronoi cell on the convex hull is unbounded, and a map made of
  half-open polygons is no use. The fix is `delaunay/bounded-cells`'
  fix -- ring the points with sentinel sites so every real site becomes
  interior -- with one change: those sentinels sit a little under one
  cell-spacing outside the bounds rather than three diagonals away. The
  far ring is right when you are going to clip the result to the
  rectangle anyway; here the boundary corners *are* the edge of the map
  and want to land near it rather than out in the middle distance.

  Sentinel centers are dropped from the finished mesh. What is left of
  them is `:border?`, on the real centers that neighbored one and on the
  corners that touched one, which is the map's edge and the thing the
  ocean is flooded from.

  ## Cost

  `delaunay/triangulate` indexes its triangles by circumcircle, so it is
  near enough linear and no longer the ceiling it was: 23ms for a
  thousand points and 136ms for four thousand on the JVM, against 83ms
  and 1.2 seconds when it scanned every triangle for every point. In a
  browser the gap is wider still. What is left is this namespace's own
  assembly, which is one pass over the triangles and is now the larger
  half of building a mesh."
  (:require [allgo.geometry.delaunay :as delaunay]
            [allgo.math :as am]
            [clojure.math :as math]))

(defn- angle-around [[ox oy] [x y]] (math/atan2 (- y oy) (- x ox)))

(defn- centroid
  "Area centroid of a ring, or its mean vertex if it has no area."
  [ring]
  (let [pts (vec ring)
        n (count pts)
        [a cx cy] (loop [i 0 a 0.0 cx 0.0 cy 0.0]
                    (if (= i n)
                      [a cx cy]
                      (let [[x0 y0] (pts i)
                            [x1 y1] (pts (mod (inc i) n))
                            c (- (* x0 y1) (* x1 y0))]
                        (recur (inc i) (+ a c)
                               (+ cx (* (+ x0 x1) c))
                               (+ cy (* (+ y0 y1) c))))))]
    (if (< (abs a) 1e-12)
      [(/ (reduce + (map first pts)) n) (/ (reduce + (map second pts)) n)]
      [(/ cx (* 3.0 a)) (/ cy (* 3.0 a))])))

(defn spacing
  "The distance between neighboring sites you would expect from `n`
  points spread over `bounds`."
  [[[x0 y0] [x1 y1]] n]
  (math/sqrt (/ (* (- x1 x0) (- y1 y0)) (max 1 (long n)))))

(defn boundary-ring
  "Sentinel sites ringing `bounds`, one `gap` outside it and spaced `gap`
  apart.

  Close in, unlike `delaunay/bounded-cells`' ring: these decide where the
  outermost real corners land, and corners are the edge of the map."
  [[[x0 y0] [x1 y1]] gap]
  (let [gap (double gap)
        ax (- x0 gap) ay (- y0 gap)
        bx (+ x1 gap) by (+ y1 gap)
        along (fn [a b] (let [steps (max 2 (long (math/ceil (/ (- b a) gap))))]
                          (map #(+ a (* (- b a) (/ (double %) steps))) (range (inc steps)))))]
    (distinct
     (concat (for [x (along ax bx)] [x ay])
             (for [x (along ax bx)] [x by])
             (for [y (along ay by)] [ax y])
             (for [y (along ay by)] [bx y])))))

(defn- corner-key
  "Circumcenters that agree to within a whisker are the same corner.

  Four cocircular sites give two triangles with the same circumcenter,
  and floating point will not make them equal. Left as two corners, the
  cell around them has a zero-length edge and the traversals get a
  neighbor that is really themselves."
  [[x y] eps]
  [(math/round (/ (double x) eps)) (math/round (/ (double y) eps))])

(defn mesh
  "The dual graph of `points` within `bounds` (`[[x0 y0] [x1 y1]]`).

  Returns `{:centers [...] :corners [...] :edges [...] :bounds ...}`, each
  vector indexed by id:

      center  {:id :point :neighbors :borders :corners :border?}
      corner  {:id :point :touches :protrudes :adjacent :border?}
      edge    {:id :centers [c0 c1] :corners [v0 v1] :midpoint}

  An edge's `:centers` holds `nil` where the cell on that side was a
  sentinel, so an edge on the rim of the map reads as `[id nil]`.

  A center's `:corners` are ordered counter-clockwise, so they are the
  polygon to draw it as. `:border?` marks the edge of the map."
  [points bounds]
  (let [pts (vec points)
        n (count pts)
        gap (spacing bounds n)
        ring (vec (boundary-ring bounds gap))
        all (into pts ring)
        site-id (into {} (map-indexed (fn [i p] [p i])) all)
        sentinel? (fn [^long id] (>= id n))
        eps (* 1e-7 (let [[[x0 y0] [x1 y1]] bounds]
                      (max (- x1 x0) (- y1 y0))))
        tris (vec (delaunay/triangulate all))
        n-tris (count tris)
        ;; One corner per triangle, merged where circumcenters coincide.
        ;;
        ;; Transients throughout this section, and not for tidiness. Built
        ;; with `update-in` over persistent maps of sets it cost about as
        ;; much as the triangulation it was consuming -- some 1.1 seconds
        ;; for seven hundred points in a browser, on top of the 1.1 the
        ;; triangulation already takes. There are only a few thousand
        ;; triangles; it is the four or five map rewrites apiece that
        ;; added up.
        [corner-pts corner-of]
        (loop [i 0 pts (transient []) by-key (transient {}) cof (transient [])]
          (if (= i n-tris)
            [(persistent! pts) (persistent! cof)]
            (let [c (get-in (tris i) [:circle :center])
                  k (corner-key c eps)]
              (if-let [existing (get by-key k)]
                (recur (inc i) pts by-key (conj! cof existing))
                (let [id (count pts)]
                  (recur (inc i) (conj! pts c) (assoc! by-key k id) (conj! cof id)))))))
        ;; Which sites each corner touches, and which corners sit on each
        ;; Delaunay edge. A Delaunay edge is shared by two triangles, and
        ;; those two circumcenters are the Voronoi edge's endpoints.
        [touches on-edge]
        (loop [i 0 tch (transient {}) oe (transient {})]
          (if (= i n-tris)
            [(persistent! tch) (persistent! oe)]
            (let [cid (corner-of i)
                  [p0 p1 p2] (:points (tris i))
                  a (site-id p0) b (site-id p1) c (site-id p2)
                  add (fn [oe u v]
                        (let [k (if (< (long u) (long v)) [u v] [v u])]
                          (assoc! oe k (conj (get oe k []) cid))))]
              (recur (inc i)
                     (assoc! tch cid (conj (get tch cid []) a b c))
                     (-> oe (add a b) (add b c) (add c a))))))
        ;; Keep the Voronoi edges that are real: two distinct corners, and
        ;; at least one real center to be a boundary of.
        raw-edges (into []
                        (comp (map (fn [[k cs]] [k (vec (distinct cs))]))
                              (filter (fn [[[a b] cs]]
                                        (and (= 2 (count cs))
                                             (or (not (sentinel? a)) (not (sentinel? b))))))
                              (map (fn [[[a b] cs]]
                                     {:centers [a b] :corners cs})))
                        on-edge)
        ;; Real centers only; a center that neighbored a sentinel is the
        ;; map's edge.
        border-centers (into #{}
                             (mapcat (fn [{[a b] :centers}]
                                       (cond (sentinel? a) [b]
                                             (sentinel? b) [a]
                                             :else [])))
                             raw-edges)
        border-corners (into #{}
                             (mapcat (fn [[cid sids]]
                                       (when (some sentinel? sids) [cid])))
                             touches)
        ;; Sentinel centers are gone from the finished mesh, so an edge
        ;; that had one keeps `nil` in its place rather than an id that
        ;; indexes nothing. An edge on the rim of the map has one side and
        ;; no other, and a renderer walking `:centers` would otherwise
        ;; reach past the end of the vector to find out.
        edges (into [] (map-indexed
                        (fn [i e]
                          (let [[v0 v1] (:corners e)
                                [x0 y0] (corner-pts v0)
                                [x1 y1] (corner-pts v1)
                                [a b] (:centers e)]
                            (assoc e
                                   :id i
                                   :centers [(when-not (sentinel? a) a)
                                             (when-not (sentinel? b) b)]
                                   :midpoint [(* 0.5 (+ x0 x1)) (* 0.5 (+ y0 y1))]))))
                    raw-edges)
        ;; Fan out the adjacency the passes will walk.
        center-edges (reduce (fn [m {:keys [id centers]}]
                               (reduce (fn [m c] (if (nil? c) m (update m c (fnil conj []) id)))
                                       m centers))
                             {} edges)
        corner-edges (reduce (fn [m {:keys [id corners]}]
                               (reduce (fn [m c] (update m c (fnil conj []) id)) m corners))
                             {} edges)
        centers (into []
                      (map-indexed
                       (fn [id p]
                         (let [eids (vec (center-edges id))
                               nbrs (into [] (comp (map #(:centers (edges %)))
                                                   (map (fn [[a b]] (if (= a id) b a)))
                                                   (remove nil?)
                                                   (distinct))
                                          eids)
                               cs (into [] (comp (map #(:corners (edges %)))
                                                 cat
                                                 (distinct))
                                        eids)]
                           {:id id
                            :point p
                            :neighbors nbrs
                            :borders eids
                            :corners (vec (sort-by #(angle-around p (corner-pts %)) cs))
                            :border? (contains? border-centers id)})))
                      pts)
        corners (into []
                      (map-indexed
                       (fn [id p]
                         (let [eids (vec (corner-edges id))]
                           {:id id
                            :point p
                            :touches (into [] (comp (remove sentinel?) (distinct))
                                           (get touches id []))
                            :protrudes eids
                            :adjacent (into [] (comp (map #(:corners (edges %)))
                                                     (map (fn [[a b]] (if (= a id) b a)))
                                                     (distinct))
                                            eids)
                            :border? (contains? border-corners id)})))
                      corner-pts)]
    {:centers centers :corners corners :edges edges :bounds bounds}))

(defn polygon
  "The ring of points to draw `center` as."
  [{:keys [corners]} center]
  (mapv #(:point (corners %)) (:corners center)))

(defn relax
  "One round of Lloyd relaxation: every point moves to the center of area
  of its cell.

  Random points clump, and clumped points make cells of wildly different
  sizes -- which shows, because a polygon map is drawn as its cells. Two
  or three rounds of this is the usual dose: enough to even the sizes out,
  few enough that the result still looks organic rather than like a
  honeycomb.

  Returns points, not a mesh, since it is the points that have moved."
  [{:keys [centers bounds] :as m}]
  (let [[[x0 y0] [x1 y1]] bounds]
    (mapv (fn [c]
            (let [ring (polygon m c)]
              (if (< (count ring) 3)
                (:point c)
                (let [[cx cy] (centroid ring)]
                  [(am/clamp cx x0 x1) (am/clamp cy y0 y1)]))))
          centers)))

(defn relaxed
  "A mesh over `points` after `rounds` of `relax`."
  [points bounds rounds]
  (loop [pts (vec points) i 0]
    (let [m (mesh pts bounds)]
      (if (>= i (long rounds))
        m
        (recur (relax m) (inc i))))))

;; ---------------------------------------------------------------------------
;; Noisy edges
;;
;; A Voronoi cell drawn as its corners is a polygon, and a map made of
;; them looks like one: the coastline is a sequence of straight lines
;; meeting at angles, which is the one thing that says "generated" no
;; matter how good the biomes are. Amit Patel's fix is to draw each
;; boundary as a wandering path rather than a segment, and it is worth
;; being clear that this changes nothing about the model -- the cells,
;; the adjacency, the elevation and the rivers are all exactly as they
;; were. It is a change to where the ink goes.
;;
;; The path for an edge is stored *on the edge*, which is the whole trick.
;; Both cells sharing it draw the same points, so the two still tile with
;; no gap and no overlap; noise generated per cell would tear them apart
;; along every boundary.
;;
;; And each path is confined to the quadrilateral made by the edge's two
;; corners and the two cell centers either side of it. A path that
;; wandered outside would cross into a cell that is not one of the two it
;; separates, and a coastline would start overlapping land two cells
;; inland. The quad is the guarantee that it cannot.

(defn- lerp [[ax ay] [bx by] t]
  [(+ ax (* (double t) (- bx ax)))
   (+ ay (* (double t) (- by ay)))])

(defn- dist [[ax ay] [bx by]]
  (let [dx (- (double bx) ax) dy (- (double by) ay)]
    (math/sqrt (+ (* dx dx) (* dy dy)))))

(defn- subdivide
  "The points strictly between `a` and `c`, wandering inside quad
  `a b c d`.

  Recursive midpoint displacement, but displaced along the quad rather
  than along a perpendicular: cut the quad into two smaller quads that
  meet at an interior point `h`, keep `h`, and recurse into both. Because
  every child quad is inside its parent, every point produced is inside
  the original -- which is the containment the whole scheme relies on,
  and it holds by construction rather than by clamping."
  [a b c d min-length rng]
  (if (or (< (dist a c) min-length) (< (dist b d) min-length))
    []
    (let [;; Where the quad gets cut. Away from the ends, or the wiggle
          ;; degenerates into a kink at a corner.
          p (+ 0.2 (* 0.6 (rng)))
          q (+ 0.2 (* 0.6 (rng)))
          e (lerp a d p)
          f (lerp b c p)
          g (lerp a b q)
          i (lerp d c q)
          h (lerp e f q)
          ;; Allowed slightly past the ends, so the child quads are not
          ;; strictly nested and the path does not visibly taper toward
          ;; the middle of every edge.
          s (+ 0.6 (* 0.8 (rng)))
          t (+ 0.6 (* 0.8 (rng)))]
      (-> []
          (into (subdivide a (lerp b g s) h (lerp d e t) min-length rng))
          (conj h)
          (into (subdivide h (lerp c f s) c (lerp d i t) min-length rng))))))

(defn- edge-path
  "One edge's wandering path, from its first corner to its second."
  [{:keys [centers corners]} {:keys [corners' centers' midpoint]} min-length rng]
  (let [[v0 v1] corners'
        [d0 d1] centers'
        p0 (:point (corners v0))
        p1 (:point (corners v1))]
    (if (or (nil? d0) (nil? d1))
      ;; An edge on the rim of the map has a cell on one side only, so
      ;; there is no quad to stay inside. Left straight, which is right:
      ;; it is the edge of the map rather than a feature of it.
      [p0 p1]
      (let [c0 (:point (centers d0))
            c1 (:point (centers d1))
            ;; Halfway to each center: the quad is the cell pair's shared
            ;; neighborhood, not the cells themselves.
            half (fn [v c] (lerp v c 0.5))
            first-half (into [p0] (conj (subdivide p0 (half p0 c0) midpoint (half p0 c1)
                                                   min-length rng)
                                        midpoint))
            second-half (into [p1] (conj (subdivide p1 (half p1 c1) midpoint (half p1 c0)
                                                    min-length rng)
                                         midpoint))]
        ;; Both halves run outward from their corner to the midpoint, so
        ;; the second is reversed and its duplicate midpoint dropped.
        (into first-half (rest (reverse second-half)))))))

(defn noisy-edges
  "Gives every edge a `:path`: the points to draw it as, from its first
  corner to its second.

  `:min-length` is where the subdivision stops, and it is the only knob
  that matters -- large and the edges stay nearly straight, small and the
  point count doubles for detail finer than anyone will see. The default
  is a fraction of the distance between neighboring sites, so it follows
  the mesh rather than the units it happens to be in.

  Costs one pass over the edges and leaves everything else untouched, so
  it can be applied at any point after the mesh is built."
  ([mesh] (noisy-edges mesh {}))
  ([{:keys [centers edges bounds] :as mesh} {:keys [min-length rng] :or {rng rand}}]
   (let [min-length (double (or min-length (* 0.12 (spacing bounds (count centers)))))]
     (assoc mesh :edges
            (mapv (fn [e]
                    (assoc e :path (edge-path mesh
                                              {:corners' (:corners e)
                                               :centers' (:centers e)
                                               :midpoint (:midpoint e)}
                                              min-length rng)))
                  edges)))))

(defn- edge-joining
  "The edge of `center` whose two corners are `a` and `b`."
  [{:keys [edges]} center a b]
  (some (fn [eid]
          (let [[v0 v1] (:corners (edges eid))]
            (when (or (and (= v0 a) (= v1 b)) (and (= v0 b) (= v1 a))) eid)))
        (:borders center)))

(defn noisy-polygon
  "The ring to draw `center` as, following each edge's `:path`.

  Falls back to `polygon` where `noisy-edges` has not been applied, so a
  renderer can call this unconditionally."
  [{:keys [corners edges] :as mesh} center]
  (let [cs (:corners center)]
    (if (or (< (count cs) 3) (nil? (:path (edges (first (:borders center))))))
      (polygon mesh center)
      (into []
            (mapcat (fn [[a b]]
                      (if-let [eid (edge-joining mesh center a b)]
                        (let [{:keys [path] [v0 _] :corners} (edges eid)]
                          ;; Each segment drops its last point; the next
                          ;; one starts there.
                          (butlast (if (= v0 a) path (reverse path))))
                        [(:point (corners a))])))
            (partition 2 1 (conj (vec cs) (first cs)))))))
