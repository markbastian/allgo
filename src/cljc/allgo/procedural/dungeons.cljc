(ns allgo.procedural.dungeons
  "Dungeon generation after the algorithm TinyKeep's developer described and
  A. Adonaac wrote up for Gamasutra:

    http://www.gamasutra.com/blogs/AAdonaac/20150903/252889/Procedural_Dungeon_Generation_Algorithm.php
    https://web.archive.org/web/20210823235922/https://gamasutra.com/blogs/AAdonaac/20150903/252889/Procedural_Dungeon_Generation_Algorithm.php

  The pipeline, each step a function you can call on its own:

    1. `scatter`    -- strew rooms of random size over a circle (or ellipse)
    2. `separate`   -- push them apart until none overlap
    3. `main-rooms` -- keep the big ones as hubs
    4. `room-graph` -- Delaunay-triangulate the hubs
    5. `connect`    -- take the minimum spanning tree (`allgo.graph`), then
                       add a few edges back so the dungeon has loops rather
                       than a single path
    6. `corridors`  -- turn each edge into a straight or L-shaped run
    7. `absorb`     -- promote the small rooms a corridor passes through
    8. `rasterize`  -- flatten the whole thing onto a tile grid

  `generate` runs all eight and hands back the rooms, the graph, the
  corridor segments and the grid.

  Everything is in world units. `:tile-size` is the grid those units snap
  to, so a tile size of 4 means every position and extent is a multiple of
  4; grid keys, by contrast, are tile indices."
  (:require [allgo.geometry.delaunay :as delaunay]
            [allgo.graph :as graph]
            [clojure.math :as math]))

(def default-config
  "Adonaac's numbers where he gives them: tile size 4, a width and height
  mean of 24, hubs at 1.25x the mean. He preferred adding 8-10% of the
  edges back where TinyKeep used 15%, so the default sits in his range."
  {:room-count       150
   :tile-size        4
   :radius           120.0
   :ellipse          nil
   :width-mean       24.0
   :width-sd         8.0
   :height-mean      24.0
   :height-sd        8.0
   :min-size         8.0
   :main-threshold   1.25
   :min-hubs         3
   :extra-edge-ratio 0.10
   :corridor-width   3
   :max-iterations   400})

;; ---------------------------------------------------------------------------
;; Random helpers

(defn roundm
  "Rounds `n` up to the next multiple of `m`, Adonaac's `roundm`. Everything
  that lands on the map goes through this, so rooms stay on the tile grid."
  [n m]
  (int (* (math/floor (/ (dec (+ n m)) m)) m)))

(defn random-point
  "A point drawn uniformly from the disc of radius `radius`.

  Two uniform draws are summed and the result folded back on itself past 1,
  which turns the flat distribution into a triangular one. That extra
  weight toward the rim is exactly what cancels the crowding you would
  otherwise get from sampling the radius directly."
  [radius]
  (let [t (* 2.0 math/PI (rand))
        u (+ (rand) (rand))
        r (if (> u 1) (- 2 u) u)]
    [(* radius r (math/cos t)) (* radius r (math/sin t))]))

(defn random-int-point
  "`random-point`, snapped to the tile grid."
  [radius tile-size]
  (let [t (* 2.0 math/PI (rand))
        u (+ (rand) (rand))
        r (if (> u 1) (- 2 u) u)]
    [(roundm (* radius r (math/cos t)) tile-size)
     (roundm (* radius r (math/sin t)) tile-size)]))

(defn random-ellipse-point
  "A point in the ellipse `ellipse-width` by `ellipse-height`, snapped to
  the tile grid.

  Scattering into a wide, flat ellipse instead of a circle is Adonaac's fix
  for rooms that are much wider than they are tall: separating those on a
  circle resolves most collisions vertically and grows a dungeon far taller
  than it is wide."
  [ellipse-width ellipse-height tile-size]
  (let [t (* 2.0 math/PI (rand))
        u (+ (rand) (rand))
        r (if (> u 1) (- 2 u) u)]
    [(roundm (/ (* ellipse-width r (math/cos t)) 2.0) tile-size)
     (roundm (/ (* ellipse-height r (math/sin t)) 2.0) tile-size)]))

(defn gaussian
  "One draw from a normal distribution, by the Box-Muller transform.

  Room sizes are normal rather than uniform because it gives you a mean and
  a spread to tune separately: the spread is what decides whether the hub
  rooms stand out from the crowd enough to be worth singling out."
  [mean sd]
  (let [u1 (max 1e-12 (rand))
        u2 (rand)]
    (+ mean (* sd (math/sqrt (* -2.0 (math/log u1))) (math/cos (* 2.0 math/PI u2))))))

(defn random-room
  "A room with a uniformly random size, centred somewhere in the disc of
  radius `radius`. Kept for callers that want one loose room; `scatter`
  is what the pipeline uses."
  [radius max-width max-height]
  {:center (random-point radius)
   :width  (* (rand) max-width)
   :height (* (rand) max-height)})

(defn center
  "Centroid of a collection of rooms."
  [rooms]
  (->> rooms
       (map :center)
       (apply map +)
       (mapv #(/ % (count rooms)))))

;; ---------------------------------------------------------------------------
;; Rectangles

(defn bounds
  "`[min-x min-y max-x max-y]` of a room."
  [{:keys [center width height]}]
  (let [[cx cy] center
        hw      (/ width 2.0)
        hh      (/ height 2.0)]
    [(- cx hw) (- cy hh) (+ cx hw) (+ cy hh)]))

(defn overlaps?
  "Whether two rooms share any interior area."
  [a b]
  (let [[ax0 ay0 ax1 ay1] (bounds a)
        [bx0 by0 bx1 by1] (bounds b)]
    (and (< ax0 bx1) (< bx0 ax1) (< ay0 by1) (< by0 ay1))))

;; ---------------------------------------------------------------------------
;; 1. Scatter

(defn- snap-size
  "Rounds an extent to an *even* number of tiles, to nearest.

  Even, because rooms are positioned by their centre: a room an odd number
  of tiles wide has its edges half a tile off the grid, and two such rooms
  placed flush straddle a shared column of tiles and bleed into each other
  when rasterized. An even extent has a whole-tile half-extent, so
  grid-aligned centres give grid-aligned edges.

  To nearest rather than `roundm`'s round-up, because the threshold that
  picks out hub rooms is a multiple of the mean size. Rounding up inflates
  every draw, which lifts the mean, which lifts the threshold -- and it
  rises faster than the sizes do, so almost nothing clears it and the
  dungeon ends up with one or two hubs instead of a dozen."
  [v min-size tile-size]
  (let [step  (* 2 tile-size)
        round (fn [x] (* step (long (math/round (/ (double x) step)))))]
    (max (round min-size) (round v) step)))

(defn scatter
  "Step 1. `room-count` rooms with normally distributed extents, centred on
  points drawn from the scatter region -- the disc of radius `:radius`, or
  the ellipse `:ellipse` when one is given. Sizes and positions are both
  snapped to the tile grid."
  [{:keys [room-count radius ellipse tile-size
           width-mean width-sd height-mean height-sd min-size]}]
  (mapv (fn [id]
          {:id     id
           :center (if ellipse
                     (random-ellipse-point (first ellipse) (second ellipse) tile-size)
                     (random-int-point radius tile-size))
           :width  (snap-size (gaussian width-mean width-sd) min-size tile-size)
           :height (snap-size (gaussian height-mean height-sd) min-size tile-size)})
        (range room-count)))

;; ---------------------------------------------------------------------------
;; 2. Separate

(def ^:private touching-tolerance
  "Rooms pushed flush against each other still overlap by a whisker, because
  their edges come out of a subtraction and the corrections converge
  geometrically rather than landing on zero. Left at `pos?`, that residue
  reads as a live collision and the solver sweeps forever over a layout
  that has long since settled. A micro-unit is far below the tile grid
  everything ends up snapped to."
  1e-6)

(defn- push-apart
  "The shortest translation that would lift room `a` clear of room `b`, or
  nil when they are already disjoint. Always along one axis: the one they
  overlap on least, which is what keeps rooms sliding apart into a layout
  rather than jittering diagonally."
  [a b]
  (let [[ax ay] (:center a)
        [bx by] (:center b)
        dx      (- ax bx)
        dy      (- ay by)
        ox      (- (/ (+ (:width a) (:width b)) 2.0) (abs dx))
        oy      (- (/ (+ (:height a) (:height b)) 2.0) (abs dy))]
    (when (and (> ox touching-tolerance) (> oy touching-tolerance))
      ;; Exactly concentric rooms have no direction to separate along, so
      ;; break the tie on id rather than letting both push the same way.
      (let [sx (cond (pos? dx) 1.0 (neg? dx) -1.0 (< (:id a) (:id b)) -1.0 :else 1.0)
            sy (cond (pos? dy) 1.0 (neg? dy) -1.0 (< (:id a) (:id b)) -1.0 :else 1.0)]
        (if (< ox oy) [(* sx ox) 0.0] [0.0 (* sy oy)])))))

(defn- nudge [room dx dy]
  (update room :center (fn [[x y]] [(+ x dx) (+ y dy)])))

(defn- separation-pass
  "One sweep over every pair, splitting each overlap between the two rooms.
  Returns `[rooms moved?]`. Later pairs see the moves made by earlier ones,
  which converges in far fewer sweeps than collecting all the corrections
  and applying them at the end."
  [rooms]
  (let [n (count rooms)]
    (loop [i 0 rooms rooms moved? false]
      (if (>= i n)
        [rooms moved?]
        (let [[rooms moved?]
              (loop [j (inc i) rooms rooms moved? moved?]
                (if (>= j n)
                  [rooms moved?]
                  (if-let [[mx my] (push-apart (rooms i) (rooms j))]
                    (recur (inc j)
                           (-> rooms
                               (update i nudge (* 0.5 mx) (* 0.5 my))
                               (update j nudge (* -0.5 mx) (* -0.5 my)))
                           true)
                    (recur (inc j) rooms moved?))))]
          (recur (inc i) rooms moved?))))))

(defn- snap-center [room tile-size]
  ;; To nearest, not up: `roundm` would shove every room the same direction
  ;; by up to a tile, and it is the *differences* in that shove that open
  ;; new overlaps. Nearest halves the worst-case displacement.
  (update room :center
          (fn [[x y]] [(* tile-size (math/round (/ x (double tile-size))))
                       (* tile-size (math/round (/ y (double tile-size))))])))

(defn- grid-separation-pass
  "A separation sweep that moves rooms in whole tiles.

  Once positions and extents are both multiples of the tile size, overlaps
  are too, so a correction rounded up to a whole tile clears the overlap
  outright and leaves the room on the grid. Unlike the continuous sweep
  this terminates exactly rather than converging on zero from above."
  [rooms tile-size]
  (let [n (count rooms)
        quantize (fn [v] (* tile-size (math/ceil (/ (/ (abs v) 2.0) tile-size))))]
    (loop [i 0 rooms rooms moved? false]
      (if (>= i n)
        [rooms moved?]
        (let [[rooms moved?]
              (loop [j (inc i) rooms rooms moved? moved?]
                (if (>= j n)
                  [rooms moved?]
                  (if-let [[mx my] (push-apart (rooms i) (rooms j))]
                    (let [dx (if (zero? mx) 0.0 (* (math/signum mx) (quantize mx)))
                          dy (if (zero? my) 0.0 (* (math/signum my) (quantize my)))]
                      (recur (inc j)
                             (-> rooms
                                 (update i nudge dx dy)
                                 (update j nudge (- dx) (- dy)))
                             true))
                    (recur (inc j) rooms moved?))))]
          (recur (inc i) rooms moved?))))))

(defn separate
  "Step 2. Sweeps until nothing overlaps, then settles the result onto the
  tile grid.

  Adonaac reached for a physics engine and let the bodies fall asleep;
  TinyKeep used separation steering. This does the same job directly --
  resolve every overlap along its shallowest axis, repeat -- which avoids a
  dependency and, unlike steering, ends with rooms actually disjoint rather
  than merely spread out.

  Two phases. The continuous one spreads the pile out and is what gives the
  layout its organic look, but its corrections shrink geometrically and
  never quite reach zero. Snapping that to the grid then reopens real
  overlaps, up to half a tile deep, because each room rounds a different
  distance -- the mismatch Adonaac notes and tolerates. The second phase
  works in whole tiles from the snapped positions and closes them, so the
  rooms that come out are grid-aligned *and* disjoint."
  [rooms {:keys [tile-size max-iterations]}]
  (let [spread  (loop [rooms (vec rooms) i 0]
                  (if (>= i max-iterations)
                    rooms
                    (let [[rooms moved?] (separation-pass rooms)]
                      (if moved? (recur rooms (inc i)) rooms))))
        snapped (mapv #(snap-center % tile-size) spread)]
    (loop [rooms snapped i 0]
      (if (>= i max-iterations)
        rooms
        (let [[rooms moved?] (grid-separation-pass rooms tile-size)]
          (if moved? (recur rooms (inc i)) rooms))))))

;; ---------------------------------------------------------------------------
;; 3. Main rooms

(defn main-rooms
  "Step 3. Tags each room `:main` or `:minor`. A room is a hub when both its
  width and its height clear `:main-threshold` times the mean -- 1.25x in
  Adonaac's write-up, so with means of 24 a hub must be over 30 on a side.

  `:min-hubs` is a floor the original does not have. The threshold is
  measured against the mean of this particular batch of rooms, and at small
  room counts that mean is noisy enough that a run can select one hub, or
  none -- which yields no edges, no corridors and an empty map. When too
  few clear the bar the biggest rooms are taken instead, so `generate`
  always returns something you could walk around in."
  [rooms {:keys [main-threshold min-hubs]}]
  (let [n       (count rooms)
        mw      (* main-threshold (/ (reduce + (map :width rooms)) n))
        mh      (* main-threshold (/ (reduce + (map :height rooms)) n))
        over    (fn [{:keys [width height]}] (and (>= width mw) (>= height mh)))
        floor   (min (or min-hubs 0) n)
        chosen  (let [qualified (filter over rooms)]
                  (if (>= (count qualified) floor)
                    (set (map :id qualified))
                    (->> rooms
                         (sort-by (fn [{:keys [width height]}] (- (* width height))))
                         (take floor)
                         (map :id)
                         set)))]
    (mapv (fn [room] (assoc room :kind (if (chosen (:id room)) :main :minor))) rooms)))

(defn hubs [rooms] (filterv (comp #{:main} :kind) rooms))

;; ---------------------------------------------------------------------------
;; 4. Delaunay graph

(defn room-graph
  "Step 4. Delaunay-triangulates the hub centres and reads the triangle
  edges off as a graph. Returns `{:nodes #{id} :edges #{#{id id}}}`.

  Delaunay is the right triangulation here because it avoids sliver
  triangles, so the edges it produces connect rooms to their genuine
  neighbours rather than skipping across the map."
  [rooms]
  (let [;; Two hubs that snapped onto the same centre would be
        ;; indistinguishable as triangulation sites, so only the first can
        ;; take part.
        by-center (reduce (fn [m {:keys [id center]}]
                            (cond-> m (not (contains? m center)) (assoc center id)))
                          {}
                          rooms)
        sites     (keys by-center)]
    (if (< (count sites) 3)
      ;; Too few points to triangulate: connect what there is directly.
      {:nodes (set (vals by-center))
       :edges (graph/complete-graph (vals by-center))}
      (let [nodes (set (vals by-center))
            edges (into #{}
                        (comp (mapcat (fn [{[a b c] :points}] [[a b] [b c] [c a]]))
                              (keep (fn [[p q]]
                                      (let [i (by-center p) j (by-center q)]
                                        (when (and i j (not= i j)) #{i j})))))
                        (delaunay/triangulate sites))]
        {:nodes nodes
         ;; Collinear hubs have no circumcircle, so the triangulation
         ;; degenerates and can hand back no triangles at all -- three
         ;; rooms in a row produce an empty graph, no corridors, and a
         ;; dungeon you cannot cross. Falling back to the complete graph
         ;; costs nothing in the normal case: the Euclidean minimum
         ;; spanning tree is a subgraph of the Delaunay triangulation, so
         ;; `connect` picks the same edges either way.
         :edges (if (graph/connected-graph? nodes edges)
                  edges
                  (graph/complete-graph nodes))}))))

;; ---------------------------------------------------------------------------
;; 5. Minimum spanning tree, plus a few edges back

(defn connect
  "Step 5. The spanning tree, plus `:extra-edge-ratio` of the edges it
  discarded put back.

  The tree alone guarantees every hub is reachable and nothing is stranded,
  but it is a tree: exactly one path between any two rooms, which plays as
  a corridor you walk down and back. Restoring a fraction of the discarded
  Delaunay edges buys loops and alternate routes. TinyKeep used 15%;
  Adonaac preferred 8-10%."
  [rooms {:keys [nodes edges]} {:keys [extra-edge-ratio]}]
  (let [by-id  (into {} (map (juxt :id identity)) rooms)
        weight (fn [e]
                 (let [[u v] (vec e)
                       [ux uy] (:center (by-id u))
                       [vx vy] (:center (by-id v))]
                   (math/hypot (- vx ux) (- vy uy))))
        tree   (graph/minimum-spanning-tree nodes edges weight)
        spare  (vec (graph/edges-not-in edges tree))
        extra  (int (math/round (* extra-edge-ratio (count spare))))]
    (into tree (take extra (shuffle spare)))))

;; ---------------------------------------------------------------------------
;; 6. Corridors

(defn- spans-x? [room x]
  (let [[x0 _ x1 _] (bounds room)] (and (>= x x0) (<= x x1))))

(defn- spans-y? [room y]
  (let [[_ y0 _ y1] (bounds room)] (and (>= y y0) (<= y y1))))

(defn corridor-segments
  "The run between two rooms, as one or two axis-aligned segments.

  Adonaac's test: take the midpoint of the two centres. If its x falls
  inside both rooms they overlap in a vertical band, so a single vertical
  run down that band joins them; likewise for y. Failing both, an L --
  out along x from the first centre, then along y into the second."
  [a b]
  (let [[ax ay] (:center a)
        [bx by] (:center b)
        mx      (/ (+ ax bx) 2.0)
        my      (/ (+ ay by) 2.0)]
    (cond
      (and (spans-x? a mx) (spans-x? b mx)) [[[mx ay] [mx by]]]
      (and (spans-y? a my) (spans-y? b my)) [[[ax my] [bx my]]]
      :else                                 [[[ax ay] [bx ay]] [[bx ay] [bx by]]])))

(defn corridors
  "Step 6. Every connecting edge turned into segments."
  [rooms edges]
  (let [by-id (into {} (map (juxt :id identity)) rooms)]
    (into [] (mapcat (fn [e]
                       (let [[u v] (vec e)]
                         (corridor-segments (by-id u) (by-id v)))))
          edges)))

;; ---------------------------------------------------------------------------
;; 7. Absorb the rooms a corridor runs through

(defn- segment-bounds
  "A segment fattened to the corridor width, as `[min-x min-y max-x max-y]`."
  [[[x1 y1] [x2 y2]] tile-size corridor-width]
  (let [half (* tile-size (/ (dec corridor-width) 2.0))]
    [(- (min x1 x2) half) (- (min y1 y2) half)
     (+ (max x1 x2) half) (+ (max y1 y2) half)]))

(defn absorb
  "Step 7. Any `:minor` room a corridor passes through becomes a
  `:hallway` room -- the irregular chambers strung along the runs, which is
  what stops corridors reading as bare lines. The rest stay `:minor` and
  are dropped from the finished map."
  [rooms segments {:keys [tile-size corridor-width]}]
  (let [boxes (mapv #(segment-bounds % tile-size corridor-width) segments)
        hits? (fn [room]
                (let [[rx0 ry0 rx1 ry1] (bounds room)]
                  (boolean
                   (some (fn [[sx0 sy0 sx1 sy1]]
                           (and (< rx0 sx1) (< sx0 rx1) (< ry0 sy1) (< sy0 ry1)))
                         boxes))))]
    (mapv (fn [room]
            (cond-> room (and (= :minor (:kind room)) (hits? room))
                    (assoc :kind :hallway)))
          rooms)))

;; ---------------------------------------------------------------------------
;; 8. Rasterize

(defn- tile-range [lo hi tile-size]
  (range (long (math/floor (/ lo tile-size)))
         (long (math/ceil (/ hi tile-size)))))

(defn- room-tiles [room tile-size]
  (let [[x0 y0 x1 y1] (bounds room)]
    (for [tx (tile-range x0 x1 tile-size)
          ty (tile-range y0 y1 tile-size)]
      [tx ty])))

(defn- segment-tiles [segment tile-size corridor-width]
  (let [[x0 y0 x1 y1] (segment-bounds segment tile-size corridor-width)]
    (for [tx (tile-range x0 (+ x1 tile-size) tile-size)
          ty (tile-range y0 (+ y1 tile-size) tile-size)]
      [tx ty])))

(defn rasterize
  "Step 8. Flattens rooms and corridors onto a tile grid: a map from
  `[tile-x tile-y]` to `:main`, `:hallway` or `:corridor`.

  Corridors go down first so a room tile always wins over the run passing
  through it, and hubs win over hallway rooms."
  [rooms segments {:keys [tile-size corridor-width]}]
  (let [put (fn [grid kind tiles] (reduce #(assoc %1 %2 kind) grid tiles))]
    (as-> {} grid
      (put grid :corridor (mapcat #(segment-tiles % tile-size corridor-width) segments))
      (put grid :hallway (mapcat #(room-tiles % tile-size)
                                 (filter (comp #{:hallway} :kind) rooms)))
      (put grid :main (mapcat #(room-tiles % tile-size)
                              (filter (comp #{:main} :kind) rooms))))))

;; ---------------------------------------------------------------------------
;; The whole pipeline

(defn generate
  "Runs every step and returns the finished dungeon:

    :rooms     every room, each `{:id :center :width :height :kind}`, where
               `:kind` is `:main`, `:hallway` or the `:minor` leftovers
    :graph     `{:nodes :edges}` -- the full Delaunay graph of the hubs
    :edges     the subset actually connected: spanning tree plus loops
    :corridors the segments those edges became
    :grid      `{[tile-x tile-y] :main|:hallway|:corridor}`
    :config    the settings used

  Takes any subset of `default-config`."
  ([] (generate {}))
  ([config]
   (let [config    (merge default-config config)
         rooms     (-> (scatter config)
                       (separate config)
                       (main-rooms config))
         graph     (room-graph (hubs rooms))
         edges     (connect rooms graph config)
         segments  (corridors rooms edges)
         rooms     (absorb rooms segments config)]
     {:rooms     rooms
      :graph     graph
      :edges     edges
      :corridors segments
      :grid      (rasterize rooms segments config)
      :config    config})))

(comment
  (let [{:keys [rooms edges grid]} (generate {:room-count 60 :radius 80.0})]
    {:rooms     (frequencies (map :kind rooms))
     :edges     (count edges)
     :grid      (count grid)})

  (->> (repeatedly 20 #(random-room 20 4 4))
       center))
