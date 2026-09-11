(ns allgo.geometry.voronoi3d
  "Voronoi diagrams in three dimensions, by half-space clipping.

  Each cell starts as the bounding box and is cut by the perpendicular
  bisector plane between its site and every other site. What survives is the
  set of points closer to this site than to any other -- the Voronoi cell,
  bounded and convex by construction.

  This is deliberately not the dual of a Delaunay tetrahedralization. That
  route needs a convex hull per cell (the 2D angular-sort shortcut has no
  3D analogue), and 3D Bowyer-Watson is fragile near cospherical inputs:
  roundoff corrupts the cavity so boundary extraction yields a non-manifold
  mesh. Clipping degrades gracefully instead -- each cut is independent, so
  error cannot cascade into bad topology. Delaunay adjacency is still
  recovered, since two sites neighbour exactly when their cells share a face.

  A cell is `{:faces [[[x y z] ...] ...] :neighbours #{site ...}}`."
  (:require [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

(def ^:private eps 1e-9)

(defn- centroid [pts]
  (mapv #(/ % (count pts)) (reduce (fn [a b] (mapv + a b)) pts)))

;; A plane is {:n normal :d offset}; a point is inside when (v/dot n p) <= d.

(defn- bisector
  "The half-space of points at least as close to `p` as to `q`."
  [p q]
  {:n (v/sub q p) :d (/ (- (v/length-squared q) (v/length-squared p)) 2.0)})

(defn- crossing [a b sa sb]
  (let [t (/ sa (- sa sb))]
    (mapv (fn [x y] (+ x (* t (- y x)))) a b)))

(defn- clip-face
  "Clip convex polygon `pts` to the inside of the plane `n.x = d`, returning
  `[surviving-loop points-on-the-cut]`. The plane arrives as primitives and
  offsets are computed in one pass, so the loop does no map lookups and one
  dot product per vertex rather than two."
  [pts nx ny nz d]
  (let [n    (count pts)
        offs (mapv (fn [[x y z]] (- (+ (* nx x) (* ny y) (* nz z)) d)) pts)]
    (loop [i 0, out (transient []), cuts (transient [])]
      (if (= i n)
        [(persistent! out) (persistent! cuts)]
        (let [j  (rem (inc i) n)
              sa (nth offs i)
              sb (nth offs j)
              a  (nth pts i)
              b  (nth pts j)
              ia (<= sa eps)
              ib (<= sb eps)]
          (cond
            (and ia ib) (recur (inc i) (conj! out b) cuts)
            ia          (let [x (crossing a b sa sb)]
                          (recur (inc i) (conj! out x) (conj! cuts x)))
            ib          (let [x (crossing a b sa sb)]
                          (recur (inc i) (-> out (conj! x) (conj! b)) (conj! cuts x)))
            :else       (recur (inc i) out cuts)))))))

(defn- distinct-points [pts]
  (reduce (fn [acc [x y z :as p]]
            (if (some (fn [[ax ay az]]
                        (let [dx (- x ax) dy (- y ay) dz (- z az)]
                          (< (+ (* dx dx) (* dy dy) (* dz dz)) 1e-14)))
                      acc)
              acc
              (conj acc p)))
          []
          pts))

(defn- cap-loop
  "The new face sealing a cut: the cut points wound in order around their
  centroid, in a basis spanning the cutting plane. Angles are precomputed --
  sort-by calls its keyfn on every comparison, which would mean several
  times more atan2 calls than there are points."
  [cuts normal]
  (let [pts (distinct-points cuts)]
    (when (>= (count pts) 3)
      (let [c (centroid pts)
            u (v/normalize (v/sub (first pts) c))
            v (v/cross (v/normalize normal) u)]
        (->> pts
             (mapv (fn [p] (let [w (v/sub p c)] [(math/atan2 (v/dot w v) (v/dot w u)) p])))
             (sort-by first)
             (mapv second))))))

;; A face carries the site whose bisector created it (nil for the bounding
;; box), so the surviving faces name the cell's true Voronoi neighbours.

(defn separating-plane?
  "True when the plane `n.x = d` is a separating axis for the convex
  polyhedron `faces`: every vertex lies inside it, so it cannot cut.

  Exact, and tighter than a radius bound, which being omnidirectional has to
  assume the worst direction. Offered because that is a genuinely useful
  predicate, but deliberately not used by `cell`: measured in the clip loop
  it is a wash, because the radius bound already rejects distant sites, so
  almost everything that reaches a clip really does cut and the pre-test
  only adds a pass. Worth reaching for when testing planes against a cell
  you are *not* about to clip anyway."
  [faces nx ny nz d]
  (every? (fn [f]
            (every? (fn [[x y z]] (<= (- (+ (* nx x) (* ny y) (* nz z)) d) eps))
                    (:pts f)))
          faces))

(defn- clip-cell
  "Cut convex polyhedron `faces` by `plane`. Returns nil when the plane does
  not reach the cell, so callers can tell a real cut from a no-op."
  [faces {[nx ny nz] :n d :d} owner]
  (let [results (mapv (fn [f] (clip-face (:pts f) nx ny nz d)) faces)
        cuts    (into [] (mapcat second) results)]
    (when-let [cap (cap-loop cuts [nx ny nz])]
      (-> (into [] (keep-indexed (fn [i [loop* _]]
                                   (when (>= (count loop*) 3)
                                     (assoc (nth faces i) :pts loop*))))
                results)
          (conj {:pts cap :site owner})))))

(defn box-faces
  "The six faces of `[[x0 y0 z0] [x1 y1 z1]]`, each wound consistently."
  [[[x0 y0 z0] [x1 y1 z1]]]
  [[[x0 y0 z0] [x0 y1 z0] [x0 y1 z1] [x0 y0 z1]]
   [[x1 y0 z0] [x1 y0 z1] [x1 y1 z1] [x1 y1 z0]]
   [[x0 y0 z0] [x0 y0 z1] [x1 y0 z1] [x1 y0 z0]]
   [[x0 y1 z0] [x1 y1 z0] [x1 y1 z1] [x0 y1 z1]]
   [[x0 y0 z0] [x1 y0 z0] [x1 y1 z0] [x0 y1 z0]]
   [[x0 y0 z1] [x0 y1 z1] [x1 y1 z1] [x1 y0 z1]]])

(defn- farthest-sq [[sx sy sz] faces]
  (reduce (fn [m f]
            (reduce (fn [m [x y z]]
                      (let [dx (- x sx) dy (- y sy) dz (- z sz)]
                        (max m (+ (* dx dx) (* dy dy) (* dz dz)))))
                    m (:pts f)))
          0.0
          faces))

(defn- finish [faces]
  {:faces (mapv :pts faces) :neighbours (into #{} (keep :site) faces)})

(defn cell
  "The Voronoi cell of `site` against `others`, clipped to `bounds`.

  `others` is walked nearest-first so the loop can stop early: a bisector
  sits half-way to its site, so once that half-distance exceeds the cell's
  farthest vertex, neither it nor anything beyond can cut the cell.

  `seed` is an optional collection of sites believed to be neighbours -- the
  previous frame's, typically. Clipping those first collapses the radius
  bound immediately, so the walk terminates far sooner. A wrong seed costs
  a little work but cannot give a wrong answer: the ordered walk still
  considers every site the bound does not exclude."
  ([site others bounds] (cell site others bounds nil))
  ([site others bounds seed]
   (let [start  (mapv (fn [pts] {:pts pts :site nil}) (box-faces bounds))
         seeded (reduce (fn [fs o]
                          (or (clip-cell fs (bisector site o) o) fs))
                        start
                        seed)
         known   (set seed)
         ordered (->> others
                      (mapv (fn [o] [(v/length-squared (v/sub o site)) o]))
                      (sort-by first)
                      (mapv second))]
     (loop [faces     seeded
            r2        (farthest-sq site seeded)
            remaining ordered]
       (if-let [other (first remaining)]
         (if (> (v/length-squared (v/sub other site)) (* 4.0 r2))
           (finish faces)
           (if (contains? known other)
             (recur faces r2 (rest remaining))
             (if-let [cut (clip-cell faces (bisector site other) other)]
               (recur cut (farthest-sq site cut) (rest remaining))
               (recur faces r2 (rest remaining)))))
         (finish faces))))))

(defn- build
  [sites bounds hints indexed?]
  (let [sites (vec sites)
        n     (count sites)
        index (when indexed?
                (persistent! (reduce (fn [m i] (assoc! m (nth sites i) i))
                                     (transient {}) (range n))))]
    (persistent!
     (reduce (fn [acc i]
               (let [site (nth sites i)
                     seed (when hints
                            (into [] (comp (keep #(nth sites % nil)) (remove #{site}))
                                  (get hints i)))
                     c    (cell site
                                (concat (subvec sites 0 i) (subvec sites (inc i)))
                                bounds
                                seed)]
                 (assoc! acc site
                         (cond-> c
                           indexed? (assoc :neighbour-idx
                                           (into #{} (keep index) (:neighbours c)))))))
             (transient {})
             (range n)))))

(defn diagram
  "Voronoi diagram of `sites` clipped to `bounds` (`[[x0 y0 z0] [x1 y1 z1]]`):
  a map from each site to its cell. The cells tile `bounds` exactly.

  The three-argument form additionally tags each cell with `:neighbour-idx`,
  its neighbours as indices into `sites`, and accepts the previous frame's
  indices (see `neighbour-hints`, and pass nil on the first frame) so each
  cell starts from its likely neighbours instead of rediscovering them.
  Positions move every frame, so an index is the only stable handle.

  Seeding is measured as a wash below a few hundred sites and is off in the
  two-argument form for that reason. The early exit has to clear every site
  within twice a cell's reach, which is around 58 sites whatever the total
  -- so under a few hundred there is no locality for a seed to exploit, and
  clipping known neighbours first only reorders work rather than avoiding
  it. Past that point the constant stops dominating and seeding pays."
  ([sites bounds] (build sites bounds nil false))
  ([sites bounds hints] (build sites bounds hints true)))

(defn neighbour-hints
  "Per-site neighbour indices from `diagram`, in `sites` order, ready to pass
  back as the `hints` argument on the next frame. Requires a diagram built
  with the three-argument form."
  [diagram sites]
  (mapv #(:neighbour-idx (get diagram %)) sites))

(defn edges
  "Delaunay edges implied by `diagram`: the unordered site pairs whose cells
  share a face."
  [diagram]
  (into #{} (for [[site {:keys [neighbours]}] diagram
                  other neighbours]
              #{site other})))
