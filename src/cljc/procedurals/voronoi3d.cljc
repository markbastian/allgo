(ns procedurals.voronoi3d
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
  (:require [clojure.math :as math]))

(def ^:private eps 1e-9)

(defn- v- [[ax ay az] [bx by bz]] [(- ax bx) (- ay by) (- az bz)])
(defn- dot [[ax ay az] [bx by bz]] (+ (* ax bx) (* ay by) (* az bz)))
(defn- cross [[ax ay az] [bx by bz]]
  [(- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))])

(defn- norm-sq [v] (dot v v))

(defn- normalize [v]
  (let [m (math/sqrt (norm-sq v))]
    (if (zero? m) v (mapv #(/ % m) v))))

(defn- centroid [pts]
  (mapv #(/ % (count pts)) (reduce (fn [a b] (mapv + a b)) pts)))

;; A plane is {:n normal :d offset}; a point is inside when (dot n p) <= d.

(defn- bisector
  "The half-space of points at least as close to `p` as to `q`."
  [p q]
  {:n (v- q p) :d (/ (- (norm-sq q) (norm-sq p)) 2.0)})

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
            u (normalize (v- (first pts) c))
            v (cross (normalize normal) u)]
        (->> pts
             (mapv (fn [p] (let [w (v- p c)] [(math/atan2 (dot w v) (dot w u)) p])))
             (sort-by first)
             (mapv second))))))

;; A face carries the site whose bisector created it (nil for the bounding
;; box), so the surviving faces name the cell's true Voronoi neighbours.

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

(defn cell
  "The Voronoi cell of `site` against `others`, clipped to `bounds`.

  `others` is walked nearest-first so the loop can stop early: a bisector
  sits half-way to its site, so once that half-distance exceeds the cell's
  farthest vertex, neither it nor anything beyond it can cut the cell. This
  is what keeps the cost per cell near-constant rather than O(n)."
  [site others bounds]
  ;; Keys precomputed: sort-by calls its keyfn on every comparison, and each
  ;; call here would allocate a difference vector.
  (let [ordered (->> others
                     (mapv (fn [o] [(norm-sq (v- o site)) o]))
                     (sort-by first)
                     (mapv second))]
    (loop [faces     (mapv (fn [pts] {:pts pts :site nil}) (box-faces bounds))
           r2        (farthest-sq site (mapv (fn [pts] {:pts pts}) (box-faces bounds)))
           remaining ordered]
      (if-let [other (first remaining)]
        (if (> (norm-sq (v- other site)) (* 4.0 r2))
          {:faces (mapv :pts faces) :neighbours (into #{} (keep :site) faces)}
          (if-let [cut (clip-cell faces (bisector site other) other)]
            (recur cut (farthest-sq site cut) (rest remaining))
            (recur faces r2 (rest remaining))))
        {:faces (mapv :pts faces) :neighbours (into #{} (keep :site) faces)}))))

(defn diagram
  "Voronoi diagram of `sites` clipped to `bounds` (`[[x0 y0 z0] [x1 y1 z1]]`):
  a map from each site to its cell. The cells tile `bounds` exactly."
  [sites bounds]
  (let [sites (vec sites)]
    (persistent!
     (reduce (fn [acc i]
               (assoc! acc (nth sites i)
                       (cell (nth sites i)
                             (concat (subvec sites 0 i) (subvec sites (inc i)))
                             bounds)))
             (transient {})
             (range (count sites))))))

(defn edges
  "Delaunay edges implied by `diagram`: the unordered site pairs whose cells
  share a face."
  [diagram]
  (into #{} (for [[site {:keys [neighbours]}] diagram
                  other neighbours]
              #{site other})))
