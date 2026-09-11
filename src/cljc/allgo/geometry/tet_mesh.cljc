(ns allgo.geometry.tet-mesh
  "Tetrahedral meshes: the volume elements a soft body is simulated on.

  A mesh is plain data,

    {:verts           [x y z x y z ...]   ; 3 per vertex
     :tet-ids         [a b c d ...]       ; 4 per tetrahedron
     :edge-ids        [a b ...]           ; 2 per unique edge
     :surface-tri-ids [a b c ...]}        ; 3 per boundary triangle

  flat rather than nested because that is what both the solver and a
  graphics buffer want, and converting between the two forms at every
  frame would cost more than the simulation.

  Only `:verts` and `:tet-ids` carry information: `edges` and
  `surface-triangles` derive the rest, so any source of tetrahedra gets
  the edges its distance constraints need and the skin to draw for free."
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Geometry

(defn vertex
  "`[x y z]` of vertex `i`."
  [verts i]
  (let [b (* 3 i)]
    [(nth verts b) (nth verts (+ b 1)) (nth verts (+ b 2))]))

(defn- v- [[ax ay az] [bx by bz]] [(- ax bx) (- ay by) (- az bz)])

(defn- cross [[ax ay az] [bx by bz]]
  [(- (* ay bz) (* az by))
   (- (* az bx) (* ax bz))
   (- (* ax by) (* ay bx))])

(defn- dot [[ax ay az] [bx by bz]] (+ (* ax bx) (* ay by) (* az bz)))

(defn tet-volume
  "Signed volume of the tetrahedron on vertices `a b c d`.

  Signed, and the sign is load-bearing: it is positive exactly when the
  vertices are wound consistently, so it doubles as the orientation test
  used to fix up a tetrahedron built the wrong way round."
  [verts a b c d]
  (let [p0 (vertex verts a)]
    (/ (dot (cross (v- (vertex verts b) p0) (v- (vertex verts c) p0))
            (v- (vertex verts d) p0))
       6.0)))

(defn- oriented
  "`[a b c d]` wound so the volume is positive."
  [verts [a b c d]]
  (if (neg? (tet-volume verts a b c d)) [b a c d] [a b c d]))

;; ---------------------------------------------------------------------------
;; Derived topology

(def ^:private tet-edge-pairs [[0 1] [0 2] [0 3] [1 2] [1 3] [2 3]])

(defn edges
  "Every distinct edge of the mesh, as a flat `[a b ...]`.

  Interior edges are shared by many tetrahedra and must appear once:
  solving the same distance constraint several times in a pass would make
  those edges stiffer than the rest purely as an artefact of connectivity."
  [tet-ids]
  (->> (partition 4 tet-ids)
       (mapcat (fn [tet]
                 (let [tet (vec tet)]
                   (map (fn [[i j]] (sort [(tet i) (tet j)])) tet-edge-pairs))))
       distinct
       (mapcat identity)
       vec))

(def ^:private tet-face-order
  "The four faces of a tetrahedron, each wound to face away from the
  vertex it omits."
  [[1 3 2] [0 2 3] [0 3 1] [0 1 2]])

(defn surface-triangles
  "The boundary of the mesh, as a flat `[a b c ...]` wound outward.

  A face between two tetrahedra belongs to both; a face on the skin
  belongs to one. So counting how often each face appears -- as a set of
  vertices, since neighbours meet it with opposite winding -- separates
  the surface from the interior without any geometric test."
  [tet-ids]
  (let [faces (for [tet  (partition 4 tet-ids)
                    :let [tet (vec tet)]
                    idx  tet-face-order]
                (mapv tet idx))]
    (->> faces
         (group-by set)
         vals
         (filter #(= 1 (count %)))
         (mapcat first)
         vec)))

(defn complete
  "Fills in `:edge-ids` and `:surface-tri-ids` from `:verts` and
  `:tet-ids`."
  [{:keys [tet-ids] :as mesh}]
  (assoc mesh
         :edge-ids (edges tet-ids)
         :surface-tri-ids (surface-triangles tet-ids)))

;; ---------------------------------------------------------------------------
;; Generators

(def ^:private kuhn-paths
  "The six corner-walks that cut a cube into tetrahedra.

  Each walks from corner 000 to corner 111 changing one axis at a time,
  and the six orderings of the three axes give six tetrahedra of equal
  volume. This is Kuhn's subdivision, chosen because it is *conforming*:
  neighbouring cubes cut their shared face the same way, so the mesh has
  no cracks. Cutting a cube into five tetrahedra is more economical but
  only conforms if alternate cubes are mirrored."
  (for [[a b c] [[0 1 2] [0 2 1] [1 0 2] [1 2 0] [2 0 1] [2 1 0]]]
    (reductions (fn [corner axis] (update corner axis inc)) [0 0 0] [a b c])))

(defn lattice-box
  "A box of `nx` by `ny` by `nz` cells, each cut into six tetrahedra.

  `size` is the edge length of one cell; the box is centred on the origin
  in x and z and sits with its base at y = 0, which is where a soft body
  dropped on a floor wants to start."
  ([nx ny nz] (lattice-box nx ny nz 1.0))
  ([nx ny nz size]
   (let [vid   (fn [i j k] (+ i (* (inc nx) (+ j (* (inc ny) k)))))
         verts (vec (for [k (range (inc nz))
                          j (range (inc ny))
                          i (range (inc nx))
                          v [(* size (- i (/ nx 2.0)))
                             (* size j)
                             (* size (- k (/ nz 2.0)))]]
                      (double v)))
         tets  (vec (for [k    (range nz)
                          j    (range ny)
                          i    (range nx)
                          path kuhn-paths
                          :let [ids (mapv (fn [[di dj dk]] (vid (+ i di) (+ j dj) (+ k dk)))
                                          path)]
                          id   (oriented verts ids)]
                      id))]
     (complete {:verts verts :tet-ids tets}))))

(defn translate
  "Shifts every vertex."
  [mesh [dx dy dz]]
  (update mesh :verts
          (fn [verts]
            (vec (map-indexed (fn [i v] (+ v (case (mod i 3) 0 dx 1 dy 2 dz)))
                              verts)))))

(defn deform
  "Moves every vertex through `f`, a function of `[x y z]`.

  Tetrahedra can invert under a deformation, so orientation is restored
  afterwards -- a tetrahedron with negative volume would be inside out,
  and the solver would work to keep it that way."
  [mesh f]
  (let [verts (vec (mapcat f (partition 3 (:verts mesh))))]
    (assoc mesh
           :verts verts
           :tet-ids (vec (mapcat #(oriented verts (vec %)) (partition 4 (:tet-ids mesh)))))))

(defn volume
  "Total volume of the mesh."
  [{:keys [verts tet-ids]}]
  (reduce + (map (fn [[a b c d]] (tet-volume verts a b c d)) (partition 4 tet-ids))))

(defn bounds
  "`[[min-x min-y min-z] [max-x max-y max-z]]`."
  [{:keys [verts]}]
  (let [pts (partition 3 verts)]
    [(mapv #(apply min (map (fn [p] (nth p %)) pts)) [0 1 2])
     (mapv #(apply max (map (fn [p] (nth p %)) pts)) [0 1 2])]))

(defn sphere-deformation
  "A deformation that pushes the corners of a box out onto a ball, for a
  rounder body than a lattice gives on its own."
  [center radius]
  (let [[cx cy cz] center]
    (fn [[x y z]]
      (let [[dx dy dz] [(- x cx) (- y cy) (- z cz)]
            len        (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz)))]
        (if (zero? len)
          [x y z]
          ;; Chebyshev distance to Euclidean: the box's own surface maps
          ;; onto the sphere and the interior is carried along with it.
          (let [box   (max (abs dx) (abs dy) (abs dz))
                scale (/ (* radius box) (* len (max radius 1e-9)))]
            [(+ cx (* dx scale)) (+ cy (* dy scale)) (+ cz (* dz scale))]))))))
