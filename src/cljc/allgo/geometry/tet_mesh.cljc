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
  (:require [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Geometry

(defn vertex
  "`[x y z]` of vertex `i`."
  [verts i]
  (let [b (* 3 i)]
    [(nth verts b) (nth verts (+ b 1)) (nth verts (+ b 2))]))

(defn tet-volume
  "Signed volume of the tetrahedron on vertices `a b c d`.

  Signed, and the sign is load-bearing: it is positive exactly when the
  vertices are wound consistently, so it doubles as the orientation test
  used to fix up a tetrahedron built the wrong way round."
  [verts a b c d]
  (let [p0 (vertex verts a)]
    (/ (v/dot (v/cross (v/sub (vertex verts b) p0) (v/sub (vertex verts c) p0))
              (v/sub (vertex verts d) p0))
       6.0)))

(defn- oriented
  "`[a b c d]` wound so the volume is positive."
  [verts [a b c d]]
  (if (neg? (tet-volume verts a b c d)) [b a c d] [a b c d]))

(defn barycentric
  "Where point `p` sits relative to the tetrahedron `a b c d`, as the four
  weights that reconstruct it: `p = b0*a + b1*b + b2*c + b3*d`, summing to
  one.

  All four non-negative means the point is inside. A negative weight says
  which face it is outside of, and how far in units of the tetrahedron's
  own size -- which is what lets a point that no tetrahedron contains
  still be attached to the one it is least outside.

  Taking `d` as the origin, the other three edges form a matrix whose
  inverse carries `p - d` into the first three weights; the fourth is
  whatever is left of one. Returns nil for a degenerate tetrahedron,
  whose matrix has no inverse."
  [verts a b c d [px py pz]]
  (let [[dx dy dz] (vertex verts d)
        [a11 a21 a31] (v/sub (vertex verts a) [dx dy dz])
        [a12 a22 a32] (v/sub (vertex verts b) [dx dy dz])
        [a13 a23 a33] (v/sub (vertex verts c) [dx dy dz])
        det (+ (* a11 (- (* a22 a33) (* a23 a32)))
               (- (* a12 (- (* a21 a33) (* a23 a31))))
               (* a13 (- (* a21 a32) (* a22 a31))))]
    (when-not (zero? det)
      (let [inv (/ 1.0 det)
            rx (- px dx) ry (- py dy) rz (- pz dz)
            b0 (* inv (+ (* (- (* a22 a33) (* a23 a32)) rx)
                         (* (- (* a13 a32) (* a12 a33)) ry)
                         (* (- (* a12 a23) (* a13 a22)) rz)))
            b1 (* inv (+ (* (- (* a23 a31) (* a21 a33)) rx)
                         (* (- (* a11 a33) (* a13 a31)) ry)
                         (* (- (* a13 a21) (* a11 a23)) rz)))
            b2 (* inv (+ (* (- (* a21 a32) (* a22 a31)) rx)
                         (* (- (* a12 a31) (* a11 a32)) ry)
                         (* (- (* a11 a22) (* a12 a21)) rz)))]
        [b0 b1 b2 (- 1.0 b0 b1 b2)]))))

(defn inside?
  "Whether `p` is within the tetrahedron, to a tolerance in barycentric
  units."
  ([verts a b c d p] (inside? verts a b c d p 0.0))
  ([verts a b c d p tolerance]
   (when-let [bary (barycentric verts a b c d p)]
     (every? #(>= % (- tolerance)) bary))))

;; ---------------------------------------------------------------------------
;; Derived topology

(def ^:private tet-edge-pairs [[0 1] [0 2] [0 3] [1 2] [1 3] [2 3]])

(defn edges
  "Every distinct edge of the mesh, as a flat `[a b ...]`.

  Interior edges are shared by many tetrahedra and must appear once:
  solving the same distance constraint several times in a pass would make
  those edges stiffer than the rest purely as an artifact of connectivity."
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
  vertices, since neighbors meet it with opposite winding -- separates
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
  neighboring cubes cut their shared face the same way, so the mesh has
  no cracks. Cutting a cube into five tetrahedra is more economical but
  only conforms if alternate cubes are mirrored."
  (for [[a b c] [[0 1 2] [0 2 1] [1 0 2] [1 2 0] [2 0 1] [2 1 0]]]
    (reductions (fn [corner axis] (update corner axis inc)) [0 0 0] [a b c])))

(defn lattice-box
  "A box of `nx` by `ny` by `nz` cells, each cut into six tetrahedra.

  `size` is the edge length of one cell; the box is centered on the origin
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
            ;; `case` on a boxed index compiles to a hash lookup; the
            ;; three-way `cond` on longs does not.
            (vec (map-indexed (fn [i v]
                                (+ v (let [axis (long (mod i 3))]
                                       (cond (zero? axis) dx (= 1 axis) dy :else dz))))
                              verts)))))

(defn deform
  "Moves every vertex through `f`, a function of `[x y z]`.

  Tetrahedra can invert under a deformation, so orientation is restored
  afterward -- a tetrahedron with negative volume would be inside out,
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
