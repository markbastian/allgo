(ns allgo.geometry.tri-mesh
  "Triangle meshes: the surface elements cloth is simulated on, and the
  skin a solid is drawn with.

  A mesh is plain data, laid out flat for the same reasons as
  `allgo.geometry.tet-mesh`:

    {:verts    [x y z x y z ...]   ; 3 per vertex
     :tri-ids  [a b c ...]         ; 3 per triangle
     :edge-ids [a b ...]           ; 2 per unique edge
     :bend-ids [a b ...]}          ; 2 per pair of adjacent triangles

  Only `:verts` and `:tri-ids` carry information; `edges` and
  `bending-edges` derive the rest. `:bend-ids` is the interesting one: for
  every pair of triangles sharing an edge it names the two vertices
  *opposite* that edge. The distance between those two is what a fold
  changes, so holding it holds the fold -- which means bending needs no
  constraint type of its own, just a distance constraint over a different
  list. Cloth is then stretching and bending at two different compliances
  over the same solver."
  (:require [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

(defn vertex
  "`[x y z]` of vertex `i`."
  [verts i]
  (let [b (* 3 i)]
    [(nth verts b) (nth verts (+ b 1)) (nth verts (+ b 2))]))

(defn triangle-area [verts a b c]
  (let [p (vertex verts a)]
    (* 0.5 (math/sqrt (reduce + (map #(* % %)
                                     (v/cross (v/sub (vertex verts b) p)
                                              (v/sub (vertex verts c) p))))))))

(defn area
  "Total surface area."
  [{:keys [verts tri-ids]}]
  (reduce + (map (fn [[a b c]] (triangle-area verts a b c)) (partition 3 tri-ids))))

;; ---------------------------------------------------------------------------
;; Topology

(defn- half-edges
  "Every directed edge as `[{:key #{a b} :slot n :tri t :corner j}]`.

  A triangle's three edges are numbered `3t + j`, which is what lets a
  neighbour be recorded as a single number and the opposite vertex
  recovered from it."
  [tri-ids]
  (let [tris (vec tri-ids)]
    (for [t (range (quot (count tris) 3))
          j (range 3)
          :let [a (tris (+ (* 3 t) j))
                b (tris (+ (* 3 t) (mod (inc j) 3)))]]
      {:key #{a b} :slot (+ (* 3 t) j) :tri t :corner j})))

(defn triangle-neighbours
  "For each of the `3 * triangles` edge slots, the slot of the triangle
  across it, or -1 for an edge on the boundary."
  [tri-ids]
  (let [slots (half-edges tri-ids)
        n     (count slots)
        by-edge (group-by :key slots)]
    (reduce (fn [acc [_ shared]]
              (if (= 2 (count shared))
                (let [[x y] shared]
                  (-> acc (assoc (:slot x) (:slot y)) (assoc (:slot y) (:slot x))))
                acc))
            (vec (repeat n -1))
            by-edge)))

(defn edges
  "Every distinct edge, as a flat `[a b ...]`.

  An interior edge belongs to two triangles and must appear once, or the
  constraint over it is solved twice a pass and that edge comes out
  stiffer than its neighbours for no reason but connectivity."
  [tri-ids]
  (->> (half-edges tri-ids)
       (map :key)
       distinct
       (mapcat #(sort %))
       vec))

(defn bending-edges
  "For every pair of triangles sharing an edge, the two vertices opposite
  it, as a flat `[a b ...]`.

  Folding along the shared edge is exactly what moves those two vertices
  together or apart, so a distance constraint between them resists the
  fold. It is a cheap stand-in for a real dihedral-angle constraint, and
  the reason cloth needs no second constraint type."
  [tri-ids]
  (let [tris      (vec tri-ids)
        neighbour (triangle-neighbours tri-ids)
        opposite  (fn [slot]
                    ;; The corner a slot's edge does not touch.
                    (let [t (quot slot 3) j (mod slot 3)]
                      (tris (+ (* 3 t) (mod (+ j 2) 3)))))]
    (->> (range (count neighbour))
         (keep (fn [slot]
                 (let [other (neighbour slot)]
                   ;; Each shared edge is seen from both sides; keep one.
                   (when (and (>= other 0) (< slot other))
                     [(opposite slot) (opposite other)]))))
         (mapcat identity)
         vec)))

(defn boundary-edges
  "The edges belonging to only one triangle -- the outline of an open
  surface, and empty for a closed one."
  [tri-ids]
  (->> (half-edges tri-ids)
       (group-by :key)
       vals
       (filter #(= 1 (count %)))
       (mapcat #(sort (:key (first %))))
       vec))

(defn complete
  "Fills in `:edge-ids` and `:bend-ids` from `:verts` and `:tri-ids`."
  [{:keys [tri-ids] :as mesh}]
  (assoc mesh
         :edge-ids (edges tri-ids)
         :bend-ids (bending-edges tri-ids)))

;; ---------------------------------------------------------------------------
;; Generators

(defn grid
  "A `nx` by `nz` grid of cells in the xz plane at height `y`, each split
  into two triangles.

  The diagonal alternates between cells. A grid whose diagonals all run
  the same way is stiffer along one bias than the other and sags into a
  lopsided shape; alternating them keeps it even."
  ([nx nz] (grid nx nz 1.0 0.0))
  ([nx nz size y]
   (let [vid   (fn [i j] (+ i (* (inc nx) j)))
         verts (vec (for [j (range (inc nz))
                          i (range (inc nx))
                          v [(* size (- i (/ nx 2.0)))
                             y
                             (* size (- j (/ nz 2.0)))]]
                      (double v)))
         tris  (vec (for [j (range nz)
                          i (range nx)
                          id (if (even? (+ i j))
                               [(vid i j) (vid (inc i) j) (vid (inc i) (inc j))
                                (vid i j) (vid (inc i) (inc j)) (vid i (inc j))]
                               [(vid i j) (vid (inc i) j) (vid i (inc j))
                                (vid (inc i) j) (vid (inc i) (inc j)) (vid i (inc j))])]
                      id))]
     (complete {:verts verts :tri-ids tris}))))

(def ^:private icosahedron
  "The twelve vertices and twenty faces, wound counter-clockwise seen from
  outside. Three golden rectangles at right angles to each other, which is
  the neatest way to write the thing down."
  (let [t (/ (inc (math/sqrt 5.0)) 2.0)]
    {:verts [[-1.0 t 0.0] [1.0 t 0.0] [-1.0 (- t) 0.0] [1.0 (- t) 0.0]
             [0.0 -1.0 t] [0.0 1.0 t] [0.0 -1.0 (- t)] [0.0 1.0 (- t)]
             [t 0.0 -1.0] [t 0.0 1.0] [(- t) 0.0 -1.0] [(- t) 0.0 1.0]]
     :faces [[0 11 5] [0 5 1] [0 1 7] [0 7 10] [0 10 11]
             [1 5 9] [5 11 4] [11 10 2] [10 7 6] [7 1 8]
             [3 9 4] [3 4 2] [3 2 6] [3 6 8] [3 8 9]
             [4 9 5] [2 4 11] [6 2 10] [8 6 7] [9 8 1]]}))

(defn geodesic
  "A sphere of `radius` made by subdividing an icosahedron `n` times.

  `n` = 0 is the icosahedron itself; each level splits every triangle into
  four and pushes the three new vertices out onto the sphere, so the mesh
  has `10 * 4^n + 2` vertices and `20 * 4^n` faces.

  The reason to build one rather than to divide latitude and longitude:
  every triangle is nearly the same size and nearly equilateral, and there
  are no poles where the parameterisation piles up and no seam where it
  wraps. That matters for anything sampled per vertex -- a procedural
  planet spends the same effort per unit of surface everywhere, rather
  than lavishing it on two points nobody is looking at.

  Only `:verts` and `:tri-ids`; run it through `complete` if you want the
  edge lists, which a mesh that is only being drawn does not need."
  ([n] (geodesic n 1.0))
  ([n radius]
   (let [radius (double radius)
         unit (fn [[x y z]]
                (let [l (math/sqrt (+ (* x x) (* y y) (* z z)))]
                  [(/ x l) (/ y l) (/ z l)]))
         {:keys [verts faces]} icosahedron]
     (loop [level 0
            verts (mapv unit verts)
            faces faces]
       (if (>= level (long n))
         {:verts (vec (mapcat (fn [[x y z]] [(* radius x) (* radius y) (* radius z)]) verts))
          :tri-ids (vec (mapcat identity faces))}
         ;; Split every edge once, not once per face that uses it: the
         ;; cache is what keeps the mesh welded, and a mesh that is not
         ;; welded has seams down every original edge.
         (let [[verts' mids]
               (reduce (fn [[vs cache] [a b]]
                         (if (contains? cache [a b])
                           [vs cache]
                           (let [i (count vs)
                                 m (unit (mapv + (nth vs a) (nth vs b)))]
                             [(conj vs m) (assoc cache [a b] i [b a] i)])))
                       [verts {}]
                       (mapcat (fn [[a b c]] [[(min a b) (max a b)]
                                              [(min b c) (max b c)]
                                              [(min c a) (max c a)]])
                               faces))
               faces' (vec (mapcat (fn [[a b c]]
                                     (let [ab (mids [a b]) bc (mids [b c]) ca (mids [c a])]
                                       [[a ab ca] [b bc ab] [c ca bc] [ab bc ca]]))
                                   faces))]
           (recur (inc level) verts' faces')))))))

(defn mean-edge-length
  "The average length of a triangle side.

  How finely the mesh samples whatever is being drawn on it, which is
  what anything procedural needs to know before it decides how much
  detail to generate -- see `allgo.procedural.fractal/octaves-for`.
  Interior edges are counted once per triangle that uses them, which
  cannot move the mean of a mesh whose triangles are all much the same
  size, and this is only ever asked of meshes like that."
  ^double [{:keys [verts tri-ids]}]
  (let [tris (partition 3 tri-ids)
        d (fn ^double [a b]
            (let [p (vertex verts a) q (vertex verts b)]
              (v/length (v/sub p q))))
        total (reduce (fn [^double acc [a b c]]
                        (+ acc (d a b) (d b c) (d c a)))
                      0.0
                      tris)]
    (if (zero? (count tris)) 0.0 (/ total (* 3.0 (count tris))))))

(defn translate [mesh [dx dy dz]]
  (update mesh :verts
          (fn [verts]
            (vec (map-indexed (fn [i v]
                                (+ v (let [axis (long (mod i 3))]
                                       (cond (zero? axis) dx (= 1 axis) dy :else dz))))
                              verts)))))

(defn deform
  "Moves every vertex through `f`, a function of `[x y z]`."
  [mesh f]
  (assoc mesh :verts (vec (mapcat f (partition 3 (:verts mesh))))))

(defn bounds
  "`[[min-x min-y min-z] [max-x max-y max-z]]`."
  [{:keys [verts]}]
  (let [pts (partition 3 verts)]
    [(mapv #(apply min (map (fn [p] (nth p %)) pts)) [0 1 2])
     (mapv #(apply max (map (fn [p] (nth p %)) pts)) [0 1 2])]))

(defn corner-vertices
  "The indices of the four extreme corners in the xz plane, which is what
  a hanging sheet is usually pinned by."
  [{:keys [verts] :as mesh}]
  (let [[lo hi] (bounds mesh)
        pts     (vec (partition 3 verts))
        nearest (fn [tx tz]
                  (apply min-key
                         (fn [i] (let [[x _ z] (pts i)]
                                   (+ (* (- x tx) (- x tx)) (* (- z tz) (- z tz)))))
                         (range (count pts))))]
    [(nearest (lo 0) (lo 2)) (nearest (hi 0) (lo 2))
     (nearest (lo 0) (hi 2)) (nearest (hi 0) (hi 2))]))
