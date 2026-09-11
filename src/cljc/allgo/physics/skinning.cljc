(ns allgo.physics.skinning
  "Driving a detailed visual mesh from a coarse simulation mesh, after Ten
  Minute Physics 12.

  A soft body's cost is in its tetrahedra, and a surface detailed enough to
  look like anything has far more vertices than a body needs elements. So
  simulate a coarse cage and carry the detailed surface along inside it:
  each visual vertex is bound once to the tetrahedron containing it, as the
  four weights that reconstruct its position from that tetrahedron's
  corners, and thereafter follows for four multiplies and three adds.

  The saving is in the ratio. A cage of a few hundred tetrahedra runs at
  the same speed whether the surface it carries has a thousand vertices or
  fifty thousand, because the surface is not simulated at all -- which is
  where the tutorial's hundredfold comes from.

  Two details make it work on real meshes:

  Binding runs tetrahedron-first, not vertex-first. Each tetrahedron
  queries an `allgo.spatial.hash` built over the visual vertices, rather
  than each vertex searching for a tetrahedron. Both find the same pairs,
  but this way the matrix inverse that turns a position into weights is
  computed once per tetrahedron instead of once per test.

  A vertex no tetrahedron contains is still bound. A coarse cage does not
  cover a detailed surface exactly, and vertices near the skin fall
  outside. The weights for such a point come out negative, and the size of
  the most negative one says how far outside it is in units of the
  tetrahedron's own extent -- so keeping the tetrahedron that minimises it
  attaches every vertex to the element it is least outside of, and the
  extrapolation carries it along just as well."
  (:require [allgo.geometry.tet-mesh :as tet]
            [allgo.spatial.hash :as spatial]
            [clojure.math :as math]))

(defn- f64
  ([n] #?(:clj (double-array n) :cljs (js/Float64Array. n)))
  ([_n coll] #?(:clj (double-array (map double coll))
                :cljs (js/Float64Array. (into-array (map double coll))))))

(defn- i32 [n] #?(:clj (int-array n) :cljs (js/Int32Array. n)))

;; ---------------------------------------------------------------------------

(defn- tet-radius
  "Centre and bounding radius of one tetrahedron, as `[cx cy cz r]`."
  [^doubles pos ids base]
  (let [ids (vec ids)
        cx (/ (reduce + (map #(aget pos (* 3 (ids (+ base %)))) (range 4))) 4.0)
        cy (/ (reduce + (map #(aget pos (inc (* 3 (ids (+ base %))))) (range 4))) 4.0)
        cz (/ (reduce + (map #(aget pos (+ 2 (* 3 (ids (+ base %))))) (range 4))) 4.0)
        r  (reduce max 0.0
                   (for [j (range 4)
                         :let [p (* 3 (ids (+ base j)))
                               dx (- (aget pos p) cx)
                               dy (- (aget pos (+ p 1)) cy)
                               dz (- (aget pos (+ p 2)) cz)]]
                     (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz)))))]
    [cx cy cz r]))

(defn bind
  "Binds every vertex of `vis-verts` to a tetrahedron of `mesh`.

  Returns `{:n :tet :bary :bound}` -- one tetrahedron index and three
  weights per visual vertex, with the fourth weight implied by the three
  summing to one. `:bound` counts the vertices that landed strictly
  inside; the rest are attached to the nearest tetrahedron and carried by
  extrapolation.

  `:border` widens each tetrahedron's search by a fraction of its own
  radius, so vertices just outside the cage are still found. Binding is
  done once, against the rest pose."
  ([mesh vis-verts] (bind mesh vis-verts {}))
  ([{:keys [verts tet-ids]} vis-verts {:keys [border] :or {border 0.25}}]
   (let [^doubles vis (f64 (count vis-verts) vis-verts)
         ^doubles pos (f64 (count verts) verts)
         ids      (vec tet-ids)
         n-vis    (quot (alength vis) 3)
         n-tets   (quot (count ids) 4)
         radii    (mapv #(tet-radius pos ids (* 4 %)) (range n-tets))
         ;; Cells the size of a typical element: a tetrahedron's query then
         ;; sweeps a couple of cells rather than hundreds or one.
         spacing  (max 1e-9 (/ (reduce + (map peek radii)) (max 1 n-tets)))
         hash     (spatial/spatial-hash spacing (max 1 n-vis))
         ^ints tet-of (i32 n-vis)
         ^doubles bary (f64 (* 3 n-vis))
         ;; How far outside its tetrahedron each vertex is, in barycentric
         ;; units. Zero means inside, and once a vertex is inside one
         ;; tetrahedron there is nothing better to find.
         ^doubles outside (f64 n-vis)]
     (dotimes [i n-vis] (aset tet-of i -1))
     (dotimes [i n-vis] (aset outside i ##Inf))
     (spatial/rebuild! hash vis n-vis)

     (dotimes [t n-tets]
       (let [base (* 4 t)
             [cx cy cz r] (radii t)
             reach (* r (+ 1.0 border))
             found (spatial/query-point! hash cx cy cz reach)]
         (when (pos? found)
           (let [a (ids base) b (ids (+ base 1)) c (ids (+ base 2)) d (ids (+ base 3))]
             (dotimes [k found]
               (let [id (spatial/neighbour hash k)]
                 ;; Already inside something: nothing can beat that.
                 (when (pos? (aget outside id))
                   (let [p  (* 3 id)
                         vx (aget vis p) vy (aget vis (+ p 1)) vz (aget vis (+ p 2))
                         dx (- vx cx) dy (- vy cy) dz (- vz cz)]
                     (when (<= (+ (* dx dx) (* dy dy) (* dz dz)) (* reach reach))
                       (when-let [weights (tet/barycentric verts a b c d [vx vy vz])]
                         ;; Out of a vector, so boxed; coerced once here
                         ;; rather than reflecting on every store.
                         (let [b0 (double (nth weights 0))
                               b1 (double (nth weights 1))
                               b2 (double (nth weights 2))
                               b3 (double (nth weights 3))
                               worst (max 0.0 (- b0) (- b1) (- b2) (- b3))]
                           (when (< worst (aget outside id))
                             (aset outside id worst)
                             (aset tet-of id t)
                             (aset bary (* 3 id) b0)
                             (aset bary (inc (* 3 id)) b1)
                             (aset bary (+ 2 (* 3 id)) b2)))))))))))))

     ;; Anything the sweep did not reach is bound to the tetrahedron whose
     ;; centre is nearest, whatever the distance. A surface much larger
     ;; than its cage has vertices no element comes close to, and leaving
     ;; those unbound strands them at the origin -- which reads as the mesh
     ;; tearing. Extrapolating from a distant element is a stretch, but it
     ;; is a stretch that moves with the body. There are usually a handful,
     ;; so the linear scan over elements costs nothing.
     (when (pos? n-tets)
       (dotimes [i n-vis]
         (when (neg? (aget tet-of i))
           (let [p  (* 3 i)
                 vx (aget vis p) vy (aget vis (+ p 1)) vz (aget vis (+ p 2))
                 nearest (reduce (fn [[_ best-d :as best] t]
                                   (let [[cx cy cz _] (radii t)
                                         dx (- vx cx) dy (- vy cy) dz (- vz cz)
                                         d  (+ (* dx dx) (* dy dy) (* dz dz))]
                                     (if (< d best-d) [t d] best)))
                                 [-1 ##Inf]
                                 (range n-tets))
                 ;; Out of a vector, so boxed.
                 t (long (first nearest))
                 base (* 4 t)]
             (when-let [weights (tet/barycentric verts (ids base) (ids (+ base 1))
                                                 (ids (+ base 2)) (ids (+ base 3))
                                                 [vx vy vz])]
               (aset tet-of i (int t))
               (aset bary (* 3 i) (double (nth weights 0)))
               (aset bary (inc (* 3 i)) (double (nth weights 1)))
               (aset bary (+ 2 (* 3 i)) (double (nth weights 2))))))))

     {:n     n-vis
      :tet   tet-of
      :bary  bary
      ;; Vertices that landed strictly inside some element, as against
      ;; those carried by extrapolation from the nearest.
      :bound (count (filter #(zero? (aget outside %)) (range n-vis)))
      :unbound (count (filter #(neg? (aget tet-of %)) (range n-vis)))})))

(defn skin!
  "Writes the visual vertex positions implied by the current simulation
  positions into `out`.

  The whole per-frame cost of the technique: one weighted sum of four
  points per visual vertex, with no search and no simulation."
  [{:keys [n ^ints tet ^doubles bary]} tet-ids ^doubles pos ^doubles out]
  (let [ids (vec tet-ids)]
    (dotimes [i n]
      (let [t (aget tet i)]
        (when (>= t 0)
          (let [base (* 4 t)
                b0 (aget bary (* 3 i))
                b1 (aget bary (inc (* 3 i)))
                b2 (aget bary (+ 2 (* 3 i)))
                b3 (- 1.0 b0 b1 b2)
                p0 (* 3 (ids base))
                p1 (* 3 (ids (+ base 1)))
                p2 (* 3 (ids (+ base 2)))
                p3 (* 3 (ids (+ base 3)))
                o  (* 3 i)]
            (aset out o (+ (* b0 (aget pos p0)) (* b1 (aget pos p1))
                           (* b2 (aget pos p2)) (* b3 (aget pos p3))))
            (aset out (+ o 1) (+ (* b0 (aget pos (+ p0 1))) (* b1 (aget pos (+ p1 1)))
                                 (* b2 (aget pos (+ p2 1))) (* b3 (aget pos (+ p3 1)))))
            (aset out (+ o 2) (+ (* b0 (aget pos (+ p0 2))) (* b1 (aget pos (+ p1 2)))
                                 (* b2 (aget pos (+ p2 2))) (* b3 (aget pos (+ p3 2)))))))))
    out))

(defn skinned-positions
  "`skin!` into a fresh array, for callers that do not keep one."
  [skinning tet-ids pos]
  (let [^doubles out (f64 (* 3 (:n skinning)))]
    (skin! skinning tet-ids pos out)))
