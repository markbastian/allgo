(ns allgo.physics.xpbd
  "Extended Position Based Dynamics, following Matthias Müller's Ten Minute
  Physics:

    https://matthias-research.github.io/pages/tenMinutePhysics/

  Ordinary physics integrates forces into velocities into positions.
  Position based dynamics works the other way: guess where the particles
  go, then repeatedly *move them* until the constraints hold, and read the
  velocity back off how far each one actually traveled. Nothing can
  explode, because nothing is ever integrated -- the worst a bad step can
  do is move a particle somewhere and have the next projection move it
  back.

  What XPBD adds to that is a physical meaning for stiffness. Plain PBD's
  stiffness depends on how many iterations you run and what the timestep
  is, so a body stiffened by running more iterations is not stiffer in any
  units you could measure. XPBD gives each constraint a *compliance* --
  the inverse of stiffness, in meters per newton -- and carries a Lagrange
  multiplier that makes the result converge to the same material whatever
  the iteration count. Compliance 0 is perfectly rigid.

  The engine is a body of particles plus constraint blocks:

    (-> (body-from-mesh mesh)
        (add-constraint (distance-constraint mesh 0.0))
        (add-constraint (volume-constraint mesh 0.0))
        (step! world))

  A constraint block holds many constraints of one kind, not one each: the
  solver walks a flat array rather than a list of objects, which is what
  keeps a few thousand tetrahedra inside a frame.

  Positions live in flat `[x y z x y z ...]` arrays, mutated in place. It
  is not the usual Clojure bargain, but a soft body is a few thousand
  constraints projected ten times a frame, and persistent vectors lose
  that by more than an order of magnitude. The mutation is confined to
  this namespace: a body is opaque, and `positions` copies out."
  (:require [allgo.array :as a]
            [allgo.spatial.hash :as spatial]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Flat arrays, portably

(defn- f64
  ([n] #?(:clj (double-array n) :cljs (js/Float64Array. n)))
  ([_n coll] #?(:clj (double-array (map double coll))
                :cljs (js/Float64Array. (into-array (map double coll))))))

;; ---------------------------------------------------------------------------
;; Bodies

(defn body
  "A body of particles from flat `verts`, all initially at rest.

  `:inv-mass` is inverse mass, not mass, because that is what the solver
  divides by and because zero is then the natural way to say *immovable*
  -- which is how pinning, grabbing and static geometry are all expressed."
  [verts]
  (let [n (quot (count verts) 3)]
    {:n           n
     :pos         (f64 (* 3 n) verts)
     :prev        (f64 (* 3 n) verts)
     :vel         (f64 (* 3 n))
     :inv-mass    (f64 n)
     :constraints []}))

(defn add-constraint [body constraint]
  (update body :constraints conj constraint))

(defn positions
  "The particle positions, copied out as a flat vector.

  Built element by element rather than with `vec`: ClojureScript's
  sequence functions do not reach into a typed array, so `vec` and `count`
  both fail on one. Clojure is happy with either, which is exactly how
  that asymmetry gets missed."
  [{:keys [pos]}]
  (let [^doubles pos pos
        n (alength pos)]
    (loop [i 0 out (transient [])]
      (if (= i n)
        (persistent! out)
        (recur (inc i) (conj! out (aget pos i)))))))

(defn particle
  "`[x y z]` of particle `i`."
  [body i]
  (let [^doubles pos (:pos body)
        b (* 3 i)]
    [(aget pos b) (aget pos (+ b 1)) (aget pos (+ b 2))]))

(defn set-particle!
  "Moves particle `i`, leaving its velocity alone."
  [body i [x y z]]
  (let [^doubles pos (:pos body)
        b (* 3 i)]
    (aset pos b (double x))
    (aset pos (+ b 1) (double y))
    (aset pos (+ b 2) (double z)))
  nil)

(defn inv-mass [body i] (aget ^doubles (:inv-mass body) i))

(defn set-inv-mass!
  "Zero makes a particle immovable; that is how pinning and grabbing work."
  [body i w]
  (aset ^doubles (:inv-mass body) i (double w))
  nil)

(defn set-velocity! [body i [x y z]]
  (let [^doubles vel (:vel body)
        b (* 3 i)]
    (aset vel b (double x))
    (aset vel (+ b 1) (double y))
    (aset vel (+ b 2) (double z)))
  nil)

;; ---------------------------------------------------------------------------
;; Mass from the mesh

(defn- tet-volume-at [^doubles pos a b c d]
  (let [a3 (* 3 a) b3 (* 3 b) c3 (* 3 c) d3 (* 3 d)
        x0 (aget pos a3) y0 (aget pos (+ a3 1)) z0 (aget pos (+ a3 2))
        e1x (- (aget pos b3) x0) e1y (- (aget pos (+ b3 1)) y0) e1z (- (aget pos (+ b3 2)) z0)
        e2x (- (aget pos c3) x0) e2y (- (aget pos (+ c3 1)) y0) e2z (- (aget pos (+ c3 2)) z0)
        e3x (- (aget pos d3) x0) e3y (- (aget pos (+ d3 1)) y0) e3z (- (aget pos (+ d3 2)) z0)
        cx (- (* e1y e2z) (* e1z e2y))
        cy (- (* e1z e2x) (* e1x e2z))
        cz (- (* e1x e2y) (* e1y e2x))]
    (/ (+ (* cx e3x) (* cy e3y) (* cz e3z)) 6.0)))

(defn distribute-mass!
  "Gives each particle the mass of the tetrahedra it belongs to.

  Mass follows the mesh rather than being spread evenly over the
  particles, so a body whose elements differ in size still behaves like
  uniform material.

  Masses are summed and inverted at the end, rather than inverse masses
  being summed as the tutorial does. Summing inverse masses combines them
  harmonically, which inverts the result: a particle sharing eight
  elements comes out four times *lighter* than one sharing two, when it
  owns four times the material and should be four times heavier. It goes
  unnoticed on a uniform mesh because the interior is then uniform either
  way and only the boundary differs."
  [body tet-ids density]
  (let [^doubles pos (:pos body)
        ^doubles inv-mass (:inv-mass body)
        tets (quot (count tet-ids) 4)
        ids  (vec tet-ids)]
    (dotimes [i (alength inv-mass)] (aset inv-mass i 0.0))
    (dotimes [t tets]
      (let [b   (* 4 t)
            vol (tet-volume-at pos (ids b) (ids (+ b 1)) (ids (+ b 2)) (ids (+ b 3)))
            ;; The quarter share of one tetrahedron's mass. Accumulated as
            ;; mass, not as its inverse -- see above.
            m   (if (pos? vol) (* density (/ vol 4.0)) 0.0)]
        (dotimes [j 4]
          (let [id (ids (+ b j))]
            (aset inv-mass id (+ (aget inv-mass id) m))))))
    (dotimes [i (alength inv-mass)]
      (let [m (aget inv-mass i)]
        (aset inv-mass i (if (pos? m) (/ 1.0 m) 0.0))))
    body))

;; ---------------------------------------------------------------------------
;; Constraints

(defprotocol Constraint
  (prepare! [constraint body]
    "Once per frame, before any substep.

  Where collision detection belongs. Finding which particles are near each
  other is far more expensive than resolving them, and over a frame they
  barely move -- so detect once against a radius widened by how far
  anything could travel, then resolve every substep against that list.")
  (reset-multipliers! [constraint]
    "Clears the accumulated Lagrange multipliers, once per substep.")
  (project! [constraint body dt]
    "One Gauss-Seidel pass: move the particles to reduce the violation."))

(defn- compliance-over [compliance dt]
  ;; alpha-tilde: compliance scaled into the timestep, the term that makes
  ;; stiffness independent of dt and of the iteration count.
  (/ compliance (* dt dt)))

(defrecord DistanceConstraints [^ints ids ^doubles rest-lengths ^doubles lambda compliance]
  Constraint
  (prepare! [_ _body] nil)
  (reset-multipliers! [_]
    (dotimes [i (alength lambda)] (aset lambda i 0.0)))
  (project! [_ body dt]
    (let [^doubles pos (:pos body)
          ^doubles inv-mass (:inv-mass body)
          alpha (compliance-over compliance dt)
          n     (alength rest-lengths)]
      (dotimes [i n]
        (let [i0 (aget ids (* 2 i))
              i1 (aget ids (inc (* 2 i)))
              w0 (aget inv-mass i0)
              w1 (aget inv-mass i1)
              w  (+ w0 w1)]
          (when (pos? w)
            (let [a  (* 3 i0) b (* 3 i1)
                  dx (- (aget pos a) (aget pos b))
                  dy (- (aget pos (+ a 1)) (aget pos (+ b 1)))
                  dz (- (aget pos (+ a 2)) (aget pos (+ b 2)))
                  len (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz)))]
              (when (pos? len)
                (let [nx (/ dx len) ny (/ dy len) nz (/ dz len)
                      c  (- len (aget rest-lengths i))
                      dl (/ (- (- c) (* alpha (aget lambda i))) (+ w alpha))]
                  (aset lambda i (+ (aget lambda i) dl))
                  (aset pos a (+ (aget pos a) (* nx dl w0)))
                  (aset pos (+ a 1) (+ (aget pos (+ a 1)) (* ny dl w0)))
                  (aset pos (+ a 2) (+ (aget pos (+ a 2)) (* nz dl w0)))
                  (aset pos b (- (aget pos b) (* nx dl w1)))
                  (aset pos (+ b 1) (- (aget pos (+ b 1)) (* ny dl w1)))
                  (aset pos (+ b 2) (- (aget pos (+ b 2)) (* nz dl w1))))))))))))

;; Which three vertices span the face opposite each vertex of a
;; tetrahedron, wound so the gradient points outward.
(def ^:private ^"[[I" volume-face-order
  #?(:clj (into-array (map int-array [[1 3 2] [0 2 3] [0 3 1] [0 1 2]]))
     :cljs (into-array (map into-array [[1 3 2] [0 2 3] [0 3 1] [0 1 2]]))))

(defrecord VolumeConstraints [^ints ids ^doubles rest-volumes ^doubles lambda compliance]
  Constraint
  (prepare! [_ _body] nil)
  (reset-multipliers! [_]
    (dotimes [i (alength lambda)] (aset lambda i 0.0)))
  (project! [_ body dt]
    (let [^doubles pos (:pos body)
          ^doubles inv-mass (:inv-mass body)
          alpha (compliance-over compliance dt)
          n     (alength rest-volumes)
          ^doubles grads (f64 12)]
      (dotimes [t n]
        ;; The gradient of a tetrahedron's volume with respect to one
        ;; vertex is the area vector of the opposite face -- moving a
        ;; vertex perpendicular to the face it faces is the only motion
        ;; that changes the volume.
        (let [base (* 4 t)
              w (loop [j 0 w 0.0]
                  (if (= j 4)
                    w
                    (let [^ints face (aget volume-face-order j)
                          p0 (* 3 (aget ids (+ base (aget face 0))))
                          p1 (* 3 (aget ids (+ base (aget face 1))))
                          p2 (* 3 (aget ids (+ base (aget face 2))))
                          e1x (- (aget pos p1) (aget pos p0))
                          e1y (- (aget pos (+ p1 1)) (aget pos (+ p0 1)))
                          e1z (- (aget pos (+ p1 2)) (aget pos (+ p0 2)))
                          e2x (- (aget pos p2) (aget pos p0))
                          e2y (- (aget pos (+ p2 1)) (aget pos (+ p0 1)))
                          e2z (- (aget pos (+ p2 2)) (aget pos (+ p0 2)))
                          gx (/ (- (* e1y e2z) (* e1z e2y)) 6.0)
                          gy (/ (- (* e1z e2x) (* e1x e2z)) 6.0)
                          gz (/ (- (* e1x e2y) (* e1y e2x)) 6.0)
                          g3 (* 3 j)]
                      (aset grads g3 gx)
                      (aset grads (+ g3 1) gy)
                      (aset grads (+ g3 2) gz)
                      (recur (inc j)
                             (+ w (* (aget inv-mass (aget ids (+ base j)))
                                     (+ (* gx gx) (* gy gy) (* gz gz))))))))]
          (when (pos? w)
            (let [vol (tet-volume-at pos (aget ids base) (aget ids (+ base 1))
                                     (aget ids (+ base 2)) (aget ids (+ base 3)))
                  c   (- vol (aget rest-volumes t))
                  dl  (/ (- (- c) (* alpha (aget lambda t))) (+ w alpha))]
              (aset lambda t (+ (aget lambda t) dl))
              (dotimes [j 4]
                (let [id (aget ids (+ base j))
                      wi (aget inv-mass id)
                      p  (* 3 id)
                      g3 (* 3 j)]
                  (aset pos p (+ (aget pos p) (* (aget grads g3) dl wi)))
                  (aset pos (+ p 1) (+ (aget pos (+ p 1)) (* (aget grads (+ g3 1)) dl wi)))
                  (aset pos (+ p 2) (+ (aget pos (+ p 2)) (* (aget grads (+ g3 2)) dl wi))))))))))))

(defrecord SphereCollision [cx cy cz radius friction]
  Constraint
  (prepare! [_ _body] nil)
  ;; Non-penetration is an inequality: it does nothing until it is
  ;; violated, and there is no compliance to accumulate against, so there
  ;; is no multiplier to keep.
  (reset-multipliers! [_] nil)
  (project! [_ body _dt]
    (let [^doubles pos (:pos body)
          ^doubles prev (:prev body)
          ^doubles inv-mass (:inv-mass body)
          n  (long (:n body))
          r  (double radius)
          mu (double friction)
          cx (double cx) cy (double cy) cz (double cz)]
      (dotimes [i n]
        (when (pos? (aget inv-mass i))
          (let [b  (* 3 i)
                dx (- (aget pos b) cx)
                dy (- (aget pos (+ b 1)) cy)
                dz (- (aget pos (+ b 2)) cz)
                d  (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz)))]
            (when (and (pos? d) (< d r))
              ;; Straight out to the surface along the radius, the shortest
              ;; move that resolves it.
              (let [s  (/ r d)
                    px (+ cx (* dx s))
                    py (+ cy (* dy s))
                    pz (+ cz (* dz s))
                    ;; Friction, the position-based way: of the distance
                    ;; traveled while in contact, undo some of the part
                    ;; running along the surface, leaving the normal part
                    ;; alone -- that part is the contact. Without it a sheet
                    ;; dropped on a ball slides off and pools on the floor,
                    ;; which is correct for a frictionless ball and useless
                    ;; as a drape.
                    nx (/ dx d) ny (/ dy d) nz (/ dz d)
                    mx (- px (aget prev b))
                    my (- py (aget prev (+ b 1)))
                    mz (- pz (aget prev (+ b 2)))
                    along (+ (* mx nx) (* my ny) (* mz nz))
                    tx (* mu (- mx (* along nx)))
                    ty (* mu (- my (* along ny)))
                    tz (* mu (- mz (* along nz)))]
                (aset pos b (- px tx))
                (aset pos (+ b 1) (- py ty))
                (aset pos (+ b 2) (- pz tz))))))))))

(defn sphere-constraint
  "Keeps every particle outside a sphere.

  A constraint rather than something applied once a frame, so it is
  enforced on every substep. Resolved per frame instead, a sheet moving at
  any speed passes straight through: the solver takes ten steps between
  checks, and the obstacle is only a few steps thick.

  `friction` from 0 to 1 is how much of the sliding along the surface is
  undone each step. At 0 a sheet dropped on a ball slides off and pools on
  the floor, which is right and dull; a little friction is what makes it
  drape."
  ([center radius] (sphere-constraint center radius 0.0))
  ([[cx cy cz] radius friction]
   (->SphereCollision (double cx) (double cy) (double cz) (double radius)
                      (double friction))))

(defrecord SelfCollision [^doubles rest-pos thickness friction hash state]
  Constraint
  (prepare! [_ body]
    ;; Detection once a frame, against a radius widened by how far a
    ;; particle can travel in one frame at the speed limit. Resolving every
    ;; substep against a list built once is the whole reason this is
    ;; affordable.
    (reset! state (spatial/adjacency hash (:pos body) (:n body)
                                     (* 2.0 (double thickness)))))
  (reset-multipliers! [_] nil)
  (project! [_ body _dt]
    (when-let [{:keys [^ints starts ^ints ids]} @state]
      (let [^doubles pos (:pos body)
            ^doubles prev (:prev body)
            ^doubles inv-mass (:inv-mass body)
            n  (long (:n body))
            th (double thickness)
            th2 (* th th)
            mu (double friction)]
        (dotimes [i n]
          (when (pos? (aget inv-mass i))
            (let [a (* 3 i)]
              (loop [k (aget starts i)]
                (when (< k (aget starts (inc i)))
                  (let [j (aget ids k)]
                    (when (pos? (aget inv-mass j))
                      (let [b  (* 3 j)
                            dx (- (aget pos b) (aget pos a))
                            dy (- (aget pos (+ b 1)) (aget pos (+ a 1)))
                            dz (- (aget pos (+ b 2)) (aget pos (+ a 2)))
                            d2 (+ (* dx dx) (* dy dy) (* dz dz))]
                        (when (and (pos? d2) (< d2 th2))
                          ;; Two particles that were always this close --
                          ;; neighbors in the sheet -- must not be pushed
                          ;; apart, or the cloth inflates. Their rest
                          ;; separation is the floor instead of the
                          ;; thickness.
                          (let [rx (- (aget rest-pos b) (aget rest-pos a))
                                ry (- (aget rest-pos (+ b 1)) (aget rest-pos (+ a 1)))
                                rz (- (aget rest-pos (+ b 2)) (aget rest-pos (+ a 2)))
                                r2 (+ (* rx rx) (* ry ry) (* rz rz))]
                            (when (<= d2 r2)
                              (let [min-d (if (< r2 th2) (math/sqrt r2) th)
                                    d     (math/sqrt d2)
                                    s     (/ (- min-d d) d 2.0)]
                                (aset pos a (- (aget pos a) (* dx s)))
                                (aset pos (+ a 1) (- (aget pos (+ a 1)) (* dy s)))
                                (aset pos (+ a 2) (- (aget pos (+ a 2)) (* dz s)))
                                (aset pos b (+ (aget pos b) (* dx s)))
                                (aset pos (+ b 1) (+ (aget pos (+ b 1)) (* dy s)))
                                (aset pos (+ b 2) (+ (aget pos (+ b 2)) (* dz s)))
                                (when (pos? mu)
                                  ;; Friction between the two layers: pull
                                  ;; both toward their average motion, so
                                  ;; folds hold instead of sliding apart.
                                  (let [vax (- (aget pos a) (aget prev a))
                                        vay (- (aget pos (+ a 1)) (aget prev (+ a 1)))
                                        vaz (- (aget pos (+ a 2)) (aget prev (+ a 2)))
                                        vbx (- (aget pos b) (aget prev b))
                                        vby (- (aget pos (+ b 1)) (aget prev (+ b 1)))
                                        vbz (- (aget pos (+ b 2)) (aget prev (+ b 2)))
                                        mx (* 0.5 (+ vax vbx))
                                        my (* 0.5 (+ vay vby))
                                        mz (* 0.5 (+ vaz vbz))]
                                    (aset pos a (+ (aget pos a) (* mu (- mx vax))))
                                    (aset pos (+ a 1) (+ (aget pos (+ a 1)) (* mu (- my vay))))
                                    (aset pos (+ a 2) (+ (aget pos (+ a 2)) (* mu (- mz vaz))))
                                    (aset pos b (+ (aget pos b) (* mu (- mx vbx))))
                                    (aset pos (+ b 1) (+ (aget pos (+ b 1)) (* mu (- my vby))))
                                    (aset pos (+ b 2) (+ (aget pos (+ b 2)) (* mu (- mz vbz))))))))))))
                    (recur (inc k))))))))))))

(defn self-collision-constraint
  "Keeps a body's own particles `thickness` apart, so cloth cannot pass
  through itself.

  `thickness` is the cloth's notional thickness, and should be somewhat
  less than the spacing between particles -- larger and the sheet holds
  itself rigid, since every neighbor is already inside it.

  Two particles that started out closer than the thickness are held at
  their *rest* separation rather than pushed to the thickness. Without
  that the constraint fights the sheet's own structure and inflates it,
  because a fine mesh's neighbors are closer together than the cloth is
  thick.

  Pair it with `:max-velocity` on the world -- `speed-limit` gives the
  right value -- or a fast-moving fold passes through itself between
  detections."
  ([body thickness] (self-collision-constraint body thickness 0.0))
  ([body thickness friction]
   (->SelfCollision (f64 (alength ^doubles (:pos body)) (seq (:pos body)))
                    (double thickness)
                    (double friction)
                    (spatial/spatial-hash (* 2.0 (double thickness)) (max 1 (:n body)))
                    (atom nil))))

(defn speed-limit
  "How fast a particle may travel and still not cross `thickness` in one
  frame. A fifth of the thickness per step, as the tutorial has it."
  [thickness dt]
  (/ (* 0.2 (double thickness)) (double dt)))

(defn distance-constraint
  "Holds pairs of particles at the distance they start out at.

  Given a mesh it uses `:edge-ids`, which resists stretching and shearing.
  Given an explicit id list it holds whatever pairs you name -- which is
  all a bending constraint is: the distance between the two vertices
  opposite a shared edge, at its own compliance. Cloth is this constraint
  twice over."
  ([mesh] (distance-constraint mesh 0.0))
  ([{:keys [verts edge-ids]} compliance]
   (distance-constraint verts edge-ids compliance))
  ([verts ids compliance]
   (let [n    (quot (count ids) 2)
         v    (vec verts)
         rest (f64 n (for [i (range n)]
                       (let [a (* 3 (nth ids (* 2 i)))
                             b (* 3 (nth ids (inc (* 2 i))))
                             dx (- (v a) (v b))
                             dy (- (v (+ a 1)) (v (+ b 1)))
                             dz (- (v (+ a 2)) (v (+ b 2)))]
                         (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz))))))]
     (->DistanceConstraints (a/i32 ids) rest (f64 n) compliance))))

(defn volume-constraint
  "Holds every tetrahedron at the volume it starts with. This is what
  makes the body resist being squashed rather than merely folded: edge
  constraints alone leave a lattice free to collapse like a hinge."
  ([mesh] (volume-constraint mesh 0.0))
  ([{:keys [verts tet-ids]} compliance]
   (let [n   (quot (count tet-ids) 4)
         pos (f64 (count verts) verts)
         rest (f64 n (for [t (range n)]
                       (let [b (* 4 t)]
                         (tet-volume-at pos (nth tet-ids b) (nth tet-ids (+ b 1))
                                        (nth tet-ids (+ b 2)) (nth tet-ids (+ b 3))))))]
     (->VolumeConstraints (a/i32 tet-ids) rest (f64 n) compliance))))

;; ---------------------------------------------------------------------------
;; The step

(def default-world
  {:gravity   [0.0 -10.0 0.0]
   :dt        (/ 1.0 60.0)
   :substeps  10
   :iterations 1
   :floor     0.0
   :damping   0.0
   ;; No limit by default; self-collision sets one from the thickness.
   :max-velocity nil})

(defn- pre-solve!
  "Integrate velocity, guess a new position, and stop anything below the
  floor.

  The floor is handled by putting the particle back where it came from and
  only then clamping its height, which cancels the tangential motion as
  well: that is where the friction of a body landing comes from."
  [body dt [gx gy gz] floor damping max-velocity]
  (let [n (:n body)
        ^doubles pos (:pos body)
        ^doubles prev (:prev body)
        ^doubles vel (:vel body)
        ^doubles inv-mass (:inv-mass body)
        decay (max 0.0 (- 1.0 (* damping dt)))
        vmax  (double (or max-velocity ##Inf))]
    (dotimes [i n]
      (when (pos? (aget inv-mass i))
        (let [b (* 3 i)
              vx (* decay (+ (aget vel b) (* gx dt)))
              vy (* decay (+ (aget vel (+ b 1)) (* gy dt)))
              vz (* decay (+ (aget vel (+ b 2)) (* gz dt)))
              ;; A speed limit, so nothing can cross a thin obstacle -- or
              ;; a sheet's own thickness -- inside one step. Scaled rather
              ;; than clipped per axis, which would change direction.
              speed (math/sqrt (+ (* vx vx) (* vy vy) (* vz vz)))
              k  (if (> speed vmax) (/ vmax speed) 1.0)
              vx (* vx k) vy (* vy k) vz (* vz k)]
          (aset vel b vx)
          (aset vel (+ b 1) vy)
          (aset vel (+ b 2) vz)
          (aset prev b (aget pos b))
          (aset prev (+ b 1) (aget pos (+ b 1)))
          (aset prev (+ b 2) (aget pos (+ b 2)))
          (aset pos b (+ (aget pos b) (* vx dt)))
          (aset pos (+ b 1) (+ (aget pos (+ b 1)) (* vy dt)))
          (aset pos (+ b 2) (+ (aget pos (+ b 2)) (* vz dt)))
          (when (and floor (< (aget pos (+ b 1)) floor))
            (aset pos b (aget prev b))
            (aset pos (+ b 2) (aget prev (+ b 2)))
            (aset pos (+ b 1) (double floor))))))))

(defn- post-solve!
  "Read velocity back off the distance actually traveled. This is the
  move that makes the scheme unconditionally stable: velocity is a
  *result* of where the solver put the particle, so it can never disagree
  with the positions or run away."
  [body dt]
  (let [n (:n body)
        ^doubles pos (:pos body)
        ^doubles prev (:prev body)
        ^doubles vel (:vel body)
        ^doubles inv-mass (:inv-mass body)
        inv-dt (/ 1.0 dt)]
    (dotimes [i n]
      (when (pos? (aget inv-mass i))
        (let [b (* 3 i)]
          (aset vel b (* (- (aget pos b) (aget prev b)) inv-dt))
          (aset vel (+ b 1) (* (- (aget pos (+ b 1)) (aget prev (+ b 1))) inv-dt))
          (aset vel (+ b 2) (* (- (aget pos (+ b 2)) (aget prev (+ b 2))) inv-dt)))))))

(defn substep!
  "One substep: predict, project, and read the velocity back.

  Merges `default-world`, so a partial world is enough here as it is for
  `step!`."
  [{:keys [constraints] :as body} world dt]
  (let [{:keys [gravity floor damping iterations max-velocity]} (merge default-world world)]
    (pre-solve! body dt gravity floor damping max-velocity)
    (doseq [c constraints] (reset-multipliers! c))
    (dotimes [_ (or iterations 1)]
      (doseq [c constraints] (project! c body dt)))
    (post-solve! body dt))
  body)

(defn step!
  "Advance the body by one frame.

  The frame is cut into `:substeps` and each is solved once, rather than
  the frame being solved `:substeps` times. That is the whole trick behind
  the claim that this cannot break: halving the timestep shrinks the
  distance a particle can travel into a bad configuration, while extra
  iterations only polish a bad guess. Ten substeps of one iteration beat
  one step of ten iterations on both stability and stiffness, for the same
  work."
  ([body] (step! body default-world))
  ([body world]
   (let [{:keys [dt substeps] :as world} (merge default-world world)
         sdt (/ dt substeps)]
     ;; Collision detection runs once here, not once per substep.
     (doseq [c (:constraints body)] (prepare! c body))
     (dotimes [_ substeps] (substep! body world sdt))
     body)))

;; ---------------------------------------------------------------------------
;; Interaction

(defn grab!
  "Pins the particle nearest `point`, returning `[body id]`.

  Its inverse mass is set to zero so the solver treats it as immovable and
  every constraint it touches resolves against it, which is what lets a
  body be dragged around by one vertex without tearing."
  [{:keys [n] :as body} [x y z]]
  (let [id (loop [i 0 best -1 best-d ##Inf]
             (if (= i n)
               best
               (let [[px py pz] (particle body i)
                     d (+ (* (- px x) (- px x)) (* (- py y) (- py y)) (* (- pz z) (- pz z)))]
                 (if (< d best-d) (recur (inc i) i d) (recur (inc i) best best-d)))))]
    (when (>= id 0)
      (let [w (inv-mass body id)]
        (set-inv-mass! body id 0.0)
        [(assoc body :grabbed {:id id :inv-mass w}) id]))))

(defn move-grabbed! [{:keys [grabbed] :as body} point]
  (when grabbed (set-particle! body (:id grabbed) point))
  body)

(defn release!
  "Gives the grabbed particle its mass back, and the velocity it was
  released with."
  ([body] (release! body [0.0 0.0 0.0]))
  ([{:keys [grabbed] :as body} velocity]
   (when grabbed
     (set-inv-mass! body (:id grabbed) (:inv-mass grabbed))
     (set-velocity! body (:id grabbed) velocity))
   (dissoc body :grabbed)))

;; ---------------------------------------------------------------------------

(defn distribute-area-mass!
  "Gives each particle a third of the area of every triangle it belongs
  to.

  The surface equivalent of `distribute-mass!`: a sheet has no volume, so
  its mass follows area instead. Spreading mass evenly over the particles
  would make a finely tessellated region heavier than a coarse one of the
  same size, and the sheet would sag wherever it happened to be detailed.

  Masses are summed and inverted at the end, for the reason given on
  `distribute-mass!`."
  [body tri-ids density]
  (let [^doubles pos (:pos body)
        ^doubles inv-mass (:inv-mass body)
        ids  (vec tri-ids)
        tris (quot (count ids) 3)]
    (dotimes [i (alength inv-mass)] (aset inv-mass i 0.0))
    (dotimes [t tris]
      (let [b  (* 3 t)
            i0 (* 3 (ids b)) i1 (* 3 (ids (+ b 1))) i2 (* 3 (ids (+ b 2)))
            e1x (- (aget pos i1) (aget pos i0))
            e1y (- (aget pos (+ i1 1)) (aget pos (+ i0 1)))
            e1z (- (aget pos (+ i1 2)) (aget pos (+ i0 2)))
            e2x (- (aget pos i2) (aget pos i0))
            e2y (- (aget pos (+ i2 1)) (aget pos (+ i0 1)))
            e2z (- (aget pos (+ i2 2)) (aget pos (+ i0 2)))
            cx (- (* e1y e2z) (* e1z e2y))
            cy (- (* e1z e2x) (* e1x e2z))
            cz (- (* e1x e2y) (* e1y e2x))
            a  (* 0.5 (math/sqrt (+ (* cx cx) (* cy cy) (* cz cz))))
            m  (if (pos? a) (* density (/ a 3.0)) 0.0)]
        (dotimes [j 3]
          (let [id (ids (+ b j))]
            (aset inv-mass id (+ (aget inv-mass id) m))))))
    (dotimes [i (alength inv-mass)]
      (let [m (aget inv-mass i)]
        (aset inv-mass i (if (pos? m) (/ 1.0 m) 0.0))))
    body))

(defn cloth
  "A body from a triangle mesh, with stretching and bending constraints.

  Both are distance constraints; what separates them is the pairs they
  hold and the compliance they hold them at. Stretching is usually rigid,
  since cloth barely stretches, while bending is where the character of a
  fabric lives: near zero is stiff like card, larger is limp like silk."
  ([mesh] (cloth mesh {}))
  ([{:keys [verts tri-ids edge-ids bend-ids]}
    {:keys [density stretch-compliance bend-compliance]
     :or   {density 1.0 stretch-compliance 0.0 bend-compliance 1.0}}]
   (cond-> (-> (body verts)
               (distribute-area-mass! tri-ids density)
               (add-constraint (distance-constraint verts edge-ids stretch-compliance)))
     (seq bend-ids)
     (add-constraint (distance-constraint verts bend-ids bend-compliance)))))

(defn pin!
  "Makes particles immovable, which is how a sheet is hung."
  [body ids]
  (doseq [i ids] (set-inv-mass! body i 0.0))
  body)

(defn soft-body
  "A body from a tetrahedral mesh with the two constraints a solid needs.

  `:edge-compliance` gives the material its give; `:volume-compliance`
  governs how much it may be squashed. Both are 0 for something rigid."
  ([mesh] (soft-body mesh {}))
  ([mesh {:keys [density edge-compliance volume-compliance]
          :or   {density 1000.0 edge-compliance 0.0 volume-compliance 0.0}}]
   (-> (body (:verts mesh))
       (distribute-mass! (:tet-ids mesh) density)
       (add-constraint (distance-constraint mesh edge-compliance))
       (add-constraint (volume-constraint mesh volume-compliance)))))

(defn mesh-volume
  "Current total volume of the tetrahedra, for checking how well the
  volume constraint is holding."
  [body tet-ids]
  (let [^doubles pos (:pos body)
        ids (vec tet-ids)]
    (reduce + (for [t (range (quot (count ids) 4))
                    :let [b (* 4 t)]]
                (tet-volume-at pos (ids b) (ids (+ b 1)) (ids (+ b 2)) (ids (+ b 3)))))))
