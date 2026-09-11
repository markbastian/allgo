(ns allgo.physics.height-field
  "Height-field water: a surface stored as one height per column, and the
  bodies that float in it. After Ten Minute Physics 20.

  `allgo.physics.fluid` and `allgo.physics.flip` simulate a volume. This
  does not. The water is a single-valued function of the plane -- one
  height per column -- which throws away everything that is not a surface:
  no overturning, no breaking waves, no splash that leaves the water. What
  is left is the wave equation on a grid, which is cheap enough to run a
  whole lake at interactive rates, and looks right for anything that does
  not break.

  Two pieces, usable apart:

    the surface   heights and vertical velocities, stepped by the wave
                  equation. `step-surface!` alone gives ripples.
    the coupling  a body displaces water, and the water pushes back with
                  the weight of what was displaced. That is all buoyancy
                  is here, and it is what `couple!` does in both
                  directions at once.

  A body joins in by implementing `Floating`, which asks only how far the
  body reaches above and below each column it covers. `Sphere` does that
  with a circle's chord; a box or a hull would do it just as well, and
  nothing in the coupling knows the difference.

  The wave speed is capped rather than trusted: information cannot cross
  more than a cell per step without the grid going unstable, so the
  requested speed is clamped to the CFL limit for the step actually being
  taken."
  (:require [clojure.math :as math]))

(defn- f32 [n] #?(:clj (float-array n) :cljs (js/Float32Array. n)))

(defn- fill! [^floats a v] (dotimes [i (alength a)] (aset a i (float v))) a)

(defn- copy! [^floats dst ^floats src]
  (dotimes [i (alength dst)] (aset dst i (aget src i)))
  dst)

;; ---------------------------------------------------------------------------
;; The surface

(defn surface
  "A grid of water columns spanning `size-x` by `size-z`, `depth` deep.

  The grid is centred on the origin, because the bodies that float in it
  are placed in world coordinates and it is one less thing to convert."
  [{:keys [size-x size-z spacing depth]
    :or   {spacing 0.02 depth 0.8}}]
  (let [nx (inc (long (math/floor (/ size-x spacing))))
        nz (inc (long (math/floor (/ size-z spacing))))
        n  (* nx nz)]
    {:nx nx :nz nz :n n
     :spacing (double spacing) :depth (double depth)
     :size-x (double size-x) :size-z (double size-z)
     :cx (long (math/floor (/ nx 2.0)))
     :cz (long (math/floor (/ nz 2.0)))
     :heights     (fill! (f32 n) depth)
     :velocities  (f32 n)
     ;; How much body is in each column now, and how much was last step.
     ;; Only the difference matters: water is pushed aside by a body
     ;; *arriving*, not by one that is already there.
     :body-heights      (f32 n)
     :prev-body-heights (f32 n)}))

(defn index [{:keys [nz]} i j] (+ (* (long i) (long nz)) (long j)))

(defn height
  "The water height of column `i j`."
  [{:keys [^floats heights] :as s} i j]
  (aget heights (index s i j)))

(defn set-height! [{:keys [^floats heights] :as s} i j v]
  (aset heights (index s i j) (float v))
  s)

(defn column-x [{:keys [cx spacing]} i] (* (- (long i) (long cx)) (double spacing)))
(defn column-z [{:keys [cz spacing]} j] (* (- (long j) (long cz)) (double spacing)))

(defn nearest-column
  "The column containing world point `x z`, clamped to the grid."
  [{:keys [nx nz cx cz spacing]} x z]
  [(min (max 0 (+ (long cx) (long (math/floor (/ (double x) (double spacing))))))
        (dec (long nx)))
   (min (max 0 (+ (long cz) (long (math/floor (/ (double z) (double spacing))))))
        (dec (long nz)))])

(defn sample
  "The water height at world point `x z`, interpolated between columns.

  Interpolating matters for anything reading the surface at a point that
  is not a column -- a boat's waterline, a camera, a reflection -- because
  the nearest column alone makes a smooth surface look like stairs."
  [{:keys [nx nz cx cz spacing] :as s} x z]
  (let [sp (double spacing)
        fx (+ (double cx) (/ (double x) sp))
        fz (+ (double cz) (/ (double z) sp))
        i0 (min (max 0 (long (math/floor fx))) (- (long nx) 2))
        j0 (min (max 0 (long (math/floor fz))) (- (long nz) 2))
        tx (min 1.0 (max 0.0 (- fx i0)))
        tz (min 1.0 (max 0.0 (- fz j0)))]
    (+ (* (- 1.0 tx) (- 1.0 tz) (height s i0 j0))
       (* tx (- 1.0 tz) (height s (inc i0) j0))
       (* (- 1.0 tx) tz (height s i0 (inc j0)))
       (* tx tz (height s (inc i0) (inc j0))))))

(defn total-volume
  "The water's volume, which nothing but the coupling should change."
  [{:keys [^floats heights spacing n]}]
  (let [a (* (double spacing) (double spacing))]
    (loop [i 0 sum 0.0]
      (if (= i (long n)) (* sum a) (recur (inc i) (+ sum (aget heights i)))))))

(defn cfl-wave-speed
  "The fastest a wave may travel on this grid at this step size.

  A wave that crosses more than half a cell per step outruns the
  stencil, which can only see its neighbours, and the grid answers by
  oscillating harder every step until it overflows."
  [{:keys [spacing]} dt]
  (* 0.5 (/ (double spacing) (double dt))))

(defn splash!
  "Drops a disc of water `amount` deep and `radius` wide at `x z`.

  Nothing in the simulation needs this; it is how you poke the surface to
  see what it does."
  [{:keys [nx nz ^floats heights] :as s} x z radius amount]
  (let [[ci cj] (nearest-column s x z)
        span (long (math/ceil (/ (double radius) (:spacing s))))
        r2   (* (double radius) (double radius))]
    (doseq [i (range (max 0 (- ci span)) (min (long nx) (+ ci span 1)))
            j (range (max 0 (- cj span)) (min (long nz) (+ cj span 1)))]
      (let [dx (- (column-x s i) (double x))
            dz (- (column-z s j) (double z))
            d2 (+ (* dx dx) (* dz dz))]
        (when (< d2 r2)
          ;; Tapered, so the splash does not start life as a cylinder with
          ;; a discontinuity around its rim.
          (let [k (- 1.0 (/ d2 r2))
                id (index s i j)]
            (aset heights id (float (+ (aget heights id) (* (double amount) k))))))))
    s))

(def default-surface-world
  {:wave-speed  2.0
   :pos-damping 1.0
   :vel-damping 0.3})

(defn step-surface!
  "One step of the wave equation over the height field.

  Each column is accelerated by how far it sits below the average of its
  four neighbours, which is the discrete Laplacian and so the wave
  equation. The two dampings do different jobs: `:vel-damping` bleeds
  energy out, so the water eventually goes flat, while `:pos-damping`
  pulls each column toward its neighbours' average, which smooths the
  grid-scale jitter the stencil cannot resolve."
  ([s dt] (step-surface! s dt {}))
  ([{:keys [nx nz spacing ^floats heights ^floats velocities] :as s} dt world]
   (let [{:keys [wave-speed pos-damping vel-damping]}
         (merge default-surface-world world)
         nx (long nx) nz (long nz)
         dt (double dt)
         sp (double spacing)
         ;; Clamped, not trusted -- and clamped locally, so that asking for
         ;; a fast wave at a large step does not quietly lower the wave
         ;; speed stored on the surface for every step after it.
         c-speed (min (double wave-speed) (cfl-wave-speed s dt))
         c  (/ (* c-speed c-speed) (* sp sp))
         pd (min (* (double pos-damping) dt) 1.0)
         vd (max 0.0 (- 1.0 (* (double vel-damping) dt)))]
     (dotimes [i nx]
       (dotimes [j nz]
         (let [id (+ (* i nz) j)
               h  (aget heights id)
               ;; A missing neighbour reads as this column's own height,
               ;; which is a wall that reflects rather than absorbs.
               sum (+ (if (> i 0) (aget heights (- id nz)) h)
                      (if (< i (dec nx)) (aget heights (+ id nz)) h)
                      (if (> j 0) (aget heights (dec id)) h)
                      (if (< j (dec nz)) (aget heights (inc id)) h))]
           (aset velocities id
                 (float (+ (aget velocities id) (* dt c (- sum (* 4.0 h))))))
           (aset heights id
                 (float (+ h (* (- (* 0.25 sum) h) pd)))))))
     (dotimes [id (* nx nz)]
       (aset velocities id (float (* (aget velocities id) vd)))
       (aset heights id (float (+ (aget heights id) (* (aget velocities id) dt)))))
     s)))

;; ---------------------------------------------------------------------------
;; Bodies

(defprotocol Floating
  "Something that displaces water in a height field.

  The coupling needs to know two things and nothing else: where in the
  plane the body might be, and how thick it is over a given column. Any
  shape that can answer those floats."
  (footprint [body]
    "`[x0 z0 x1 z1]`, the body's extent in the plane. Only a bound; the
     coupling still asks each column inside it.")
  (half-height-at [body x z]
    "Half the body's vertical extent over the column at `x z`, or 0 where
     the column misses the body entirely.")
  (centre-y [body] "The body's centre height.")
  (body-volume [body] "The volume displaced when the body is fully under.")
  (add-vertical-velocity [body dv] "The body, moving `dv` faster upward."))

(defrecord Sphere [pos vel radius mass restitution]
  Floating
  (footprint [_]
    (let [[x _ z] pos r (double radius)]
      [(- x r) (- z r) (+ x r) (+ z r)]))
  (half-height-at [_ x z]
    (let [[px _ pz] pos
          dx (- (double px) (double x))
          dz (- (double pz) (double z))
          d2 (+ (* dx dx) (* dz dz))
          r2 (* (double radius) (double radius))]
      ;; Half the chord of the circle cut at this distance from centre.
      (if (< d2 r2) (math/sqrt (- r2 d2)) 0.0)))
  (centre-y [_] (nth pos 1))
  (body-volume [_] (* 4.0 (/ math/PI 3.0) radius radius radius))
  (add-vertical-velocity [b dv] (update-in b [:vel 1] + dv)))

(defn sphere
  "A floating sphere. `density` is relative to the water, so 1.0 is
  neutrally buoyant, below floats and above sinks."
  [{:keys [pos radius density restitution]
    :or   {density 1.0 restitution 0.1}}]
  (map->Sphere
   {:pos (vec pos)
    :vel [0.0 0.0 0.0]
    :radius (double radius)
    :mass (* 4.0 (/ math/PI 3.0) radius radius radius (double density))
    :restitution (double restitution)}))

(defn volume
  "What the body displaces when fully submerged."
  [body]
  (body-volume body))

;; ---------------------------------------------------------------------------
;; Coupling

(def default-coupling
  {:alpha     0.5
   :smoothing 2
   ;; Per second, and scaled by how much of the body is actually under --
   ;; see `couple!`. The reference's per-column damping works out near 8
   ;; for a body of this size at a 30Hz step; a little more than that is
   ;; wanted here because the drag is no longer applied to a body that has
   ;; left the water. It does not affect where a body finally floats, only
   ;; how long it takes to get there -- so raise it for a body much
   ;; lighter than water, which otherwise bobs for a very long time.
   :drag      12.0})

(defn- smooth-body-heights!
  "Spreads the displaced volume over neighbouring columns.

  A sphere's footprint has a hard rim -- the chord goes to zero at the
  edge but the column either contains it or does not -- and without this
  the water inherits that edge as a ring of grid-aligned spikes."
  [{:keys [nx nz ^floats body-heights]} iterations]
  (let [nx (long nx) nz (long nz)]
    (dotimes [_ (long iterations)]
      (dotimes [i nx]
        (dotimes [j nz]
          (let [id (+ (* i nz) j)
                n  (+ (if (and (> i 0) (< i (dec nx))) 2 1)
                      (if (and (> j 0) (< j (dec nz))) 2 1))
                sum (+ (if (> i 0) (aget body-heights (- id nz)) 0.0)
                       (if (< i (dec nx)) (aget body-heights (+ id nz)) 0.0)
                       (if (> j 0) (aget body-heights (dec id)) 0.0)
                       (if (< j (dec nz)) (aget body-heights (inc id)) 0.0))]
            (aset body-heights id (float (/ sum n)))))))))

(defn couple!
  "Displacement, both ways.

  Each column under a body gets the depth of body submerged in it. That
  depth times the column's area is displaced water, and its weight is the
  buoyant force on the body -- Archimedes, arrived at by measuring rather
  than by any special case. The same depths are written into the height
  field, so the water rises where the body went in.

  Only the *change* in displacement is added to the water. A body sitting
  still has already pushed its water aside and must not keep pushing it.

  Returns the bodies with buoyancy applied; the surface is mutated."
  ([s bodies dt gravity] (couple! s bodies dt gravity {}))
  ([{:keys [nx nz spacing ^floats heights ^floats body-heights
            ^floats prev-body-heights] :as s}
    bodies dt gravity world]
   (let [{:keys [alpha smoothing drag]} (merge default-coupling world)
         nx (long nx) nz (long nz)
         sp (double spacing)
         area (* sp sp)
         dt (double dt)
         g  (double gravity)
         ;; Last step's displacement, kept so only the difference moves water.
         _ (copy! prev-body-heights body-heights)
         _ (fill! body-heights 0.0)
         bodies
         (mapv
          (fn [body]
            (let [[bx0 bz0 bx1 bz1] (footprint body)
                  [i0 j0] (nearest-column s bx0 bz0)
                  [i1 j1] (nearest-column s bx1 bz1)
                  cy (double (centre-y body))
                  force
                  (loop [i (long i0) force 0.0]
                    (if (> i (long i1))
                      force
                      (recur
                       (inc i)
                       (loop [j (long j0) force force]
                         (if (> j (long j1))
                           force
                           (let [hh (double (half-height-at body (column-x s i) (column-z s j)))]
                             (if (pos? hh)
                               (let [id (+ (* i nz) j)
                                     water (aget heights id)
                                     ;; The body's extent, clipped to the
                                     ;; water it is actually inside: below
                                     ;; the surface and above the floor.
                                     lo (max (- cy hh) 0.0)
                                     hi (min (+ cy hh) water)
                                     depth (max (- hi lo) 0.0)]
                                 (if (pos? depth)
                                   (do (aset body-heights id
                                             (float (+ (aget body-heights id) depth)))
                                       (recur (inc j) (+ force (* depth area (- g)))))
                                   (recur (inc j) force)))
                               (recur (inc j) force))))))))
                  ;; The reference applies the buoyant force column by
                  ;; column and its velocity drag along with it, so a body
                  ;; covering four hundred columns is damped four hundred
                  ;; times a step and how much drag it feels depends on
                  ;; the grid resolution. Summing the force first and
                  ;; damping once makes drag a property of the water --
                  ;; but it has to be scaled by how much of the body is
                  ;; under, or a body in mid-air is dragged by water it is
                  ;; nowhere near, and a barely-dipped one is dragged as
                  ;; hard as a submerged one.
                  displaced (/ force (- g))
                  submerged (min 1.0 (/ displaced (max 1e-12 (double (body-volume body)))))
                  damped (max 0.0 (- 1.0 (* (double drag) dt submerged)))]
              ;; Drag first, then the impulse. The other order damps the
              ;; buoyant impulse but not the gravity added afterwards, and
              ;; that asymmetry moves the resting depth: a body has to sit
              ;; deeper, by a factor of 1/damped, before the two balance.
              ;; Archimedes then depends on the drag coefficient, which is
              ;; not a thing Archimedes depends on.
              (-> body
                  (update :vel (fn [[vx vy vz]]
                                 [(* vx damped) (* vy damped) (* vz damped)]))
                  (add-vertical-velocity (/ (* dt force) (double (:mass body)))))))
          bodies)]
     (smooth-body-heights! s smoothing)
     (dotimes [id (* nx nz)]
       (let [change (- (aget body-heights id) (aget prev-body-heights id))]
         (aset heights id (float (+ (aget heights id) (* (double alpha) change))))))
     bodies)))

;; ---------------------------------------------------------------------------
;; Bodies in a tank

(def default-world
  {:dt        (/ 1.0 30.0)
   :gravity   -10.0
   :tank      {:size [2.5 1.0 3.0] :border 0.03}
   :alpha     0.5
   :smoothing 2
   :drag      12.0
   :wave-speed 2.0
   :pos-damping 1.0
   :vel-damping 0.3})

(defn- clamp-to-tank
  "Keeps a body inside the tank, bouncing off the walls and the floor."
  [{:keys [radius restitution] :as body} {:keys [size border]}]
  (let [[sx _ sz] size
        r (double radius) e (double restitution) b (double border)
        wx (- (* 0.5 sx) r (* 0.5 b))
        wz (- (* 0.5 sz) r (* 0.5 b))
        [px py pz] (:pos body)
        [vx vy vz] (:vel body)
        [px vx] (cond (< px (- wx)) [(- wx) (* (- e) vx)]
                      (> px wx)     [wx (* (- e) vx)]
                      :else         [px vx])
        [pz vz] (cond (< pz (- wz)) [(- wz) (* (- e) vz)]
                      (> pz wz)     [wz (* (- e) vz)]
                      :else         [pz vz])
        [py vy] (if (< py r) [r (* (- e) vy)] [py vy])]
    (assoc body :pos [px py pz] :vel [vx vy vz])))

(defn step-body
  "Gravity, then motion, then the tank."
  [body dt gravity tank]
  (let [dt (double dt)
        [vx vy vz] (:vel body)
        vy (+ vy (* dt (double gravity)))
        [px py pz] (:pos body)]
    (-> body
        (assoc :vel [vx vy vz]
               :pos [(+ px (* dt vx)) (+ py (* dt vy)) (+ pz (* dt vz))])
        (clamp-to-tank tank))))

(defn- collide-pair
  "Separates two overlapping bodies and exchanges momentum along the line
  between them."
  [a b]
  (let [[ax ay az] (:pos a) [bx by bz] (:pos b)
        dx (- bx ax) dy (- by ay) dz (- bz az)
        d  (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz)))
        min-d (+ (double (:radius a)) (double (:radius b)))]
    (if (or (zero? d) (>= d min-d))
      [a b]
      (let [nx (/ dx d) ny (/ dy d) nz (/ dz d)
            corr (/ (- min-d d) 2.0)
            a (update a :pos (fn [[x y z]]
                               [(- x (* nx corr)) (- y (* ny corr)) (- z (* nz corr))]))
            b (update b :pos (fn [[x y z]]
                               [(+ x (* nx corr)) (+ y (* ny corr)) (+ z (* nz corr))]))
            [avx avy avz] (:vel a) [bvx bvy bvz] (:vel b)
            v1 (+ (* avx nx) (* avy ny) (* avz nz))
            v2 (+ (* bvx nx) (* bvy ny) (* bvz nz))
            m1 (double (:mass a)) m2 (double (:mass b))
            e  (double (:restitution a))
            nv1 (/ (- (+ (* m1 v1) (* m2 v2)) (* m2 (- v1 v2) e)) (+ m1 m2))
            nv2 (/ (- (+ (* m1 v1) (* m2 v2)) (* m1 (- v2 v1) e)) (+ m1 m2))
            da (- nv1 v1) db (- nv2 v2)]
        [(assoc a :vel [(+ avx (* nx da)) (+ avy (* ny da)) (+ avz (* nz da))])
         (assoc b :vel [(+ bvx (* nx db)) (+ bvy (* ny db)) (+ bvz (* nz db))])]))))

(defn resolve-collisions
  "Every pair, once."
  [bodies]
  (let [n (count bodies)]
    (loop [i 0 bs (vec bodies)]
      (if (>= i n)
        bs
        (recur (inc i)
               (loop [j (inc i) bs bs]
                 (if (>= j n)
                   bs
                   (let [[a b] (collide-pair (bs i) (bs j))]
                     (recur (inc j) (-> bs (assoc i a) (assoc j b)))))))))))

(defn step!
  "One frame of water and everything floating in it.

  Coupling first, so the buoyant force is worked out against the surface
  the bodies are actually sitting in, then the surface, then the bodies.

  Returns `{:surface :bodies}`; the surface is mutated in place and the
  bodies are replaced."
  ([world] (step! world {}))
  ([{:keys [surface bodies]} world]
   (let [{:keys [dt gravity tank] :as w} (merge default-world world)
         bodies (couple! surface bodies dt gravity w)]
     (step-surface! surface dt w)
     {:surface surface
      :bodies (-> (mapv #(step-body % dt gravity tank) bodies)
                  resolve-collisions)})))
