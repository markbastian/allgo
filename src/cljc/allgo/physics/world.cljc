(ns allgo.physics.world
  "A scene with both kinds of body in it.

  `allgo.physics.solver` simulates loose rigid bodies and is the right
  thing for a wall of bricks. `allgo.physics.articulated` simulates a
  jointed figure in reduced coordinates and is the right thing for a
  ragdoll. Until this namespace existed the two had never met: a ragdoll
  could fall down stairs, because stairs do not move, and could not
  knock a brick off a wall.

  ## Why they cannot simply be joined

  `allgo.physics.solver`'s sweep is written around flat per-body arrays
  of velocity, inverse mass and inverse inertia. That is most of why it
  is fast and it is not an accident. An articulated model has no `1/m`
  to put in one: its answer to *how much velocity does a push here buy*
  is a walk through a tree, and an impulse on it changes a generalised
  velocity rather than one body's own.

  So this is not a connector between two solvers. The split that works
  is per *contact*, not per solver -- a contact knows what is on each
  side of it, and only the sides that are articulated need the walk.

  ## What a participant has to answer

  Three questions, and every kind of body can answer them:

      how fast is the surface at `p` moving along `d`
      how much velocity does a unit impulse there buy
      here is an impulse, take it

  A static body answers nothing, nothing, and does nothing. A rigid body
  answers from its velocity and its inverse inertia. An articulated link
  answers from `allgo.physics.articulated/response-at`, which walks its
  tree. Once each side of a contact has answered, the sequential impulse
  that follows is the same arithmetic whatever asked it, and one
  Gauss-Seidel sweep covers all of them together -- which is what makes
  the answer consistent rather than two solvers arguing over the same
  brick.

  The expensive half of an articulated answer depends only on where its
  joints are, so it is built once a step by
  `allgo.physics.articulated/configuration` and every contact on that
  model shares it. Per iteration what remains is dot products.

  ## What this is not

  It is a sequential impulse solver and only that. There is no
  substepping, no XPBD, no sleeping, and warm starting lasts a step
  rather than crossing between them. `allgo.physics.solver` has all of
  those and is still the thing to reach for when the scene is only
  bricks; this one is for a scene that mixes kinds, and it says so
  rather than pretending to replace it.

  Two articulated models cannot yet touch each other. One can touch
  itself, and either can touch anything rigid."
  (:require [allgo.array :as a]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.contact :as contact]
            [allgo.physics.rigid :as rigid]))

(def defaults
  {:gravity [0.0 -9.81 0.0]
   :self-collide? true
   :iterations 8
   :friction 0.6
   :restitution 0.0
   :slop 0.005
   :bias-factor 0.2
   :max-push-speed 3.0})

(defn world
  "A scene of loose `bodies` and jointed `models`.

  Each model is `{:model m :pose p}` -- the description and the state
  `allgo.physics.articulated` keeps apart, kept apart here too."
  ([] (world nil nil nil))
  ([bodies models] (world bodies models nil))
  ([bodies models opts]
   (merge defaults opts {:bodies (vec bodies) :models (vec models)})))

;; ---------------------------------------------------------------------------
;; Who is on each side of a contact

;; A side is `[:rigid i]` or `[:link m i]`, where `m` indexes `:models`
;; and `i` is the link within it -- `-1` being that model's free root,
;; as everywhere else in `allgo.physics.articulated`.

(defn- rigid-side? [side] (= :rigid (first side)))

(defn- same-model? [a b]
  (and (not (rigid-side? a)) (not (rigid-side? b)) (= (second a) (second b))))

;; ---------------------------------------------------------------------------
;; Finding the contacts

(defn- model-bodies
  "Every shaped part of every model, as rigid bodies with the side that
  produced them."
  [models configs]
  (into []
        (for [[mi {:keys [model pose]}] (map-indexed vector models)
              {:keys [link body]} (ab/collision-bodies model (:q pose) (:base pose)
                                                       (:frames (nth configs mi)))]
          {:side [:link mi link] :body body})))

(defn- pair-contacts
  "Contacts between two bodies, with the normal pointing the way `a` has
  to be pushed."
  [side-a body-a side-b body-b]
  (mapv (fn [c] {:a side-a :b side-b
                 :point (:point c)
                 :normal (v/negate (:normal c))
                 :depth (:depth c)})
        (contact/between 0 1 body-a body-b)))

(defn contacts
  "Every contact in the scene: rigid against rigid, rigid against a
  model's parts, and a model against itself.

  The first of those goes through `allgo.physics.contact/all`, which has
  the broad phase; the rest are pair tests, because a model has few
  enough parts that sorting them would cost more than trying them."
  [{:keys [bodies models self-collide?] :or {self-collide? true}} configs]
  (let [parts (model-bodies models configs)]
    (-> []
        ;; Rigid against rigid.
        (into (map (fn [c] {:a [:rigid (:a c)] :b [:rigid (:b c)]
                            :point (:point c)
                            :normal (v/negate (:normal c))
                            :depth (:depth c)}))
              (contact/all bodies))
        ;; Every model's parts against every rigid body.
        (into (for [{:keys [side body]} parts
                    [bi rb] (map-indexed vector bodies)
                    c (pair-contacts side body [:rigid bi] rb)]
                c))
        ;; And each model against itself, unless told not to.
        (into (when self-collide?
                (for [[mi {:keys [model pose]}] (map-indexed vector models)
                      c (ab/self-contacts model (:q pose) (:base pose)
                                          (:frames (nth configs mi)))]
                  {:a [:link mi (:link c)] :b [:link mi (:other c)]
                   :point (:point c) :normal (:normal c) :depth (:depth c)}))))))

;; ---------------------------------------------------------------------------
;; What each side of a contact answers

;; A prepared side is a map the sweep reads every iteration, and nothing
;; in it is computed twice. For a rigid body that is the lever arm, the
;; inverse mass, and the angular response `I^-1 (r x d)` already turned
;; into the world -- `allgo.physics.solver` recomputes those cross
;; products every sweep, and there is no need to. For a link it is the
;; generalised force and the response to it, which came from one walk of
;; its tree.

(defn- prepare-rigid
  [bodies i point dir]
  (let [b (nth bodies i)]
    (if (rigid/static? b)
      {:kind :static}
      (let [r (v/sub point (:pos b))
            ;; I^-1 (r x d), in the world, which is both the angular
            ;; share of the effective mass and the spin an impulse adds.
            iw (q/rotate (:rot b)
                         (v/mul (:inv-inertia b)
                                (q/rotate (:inv-rot b) (v/cross r dir))))]
        {:kind :rigid :i i :r r :dir dir
         :im (double (:inv-mass b))
         :iw iw
         :w (+ (double (:inv-mass b)) (v/dot (v/cross r dir) iw))}))))

(defn- prepare-side
  [{:keys [bodies models]} configs side point dir]
  (if (rigid-side? side)
    (prepare-rigid bodies (second side) point dir)
    (let [[_ mi li] side
          {:keys [model]} (nth models mi)
          {:keys [force delta mass] :as _r} (let [r (ab/response-at model (nth configs mi)
                                                                    li point dir)]
                                              {:force (:g r) :delta (:delta r) :mass (:m r)})]
      {:kind :link :mi mi :force force :delta delta
       :w (if (pos? (double mass)) (/ 1.0 (double mass)) 0.0)})))

(defn- prepare-pair
  "Both sides at once, for a contact between two parts of one model.

  Asking each separately would be wrong, not merely slower: the two
  responses interact through everything the pair have in common, which
  is at least the root."
  [{:keys [models]} configs a b point dir]
  (let [[_ mi ai] a
        [_ _ bi] b
        {:keys [model]} (nth models mi)
        r (ab/pair-response-at model (nth configs mi) ai bi point dir)]
    {:kind :pair :mi mi :force (:g r) :delta (:delta r)
     :w (if (pos? (double (:m r))) (/ 1.0 (double (:m r))) 0.0)}))

;; ---------------------------------------------------------------------------
;; The sweep

(defn- body-arrays
  "Velocities out of the bodies and into two flat arrays."
  [bodies]
  (let [n (count bodies)
        vel (a/f64 (* 3 n))
        omega (a/f64 (* 3 n))]
    (dotimes [i n]
      (let [b (nth bodies i)
            [vx vy vz] (:vel b)
            [wx wy wz] (:omega b)]
        (aset vel (* 3 i) (double vx))
        (aset vel (+ (* 3 i) 1) (double vy))
        (aset vel (+ (* 3 i) 2) (double vz))
        (aset omega (* 3 i) (double wx))
        (aset omega (+ (* 3 i) 1) (double wy))
        (aset omega (+ (* 3 i) 2) (double wz))))
    [vel omega]))

(defn- write-back [bodies ^doubles vel ^doubles omega]
  (mapv (fn [i b]
          (if (rigid/static? b)
            b
            (assoc b
                   :vel [(aget vel (* 3 i)) (aget vel (+ (* 3 i) 1)) (aget vel (+ (* 3 i) 2))]
                   :omega [(aget omega (* 3 i)) (aget omega (+ (* 3 i) 1))
                           (aget omega (+ (* 3 i) 2))])))
        (range (count bodies))
        bodies))

(defn- side-speed
  "How fast this side's surface is moving along the contact direction."
  ^double [side ^doubles vel ^doubles omega us]
  (case (:kind side)
    :static 0.0
    :rigid (let [i (long (:i side))
                 [rx ry rz] (:r side)
                 [dx dy dz] (:dir side)
                 b (* 3 i)
                 wx (aget omega b) wy (aget omega (+ b 1)) wz (aget omega (+ b 2))]
             (+ (* dx (+ (aget vel b) (- (* wy (double rz)) (* wz (double ry)))))
                (* dy (+ (aget vel (+ b 1)) (- (* wz (double rx)) (* wx (double rz)))))
                (* dz (+ (aget vel (+ b 2)) (- (* wx (double ry)) (* wy (double rx)))))))
    (:link :pair) (let [^doubles f (:force side)
                        ^doubles u (nth us (:mi side))]
                    (loop [k 0 acc 0.0]
                      (if (= k (alength f))
                        acc
                        (recur (inc k) (+ acc (* (aget f k) (aget u k)))))))))

(defn- push!
  "Give this side an impulse of `lambda` along the contact direction.

  `sign` is +1 for the side the normal points towards and -1 for the
  other, which is the whole of Newton's third law as this solver needs
  to know it."
  [side ^doubles vel ^doubles omega us lambda sign]
  (let [s (* (double lambda) (double sign))]
    (case (:kind side)
      :static nil
      :rigid (let [i (long (:i side))
                   b (* 3 i)
                   [dx dy dz] (:dir side)
                   [ix iy iz] (:iw side)
                   im (double (:im side))]
               (aset vel b (+ (aget vel b) (* s im (double dx))))
               (aset vel (+ b 1) (+ (aget vel (+ b 1)) (* s im (double dy))))
               (aset vel (+ b 2) (+ (aget vel (+ b 2)) (* s im (double dz))))
               (aset omega b (+ (aget omega b) (* s (double ix))))
               (aset omega (+ b 1) (+ (aget omega (+ b 1)) (* s (double iy))))
               (aset omega (+ b 2) (+ (aget omega (+ b 2)) (* s (double iz)))))
      (:link :pair) (let [^doubles d (:delta side)
                          ^doubles u (nth us (:mi side))]
                      (dotimes [k (alength d)]
                        (aset u k (+ (aget u k) (* s (aget d k)))))))
    nil))

(defn- tangents [n]
  (let [a (if (< (abs (double (nth n 0))) 0.9) [1.0 0.0 0.0] [0.0 1.0 0.0])
        t1 (v/normalize (v/cross n a))]
    [t1 (v/cross n t1)]))

(defn- prepare
  "One contact, ready to be iterated: three directions, each with what
  both sides answer, and the bias that pushes an overlap apart."
  [w configs c dt opts]
  (let [dt (double dt)
        {:keys [slop bias-factor max-push-speed]} opts
        n (v/normalize (:normal c))
        [t1 t2] (tangents n)
        pair? (same-model? (:a c) (:b c))
        along (fn [d]
                (if pair?
                  (let [s (prepare-pair w configs (:a c) (:b c) (:point c) d)]
                    {:pair s :w (double (:w s))})
                  (let [sa (prepare-side w configs (:a c) (:point c) d)
                        sb (prepare-side w configs (:b c) (:point c) d)]
                    {:a sa :b sb :w (+ (double (or (:w sa) 0.0))
                                       (double (or (:w sb) 0.0)))})))]
    {:dirs [(along n) (along t1) (along t2)]
     :bias (min (double max-push-speed)
                (/ (* (double bias-factor)
                      (max 0.0 (- (double (:depth c)) (double slop))))
                   dt))}))

(defn- speed-of ^double [row ^doubles vel ^doubles omega us]
  (if-let [p (:pair row)]
    (side-speed p vel omega us)
    (- (side-speed (:a row) vel omega us)
       (side-speed (:b row) vel omega us))))

(defn- apply-to! [row ^doubles vel ^doubles omega us lambda]
  (if-let [p (:pair row)]
    (push! p vel omega us lambda 1.0)
    (do (push! (:a row) vel omega us lambda 1.0)
        (push! (:b row) vel omega us lambda -1.0))))

(defn- solve!
  "Sequential impulse over every contact in the scene, whatever is on
  each side of it.

  One sweep, not one per kind. That is the point of the whole namespace:
  a brick being pushed by a ragdoll's hand and by the brick beneath it
  has to see both pushes in the same iteration, or the two answers are
  computed against a velocity neither of them ends up with."
  [prepared ^doubles vel ^doubles omega us opts]
  (let [{:keys [iterations friction restitution]} opts
        k (count prepared)
        approach (a/f64 (mapv (fn [row]
                                (speed-of (nth (:dirs row) 0) vel omega us))
                              prepared))
        acc (a/f64 (* 3 k))]
    (dotimes [_ (long iterations)]
      (dotimes [i k]
        (let [row (nth prepared i)
              dirs (:dirs row)
              b (* 3 i)
              va (aget approach i)
              target (max (double (:bias row))
                          (if (< va -0.5) (* (- (double restitution)) va) 0.0))
              nrow (nth dirs 0)
              wn (double (:w nrow))
              vn (speed-of nrow vel omega us)
              old (aget acc b)
              nw (max 0.0 (+ old (if (> wn 1e-12) (/ (- target vn) wn) 0.0)))]
          (apply-to! nrow vel omega us (- nw old))
          (aset acc b nw)
          (let [limit (* (double friction) nw)]
            (dotimes [t 2]
              (let [row' (nth dirs (inc t))
                    wt (double (:w row'))
                    vt (speed-of row' vel omega us)
                    idx (+ b 1 t)
                    o (aget acc idx)
                    a' (min limit (max (- limit)
                                       (+ o (if (> wt 1e-12) (/ (- vt) wt) 0.0))))]
                (apply-to! row' vel omega us (- a' o))
                (aset acc idx a')))))))
    nil))

;; ---------------------------------------------------------------------------
;; A step

(defn step
  "One step of the whole scene.

  Velocities first for everybody, then one contact sweep over all of
  them together, then positions. The order is the same one the other
  two solvers use and for the same reason: a body is stopped before it
  moves into something, rather than pulled back out afterwards."
  [{:keys [bodies models gravity] :as w} dt]
  (let [dt (double dt)
        ;; Gravity on the loose bodies. `rigid/integrate` moves them as
        ;; well, so only the velocity half is wanted here.
        bodies (mapv (fn [b]
                       (if (rigid/static? b)
                         b
                         (update b :vel v/add-scaled gravity dt)))
                     bodies)
        ;; And the jointed ones, whose accelerations come from the whole
        ;; tree at once.
        models (mapv (fn [{:keys [model pose] :as m}]
                       (let [root (ab/base model)
                             {:keys [qdd base-acc]}
                             (ab/forward-dynamics model (:q pose) (:qd pose)
                                                  (vec (repeat (ab/dof model) 0.0))
                                                  (cond-> {:gravity gravity}
                                                    root (assoc :base (:base pose))))
                             accel (if root (vec (concat base-acc qdd)) (vec qdd))]
                         (assoc m :pose
                                (ab/with-velocity
                                  model pose
                                  (mapv (fn [a c] (+ (double a) (* dt (double c))))
                                        (ab/velocity model pose) accel)))))
                     models)
        w (assoc w :bodies bodies :models models)
        configs (mapv (fn [{:keys [model pose]}]
                        (ab/configuration model (:q pose) (:base pose)))
                      models)
        cs (contacts w configs)
        [vel omega] (body-arrays bodies)
        us (mapv (fn [{:keys [model pose]}] (a/f64 (ab/velocity model pose))) models)]
    (when (seq cs)
      (solve! (mapv #(prepare w configs % dt w) cs) vel omega us w))
    (let [bodies (write-back bodies vel omega)
          models (mapv (fn [m i]
                         (update m :pose
                                 #(ab/with-velocity (:model m) % (vec (seq ^doubles (nth us i))))))
                       models
                       (range (count models)))]
      (assoc w
             ;; Positions last, from whatever the sweep left behind.
             :bodies (mapv (fn [b]
                             (if (rigid/static? b)
                               b
                               (rigid/integrate b dt [0.0 0.0 0.0])))
                           bodies)
             ;; Only moved, not accelerated again: everything that was
             ;; going to change a velocity has already changed it.
             :models (mapv (fn [{:keys [model pose] :as m}]
                             (assoc m :pose (ab/advance pose model dt)))
                           models)))))
