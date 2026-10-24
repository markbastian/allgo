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
  is a walk through a tree, and an impulse on it changes a generalized
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

  There is no substepping and no XPBD. By default (`:solver :convex`)
  every contact, pin and joint limit is solved at once by Newton's
  method (`allgo.physics.convex`, `solve-convex!`): a column of thirty
  cubes stands where the sequential impulse holds six.
  `:solver :sequential-impulse` sweeps them instead, rigid with a
  Baumgarte push unless `:soft?` asks for springs; it is cheaper (the
  ragdoll scene 1.4 ms a step against 2.4), and it is what the
  motorcycle was tuned under and stays on (`motorcycle/world-defaults`).
  Bodies and models sleep, by islands,
  as `allgo.physics.solver`'s do (see Sleeping, below). `allgo.physics.solver`
  has all of those and is still the thing to reach for when the scene is
  only bricks; this one is for a scene that mixes kinds, and it says so
  rather than pretending to replace it.

  ## Warm starting

  Each contact remembers the impulse it ended a step with, keyed as
  `allgo.physics.solver` keys them, and starts the next step from it;
  so does each pin. A joint limit does not: its rows carry nothing that
  says which limit they are, and they are few.

  A contact the keys do not recognize takes the impulse of last step's
  nearest unclaimed point on the same pair, within `:recall-reach`
  (`recall`). The box test renames the points of a swaying brick as its
  corners cross the edges of the one below, and here, with its rigid
  contacts, a column sways enough for that to matter: matched by key
  alone it loses six bricks, by nearness too it stands them, and at 32
  iterations ten where it stood eight. `allgo.physics.solver`'s columns
  were no taller for it, and it matches by key alone.

  Two articulated models touch each other through the same sweep, and
  can be pinned together; one can touch itself, and either can touch
  anything rigid."
  (:require [allgo.array :as a]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.math :as am]
            [allgo.physics.articulated :as ab]
            [allgo.physics.contact :as contact]
            [allgo.physics.convex :as convex]
            [allgo.physics.rigid :as rigid]
            [clojure.math :as math]))

(def defaults
  {:gravity [0.0 -9.81 0.0]
   :self-collide? true
   ;; `:convex` solves every contact, pin and limit at once by Newton's
   ;; method; `:sequential-impulse` sweeps them `iterations` times.
   :solver :convex
   :iterations 8
   :warm-start? true
   ;; A contact whose key is new takes the impulse of last step's nearest
   ;; unclaimed point on the same pair within this many meters: see
   ;; `recall`. Nil matches by key alone.
   :recall-reach 0.01
   :friction 0.6
   :restitution 0.0
   :slop 0.005
   :bias-factor 0.2
   ;; Soft contacts, as `allgo.physics.solver` has them: an overlap is a
   ;; damped spring of `contact-hertz` and damping ratio `contact-damping`
   ;; instead of being pushed out at `bias-factor` of it a step. They push
   ;; a deep overlap out about a third as fast and rest a millimeter
   ;; deeper. They do not stand a taller column: columns here fall for
   ;; want of iterations, not of softness, and the solver's own
   ;; sequential impulse, soft all along, falls at the same height. So
   ;; they are there to be asked for, not the default.
   :soft? false
   :contact-hertz 30.0
   :contact-damping 10.0
   :max-push-speed 3.0
   ;; Sleeping, as `allgo.physics.solver` has it: an island still for
   ;; `time-to-sleep` stops being simulated, and wakes when something
   ;; awake touches it.
   :allow-sleep? true
   :time-to-sleep 0.5
   :linear-sleep-speed 0.01
   :angular-sleep-speed 0.035})

(defn world
  "A scene of loose `bodies` and jointed `models`.

  Each model is `{:model m :pose p}` -- the description and the state
  `allgo.physics.articulated` keeps apart, kept apart here too. A model
  may also carry `:tau`, the generalized forces on its joints for the
  next step, one per degree of freedom: a motor's torque, a spring's
  push. It stays until replaced, so a caller driving the model sets it
  before each step from wherever the model has got to.

  `:pins` in `opts` joins two bodies at a point: each is
  `{:a side :pa local :b side :pb local}`, where a side is `[:rigid i]`,
  `[:link model-index link-index]` or `[:static]` and the points are in
  that body's own frame. A pin given a `:break-force` lets go when the
  force through it exceeds that, and moves to `:broken-pins`. Given a
  `:break-time` too, it lets go only when the force averaged over that
  many seconds does: a rigid pin passes on any sudden change of velocity
  in a single step, as a force that is enormous and gone again, and a
  pin standing for something with give in it -- a hand on a grip --
  should not break on that. Each pin carries its `:load`, the force it
  is judged on.

  Joint limits on the models are solved in the same sweep as the
  contacts and pins."
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

(defn- asleep?
  "Whether the body on `side` is asleep: a loose body, or a whole model."
  [{:keys [bodies models]} side]
  (case (first side)
    :rigid (rigid/sleeping? (nth bodies (second side)))
    :link (true? (:sleeping? (nth models (second side))))
    false))

(defn- awake?
  "Whether the body on `side` moves this step: not static, not asleep."
  [{:keys [bodies models]} side]
  (case (first side)
    :rigid (not (rigid/inert? (nth bodies (second side))))
    :link (not (:sleeping? (nth models (second side))))
    false))

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
                 :depth (:depth c)
                 :id (:id c)})
        (contact/between 0 1 body-a body-b)))

(defn- reach
  "How far a body's shape extends from its center, at most: the radius of
  a ball that holds it, which is all the test below needs."
  ^double [{:keys [shape radius size major minor]}]
  (case shape
    :ball (double radius)
    :torus (+ (double major) (double minor))
    (* 0.5 (math/sqrt (reduce + (map #(* (double %) (double %)) size))))))

(defn- near?
  "Whether two bodies are close enough for a contact to be worth asking
  about -- their bounding balls, with the speculative margin to spare."
  [a b]
  (< (v/distance (:pos a) (:pos b)) (+ (reach a) (reach b) 0.1)))

(defn contacts
  "Every contact in the scene: rigid against rigid, rigid against a
  model's parts, one model's parts against another's, and a model
  against itself.

  The first of those goes through `allgo.physics.contact/all`, which has
  the broad phase; the rest are pair tests, because a model has few
  enough parts that sorting them would cost more than trying them.
  Between two models the pairs are thrown out by their bounding balls
  first, since most of two bodies' parts are nowhere near each other.

  `:collide-models?` false leaves different models passing through one
  another, which is what two bodies pinned together in an overlapping
  pose want -- a rider's knees tucked in against a tank."
  [{:keys [bodies models self-collide? collide-models?]
    :or {self-collide? true collide-models? true} :as w} configs]
  (let [parts (model-bodies models configs)]
    (-> []
        ;; Rigid against rigid.
        (into (map (fn [c] {:a [:rigid (:a c)] :b [:rigid (:b c)]
                            :point (:point c)
                            :normal (v/negate (:normal c))
                            :depth (:depth c)
                            :id (:id c)}))
              (contact/all bodies))
        ;; Every model's parts against every rigid body.
        (into (for [{:keys [side body]} parts
                    [bi rb] (map-indexed vector bodies)
                    ;; A course is dozens of boxes and a part is near
                    ;; two or three of them; the rest are ruled out by
                    ;; their bounding balls before any exact test.
                    :when (near? body rb)
                    ;; two bodies asleep, or asleep and static, have
                    ;; nothing to say to each other
                    :when (or (awake? w side) (awake? w [:rigid bi]))
                    c (pair-contacts side body [:rigid bi] rb)]
                c))
        ;; One model's parts against another's. Each side answers for its
        ;; own tree, exactly as it does against a rigid body: the two
        ;; models share nothing, so neither's response depends on the
        ;; other's.
        (into (when collide-models?
                (for [{side-a :side body-a :body} parts
                      {side-b :side body-b :body} parts
                      :when (< (long (second side-a)) (long (second side-b)))
                      :when (near? body-a body-b)
                      :when (or (awake? w side-a) (awake? w side-b))
                      c (pair-contacts side-a body-a side-b body-b)]
                  c)))
        ;; And each model against itself, unless told not to.
        (into (when self-collide?
                (for [[mi {:keys [model pose sleeping?]}] (map-indexed vector models)
                      :when (not sleeping?)
                      c (ab/self-contacts model (:q pose) (:base pose)
                                          (:frames (nth configs mi)))]
                  {:a [:link mi (:link c)] :b [:link mi (:other c)]
                   :point (:point c) :normal (:normal c) :depth (:depth c)
                   :id (:id c)}))))))

;; ---------------------------------------------------------------------------
;; What each side of a contact answers

;; A prepared side is a map the sweep reads every iteration, and nothing
;; in it is computed twice. For a rigid body that is the lever arm, the
;; inverse mass, and the angular response `I^-1 (r x d)` already turned
;; into the world -- `allgo.physics.solver` recomputes those cross
;; products every sweep, and there is no need to. For a link it is the
;; generalized force and the response to it, which came from one walk of
;; its tree.

(defn- turn!
  "`q/rotate` on arrays: `v` (at 3, in `p`) turned by the quaternion at
  `qo` in `q`, into `p` at 3 -- the same arithmetic in the same order."
  [^doubles q qo ^doubles p]
  (let [o (long qo)
        qx (aget q o) qy (aget q (+ o 1)) qz (aget q (+ o 2)) qw (aget q (+ o 3))
        vx (aget p 3) vy (aget p 4) vz (aget p 5)
        tx (* (- (* qy vz) (* qz vy)) 2.0)
        ty (* (- (* qz vx) (* qx vz)) 2.0)
        tz (* (- (* qx vy) (* qy vx)) 2.0)]
    (aset p 3 (+ (+ vx (* tx qw)) (- (* qy tz) (* qz ty))))
    (aset p 4 (+ (+ vy (* ty qw)) (- (* qz tx) (* qx tz))))
    (aset p 5 (+ (+ vz (* tz qw)) (- (* qx ty) (* qy tx))))))

(defn- prepare-rigid
  [bodies i point dir]
  (let [b (nth bodies i)]
    (if (rigid/inert? b)
      {:kind :static}
      ;; In primitives, with the arithmetic of `v/sub`, `v/cross`,
      ;; `q/rotate`, `v/mul` and `v/dot` it replaced, in their order: as
      ;; vectors this was a fifth of a sequential impulse step in a
      ;; browser, where reading a vector is a generic `nth`.
      (let [[px py pz] point [bx by bz] (:pos b) [dx dy dz] dir [ix iy iz] (:inv-inertia b)
            dx (double dx) dy (double dy) dz (double dz)
            rx (- (double px) (double bx)) ry (- (double py) (double by)) rz (- (double pz) (double bz))
            cx (- (* ry dz) (* rz dy)) cy (- (* rz dx) (* rx dz)) cz (- (* rx dy) (* ry dx))
            ;; q, the inverse rotation then the rotation; p, the working
            ;; vector at 3 and r and dir around it once it is done
            q (a/f64 (concat (:inv-rot b) (:rot b)))
            p (a/f64 9)
            _ (do (aset p 3 cx) (aset p 4 cy) (aset p 5 cz)
                  (turn! q 0 p)
                  (aset p 3 (* (double ix) (aget p 3)))
                  (aset p 4 (* (double iy) (aget p 4)))
                  (aset p 5 (* (double iz) (aget p 5)))
                  (turn! q 4 p))
            ;; I^-1 (r x d), in the world, which is both the angular
            ;; share of the effective mass and the spin an impulse adds
            wx (aget p 3) wy (aget p 4) wz (aget p 5)
            im (double (:inv-mass b))]
        ;; packed for the sweep to read every iteration: r, dir, I^-1 (r x d)
        (aset p 0 rx) (aset p 1 ry) (aset p 2 rz)
        (aset p 3 dx) (aset p 4 dy) (aset p 5 dz)
        (aset p 6 wx) (aset p 7 wy) (aset p 8 wz)
        {:kind :rigid :i i :r [rx ry rz] :dir dir :im im :p p
         :w (+ im (+ (* cx wx) (* cy wy) (* cz wz)))}))))

(defn- prepare-side
  [{:keys [bodies models] :as w} configs side point dir]
  (cond
    (= :static (first side)) {:kind :static}
    ;; asleep, a model is as fixed as the ground until something wakes it
    (asleep? w side) {:kind :static}
    (rigid-side? side)
    (prepare-rigid bodies (second side) point dir)
    ;; The convex solve holds each model's mass matrix itself and asks
    ;; only for the row; the response is most of the cost of a side.
    (= :convex (:solver w))
    (let [[_ mi li] side]
      {:kind :link :mi mi :force (ab/force-at (:model (nth models mi)) (nth configs mi) li point dir)})
    :else
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
  [{:keys [models] :as w} configs a b point dir]
  (let [[_ mi ai] a
        [_ _ bi] b
        {:keys [model]} (nth models mi)]
    (if (= :convex (:solver w))
      {:kind :pair :mi mi :force (ab/pair-force-at model (nth configs mi) ai bi point dir)}
      (let [r (ab/pair-response-at model (nth configs mi) ai bi point dir)]
        {:kind :pair :mi mi :force (:g r) :delta (:delta r)
         :w (if (pos? (double (:m r))) (/ 1.0 (double (:m r))) 0.0)}))))

;; ---------------------------------------------------------------------------
;; The sweep

;; ---------------------------------------------------------------------------
;; Pins and limits

(defn- side-frame
  "Where the body on this side of a pin is now."
  [{:keys [bodies models]} configs side]
  (case (first side)
    :static {:pos v/zero :rot q/identity-q}
    :rigid (let [b (nth bodies (second side))] {:pos (:pos b) :rot (:rot b)})
    :link (let [[_ mi li] side
                {:keys [model pose]} (nth models mi)]
            (ab/frame-of model (:q pose) (:base pose) li (:frames (nth configs mi))))))

(defn pin-points
  "Both ends of `pin`, in the world, as `[pa pb]`."
  [w configs {:keys [a b pa pb]}]
  (let [fa (side-frame w configs a)
        fb (side-frame w configs b)]
    [(v/add (:pos fa) (q/rotate (:rot fa) pa))
     (v/add (:pos fb) (q/rotate (:rot fb) pb))]))

(def ^:private axes [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]])

(defn- prepare-pin
  "A pin as three rows, one per world axis, each asking the two ends to
  close whatever gap has opened along it this step.

  Bilateral, where a contact is one-sided: a pin pulls as readily as it
  pushes, so its impulse is never clamped. That is what makes it a joint
  rather than a surface -- and what lets it close a loop that a tree of
  reduced coordinates cannot, a rider's two hands on one handlebar."
  [w configs pin dt opts]
  (let [{:keys [a b]} pin
        [pa pb] (pin-points w configs pin)
        gap (v/sub pb pa)
        dt (double dt)
        bias-factor (double (or (:pin-bias opts) 0.3))
        max-push (* 2.0 (double (:max-push-speed opts)))
        pair? (and (same-model? a b) (not (asleep? w a)))
        mid (v/scale (v/add pa pb) 0.5)]
    {:rows (mapv (fn [d]
                   (if pair?
                     (let [s (prepare-pair w configs a b mid d)]
                       {:pair s :w (double (or (:w s) 0.0))})
                     (let [sa (prepare-side w configs a pa d)
                           sb (prepare-side w configs b pb d)]
                       {:a sa :b sb :w (+ (double (or (:w sa) 0.0))
                                          (double (or (:w sb) 0.0)))})))
                 axes)
     :bias (mapv (fn [d]
                   (am/clamp (/ (* bias-factor (v/dot gap d)) dt) (- max-push) max-push))
                 axes)}))

(defn- prepare-limits
  "Every model's joints that are past their limits, as one-sided rows on
  that model's generalized velocity."
  [{:keys [models] :as opts} configs dt]
  (into []
        (for [[mi {:keys [model]}] (map-indexed vector models)
              :when (and (ab/limited? model) (not (:sleeping? (nth models mi))))
              row (ab/limit-constraints model (nth configs mi) dt opts)]
          (assoc row :mi mi))))

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
                 ^doubles p (:p side)
                 rx (aget p 0) ry (aget p 1) rz (aget p 2)
                 dx (aget p 3) dy (aget p 4) dz (aget p 5)
                 b (* 3 i)
                 wx (aget omega b) wy (aget omega (+ b 1)) wz (aget omega (+ b 2))]
             (+ (* dx (+ (aget vel b) (- (* wy rz) (* wz ry))))
                (* dy (+ (aget vel (+ b 1)) (- (* wz rx) (* wx rz))))
                (* dz (+ (aget vel (+ b 2)) (- (* wx ry) (* wy rx))))))
    (:link :pair) (let [^doubles f (:force side)
                        ^doubles u (nth us (:mi side))]
                    (loop [k 0 acc 0.0]
                      (if (== k (alength f))
                        acc
                        (recur (inc k) (+ acc (* (aget f k) (aget u k)))))))))

(defn- push!
  "Give this side an impulse of `lambda` along the contact direction.

  `sign` is +1 for the side the normal points toward and -1 for the
  other, which is the whole of Newton's third law as this solver needs
  to know it."
  [side ^doubles vel ^doubles omega us lambda sign]
  (let [s (* (double lambda) (double sign))]
    (case (:kind side)
      :static nil
      :rigid (let [i (long (:i side))
                   b (* 3 i)
                   ^doubles p (:p side)
                   im (double (:im side))]
               (aset vel b (+ (aget vel b) (* s im (aget p 3))))
               (aset vel (+ b 1) (+ (aget vel (+ b 1)) (* s im (aget p 4))))
               (aset vel (+ b 2) (+ (aget vel (+ b 2)) (* s im (aget p 5))))
               (aset omega b (+ (aget omega b) (* s (aget p 6))))
               (aset omega (+ b 1) (+ (aget omega (+ b 1)) (* s (aget p 7))))
               (aset omega (+ b 2) (+ (aget omega (+ b 2)) (* s (aget p 8)))))
      (:link :pair) (let [^doubles d (:delta side)
                          ^doubles u (nth us (:mi side))]
                      (dotimes [k (alength d)]
                        (aset u k (+ (aget u k) (* s (aget d k)))))))
    nil))

(defn- tangents [n]
  (let [a (if (< (abs (double (nth n 0))) 0.9) [1.0 0.0 0.0] [0.0 1.0 0.0])
        t1 (v/normalize (v/cross n a))]
    [t1 (v/cross n t1)]))

(defn- softness
  "The soft contact's three numbers for a step of `dt`, as
  `allgo.physics.solver` has them (Catto's soft constraints): the rate an
  overlap is pushed out at, how much of the ideal impulse to apply, and
  how much of what has accumulated to give back. The frequency is capped
  at a quarter of the step rate, above which the spring would ring
  between steps rather than damp."
  [{:keys [contact-hertz contact-damping]} dt]
  (let [dt (double dt)
        hertz (min (double contact-hertz) (* 0.25 (/ 1.0 dt)))
        zeta (double contact-damping)
        w (* 2.0 math/PI hertz)
        c (* dt w (+ (* 2.0 zeta) (* dt w)))]
    {:bias-rate (/ w (+ (* 2.0 zeta) (* dt w)))
     :mass-scale (/ c (+ 1.0 c))
     :impulse-scale (/ 1.0 (+ 1.0 c))}))

(defn- prepare
  "One contact, ready to be iterated: three directions, each with what
  both sides answer, and the bias that pushes an overlap apart -- with,
  for a soft contact, the scales the normal is solved with."
  [w configs c dt opts]
  (let [dt (double dt)
        {:keys [slop bias-factor max-push-speed soft?]} opts
        n (v/normalize (:normal c))
        [t1 t2] (tangents n)
        pair? (same-model? (:a c) (:b c))
        along (fn [d]
                (if pair?
                  (let [s (prepare-pair w configs (:a c) (:b c) (:point c) d)]
                    {:pair s :w (double (or (:w s) 0.0))})
                  (let [sa (prepare-side w configs (:a c) (:point c) d)
                        sb (prepare-side w configs (:b c) (:point c) d)]
                    {:a sa :b sb :w (+ (double (or (:w sa) 0.0))
                                       (double (or (:w sb) 0.0)))})))]
    (merge
     {:dirs [(along n) (along t1) (along t2)]
      :mass-scale 1.0 :impulse-scale 0.0}
     (if soft?
       ;; Shallower than the slop, the contact is rigid and lets the two
       ;; close only as far as the slop; deeper, it is the spring.
       (let [over (- (double (:depth c)) (double slop))]
         (if (pos? over)
           (let [{:keys [bias-rate mass-scale impulse-scale]} (or (:softness opts) (softness opts dt))]
             {:bias (min (double max-push-speed) (* (double bias-rate) over))
              :mass-scale mass-scale :impulse-scale impulse-scale})
           {:bias (/ over dt)}))
       {:bias (min (double max-push-speed)
                   (/ (* (double bias-factor)
                         (max 0.0 (- (double (:depth c)) (double slop))))
                      dt))}))))

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
  computed against a velocity neither of them ends up with.

  `warm` is the impulse each contact ended last step with, as
  `[normal t1 t2]` or nil, and `pin-warm` each pin's along the world
  axes. Both are applied before the first iteration and the sweep
  corrects from there. Returns what each contact and each pin carried
  over the step, as `{:contacts :pins}` in the same shapes."
  [prepared limits pins ^doubles vel ^doubles omega us opts warm pin-warm]
  (let [{:keys [iterations friction restitution]} opts
        k (count prepared)
        nl (count limits)
        np (count pins)
        lacc (a/f64 (max 1 nl))
        pacc (a/f64 (max 1 (* 3 np)))
        ;; Measured before the warm start, which would otherwise hide
        ;; the approach that restitution answers.
        approach (a/f64 (mapv (fn [row]
                                (speed-of (nth (:dirs row) 0) vel omega us))
                              prepared))
        acc (a/f64 (* 3 k))]
    (dotimes [i k]
      (when-let [was (nth warm i)]
        (let [dirs (:dirs (nth prepared i))]
          (dotimes [t 3]
            (let [p (double (nth was t))]
              (when-not (zero? p)
                (apply-to! (nth dirs t) vel omega us p)
                (aset acc (+ (* 3 i) t) p)))))))
    (dotimes [i np]
      (when-let [was (nth pin-warm i)]
        (let [rows (:rows (nth pins i))]
          (dotimes [r 3]
            (let [p (double (nth was r))]
              (when (and (not (zero? p)) (> (double (:w (nth rows r))) 1e-12))
                (apply-to! (nth rows r) vel omega us p)
                (aset pacc (+ (* 3 i) r) p)))))))
    (dotimes [_ (long iterations)]
      ;; Limits, then pins, then contacts: a knee folded backward is a
      ;; worse thing to see than a hand a millimeter off the bars, and
      ;; both are worse than a foot a millimeter into the floor.
      (dotimes [i nl]
        (let [{:keys [m bias mi]} (nth limits i)
              ^doubles g (:g (nth limits i))
              ^doubles d (:delta (nth limits i))
              ^doubles u (nth us mi)
              vn (loop [j 0 acc 0.0]
                   (if (== j (alength g)) acc (recur (inc j) (+ acc (* (aget g j) (aget u j))))))
              old (aget lacc i)
              nw (max 0.0 (+ old (* (double m) (- (double bias) vn))))
              dl (- nw old)]
          (dotimes [j (alength d)]
            (aset u j (+ (aget u j) (* dl (aget d j)))))
          (aset lacc i nw)))
      (dotimes [i np]
        (let [{:keys [rows bias]} (nth pins i)]
          (dotimes [r 3]
            (let [row (nth rows r)
                  wr (double (:w row))]
              (when (> wr 1e-12)
                (let [lambda (/ (- (double (nth bias r)) (speed-of row vel omega us)) wr)
                      idx (+ (* 3 i) r)]
                  (apply-to! row vel omega us lambda)
                  (aset pacc idx (+ (aget pacc idx) lambda))))))))
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
              ;; a soft normal gives back some of what it has
              ;; accumulated each time, which is what keeps it from
              ;; storing energy; a hard one has scales 1 and 0
              nw (max 0.0 (+ old (if (> wn 1e-12)
                                   (- (/ (* (double (:mass-scale row)) (- target vn)) wn)
                                      (* (double (:impulse-scale row)) old))
                                   0.0)))]
          (apply-to! nrow vel omega us (- nw old))
          (aset acc b nw)
          (let [limit (* (double friction) nw)]
            (dotimes [t 2]
              (let [row' (nth dirs (inc t))
                    wt (double (:w row'))
                    vt (speed-of row' vel omega us)
                    idx (+ b 1 t)
                    o (aget acc idx)
                    a' (am/clamp (+ o (if (> wt 1e-12) (/ (- vt) wt) 0.0)) (- limit) limit)]
                (apply-to! row' vel omega us (- a' o))
                (aset acc idx a')))))))
    {:contacts (mapv (fn [i] [(aget acc (* 3 i)) (aget acc (+ (* 3 i) 1))
                              (aget acc (+ (* 3 i) 2))])
                     (range k))
     :pins (mapv (fn [i] [(aget pacc (* 3 i)) (aget pacc (+ (* 3 i) 1))
                          (aget pacc (+ (* 3 i) 2))])
                 (range np))}))

;; ---------------------------------------------------------------------------
;; Every constraint at once

(defn- rigid-block
  "A loose body's six numbers: its mass and inverse as dense 6x6s, the
  world inertia beside `m I`."
  [{:keys [inv-mass inv-inertia rot]}]
  (let [[[m00 m01 m02] [m10 m11 m12] [m20 m21 m22]] (q/to-matrix rot)
        r [[m00 m01 m02] [m10 m11 m12] [m20 m21 m22]]
        world (fn [diag] (vec (for [i (range 3) j (range 3)]
                                (reduce + (map (fn [k] (* (double (get-in r [i k])) (double (nth diag k))
                                                          (double (get-in r [j k]))))
                                               (range 3))))))
        im (double inv-mass)
        inertia (world (mapv #(/ 1.0 (double %)) inv-inertia))
        inverse (world (vec inv-inertia))
        six (fn [lin ang] (vec (for [i (range 6) j (range 6)]
                                 (cond (and (< i 3) (< j 3)) (if (== i j) lin 0.0)
                                       (and (>= i 3) (>= j 3)) (nth ang (+ (* 3 (- i 3)) (- j 3)))
                                       :else 0.0))))]
    [(six (/ 1.0 im) inertia) (six im inverse)]))

(defn- side-row
  "What the body on a prepared side is, as a block, and its stretch of a
  row: a loose body's direction and lever arm crossed with it, a model's
  generalized force -- nil for anything that does not move."
  [side]
  (case (:kind side)
    :static nil
    :rigid (let [[rx ry rz] (:r side) [dx dy dz] (:dir side)]
             [[:rigid (:i side)]
              [dx dy dz (- (* ry dz) (* rz dy)) (- (* rz dx) (* rx dz)) (- (* rx dy) (* ry dx))]])
    (:link :pair) [[:model (:mi side)] (vec (:force side))]))

(defn- solve-convex!
  "`solve!`'s job, by `allgo.physics.convex`: every contact, pin and joint
  limit in the scene as one convex problem, solved by Newton's method --
  what MuJoCo does, and what lets a column stand as tall as it should.

  The blocks are whatever the rows touch that is awake: a loose body's
  six numbers, or a model's whole generalized velocity with its mass
  matrix, `articulated/generalized-mass-matrix`. Contacts are
  regularized and asked to push out an overlap as `allgo.physics.solver`'s
  convex solver asks; a joint limit is a one-sided row asking for its
  bias, a pin three two-sided ones. Writes the velocities back into
  `vel`, `omega` and `us`, and returns what each contact and pin carried,
  as `solve!` does, and the Newton factors as `:convex-factors`: given
  back next step while the blocks are the same bodies, they precondition
  its Newton steps."
  [w configs cs prepared limits pins ^doubles vel ^doubles omega us dt warm pin-warm]
  (let [dt (double dt)
        {:keys [bodies models friction restitution max-push-speed]} w
        beta (double (or (:convex-beta w) 0.1))
        tau (double (or (:convex-tau w) 1.0))
        e (double restitution)
        ;; Its own slop, and a small one. The slop is there to keep an
        ;; impulse solver from chasing an overlap that rounding puts back,
        ;; which this does not do; and a very stiff contact still gives a
        ;; little under load each step, so every contact in a stack ends
        ;; up the whole slop deep -- five millimeters a course, and a
        ;; column of thirty fifteen centimeters short, at the slop the
        ;; impulse solve wants.
        slop (double (or (:convex-slop w) 0.001))
        regularize (fn [wn] (/ (* beta beta (max (double wn) 1e-9)) (* 4.0 math/PI math/PI)))
        ;; blocks, made as the rows reach them
        blocks (volatile! {})
        ;; each model's inverse mass matrix, flat, as its block was made
        minvs (volatile! {})
        builder (volatile! (convex/builder))
        block! (fn [key]
                 (if-let [i (get @blocks key)]
                   i
                   (let [[kind i] key
                         [n mass minv vstar]
                         (case kind
                           :rigid (let [b (nth bodies i)]
                                    (when-not (rigid/inert? b)
                                      (let [[m mi] (rigid-block b) i3 (* 3 (long i))]
                                        [6 m mi (into (mapv #(aget vel (+ i3 %)) (range 3))
                                                      (mapv #(aget omega (+ i3 %)) (range 3)))])))
                           :model (let [{:keys [model pose sleeping?]} (nth models i)]
                                    (when-not sleeping?
                                      (let [mass (ab/generalized-mass-matrix model (:q pose) (:ls (nth configs i)))
                                            ^doubles u (nth us i)
                                            n (alength u)
                                            minv (convex/spd-inverse mass n)]
                                        (vswap! minvs assoc i minv)
                                        [n (vec mass) (vec minv) (vec u)]))))]
                     (if-not n
                       (do (vswap! blocks assoc key -1) -1)
                       (let [[b' idx] (convex/add-block @builder n mass minv vstar)]
                         (vreset! builder b')
                         (vswap! blocks assoc key idx)
                         idx)))))
        ;; a prepared row as `[ka kb against-a against-b]`, the second side
        ;; pushed the other way
        row-of (fn [row]
                 (if-let [p (:pair row)]
                   (let [[key s] (side-row p)] [(block! key) -1 s nil])
                   (let [[ka sa] (side-row (:a row)) [kb sb] (side-row (:b row))
                         ia (if ka (block! ka) -1) ib (if kb (block! kb) -1)]
                     [ia ib (when (>= (long ia) 0) sa) (when (>= (long ib) 0) (mapv - sb))])))
        ;; a row's inverse effective mass, what its regularization is
        ;; scaled by: a model side's is g' M^-1 g, from the matrix the
        ;; block was made with, since its side was prepared without one
        side-w (fn [side]
                 (case (:kind side)
                   (:link :pair)
                   (let [mi (:mi side)]
                     (if (neg? (long (block! [:model mi])))
                       0.0
                       (let [^doubles g (:force side) ^doubles m (get @minvs mi) n (alength g)
                             ;; row r of M^-1 against g
                             mg (fn [r] (let [o (* (long r) n)]
                                          (loop [c 0 x 0.0]
                                            (if (== c n) x (recur (inc c) (+ x (* (aget m (+ o c)) (aget g c))))))))]
                         (loop [r 0 acc 0.0]
                           (if (== r n) acc (recur (inc r) (+ acc (* (aget g r) (double (mg r))))))))))
                   (double (or (:w side) 0.0))))
        row-w (fn [row] (if-let [p (:pair row)] (side-w p) (+ (double (side-w (:a row))) (double (side-w (:b row))))))
        add! (fn [kind rows rn rt vhat mu]
               (let [[ka kb] (first rows)]
                 (vswap! builder convex/add-constraint kind ka kb
                         (mapv (fn [[_ _ sa sb]] [sa sb]) rows) rn rt vhat mu)))
        live? (fn [[ka kb]] (or (>= (long ka) 0) (>= (long kb) 0)))
        ;; contacts, then limits, then pins' axes: each constraint's place,
        ;; or nil for one with nothing awake on either side
        contact-at (vec (for [[c row] (map vector cs prepared)]
                          (let [rows (mapv row-of (:dirs row))
                                normal (nth (:dirs row) 0)]
                            (when (live? (first rows))
                              (let [wn (double (row-w normal))
                                    over (- (double (:depth c)) slop)
                                    push (if (pos? over) (min (double max-push-speed) (/ over (* tau dt))) (/ over dt))
                                    va (double (speed-of normal vel omega us))
                                    vhat (if (and (pos? e) (< va -0.5)) (max push (* (- e) va)) push)]
                                (add! convex/contact rows (regularize wn) (* 1e-3 (max wn 1e-9)) vhat friction)
                                (dec (count (:kind @builder))))))))
        _ (doseq [{:keys [mi g m bias]} limits]
            (let [row [(block! [:model mi]) -1 (vec g) nil]]
              (when (live? row)
                (add! convex/unilateral [row] (regularize (if (pos? (double m)) (/ 1.0 (double m)) 0.0)) 0.0 bias 0.0))))
        pin-at (vec (for [{:keys [rows bias]} pins]
                      (vec (for [r (range 3)]
                             (let [row (nth rows r) built (row-of row) wr (double (row-w row))]
                               (when (and (live? built) (> wr 1e-12))
                                 (add! convex/bilateral [built] (regularize wr) 0.0 (nth bias r) 0.0)
                                 (dec (count (:kind @builder)))))))))
        problem (convex/build @builder)
        nk (long (:nk problem))
        ;; last step's impulses, for the warm start
        last-impulses (let [g (a/f64 (* 3 nk))]
                        (doseq [[k was] (map vector contact-at warm) :when (and k was)]
                          (dotimes [t 3] (aset g (+ (* 3 (long k)) t) (double (nth was t)))))
                        (doseq [[ks was] (map vector pin-at pin-warm) :when was
                                [r k] (map-indexed vector ks) :when k]
                          (aset g (* 3 (long k)) (double (nth was r))))
                        g)
        warm-v (convex/pushed problem last-impulses)
        v0 (when (and (:warm-start? w) (< (double (convex/cost problem warm-v))
                                          (double (convex/cost problem (:vstar problem)))))
             warm-v)
        placed @blocks
        kept (:convex-factors w)
        {:keys [^doubles v ^doubles impulses factors]}
        (convex/solve problem v0 (if (= placed (:blocks kept)) {:factors (:factors kept)} {}))
        ^ints boff (:boff problem)]
    ;; the velocities back where the rest of the step reads them
    (doseq [[[kind i] b] placed :when (>= (long b) 0)]
      (let [o (aget boff b)]
        (case kind
          :rigid (let [i3 (* 3 (long i))]
                   (dotimes [k 3]
                     (aset vel (+ i3 k) (aget v (+ o k)))
                     (aset omega (+ i3 k) (aget v (+ o 3 k)))))
          :model (let [^doubles u (nth us i)]
                   (dotimes [k (alength u)] (aset u k (aget v (+ o k))))))))
    {:contacts (mapv (fn [k] (if k
                               (let [k3 (* 3 (long k))]
                                 [(aget impulses k3) (aget impulses (+ k3 1)) (aget impulses (+ k3 2))])
                               [0.0 0.0 0.0]))
                     contact-at)
     :pins (mapv (fn [ks] (mapv (fn [k] (if k (aget impulses (* 3 (long k))) 0.0)) ks)) pin-at)
     :convex-factors {:blocks placed :factors factors}}))

;; ---------------------------------------------------------------------------
;; A step

(defn- contact-keys
  "What counts as the same contact as last step, for each of `cs`.

  The two sides and the feature id `allgo.physics.contact` stamps on the
  point, as `allgo.physics.solver` keys them -- and a count, because not
  every shape tells its points apart: a tire against the road calls all
  of its points `:torus-box`, and without the count they would share one
  remembered impulse."
  [cs]
  (first (reduce (fn [[out seen] {:keys [a b id]}]
                   (let [k [a b id]
                         n (long (get seen k 0))]
                     [(conj out (conj k n)) (assoc seen k (inc n))]))
                 [[] {}]
                 cs)))

;; ---------------------------------------------------------------------------
;; Sleeping
;;
;; As `allgo.physics.solver` sleeps its bodies, by islands -- everything
;; joined by a contact or a pin sleeps together or not at all -- with a
;; whole model as one member of an island. A model sleeps when its base and
;; every joint are slow, and never while anything drives it: a `:tau` with
;; anything in it is somebody holding the controls.
;;
;; Asleep, a body is left out of the step entirely and answers every
;; contact as the ground does. When an island goes to sleep its members
;; are stamped with an island id, and anything awake that touches one of
;; them wakes them all at once: nothing finds contacts between two
;; sleeping bodies, so there would be no other way to know that the brick
;; under the one struck was holding it up.

(defn- driven? [{:keys [tau]}]
  (boolean (some #(not (zero? (double %))) tau)))

(defn- who
  "What waking the body on `side` wakes: its island, or the body alone
  if it fell asleep outside one -- placed asleep with `rigid/asleep`."
  [{:keys [bodies models]} side]
  (let [x (if (rigid-side? side) (nth bodies (second side)) (nth models (second side)))]
    (cond (:island x) [:island (:island x)]
          (rigid-side? side) side
          :else [:link (second side) -1])))

(defn- wake-all
  "Wakes everything asleep that `disturbed` names: by island, or by side."
  [w disturbed]
  (let [hit? (fn [x side] (or (disturbed side) (some-> (:island x) (#(disturbed [:island %])))))]
    (assoc w
           :bodies (vec (map-indexed (fn [i b]
                                       (if (and (rigid/sleeping? b) (hit? b [:rigid i]))
                                         (dissoc (rigid/wake b) :island)
                                         b))
                                     (:bodies w)))
           :models (vec (map-indexed (fn [mi m]
                                       (if (and (:sleeping? m) (hit? m [:link mi -1]))
                                         (-> m (assoc :sleeping? false :sleep-time 0.0) (dissoc :island))
                                         m))
                                     (:models w))))))

(defn wake
  "`w` with whatever is on `side` awake, and the rest of its island.

  Anything that changes a body from outside the step -- a push, a new
  velocity, a new place -- has to wake it, or a sleeping body ignores
  it until something bumps into it. A side is `[:rigid i]` or
  `[:link model-index link-index]`, as for pins."
  [w side]
  (if (asleep? w side) (wake-all w #{(who w side)}) w))

(defn- rouse
  "Wakes whatever something awake has come to touch or is pinned to,
  with the rest of its island, and any sleeping model given a torque.
  Returns `[world woke?]`."
  [w cs]
  (let [disturbed (into #{}
                        (mapcat (fn [{:keys [a b]}]
                                  (cond
                                    (and (asleep? w a) (awake? w b)) [(who w a)]
                                    (and (asleep? w b) (awake? w a)) [(who w b)]
                                    :else [])))
                        (concat cs (:pins w)))
        disturbed (into disturbed
                        (keep-indexed (fn [mi m]
                                        (when (and (:sleeping? m) (driven? m))
                                          (who w [:link mi -1]))))
                        (:models w))]
    (if (empty? disturbed)
      [w false]
      [(wake-all w disturbed) true])))

(defn- still-model?
  "A model at rest: its base turning and moving slower than the
  thresholds, and every joint turning slower than the angular one."
  [{:keys [model pose]} lin ang]
  (let [u (ab/velocity model pose)
        base? (ab/base model)]
    (every? true? (map-indexed (fn [i x]
                                 (< (abs (double x)) (if (and base? (<= 3 i 5)) lin ang)))
                               u))))

(defn- settle
  "Runs each awake body's sleep clock, and puts to sleep every island
  whose members have all been still for `:time-to-sleep`."
  [{:keys [bodies models pins] :as w} cs dt]
  (let [dt (double dt)
        lin (double (:linear-sleep-speed w))
        ang (double (:angular-sleep-speed w))
        wait (double (:time-to-sleep w))
        tick (fn [x still?] (assoc x :sleep-time (if still? (+ (double (or (:sleep-time x) 0.0)) dt) 0.0)))
        bodies (mapv (fn [b]
                       (if (rigid/inert? b)
                         b
                         (tick b (and (< (v/length (:vel b)) lin) (< (v/length (:omega b)) ang)))))
                     bodies)
        models (mapv (fn [m]
                       (if (:sleeping? m)
                         m
                         (tick m (and (not (driven? m)) (still-model? m lin ang)))))
                     models)
        w (assoc w :bodies bodies :models models)
        nb (count bodies)
        node (fn [side] (when (awake? w side)
                          (if (rigid-side? side) (long (second side)) (+ nb (long (second side))))))
        clock (fn [^long i] (double (or (:sleep-time (if (< i nb) (nth bodies i) (nth models (- i nb)))) 0.0)))
        n (+ nb (count models))]
    ;; An island sleeps on the clock of its youngest member, so until
    ;; some body has been still long enough no island can.
    (if-not (some #(and (node %) (>= (clock (node %)) wait))
                  (concat (map (fn [i] [:rigid i]) (range nb)) (map (fn [i] [:link i -1]) (range (count models)))))
      w
      (let [^ints parent (a/i32 (range n))
            find (fn find [^long i]
                   (let [p (aget parent i)]
                     (if (== p i) i (let [r (find p)] (aset parent i (int r)) r))))]
        (doseq [{:keys [a b]} (concat cs pins)]
          (let [x (node a) y (node b)]
            (when (and x y)
              (let [rx (find x) ry (find y)]
                (when-not (== rx ry) (aset parent (int rx) (int ry)))))))
        (let [members (group-by find (filter #(node (if (< (long %) nb) [:rigid %] [:link (- (long %) nb) -1]))
                                             (range n)))
              sleepy (filter (fn [[_ ms]] (every? #(>= (clock %) wait) ms)) members)
              base (long (or (:next-island w) 0))]
          (reduce (fn [w [k [_ ms]]]
                    (let [id (+ base (long k))]
                      (reduce (fn [w i]
                                (if (< (long i) nb)
                                  (update-in w [:bodies i] assoc :sleeping? true :island id
                                             :vel v/zero :omega v/zero)
                                  (update-in w [:models (- (long i) nb)]
                                             (fn [{:keys [model] :as m}]
                                               (-> m
                                                   (assoc :sleeping? true :island id)
                                                   (update :pose #(ab/with-velocity model %
                                                                    (vec (repeat (count (ab/velocity model %)) 0.0)))))))))
                              w ms)))
                  (update w :next-island #(+ (long (or % 0)) (count sleepy)))
                  (map-indexed vector sleepy)))))))

(defn- sleepers
  "The entries of `m`, keyed by contact, whose two sides are both asleep
  or fixed: contacts not solved this step, whose memory is kept for when
  they wake."
  [w m]
  (into {} (filter (fn [[k _]] (not (or (awake? w (nth k 0)) (awake? w (nth k 1)))))) m))

(defn- recall
  "Last step's impulse for each of `cs`, keyed `ks`: by key, or, for a
  contact whose key is new, from the nearest point of the same pair that
  no contact claimed by key, within `reach` meters.

  A key is the name the box test gives a point, after the corners and
  cuts that made it, and a brick swaying over the one below it renames
  its points as its corners cross the edges beneath: a corner inside
  becomes a cut and back again, and the impulse it carried is lost just
  as the sway most needs it. The point itself has hardly moved. Matching
  by nearness has no edges for it to slide across, which rounding to a
  grid did."
  [cs ks remembered points reach]
  (let [reach (double reach)
        claimed (into #{} (filter #(contains? remembered %)) ks)
        spare (group-by #(subvec % 0 2) (remove claimed (keys remembered)))]
    (first
     (reduce (fn [[out taken] [c k]]
               (if (claimed k)
                 [(conj out (get remembered k)) taken]
                 (let [p (:point c)
                       best (reduce (fn [best o]
                                      (let [d (v/distance p (get points o))]
                                        (if (and (< d reach) (not (taken o))
                                                 (or (nil? best) (< d (double (second best)))))
                                          [o d]
                                          best)))
                                    nil
                                    (get spare (subvec k 0 2)))]
                   (if best
                     [(conj out (get remembered (first best))) (conj taken (first best))]
                     [(conj out nil) taken]))))
             [[] #{}]
             (map vector cs ks)))))

(defn step
  "One step of the whole scene.

  Velocities first for everybody, then one contact sweep over all of
  them together, then positions. The order is the same one the other
  two solvers use and for the same reason: a body is stopped before it
  moves into something, rather than pulled back out afterward."
  [w dt]
  (let [dt (double dt)
        ;; Where everything is does not change until the end of the step,
        ;; so neither do the joints' frames or the contacts.
        configs (mapv (fn [{:keys [model pose]}]
                        (ab/configuration model (:q pose) (:base pose)))
                      (:models w))
        ;; Waking comes first, so that whatever is woken is solved in the
        ;; same step that reached it. Nothing found contacts among the
        ;; sleepers, so once some wake their contacts are found again --
        ;; or a column woken at the top would fall through itself for a
        ;; step -- and those may wake more.
        [w cs] (loop [w w cs (contacts w configs)]
                 (let [[w' woke?] (rouse w cs)]
                   (if woke? (recur w' (contacts w' configs)) [w cs])))
        {:keys [bodies models gravity]} w
        ;; Gravity on the loose bodies. `rigid/integrate` moves them as
        ;; well, so only the velocity half is wanted here.
        bodies (mapv (fn [b]
                       (if (rigid/inert? b)
                         b
                         (update b :vel v/add-scaled gravity dt)))
                     bodies)
        ;; And the jointed ones, whose accelerations come from the whole
        ;; tree at once.
        models (mapv (fn [{:keys [model pose tau sleeping?] :as m}]
                       (if sleeping?
                         m
                         (let [root (ab/base model)
                               {:keys [qdd base-acc]}
                               (ab/forward-dynamics model (:q pose) (:qd pose)
                                                    (or tau (vec (repeat (ab/dof model) 0.0)))
                                                    (cond-> {:gravity gravity}
                                                      root (assoc :base (:base pose))))
                               accel (if root (vec (concat base-acc qdd)) (vec qdd))]
                           (assoc m :pose
                                  (ab/with-velocity
                                    model pose
                                    (mapv (fn [a c] (+ (double a) (* dt (double c))))
                                          (ab/velocity model pose) accel))))))
                     models)
        w (assoc w :bodies bodies :models models)
        ;; Two sleepers, or a sleeper and the ground, have nothing to
        ;; solve; what they pushed with is kept for when they wake.
        cs (filterv #(or (awake? w (:a %)) (awake? w (:b %))) cs)
        limits (prepare-limits w configs dt)
        pins (mapv #(prepare-pin w configs % dt w) (:pins w))
        [vel omega] (body-arrays bodies)
        us (mapv (fn [{:keys [model pose]}] (a/f64 (ab/velocity model pose))) models)
        warm? (:warm-start? w)
        ks (contact-keys cs)
        remembered (:contact-impulses w)
        solved (if (or (seq cs) (seq limits) (seq pins))
                 (let [prepared (let [opts (cond-> w (:soft? w) (assoc :softness (softness w dt)))]
                                  (mapv #(prepare w configs % dt opts) cs))
                       warm (if warm?
                              (if-let [reach (:recall-reach w)]
                                (recall cs ks remembered (:contact-points w) reach)
                                (mapv #(get remembered %) ks))
                              (mapv (constantly nil) ks))
                       pin-warm (mapv #(when warm? (:impulse %)) (:pins w))]
                   (if (= :convex (:solver w))
                     (solve-convex! w configs cs prepared limits pins vel omega us dt warm pin-warm)
                     (solve! prepared limits pins vel omega us w warm pin-warm)))
                 {:contacts [] :pins []})
        carried (mapv (fn [[x y z]]
                        (math/sqrt (+ (* (double x) (double x)) (* (double y) (double y))
                                      (* (double z) (double z)))))
                      (:pins solved))
        ;; A pin pulled harder than it can hold lets go, and stays gone.
        ;; With a `:break-time`, harder for that long: the force is
        ;; averaged over it first, and the average is what is held.
        loaded (mapv (fn [pin impulse]
                       (let [f (/ (double impulse) (double dt))
                             tau (double (or (:break-time pin) 0.0))]
                         (if (pos? tau)
                           (let [load (double (or (:load pin) 0.0))]
                             (assoc pin :load (+ load (* (- f load) (min 1.0 (/ (double dt) tau))))))
                           (assoc pin :load f))))
                     (map (fn [pin impulse] (assoc pin :impulse impulse))
                          (:pins w) (:pins solved))
                     carried)
        held? (mapv (fn [pin]
                      (let [limit (:break-force pin)]
                        (or (nil? limit) (<= (double (:load pin)) (double limit)))))
                    loaded)
        w (assoc w :pins loaded)
        w (if (every? true? held?)
            w
            (assoc w
                   :pins (vec (keep-indexed (fn [i p] (when (nth held? i) p)) (:pins w)))
                   :broken-pins (into (vec (:broken-pins w))
                                      (keep-indexed (fn [i p] (when-not (nth held? i) p))
                                                    (:pins w)))))
        bodies (write-back bodies vel omega)
        models (mapv (fn [m i]
                       (update m :pose
                               #(ab/with-velocity (:model m) % (vec (seq ^doubles (nth us i))))))
                     models
                     (range (count models)))]
    (cond->
     (assoc w
            :pin-impulses carried
            :convex-factors (:convex-factors solved)
            :contact-impulses (into (sleepers w remembered) (zipmap ks (:contacts solved)))
            :contact-points (into (sleepers w (:contact-points w)) (zipmap ks (map :point cs)))
            ;; Positions last, from whatever the sweep left behind.
            :bodies (mapv (fn [b]
                            (if (rigid/inert? b)
                              b
                              (rigid/integrate b dt [0.0 0.0 0.0])))
                          bodies)
            ;; Only moved, not accelerated again: everything that was
            ;; going to change a velocity has already changed it.
            :models (mapv (fn [{:keys [model pose sleeping?] :as m}]
                            (if sleeping? m (assoc m :pose (ab/advance pose model dt))))
                          models))
      (:allow-sleep? w) (settle cs dt))))
