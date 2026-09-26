(ns allgo.physics.solver
  "Three ways to make contacts hold, in one engine, so they can be
  compared on the same pile of bricks.

  `allgo.physics.contact` says where bodies touch. What to *do* about it
  is a separate question with more than one respectable answer, and the
  three here are the ones in use:

    :sequential-impulse   the classic. Work in velocities: for each
                          contact compute the impulse that would stop the
                          two surfaces approaching, apply it, move on, and
                          go round again. Erin Catto's Box2D and every
                          engine descended from it.
    :tgs                  the same impulses, but the bodies move *during*
                          the solve. Split the step into substeps and
                          advance positions between them, so each round
                          is solved against where the bodies are rather
                          than where they were. PhysX's default.
    :xpbd                 do not solve velocities at all. Move the bodies
                          until they no longer overlap, then read the
                          velocity back off how far each one traveled.
                          `allgo.physics.rigid` and Müller's papers.

  ## What actually separates them

  All three are Gauss-Seidel: they visit constraints one at a time and
  each one sees the effect of the last. The differences are *when the
  geometry is evaluated* and *what quantity is being corrected*.

  Sequential impulse linearizes once, at the top of the step, and every
  iteration after that works against stale lever arms. That is fine when
  bodies barely move in a step and wrong when they do -- a fast spin, or
  a deep stack settling a long way. It is also the cheapest per
  iteration, because nothing has to be recomputed.

  TGS re-evaluates. Anchors are kept in each body's own frame, so after
  the bodies move the contact's world position and separation can be
  recovered without running collision detection again. That is the whole
  trick, and it is why TGS holds a tall stack that sequential impulse
  lets sag: the twentieth iteration is solving the real problem rather
  than a linearization of the problem as it was twenty iterations ago.

  XPBD corrects positions instead of velocities, which cannot inject
  energy -- the worst a bad correction can do is put a body somewhere and
  have the next one move it back. It pays for that by needing the
  velocity pass afterward to get restitution and friction, and by being
  stiffer to tune.

  ## Which one to use

  Measured, on the demo's own wall -- a running bond four bricks wide,
  gripping at 0.55, left alone for twenty seconds. The count is the
  bricks still where they were laid; one short of the total is a full
  marks, because every course of a running bond ends in a brick
  overhanging by exactly half, and that one falls off on its own.

      courses            4      6      8     12     16
      sequential impulse 15/16  23/24  31/32  10/48   0/64
      tgs                15/16  23/24  31/32  47/48  63/64
      xpbd               15/16  23/24  31/32   2/48   0/64

  **Up to about eight courses, take any of them.** They all hold, and
  the choice is about character rather than whether the wall stands:
  sequential impulse is the cheapest step, XPBD settles the deadest,
  TGS is the most controllable on impact.

  **Past eight courses, use TGS.** It is the only one of the three that
  substeps, and holding a tall stack is exactly what substepping buys:
  the bodies move between iterations and the geometry is measured again,
  so support propagates from the floor to the top *within* the step.

  Sequential impulse cannot get there and more iterations will not take
  it -- sixteen and sixty-four both leave nineteen bricks of a sixteen
  course wall standing, which is the same answer at four times the cost.
  It linearizes once, at the top of the step, so its later iterations are
  solving lever arms and overlaps measured before anything moved. That is
  not a bug to be fixed; it is what sequential impulse *is*, and it is
  why TGS exists.

  XPBD stops at about the same height and for a different reason. It
  reads velocity back off the position correction, at the substep rate,
  so a correction that does not settle becomes a velocity that does not
  settle. Past what it can hold, it now comes down rather than taking
  off -- the correction is bounded by `max-push-speed`, the same bound
  the impulse solvers have always had -- but it does come down.

  So: **tall stacks want TGS**. Shallow piles, ragdolls, anything thrown
  around, and anything where a scene settling dead matters more than a
  scene standing tall, can have whichever suits.

  ## The state, and why it is in arrays

  Bodies come in and go out as the maps `allgo.physics.rigid` makes. In
  between, the solve runs on flat arrays: velocity, angular velocity, and
  the world-space inverse inertia tensor, one row per body. A wall of two
  hundred bricks is six hundred contacts, and eight iterations over them
  is five thousand contact solves a step -- persistent vectors lose that
  by an order of magnitude, and the whole ask here was that it be fast.

  Warm starting is the other reason contacts are arrays: the impulse a
  contact needed last step is a very good guess at what it needs this
  step, and applying it before the first iteration is most of what lets a
  stack stand at a handful of iterations instead of fifty."
  (:require [allgo.array :as a]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.contact :as contact]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.toi :as toi]
            [clojure.math :as math]))

(def defaults
  {:solver :sequential-impulse
   :gravity [0.0 -10.0 0.0]
   :iterations 8
   ;; Only :tgs and :xpbd use these; sequential impulse takes the step
   ;; whole, which is the point of comparing them.
   :substeps 4
   :friction 0.5
   :restitution 0.0
   ;; Contacts shallower than this are left alone. Without a little slop
   ;; a resting stack jitters forever, chasing a penetration that the
   ;; next step's rounding puts straight back.
   :slop 0.005
   ;; Contacts are springs, not rules. A rigid contact solved with a
   ;; Baumgarte bias pushes overlapping bodies apart at `bias/dt` times
   ;; the overlap, which for a stack too deep to converge in the
   ;; iterations it has is enough to throw the bricks off each other: they
   ;; leave, they lose the contact and its remembered impulse, they fall
   ;; back deeper, and it pushes harder. A twenty brick column pogoed
   ;; itself apart in three seconds that way.
   ;;
   ;; A spring of a stated frequency and damping cannot do that. The
   ;; pushout is bounded by construction, the damping ratio says how it
   ;; comes to rest, and the constraint forgets a little of its
   ;; accumulated impulse each iteration -- which is what makes it unable
   ;; to pump. Erin Catto's soft constraints; Box2D v3 and PhysX's TGS
   ;; Soft are both this.
   :contact-hertz 30.0
   :contact-damping 10.0
   ;; However deep the overlap, it is never pushed out faster than this.
   :max-push-speed 3.0
   ;; Below this closing speed a contact is treated as resting and does
   ;; not bounce, however elastic the material.
   :restitution-threshold 1.0
   :warm-start? true
   ;; Sleeping. A body still for `time-to-sleep` stops being simulated at
   ;; all, and wakes when something awake comes to touch it.
   :allow-sleep? true
   :time-to-sleep 0.5
   :linear-sleep-speed 0.01
   :angular-sleep-speed 0.035})

(defn- awake?
  "A body the step should move: not static, not asleep."
  [b]
  (not (rigid/inert? b)))

(defn- dynamic?
  "Anything that can move, whether or not it is moving now.

  Which is what islands are built from: a sleeping brick is still part of
  the group it settled with, and has to wake with it."
  [b]
  (not (rigid/static? b)))

;; ---------------------------------------------------------------------------
;; The solve's own view of the bodies

(defn- inverse-inertia-world
  "`R * diag(invI) * R^T`, written into `out` at row `i`.

  The body's inertia is diagonal in its own frame and is not diagonal in
  the world's, so this is needed once per body per step and then used by
  every contact on it."
  [^doubles out ^long i rot inv-inertia]
  (let [[[m00 m01 m02] [m10 m11 m12] [m20 m21 m22]] (q/to-matrix rot)
        [ix iy iz] inv-inertia
        ix (double ix) iy (double iy) iz (double iz)
        b (* i 9)]
    ;; Columns of R scaled by the diagonal, times R^T.
    (aset out b (+ (* m00 m00 ix) (* m01 m01 iy) (* m02 m02 iz)))
    (aset out (+ b 1) (+ (* m00 m10 ix) (* m01 m11 iy) (* m02 m12 iz)))
    (aset out (+ b 2) (+ (* m00 m20 ix) (* m01 m21 iy) (* m02 m22 iz)))
    (aset out (+ b 3) (aget out (+ b 1)))
    (aset out (+ b 4) (+ (* m10 m10 ix) (* m11 m11 iy) (* m12 m12 iz)))
    (aset out (+ b 5) (+ (* m10 m20 ix) (* m11 m21 iy) (* m12 m22 iz)))
    (aset out (+ b 6) (aget out (+ b 2)))
    (aset out (+ b 7) (aget out (+ b 5)))
    (aset out (+ b 8) (+ (* m20 m20 ix) (* m21 m21 iy) (* m22 m22 iz)))
    out))

(defn- tangents
  "Two directions perpendicular to `n` and to each other.

  Friction acts in the contact plane, and the plane needs a basis. Which
  basis does not matter as long as it is consistent within a step."
  [[nx ny nz]]
  (let [t (if (< (abs (double nx)) 0.57735)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0])
        t1 (v/normalize (v/cross [nx ny nz] t))]
    [t1 (v/cross [nx ny nz] t1)]))

;; ---------------------------------------------------------------------------
;; Impulses

(defn- ang-vel [^doubles w ^long i]
  [(aget w (* i 3)) (aget w (+ (* i 3) 1)) (aget w (+ (* i 3) 2))])

(defn- lin-vel [^doubles vel ^long i]
  [(aget vel (* i 3)) (aget vel (+ (* i 3) 1)) (aget vel (+ (* i 3) 2))])

(defn- mul-inertia [^doubles ii ^long i [x y z]]
  (let [b (* i 9)]
    [(+ (* (aget ii b) x) (* (aget ii (+ b 1)) y) (* (aget ii (+ b 2)) z))
     (+ (* (aget ii (+ b 3)) x) (* (aget ii (+ b 4)) y) (* (aget ii (+ b 5)) z))
     (+ (* (aget ii (+ b 6)) x) (* (aget ii (+ b 7)) y) (* (aget ii (+ b 8)) z))]))

(defn- apply-impulse!
  "Adds impulse `p` at lever arm `r` to body `i`."
  ;; No primitive hints on `i` or `sign`: a fn taking primitives is
  ;; limited to four arguments, and array hints do not count.
  [^doubles vel ^doubles omega ^doubles ii ^doubles inv-mass i r p sign]
  (let [i (long i)
        im (aget inv-mass i)]
    (when (pos? im)
      (let [s (double sign)
            [px py pz] p
            b3 (* i 3)]
        (aset vel b3 (+ (aget vel b3) (* s px im)))
        (aset vel (+ b3 1) (+ (aget vel (+ b3 1)) (* s py im)))
        (aset vel (+ b3 2) (+ (aget vel (+ b3 2)) (* s pz im)))
        (let [[dx dy dz] (mul-inertia ii i (v/cross r (v/scale p s)))]
          (aset omega b3 (+ (aget omega b3) dx))
          (aset omega (+ b3 1) (+ (aget omega (+ b3 1)) dy))
          (aset omega (+ b3 2) (+ (aget omega (+ b3 2)) dz)))))))

(defn- effective-mass
  "The inverse of how much the pair gives to an impulse along `n`."
  [^doubles ii ^doubles inv-mass ia ib ra rb n]
  (let [ia (long ia)
        ib (long ib)
        k (+ (aget inv-mass ia) (aget inv-mass ib)
             (v/dot n (v/cross (mul-inertia ii ia (v/cross ra n)) ra))
             (v/dot n (v/cross (mul-inertia ii ib (v/cross rb n)) rb)))]
    (if (> k 1e-12) (/ 1.0 k) 0.0)))

(defn- relative-velocity
  [^doubles vel ^doubles omega ia ib ra rb]
  (let [ia (long ia) ib (long ib)]
    (v/sub (v/add (lin-vel vel ib) (v/cross (ang-vel omega ib) rb))
           (v/add (lin-vel vel ia) (v/cross (ang-vel omega ia) ra)))))

;; ---------------------------------------------------------------------------
;; The world

(defn world
  "A simulation over `bodies`, which are `allgo.physics.rigid` maps.

  `:broad` is the broad phase's kept state, and the reason the world is
  threaded through `step` rather than rebuilt: sweep and prune is cheap
  because the order it sorted the bodies into last step is nearly right
  this step, and that is only true if it is the same one."
  ([bodies] (world bodies {}))
  ([bodies opts]
   (merge defaults opts {:bodies (vec bodies)
                         :contacts []
                         :broad (contact/broad-phase (count bodies))})))

;; ---------------------------------------------------------------------------
;; Pose arithmetic on the flat arrays
;;
;; The same operations `allgo.geometry.quaternion` and
;; `allgo.physics.rigid` already have, written to read and write array
;; slots instead of returning vectors. They are here rather than there
;; because there is where the readable version belongs: a body is a map
;; and a rotation is four numbers, and that is the right trade for
;; everything except the innermost loop of a solve that runs two and a
;; half thousand times a pass.
;;
;; Four arguments and one primitive hint apiece is not a style: a fn
;; taking primitives is limited to four arguments, and array hints do
;; not count toward it. Vectors come in and out through `^doubles`
;; scratch of length three, allocated once by the caller.

(defn- qrot!
  "Turns the vector in `v` by the quaternion at `i`, into `out`.

  `q/rotate`'s two cross products, in place. `out` may be `v`."
  [^doubles quats ^long i ^doubles v ^doubles out]
  (let [b (* i 4)
        qx (aget quats b) qy (aget quats (+ b 1))
        qz (aget quats (+ b 2)) qw (aget quats (+ b 3))
        vx (aget v 0) vy (aget v 1) vz (aget v 2)
        tx (* 2.0 (- (* qy vz) (* qz vy)))
        ty (* 2.0 (- (* qz vx) (* qx vz)))
        tz (* 2.0 (- (* qx vy) (* qy vx)))]
    (aset out 0 (+ vx (* qw tx) (- (* qy tz) (* qz ty))))
    (aset out 1 (+ vy (* qw ty) (- (* qz tx) (* qx tz))))
    (aset out 2 (+ vz (* qw tz) (- (* qx ty) (* qy tx))))
    out))

(defn- inertia-mul!
  "`invI * v` in the body's own frame, where the inertia is diagonal."
  [^doubles inv-inertia ^long i ^doubles v ^doubles out]
  (let [b (* i 3)]
    (aset out 0 (* (aget inv-inertia b) (aget v 0)))
    (aset out 1 (* (aget inv-inertia (+ b 1)) (aget v 1)))
    (aset out 2 (* (aget inv-inertia (+ b 2)) (aget v 2)))
    out))

(defn- angular-mass
  "How much the body at `i` gives to a correction along `n` acting at the
  lever arm in `r`.

  `rigid/inverse-mass`'s angular half: the lever arm crossed with the
  direction, taken into the body's frame, weighted by the inertia
  there. `scratch` is written over."
  ^double [arrays ^long i ^doubles r-cross-n ^doubles scratch]
  (let [^doubles inv-rot (:inv-rot arrays)
        ^doubles inv-inertia (:inv-inertia arrays)
        _ (qrot! inv-rot i r-cross-n scratch)
        b (* i 3)
        rx (aget scratch 0) ry (aget scratch 1) rz (aget scratch 2)]
    (+ (* rx rx (aget inv-inertia b))
       (* ry ry (aget inv-inertia (+ b 1)))
       (* rz rz (aget inv-inertia (+ b 2))))))

(defn- spin!
  "Turns the body at `i` by the small rotation in `w`.

  `rot += 1/2 * (w, 0) * rot`, renormalized, with the inverse
  orientation kept in step -- `rigid/with-rotation` and the angular half
  of `rigid/apply-impulse` together. For a unit quaternion the inverse
  *is* the conjugate, which is what `q/inverse` reduces to once the
  normalize above it has run."
  [^doubles rot ^doubles inv-rot ^long i ^doubles w]
  (let [b (* i 4)
        rx (aget rot b) ry (aget rot (+ b 1))
        rz (aget rot (+ b 2)) rw (aget rot (+ b 3))
        ox (aget w 0) oy (aget w 1) oz (aget w 2)
        ;; (ox oy oz 0) * rot
        dx (+ (* ox rw) (* oy rz) (- (* oz ry)))
        dy (+ (- (* ox rz)) (* oy rw) (* oz rx))
        dz (+ (* ox ry) (- (* oy rx)) (* oz rw))
        dw (- (+ (* ox rx) (* oy ry) (* oz rz)))
        nx (+ rx (* 0.5 dx)) ny (+ ry (* 0.5 dy))
        nz (+ rz (* 0.5 dz)) nw (+ rw (* 0.5 dw))
        l (math/sqrt (+ (* nx nx) (* ny ny) (* nz nz) (* nw nw)))]
    (if (zero? l)
      (do (aset rot b 0.0) (aset rot (+ b 1) 0.0)
          (aset rot (+ b 2) 0.0) (aset rot (+ b 3) 1.0)
          (aset inv-rot b 0.0) (aset inv-rot (+ b 1) 0.0)
          (aset inv-rot (+ b 2) 0.0) (aset inv-rot (+ b 3) 1.0))
      (let [ux (/ nx l) uy (/ ny l) uz (/ nz l) uw (/ nw l)]
        (aset rot b ux) (aset rot (+ b 1) uy)
        (aset rot (+ b 2) uz) (aset rot (+ b 3) uw)
        (aset inv-rot b (- ux)) (aset inv-rot (+ b 1) (- uy))
        (aset inv-rot (+ b 2) (- uz)) (aset inv-rot (+ b 3) uw)))))

(defn- anchor!
  "Where the contact anchor at `k` on the body at `i` now is, into `out`.

  The point is kept in the body's frame, so this is `rigid/local->world`
  and is what recovers the contact geometry after the bodies have moved
  without running collision detection again."
  [arrays ^long i ^doubles local ^doubles out]
  (let [^doubles rot (:rot arrays)
        ^doubles pos (:pos arrays)
        i3 (* i 3)]
    (qrot! rot i local out)
    (aset out 0 (+ (aget out 0) (aget pos i3)))
    (aset out 1 (+ (aget out 1) (aget pos (+ i3 1))))
    (aset out 2 (+ (aget out 2) (aget pos (+ i3 2))))
    out))

(defn- load3!
  "Copies three numbers out of a packed array at row `k`."
  [^doubles src ^long k ^doubles out]
  (let [k3 (* k 3)]
    (aset out 0 (aget src k3))
    (aset out 1 (aget src (+ k3 1)))
    (aset out 2 (aget src (+ k3 2)))
    out))

;; ---------------------------------------------------------------------------
;; The contacts' own arrays

(defn- contact-key
  "What counts as the same contact as last step, for warm starting.

  The pair, and the feature id `allgo.physics.contact` stamps on the
  point -- which corner of which face is pressed into which face. It used
  to be the contact point rounded to two centimeters, and that is a fine
  key for a stack that is already still and a bad one for a stack that is
  moving: the points slide across the rounding and the match is lost
  exactly when the impulse history is most needed. A twenty brick column
  matched a quarter of its contacts that way and fell over; on ids it
  matches all of them."
  [c]
  [(:a c) (:b c) (:id c)])

(defn- prepare
  "Turns a frame's contacts into the solver's own arrays.

  Anchors are stored in each body's frame as well as the world's. TGS and
  XPBD move the bodies mid-solve and need to know where the contact went;
  sequential impulse does not and ignores them.

  One pass, writing into the arrays directly. It used to be a dozen
  `mapcat`s over the contacts, one per array: a dozen lazy sequences of
  boxed doubles, `tangents` called twice per contact for the two halves
  of the basis it computes in one go, and a `nth` into the body vector
  four times over. That came to 5.4ms of a 32ms step on a 16x16 wall,
  which is more than the constraint solve it was preparing for."
  [arrays contacts previous]
  (let [contacts (vec contacts)
        n (count contacts)
        ;; Last step's impulses for the same contact, which is the guess
        ;; warm starting rests on.
        old (reduce (fn [m c] (assoc m (contact-key c) c)) {} previous)
        ^doubles pos (:pos arrays)
        ^doubles inv-rot (:inv-rot arrays)
        ^ints a-arr (a/i32 n) ^ints b-arr (a/i32 n)
        normal (a/f64 (* n 3)) t1-arr (a/f64 (* n 3)) t2-arr (a/f64 (* n 3))
        depth (a/f64 n) depth0 (a/f64 n)
        ra (a/f64 (* n 3)) rb (a/f64 (* n 3))
        ral (a/f64 (* n 3)) rbl (a/f64 (* n 3))
        pn (a/f64 n) p1 (a/f64 n) p2 (a/f64 n)
        t0 (a/f64 3) s0 (a/f64 3)]
    (dotimes [k n]
      (let [c (nth contacts k)
            ia (long (:a c)) ib (long (:b c))
            nrm (:normal c)
            [nx ny nz] nrm
            [qx qy qz] (:point c)
            ;; One call, both directions. Friction needs a basis for the
            ;; contact plane and `tangents` builds the whole of it.
            [[ax ay az] [bx by bz]] (tangents nrm)
            k3 (* k 3) ia3 (* ia 3) ib3 (* ib 3)]
        (aset a-arr k (int ia))
        (aset b-arr k (int ib))
        (aset normal k3 (double nx))
        (aset normal (+ k3 1) (double ny))
        (aset normal (+ k3 2) (double nz))
        (aset t1-arr k3 (double ax))
        (aset t1-arr (+ k3 1) (double ay))
        (aset t1-arr (+ k3 2) (double az))
        (aset t2-arr k3 (double bx))
        (aset t2-arr (+ k3 1) (double by))
        (aset t2-arr (+ k3 2) (double bz))
        (aset depth k (double (:depth c)))
        ;; The overlap as collision detection found it, kept so that
        ;; `refresh-anchors!` can add the drift to it rather than to a
        ;; value the last substep already moved.
        (aset depth0 k (double (:depth c)))
        ;; The lever arm in the world is the point less the center, and
        ;; the anchor in the body's own frame is that turned back by the
        ;; body's orientation -- so the world arm is what `world->local`
        ;; wanted anyway and is computed once for both.
        (let [rx (- (double qx) (aget pos ia3))
              ry (- (double qy) (aget pos (+ ia3 1)))
              rz (- (double qz) (aget pos (+ ia3 2)))]
          (aset ra k3 rx) (aset ra (+ k3 1) ry) (aset ra (+ k3 2) rz)
          (aset t0 0 rx) (aset t0 1 ry) (aset t0 2 rz)
          (qrot! inv-rot ia t0 s0)
          (aset ral k3 (aget s0 0))
          (aset ral (+ k3 1) (aget s0 1))
          (aset ral (+ k3 2) (aget s0 2)))
        (let [rx (- (double qx) (aget pos ib3))
              ry (- (double qy) (aget pos (+ ib3 1)))
              rz (- (double qz) (aget pos (+ ib3 2)))]
          (aset rb k3 rx) (aset rb (+ k3 1) ry) (aset rb (+ k3 2) rz)
          (aset t0 0 rx) (aset t0 1 ry) (aset t0 2 rz)
          (qrot! inv-rot ib t0 s0)
          (aset rbl k3 (aget s0 0))
          (aset rbl (+ k3 1) (aget s0 1))
          (aset rbl (+ k3 2) (aget s0 2)))
        (when-let [was (get old (contact-key c))]
          (aset pn k (double (:pn was)))
          ;; The tangential impulses come back too, now that the match is
          ;; exact. They were being dropped -- stored by `contact-state`
          ;; every step and never read -- so friction began each step
          ;; from nothing and had to be rediscovered in the iterations it
          ;; had left, which is what let a stack creep sideways while it
          ;; stood.
          (aset p1 k (double (:p1 was)))
          (aset p2 k (double (:p2 was))))))
    {:n n
     :a a-arr :b b-arr
     :normal normal :t1 t1-arr :t2 t2-arr
     :depth depth :depth0 depth0
     :ra ra :rb rb :ra-local ral :rb-local rbl
     :pn pn :p1 p1 :p2 p2
     ;; What the XPBD position solve had to push with to keep this contact
     ;; apart, summed over the last substep. Only that solver fills it,
     ;; and it is what bounds friction there.
     :lambda (a/f64 n)
     ;; The closing speed as the step began, which is what restitution is
     ;; measured against -- after an iteration or two it is gone.
     :approach (a/f64 n)
     :contacts contacts}))

(defn- body-arrays
  "The bodies as flat arrays, in one pass over the maps.

  Velocity and world inertia are what the impulse solvers work in. The
  pose is here too -- position, orientation, and the inverse orientation
  that takes a world direction into the body's frame -- because the
  position solve moves the bodies *during* the solve, and doing that
  through the body vector means rebuilding a 257 element vector of maps
  once per contact. On a sixteen course wall that was twenty thousand
  rebuilds a substep and half of XPBD's step.

  `inv-inertia` is the body frame diagonal rather than the world tensor
  in `ii`, because the position solve turns the bodies as it goes and a
  world tensor built at the top of the substep would be stale by the
  second contact.

  `movable` is what a correction may touch: anything not static. Not
  `awake?` -- a sleeping body is still pushed out of an overlap, and
  only its velocity read-back is skipped."
  [bodies]
  (let [n (count bodies)
        vel (a/f64 (* n 3))
        omega (a/f64 (* n 3))
        inv-mass (a/f64 n)
        ii (a/f64 (* n 9))
        pos (a/f64 (* n 3))
        rot (a/f64 (* n 4))
        inv-rot (a/f64 (* n 4))
        inv-inertia (a/f64 (* n 3))
        ^ints movable (a/i32 n)]
    (dotimes [i n]
      (let [b (nth bodies i)
            [vx vy vz] (:vel b)
            [wx wy wz] (:omega b)
            [px py pz] (:pos b)
            [rx ry rz rw] (:rot b)
            [qx qy qz qw] (:inv-rot b)
            [ax ay az] (:inv-inertia b)
            i3 (* i 3)
            i4 (* i 4)]
        (aset vel i3 (double vx))
        (aset vel (+ i3 1) (double vy))
        (aset vel (+ i3 2) (double vz))
        (aset omega i3 (double wx))
        (aset omega (+ i3 1) (double wy))
        (aset omega (+ i3 2) (double wz))
        (aset pos i3 (double px))
        (aset pos (+ i3 1) (double py))
        (aset pos (+ i3 2) (double pz))
        (aset rot i4 (double rx))
        (aset rot (+ i4 1) (double ry))
        (aset rot (+ i4 2) (double rz))
        (aset rot (+ i4 3) (double rw))
        (aset inv-rot i4 (double qx))
        (aset inv-rot (+ i4 1) (double qy))
        (aset inv-rot (+ i4 2) (double qz))
        (aset inv-rot (+ i4 3) (double qw))
        (aset inv-inertia i3 (double ax))
        (aset inv-inertia (+ i3 1) (double ay))
        (aset inv-inertia (+ i3 2) (double az))
        (aset inv-mass i (double (:inv-mass b)))
        (aset movable i (if (rigid/static? b) 0 1))
        (inverse-inertia-world ii i (:rot b) (:inv-inertia b))))
    {:n n :vel vel :omega omega :inv-mass inv-mass :ii ii
     :pos pos :rot rot :inv-rot inv-rot :inv-inertia inv-inertia
     :movable movable}))

(defn- write-poses
  "Puts the moved poses back on the bodies.

  `prev-pos` and `prev-rot` are left alone: XPBD reads its velocity off
  the difference between them and where the substep ended, so the pose
  the substep started from has to survive the solve."
  [bodies {:keys [pos rot inv-rot movable]}]
  (let [^doubles pos pos ^doubles rot rot ^doubles inv-rot inv-rot
        ^ints movable movable]
    (mapv (fn [i b]
            (let [i (long i)]
              (if (zero? (aget movable i))
                b
                (let [i3 (* i 3) i4 (* i 4)]
                  (assoc b
                         :pos [(aget pos i3) (aget pos (+ i3 1)) (aget pos (+ i3 2))]
                         :rot [(aget rot i4) (aget rot (+ i4 1))
                               (aget rot (+ i4 2)) (aget rot (+ i4 3))]
                         :inv-rot [(aget inv-rot i4) (aget inv-rot (+ i4 1))
                                   (aget inv-rot (+ i4 2)) (aget inv-rot (+ i4 3))])))))
          (range (count bodies))
          bodies)))

(defn- vec3-at [^doubles arr ^long i]
  [(aget arr (* i 3)) (aget arr (+ (* i 3) 1)) (aget arr (+ (* i 3) 2))])

;; ---------------------------------------------------------------------------
;; The velocity solve, shared by sequential impulse and TGS

(defn- solve-velocities!
  "One Gauss-Seidel sweep over the contacts.

  Written out in primitive doubles rather than through
  `allgo.geometry.vec3`, and that is not a style choice. Each contact
  solve is a couple of dozen small vector operations, and in a browser
  every one of them allocates a three-element persistent vector: a wall
  of seventy bricks came to 308 contacts, and eight sweeps over them took
  343ms -- 139 microseconds to do about a hundred floating point
  operations. The same arithmetic on locals is a fraction of that. The
  JVM's escape analysis hides most of the cost; JavaScript's does not,
  and this has to run in both."
  [{:keys [vel omega ii inv-mass]} cs dt opts]
  (let [{:keys [friction slop restitution]} opts
        ^doubles vel vel ^doubles omega omega ^doubles ii ii ^doubles inv-mass inv-mass
        ^doubles normal (:normal cs)
        ^doubles t1 (:t1 cs)
        ^doubles t2 (:t2 cs)
        ^doubles depth (:depth cs)
        ^doubles ra-arr (:ra cs)
        ^doubles rb-arr (:rb cs)
        ^doubles pn (:pn cs)
        ^doubles p1 (:p1 cs)
        ^doubles p2 (:p2 cs)
        ^doubles approach (:approach cs)
        ^ints ia-arr (:a cs)
        ^ints ib-arr (:b cs)
        dt (double dt)
        friction (double friction)
        slop (double slop)
        restitution (double restitution)
        threshold (double (:restitution-threshold opts))
        relax? (boolean (:relax? opts))
        bounce? (boolean (:bounce? opts))
        ;; The spring, as three numbers the solve can use directly. A
        ;; contact cannot be stiffer than the step can represent, so the
        ;; frequency is capped at a quarter of the step rate -- above that
        ;; the spring oscillates between steps instead of damping.
        hertz (min (double (:contact-hertz opts)) (* 0.25 (/ 1.0 dt)))
        zeta (double (:contact-damping opts))
        w (* 2.0 math/PI hertz)
        c (* dt w (+ (* 2.0 zeta) (* dt w)))
        bias-rate (/ w (+ (* 2.0 zeta) (* dt w)))
        ;; How much of the ideal impulse to apply, and how much of what
        ;; has accumulated to give back. The second is the whole trick:
        ;; an accumulated impulse that decays cannot store energy.
        mass-scale (/ c (+ 1.0 c))
        impulse-scale (/ 1.0 (+ 1.0 c))
        max-push (double (:max-push-speed opts))
        n-contacts (long (:n cs))]
    (dotimes [k n-contacts]
      (let [ia (aget ia-arr k) ib (aget ib-arr k)
            k3 (* 3 k)
            rax (aget ra-arr k3) ray (aget ra-arr (+ k3 1)) raz (aget ra-arr (+ k3 2))
            rbx (aget rb-arr k3) rby (aget rb-arr (+ k3 1)) rbz (aget rb-arr (+ k3 2))
            ia3 (* 3 ia) ib3 (* 3 ib) ia9 (* 9 ia) ib9 (* 9 ib)
            ima (aget inv-mass ia) imb (aget inv-mass ib)
            ;; Positive is a gap, negative is an overlap -- the sign
            ;; convention the speculative margin brings with it.
            sep (- slop (aget depth k))
            gap? (pos? sep)
            ;; A gap asks only that the two not close it faster than this
            ;; step can afford, so they meet at the surface instead of
            ;; inside it, and it asks that rigidly -- there is nothing
            ;; soft about arriving. An overlap is the spring.
            b-term (if gap?
                     (- (/ sep dt))
                     (min (* bias-rate (- sep)) max-push))
            ms (if gap? 1.0 mass-scale)
            is (if gap? 0.0 impulse-scale)
            va (aget approach k)
            ;; A gap only bounces if the surfaces actually reach each
            ;; other before the step is out. Without the test a ball
            ;; falling toward a floor it will not touch this step is
            ;; still handed a restitution target, and since the margin
            ;; is as wide as the body can travel there is always such a
            ;; step -- the ball bounces off nothing, short of the floor.
            arrives? (or (not gap?) (>= (* (- va) dt) sep))
            r-term (if (and arrives? (< va (- threshold)))
                     (* (- restitution) va)
                     0.0)
            ;; The relax pass asks for one thing only: that the surfaces
            ;; stop closing. No push, no bounce, nothing soft. It runs
            ;; after the bodies have moved and its job is to take back the
            ;; velocity the push put in, which is the cheap and standard
            ;; alternative to carrying that push in a pseudo velocity.
            ;; One inequality, four right-hand sides.
            ;;
            ;; A gap comes first and answers for every pass, which is
            ;; the part that is easy to get wrong. `Do not close this
            ;; gap faster than the step allows` is true of the solve,
            ;; of the relax that follows it, and of a gap that is never
            ;; reached at all -- it is a statement about geometry, not
            ;; about which pass is asking. Letting relax override it
            ;; with its usual zero stops a body dead the moment
            ;; anything comes within the margin, which with a margin as
            ;; wide as a step's travel means stopping in mid-air a
            ;; tenth of a meter short of the floor.
            ;;
            ;; Otherwise: the solve pushes overlap out and adds the
            ;; rebound, the relax pass takes back what the push added
            ;; and asks for nothing else, and the bounce pass asks for
            ;; the rebound alone once the other two have finished.
            normal-target (cond
                            bounce? r-term
                            gap? (if (pos? r-term) r-term b-term)
                            relax? 0.0
                            :else (+ b-term r-term))
            ms (if (or relax? bounce? gap?) 1.0 ms)
            is (if (or relax? bounce? gap?) 0.0 is)
            ;; The bounce pass is asked only of contacts that actually
            ;; pushed: `pn` above zero is the evidence that the
            ;; surfaces met rather than merely came near.
            skip? (and bounce? (or (zero? r-term) (<= (aget pn k) 0.0)))]
        ;; Three directions in turn -- the normal, then two tangents --
        ;; rather than a closure called three times. The closure was
        ;; allocated per contact per sweep, two and a half thousand of
        ;; them a step, and in a browser that cost more than the
        ;; arithmetic it wrapped.
        (loop [dir 0 limit 0.0]
          (when (and (not skip?) (< dir (if bounce? 1 3)))
            (let [^doubles darr (case dir 0 normal 1 t1 t2)
                  ^doubles acc (case dir 0 pn 1 p1 p2)
                  dx (aget darr k3) dy (aget darr (+ k3 1)) dz (aget darr (+ k3 2))
                  normal? (zero? dir)
                  lo (if normal? 0.0 (- limit))
                  hi (if normal? ##Inf limit)
                  target (if normal? normal-target 0.0)
                  rvx (- (+ (aget vel ib3)
                            (- (* (aget omega (+ ib3 1)) rbz) (* (aget omega (+ ib3 2)) rby)))
                         (+ (aget vel ia3)
                            (- (* (aget omega (+ ia3 1)) raz) (* (aget omega (+ ia3 2)) ray))))
                  rvy (- (+ (aget vel (+ ib3 1))
                            (- (* (aget omega (+ ib3 2)) rbx) (* (aget omega ib3) rbz)))
                         (+ (aget vel (+ ia3 1))
                            (- (* (aget omega (+ ia3 2)) rax) (* (aget omega ia3) raz))))
                  rvz (- (+ (aget vel (+ ib3 2))
                            (- (* (aget omega ib3) rby) (* (aget omega (+ ib3 1)) rbx)))
                         (+ (aget vel (+ ia3 2))
                            (- (* (aget omega ia3) ray) (* (aget omega (+ ia3 1)) rax))))
                  vd (+ (* rvx dx) (* rvy dy) (* rvz dz))
                  ;; r x d, through the inverse inertia, and back across r
                  ;; -- the first is the angular impulse, the pair is the
                  ;; angular share of the effective mass.
                  cax (- (* ray dz) (* raz dy))
                  cay (- (* raz dx) (* rax dz))
                  caz (- (* rax dy) (* ray dx))
                  iax (+ (* (aget ii ia9) cax) (* (aget ii (+ ia9 1)) cay) (* (aget ii (+ ia9 2)) caz))
                  iay (+ (* (aget ii (+ ia9 3)) cax) (* (aget ii (+ ia9 4)) cay) (* (aget ii (+ ia9 5)) caz))
                  iaz (+ (* (aget ii (+ ia9 6)) cax) (* (aget ii (+ ia9 7)) cay) (* (aget ii (+ ia9 8)) caz))
                  ka (+ (* dx (- (* iay raz) (* iaz ray)))
                        (* dy (- (* iaz rax) (* iax raz)))
                        (* dz (- (* iax ray) (* iay rax))))
                  cbx (- (* rby dz) (* rbz dy))
                  cby (- (* rbz dx) (* rbx dz))
                  cbz (- (* rbx dy) (* rby dx))
                  ibx (+ (* (aget ii ib9) cbx) (* (aget ii (+ ib9 1)) cby) (* (aget ii (+ ib9 2)) cbz))
                  iby (+ (* (aget ii (+ ib9 3)) cbx) (* (aget ii (+ ib9 4)) cby) (* (aget ii (+ ib9 5)) cbz))
                  ibz (+ (* (aget ii (+ ib9 6)) cbx) (* (aget ii (+ ib9 7)) cby) (* (aget ii (+ ib9 8)) cbz))
                  kb (+ (* dx (- (* iby rbz) (* ibz rby)))
                        (* dy (- (* ibz rbx) (* ibx rbz)))
                        (* dz (- (* ibx rby) (* iby rbx))))
                  keff (+ ima imb ka kb)
                  m (if (> keff 1e-12) (/ 1.0 keff) 0.0)
                  old (aget acc k)
                  ;; The normal is a spring and gives some of its
                  ;; accumulated impulse back each time; the tangents are
                  ;; rigid, and `mass-scale` 1 with `impulse-scale` 0 is
                  ;; exactly the hard solve they had before.
                  lambda (if normal?
                           (- (* m ms (- target vd)) (* is old))
                           (* m (- target vd)))
                  nw (min hi (max lo (+ old lambda)))
                  d (- nw old)]
              (aset acc k nw)
              (when (pos? ima)
                (aset vel ia3 (- (aget vel ia3) (* d dx ima)))
                (aset vel (+ ia3 1) (- (aget vel (+ ia3 1)) (* d dy ima)))
                (aset vel (+ ia3 2) (- (aget vel (+ ia3 2)) (* d dz ima)))
                (aset omega ia3 (- (aget omega ia3) (* d iax)))
                (aset omega (+ ia3 1) (- (aget omega (+ ia3 1)) (* d iay)))
                (aset omega (+ ia3 2) (- (aget omega (+ ia3 2)) (* d iaz))))
              (when (pos? imb)
                (aset vel ib3 (+ (aget vel ib3) (* d dx imb)))
                (aset vel (+ ib3 1) (+ (aget vel (+ ib3 1)) (* d dy imb)))
                (aset vel (+ ib3 2) (+ (aget vel (+ ib3 2)) (* d dz imb)))
                (aset omega ib3 (+ (aget omega ib3) (* d ibx)))
                (aset omega (+ ib3 1) (+ (aget omega (+ ib3 1)) (* d iby)))
                (aset omega (+ ib3 2) (+ (aget omega (+ ib3 2)) (* d ibz))))
              ;; Friction is bounded by the normal force it rides on, so
              ;; the normal pass hands its impulse to the two after it.
              (recur (inc dir) (if normal? (* friction nw) limit)))))))))

(defn- solve-xpbd-velocities!
  "The velocity pass XPBD needs after moving the bodies.

  Reading velocity back off the position change is what makes XPBD
  stable, and it is also what makes this pass necessary: the correction
  that pushed two bodies apart shows up in that reading as separating
  velocity, so a body that was resolved out of a deep overlap leaves the
  position solve traveling fast. A dropped box came off the floor at
  fifteen meters a second before this existed.

  So the normal velocity is *set* rather than pushed at: to whatever
  restitution asks of the approach speed measured before the solve, and
  to nothing at all when the material does not bounce. Unlike the
  impulse solvers there is no non-negative clamp, because the correction
  wanted here is usually negative -- it is taking energy the position
  solve put in."
  [{:keys [vel omega ii inv-mass]} cs h opts]
  (let [^doubles normal (:normal cs)
        ^doubles t1 (:t1 cs) ^doubles t2 (:t2 cs)
        ^doubles ra-arr (:ra cs) ^doubles rb-arr (:rb cs)
        ^doubles approach (:approach cs)
        ^doubles lambda (:lambda cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)
        h (double h)
        e (double (:restitution opts))
        threshold (double (:restitution-threshold opts))
        friction (double (:friction opts))]
    (dotimes [k (long (:n cs))]
      ;; Only contacts the position solve actually pushed on. A pair that
      ;; collision detection reported at the top of the frame but that is
      ;; apart by the time this runs has nothing to correct -- and since
      ;; the normal velocity here is *set* rather than pushed at, touching
      ;; one would brake a body that is legitimately flying away from it.
      ;; That is what killed the bounce: the substep after the impact, the
      ;; ball was leaving at four meters a second and this pass set it to
      ;; nothing.
      (when (pos? (aget lambda k))
        (let [ia (aget ia-arr k) ib (aget ib-arr k)
              ra (vec3-at ra-arr k) rb (vec3-at rb-arr k)
              n (vec3-at normal k)
              pre (aget approach k)
              rv (relative-velocity vel omega ia ib ra rb)
              vn (v/dot rv n)
              target (if (< pre (- threshold)) (* (- e) pre) 0.0)
              dv (- target vn)
              m (effective-mass ii inv-mass ia ib ra rb n)
              p (v/scale n (* m dv))]
          (apply-impulse! vel omega ii inv-mass ia ra p -1.0)
          (apply-impulse! vel omega ii inv-mass ib rb p 1.0)
          ;; Friction against the normal force that actually held this
          ;; contact -- the multiplier the position solve needed, over the
          ;; substep: `mu * lambda / h`.
          ;;
          ;; It used to be bounded by the approach speed instead, and that
          ;; is zero for anything at rest, so a settled stack had no
          ;; friction at all. Nothing then removed the sideways and
          ;; angular velocity the position solve hands back -- it reads
          ;; velocity off the correction at the substep rate, so a
          ;; millimeter of pushout is a quarter of a meter a second -- and
          ;; a six brick column wound itself from 6mm/s to 29m/s in five
          ;; seconds.
          (let [limit (* friction (/ (abs (aget lambda k)) h))
                rv (relative-velocity vel omega ia ib ra rb)]
            (doseq [tarr [t1 t2]]
              (let [t (vec3-at ^doubles tarr k)
                    vt (v/dot rv t)
                    mt (effective-mass ii inv-mass ia ib ra rb t)
                    j (max (- limit) (min limit (* mt (- vt))))
                    p (v/scale t j)]
                (apply-impulse! vel omega ii inv-mass ia ra p -1.0)
                (apply-impulse! vel omega ii inv-mass ib rb p 1.0)))))))))

(defn- warm-start! [{:keys [vel omega ii inv-mass]} cs]
  (let [^doubles normal (:normal cs)
        ^doubles t1 (:t1 cs) ^doubles t2 (:t2 cs)
        ^doubles ra-arr (:ra cs) ^doubles rb-arr (:rb cs)
        ^doubles pn (:pn cs) ^doubles p1 (:p1 cs) ^doubles p2 (:p2 cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)]
    (dotimes [k (long (:n cs))]
      (let [ia (aget ia-arr k) ib (aget ib-arr k)
            ra (vec3-at ra-arr k) rb (vec3-at rb-arr k)
            p (v/add (v/scale (vec3-at normal k) (aget pn k))
                     (v/add (v/scale (vec3-at t1 k) (aget p1 k))
                            (v/scale (vec3-at t2 k) (aget p2 k))))]
        (apply-impulse! vel omega ii inv-mass ia ra p -1.0)
        (apply-impulse! vel omega ii inv-mass ib rb p 1.0)))))

(defn- record-approach! [{:keys [vel omega]} cs]
  (let [^doubles normal (:normal cs)
        ^doubles ra-arr (:ra cs) ^doubles rb-arr (:rb cs)
        ^doubles approach (:approach cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)]
    (dotimes [k (long (:n cs))]
      (aset approach k
            (v/dot (relative-velocity vel omega (aget ia-arr k) (aget ib-arr k)
                                      (vec3-at ra-arr k) (vec3-at rb-arr k))
                   (vec3-at normal k))))))

(defn- refresh-anchors!
  "Recomputes lever arms and overlap from where the bodies now are.

  This is TGS. Sequential impulse never calls it, and that single
  difference is the whole of what TGS buys: its later iterations are
  solving the problem as it stands rather than as it was linearized at
  the top of the step.

  It runs twice a substep, over every contact, which on a sixteen course
  wall is eighteen thousand times a step -- so it reads the poses out of
  `body-arrays` rather than the body maps. Through the maps it was a
  `nth` into a persistent vector, four keyword lookups and three
  allocated vectors per contact, and a quarter of the whole TGS step."
  [arrays cs]
  (let [^doubles ra-arr (:ra cs) ^doubles rb-arr (:rb cs)
        ^doubles ral (:ra-local cs) ^doubles rbl (:rb-local cs)
        ^doubles depth (:depth cs) ^doubles normal (:normal cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)
        ^doubles depth0 (:depth0 cs)
        ^doubles pos (:pos arrays)
        pa (a/f64 3) pb (a/f64 3) t0 (a/f64 3)]
    (dotimes [k (long (:n cs))]
      (let [ia (aget ia-arr k) ib (aget ib-arr k)
            _ (anchor! arrays ia (load3! ral k t0) pa)
            _ (anchor! arrays ib (load3! rbl k t0) pb)
            k3 (* k 3) ia3 (* ia 3) ib3 (* ib 3)]
        (aset ra-arr k3 (- (aget pa 0) (aget pos ia3)))
        (aset ra-arr (+ k3 1) (- (aget pa 1) (aget pos (+ ia3 1))))
        (aset ra-arr (+ k3 2) (- (aget pa 2) (aget pos (+ ia3 2))))
        (aset rb-arr k3 (- (aget pb 0) (aget pos ib3)))
        (aset rb-arr (+ k3 1) (- (aget pb 1) (aget pos (+ ib3 1))))
        (aset rb-arr (+ k3 2) (- (aget pb 2) (aget pos (+ ib3 2))))
        ;; How far the two anchors have drifted along the normal since
        ;; the contact was found, added to the overlap it had then.
        (aset depth k (+ (aget depth0 k)
                         (* (- (aget pa 0) (aget pb 0)) (aget normal k3))
                         (* (- (aget pa 1) (aget pb 1)) (aget normal (+ k3 1)))
                         (* (- (aget pa 2) (aget pb 2)) (aget normal (+ k3 2)))))))))

;; ---------------------------------------------------------------------------
;; Moving the bodies

(defn- accelerate!
  "Gravity, over `dt`, on everything that is awake and can move."
  [{:keys [vel]} bodies gravity dt]
  (let [^doubles vel vel
        [gx gy gz] gravity
        dt (double dt)]
    (dotimes [i (count bodies)]
      (when (awake? (nth bodies i))
        (let [b3 (* i 3)]
          (aset vel b3 (+ (aget vel b3) (* (double gx) dt)))
          (aset vel (+ b3 1) (+ (aget vel (+ b3 1)) (* (double gy) dt)))
          (aset vel (+ b3 2) (+ (aget vel (+ b3 2)) (* (double gz) dt))))))))

(defn- write-back
  "Puts the solved velocities back on the bodies."
  [bodies {:keys [vel omega]}]
  (mapv (fn [i b]
          (if-not (awake? b)
            b
            (assoc b :vel (vec3-at ^doubles vel i) :omega (vec3-at ^doubles omega i))))
        (range (count bodies))
        bodies))

(defn- advance
  "Moves and turns every body by its current velocity."
  [bodies dt]
  (let [dt (double dt)]
    (mapv (fn [{:keys [pos rot vel omega damping angular-damping] :as b}]
            (if-not (awake? b)
              b
              (let [lin (max 0.0 (- 1.0 (* (double (or damping 0.0)) dt)))
                    ang (max 0.0 (- 1.0 (* (double (or angular-damping 0.0)) dt)))
                    vel (v/scale vel lin)
                    omega (v/scale omega ang)
                    [dx dy dz dw] (q/mul [(omega 0) (omega 1) (omega 2) 0.0] rot)
                    [rx ry rz rw] rot
                    h (* 0.5 dt)]
                (-> b
                    (assoc :vel vel :omega omega
                           :prev-pos pos :prev-rot rot
                           :pos (v/add-scaled pos vel dt)
                           :rot (q/normalize [(+ rx (* h dx)) (+ ry (* h dy))
                                              (+ rz (* h dz)) (+ rw (* h dw))]))
                    (as-> b' (assoc b' :inv-rot (q/inverse (:rot b'))))))))
          bodies)))

;; ---------------------------------------------------------------------------
;; The three

(defn- contact-state
  "The contacts with their accumulated impulses, for next step to guess
  from.

  Warm starting is most of why a stack stands at eight iterations rather
  than fifty: the force a contact needed last step is a very good
  estimate of what it needs this step, and starting from it means the
  iterations refine an answer instead of discovering one."
  [cs]
  (let [^doubles pn (:pn cs) ^doubles p1 (:p1 cs) ^doubles p2 (:p2 cs)]
    (mapv (fn [k c] (assoc c :pn (aget pn (long k)) :p1 (aget p1 (long k)) :p2 (aget p2 (long k))))
          (range (:n cs))
          (:contacts cs))))

(defn- project-contacts!
  "One XPBD pass: push overlapping bodies apart, positions only.

  This is `allgo.physics.rigid/correct` -- share a correction between
  two bodies by how much each gives at the point it acts -- written
  against `body-arrays` rather than the body vector. It is the same
  arithmetic in the same order, and the reason for the copy is that
  `correct` returns a new `:bodies`: a persistent vector of 257 maps,
  rebuilt twice for every contact that pushed. At 2,700 contacts and
  eight passes that was three quarters of XPBD's step spent allocating
  bodies it was about to throw away. Here the pass moves the poses in
  place and the body vector is rebuilt once, at the end of the substep.

  The angular half of each correction is relaxed by half, which
  `allgo.physics.joint` also asks for and for the same reason. A brick in
  a wall carries eight or ten contact points at once and each one turns
  it the whole way on its own; they fight, the substep ends somewhere
  none of them asked for, and that leftover displacement comes back as
  velocity divided by `h`. Relaxed, a twelve course wall settles dead;
  unrelaxed it was doing twenty meters a second within a second, and
  solving it harder -- more passes, more substeps -- made it worse, which
  is what says the trouble is the read-back and not convergence."
  [arrays cs slop max-push]
  (let [^doubles normal (:normal cs)
        ^doubles ral (:ra-local cs) ^doubles rbl (:rb-local cs)
        ^doubles depth0 (:depth0 cs)
        ^doubles lambda (:lambda cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)
        ^doubles pos (:pos arrays)
        ^doubles rot (:rot arrays) ^doubles inv-rot (:inv-rot arrays)
        ^doubles inv-mass (:inv-mass arrays)
        ^ints movable (:movable arrays)
        slop (double slop) max-push (double max-push)
        ;; Scratch, allocated once for the whole pass rather than per
        ;; contact. `pa` and `pb` have to survive the correction that
        ;; follows them; the rest are working room.
        pa (a/f64 3) pb (a/f64 3) t0 (a/f64 3) t1 (a/f64 3) t2 (a/f64 3)]
    (dotimes [k (long (:n cs))]
      (let [ia (aget ia-arr k) ib (aget ib-arr k)
            ma (aget movable ia) mb (aget movable ib)]
        (when-not (and (zero? ma) (zero? mb))
          (let [_ (anchor! arrays ia (load3! ral k t0) pa)
                _ (anchor! arrays ib (load3! rbl k t0) pb)
                k3 (* k 3)
                nx (aget normal k3) ny (aget normal (+ k3 1)) nz (aget normal (+ k3 2))
                pen (+ (aget depth0 k)
                       (* (- (aget pa 0) (aget pb 0)) nx)
                       (* (- (aget pa 1) (aget pb 1)) ny)
                       (* (- (aget pa 2) (aget pb 2)) nz))
                ;; Bounded, for the same reason and by the same number
                ;; the impulse solvers use: however deep the overlap, it
                ;; is not pushed out faster than `max-push-speed`. A deep
                ;; overlap takes several substeps to clear instead of
                ;; one, which is slower and is not 60 meters a second.
                c (min (- pen slop) max-push)]
            (when (pos? c)
              ;; The correction is `-c` along the normal, so the unit
              ;; direction the pair's compliance is measured along is the
              ;; normal reversed. Which way it points makes no difference
              ;; to the effective mass -- that is quadratic in it -- but
              ;; it does to the impulse.
              (let [dx (- nx) dy (- ny) dz (- nz)
                    ia3 (* ia 3) ib3 (* ib 3)
                    ;; (at - pos) x d, the lever arm crossed with the
                    ;; direction, for each body.
                    ax (- (aget pa 0) (aget pos ia3))
                    ay (- (aget pa 1) (aget pos (+ ia3 1)))
                    az (- (aget pa 2) (aget pos (+ ia3 2)))
                    bx (- (aget pb 0) (aget pos ib3))
                    by (- (aget pb 1) (aget pos (+ ib3 1)))
                    bz (- (aget pb 2) (aget pos (+ ib3 2)))
                    wa (if (zero? ma)
                         0.0
                         (do (aset t1 0 (- (* ay dz) (* az dy)))
                             (aset t1 1 (- (* az dx) (* ax dz)))
                             (aset t1 2 (- (* ax dy) (* ay dx)))
                             (+ (aget inv-mass ia) (angular-mass arrays ia t1 t2))))
                    wb (if (zero? mb)
                         0.0
                         (do (aset t1 0 (- (* by dz) (* bz dy)))
                             (aset t1 1 (- (* bz dx) (* bx dz)))
                             (aset t1 2 (- (* bx dy) (* by dx)))
                             (+ (aget inv-mass ib) (angular-mass arrays ib t1 t2))))
                    w (+ wa wb)]
                (when (pos? w)
                  ;; `lambda` is `-c/w`, and the impulse is the direction
                  ;; scaled by its negation. What the velocity pass wants
                  ;; recorded is the magnitude, which is `c/w` -- the
                  ;; force over dt squared, times dt squared.
                  (let [mag (/ c w)
                        px (* dx mag) py (* dy mag) pz (* dz mag)]
                    (aset lambda k (+ (aget lambda k) mag))
                    (when-not (zero? ma)
                      (let [im (aget inv-mass ia)]
                        (aset pos ia3 (+ (aget pos ia3) (* px im)))
                        (aset pos (+ ia3 1) (+ (aget pos (+ ia3 1)) (* py im)))
                        (aset pos (+ ia3 2) (+ (aget pos (+ ia3 2)) (* pz im))))
                      ;; The lever arm is measured from the moved center,
                      ;; as the reference measures it. It makes no
                      ;; difference: the center moved along `p`, and a
                      ;; vector crossed with something parallel to itself
                      ;; is zero.
                      (let [rx (- (aget pa 0) (aget pos ia3))
                            ry (- (aget pa 1) (aget pos (+ ia3 1)))
                            rz (- (aget pa 2) (aget pos (+ ia3 2)))]
                        (aset t1 0 (- (* ry pz) (* rz py)))
                        (aset t1 1 (- (* rz px) (* rx pz)))
                        (aset t1 2 (- (* rx py) (* ry px)))
                        (qrot! inv-rot ia t1 t2)
                        (inertia-mul! (:inv-inertia arrays) ia t2 t1)
                        (qrot! rot ia t1 t2)
                        ;; Halved, which `allgo.physics.joint` also asks
                        ;; for and for the same reason -- see above.
                        (aset t2 0 (* 0.5 (aget t2 0)))
                        (aset t2 1 (* 0.5 (aget t2 1)))
                        (aset t2 2 (* 0.5 (aget t2 2)))
                        (spin! rot inv-rot ia t2)))
                    (when-not (zero? mb)
                      (let [im (aget inv-mass ib)]
                        (aset pos ib3 (- (aget pos ib3) (* px im)))
                        (aset pos (+ ib3 1) (- (aget pos (+ ib3 1)) (* py im)))
                        (aset pos (+ ib3 2) (- (aget pos (+ ib3 2)) (* pz im))))
                      (let [rx (- (aget pb 0) (aget pos ib3))
                            ry (- (aget pb 1) (aget pos (+ ib3 1)))
                            rz (- (aget pb 2) (aget pos (+ ib3 2)))
                            qx (- px) qy (- py) qz (- pz)]
                        (aset t1 0 (- (* ry qz) (* rz qy)))
                        (aset t1 1 (- (* rz qx) (* rx qz)))
                        (aset t1 2 (- (* rx qy) (* ry qx)))
                        (qrot! inv-rot ib t1 t2)
                        (inertia-mul! (:inv-inertia arrays) ib t2 t1)
                        (qrot! rot ib t1 t2)
                        (aset t2 0 (* 0.5 (aget t2 0)))
                        (aset t2 1 (* 0.5 (aget t2 1)))
                        (aset t2 2 (* 0.5 (aget t2 2)))
                        (spin! rot inv-rot ib t2)))))))))))
    arrays))

(defn- roots
  "Union-find over the contacts: which dynamic bodies move together.

  Static bodies join nothing. The floor touches every brick in the scene,
  and an island that included it would be the whole scene -- one group
  that can only sleep when the last thing in it has stopped, which on a
  wall with a ball still rolling somewhere is never."
  [bodies contacts]
  (let [n (count bodies)
        ^ints parent (a/i32 (range n))
        find (fn find [^long i]
               (let [p (aget parent i)]
                 (if (= p i)
                   i
                   (let [r (find p)]
                     (aset parent i (int r))
                     r))))]
    (doseq [c contacts]
      (let [a (long (:a c)) b (long (:b c))]
        (when (and (dynamic? (nth bodies a)) (dynamic? (nth bodies b)))
          (let [ra (find a) rb (find b)]
            (when-not (= ra rb) (aset parent (int ra) (int rb)))))))
    (mapv (fn [i] (if (dynamic? (nth bodies i)) (find i) -1)) (range n))))

(defn- settle
  "Sleeps whole islands, and wakes any island something has disturbed.

  Two halves, and both are necessary. A body under the speed thresholds
  accumulates time on its own clock; one over them puts its clock back to
  zero. But a body cannot sleep on its own account -- it sleeps when
  every body it is touching has also been still long enough, which is
  what the union-find is for. A brick resting on a brick that is still
  rolling has to stay awake, and it looks perfectly still while doing it.

  This is what makes a settled scene actually settle. A stack that has
  stopped is not nearly stopped: it is carrying the solver's own residual,
  a millimeter a second of it, and an unstable arrangement amplifies that
  until it falls over. Sleeping is the only thing that takes the residual
  to zero, and it is why a wall in any shipping engine stands overnight."
  [bodies contacts dt opts]
  (if-not (:allow-sleep? opts)
    bodies
    (let [dt (double dt)
          lin (double (:linear-sleep-speed opts))
          ang (double (:angular-sleep-speed opts))
          wait (double (:time-to-sleep opts))
          still? (fn [b] (and (< (v/length (:vel b)) lin)
                              (< (v/length (:omega b)) ang)))
          ;; Each body's own clock first.
          ticked (mapv (fn [b]
                         (if (rigid/static? b)
                           b
                           (assoc b :sleep-time
                                  (if (still? b)
                                    (+ (double (or (:sleep-time b) 0.0)) dt)
                                    0.0))))
                       bodies)
          ;; An island sleeps on the clock of whichever member has been
          ;; still for the least time, so if nothing has reached the
          ;; threshold then no island can, and the whole of the rest of
          ;; this -- the union-find and the pass over it -- is work with a
          ;; known answer. A scene in motion never gets past here.
          group (when (some #(>= (double (or (:sleep-time %) 0.0)) wait) ticked)
                  (roots ticked contacts))
          ;; An island sleeps on the clock of whichever member has been
          ;; still for the least time.
          youngest (when group
                     (reduce (fn [m i]
                               (let [r (nth group i)]
                                 (if (neg? (long r))
                                   m
                                   (assoc m r (min (double (get m r ##Inf))
                                                   (double (or (:sleep-time (nth ticked i)) 0.0)))))))
                             {}
                             (range (count ticked))))]
      (if-not group
        ticked
        (mapv (fn [i b]
                (let [r (nth group i)]
                  (if (neg? (long r))
                    b
                    (let [asleep? (>= (double (get youngest r 0.0)) wait)]
                      (if asleep?
                        (assoc b :sleeping? true :vel v/zero :omega v/zero)
                        (assoc b :sleeping? false))))))
              (range (count ticked))
              ticked)))))

(defn- rouse
  "Wakes anything asleep that an awake body has come to touch.

  Run before the solve, on this step's contacts, so a sleeping stack is
  awake and solved in the same step the ball reaches it rather than the
  one after -- a step late is a ball halfway through a wall."
  [bodies contacts]
  (let [disturbed (into #{}
                        (mapcat (fn [c]
                                  (let [a (nth bodies (:a c)) b (nth bodies (:b c))]
                                    (cond
                                      (and (rigid/sleeping? a) (awake? b)) [(:a c)]
                                      (and (rigid/sleeping? b) (awake? a)) [(:b c)]
                                      :else []))))
                        contacts)]
    (if (empty? disturbed)
      bodies
      (reduce (fn [bs i] (update bs i rigid/wake)) bodies disturbed))))

(defn- step-sequential-impulse
  "Linearize once, solve velocities, then move.

  Every iteration works against the lever arms and overlaps measured at
  the top of the step. Cheapest per iteration, and the one that sags when
  a stack is tall or a body is turning quickly."
  [{:keys [bodies gravity iterations warm-start?] :as w} dt]
  (let [contacts (contact/all bodies (:broad w) dt)
        bodies (rouse bodies contacts)
        arrays (body-arrays bodies)
        _ (accelerate! arrays bodies gravity dt)
        cs (prepare arrays contacts (:contacts w))]
    (record-approach! arrays cs)
    (when warm-start? (warm-start! arrays cs))
    (dotimes [_ (long iterations)]
      (solve-velocities! arrays cs dt w))
    (assoc w
           :bodies (-> bodies (write-back arrays) (advance dt)
                       (settle contacts dt w))
           :contacts (contact-state cs))))

(defn- step-tgs
  "Small steps: several short steps, each solved and each moved.

  The shape is the one Box2D v3 and PhysX's TGS Soft use, and Macklin's
  *Small Steps in Physics Simulation* is the argument for it -- several
  small steps solved once beats one step solved many times, because every
  substep re-integrates and the constraint is solved against where the
  bodies now are rather than where they were.

  Per substep, in this order and for these reasons:

    integrate velocities   gravity, over the substep
    refresh anchors        lever arms and separation from where the
                           bodies are now. Collision detection itself
                           runs once for the whole step -- it costs more
                           than the entire solve -- but the anchors are
                           kept in body space, so this recovers the
                           geometry without it. That is the whole trick.
    warm start             apply the impulses this contact needed last
                           time. Every substep, not just the first: the
                           accumulated impulse is a force over a substep,
                           and a force that held the stack up a moment
                           ago is the best guess there is for now.
    solve, with bias       the soft contact, pushing overlap out and
                           holding gaps open
    integrate positions    move by the velocity just solved
    relax, without bias    solve again with no push and no bounce, which
                           takes back the velocity the push added. This
                           is what stops a settled stack from breathing,
                           and it is why no pseudo velocity is needed.

  Restitution is measured once, against the closing speed as the step
  began -- after a substep of solving it is gone. Applying it is a pass
  of its own, after the last substep, and that is not tidiness. The
  relax pass exists to remove separating velocity that the bias put in,
  and it cannot tell that from a rebound: a bouncy ball solved with
  restitution inside the substep left the floor and had the leaving
  speed taken straight back off it by the relax that followed. TGS did
  not bounce at all, at any coefficient. So the substeps run dead and
  the rebound is applied once at the end, to the contacts that actually
  pushed, which is what Box2D v3 does and for this reason."
  [{:keys [bodies gravity iterations substeps warm-start?] :as w} dt]
  (let [substeps (max 1 (long substeps))
        h (/ (double dt) substeps)
        per (max 1 (quot (long iterations) substeps))
        ;; No bounce inside the substeps -- see above.
        solve-opts (assoc w :restitution 0.0)
        relax-opts (assoc w :relax? true)
        bounce-opts (assoc w :bounce? true)
        contacts (contact/all bodies (:broad w) dt)
        bodies (rouse bodies contacts)]
    (loop [bodies bodies
           cs (prepare (body-arrays bodies) contacts (:contacts w))
           n substeps
           first? true]
      (if (zero? n)
        (let [bounced (let [arrays (body-arrays bodies)]
                        (refresh-anchors! arrays cs)
                        (solve-velocities! arrays cs h bounce-opts)
                        (write-back bodies arrays))]
          (assoc w
                 :bodies (settle bounced contacts dt w)
                 :contacts (contact-state cs)))
        (let [arrays (body-arrays bodies)
              _ (accelerate! arrays bodies gravity h)
              _ (refresh-anchors! arrays cs)
              _ (when first? (record-approach! arrays cs))
              _ (when warm-start? (warm-start! arrays cs))
              _ (dotimes [_ per] (solve-velocities! arrays cs h solve-opts))
              moved (-> bodies (write-back arrays) (advance h))
              ;; The relax pass sees the bodies where the substep left
              ;; them, so the anchors are measured again first.
              relaxed (let [arrays (body-arrays moved)]
                        (refresh-anchors! arrays cs)
                        (solve-velocities! arrays cs h relax-opts)
                        (write-back moved arrays))]
          (recur relaxed cs (dec n) false))))))

(def ^:private max-slices
  "The most substeps XPBD will cut a step into.

  Conservative advancement shortens a substep to the moment of impact,
  and a scene contrived enough -- several fast things arriving at
  different times -- could ask for that over and over. Sixty-four
  bounds the work, and the step it leaves is short enough that anything
  still stepping over something at that size was going to need a
  different engine anyway."
  64)

(defn- thinnest
  "A body's narrowest dimension -- what it could be stepped over."
  ^double [b]
  (if (= :ball (:shape b))
    (* 2.0 (double (:radius b)))
    (reduce min (map double (:size b)))))

(defn- outrunning?
  "Whether any awake body would cross more than half its own thinnest
  dimension in a substep of `h`.

  The question that decides whether a time of impact is worth asking
  for. Nothing in an ordinary scene answers yes, so an ordinary scene
  pays one length and one comparison per body and nothing else."
  [bodies ^double h]
  (boolean (some (fn [b]
                   (and (awake? b)
                        (> (* (v/length (:vel b)) h) (* 0.5 (thinnest b)))))
                 bodies)))

(defn- safe-slice
  "The longest a substep may be without anything stepping over anything,
  or nil if nothing is going fast enough for that to be in question.

  Conservative advancement, over the pairs collision detection already
  offered. A time of impact of zero is a pair that is touching now,
  which is the ordinary case and says nothing about how long the next
  substep may be -- it is the pairs that are *about* to touch that set
  the limit."
  [bodies contacts h remaining slop]
  (let [h (double h) remaining (double remaining) slop (double slop)]
    (when (outrunning? bodies h)
      (let [pairs (into #{} (map (juxt :a :b)) contacts)]
        (reduce (fn [best [i j]]
                  (let [a (nth bodies i) b (nth bodies j)]
                    (if (and (rigid/inert? a) (rigid/inert? b))
                      best
                      (let [{:keys [t closing]} (toi/impact a b remaining)]
                        (if (and t (> (double t) 1e-9))
                        ;; Past the moment of impact, not onto it. A
                        ;; position solver corrects `overlap - slop` and
                        ;; has nothing to do at an overlap of zero, so
                        ;; arriving exactly on contact leaves the body
                        ;; traveling at the speed it arrived with and
                        ;; the next slice carries it through. Twice the
                        ;; slop is the shallowest overlap that is
                        ;; certainly worth solving.
                          (let [t' (if (> (double closing) 1e-9)
                                     (+ (double t) (/ (* 2.0 slop) (double closing)))
                                     (double t))]
                            (if (or (nil? best) (< t' (double best))) t' best))
                          best)))))
                nil
                pairs)))))

(defn- step-xpbd
  "Do not solve velocities: move the bodies until they no longer overlap,
  then read the velocity back off how far each traveled.

  Nothing here can inject energy, which is the appeal. The cost is a
  velocity pass afterward for restitution and friction, which positions
  alone cannot express."
  [{:keys [bodies gravity substeps iterations] :as w} dt]
  (let [substeps (max 1 (long substeps))
        h (/ (double dt) substeps)
        contacts (contact/all bodies (:broad w) dt)
        bodies (rouse bodies contacts)
        cs (prepare (body-arrays bodies) contacts (:contacts w))
        slop (double (:slop w))
        ;; The furthest a contact may push in one substep is bounded --
        ;; the impulse solvers have had this since they were written and
        ;; this one did not, which is the whole of why it explodes: it
        ;; reads velocity back off the correction, over the substep, so
        ;; an unbounded correction is an unbounded velocity. It is bound
        ;; per slice now, since the slices are no longer all one length.
        ^doubles lambda (:lambda cs)
        passes (max 1 (quot (long iterations) substeps))]
    (loop [bodies bodies remaining (double dt) slices 0]
      (if (or (<= remaining 1e-9) (>= slices max-slices))
        (assoc w
               :bodies (settle bodies contacts dt w)
               :contacts (contact-state cs))
        (let [;; How long this substep may be. Ordinarily the nominal
              ;; one; shorter when something is about to be stepped
              ;; over, which is the whole of what conservative
              ;; advancement is for here.
              h (max (/ (double dt) max-slices)
                     (min h remaining
                          (or (safe-slice bodies contacts h remaining slop) h)))
              push-limit (* (double (:max-push-speed w)) h)
              ;; Only this substep's multipliers are wanted, so they start
              ;; again each time round.
              _ (dotimes [k (long (:n cs))] (aset lambda k 0.0))
              bodies (mapv #(if (awake? %) (rigid/integrate % h gravity) %) bodies)
              ;; One read of the bodies for the whole position solve.
              ;; The passes move the poses in these arrays and the body
              ;; vector is rebuilt once, below.
              posed (body-arrays bodies)
              _ (refresh-anchors! posed cs)
              ;; How fast the surfaces were closing before the solve.
              ;; Restitution is measured against this; after the position
              ;; solve it is whatever the pushout left behind.
              _ (record-approach! posed cs)
              _ (dotimes [_ (long passes)] (project-contacts! posed cs slop push-limit))
              bodies (write-poses bodies posed)
              bodies (mapv #(if (awake? %) (rigid/update-velocities % h) %) bodies)
              ;; The velocity pass belongs *inside* the substep, not once
              ;; at the end of the frame. It is the only thing that takes
              ;; energy back out -- restitution and friction both -- and
              ;; running it once per frame left three substeps' worth of
              ;; sideways motion to accumulate unopposed.
              arrays (body-arrays bodies)
              _ (refresh-anchors! arrays cs)
              _ (solve-xpbd-velocities! arrays cs h w)]
          (recur (write-back bodies arrays) (- remaining h) (inc slices)))))))

(defn step
  "One step of the world, by whichever `:solver` it carries."
  ([w] (step w (/ 1.0 60.0)))
  ([w dt]
   (case (:solver w)
     :tgs (step-tgs w dt)
     :xpbd (step-xpbd w dt)
     (step-sequential-impulse w dt))))
