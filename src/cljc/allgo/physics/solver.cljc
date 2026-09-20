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
                          velocity back off how far each one travelled.
                          `allgo.physics.rigid` and Müller's papers.

  ## What actually separates them

  All three are Gauss-Seidel: they visit constraints one at a time and
  each one sees the effect of the last. The differences are *when the
  geometry is evaluated* and *what quantity is being corrected*.

  Sequential impulse linearises once, at the top of the step, and every
  iteration after that works against stale lever arms. That is fine when
  bodies barely move in a step and wrong when they do -- a fast spin, or
  a deep stack settling a long way. It is also the cheapest per
  iteration, because nothing has to be recomputed.

  TGS re-evaluates. Anchors are kept in each body's own frame, so after
  the bodies move the contact's world position and separation can be
  recovered without running collision detection again. That is the whole
  trick, and it is why TGS holds a tall stack that sequential impulse
  lets sag: the twentieth iteration is solving the real problem rather
  than a linearisation of the problem as it was twenty iterations ago.

  XPBD corrects positions instead of velocities, which cannot inject
  energy -- the worst a bad correction can do is put a body somewhere and
  have the next one move it back. It pays for that by needing the
  velocity pass afterwards to get restitution and friction, and by being
  stiffer to tune.

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
   :warm-start? true})

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
  "A simulation over `bodies`, which are `allgo.physics.rigid` maps."
  ([bodies] (world bodies {}))
  ([bodies opts] (merge defaults opts {:bodies (vec bodies) :contacts []})))

(defn- contact-key
  "What counts as the same contact as last step, for warm starting.

  The pair, and the feature id `allgo.physics.contact` stamps on the
  point -- which corner of which face is pressed into which face. It used
  to be the contact point rounded to two centimetres, and that is a fine
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
  sequential impulse does not and ignores them."
  [bodies contacts previous]
  (let [n (count contacts)
        ;; Last step's impulses for the same contact, which is the guess
        ;; warm starting rests on.
        old (reduce (fn [m c] (assoc m (contact-key c) [(:pn c) (:p1 c) (:p2 c)]))
                    {}
                    previous)
        remembered (mapv (fn [c] (or (get old (contact-key c)) [0.0 0.0 0.0])) contacts)]
    {:n n
     :a (a/i32 (map :a contacts))
     :b (a/i32 (map :b contacts))
     :normal (a/f64 (mapcat :normal contacts))
     :t1 (a/f64 (mapcat #(first (tangents (:normal %))) contacts))
     :t2 (a/f64 (mapcat #(second (tangents (:normal %))) contacts))
     :depth (a/f64 (map :depth contacts))
     ;; The overlap as collision detection found it, kept so that
     ;; `refresh-anchors!` can add the drift to it rather than to a value
     ;; the last substep already moved.
     :depth0 (a/f64 (map :depth contacts))
     :ra (a/f64 (mapcat (fn [c] (v/sub (:point c) (:pos (nth bodies (:a c))))) contacts))
     :rb (a/f64 (mapcat (fn [c] (v/sub (:point c) (:pos (nth bodies (:b c))))) contacts))
     :ra-local (a/f64 (mapcat (fn [c] (rigid/world->local (nth bodies (:a c)) (:point c)))
                              contacts))
     :rb-local (a/f64 (mapcat (fn [c] (rigid/world->local (nth bodies (:b c)) (:point c)))
                              contacts))
     :pn (a/f64 (map #(nth % 0) remembered))
     ;; The tangential impulses come back too, now that the match is
     ;; exact. They were being dropped -- stored by `contact-state` every
     ;; step and never read -- so friction began each step from nothing
     ;; and had to be rediscovered in the iterations it had left, which
     ;; is what let a stack creep sideways while it stood.
     :p1 (a/f64 (map #(nth % 1) remembered))
     :p2 (a/f64 (map #(nth % 2) remembered))
     ;; What the XPBD position solve had to push with to keep this contact
     ;; apart, summed over the last substep. Only that solver fills it,
     ;; and it is what bounds friction there.
     :lambda (a/f64 n)
     ;; The closing speed as the step began, which is what restitution is
     ;; measured against -- after an iteration or two it is gone.
     :approach (a/f64 n)
     :contacts (vec contacts)}))

(defn- body-arrays [bodies]
  (let [n (count bodies)
        vel (a/f64 (* n 3))
        omega (a/f64 (* n 3))
        inv-mass (a/f64 n)
        ii (a/f64 (* n 9))]
    (dotimes [i n]
      (let [b (nth bodies i)
            [vx vy vz] (:vel b)
            [wx wy wz] (:omega b)]
        (aset vel (* i 3) (double vx))
        (aset vel (+ (* i 3) 1) (double vy))
        (aset vel (+ (* i 3) 2) (double vz))
        (aset omega (* i 3) (double wx))
        (aset omega (+ (* i 3) 1) (double wy))
        (aset omega (+ (* i 3) 2) (double wz))
        (aset inv-mass i (double (:inv-mass b)))
        (inverse-inertia-world ii i (:rot b) (:inv-inertia b))))
    {:vel vel :omega omega :inv-mass inv-mass :ii ii}))

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
            pen (- (aget depth k) slop)
            b-term (if (pos? pen) (min (* bias-rate pen) max-push) 0.0)
            va (aget approach k)
            r-term (if (< va (- threshold)) (* (- restitution) va) 0.0)
            normal-target (+ b-term r-term)]
        ;; Three directions in turn -- the normal, then two tangents --
        ;; rather than a closure called three times. The closure was
        ;; allocated per contact per sweep, two and a half thousand of
        ;; them a step, and in a browser that cost more than the
        ;; arithmetic it wrapped.
        (loop [dir 0 limit 0.0]
          (when (< dir 3)
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
                           (- (* m mass-scale (- target vd)) (* impulse-scale old))
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
  position solve travelling fast. A dropped box came off the floor at
  fifteen metres a second before this existed.

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
      ;; ball was leaving at four metres a second and this pass set it to
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
          ;; millimetre of pushout is a quarter of a metre a second -- and
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
  solving the problem as it stands rather than as it was linearised at
  the top of the step."
  [bodies cs]
  (let [^doubles ra-arr (:ra cs) ^doubles rb-arr (:rb cs)
        ^doubles ral (:ra-local cs) ^doubles rbl (:rb-local cs)
        ^doubles depth (:depth cs) ^doubles normal (:normal cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)
        ^doubles depth0 (:depth0 cs)]
    (dotimes [k (long (:n cs))]
      (let [a (nth bodies (aget ia-arr k))
            b (nth bodies (aget ib-arr k))
            pa (rigid/local->world a (vec3-at ral k))
            pb (rigid/local->world b (vec3-at rbl k))
            n (vec3-at normal k)
            ra (v/sub pa (:pos a))
            rb (v/sub pb (:pos b))]
        (dotimes [c 3]
          (aset ra-arr (+ (* k 3) c) (double (nth ra c)))
          (aset rb-arr (+ (* k 3) c) (double (nth rb c))))
        ;; How far the two anchors have drifted along the normal since
        ;; the contact was found, added to the overlap it had then.
        (aset depth k (+ (aget depth0 k) (v/dot (v/sub pa pb) n)))))))

;; ---------------------------------------------------------------------------
;; Moving the bodies

(defn- accelerate! [{:keys [vel]} bodies gravity dt]
  (let [^doubles vel vel
        [gx gy gz] gravity
        dt (double dt)]
    (dotimes [i (count bodies)]
      (when-not (rigid/static? (nth bodies i))
        (let [b3 (* i 3)]
          (aset vel b3 (+ (aget vel b3) (* (double gx) dt)))
          (aset vel (+ b3 1) (+ (aget vel (+ b3 1)) (* (double gy) dt)))
          (aset vel (+ b3 2) (+ (aget vel (+ b3 2)) (* (double gz) dt))))))))

(defn- write-back
  "Puts the solved velocities back on the bodies."
  [bodies {:keys [vel omega]}]
  (mapv (fn [i b]
          (if (rigid/static? b)
            b
            (assoc b :vel (vec3-at ^doubles vel i) :omega (vec3-at ^doubles omega i))))
        (range (count bodies))
        bodies))

(defn- advance
  "Moves and turns every body by its current velocity."
  [bodies dt]
  (let [dt (double dt)]
    (mapv (fn [{:keys [pos rot vel omega damping angular-damping] :as b}]
            (if (rigid/static? b)
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

(defn- project-contacts
  "One XPBD pass: push overlapping bodies apart, positions only.

  `allgo.physics.rigid/correct` already knows how to share a correction
  between two bodies by how much each gives at the point it acts -- which
  is the same question a contact asks, so there is nothing to add here
  but the error vector, and the multiplier it needed on the way back.

  The angular half of each correction is relaxed by half, which
  `allgo.physics.joint` also asks for and for the same reason. A brick in
  a wall carries eight or ten contact points at once and each one turns
  it the whole way on its own; they fight, the substep ends somewhere
  none of them asked for, and that leftover displacement comes back as
  velocity divided by `h`. Relaxed, a twelve course wall settles dead;
  unrelaxed it was doing twenty metres a second within a second, and
  solving it harder -- more passes, more substeps -- made it worse, which
  is what says the trouble is the read-back and not convergence."
  [bodies cs slop dt]
  (let [^doubles normal (:normal cs)
        ^doubles ral (:ra-local cs) ^doubles rbl (:rb-local cs)
        ^doubles depth0 (:depth0 cs)
        ^doubles lambda (:lambda cs)
        ^ints ia-arr (:a cs) ^ints ib-arr (:b cs)
        dt (double dt)]
    (reduce (fn [bs k]
              (let [ia (aget ia-arr k) ib (aget ib-arr k)
                    a (nth bs ia) b (nth bs ib)
                    pa (rigid/local->world a (vec3-at ral k))
                    pb (rigid/local->world b (vec3-at rbl k))
                    n (vec3-at normal k)
                    pen (+ (aget depth0 k) (v/dot (v/sub pa pb) n))
                    c (- pen (double slop))]
                (if (pos? c)
                  (let [{:keys [bodies force]}
                        (rigid/correct bs {:a ia :b ib
                                           :corr (v/scale n (- c))
                                           :at pa :other-at pb :dt dt
                                           :angular-relaxation 0.5})]
                    ;; `correct` reports the force; the multiplier behind
                    ;; it is that over dt squared, and summing it is how
                    ;; the velocity pass learns how hard this contact was
                    ;; pushing.
                    (aset lambda k (+ (aget lambda k)
                                      (abs (* (double force) dt dt))))
                    bodies)
                  bs)))
            bodies
            (range (long (:n cs))))))

(defn- step-sequential-impulse
  "Linearise once, solve velocities, then move.

  Every iteration works against the lever arms and overlaps measured at
  the top of the step. Cheapest per iteration, and the one that sags when
  a stack is tall or a body is turning quickly."
  [{:keys [bodies gravity iterations warm-start?] :as w} dt]
  (let [arrays (body-arrays bodies)
        _ (accelerate! arrays bodies gravity dt)
        cs (prepare bodies (contact/all bodies) (:contacts w))]
    (record-approach! arrays cs)
    (when warm-start? (warm-start! arrays cs))
    (dotimes [_ (long iterations)]
      (solve-velocities! arrays cs dt w))
    (assoc w
           :bodies (-> bodies (write-back arrays) (advance dt))
           :contacts (contact-state cs))))

(defn- step-tgs
  "Solve, move a little, re-measure, solve again.

  Collision detection still runs once per step -- it is far too expensive
  to repeat -- but the anchors are kept in body space, so after each
  substep the contact's lever arms and overlap can be recovered from
  where the bodies now are. That is what the later iterations are solving
  against, and it is the whole difference from sequential impulse."
  [{:keys [bodies gravity iterations substeps warm-start?] :as w} dt]
  (let [substeps (max 1 (long substeps))
        h (/ (double dt) substeps)
        per (max 1 (quot (long iterations) substeps))
        contacts (contact/all bodies)]
    (loop [bodies bodies cs (prepare bodies contacts (:contacts w)) n substeps first? true]
      (if (zero? n)
        (assoc w :bodies bodies :contacts (contact-state cs))
        (let [arrays (body-arrays bodies)]
          (accelerate! arrays bodies gravity h)
          (refresh-anchors! bodies cs)
          (when first? (record-approach! arrays cs))
          (when (and warm-start? first?) (warm-start! arrays cs))
          (dotimes [_ per] (solve-velocities! arrays cs h w))
          (recur (-> bodies (write-back arrays) (advance h)) cs (dec n) false))))))

(defn- step-xpbd
  "Do not solve velocities: move the bodies until they no longer overlap,
  then read the velocity back off how far each travelled.

  Nothing here can inject energy, which is the appeal. The cost is a
  velocity pass afterwards for restitution and friction, which positions
  alone cannot express."
  [{:keys [bodies gravity substeps iterations] :as w} dt]
  (let [substeps (max 1 (long substeps))
        h (/ (double dt) substeps)
        contacts (contact/all bodies)
        cs (prepare bodies contacts (:contacts w))
        slop (double (:slop w))
        ^doubles lambda (:lambda cs)
        passes (max 1 (quot (long iterations) substeps))]
    (loop [bodies bodies n substeps]
      (if (zero? n)
        (assoc w :bodies bodies :contacts (contact-state cs))
        (let [;; Only this substep's multipliers are wanted, so they start
              ;; again each time round.
              _ (dotimes [k (long (:n cs))] (aset lambda k 0.0))
              bodies (mapv #(rigid/integrate % h gravity) bodies)
              _ (refresh-anchors! bodies cs)
              ;; How fast the surfaces were closing before the solve.
              ;; Restitution is measured against this; after the position
              ;; solve it is whatever the pushout left behind.
              _ (record-approach! (body-arrays bodies) cs)
              bodies (reduce (fn [bs _] (project-contacts bs cs slop h))
                             bodies
                             (range passes))
              bodies (mapv #(rigid/update-velocities % h) bodies)
              ;; The velocity pass belongs *inside* the substep, not once
              ;; at the end of the frame. It is the only thing that takes
              ;; energy back out -- restitution and friction both -- and
              ;; running it once per frame left three substeps' worth of
              ;; sideways motion to accumulate unopposed.
              arrays (body-arrays bodies)
              _ (refresh-anchors! bodies cs)
              _ (solve-xpbd-velocities! arrays cs h w)]
          (recur (write-back bodies arrays) (dec n)))))))

(defn step
  "One step of the world, by whichever `:solver` it carries."
  ([w] (step w (/ 1.0 60.0)))
  ([w dt]
   (case (:solver w)
     :tgs (step-tgs w dt)
     :xpbd (step-xpbd w dt)
     (step-sequential-impulse w dt))))
