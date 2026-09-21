(ns allgo.physics.articulated
  "Articulated bodies in reduced coordinates, after Featherstone.

  Everywhere else in this library a joint is a *constraint*: two free
  rigid bodies, each with six degrees of freedom, and a solver that
  pushes them back together every step until they behave as though they
  were pinned. `allgo.physics.joint` does that, and for a pile of loose
  parts it is the right shape -- contacts and joints go through the same
  solver and nothing needs to know which is which.

  It is the wrong shape for an arm. A six-link arm has six degrees of
  freedom and the constraint formulation gives it thirty-six, then spends
  its iterations taking thirty of them away again. They never quite go.
  Two links a metre long hanging off a fixed base, the lower one heavier
  than the upper, worst separation between a joint's two anchors over
  two seconds -- anchors that are supposed to be the same point:

      mass ratio   5 substeps   10      20
               1      0.00475   0.00126   0.00029
              10      0.01633   0.00393   0.00098
             100      0.11241   0.04149   0.01001
            1000      2.75497   0.29095   0.09466

  At a thousand to one and five substeps the joint has come nearly three
  metres apart, on links a metre long. More substeps help and do not
  cure: every column still climbs with the ratio. That is the case a
  ragdoll is made of -- a hand on a forearm on an upper arm, each several
  times the mass of the one beyond it.

  In reduced coordinates there is no such number to report. A joint angle
  cannot be violated because there is nowhere for the violation to live.

  The answer is to give the arm six numbers in the first place, so that
  nothing has to be solved to keep it together and what remains is
  working out how it accelerates. That is what this does, in time linear
  in the number of links.

  What it gives up is being able to come apart or close a loop: a tree of
  joints is the whole of what it can describe.

  ## Hitting itself

  `:self-collide?` adds the model's own parts against each other, which
  is what stops a forearm passing through a thigh. A link and its
  parent are never tested -- they meet at the joint and overlap there by
  construction, so a contact between them would be permanent and would
  push the model apart from the inside -- and `:no-collide` on a link
  names any other pair to leave alone.

  The impulse is between two bodies that can both move, so what it
  meets is a *reduced* mass, smaller than either body's own and not
  larger: both ends of the push give, where a contact with the floor
  has one end that does not. Nothing outside the model is touched by
  it, which is the sharp check -- an internal impulse leaves both
  momenta exactly where they were.

  ## Shapes

  A part gains geometry by being given a `:shape` -- `:box` with a
  `:size` or `:ball` with a `:radius` -- and it sits at the part's
  centre of mass unless `:shape-pose` says otherwise. A link's frame is
  at its *joint*, not in the middle of it, so a shape left at the frame
  origin would stick out of the elbow.

  A part with no shape collides with nothing, which is what a link that
  exists only to carry a degree of freedom wants -- the middle of three
  stacked hinges standing in for a shoulder is not a thing that can be
  hit.

  ## Being hit

  It can be hit, though, which took one more thing. A contact solver
  asks a body two questions -- how much velocity does a push here buy,
  and please take this impulse -- and neither has an obvious answer for
  a body whose motion is described by joint angles.

  Both come out of the articulated body algorithm run with the velocity
  terms taken out. An impulse is a force with no duration, so nothing
  has time to contribute a velocity product and the bias forces are the
  impulse itself; the articulated inertias are the ones already built
  for this configuration. What is left is an inward walk from the body
  that was hit to the root and an outward sweep back, which is linear in
  the links -- `impulse-delta`.

  `inverse-mass-matrix` says the same thing as a whole matrix, and is
  kept because it can be checked against inverse dynamics and so carries
  the evidence over to the fast path. It is not what runs: it takes
  `n + 6` runs of the full algorithm to build, and that was, measured,
  the entire cost of a contact -- 33ms of a 34ms solve on sixteen links.

  That number is the interesting one. A hand on the end of an
  outstretched arm is light; the same hand with the arm folded against
  the chest is most of a torso. The constraint formulation gets this
  right too, eventually, by iterating; here it is one matrix product
  and it is exact.

  ## Knowing it is right

  There is not much to see in a number like `qdd`, so the checks are
  quantities that have to hold whatever the algorithm does. A hinged rod
  matches the closed form to twelve digits. The fast path matches
  `M^-1 (tau - bias)`, built out of the slow one, to about fourteen.
  And a free chain under no gravity conserves both its momenta, which is
  the strongest of them: momentum belongs to the whole system rather
  than to any link, so a sign lost anywhere in the recursion shows up
  there even when every link looks plausible on its own. It is what
  caught gravity being added in the wrong frame.

  The same measure covers being hit. An impulse on a free model must
  change its linear momentum by exactly that impulse and its angular
  momentum by exactly that impulse's moment -- which holds only if the
  Jacobian, the inverse inertia and the momentum sum, three separate
  pieces of arithmetic, all agree. They do, to about 5e-16, for a
  rotated base and an oblique push on a middle link.

  On the ground there is no exact answer to check against, so the checks
  are the ones a person would make by eye if eyes were precise: a box
  comes to rest at half its own height and stays there; nothing ends up
  buried; a slope is slid down when it is slippery and not when it
  grips; and a bounce returns to `0.2 + 0.8 e^2` from the height it was
  dropped, to within two centimetres across the range of `e`.

  ## What it costs

  A chain lying on the floor, JVM, per step:

      links         0      2      4      8     16
      contacts      4      6      8     16     52
      at first   0.71   4.92  12.19  29.85  154.1
      structure  0.67   2.03   4.24   6.87   31.35
      flat       0.45   0.86   1.36   1.78    5.21

  Thirty times, in two halves. The structural half: a contact used to
  cost a whole inverse inertia matrix, which cost `n + 6` runs of the
  articulated body algorithm, and is now three O(n) impulse responses;
  the per-configuration data is built once a step rather than once per
  matrix column; the inward walk visits only the path from the body that
  was hit to the root.

  The other half is arithmetic. Everything in the section above is the
  readable statement of the algebra and is what the tests check, and it
  is *not* what the inner loops run on -- a 6x6 as a vector of vectors
  is seven objects and thirty-six boxed doubles. The per-configuration
  cache holds flat `double` arrays instead, the passes over them
  allocate one result apiece, and the contact sweep keeps one mutable
  generalised velocity that every impulse adds into. The conversion
  happens once per link per step, where it costs nothing.

  A word of warning from having done it: the first attempt was *slower*,
  because three type hints were missing and reflective `aget` costs two
  orders of magnitude. `make reflect` is a build step for this reason
  and it is worth running before believing any measurement.

  Two more of those were found later and are worth naming, because both
  were hidden behind a sentence claiming they did not matter. `links`
  built its transforms and inertias as vectors of vectors and converted
  them -- four milliseconds of a settled ragdoll's step -- and
  `advance-positions` built the whole cache again to read ten integers
  off it. The ragdoll went from 2.16ms a step to 1.27.

  The recursion over it went the same way afterwards. `ls` and the
  articulated inertias are `deftype`s rather than maps, so a link's
  transform is a field read instead of a hash lookup -- six hundred of
  those a step on a settled ragdoll -- and the total degrees of freedom
  is read off the last link rather than walked for, because the inner
  recursion asked for it on every one of its seventy calls.

  That is worth 35% in a browser and is lost in the noise on a JVM,
  which is the expected shape: a keyword lookup costs relatively more
  in JavaScript than beside JVM floating point.

  ## Measuring any of this from a browser

  With care, and the care is not optional. These numbers were nearly
  got wrong twice.

  An automated browser tab is backgrounded, and Chrome gives a
  backgrounded tab less of a processor. Ten million square roots and
  multiplies -- no Clojure in it at all -- take 6.3ms on this JVM and
  110ms in that tab. Sixteen times. So an absolute millisecond figure
  measured there means nothing, and one was published in this
  repository before anybody checked.

  What does survive is a *ratio*. Run the calibration loop beside the
  thing being measured, in the same tab, and divide: the throttle
  cancels. That is how the 35% above was measured, and it is also how
  the ragdoll was shown to be about twice the JVM's cost rather than
  thirty times -- which is the ordinary price of Clojure in JavaScript
  and not a bug worth hunting. A day was spent hunting it.

  ## Spatial vectors

  The algorithms are short because the algebra is. A rigid body's motion
  is six numbers and so is a force on it, and Featherstone's spatial
  vectors treat each as one object rather than an angular part and a
  linear part that have to be kept in step by hand:

      motion  [wx wy wz vx vy vz]   angular velocity, then the velocity
                                    of the body-fixed point at the origin
      force   [nx ny nz fx fy fz]   moment about the origin, then force

  Written this way, Newton and Euler together are `f = I a + v x* I v`,
  one line rather than two coupled ones, and an inertia is a single 6x6.
  The linear part of a motion vector is the velocity of *the point at the
  origin*, not of the centre of mass, which is the one thing that trips
  everybody: a body spinning about a distant origin has a large linear
  part and is not going anywhere.

  Two cross products, not one, because motion and force vectors
  transform differently: `crm` for the rate of change of a motion vector
  carried along by a moving frame, `crf` for a force vector, and
  `crf(v) = -crm(v)^T`.

  ## The model

  A vector of links in topological order -- a link's parent always comes
  before it -- each `{:parent :joint :axis :origin :mass :com :inertia}`.
  `:origin` is the fixed transform from the parent's frame to this
  joint's frame, `:axis` the joint's axis in that frame, `:inertia` the
  3x3 rotational inertia about the link's own centre of mass. Parent -1
  means the root.

  That vector on its own is a chain bolted to the world. Wrapped as
  `{:base {:mass :com :inertia} :links [...]}` the root is free, which
  is what anything that has to fall over needs -- a ragdoll's pelvis
  goes nowhere useful if it is nailed to the origin.

  The base's six degrees of freedom are not coordinates. They are a
  pose and a spatial velocity, carried in the *state* rather than in
  `q`, so that nothing ever has to parameterise a rotation with three
  numbers and find the singularity in it. The base costs one 6x6 solve
  and that is all.

  ## How wide a joint is

  `:revolute` and `:prismatic` are one number, `:spherical` is three,
  `:fixed` is none. So `q` has an entry per *link* -- an angle, a
  distance, or for a ball joint a quaternion -- while `qd` and `tau` are
  flat and as long as the model has degrees of freedom, and a link knows
  where its own rates start in them.

  A ball joint is three numbers of velocity and four of configuration
  for the same reason the base is: the velocity is an angular velocity,
  which is genuinely three numbers, and the configuration is a rotation,
  which no three numbers can name without a direction in which they stop
  working. Three stacked hinges would be the alternative and it is the
  one that finds that direction -- a ragdoll tumbles, so it finds it.
  It is also cheaper this way: a ball joint's motion subspace is the
  identity on the angular half, with no intermediate frames to transform
  through.

  ## Limits

  A joint can be given a range to stay inside, and a ragdoll needs one:
  without limits a body settles with its head folded back on itself and
  its knees bent the wrong way, which is a pile of sticks rather than a
  figure.

  `:limit [lo hi]` bounds a hinge or a slider. `:cone theta` bounds how
  far a ball joint's bone may swing from where it rests, and `:twist`
  how far it may turn about that bone -- two different things, and a
  cone alone leaves a head free to face backwards.

  They are solved as one-sided constraints in the same sweep as the
  contacts, which is why they can be: a limit is a push between a link
  and its own parent where a contact is a push between a link and the
  floor, and the impulse response does not need to be told which. Limits
  go first in each sweep, because a knee that has folded backwards is a
  worse thing to look at than a foot a millimetre into the floor.

  ## Limits

  A joint can be given a range to stay inside, and a ragdoll needs one:
  without limits a body settles with its head folded back on itself and
  its knees bent the wrong way, which is a pile of sticks rather than a
  figure.

  `:limit [lo hi]` bounds a hinge or a slider. `:cone theta` bounds how
  far a ball joint's bone may swing from where it rests, and `:twist`
  how far it may turn about that bone -- two different things, and a
  cone alone leaves a head free to face backwards.

  They are solved as one-sided constraints in the same sweep as the
  contacts, which is why they can be: a limit is a push between a link
  and its own parent where a contact is a push between a link and the
  floor, and the impulse response does not need to be told which. Limits
  go first in each sweep, because a knee that has folded backwards is a
  worse thing to look at than a foot a millimetre into the floor.

  A ball joint does need its body to have inertia about every axis. A
  mathematically thin rod has none about its own length, so the three by
  three it has to invert is singular, and a joint whose accelerations
  cannot be computed does not move -- it welds itself, quietly. Limbs
  should be given the inertia of the shape they are, not of a line.

  ## What runs on it

  `inverse-dynamics` is the recursive Newton-Euler algorithm: given where
  the joints are, how fast they are moving and how fast they are being
  accelerated, what torques does that take. `forward-dynamics` is the
  articulated body algorithm: given the torques, what accelerations. Both
  are O(n) and the second is the one a simulation wants.

  Gravity arrives through neither of them. The whole recursion runs in a
  frame falling at `g`, where there is no gravity to account for, a free
  root is simply force-free and the answer needs one offset put back at
  the end. The offset is in the root's own coordinates, which is the
  thing to be careful about: left in the world's it is right only while
  the root is unrotated, and a body that falls and spins at the same time
  quietly stops conserving sideways momentum.

  Having both is also how either is trusted. The mass matrix can be read
  out of inverse dynamics a column at a time -- unit acceleration on one
  joint, none on the others, no gravity, no velocity -- and the bias by
  asking for no acceleration at all. Solve `M qdd = tau - bias` with a
  general matrix inverse and the answer must be what the articulated body
  algorithm produced without ever forming M. The two share the spatial
  algebra and nothing else, so agreeing to fourteen digits is worth
  something."
  (:require [allgo.array :as a]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.numerics.linear :as lin]
            [allgo.physics.contact :as contact]
            [allgo.physics.rigid :as rigid]))

;; ---------------------------------------------------------------------------
;; Spatial algebra

(defn skew
  "The 3x3 matrix with `skew(a) b` = `a x b`."
  [[x y z]]
  (let [x (double x) y (double y) z (double z)]
    [[0.0 (- z) y]
     [z 0.0 (- x)]
     [(- y) x 0.0]]))

(defn- blocks
  "A 6x6 from four 3x3 blocks, reading top-left, top-right, bottom-left,
  bottom-right."
  [a b c d]
  (into (mapv into a b) (mapv into c d)))

(def ^:private zero3 [[0.0 0.0 0.0] [0.0 0.0 0.0] [0.0 0.0 0.0]])

(defn transform
  "The Plucker transform taking motion vectors from a frame to one
  displaced from it by `pos` and rotated by `rot`.

  `rot` turns vectors of the parent frame into the child's orientation,
  so the coordinate transform is its inverse -- the usual place to lose a
  sign. `pos` is where the child's origin sits, in parent coordinates.

  The transpose of this takes force vectors the other way, from the child
  back to the parent, and that duality is the whole reason the inward
  pass of both algorithms is as short as it is."
  [{:keys [rot pos]}]
  (let [e (q/to-matrix (q/conjugate (or rot q/identity-q)))
        r (or pos v/zero)]
    (blocks e zero3
            (lin/mat-scale (lin/mat-mul e (skew r)) -1.0)
            e)))

(defn crm
  "Spatial cross product for motion vectors: `crm(v) m` = `v x m`."
  [[wx wy wz vx vy vz]]
  (let [w (skew [wx wy wz])]
    (blocks w zero3 (skew [vx vy vz]) w)))

(defn crf
  "Spatial cross product for force vectors, `-crm(v)^T`.

  Not the same operator and not an optimisation: a force vector is a
  different kind of thing from a motion vector -- it pairs with one to
  give power -- and carrying it along a moving frame is a different sum."
  [[wx wy wz vx vy vz]]
  (let [w (skew [wx wy wz])]
    (blocks w (skew [vx vy vz]) zero3 w)))

(defn spatial-inertia
  "The 6x6 inertia of a body of `mass` whose centre of mass is at `com`
  and whose rotational inertia about that centre is `inertia`.

  About the *frame origin*, not the centre of mass, which is what lets
  every link use its own joint frame and never think about the parallel
  axis theorem again -- the off-diagonal blocks are it."
  [mass com inertia]
  (let [m (double mass)
        c (skew com)
        ct (lin/transpose c)]
    (blocks (lin/mat-add inertia (lin/mat-scale (lin/mat-mul c ct) m))
            (lin/mat-scale c m)
            (lin/mat-scale ct m)
            (lin/mat-scale (lin/eye 3) m))))

;; ---------------------------------------------------------------------------
;; Joints

;; ---------------------------------------------------------------------------
;; The same algebra, flat
;;
;; Everything above is the readable statement of what a spatial vector is
;; and what may be done to one, and it is what the tests check. It is not
;; what the inner loops run on. A 6x6 as a vector of vectors is seven
;; objects and thirty-six boxed doubles; `mat-vec` over one allocates six
;; lazy sequences to produce a seventh vector. A settled sixteen link
;; model asks for several hundred of those a step, which was, measured,
;; most of a frame.
;;
;; So the per-configuration cache holds `double` arrays -- a 6x6 as
;; thirty-six numbers in row-major order, a spatial vector as six -- and
;; the passes over it allocate one result array apiece and nothing else.
;; The conversion happens once per link per step, where it costs nothing.

(defn- flat36
  "A 6x6 from rows into one array, row-major."
  ^doubles [m]
  (let [out (a/f64 36)]
    (dotimes [r 6]
      (let [row (nth m r)]
        (dotimes [c 6] (aset out (+ (* 6 r) c) (double (nth row c))))))
    out))

(defn- e-matrix
  "The rotation `rot` as nine numbers, transposed -- the coordinate
  transform rather than the rotation."
  ^doubles [rot]
  (let [r (or rot q/identity-q)
        [ax ay az] (q/rotate r [1.0 0.0 0.0])
        [bx by bz] (q/rotate r [0.0 1.0 0.0])
        [cx cy cz] (q/rotate r [0.0 0.0 1.0])
        o (a/f64 9)]
    ;; Columns of the rotation are where the axes go, so its transpose
    ;; -- which is what a coordinate transform wants -- has them as rows.
    (aset o 0 (double ax)) (aset o 1 (double ay)) (aset o 2 (double az))
    (aset o 3 (double bx)) (aset o 4 (double by)) (aset o 5 (double bz))
    (aset o 6 (double cx)) (aset o 7 (double cy)) (aset o 8 (double cz))
    o))

(defn- transform-flat
  "`transform`, built straight into thirty-six numbers.

  The readable one above goes through four vector-of-vector matrices and
  a `mat-mul` to say the same thing, and `links` asks for two of them
  per link per step. On a settled ragdoll that was four milliseconds of
  a forty-two millisecond step -- which is what comes of writing that a
  conversion `costs nothing` without measuring it."
  ^doubles [{:keys [rot pos]}]
  (let [^doubles e (e-matrix rot)
        [rx ry rz] (or pos v/zero)
        rx (double rx) ry (double ry) rz (double rz)
        o (a/f64 36)]
    (dotimes [r 3]
      (let [b (* 3 r)
            e0 (aget e b) e1 (aget e (+ b 1)) e2 (aget e (+ b 2))]
        ;; E in both diagonal blocks.
        (aset o (+ (* 6 r) 0) e0) (aset o (+ (* 6 r) 1) e1) (aset o (+ (* 6 r) 2) e2)
        (aset o (+ (* 6 (+ r 3)) 3) e0) (aset o (+ (* 6 (+ r 3)) 4) e1)
        (aset o (+ (* 6 (+ r 3)) 5) e2)
        ;; Bottom left is -E skew(r).
        (let [row (* 6 (+ r 3))]
          (aset o row (- (- (* e1 rz) (* e2 ry))))
          (aset o (+ row 1) (- (- (* e2 rx) (* e0 rz))))
          (aset o (+ row 2) (- (- (* e0 ry) (* e1 rx)))))))
    o))

(defn- spatial-inertia-flat
  "`spatial-inertia`, built straight into thirty-six numbers."
  ^doubles [mass com inertia]
  (let [m (double mass)
        [cx cy cz] com
        cx (double cx) cy (double cy) cz (double cz)
        o (a/f64 36)
        ;; C = skew(c); the upper left is Ic + m C C^T, and C C^T is
        ;; (c.c) I - c c^T.
        cc (+ (* cx cx) (* cy cy) (* cz cz))
        cvec [cx cy cz]]
    (dotimes [r 3]
      (dotimes [k 3]
        (aset o (+ (* 6 r) k)
              (+ (double (nth (nth inertia r) k))
                 (* m (- (if (= r k) cc 0.0)
                         (* (double (nth cvec r)) (double (nth cvec k)))))))))
    ;; m C top right, m C^T bottom left, m on the bottom right diagonal.
    (let [c [[0.0 (- cz) cy] [cz 0.0 (- cx)] [(- cy) cx 0.0]]]
      (dotimes [r 3]
        (dotimes [k 3]
          (aset o (+ (* 6 r) 3 k) (* m (double (nth (nth c r) k))))
          (aset o (+ (* 6 (+ r 3)) k) (* m (double (nth (nth c k) r)))))))
    (dotimes [r 3] (aset o (+ (* 6 (+ r 3)) 3 r) m))
    o))

(defn- mv6
  "`m x` for a flat 6x6 and a flat spatial vector."
  ^doubles [^doubles m ^doubles x]
  (let [out (a/f64 6)]
    (dotimes [r 6]
      (let [b (* 6 r)]
        (aset out r (+ (* (aget m b) (aget x 0))
                       (* (aget m (+ b 1)) (aget x 1))
                       (* (aget m (+ b 2)) (aget x 2))
                       (* (aget m (+ b 3)) (aget x 3))
                       (* (aget m (+ b 4)) (aget x 4))
                       (* (aget m (+ b 5)) (aget x 5))))))
    out))

(defn- mm6
  "`a b` for two flat 6x6s."
  ^doubles [^doubles a ^doubles b]
  (let [out (a/f64 36)]
    (dotimes [r 6]
      (dotimes [c 6]
        (let [ar (* 6 r)]
          (aset out (+ ar c)
                (+ (* (aget a ar) (aget b c))
                   (* (aget a (+ ar 1)) (aget b (+ 6 c)))
                   (* (aget a (+ ar 2)) (aget b (+ 12 c)))
                   (* (aget a (+ ar 3)) (aget b (+ 18 c)))
                   (* (aget a (+ ar 4)) (aget b (+ 24 c)))
                   (* (aget a (+ ar 5)) (aget b (+ 30 c))))))))
    out))

(defn- transpose6
  ^doubles [^doubles m]
  (let [out (a/f64 36)]
    (dotimes [r 6] (dotimes [c 6] (aset out (+ (* 6 c) r) (aget m (+ (* 6 r) c)))))
    out))

(defn- add6! ^doubles [^doubles dst ^doubles x] (dotimes [i 6] (aset dst i (+ (aget dst i) (aget x i)))) dst)
(defn- add36! ^doubles [^doubles dst ^doubles x] (dotimes [i 36] (aset dst i (+ (aget dst i) (aget x i)))) dst)

(defn- dot6 ^double [^doubles a ^doubles b]
  (+ (* (aget a 0) (aget b 0)) (* (aget a 1) (aget b 1)) (* (aget a 2) (aget b 2))
     (* (aget a 3) (aget b 3)) (* (aget a 4) (aget b 4)) (* (aget a 5) (aget b 5))))

(defn- copy6 ^doubles [^doubles x] (let [o (a/f64 6)] (dotimes [i 6] (aset o i (aget x i))) o))
(defn- copy36 ^doubles [^doubles x] (let [o (a/f64 36)] (dotimes [i 36] (aset o i (aget x i))) o))

(defn- cross-motion
  "`crm(v) x`, written out rather than built.

  The 6x6 of `crm` never has to exist: the angular half of the answer is
  `w x a` and the linear half is `u x a + w x b`, which is twelve
  multiplications against thirty-six and no matrix."
  ^doubles [^doubles v ^doubles x]
  (let [wx (aget v 0) wy (aget v 1) wz (aget v 2)
        ux (aget v 3) uy (aget v 4) uz (aget v 5)
        ax (aget x 0) ay (aget x 1) az (aget x 2)
        bx (aget x 3) by (aget x 4) bz (aget x 5)
        out (a/f64 6)]
    (aset out 0 (- (* wy az) (* wz ay)))
    (aset out 1 (- (* wz ax) (* wx az)))
    (aset out 2 (- (* wx ay) (* wy ax)))
    (aset out 3 (+ (- (* uy az) (* uz ay)) (- (* wy bz) (* wz by))))
    (aset out 4 (+ (- (* uz ax) (* ux az)) (- (* wz bx) (* wx bz))))
    (aset out 5 (+ (- (* ux ay) (* uy ax)) (- (* wx by) (* wy bx))))
    out))

(defn- cross-force
  "`crf(v) f`, likewise. Angular is `w x n + u x f`, linear is `w x f`."
  ^doubles [^doubles v ^doubles f]
  (let [wx (aget v 0) wy (aget v 1) wz (aget v 2)
        ux (aget v 3) uy (aget v 4) uz (aget v 5)
        nx (aget f 0) ny (aget f 1) nz (aget f 2)
        fx (aget f 3) fy (aget f 4) fz (aget f 5)
        out (a/f64 6)]
    (aset out 0 (+ (- (* wy nz) (* wz ny)) (- (* uy fz) (* uz fy))))
    (aset out 1 (+ (- (* wz nx) (* wx nz)) (- (* uz fx) (* ux fz))))
    (aset out 2 (+ (- (* wx ny) (* wy nx)) (- (* ux fy) (* uy fx))))
    (aset out 3 (- (* wy fz) (* wz fy)))
    (aset out 4 (- (* wz fx) (* wx fz)))
    (aset out 5 (- (* wx fy) (* wy fx)))
    out))

(defn- s-dot
  "`S^T f`, as a flat `ndof` array -- each of a joint's axes' share."
  ^doubles [^doubles s nd ^doubles f]
  (let [nd (long nd) out (a/f64 nd)]
    (dotimes [k nd]
      (let [b (* 6 k)]
        (aset out k (+ (* (aget s b) (aget f 0)) (* (aget s (+ b 1)) (aget f 1))
                       (* (aget s (+ b 2)) (aget f 2)) (* (aget s (+ b 3)) (aget f 3))
                       (* (aget s (+ b 4)) (aget f 4)) (* (aget s (+ b 5)) (aget f 5))))))
    out))

(defn- s-apply!
  "`dst += S x`, taking the joint's rates from `x` at `off`."
  ^doubles [^doubles dst ^doubles s nd x off]
  (dotimes [k (long nd)]
    (let [r (double (nth x (+ (long off) k))) b (* 6 k)]
      (when-not (zero? r)
        (dotimes [i 6] (aset dst i (+ (aget dst i) (* r (aget s (+ b i)))))))))
  dst)

(defn- s-apply-arr!
  "`dst += S x` for a flat `x` of exactly this joint's rates."
  ^doubles [^doubles dst ^doubles s nd ^doubles x]
  (dotimes [k nd]
    (let [r (aget x k) b (* 6 k)]
      (when-not (zero? r)
        (dotimes [i 6] (aset dst i (+ (aget dst i) (* r (aget s (+ b i)))))))))
  dst)

(defn- flat-d-inverse
  "The inverse of a joint's `D`, or nil where it has none.

  Nil rather than a huge number: a joint whose subtree presents no
  inertia along one of its axes has no acceleration there to compute."
  ^doubles [^doubles d nd]
  (case (long nd)
    0 nil
    1 (let [x (aget d 0)] (when (> (abs x) 1e-12) (a/f64 [(/ 1.0 x)])))
    (let [rows (mapv (fn [r] (mapv (fn [c] (aget d (+ (* nd r) c))) (range nd))) (range nd))]
      (when-let [inv (lin/inverse rows)]
        (a/f64 (mapcat identity inv))))))

(defn- minus-u-dinv-ut
  "`IA - U D^-1 U^T`, the inertia a joint presents once it has been
  allowed to give."
  ^doubles [^doubles ia ^doubles u ^doubles dinv nd]
  (let [nd (long nd) out (copy36 ia)]
    (dotimes [k nd]
      (dotimes [l nd]
        (let [w (aget dinv (+ (* nd k) l))]
          (when-not (zero? w)
            (let [bk (* 6 k) bl (* 6 l)]
              (dotimes [r 6]
                (let [uk (* w (aget u (+ bk r)))]
                  (when-not (zero? uk)
                    (dotimes [c 6]
                      (aset out (+ (* 6 r) c)
                            (- (aget out (+ (* 6 r) c)) (* uk (aget u (+ bl c))))))))))))))
    out))

(defn- dot-n
  "The dot product of two flat arrays of the same length."
  ^double [^doubles a ^doubles b]
  (let [n (alength a)]
    (loop [i 0 acc 0.0] (if (= i n) acc (recur (inc i) (+ acc (* (aget a i) (aget b i))))))))

(defn- axpy-n!
  "`dst += s x` over a whole flat array."
  ^doubles [^doubles dst ^doubles x s]
  (let [s (double s) n (alength dst)]
    (when-not (zero? s)
      (dotimes [i n] (aset dst i (+ (aget dst i) (* s (aget x i))))))
    dst))

(defn- root-solve
  "`IA0^-1 p` for the free root's six by six.

  Symmetric by construction and not quite by arithmetic, after a chain
  of congruences, so it is averaged with its own transpose before being
  factored. Cholesky wants it to be and it costs thirty-six additions."
  [^doubles ia0 ^doubles p]
  (let [rows (mapv (fn [r] (mapv (fn [c]
                                   (* 0.5 (+ (aget ia0 (+ (* 6 r) c))
                                             (aget ia0 (+ (* 6 c) r)))))
                                 (range 6)))
                   (range 6))]
    (lin/cholesky-solve rows (mapv (fn [i] (aget p (long i))) (range 6)))))

(defn- small-solve
  "`D^-1 y` for the one by one or three by three a joint presents.

  Nil `dinv` is a joint whose subtree resists nothing along one of its
  axes; there is no acceleration there to compute and inventing one is
  how a zero-inertia limb tears a model apart."
  ^doubles [^doubles dinv nd ^doubles y]
  (let [nd (long nd) out (a/f64 nd)]
    (dotimes [r nd]
      (let [b (* nd r)]
        (aset out r (double (loop [c 0 acc 0.0]
                              (if (= c nd)
                                acc
                                (recur (inc c) (+ acc (* (aget dinv (+ b c)) (aget y c))))))))))
    out))

(defn- kind-of [link] (or (:joint link) :revolute))

(defn joint-dof
  "How many numbers it takes to say how fast this joint is moving.

  One for a hinge or a slider, three for a ball, none at all for a weld.
  A joint's *coordinate* is a separate question and not always the same
  count: a ball joint's velocity is three numbers and its configuration
  is four, because three numbers cannot name a rotation without a
  direction in which they stop working."
  [link]
  (case (kind-of link)
    (:revolute :prismatic) 1
    :spherical 3
    :fixed 0))

(defn- joint-transform
  "How the child frame sits relative to the joint frame at coordinate `x`.

  For a ball joint `x` *is* the rotation, carried as a quaternion. Every
  other kind takes a single number."
  [link x]
  (case (kind-of link)
    :revolute {:rot (q/from-axis-angle (:axis link) (double x)) :pos v/zero}
    :prismatic {:rot q/identity-q :pos (v/scale (:axis link) (double x))}
    :spherical {:rot (if (and (sequential? x) (= 4 (count x)))
                       (q/normalize (vec x))
                       q/identity-q)
                :pos v/zero}
    :fixed {:rot q/identity-q :pos v/zero}))

(defn- subspace
  "The joint's motion subspace: one column per degree of freedom, each
  the spatial velocity that one unit of that rate produces, in the
  child's own coordinates.

  A ball joint's three columns are the identity on the angular half,
  which is what makes it cheaper than three stacked hinges rather than
  merely tidier -- there is no intermediate frame to transform through
  and no configuration in which two of the three axes line up."
  [link]
  (let [[ax ay az] (:axis link)]
    (case (kind-of link)
      :revolute [[(double ax) (double ay) (double az) 0.0 0.0 0.0]]
      :prismatic [[0.0 0.0 0.0 (double ax) (double ay) (double az)]]
      :spherical [[1.0 0.0 0.0 0.0 0.0 0.0]
                  [0.0 1.0 0.0 0.0 0.0 0.0]
                  [0.0 0.0 1.0 0.0 0.0 0.0]]
      :fixed [])))

(defn chain
  "The links of a model, whichever form it was given in.

  A bare vector of links is a chain bolted to the world. A map
  `{:base {:mass :com :inertia} :links [...]}` gives it a root body that
  is free to move, which is what a ragdoll's pelvis is and what anything
  that has to fall over needs."
  [model]
  (if (map? model) (:links model) model))

(defn base
  "The free root body, or nil if the model is bolted to the world."
  [model]
  (when (map? model) (:base model)))

(defn dof
  "How many numbers it takes to say how fast the model's joints are
  moving: one per hinge or slider, three per ball joint, none per weld.

  The base's six, if it has any, are not among them -- they are a pose
  and a velocity, not coordinates, precisely so that no one has to
  parameterise a rotation with three numbers."
  [model]
  (reduce + (map joint-dof (chain model))))

;; ---------------------------------------------------------------------------
;; Shared per-link setup

(deftype Link [^doubles xup ^doubles xt ^doubles s ^long ndof ^long offset
               ^doubles inertia ^long parent])

(deftype Node [^doubles ia ^doubles u dinv ^doubles ia-free])

;; `deftype` rather than a map, and it is worth saying why because a map
;; reads better. The two above are what every pass looks at once per
;; link, and a contact asks for three passes -- so on a settled ragdoll
;; they were read some six hundred times a step. A map destructured by
;; keyword is a hash lookup apiece; a `deftype` field is a property
;; read. On a JVM the difference is lost beside the floating point. In
;; JavaScript it was most of the cost of a step.

(defn- links
  "Each link's transform from its parent, motion subspace and inertia, at
  the configuration `q-vec`."
  [model q-vec]
  (let [parts (chain model)]
    (first
     (reduce
      (fn [[acc ^long off] i]
        (let [link (nth parts i)
              nd (long (joint-dof link))
              xup (mm6 (transform-flat (joint-transform link (nth q-vec i nil)))
                       (transform-flat (:origin link)))
              cols (subspace link)
              sarr (a/f64 (* 6 (max 1 nd)))]
          (dotimes [k nd]
            (dotimes [j 6] (aset sarr (+ (* 6 k) j) (double (nth (nth cols k) j)))))
          [(conj acc (Link. xup
                             ;; The transpose is kept rather than taken
                             ;; again: every inward pass, of which there
                             ;; are several a step, wants it.
                            (transpose6 xup)
                            sarr
                            nd
                             ;; Where this joint's rates start in the
                             ;; flat velocity vector. Joints are no
                             ;; longer all one number wide, so nothing
                             ;; can index by link any more.
                            off
                            (spatial-inertia-flat (:mass link) (:com link) (:inertia link))
                            (long (:parent link))))
           (+ off nd)]))
      [[] 0]
      (range (count parts))))))

(defn- total-dof
  "How many joint rates the model has, off the last link.

  Offsets accumulate in link order, so the last one plus its own width
  is the total -- which is worth an O(1) read rather than a walk,
  because the inner recursion asks for it on every call and there are
  some seventy of those in a settled ragdoll's step."
  ^long [ls]
  (let [n (count ls)]
    (if (zero? n)
      0
      (let [^Link l (nth ls (dec n))] (+ (.-offset l) (.-ndof l))))))

(defn- spatial-velocities
  "Every link's spatial velocity, in its own coordinates.

  Linear in `(v0, qd)` and used three ways because of it: to say what a
  model is doing now, and -- one unit of generalised velocity at a time
  -- to build the Jacobian of any point on it."
  [ls ^doubles v0 qd]
  (reduce (fn [acc i]
            (let [^Link l (nth ls i)
                  parent (.-parent l)
                  vp (if (neg? parent) v0 (nth acc parent))]
              (conj acc (s-apply! (mv6 (.-xup l) vp) (.-s l) (.-ndof l) qd (.-offset l)))))
          []
          (range (count ls))))

(defn- base-acceleration
  "The base's spatial acceleration, which is how gravity gets in.

  Not by adding a force to every link. A fictitious upward acceleration
  of the whole base produces exactly the same relative motion and costs
  nothing -- gravity then arrives through the same recursion as
  everything else, and no link has to be told about it."
  [gravity]
  (let [[gx gy gz] (or gravity [0.0 -9.81 0.0])]
    [0.0 0.0 0.0 (- (double gx)) (- (double gy)) (- (double gz))]))

;; ---------------------------------------------------------------------------
;; Inverse dynamics

(defn inverse-dynamics
  "The torques that produce accelerations `qdd` from state `q-vec`, `qd`.

  Recursive Newton-Euler: outward for the velocity and acceleration of
  every link, then the net force each one needs; inward summing those
  forces onto the parents, taking the joint's share off each as it
  passes."
  ([model q-vec qd qdd] (inverse-dynamics model q-vec qd qdd nil))
  ([model q-vec qd qdd {:keys [gravity]}]
   (let [ls (links model q-vec)
         n (count ls)
         ;; Outward: where every link is going, and what it takes.
         out (reduce
              (fn [acc i]
                (let [^Link l (nth ls i)
                      parent (.-parent l)
                      ^doubles inertia (.-inertia l)
                      vj (s-apply! (a/f64 6) (.-s l) (.-ndof l) qd (.-offset l))
                      vp (if (neg? parent) (a/f64 6) (:v (nth acc parent)))
                      ap (if (neg? parent)
                           (a/f64 (base-acceleration gravity))
                           (:a (nth acc parent)))
                      vi (add6! (mv6 (.-xup l) vp) vj)
                      ai (-> (mv6 (.-xup l) ap)
                             (s-apply! (.-s l) (.-ndof l) qdd (.-offset l))
                             (add6! (cross-motion vi vj)))
                      fi (add6! (mv6 inertia ai) (cross-force vi (mv6 inertia vi)))]
                  (conj acc {:v vi :a ai :f fi})))
              []
              (range n))
         ;; Inward: each link's force, plus everything its children left
         ;; on it, read off along the joint axis.
         forces (reduce (fn [fs i]
                          (let [^Link l (nth ls i)
                                fi (nth fs i)]
                            (if (neg? (.-parent l))
                              fs
                              (update fs (.-parent l) #(add6! % (mv6 (.-xt l) fi))))))
                        (mapv #(copy6 (:f %)) out)
                        (reverse (range n)))]
     ;; Each joint's share of the force on its own link, laid back into
     ;; the flat vector the caller handed its rates in.
     (reduce (fn [acc i]
               (let [^Link l (nth ls i)
                     ^doubles share (s-dot (.-s l) (.-ndof l) (nth forces i))]
                 (reduce (fn [a k] (assoc a (+ (.-offset l) k) (aget share k)))
                         acc
                         (range (.-ndof l)))))
             (vec (repeat (dof model) 0.0))
             (range n)))))

(defn mass-matrix
  "The joint-space inertia `M(q)`, a column at a time out of
  `inverse-dynamics`.

  Unit acceleration on one joint, none on any other, nothing moving and
  no gravity: what is left is that column of M. It is O(n^2) calls of an
  O(n) algorithm where the composite rigid body algorithm would be
  O(n^2) outright, and it is here to check the fast path rather than to
  be the fast path."
  [model q-vec]
  (let [n (dof model)
        zero (vec (repeat n 0.0))]
    (lin/transpose
     (mapv (fn [i]
             (inverse-dynamics model q-vec zero (assoc zero i 1.0)
                               {:gravity [0.0 0.0 0.0]}))
           (range n)))))

(defn bias-forces
  "The torques needed to hold the joints at zero acceleration: gravity,
  centrifugal and Coriolis together."
  ([model q-vec qd] (bias-forces model q-vec qd nil))
  ([model q-vec qd opts]
   (inverse-dynamics model q-vec qd (vec (repeat (dof model) 0.0)) opts)))

;; ---------------------------------------------------------------------------
;; Forward dynamics

(defn- articulated-inertias
  "Each joint's articulated inertia and the two numbers taken off it,
  plus the root's.

  `IA` is what everything beyond a joint weighs given that its own
  joints are free to move; `U = IA S` and `d = S^T U` project that onto
  the joint's own axis, and `ia-free` is `IA` with the axis divided out
  -- what the subtree presents to the link above once the joint between
  them has been allowed to give.

  None of it depends on how fast anything is going or on what is pushing
  it, only on where the joints are. That is what lets one build serve
  both the step's accelerations and every contact impulse asked about
  afterwards, and it is the difference between a contact costing O(n)
  and costing a whole inverse inertia matrix."
  [model ls]
  (let [n (count ls)
        root (base model)
        ;; Filled in from the tips down. Every `ia` starts as the link's
        ;; own and is added to in place by its children, so there is
        ;; nothing to thread through a reduce -- the arrays are the
        ;; accumulator.
        nodes (object-array n)
        ia0 (when root (spatial-inertia-flat (:mass root) (:com root) (:inertia root)))]
    (dotimes [k n]
      (let [^Link l (nth ls k)]
        (aset nodes k (Node. (copy36 (.-inertia l)) nil nil nil))))
    (loop [i (dec n)]
      (when (>= i 0)
        (let [^Link l (nth ls i)
              ^Node node (aget nodes i)
              ^doubles s (.-s l)
              nd (.-ndof l)
              ^doubles ia (.-ia node)
              ;; One column of U per degree of freedom, and D their
              ;; projection onto each other -- a number for a hinge, a
              ;; three by three for a ball joint.
              u (a/f64 (* 6 (max 1 nd)))
              _ (dotimes [k nd]
                  (let [col (a/f64 6)]
                    (dotimes [j 6] (aset col j (aget s (+ (* 6 k) j))))
                    (let [^doubles uk (mv6 ia col)]
                      (dotimes [j 6] (aset u (+ (* 6 k) j) (aget uk j))))))
              d (let [out (a/f64 (* (max 1 nd) (max 1 nd)))]
                  (dotimes [r nd]
                    (dotimes [c nd]
                      (let [br (* 6 r) bc (* 6 c)]
                        (aset out (+ (* nd r) c)
                              (double (loop [j 0 acc 0.0]
                                        (if (= j 6)
                                          acc
                                          (recur (inc j)
                                                 (+ acc (* (aget s (+ br j))
                                                           (aget u (+ bc j))))))))))))
                  out)
              dinv (flat-d-inverse d nd)
              ia' (if dinv (minus-u-dinv-ut ia u dinv nd) ia)
              up (mm6 (.-xt l) (mm6 ia' (.-xup l)))]
          (aset nodes i (Node. ia u dinv ia'))
          (if (neg? (.-parent l))
            (when ia0 (add36! ia0 up))
            (let [^Node up-node (aget nodes (.-parent l))]
              (add36! (.-ia up-node) up))))
        (recur (dec i))))
    {:links nodes :ia0 ia0}))

(defn forward-dynamics
  "The accelerations that torques `tau` produce, in O(n).

  Returns `{:qdd :base-acc}` -- the joint accelerations, and the base's
  spatial acceleration, which is zero for a model bolted to the world.

  Featherstone's articulated body algorithm. The idea it turns on: seen
  from one joint, everything beyond it behaves like a single rigid body
  with an inertia of its own -- an *articulated* inertia, which is what
  the subtree's own joints being free to move does to the inertia it
  presents. Compute that inward from the tips and every joint's
  acceleration follows from its parent's outward, with no matrix to
  invert bigger than one degree of freedom.

  Three passes: outward for velocities and the bias force each link's own
  motion generates, inward accumulating articulated inertias and their
  bias forces onto the parents, outward again for the accelerations.

  A free base costs one 6x6 solve and no change to any of that. The
  inward pass carries the whole body's articulated inertia and bias force
  up to the root, and a root nothing is holding is a root with no force
  on it, so its acceleration is whatever makes that sum vanish. The
  bolted case is the same equation with the answer already known."
  ([model q-vec qd tau] (forward-dynamics model q-vec qd tau nil))
  ([model q-vec qd tau {gravity :gravity root-state :base
                        ext :external root-force :base-force
                        given-ls :ls given-ai :ai}]
   (let [ls (or given-ls (links model q-vec))
         n (count ls)
         root (base model)
         ;; Everything is computed in a frame falling at `gravity`, where
         ;; there is nothing to account for and a free root is simply
         ;; force-free. The offset comes back out at the end.
         ;;
         ;; In the root's *own* coordinates, which is the part that is
         ;; easy to get wrong: every spatial quantity here is in the
         ;; frame it belongs to, so gravity has to be turned into the
         ;; base's before it can be added to anything. Left in world
         ;; coordinates it is right only while the base is unrotated,
         ;; and a falling body that is also spinning quietly stops
         ;; conserving horizontal momentum.
         a-grav (let [[gx gy gz] (if root
                                   (q/rotate (q/conjugate (or (:rot root-state)
                                                              q/identity-q))
                                             (or gravity [0.0 -9.81 0.0]))
                                   (or gravity [0.0 -9.81 0.0]))]
                  (a/f64 [0.0 0.0 0.0 (double gx) (double gy) (double gz)]))
         v0 (a/f64 (if root (vec (or (:vel root-state) (repeat 6 0.0))) (repeat 6 0.0)))
         root-i (when root
                  (spatial-inertia-flat (:mass root) (:com root) (:inertia root)))
         ;; What every joint weighs from above, which depends only on
         ;; where the joints are. Handed in when the caller has already
         ;; built it for this configuration.
         ai (or given-ai (articulated-inertias model ls))
         ;; Pass one, outward: velocity, the acceleration that velocity
         ;; alone implies, and the force needed to hold the link on it.
         pass1 (reduce
                (fn [acc i]
                  (let [^Link l (nth ls i)
                        parent (.-parent l)
                        ^doubles inertia (.-inertia l)
                        vj (s-apply! (a/f64 6) (.-s l) (.-ndof l) qd (.-offset l))
                        vp (if (neg? parent) v0 (:v (nth acc parent)))
                        vi (add6! (mv6 (.-xup l) vp) vj)
                        ;; An outside force on a link reduces the bias
                        ;; force the link needs of its parent by exactly
                        ;; itself, which is the whole of how anything
                        ;; external gets in -- a contact, a thruster, a
                        ;; hand pushing.
                        f (when ext (nth ext i nil))
                        bias (cross-force vi (mv6 inertia vi))]
                    (when f (dotimes [k 6] (aset bias k (- (aget bias k) (double (nth f k))))))
                    (conj acc {:v vi
                               :c (cross-motion vi vj)
                               :pa bias})))
                []
                (range n))
         ;; Pass two, inward: what the subtree beyond each joint looks
         ;; like to the link above it, and what the whole of it looks
         ;; like to the root.
         ;; Pass two, inward: the bias forces only. The inertias they
         ;; ride on were built once, above.
         pass2 (reduce
                (fn [{:keys [links pa0] :as acc} i]
                  (let [^Link l (nth ls i)
                        nd (.-ndof l)
                        parent (.-parent l)
                        offset (.-offset l)
                        ^doubles s (.-s l)
                        {:keys [^doubles pa ^doubles c]} (nth links i)
                        ^Node node (aget ^objects (:links ai) i)
                        ^doubles u (.-u node)
                        dinv (.-dinv node)
                        ^doubles ia-free (.-ia-free node)
                        ;; What each of this joint's axes has left over
                        ;; once the bias force has taken its share.
                        uu (let [^doubles share (s-dot s nd pa) out (a/f64 (max 1 nd))]
                             (dotimes [k nd]
                               (aset out k (- (double (nth tau (+ (long offset) k)))
                                              (aget share k))))
                             out)
                        acc (assoc-in acc [:links i :uu] uu)
                        pa' (cond-> (add6! (mv6 ia-free c) pa)
                              dinv (s-apply-arr! u nd (small-solve dinv nd uu)))
                        up-p (mv6 (.-xt l) pa')]
                    (if (neg? parent)
                      (assoc acc :pa0 (when pa0 (add6! pa0 up-p)))
                      (update-in acc [:links parent :pa] #(add6! % up-p)))))
                {:links pass1
                 :pa0 (when root
                        (let [bias (cross-force v0 (mv6 root-i v0))]
                          (when root-force
                            (dotimes [k 6]
                              (aset bias k (- (aget bias k) (double (nth root-force k))))))
                          bias))}
                (reverse (range n)))
         solved (:links pass2)
         ;; The root, if it is free: force-free in the falling frame.
         a0 (if root
              (a/f64 (mapv - (root-solve (:ia0 ai) (:pa0 pass2))))
              (let [o (a/f64 6)] (dotimes [k 6] (aset o k (- (aget a-grav k)))) o))
         ;; Pass three, outward: each joint's acceleration, then the link's.
         out (reduce
              (fn [{:keys [a] :as acc} i]
                (let [^Link l (nth ls i)
                      nd (.-ndof l)
                      ^doubles s (.-s l)
                      {:keys [^doubles c ^doubles uu]} (nth solved i)
                      ^Node node (aget ^objects (:links ai) i)
                      ^doubles u (.-u node)
                      dinv (.-dinv node)
                      ap (if (neg? (.-parent l)) a0 (nth a (.-parent l)))
                      ^doubles a' (add6! (mv6 (.-xup l) ap) c)
                      qddi (if dinv
                             (small-solve dinv nd
                                          (let [y (a/f64 (max 1 nd))]
                                            (dotimes [k nd]
                                              (let [b (* 6 k)]
                                                (aset y k (- (aget uu k)
                                                             (loop [j 0 s2 0.0]
                                                               (if (= j 6)
                                                                 s2
                                                                 (recur (inc j)
                                                                        (+ s2 (* (aget u (+ b j))
                                                                                 (aget a' j))))))))))
                                            y))
                             (a/f64 (max 1 nd)))]
                  (-> acc
                      (update :a conj (s-apply-arr! (copy6 a') s nd qddi))
                      (update :qdd into (take nd (seq qddi))))))
              {:a [] :qdd []}
              (range n))]
     {:qdd (vec (:qdd out))
      ;; Back out of the falling frame, which is where the root's own
      ;; share of gravity comes from.
      :base-acc (if root
                  (mapv (fn [i] (+ (aget ^doubles a0 (long i)) (aget ^doubles a-grav (long i))))
                        (range 6))
                  (vec (repeat 6 0.0)))})))

;; ---------------------------------------------------------------------------
;; Moving it

(declare solve-constraints contacts-with self-contacts velocity with-velocity
         poses limited?)

(defn- advance-coordinate
  "One joint's configuration, moved by its own rates.

  A hinge or a slider adds. A ball joint cannot: its configuration is a
  quaternion and its rates are an angular velocity in the *child's* own
  frame, so the step is `qj += dt/2 qj (w,0)` -- the joint's own
  quaternion on the left, where the base's world-frame spin puts it on
  the right. Getting that round the wrong way turns a ball joint into a
  thing that drifts sideways under its own rotation."
  [link x qd offset dt]
  (let [offset (long offset) dt (double dt)]
    (case (kind-of link)
      (:revolute :prismatic) (+ (double x) (* dt (double (nth qd offset))))
      :spherical (let [rot (if (and (sequential? x) (= 4 (count x)))
                             (vec x)
                             q/identity-q)
                       w [(double (nth qd offset))
                          (double (nth qd (+ offset 1)))
                          (double (nth qd (+ offset 2)))]
                       [dx dy dz dw] (q/mul rot [(w 0) (w 1) (w 2) 0.0])
                       [rx ry rz rw] rot
                       h (* 0.5 dt)]
                   (q/normalize [(+ rx (* h dx)) (+ ry (* h dy))
                                 (+ rz (* h dz)) (+ rw (* h dw))]))
      :fixed x)))

(defn- advance-positions
  "Move by whatever the velocities now are."
  [model state ^double dt]
  (let [root (base model)
        b (:base state)
        parts (chain model)
        ;; Where each joint's rates start, and nothing else. Building the
        ;; whole per-configuration cache to read ten integers off it cost
        ;; four milliseconds a step, which was a tenth of the frame.
        offsets (reductions + 0 (map joint-dof parts))
        q' (mapv (fn [i]
                   (advance-coordinate (nth parts i) (nth (:q state) i nil)
                                       (:qd state) (nth offsets i) dt))
                 (range (count parts)))]
    (if-not root
      (assoc state :q q')
      (let [v (vec (:vel b))
            rot (or (:rot b) q/identity-q)
            pos (or (:pos b) v/zero)
            ;; Both halves of the spatial velocity are in the base's own
            ;; frame, so both are turned into the world before they move
            ;; anything.
            w (q/rotate rot (subvec v 0 3))
            u (q/rotate rot (subvec v 3 6))
            ;; rot += dt/2 (omega, 0) rot, then back onto the unit
            ;; sphere -- the same first-order step `allgo.physics.rigid`
            ;; takes, and it walks off for the same reason.
            [dx dy dz dw] (q/mul [(w 0) (w 1) (w 2) 0.0] rot)
            [rx ry rz rw] rot
            h (* 0.5 dt)]
        (assoc state
               :q q'
               :base (assoc b
                            :rot (q/normalize [(+ rx (* h dx)) (+ ry (* h dy))
                                               (+ rz (* h dz)) (+ rw (* h dw))])
                            :pos (v/add-scaled pos u dt)))))))

(defn advance
  "Move by whatever the velocities already are, and change nothing else.

  `step` is accelerate, then contacts, then this. A caller that has done
  its own accelerating and its own contacts -- `allgo.physics.world`,
  which has other bodies to account for in the same sweep -- wants only
  the last of the three."
  [state model dt]
  (advance-positions model state dt))

(defn step
  "One semi-implicit Euler step of `{:q :qd}`, and `:base` if there is
  one, under `tau`.

  Velocity first, then contacts, then position. The order is the whole
  of what makes a contact hold: the velocities are corrected before
  anything moves, so a body that was about to be driven into the floor
  never is, rather than being pulled back out afterwards.

  Contacts come either ready-made as `:contacts` or, more usually, from
  `:obstacles` -- static `allgo.physics.rigid` bodies this model is to
  be generated against -- and from `:self-collide?`, which adds the
  model's own parts against each other. They are found at the
  configuration the step *starts* from, which is where the model
  actually is when the question is asked.

  The base's state is a pose and a spatial velocity in its own frame,
  not six more coordinates. Three reasons, and the third is the one that
  matters: a rotation has no three-number parameterisation without a
  singularity in it; a quaternion is the shape the rest of this library
  speaks; and a velocity in body coordinates is what the algorithms
  already produce, so nothing has to be converted on the way in."
  ([state model dt] (step state model dt nil))
  ([{:keys [q qd] :as state} model ^double dt
    {:keys [tau contacts obstacles self-collide?] :as opts}]
   (let [root (base model)
         b (:base state)
         tau (or tau (vec (repeat (dof model) 0.0)))
         ;; Where everything is, computed once and handed to every pass
         ;; that needs it. It depends only on `q` and the root's pose,
         ;; neither of which changes until the very end of the step.
         ls (links model q)
         frames (poses model q b)
         ;; What every joint weighs from above. It depends only on `q`,
         ;; so the accelerations and every contact impulse afterwards
         ;; share one build -- they were making one apiece.
         ai (articulated-inertias model ls)
         {:keys [qdd base-acc]} (forward-dynamics model q qd tau
                                                  (cond-> (assoc opts :ls ls :ai ai)
                                                    root (assoc :base b)))
         accel (if root (vec (concat base-acc qdd)) (vec qdd))
         moving (with-velocity model state
                  (mapv (fn [a c] (+ (double a) (* dt (double c))))
                        (velocity model state)
                        accel))
         cs (or contacts
                (cond-> []
                  (seq obstacles) (into (contacts-with model q b obstacles frames))
                  self-collide? (into (self-contacts model q b frames))))]
     (advance-positions model
                        ;; Limits are solved even with nothing to stand
                        ;; on. A joint folding backwards in mid-air is
                        ;; still a joint folding backwards, and a model
                        ;; with no limits at all pays only this test.
                        (if (or (seq cs) (limited? model))
                          (solve-constraints model moving (vec cs) dt
                                             (assoc opts :ls ls :ai ai :frames frames))
                          moving)
                        dt))))

;; ---------------------------------------------------------------------------
;; Where the links are

(defn poses
  "Each link's frame in world coordinates, as `{:rot :pos}`.

  Reduced coordinates describe an arm by its angles, and something has to
  turn those back into places to draw. `root` is where the base sits, for
  a model that has one; without it everything comes out relative to the
  world origin, which is where a bolted model is."
  ([model q-vec] (poses model q-vec nil))
  ([model q-vec root]
   (:out
    (reduce
     (fn [{:keys [frames] :as acc} i]
       (let [link (nth (chain model) i)
             parent (long (:parent link))
             {prot :rot ppos :pos} (if (neg? parent)
                                     {:rot (or (:rot root) q/identity-q)
                                      :pos (or (:pos root) v/zero)}
                                     (nth frames parent))
             {orot :rot opos :pos} (:origin link)
             {jrot :rot jpos :pos} (joint-transform link (nth q-vec i nil))
             rot (q/mul (q/mul prot (or orot q/identity-q)) jrot)
             pos (v/add (v/add ppos (q/rotate prot (or opos v/zero)))
                        (q/rotate (q/mul prot (or orot q/identity-q)) jpos))
             frame {:rot rot :pos pos}]
         (-> acc (update :frames conj frame) (update :out conj frame))))
     {:frames [] :out []}
     (range (count (chain model)))))))

(defn energy
  "Kinetic plus potential, for checking that nothing is being invented.

  A chain under no torque and no damping must keep this constant, and an
  error in the spatial algebra almost always shows up here before it
  shows up anywhere a person would notice."
  ([model state] (energy model state nil))
  ([model {:keys [q qd] :as state} {:keys [gravity]}]
   (let [g (or gravity [0.0 -9.81 0.0])
         root (base model)
         b (:base state)
         ls (links model q)
         fs (poses model q b)
         v0 (a/f64 (if root (or (:vel b) (repeat 6 0.0)) (repeat 6 0.0)))
         vels (spatial-velocities ls v0 qd)
         part (fn [^doubles inertia ^doubles vi rot pos mass com]
                (let [ke (* 0.5 (dot6 vi (mv6 inertia vi)))
                      c (v/add pos (q/rotate rot com))]
                  (- ke (* (double mass) (v/dot g c)))))]
     (+ (if root
          (part (flat36 (spatial-inertia (:mass root) (:com root) (:inertia root)))
                v0 (or (:rot b) q/identity-q) (or (:pos b) v/zero)
                (:mass root) (:com root))
          0.0)
        (reduce +
                (map (fn [i]
                       (let [link (nth (chain model) i)
                             {:keys [rot pos]} (nth fs i)]
                         (part (let [^Link l (nth ls i)] (.-inertia l))
                               (nth vels i) rot pos (:mass link) (:com link))))
                     (range (count ls))))))))

(defn momentum
  "The whole model's linear and angular momentum about the world origin.

  Under no torque and no gravity both are conserved exactly, whatever the
  joints are doing, and that is the strongest statement available about a
  floating chain -- it is a property of the physics rather than of any
  one link, so a sign error anywhere in the recursion breaks it."
  [model {:keys [q qd] :as state}]
  (let [root (base model)
        b (:base state)
        ls (links model q)
        fs (poses model q b)
        v0 (a/f64 (if root (or (:vel b) (repeat 6 0.0)) (repeat 6 0.0)))
        vels (spatial-velocities ls v0 qd)
        ;; Each body's spatial momentum `I v` is in its own frame and
        ;; about its own origin. Both have to be carried to the world
        ;; before they can be added up.
        one (fn [^doubles inertia ^doubles vi {:keys [rot pos]}]
              (let [^doubles h (mv6 inertia vi)
                    ang (q/rotate rot [(aget h 0) (aget h 1) (aget h 2)])
                    lin (q/rotate rot [(aget h 3) (aget h 4) (aget h 5)])]
                [lin (v/add ang (v/cross pos lin))]))
        parts (cond-> (mapv (fn [i] (one (let [^Link l (nth ls i)] (.-inertia l))
                                         (nth vels i) (nth fs i)))
                            (range (count ls)))
                root (conj (one (flat36 (spatial-inertia (:mass root) (:com root)
                                                         (:inertia root)))
                                v0
                                {:rot (or (:rot b) q/identity-q)
                                 :pos (or (:pos b) v/zero)})))]
    {:linear (reduce v/add v/zero (map first parts))
     :angular (reduce v/add v/zero (map second parts))}))

;; ---------------------------------------------------------------------------
;; Being pushed

(defn generalised-dof
  "How many numbers it takes to say how fast the whole model is moving:
  one per joint, plus six for a root that is free to move."
  [model]
  (+ (dof model) (if (base model) 6 0)))

(defn- split
  "A generalised velocity back into the base's six and the joints' rest."
  [model u]
  (if (base model)
    [(vec (take 6 u)) (vec (drop 6 u))]
    [(vec (repeat 6 0.0)) (vec u)]))

(defn frame-of
  "Where body `i` is, in world terms. `-1` is the free root itself.

  Given `frames` -- an already-computed `poses` -- it reads from those
  instead of walking the tree again, which is what a caller asking about
  several bodies at the same configuration should do."
  ([model q-vec root-state i] (frame-of model q-vec root-state i nil))
  ([model q-vec root-state i frames]
   (if (neg? (long i))
     {:rot (or (:rot root-state) q/identity-q) :pos (or (:pos root-state) v/zero)}
     (nth (or frames (poses model q-vec root-state)) i))))

(defn- velocity-at
  "The world velocity of the point of body `i` at `p`, given every body's
  spatial velocity and where they are."
  [vels ^doubles v0 {:keys [rot pos]} i p]
  (let [^doubles vi (if (neg? (long i)) v0 (nth vels i))
        w (q/rotate rot [(aget vi 0) (aget vi 1) (aget vi 2)])
        at-origin (q/rotate rot [(aget vi 3) (aget vi 4) (aget vi 5)])]
    (v/add at-origin (v/cross w (v/sub p pos)))))

(defn point-velocity
  "How fast the point of body `i` that is at world position `p` is moving.
  `-1` is the free root itself.

  A point *fixed to the body*, not a point in space: the body's own
  material carries it, so it moves with its spin as well as its travel."
  [model q-vec state i p]
  (let [ls (links model q-vec)
        root-state (:base state)
        v0 (a/f64 (or (:vel root-state) (repeat 6 0.0)))
        vels (spatial-velocities ls v0 (:qd state))]
    (velocity-at vels v0 (frame-of model q-vec root-state i) i p)))

(defn- jacobian*
  "The Jacobian, given the per-configuration data already computed."
  [model ls frame i p]
  (let [n (generalised-dof model)
        column (fn [j]
                 (let [u (assoc (vec (repeat n 0.0)) j 1.0)
                       [v0 qd] (split model u)
                       v0 (a/f64 v0)]
                   (velocity-at (spatial-velocities ls v0 qd) v0 frame i p)))]
    (lin/transpose (mapv column (range n)))))

(defn point-jacobian
  "The 3 by `generalised-dof` matrix taking a generalised velocity to the
  world velocity of the point of link `i` at `p`.

  Built a column at a time, by asking what one unit of each generalised
  velocity on its own does. That is not an approximation -- the map is
  linear, which is the whole reason a Jacobian exists -- and it costs a
  velocity recursion per column where a purpose-built one would walk the
  path from the root once. At a ragdoll's twenty-odd degrees of freedom
  the difference is not worth the second implementation to get wrong."
  ([model q-vec root-state i p]
   (jacobian* model (links model q-vec) (frame-of model q-vec root-state i) i p))
  ([model q-vec root-state i p ls frames]
   (jacobian* model ls (frame-of model q-vec root-state i frames) i p)))

(defn inverse-mass-matrix
  "`H^-1`, the generalised inverse inertia, a column at a time out of the
  articulated body algorithm.

  Nothing is moving and there is no gravity, so the accelerations a unit
  generalised force produces *are* that column of the inverse. The
  inverse rather than the matrix itself because everything asked of it
  here -- how hard is this point to push, what does an impulse do -- is
  a question about the inverse, and forming H only to factor it again
  would be work in both directions."
  ([model q-vec] (inverse-mass-matrix model q-vec (links model q-vec)))
  ([model q-vec ls]
   (let [n (generalised-dof model)
         root (base model)
         nj (dof model)
         zero-q (vec (repeat nj 0.0))
         free {:gravity [0.0 0.0 0.0] :ls ls}
         column (fn [j]
                  (let [{:keys [qdd base-acc]}
                        (if (and root (< j 6))
                          (forward-dynamics model q-vec zero-q zero-q
                                            (assoc free
                                                   :base-force (assoc (vec (repeat 6 0.0)) j 1.0)
                                                   :base {:vel (vec (repeat 6 0.0))}))
                          (forward-dynamics model q-vec zero-q
                                            (assoc zero-q (- j (if root 6 0)) 1.0)
                                            (cond-> free
                                              root (assoc :base {:vel (vec (repeat 6 0.0))}))))]
                    (if root (vec (concat base-acc qdd)) (vec qdd))))]
     (lin/transpose (mapv column (range n))))))

(defn- point-force
  "A unit force at world point `p` along `dir`, as a spatial force in the
  coordinates of the body whose `frame` is given."
  ^doubles [{:keys [rot pos]} p dir]
  (let [inv (q/conjugate rot)
        f (q/rotate inv dir)
        r (q/rotate inv (v/sub p pos))]
    (a/f64 (concat (v/cross r f) f))))

(defn- generalised-force
  "The generalised force a spatial force `f` on body `i` makes.

  A walk from the body it acts on up to the root, taking each joint's
  share along its own axis and carrying the rest to the parent. Bodies
  off that path feel nothing at all, which is why this is a walk and not
  a matrix -- a hand pushed sideways says nothing about the other arm."
  ^doubles [model ls i f]
  (let [root? (some? (base model))
        nj (total-dof ls)
        out (a/f64 (+ nj (if root? 6 0)))
        base-off (if root? 6 0)]
    (loop [j (long i) ^doubles f f]
      (if (neg? j)
        (when root? (dotimes [k 6] (aset out k (aget f k))))
        (let [^Link l (nth ls j)
              ^doubles share (s-dot (.-s l) (.-ndof l) f)]
          (dotimes [k (.-ndof l)]
            (aset out (+ base-off (.-offset l) k) (aget share k)))
          (recur (.-parent l) (mv6 (.-xt l) f)))))
    out))

(defn- ancestry
  "Body `i` and every parent above it, tips first. `-1` is the root and
  has no entry -- it is where the walk ends."
  [ls i]
  (loop [j (long i) acc []]
    (if (neg? j)
      acc
      (recur (let [^Link l (nth ls j)] (.-parent l)) (conj acc j)))))

(defn- delta-from
  "The change in generalised velocity an impulse makes, whatever kind.

  The articulated body algorithm with the velocity terms gone. An
  impulse is a force with no duration, so there is no time for a
  velocity product to contribute anything and the bias forces are the
  impulse itself; the inertias are the ones already built for this
  configuration. What is left is an inward walk and an outward sweep,
  both linear in the links.

  `pa` maps a body to the spatial impulse sitting on it and `pa0` is the
  root's, which is how a contact enters. `extra` maps a joint to a
  generalised impulse applied along its own axes, which is how a limit
  enters. Both end up in the same `uu` and neither needs the recursion
  to know which it was.

  This is what replaced forming the inverse inertia matrix. That took
  `n + 6` runs of the full algorithm to build and was, measured, the
  entire cost of a contact -- 33ms of a 34ms solve on a sixteen link
  model."
  ^doubles [model ls ai path pa pa0 extra]
  (let [n (count ls)
        root (base model)
        inward (reduce
                (fn [{:keys [pa pa0 uu] :as acc} j]
                  (let [^Link l (nth ls j)
                        nd (.-ndof l)
                        ^doubles s (.-s l)
                        ^Node node (aget ^objects (:links ai) j)
                        ^doubles u (.-u node)
                        dinv (.-dinv node)
                        ^doubles paj (or (get pa j) (a/f64 6))
                        ^doubles seed (get extra j)
                        uj (let [^doubles share (s-dot s nd paj) o (a/f64 (max 1 nd))]
                             (dotimes [k nd]
                               (aset o k (- (if seed (aget seed k) 0.0) (aget share k))))
                             o)
                        acc (assoc acc :uu (assoc uu j uj))
                        up (mv6 (.-xt l) (cond-> (copy6 paj)
                                           dinv (s-apply-arr!
                                                 u nd (small-solve dinv nd uj))))]
                    (if (neg? (.-parent l))
                      (assoc acc :pa0 (add6! pa0 up))
                      (assoc acc :pa (update pa (.-parent l)
                                             #(add6! (or % (a/f64 6)) up))))))
                {:pa pa :pa0 pa0 :uu (vec (repeat n nil))}
                path)
        dv0 (if root
              (a/f64 (mapv - (root-solve (:ia0 ai) (:pa0 inward))))
              (a/f64 6))
        outward (reduce
                 (fn [{:keys [dv] :as acc} j]
                   (let [^Link l (nth ls j)
                         nd (.-ndof l)
                         ^doubles s (.-s l)
                         ^Node node (aget ^objects (:links ai) j)
                         ^doubles u (.-u node)
                         dinv (.-dinv node)
                         a' (mv6 (.-xup l)
                                 (if (neg? (.-parent l)) dv0 (nth dv (.-parent l))))
                         ;; Bodies off the path have no bias force on
                         ;; them at all, which is what the nil says;
                         ;; their joints still accelerate, because
                         ;; everything above them moved.
                         ^doubles uj (or (nth (:uu inward) j) (a/f64 (max 1 nd)))
                         dq (if dinv
                              (small-solve dinv nd
                                           (let [y (a/f64 (max 1 nd))]
                                             (dotimes [k nd]
                                               (let [b (* 6 k)]
                                                 (aset y k
                                                       (- (aget uj k)
                                                          (loop [t 0 acc2 0.0]
                                                            (if (= t 6)
                                                              acc2
                                                              (recur (inc t)
                                                                     (+ acc2 (* (aget u (+ b t))
                                                                                (aget a' t))))))))))
                                             y))
                              (a/f64 (max 1 nd)))]
                     (-> acc
                         (update :dv conj (s-apply-arr! (copy6 a') s nd dq))
                         (update :dq conj dq))))
                 {:dv [] :dq []}
                 (range n))
        base-off (if root 6 0)
        out (a/f64 (+ (total-dof ls) base-off))]
    (when root (dotimes [k 6] (aset out k (aget ^doubles dv0 k))))
    (dotimes [j n]
      (let [^Link l (nth ls j)
            ^doubles dq (nth (:dq outward) j)]
        (dotimes [k (.-ndof l)]
          (aset out (+ base-off (.-offset l) k) (aget dq k)))))
    out))

(defn- impulse-delta
  "The change in generalised velocity from a spatial impulse `f` on body
  `i`, in that body's own coordinates. `-1` is the root.

  The articulated body algorithm again, with the velocity terms gone.
  An impulse is a force with no duration, so there is no time for a
  velocity product to contribute anything and the bias forces are the
  impulse itself; the inertias are the ones already built for this
  configuration. What is left is an inward pass and an outward one, both
  linear in the links.

  This is what replaced forming the inverse inertia matrix. That took
  `n + 6` runs of the full algorithm to build and was, measured, the
  entire cost of a contact -- 33ms of a 34ms solve on a sixteen link
  model. A contact asks about three directions, so three of these do
  instead."
  ^doubles [model ls ai i ^doubles f]
  (let [neg-f (let [o (a/f64 6)] (dotimes [k 6] (aset o k (- (aget f k)))) o)]
    (delta-from model ls ai
                (ancestry ls i)
                (if (neg? (long i)) {} {i neg-f})
                (if (neg? (long i)) (copy6 neg-f) (a/f64 6))
                {})))

(defn- joint-delta
  "The change in generalised velocity from a generalised impulse `w` on
  joint `j`'s own axes.

  What a joint limit needs, where a contact needs `impulse-delta`. The
  walk is the same one -- a limit is a push between a link and its own
  parent rather than between a link and the floor, and neither of them
  is anything the recursion has to be told about."
  ^doubles [model ls ai j ^doubles w]
  (delta-from model ls ai (ancestry ls j) {} (a/f64 6) {j w}))

(defn- pair-ancestry
  "Everything between either of two bodies and the root, children first.

  Descending index order is that order, since a link's parent always has
  a smaller index than the link. Each body appears once even where the
  two walks meet, which they do at the first common ancestor -- and
  meeting is the point: an arm hitting a thigh is felt at the pelvis,
  and the recursion has to add both shares to it before passing it on."
  [ls i j]
  (vec (sort > (distinct (concat (ancestry ls i) (ancestry ls j))))))

(defn- pair-delta
  "The change in generalised velocity from an impulse between two bodies
  of the same model -- `f` on body `i` and minus `f'` on body `j`.

  The two spatial forces are given separately because each is in its own
  body's coordinates: the same push, written twice, because a spatial
  vector only means anything alongside the frame it is in."
  ^doubles [model ls ai i ^doubles f j ^doubles f']
  (let [neg (fn [^doubles x] (let [o (a/f64 6)] (dotimes [k 6] (aset o k (- (aget x k)))) o))
        pa (cond-> {}
             (not (neg? (long i))) (assoc i (neg f))
             (not (neg? (long j))) (assoc j (copy6 f')))
        pa0 (cond-> (a/f64 6)
              (neg? (long i)) (add6! (neg f))
              (neg? (long j)) (add6! f'))]
    (delta-from model ls ai (pair-ancestry ls i j) pa pa0 {})))

(defn- pair-force
  "The generalised force an impulse between two bodies makes: what it
  does to `i` less what it does to `j`.

  The difference, because what a contact constrains is the *relative*
  velocity of the two surfaces, and a push that moves both bodies the
  same way does not change that at all."
  ^doubles [model ls i ^doubles f j ^doubles f']
  (let [ga (generalised-force model ls i f)
        gb (generalised-force model ls j f')
        n (alength ga)
        out (a/f64 n)]
    (dotimes [k n] (aset out k (- (aget ga k) (aget gb k))))
    out))

(defn- response
  "How a unit impulse at world point `p` on body `i` along `dir` is felt.

  `:g` is the generalised force it makes, so the closing speed is
  `g . u`; `:delta` is what it does to `u`; `:m` is the mass felt there.
  Zero rather than infinity when nothing can move that way -- this is
  the number an impulse gets multiplied by, and a direction that cannot
  give is a direction no impulse is worth applying."
  [model ls ai frame i p dir]
  (let [f (point-force frame p dir)
        g (generalised-force model ls i f)
        delta (impulse-delta model ls ai i f)
        w (dot-n g delta)]
    {:g g :delta delta :m (if (> w 1e-12) (/ 1.0 w) 0.0)}))

(defn- pair-response
  "How a unit impulse between two bodies of the same model is felt.

  The same three numbers as `response`, for a push on `i` and an equal
  and opposite one on `j` at the same world point.

  `m` comes out *smaller* than either body's own effective mass there,
  which is the reverse of what it sounds like it should be and is
  right. It is a reduced mass. The constraint is on how fast the two
  surfaces approach each other, and an impulse between two things that
  can both move changes that faster than the same impulse against
  something that cannot move at all -- both ends give."
  [model ls ai frame-i i frame-j j p dir]
  (let [f (point-force frame-i p dir)
        f' (point-force frame-j p dir)
        g (pair-force model ls i f j f')
        delta (pair-delta model ls ai i f j f')
        w (dot-n g delta)]
    {:g g :delta delta :m (if (> w 1e-12) (/ 1.0 w) 0.0)}))

(defn configuration
  "Everything about where a model is that a step's worth of questions
  can share: each link's transform and inertia, what every joint weighs
  from above, and where every part has ended up.

  All of it depends only on `q` and the root's pose, so it is built once
  and handed to `response-at` as many times as the contacts need. Opaque
  -- the shapes inside are this namespace's business."
  [model q-vec root-state]
  (let [ls (links model q-vec)]
    {:ls ls
     :ai (articulated-inertias model ls)
     :frames (poses model q-vec root-state)
     :q q-vec
     :root root-state}))

(defn response-at
  "How a unit impulse at world point `p` on body `i` along `dir` is felt,
  given a `configuration`.

  `{:force :delta :mass}` -- the generalised force the impulse makes, so
  that the closing speed along `dir` is `force . u`; what it does to
  `u`; and the mass it meets there. Both vectors are flat arrays as long
  as `generalised-dof`, which is raw for a public interface and is the
  point: this is what a contact solver holds per contact per direction
  and iterates over, and boxing it would undo the reason it is fast.

  `mass` is zero rather than infinite where nothing can move that way --
  it is the number an impulse gets multiplied by."
  [model cfg i p dir]
  (let [{:keys [ls ai frames q root]} cfg]
    (response model ls ai (frame-of model q root i frames) i p (v/normalize dir))))

(defn pair-response-at
  "The same, for an impulse between two bodies of the *same* model --
  `i` pushed along `dir` and `j` the other way.

  Not the same as asking twice. The two responses interact through
  everything the pair have in common, which is at least the root, and
  the mass it meets is a reduced mass rather than either body's own."
  [model cfg i j p dir]
  (let [{:keys [ls ai frames q root]} cfg]
    (pair-response model ls ai
                   (frame-of model q root i frames) i
                   (frame-of model q root j frames) j
                   p (v/normalize dir))))

(defn impulse-at
  "What a unit impulse at world point `p` on body `i`, along `dir`, does.
  `-1` is the free root itself.

  Returns `{:delta-u :effective-mass}`: the change in generalised
  velocity per unit of impulse, and the mass the impulse feels there.

  The second is `1 / (d^T J H^-1 J^T d)`, and it is the number a contact
  solver actually wants -- how much velocity a given push buys at this
  point in this direction, with the whole articulated body hanging off
  it. A hand on the end of an outstretched arm is light; the same hand
  with the arm folded against the chest is most of a torso."
  [model q-vec root-state i p dir]
  (let [ls (links model q-vec)
        {:keys [delta m]} (response model ls (articulated-inertias model ls)
                                    (frame-of model q-vec root-state i)
                                    i p (v/normalize dir))]
    {:delta-u (vec (seq ^doubles delta))
     :effective-mass (if (pos? (double m)) m ##Inf)}))

(defn pair-impulse-at
  "What a unit impulse between two of the model's own bodies does --
  pushing `i` along `dir` at world point `p` and `j` the other way.

  Returns `{:delta-u :effective-mass}`, as `impulse-at` does. The mass
  is a reduced mass and comes out smaller than either body's alone:
  both ends of this push give, where a contact with the floor has one
  end that does not."
  [model q-vec root-state i j p dir]
  (let [ls (links model q-vec)
        frames (poses model q-vec root-state)
        {:keys [delta m]} (pair-response model ls (articulated-inertias model ls)
                                         (frame-of model q-vec root-state i frames) i
                                         (frame-of model q-vec root-state j frames) j
                                         p (v/normalize dir))]
    {:delta-u (vec (seq ^doubles delta))
     :effective-mass (if (pos? (double m)) m ##Inf)}))

(defn apply-impulse
  "`state` after an impulse of `magnitude` at world point `p` on body
  `i`, along `dir`. `-1` is the free root itself.

  Positions do not move -- an impulse is instantaneous by definition --
  so only the velocities change, and they change by `H^-1 J^T d` times
  the magnitude. This is the whole of what a contact solver needs of an
  articulated body, and it is why contacts can be added to one without
  the dynamics knowing anything about them."
  [model state i p dir magnitude]
  (let [magnitude (double magnitude)
        {:keys [delta-u]} (impulse-at model (:q state) (:base state) i p dir)
        root (base model)
        n (dof model)
        scaled-delta (mapv #(* (double %) magnitude) delta-u)]
    (if root
      (-> state
          (update-in [:base :vel] #(mapv + (vec (or % (repeat 6 0.0)))
                                         (take 6 scaled-delta)))
          (update :qd #(mapv + % (drop 6 scaled-delta))))
      (update state :qd #(mapv + % (take n scaled-delta))))))

(defn apply-pair-impulse
  "`state` after an impulse of `magnitude` between two of the model's own
  bodies: `i` pushed along `dir` at `p`, `j` pushed the other way.

  Nothing outside the model is touched, so its momentum is unchanged --
  which is the sharpest check there is that the two halves of the push
  really are equal and opposite."
  [model state i j p dir magnitude]
  (let [magnitude (double magnitude)
        {:keys [delta-u]} (pair-impulse-at model (:q state) (:base state) i j p dir)
        root (base model)
        n (dof model)
        scaled (mapv #(* (double %) magnitude) delta-u)]
    (if root
      (-> state
          (update-in [:base :vel] #(mapv + (vec (or % (repeat 6 0.0))) (take 6 scaled)))
          (update :qd #(mapv + % (drop 6 scaled))))
      (update state :qd #(mapv + % (take n scaled))))))

;; ---------------------------------------------------------------------------
;; Generalised velocity

(defn velocity
  "The model's generalised velocity: the root's spatial six, then one per
  joint. A bolted model has only the joints'."
  [model state]
  (if (base model)
    (vec (concat (or (:vel (:base state)) (repeat 6 0.0)) (:qd state)))
    (vec (:qd state))))

(defn with-velocity
  "`state` moving at the generalised velocity `u`."
  [model state u]
  (if (base model)
    (-> state
        (assoc :base (assoc (:base state) :vel (vec (take 6 u))))
        (assoc :qd (vec (drop 6 u))))
    (assoc state :qd (vec u))))

;; ---------------------------------------------------------------------------
;; Shapes, and what they run into

(defn- shape-pose
  "Where a part's collision shape sits in that part's own frame.

  Centred on the centre of mass unless told otherwise, which is right
  for a limb and saves every model repeating it. A link's frame is at
  its *joint*, not in the middle of it, so a shape left at the origin
  would stick out of the elbow."
  [part]
  (or (:shape-pose part) {:rot q/identity-q :pos (:com part)}))

(defn- part-body
  "A part's collision shape as an `allgo.physics.rigid` body, placed
  where the part is now, or nil if it was never given one."
  [part frame]
  (when-let [kind (:shape part)]
    (let [{srot :rot spos :pos} (shape-pose part)
          rot (q/mul (:rot frame) (or srot q/identity-q))
          pos (v/add (:pos frame) (q/rotate (:rot frame) (or spos v/zero)))
          common {:pos pos :rot rot :density 1.0}]
      (case kind
        :box (rigid/box (assoc common :size (:size part)))
        :ball (rigid/ball (assoc common :radius (:radius part)))))))

(defn collision-bodies
  "Every shaped part of the model, as `[{:link i :body b} ...]`.

  `-1` is the free root. Parts with no `:shape` are not here: a link
  that is only there to carry a degree of freedom -- the middle of three
  stacked hinges standing in for a shoulder -- has no geometry and
  should collide with nothing."
  ([model q-vec root-state]
   (collision-bodies model q-vec root-state (poses model q-vec root-state)))
  ([model q-vec root-state frames]
   (let [parts (chain model)
         root (base model)]
     (into (if-let [b (and root (part-body root (frame-of model q-vec root-state -1)))]
             [{:link -1 :body b}]
             [])
           (keep (fn [i]
                   (when-let [b (part-body (nth parts i) (nth frames i))]
                     {:link i :body b})))
           (range (count parts))))))

(defn contacts-with
  "Every contact between the model's shapes and the static `obstacles`.

  `obstacles` are ordinary `allgo.physics.rigid` bodies -- a floor, a
  ramp, scenery. Nothing is applied back to them, which is what makes
  them static and what makes this the easy half: one side of every
  contact has no degrees of freedom to account for.

  The normal points the way the link has to be pushed, which is the
  opposite of what `allgo.physics.contact` reports for the pair. Stating
  it in the direction the solver will use it saves a negation at every
  later step and a sign error at one of them."
  ([model q-vec root-state obstacles]
   (contacts-with model q-vec root-state obstacles (poses model q-vec root-state)))
  ([model q-vec root-state obstacles frames]
   (into []
         (for [{:keys [link body]} (collision-bodies model q-vec root-state frames)
               ob obstacles
               c (contact/between 0 1 body ob)]
           {:link link
            :point (:point c)
            :normal (v/negate (:normal c))
            :depth (:depth c)}))))

;; ---------------------------------------------------------------------------
;; Solving them

(def default-contact
  "Settings for `solve-contacts`, matching `allgo.physics.solver`'s where
  they mean the same thing."
  {:iterations 8
   :friction 0.6
   :restitution 0.0
   :slop 0.005
   :bias-factor 0.2
   :max-push-speed 3.0})

(defn- tangents
  "Two unit directions across `n`, any two."
  [n]
  (let [a (if (< (abs (double (nth n 0))) 0.9) [1.0 0.0 0.0] [0.0 1.0 0.0])
        t1 (v/normalize (v/cross n a))]
    [t1 (v/cross n t1)]))

(defn limited?
  "Whether any joint has a limit to be checked at all."
  [model]
  (boolean (some #(or (:limit %) (:cone %) (:twist %)) (chain model))))

(defn- world-extent
  "Half the width of a body's world-axis box, per axis.

  Only good enough to reject a pair before the exact test is asked, and
  that is all it is for: the separating axis test costs fifteen axes of
  dot and cross products and most pairs of a body's own limbs are
  nowhere near each other."
  [b]
  (if (= :ball (:shape b))
    (let [r (double (:radius b))] [r r r])
    (let [[hx hy hz] (mapv #(* 0.5 (double %)) (:size b))
          rot (:rot b)
          ax (q/rotate rot [1.0 0.0 0.0])
          ay (q/rotate rot [0.0 1.0 0.0])
          az (q/rotate rot [0.0 0.0 1.0])]
      (mapv (fn [k] (+ (* hx (abs (double (nth ax k))))
                       (* hy (abs (double (nth ay k))))
                       (* hz (abs (double (nth az k))))))
            (range 3)))))

(defn- apart?
  "Whether two bodies' world-axis boxes miss each other."
  [a b]
  (let [ea (world-extent a) eb (world-extent b)]
    (boolean (some (fn [k]
                     (> (abs (- (double (nth (:pos a) k)) (double (nth (:pos b) k))))
                        (+ (double (nth ea k)) (double (nth eb k)))))
                   (range 3)))))

(defn- collidable?
  "Whether two of a model's own parts are allowed to touch.

  A link and its parent are not: they meet at the joint and overlap
  there by construction, so a contact between them is a permanent one
  pushing the body apart. Anything else is fair game unless a link's
  `:no-collide` says otherwise, which is how a shoulder that sits inside
  the ribcage is told to stay there."
  [parts i j]
  (let [i (long i) j (long j)
        parent-of (fn [k] (if (neg? (long k)) -2 (long (:parent (nth parts k)))))
        excluded (fn [k other]
                   (and (not (neg? (long k)))
                        (contains? (set (:no-collide (nth parts k))) other)))]
    (and (not= i j)
         (not= (parent-of i) j)
         (not= (parent-of j) i)
         (not (excluded i j))
         (not (excluded j i)))))

(defn self-contacts
  "Every contact between one of the model's shapes and another of them.

  Without this a forearm passes through a thigh, which on a ragdoll is
  the most visible thing left wrong. It is also the expensive thing: a
  model's parts are all awake and all moving, so there is no sleeping to
  lean on the way a wall of bricks does, and the pairs go as the square
  of the parts. The world-axis boxes reject most of them before the
  exact test is asked.

  Normals point the way the *first* link of each pair has to be pushed."
  [model q-vec root-state frames]
  (let [parts (chain model)
        bodies (collision-bodies model q-vec root-state frames)]
    (into []
          (for [[x y] (map vector (range) bodies)
                [x' y'] (map vector (range) bodies)
                :when (< (long x) (long x'))
                :let [{la :link ba :body} y
                      {lb :link bb :body} y']
                :when (collidable? parts la lb)
                :when (not (apart? ba bb))
                c (contact/between 0 1 ba bb)]
            {:link la :other lb
             :point (:point c)
             :normal (v/negate (:normal c))
             :depth (:depth c)}))))

(defn- limit-rows
  "The joint limits currently being pushed against, as one-sided
  constraints on the generalised velocity.

  A hinge with `:limit [lo hi]` is the easy half: the violation is a
  number and the direction to push is its own axis.

  A ball joint with `:cone theta` is a limit on how far the bone may
  swing from where it points at rest -- `:limit-axis` in the child's own
  frame, the direction of the centre of mass unless said otherwise. Turn
  it by the joint's rotation and the angle to the rest direction is the
  swing; the axis to turn about to reduce it is the cross product of the
  two, carried back into the child's frame because that is where a
  spherical joint's velocity lives. There is no limit on twist: a
  ragdoll's shoulder needs a cone and does not care.

  Only violated limits are returned. A joint inside its range is not a
  constraint and solving it as one costs a row for nothing."
  [model ls ai q dt bias-factor max-push]
  (let [dt (double dt) bias-factor (double bias-factor) max-push (double max-push)
        parts (chain model)
        root? (some? (base model))
        base-off (if root? 6 0)
        nd-total (+ base-off (total-dof ls))]
    (into []
          (keep
           (fn [j]
             (let [link (nth parts j)
                   ^Link l (nth ls j)
                   nd (.-ndof l)
                   offset (.-offset l)
                   x (nth q j nil)
                   ;; `depth` is how far past the limit, `w` the
                   ;; generalised impulse direction that comes back.
                   [depth ^doubles w]
                   (case (kind-of link)
                     (:revolute :prismatic)
                     (when-let [[lo hi] (:limit link)]
                       (let [xv (double x)]
                         (cond (< xv (double lo)) [(- (double lo) xv) (a/f64 [1.0])]
                               (> xv (double hi)) [(- xv (double hi)) (a/f64 [-1.0])]
                               :else nil)))
                     :spherical
                     (let [rest-dir (v/normalize (or (:limit-axis link) (:com link)))
                           rot (if (and (sequential? x) (= 4 (count x)))
                                 ;; Both halves of the sphere name the
                                 ;; same rotation; the twist angle read
                                 ;; off the wrong one is out by a turn.
                                 (let [r (vec x)] (if (neg? (double (nth r 3))) (mapv - r) r))
                                 q/identity-q)
                           bone (q/rotate rot rest-dir)
                           swing (Math/acos (max -1.0 (min 1.0 (v/dot bone rest-dir))))
                           cone (:cone link)
                           twist-max (:twist link)
                           ;; How far it has turned about the bone
                           ;; itself, which the cone says nothing about:
                           ;; a head can be within forty degrees of
                           ;; upright and still be facing backwards.
                           proj (v/dot [(nth rot 0) (nth rot 1) (nth rot 2)] rest-dir)
                           twist (* 2.0 (Math/atan2 proj (double (nth rot 3))))]
                       (cond
                         (and cone (> swing (double cone)))
                         (let [n (v/cross bone rest-dir)
                               len (v/length n)]
                           (when (> len 1e-9)
                             (let [axis (q/rotate (q/conjugate rot) (v/scale n (/ 1.0 len)))]
                               [(- swing (double cone))
                                (a/f64 [(nth axis 0) (nth axis 1) (nth axis 2)])])))

                         (and twist-max (> (abs twist) (double twist-max)))
                         (let [sgn (if (pos? twist) -1.0 1.0)]
                           [(- (abs twist) (double twist-max))
                            (a/f64 [(* sgn (nth rest-dir 0))
                                    (* sgn (nth rest-dir 1))
                                    (* sgn (nth rest-dir 2))])])

                         :else nil))
                     nil)]
               (when depth
                 (let [delta (joint-delta model ls ai j w)
                       g (let [o (a/f64 nd-total)]
                           (dotimes [k nd]
                             (aset o (+ base-off (long offset) k) (aget w k)))
                           o)
                       wgt (dot-n g delta)]
                   {:g g :delta delta
                    :m (if (> wgt 1e-12) (/ 1.0 wgt) 0.0)
                    :bias (min max-push (/ (* bias-factor (double depth)) dt))})))))
          (range (count parts)))))

(defn solve-constraints
  "`state` with its velocities corrected so that neither the contacts nor
  the joint limits are being driven into.

  Sequential impulse, the same as `allgo.physics.solver` runs, over the
  same kind of accumulated clamped impulses -- and the reason it can be
  the same is the point of the previous section. A contact solver only
  ever asks a body how much velocity a push buys and then pushes; it
  does not care that the answer came through a chain of joints.

  What is different is where the work goes. `H^-1` depends only on where
  the joints are, so it is built once for the step and not once per
  contact per iteration; each contact's three directions turn into a
  generalised force `g = J^T d` and a response `H^-1 g` once, and after
  that an iteration is dot products. Otherwise a ragdoll would spend its
  frame rebuilding the same matrix eighty times."
  [model state cs dt opts]
  (let [{:keys [iterations friction restitution slop bias-factor max-push-speed]}
        (merge default-contact opts)
        dt (double dt)
        q (:q state)
        root-state (:base state)
        ;; All of this depends only on where the joints are, so it is
        ;; built once for the step. Rebuilt per contact per iteration it
        ;; was most of the frame.
        ls (or (:ls opts) (links model q))
        frames (or (:frames opts) (poses model q root-state))
        ai (or (:ai opts) (articulated-inertias model ls))
        u0 (a/f64 (velocity model state))
        prep (mapv (fn [c]
                     (let [i (:link c)
                           other (:other c)
                           frame (frame-of model q root-state i frames)
                           n (v/normalize (:normal c))
                           [t1 t2] (tangents n)
                           ;; A contact against the world pushes one
                           ;; body; a contact between two of the model's
                           ;; own parts pushes both, and what it
                           ;; constrains is their relative velocity.
                           along (if other
                                   (let [frame' (frame-of model q root-state other frames)]
                                     #(pair-response model ls ai frame i frame' other
                                                     (:point c) %))
                                   #(response model ls ai frame i (:point c) %))
                           dirs {:n (along n) :t1 (along t1) :t2 (along t2)}]
                       (assoc dirs
                              ;; The closing speed as the step began.
                              ;; Restitution is measured against this and
                              ;; nothing later, for the same reason the
                              ;; rigid solvers record it once.
                              :approach (dot-n (:g (:n dirs)) u0)
                              :bias (min (double max-push-speed)
                                         (/ (* (double bias-factor)
                                               (max 0.0 (- (double (:depth c))
                                                           (double slop))))
                                            dt)))))
                   cs)
        limits (limit-rows model ls ai q dt bias-factor max-push-speed)
        k (count prep)
        nl (count limits)]
    (if (and (zero? k) (zero? nl))
      state
      ;; One mutable generalised velocity for the whole sweep, and one
      ;; flat array of accumulated impulses beside it. Every push is a
      ;; scaled add into the first and every closing speed a dot product
      ;; out of it, so the inner loop allocates nothing: eight
      ;; iterations over fifty contacts is twelve hundred of each, and
      ;; on boxed vectors that was most of what a contact cost.
      (let [^doubles u u0
            ^doubles acc (a/f64 (* 3 k))
            ^doubles lacc (a/f64 (max 1 nl))]
        (dotimes [_ (long iterations)]
          ;; Limits first. A joint being held inside its range changes
          ;; what the contacts below it are pushing against, and a knee
          ;; that has folded backwards is a worse thing to look at than
          ;; a foot a millimetre into the floor.
          (dotimes [i nl]
            (let [{:keys [m bias]} (nth limits i)
                  ^doubles gl (:g (nth limits i))
                  ^doubles dl (:delta (nth limits i))
                  vn (dot-n gl u)
                  old (aget lacc i)
                  nw (max 0.0 (+ old (* (double m) (- (double bias) vn))))]
              (axpy-n! u dl (- nw old))
              (aset lacc i nw)))
          (dotimes [i k]
            (let [{:keys [n t1 t2 approach bias]} (nth prep i)
                  b (* 3 i)
                  ;; The normal first: friction is bounded by the force
                  ;; it rides on, so it wants this sweep's answer and
                  ;; not the last one's.
                  target (max (double bias)
                              (if (< (double approach) -0.5)
                                (* (- (double restitution)) (double approach))
                                0.0))
                  ^doubles gn (:g n)
                  ^doubles dn (:delta n)
                  vn (dot-n gn u)
                  an (aget acc b)
                  an' (max 0.0 (+ an (* (double (:m n)) (- target vn))))]
              (axpy-n! u dn (- an' an))
              (aset acc b an')
              (let [limit (* (double friction) an')]
                (dotimes [t 2]
                  (let [dir (if (zero? t) t1 t2)
                        ^doubles gt (:g dir)
                        ^doubles dt' (:delta dir)
                        idx (+ b 1 t)
                        vt (dot-n gt u)
                        old (aget acc idx)
                        a' (min limit (max (- limit)
                                           (+ old (* (double (:m dir)) (- vt)))))]
                    (axpy-n! u dt' (- a' old))
                    (aset acc idx a')))))))
        (with-velocity model state (vec (seq u)))))))
