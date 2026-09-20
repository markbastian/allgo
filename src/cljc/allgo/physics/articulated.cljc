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

  Still more than it should. A chain lying on the floor, JVM, per step:

      links         0      2      4      8     16
      contacts      4      6      8     16     52
      before     0.71   4.92  12.19  29.85  154.1
      now        0.67   2.03   4.24   6.87   31.35

  Five times better and not yet enough. What was fixed was structural: a
  contact cost a whole inverse inertia matrix, which cost `n + 6` runs
  of the articulated body algorithm, and it is now three O(n) impulse
  responses. The per-configuration data -- link transforms, poses,
  articulated inertias -- is built once a step instead of once per
  matrix column, and the inward walk visits only the path from the body
  that was hit to the root.

  What is left is not structural. Every 6x6 here is a vector of vectors
  and every spatial vector a vector of six boxed doubles, so a step of a
  settled sixteen-link chain allocates on the order of fifty thousand of
  them. That is the same treatment `allgo.physics.solver` gave its
  contact solve -- primitive doubles, no allocation in the inner loop --
  and the same one `allgo.physics.contact` is still waiting for.

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
  (:require [allgo.geometry.quaternion :as q]
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

(def ^:private zero6 [0.0 0.0 0.0 0.0 0.0 0.0])

(defn- s-times
  "`S x`, the spatial vector a joint's rates make, reading them out of
  the flat velocity `x` at this joint's offset."
  [s x ^long offset]
  (reduce (fn [acc k]
            (let [r (double (nth x (+ offset k)))
                  col (nth s k)]
              (if (zero? r) acc (mapv (fn [a c] (+ (double a) (* r (double c)))) acc col))))
          zero6
          (range (count s))))

(defn- st-times
  "`S^T f`, each of a joint's axes' share of a spatial force."
  [s f]
  (mapv (fn [col] (reduce + (map * col f))) s))

(defn- d-inverse
  "The inverse of `S^T IA S`, or nil where the joint cannot move.

  Nil rather than a huge number: a joint whose subtree presents no
  inertia along one of its axes has no acceleration there to compute,
  and inventing one is how a zero-mass link tears a model apart."
  [d]
  (case (count d)
    0 []
    1 (let [x (double (get-in d [0 0]))]
        (when (> (abs x) 1e-12) [[(/ 1.0 x)]]))
    (lin/inverse d)))

(defn- outer
  "The 6x6 `a b^T`."
  [a b]
  (mapv (fn [ai] (mapv (fn [bj] (* (double ai) (double bj))) b)) a))

(defn- u-dinv-ut
  "`U D^-1 U^T`, the part of an articulated inertia the joint gives away."
  [u dinv]
  (let [n (count u)]
    (reduce (fn [m [k l]]
              (let [w (double (get-in dinv [k l]))]
                (if (zero? w)
                  m
                  (lin/mat-add m (lin/mat-scale (outer (nth u k) (nth u l)) w)))))
            (vec (repeat 6 zero6))
            (for [k (range n) l (range n)] [k l]))))

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
              xup (lin/mat-mul (transform (joint-transform link (nth q-vec i nil)))
                               (transform (:origin link)))]
          [(conj acc {:xup xup
                      ;; Kept rather than taken again: every inward
                      ;; pass, of which there are several a step, wants
                      ;; it.
                      :xt (lin/transpose xup)
                      :s (subspace link)
                      :ndof nd
                      ;; Where this joint's rates start in the flat
                      ;; velocity vector. Joints are no longer all one
                      ;; number wide, so nothing can index by link any
                      ;; more.
                      :offset off
                      :i (spatial-inertia (:mass link) (:com link) (:inertia link))
                      :parent (long (:parent link))})
           (+ off nd)]))
      [[] 0]
      (range (count parts))))))

(defn- scaled [s ^double x] (mapv #(* (double %) x) s))

(defn- spatial-velocities
  "Every link's spatial velocity, in its own coordinates.

  Linear in `(v0, qd)` and used three ways because of it: to say what a
  model is doing now, and -- one unit of generalised velocity at a time
  -- to build the Jacobian of any point on it."
  [ls v0 qd]
  (reduce (fn [acc i]
            (let [{:keys [xup s parent offset]} (nth ls i)
                  vp (if (neg? (long parent)) v0 (nth acc parent))]
              (conj acc (mapv + (lin/mat-vec xup vp) (s-times s qd offset)))))
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
                (let [{:keys [xup s parent offset]} (nth ls i)
                      inertia (:i (nth ls i))
                      vj (s-times s qd offset)
                      vp (if (neg? parent) zero6 (:v (nth acc parent)))
                      ap (if (neg? parent) (base-acceleration gravity) (:a (nth acc parent)))
                      vi (mapv + (lin/mat-vec xup vp) vj)
                      ai (mapv + (lin/mat-vec xup ap)
                               (s-times s qdd offset)
                               (lin/mat-vec (crm vi) vj))
                      fi (mapv + (lin/mat-vec inertia ai)
                               (lin/mat-vec (crf vi) (lin/mat-vec inertia vi)))]
                  (conj acc {:v vi :a ai :f fi})))
              []
              (range n))
         ;; Inward: each link's force, plus everything its children left
         ;; on it, read off along the joint axis.
         forces (reduce (fn [fs i]
                          (let [{:keys [xup parent]} (nth ls i)
                                fi (nth fs i)]
                            (if (neg? parent)
                              fs
                              (update fs parent
                                      #(mapv + % (lin/mat-vec (lin/transpose xup) fi))))))
                        (mapv :f out)
                        (reverse (range n)))]
     ;; Each joint's share of the force on its own link, laid back into
     ;; the flat vector the caller handed its rates in.
     (reduce (fn [acc i]
               (let [{:keys [s offset]} (nth ls i)]
                 (reduce (fn [a k] (assoc a (+ (long offset) k) (nth (st-times s (nth forces i)) k)))
                         acc
                         (range (count s)))))
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
        root (base model)]
    (reduce
     (fn [{:keys [links ia0]} i]
       (let [{:keys [xup xt parent s]} (nth ls i)
             ia (:ia (nth links i))
             ;; One column of U per degree of freedom, and D their
             ;; projection onto each other -- a number for a hinge, a
             ;; three by three for a ball joint.
             u (mapv #(lin/mat-vec ia %) s)
             d (mapv (fn [sk] (mapv (fn [ul] (reduce + (map * sk ul))) u)) s)
             dinv (d-inverse d)
             ia' (if dinv (lin/mat-sub ia (u-dinv-ut u dinv)) ia)
             links (-> links
                       (assoc-in [i :u] u)
                       (assoc-in [i :dinv] dinv)
                       (assoc-in [i :ia-free] ia'))
             up (lin/mat-mul xt (lin/mat-mul ia' xup))]
         (if (neg? (long parent))
           {:links links :ia0 (when ia0 (lin/mat-add ia0 up))}
           {:links (update-in links [parent :ia] #(lin/mat-add % up)) :ia0 ia0})))
     {:links (mapv (fn [l] {:ia (:i l)}) ls)
      :ia0 (when root (spatial-inertia (:mass root) (:com root) (:inertia root)))}
     (reverse (range n)))))

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
                  [0.0 0.0 0.0 (double gx) (double gy) (double gz)])
         v0 (if root (vec (or (:vel root-state) (repeat 6 0.0))) (vec (repeat 6 0.0)))
         root-i (when root (spatial-inertia (:mass root) (:com root) (:inertia root)))
         ;; What every joint weighs from above, which depends only on
         ;; where the joints are. Handed in when the caller has already
         ;; built it for this configuration.
         ai (or given-ai (articulated-inertias model ls))
         ;; Pass one, outward: velocity, the acceleration that velocity
         ;; alone implies, and the force needed to hold the link on it.
         pass1 (reduce
                (fn [acc i]
                  (let [{:keys [xup s parent offset]} (nth ls i)
                        inertia (:i (nth ls i))
                        vj (s-times s qd offset)
                        vp (if (neg? parent) v0 (:v (nth acc parent)))
                        vi (mapv + (lin/mat-vec xup vp) vj)
                        ;; An outside force on a link reduces the bias
                        ;; force the link needs of its parent by exactly
                        ;; itself, which is the whole of how anything
                        ;; external gets in -- a contact, a thruster, a
                        ;; hand pushing.
                        f (when ext (nth ext i nil))]
                    (conj acc {:v vi
                               :c (lin/mat-vec (crm vi) vj)
                               :pa (let [bias (lin/mat-vec (crf vi)
                                                           (lin/mat-vec inertia vi))]
                                     (if f (mapv - bias f) bias))})))
                []
                (range n))
         ;; Pass two, inward: what the subtree beyond each joint looks
         ;; like to the link above it, and what the whole of it looks
         ;; like to the root.
         ;; Pass two, inward: the bias forces only. The inertias they
         ;; ride on were built once, above.
         pass2 (reduce
                (fn [{:keys [links pa0] :as acc} i]
                  (let [{:keys [s parent xt offset]} (nth ls i)
                        {:keys [pa c]} (nth links i)
                        {:keys [u dinv ia-free]} (nth (:links ai) i)
                        ;; What each of this joint's axes has left over
                        ;; once the bias force has taken its share.
                        uu (mapv (fn [k tk] (- (double tk) (double k)))
                                 (st-times s pa)
                                 (map #(nth tau (+ (long offset) %)) (range (count s))))
                        acc (assoc-in acc [:links i :uu] uu)
                        pa' (cond-> (mapv + pa (lin/mat-vec ia-free c))
                              dinv (as-> x (reduce (fn [a l]
                                                     (mapv + a (scaled (nth u l)
                                                                       (nth (lin/mat-vec dinv uu) l))))
                                                   x
                                                   (range (count u)))))
                        up-p (lin/mat-vec xt pa')]
                    (if (neg? (long parent))
                      (assoc acc :pa0 (when pa0 (mapv + pa0 up-p)))
                      (update-in acc [:links parent :pa] #(mapv + % up-p)))))
                {:links pass1
                 :pa0 (when root
                        (let [bias (lin/mat-vec (crf v0) (lin/mat-vec root-i v0))]
                          (if root-force (mapv - bias root-force) bias)))}
                (reverse (range n)))
         solved (:links pass2)
         ;; The root, if it is free: force-free in the falling frame.
         a0 (if root
              (let [m (:ia0 ai)
                    ;; Symmetric by construction and not quite by
                    ;; arithmetic, after a chain of congruences. Cholesky
                    ;; wants it to be, and averaging costs nothing.
                    sym (lin/mat-scale (lin/mat-add m (lin/transpose m)) 0.5)]
                (mapv - (lin/cholesky-solve sym (:pa0 pass2))))
              (mapv - a-grav))
         ;; Pass three, outward: each joint's acceleration, then the link's.
         out (reduce
              (fn [{:keys [a] :as acc} i]
                (let [{:keys [xup s parent]} (nth ls i)
                      {:keys [c uu]} (nth solved i)
                      {:keys [u dinv]} (nth (:links ai) i)
                      ap (if (neg? (long parent)) a0 (nth a parent))
                      a' (mapv + (lin/mat-vec xup ap) c)
                      qddi (if dinv
                             (lin/mat-vec dinv
                                          (mapv (fn [k ul] (- (double k) (reduce + (map * ul a'))))
                                                uu u))
                             (vec (repeat (count s) 0.0)))]
                  (-> acc
                      (update :a conj (mapv + a' (s-times s qddi 0)))
                      (update :qdd into qddi))))
              {:a [] :qdd []}
              (range n))]
     {:qdd (:qdd out)
      ;; Back out of the falling frame, which is where the root's own
      ;; share of gravity comes from.
      :base-acc (if root (mapv + a0 a-grav) (vec (repeat 6 0.0)))})))

;; ---------------------------------------------------------------------------
;; Moving it

(declare solve-contacts contacts-with velocity with-velocity poses)

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
        ls (links model (:q state))
        q' (mapv (fn [i]
                   (advance-coordinate (nth parts i) (nth (:q state) i nil)
                                       (:qd state) (:offset (nth ls i)) dt))
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

(defn step
  "One semi-implicit Euler step of `{:q :qd}`, and `:base` if there is
  one, under `tau`.

  Velocity first, then contacts, then position. The order is the whole
  of what makes a contact hold: the velocities are corrected before
  anything moves, so a body that was about to be driven into the floor
  never is, rather than being pulled back out afterwards.

  Contacts come either ready-made as `:contacts` or, more usually, from
  `:obstacles` -- static `allgo.physics.rigid` bodies this model is to
  be generated against. They are found at the configuration the step
  *starts* from, which is where the model actually is when the question
  is asked.

  The base's state is a pose and a spatial velocity in its own frame,
  not six more coordinates. Three reasons, and the third is the one that
  matters: a rotation has no three-number parameterisation without a
  singularity in it; a quaternion is the shape the rest of this library
  speaks; and a velocity in body coordinates is what the algorithms
  already produce, so nothing has to be converted on the way in."
  ([state model dt] (step state model dt nil))
  ([{:keys [q qd] :as state} model ^double dt {:keys [tau contacts obstacles] :as opts}]
   (let [root (base model)
         b (:base state)
         tau (or tau (vec (repeat (dof model) 0.0)))
         ;; Where everything is, computed once and handed to every pass
         ;; that needs it. It depends only on `q` and the root's pose,
         ;; neither of which changes until the very end of the step.
         ls (links model q)
         frames (poses model q b)
         {:keys [qdd base-acc]} (forward-dynamics model q qd tau
                                                  (cond-> (assoc opts :ls ls)
                                                    root (assoc :base b)))
         accel (if root (vec (concat base-acc qdd)) (vec qdd))
         moving (with-velocity model state
                  (mapv (fn [a c] (+ (double a) (* dt (double c))))
                        (velocity model state)
                        accel))
         cs (or contacts
                (when (seq obstacles) (contacts-with model q b obstacles frames)))]
     (advance-positions model
                        (if (seq cs)
                          (solve-contacts model moving cs dt
                                          (assoc opts :ls ls :frames frames))
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
         v0 (if root (vec (or (:vel b) (repeat 6 0.0))) (vec (repeat 6 0.0)))
         vels (spatial-velocities ls v0 qd)
         part (fn [inertia vi rot pos mass com]
                (let [ke (* 0.5 (reduce + (map * vi (lin/mat-vec inertia vi))))
                      c (v/add pos (q/rotate rot com))]
                  (- ke (* (double mass) (v/dot g c)))))]
     (+ (if root
          (part (spatial-inertia (:mass root) (:com root) (:inertia root))
                v0 (or (:rot b) q/identity-q) (or (:pos b) v/zero)
                (:mass root) (:com root))
          0.0)
        (reduce +
                (map (fn [i]
                       (let [link (nth (chain model) i)
                             {:keys [rot pos]} (nth fs i)]
                         (part (:i (nth ls i)) (nth vels i) rot pos
                               (:mass link) (:com link))))
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
        v0 (if root (vec (or (:vel b) (repeat 6 0.0))) (vec (repeat 6 0.0)))
        vels (spatial-velocities ls v0 qd)
        ;; Each body's spatial momentum `I v` is in its own frame and
        ;; about its own origin. Both have to be carried to the world
        ;; before they can be added up.
        one (fn [inertia vi {:keys [rot pos]}]
              (let [h (lin/mat-vec inertia vi)
                    ang (q/rotate rot (subvec (vec h) 0 3))
                    lin (q/rotate rot (subvec (vec h) 3 6))]
                [lin (v/add ang (v/cross pos lin))]))
        parts (cond-> (mapv (fn [i] (one (:i (nth ls i)) (nth vels i) (nth fs i)))
                            (range (count ls)))
                root (conj (one (spatial-inertia (:mass root) (:com root) (:inertia root))
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
  [vels v0 {:keys [rot pos]} i p]
  (let [vi (vec (if (neg? (long i)) v0 (nth vels i)))
        w (q/rotate rot (subvec vi 0 3))
        at-origin (q/rotate rot (subvec vi 3 6))]
    (v/add at-origin (v/cross w (v/sub p pos)))))

(defn point-velocity
  "How fast the point of body `i` that is at world position `p` is moving.
  `-1` is the free root itself.

  A point *fixed to the body*, not a point in space: the body's own
  material carries it, so it moves with its spin as well as its travel."
  [model q-vec state i p]
  (let [ls (links model q-vec)
        root-state (:base state)
        v0 (vec (or (:vel root-state) (repeat 6 0.0)))
        vels (spatial-velocities ls v0 (:qd state))]
    (velocity-at vels v0 (frame-of model q-vec root-state i) i p)))

(defn- jacobian*
  "The Jacobian, given the per-configuration data already computed."
  [model ls frame i p]
  (let [n (generalised-dof model)
        column (fn [j]
                 (let [u (assoc (vec (repeat n 0.0)) j 1.0)
                       [v0 qd] (split model u)]
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
  [{:keys [rot pos]} p dir]
  (let [inv (q/conjugate rot)
        f (q/rotate inv dir)
        r (q/rotate inv (v/sub p pos))]
    (vec (concat (v/cross r f) f))))

(defn- generalised-force
  "The generalised force a spatial force `f` on body `i` makes.

  A walk from the body it acts on up to the root, taking each joint's
  share along its own axis and carrying the rest to the parent. Bodies
  off that path feel nothing at all, which is why this is a walk and not
  a matrix -- a hand pushed sideways says nothing about the other arm."
  [model ls i f]
  (let [nd (reduce + (map :ndof ls))
        walk (loop [j (long i) f (vec f) acc (vec (repeat nd 0.0))]
               (if (neg? j)
                 {:root f :joints acc}
                 (let [{:keys [s parent xt offset]} (nth ls j)
                       share (st-times s f)]
                   (recur (long parent)
                          (lin/mat-vec xt f)
                          (reduce (fn [a k] (assoc a (+ (long offset) k) (nth share k)))
                                  acc
                                  (range (count s)))))))]
    (if (base model)
      (vec (concat (:root walk) (:joints walk)))
      (vec (:joints walk)))))

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
  [model ls ai i f]
  (let [n (count ls)
        root (base model)
        zero6 (vec (repeat 6 0.0))
        ;; Only the bodies between the one that was hit and the root can
        ;; have any bias force on them, so the inward pass is that walk
        ;; and not a sweep of everything. For a ragdoll -- a shallow
        ;; tree, not a chain -- that is three or four links out of
        ;; sixteen. Every other link's share is zero and stays zero,
        ;; which is what the outward pass then reads.
        path (loop [j (long i) acc []]
               (if (neg? j) acc (recur (long (:parent (nth ls j))) (conj acc j))))
        inward (reduce
                (fn [{:keys [pa pa0 uu] :as acc} j]
                  (let [{:keys [s parent xt]} (nth ls j)
                        {:keys [u dinv]} (nth (:links ai) j)
                        paj (get pa j zero6)
                        uj (mapv - (st-times s paj))
                        acc (assoc acc :uu (assoc uu j uj))
                        up (lin/mat-vec
                            xt
                            (cond-> paj
                              dinv (as-> x
                                         (reduce (fn [a l]
                                                   (mapv + a (scaled (nth u l)
                                                                     (nth (lin/mat-vec dinv uj) l))))
                                                 x
                                                 (range (count u))))))]
                    (if (neg? (long parent))
                      (assoc acc :pa0 (mapv + pa0 up))
                      (assoc acc :pa (update pa parent #(mapv + (or % zero6) up))))))
                {:pa (if (neg? (long i)) {} {i (mapv - f)})
                 :pa0 (if (neg? (long i)) (mapv - f) zero6)
                 :uu (vec (repeat n nil))}
                path)
        dv0 (if root
              (let [m (:ia0 ai)
                    sym (lin/mat-scale (lin/mat-add m (lin/transpose m)) 0.5)]
                (mapv - (lin/cholesky-solve sym (:pa0 inward))))
              zero6)
        outward (reduce
                 (fn [{:keys [dv] :as acc} j]
                   (let [{:keys [xup s parent]} (nth ls j)
                         {:keys [u dinv]} (nth (:links ai) j)
                         a' (lin/mat-vec xup (if (neg? (long parent)) dv0 (nth dv parent)))
                         ;; Bodies off the path to the impulse have no
                         ;; bias force on them at all, which is what the
                         ;; nil says; their joints still accelerate,
                         ;; because everything above them moved.
                         uj (or (nth (:uu inward) j) (vec (repeat (count s) 0.0)))
                         dq (if dinv
                              (lin/mat-vec dinv
                                           (mapv (fn [k ul]
                                                   (- (double k) (reduce + (map * ul a'))))
                                                 uj u))
                              (vec (repeat (count s) 0.0)))]
                     (-> acc
                         (update :dv conj (mapv + a' (s-times s dq 0)))
                         (update :dqd into dq))))
                 {:dv [] :dqd []}
                 (range n))]
    (if root
      (vec (concat dv0 (:dqd outward)))
      (vec (:dqd outward)))))

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
        w (reduce + (map * g delta))]
    {:g g :delta delta :m (if (> w 1e-12) (/ 1.0 w) 0.0)}))

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
    {:delta-u delta
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

(defn solve-contacts
  "`state` with its velocities corrected so the contacts are not being
  driven into.

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
        u0 (velocity model state)
        prep (mapv (fn [c]
                     (let [i (:link c)
                           frame (frame-of model q root-state i frames)
                           n (v/normalize (:normal c))
                           [t1 t2] (tangents n)
                           along #(response model ls ai frame i (:point c) %)
                           dirs {:n (along n) :t1 (along t1) :t2 (along t2)}]
                       (assoc dirs
                              ;; The closing speed as the step began.
                              ;; Restitution is measured against this and
                              ;; nothing later, for the same reason the
                              ;; rigid solvers record it once.
                              :approach (reduce + (map * (:g (:n dirs)) u0))
                              :bias (min (double max-push-speed)
                                         (/ (* (double bias-factor)
                                               (max 0.0 (- (double (:depth c))
                                                           (double slop))))
                                            dt)))))
                   cs)
        k (count prep)]
    (if (zero? k)
      state
      (with-velocity
        model state
        (loop [u u0
               acc (vec (repeat k [0.0 0.0 0.0]))
               iter 0]
          (if (= iter (long iterations))
            u
            (let [[u acc]
                  (reduce
                   (fn [[u acc] i]
                     (let [{:keys [n t1 t2 approach bias]} (nth prep i)
                           [an at1 at2] (nth acc i)
                           ;; Normal first: friction is bounded by the
                           ;; force it rides on, so it needs this
                           ;; iteration's answer and not last one's.
                           target (max (double bias)
                                       (if (< (double approach) -0.5)
                                         (* (- (double restitution)) (double approach))
                                         0.0))
                           vn (reduce + (map * (:g n) u))
                           dn (* (double (:m n)) (- target vn))
                           an' (max 0.0 (+ (double an) dn))
                           u (mapv #(+ (double %1) (* (double %2) (- an' (double an))))
                                   u (:delta n))
                           limit (* (double friction) an')
                           [u at1'] (let [vt (reduce + (map * (:g t1) u))
                                          d (* (double (:m t1)) (- vt))
                                          a' (min limit (max (- limit) (+ (double at1) d)))]
                                      [(mapv #(+ (double %1)
                                                 (* (double %2) (- a' (double at1))))
                                             u (:delta t1))
                                       a'])
                           [u at2'] (let [vt (reduce + (map * (:g t2) u))
                                          d (* (double (:m t2)) (- vt))
                                          a' (min limit (max (- limit) (+ (double at2) d)))]
                                      [(mapv #(+ (double %1)
                                                 (* (double %2) (- a' (double at2))))
                                             u (:delta t2))
                                       a'])]
                       [u (assoc acc i [an' at1' at2'])]))
                   [u acc]
                   (range k))]
              (recur u acc (inc iter)))))))))
