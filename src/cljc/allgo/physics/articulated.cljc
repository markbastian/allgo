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

  ## Being hit

  It can be hit, though, which took one more thing. A contact solver
  asks a body two questions -- how much velocity does a push here buy,
  and please take this impulse -- and neither has an obvious answer for
  a body whose motion is described by joint angles.

  Both come out of `H^-1`, the generalised inverse inertia, which the
  articulated body algorithm hands over a column at a time: nothing
  moving, no gravity, one unit of generalised force, and the
  accelerations that result *are* that column. With a Jacobian for the
  contact point, `1 / (d^T J H^-1 J^T d)` is the mass felt there and
  `H^-1 J^T d` is what an impulse does to the joint rates.

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
  numbers and find the singularity in it. The joints keep their one
  number each and gain nothing from the base being free; the base gains
  one 6x6 solve, and that is the whole of what a floating base costs.

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
            [allgo.numerics.linear :as lin]))

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

(defn- joint-transform
  "How the child frame sits relative to the joint frame at coordinate `x`."
  [{:keys [joint axis]} ^double x]
  (case (or joint :revolute)
    :revolute {:rot (q/from-axis-angle axis x) :pos v/zero}
    :prismatic {:rot q/identity-q :pos (v/scale axis x)}
    :fixed {:rot q/identity-q :pos v/zero}))

(defn- subspace
  "The joint's motion subspace: the spatial velocity one unit of joint
  rate produces, in the child's own coordinates."
  [{:keys [joint axis]}]
  (let [[ax ay az] axis]
    (case (or joint :revolute)
      :revolute [(double ax) (double ay) (double az) 0.0 0.0 0.0]
      :prismatic [0.0 0.0 0.0 (double ax) (double ay) (double az)]
      :fixed [0.0 0.0 0.0 0.0 0.0 0.0])))

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
  "How many joint coordinates the model has -- one per link, since every
  joint here is a single degree of freedom. The base's six, if it has
  any, are not among them: they are a pose and a velocity, not
  coordinates, precisely so that no one has to parameterise a rotation."
  [model]
  (count (chain model)))

;; ---------------------------------------------------------------------------
;; Shared per-link setup

(defn- links
  "Each link's transform from its parent, motion subspace and inertia, at
  the configuration `q-vec`."
  [model q-vec]
  (mapv (fn [link ^double x]
          (let [s (subspace link)]
            {:xup (lin/mat-mul (transform (joint-transform link x))
                               (transform (:origin link)))
             :s s
             :i (spatial-inertia (:mass link) (:com link) (:inertia link))
             :parent (long (:parent link))}))
        (chain model)
        q-vec))

(defn- scaled [s ^double x] (mapv #(* (double %) x) s))

(defn- spatial-velocities
  "Every link's spatial velocity, in its own coordinates.

  Linear in `(v0, qd)` and used three ways because of it: to say what a
  model is doing now, and -- one unit of generalised velocity at a time
  -- to build the Jacobian of any point on it."
  [ls v0 qd]
  (reduce (fn [acc i]
            (let [{:keys [xup s parent]} (nth ls i)
                  vp (if (neg? (long parent)) v0 (nth acc parent))]
              (conj acc (mapv + (lin/mat-vec xup vp) (scaled s (nth qd i))))))
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
                (let [{:keys [xup s parent]} (nth ls i)
                      inertia (:i (nth ls i))
                      vj (scaled s (nth qd i))
                      vp (if (neg? parent) (vec (repeat 6 0.0)) (:v (nth acc parent)))
                      ap (if (neg? parent) (base-acceleration gravity) (:a (nth acc parent)))
                      vi (mapv + (lin/mat-vec xup vp) vj)
                      ai (mapv + (lin/mat-vec xup ap)
                               (scaled s (nth qdd i))
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
     (mapv (fn [i] (reduce + (map * (:s (nth ls i)) (nth forces i)))) (range n)))))

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

(defn- outer
  "The 6x6 `a b^T`."
  [a b]
  (mapv (fn [ai] (mapv (fn [bj] (* (double ai) (double bj))) b)) a))

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
                        ext :external root-force :base-force}]
   (let [ls (links model q-vec)
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
         ;; Pass one, outward: velocity, the acceleration that velocity
         ;; alone implies, and the force needed to hold the link on it.
         pass1 (reduce
                (fn [acc i]
                  (let [{:keys [xup s parent]} (nth ls i)
                        inertia (:i (nth ls i))
                        vj (scaled s (nth qd i))
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
                               :ia inertia
                               :pa (let [bias (lin/mat-vec (crf vi)
                                                           (lin/mat-vec inertia vi))]
                                     (if f (mapv - bias f) bias))})))
                []
                (range n))
         ;; Pass two, inward: what the subtree beyond each joint looks
         ;; like to the link above it, and what the whole of it looks
         ;; like to the root.
         pass2 (reduce
                (fn [{:keys [links ia0 pa0] :as acc} i]
                  (let [{:keys [xup s parent]} (nth ls i)
                        {:keys [ia pa c]} (nth links i)
                        u (lin/mat-vec ia s)
                        d (reduce + (map * s u))
                        uu (- (double (nth tau i)) (reduce + (map * s pa)))
                        acc (assoc acc :links (-> links
                                                  (assoc-in [i :u] u)
                                                  (assoc-in [i :d] d)
                                                  (assoc-in [i :uu] uu)))]
                    (if (< (abs d) 1e-12)
                      acc
                      (let [ia' (lin/mat-sub ia (lin/mat-scale (outer u u) (/ 1.0 d)))
                            pa' (mapv + pa (lin/mat-vec ia' c) (scaled u (/ uu d)))
                            xt (lin/transpose xup)
                            up-i (lin/mat-mul xt (lin/mat-mul ia' xup))
                            up-p (lin/mat-vec xt pa')]
                        (if (neg? parent)
                          (assoc acc
                                 :ia0 (when ia0 (lin/mat-add ia0 up-i))
                                 :pa0 (when pa0 (mapv + pa0 up-p)))
                          (update acc :links
                                  #(-> %
                                       (update-in [parent :ia] (fn [m] (lin/mat-add m up-i)))
                                       (update-in [parent :pa] (fn [v] (mapv + v up-p))))))))))
                {:links pass1
                 :ia0 root-i
                 :pa0 (when root
                        (let [bias (lin/mat-vec (crf v0) (lin/mat-vec root-i v0))]
                          (if root-force (mapv - bias root-force) bias)))}
                (reverse (range n)))
         solved (:links pass2)
         ;; The root, if it is free: force-free in the falling frame.
         a0 (if root
              (let [m (:ia0 pass2)
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
                      {:keys [c u d uu]} (nth solved i)
                      ap (if (neg? parent) a0 (nth a parent))
                      a' (mapv + (lin/mat-vec xup ap) c)
                      qddi (if (< (abs (double d)) 1e-12)
                             0.0
                             (/ (- (double uu) (reduce + (map * u a'))) (double d)))]
                  (-> acc
                      (update :a conj (mapv + a' (scaled s qddi)))
                      (update :qdd conj qddi))))
              {:a [] :qdd []}
              (range n))]
     {:qdd (:qdd out)
      ;; Back out of the falling frame, which is where the root's own
      ;; share of gravity comes from.
      :base-acc (if root (mapv + a0 a-grav) (vec (repeat 6 0.0)))})))

;; ---------------------------------------------------------------------------
;; Moving it

(defn step
  "One semi-implicit Euler step of `{:q :qd}`, and `:base` if there is
  one, under `tau`.

  Velocity first and then position, which is the same choice
  `allgo.physics.rigid` makes and for the same reason: it costs nothing
  and it does not pump energy into an oscillation the way explicit Euler
  does.

  The base's state is a pose and a spatial velocity in its own frame,
  not six more coordinates. Three reasons, and the third is the one that
  matters: a rotation has no three-number parameterisation without a
  singularity in it; a quaternion is the shape the rest of this library
  speaks; and a velocity in body coordinates is what the algorithms
  already produce, so nothing has to be converted on the way in."
  ([state model dt] (step state model dt nil))
  ([{:keys [q qd] :as state} model ^double dt {:keys [tau] :as opts}]
   (let [root (base model)
         b (:base state)
         tau (or tau (vec (repeat (dof model) 0.0)))
         {:keys [qdd base-acc]} (forward-dynamics model q qd tau
                                                  (cond-> opts root (assoc :base b)))
         qd' (mapv (fn [a c] (+ (double a) (* dt (double c)))) qd qdd)
         q' (mapv (fn [a c] (+ (double a) (* dt (double c)))) q qd')]
     (if-not root
       {:q q' :qd qd'}
       (let [v (vec (or (:vel b) (repeat 6 0.0)))
             v' (mapv (fn [a c] (+ (double a) (* dt (double c)))) v base-acc)
             rot (or (:rot b) q/identity-q)
             pos (or (:pos b) v/zero)
             ;; Both halves of the spatial velocity are in the base's own
             ;; frame, so both are turned into the world before they move
             ;; anything.
             w (q/rotate rot (subvec v' 0 3))
             u (q/rotate rot (subvec v' 3 6))
             ;; rot += dt/2 (omega, 0) rot, then back onto the unit
             ;; sphere -- the same first-order step `allgo.physics.rigid`
             ;; takes, and it walks off for the same reason.
             [dx dy dz dw] (q/mul [(w 0) (w 1) (w 2) 0.0] rot)
             [rx ry rz rw] rot
             h (* 0.5 dt)]
         {:q q' :qd qd'
          :base {:rot (q/normalize [(+ rx (* h dx)) (+ ry (* h dy))
                                    (+ rz (* h dz)) (+ rw (* h dw))])
                 :pos (v/add-scaled pos u dt)
                 :vel v'}})))))

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
             {jrot :rot jpos :pos} (joint-transform link (nth q-vec i))
             rot (q/mul (q/mul prot (or orot q/identity-q)) jrot)
             pos (v/add (v/add ppos (q/rotate prot (or opos v/zero)))
                        (q/rotate (q/mul prot (or orot q/identity-q)) jpos))
             frame {:rot rot :pos pos}]
         (-> acc (update :frames conj frame) (update :out conj frame))))
     {:frames [] :out []}
     (range (dof model))))))

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
  "Where body `i` is, in world terms. `-1` is the free root itself."
  [model q-vec root-state i]
  (if (neg? (long i))
    {:rot (or (:rot root-state) q/identity-q) :pos (or (:pos root-state) v/zero)}
    (nth (poses model q-vec root-state) i)))

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

(defn point-jacobian
  "The 3 by `generalised-dof` matrix taking a generalised velocity to the
  world velocity of the point of link `i` at `p`.

  Built a column at a time, by asking what one unit of each generalised
  velocity on its own does. That is not an approximation -- the map is
  linear, which is the whole reason a Jacobian exists -- and it costs a
  velocity recursion per column where a purpose-built one would walk the
  path from the root once. At a ragdoll's twenty-odd degrees of freedom
  the difference is not worth the second implementation to get wrong."
  [model q-vec root-state i p]
  (let [ls (links model q-vec)
        n (generalised-dof model)
        frame (frame-of model q-vec root-state i)
        column (fn [j]
                 (let [u (assoc (vec (repeat n 0.0)) j 1.0)
                       [v0 qd] (split model u)]
                   (velocity-at (spatial-velocities ls v0 qd) v0 frame i p)))]
    (lin/transpose (mapv column (range n)))))

(defn inverse-mass-matrix
  "`H^-1`, the generalised inverse inertia, a column at a time out of the
  articulated body algorithm.

  Nothing is moving and there is no gravity, so the accelerations a unit
  generalised force produces *are* that column of the inverse. The
  inverse rather than the matrix itself because everything asked of it
  here -- how hard is this point to push, what does an impulse do -- is
  a question about the inverse, and forming H only to factor it again
  would be work in both directions."
  [model q-vec]
  (let [n (generalised-dof model)
        root (base model)
        nj (dof model)
        zero-q (vec (repeat nj 0.0))
        free {:gravity [0.0 0.0 0.0]}
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
    (lin/transpose (mapv column (range n)))))

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
  (let [j (point-jacobian model q-vec root-state i p)
        hinv (inverse-mass-matrix model q-vec)
        ;; J^T d: the generalised force a unit impulse along `dir` makes.
        jtd (lin/mat-vec (lin/transpose j) (vec dir))
        delta (lin/mat-vec hinv jtd)
        w (reduce + (map * jtd delta))]
    {:delta-u delta
     :effective-mass (if (> w 1e-12) (/ 1.0 w) ##Inf)}))

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
