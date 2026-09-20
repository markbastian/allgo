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

  What it gives up is everything the constraint formulation was good at.
  A reduced-coordinate chain cannot be hit by a brick, cannot be taken
  apart, and cannot close a loop -- a tree of joints is the whole of what
  it can describe. Contacts are the interesting half of that and are not
  here yet.

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
  is the base, which is fixed: a floating base is a six degree of freedom
  joint to the world and is not here yet.

  ## What runs on it

  `inverse-dynamics` is the recursive Newton-Euler algorithm: given where
  the joints are, how fast they are moving and how fast they are being
  accelerated, what torques does that take. `forward-dynamics` is the
  articulated body algorithm: given the torques, what accelerations. Both
  are O(n) and the second is the one a simulation wants.

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

(defn dof
  "How many coordinates the model has -- one per link, since every joint
  here is a single degree of freedom."
  [model]
  (count model))

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
        model
        q-vec))

(defn- scaled [s ^double x] (mapv #(* (double %) x) s))

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
  "The accelerations `qdd` that torques `tau` produce, in O(n).

  Featherstone's articulated body algorithm. The idea it turns on: seen
  from one joint, everything beyond it behaves like a single rigid body
  with an inertia of its own -- an *articulated* inertia, which is what
  the subtree's own joints being free to move does to the inertia it
  presents. Compute that inward from the tips and every joint's
  acceleration follows from its parent's outward, with no matrix to
  invert bigger than one degree of freedom.

  Three passes: outward for velocities and the bias force each link's own
  motion generates, inward accumulating articulated inertias and their
  bias forces onto the parents, outward again for the accelerations."
  ([model q-vec qd tau] (forward-dynamics model q-vec qd tau nil))
  ([model q-vec qd tau {:keys [gravity]}]
   (let [ls (links model q-vec)
         n (count ls)
         ;; Pass one, outward: velocity, the acceleration that velocity
         ;; alone implies, and the force needed to hold the link on it.
         pass1 (reduce
                (fn [acc i]
                  (let [{:keys [xup s parent]} (nth ls i)
                        inertia (:i (nth ls i))
                        vj (scaled s (nth qd i))
                        vp (if (neg? parent) (vec (repeat 6 0.0)) (:v (nth acc parent)))
                        vi (mapv + (lin/mat-vec xup vp) vj)]
                    (conj acc {:v vi
                               :c (lin/mat-vec (crm vi) vj)
                               :ia inertia
                               :pa (lin/mat-vec (crf vi) (lin/mat-vec inertia vi))})))
                []
                (range n))
         ;; Pass two, inward: what the subtree beyond each joint looks
         ;; like to the link above it.
         pass2 (reduce
                (fn [acc i]
                  (let [{:keys [xup s parent]} (nth ls i)
                        {:keys [ia pa c]} (nth acc i)
                        u (lin/mat-vec ia s)
                        d (reduce + (map * s u))
                        uu (- (double (nth tau i)) (reduce + (map * s pa)))
                        acc (assoc-in acc [i :u] u)
                        acc (assoc-in acc [i :d] d)
                        acc (assoc-in acc [i :uu] uu)]
                    (if (or (neg? parent) (< (abs d) 1e-12))
                      acc
                      (let [ia' (lin/mat-sub ia (lin/mat-scale (outer u u) (/ 1.0 d)))
                            pa' (mapv + pa (lin/mat-vec ia' c) (scaled u (/ uu d)))
                            xt (lin/transpose xup)]
                        (-> acc
                            (update-in [parent :ia]
                                       #(lin/mat-add % (lin/mat-mul xt (lin/mat-mul ia' xup))))
                            (update-in [parent :pa]
                                       #(mapv + % (lin/mat-vec xt pa'))))))))
                pass1
                (reverse (range n)))]
     ;; Pass three, outward: each joint's acceleration, then the link's.
     (:qdd
      (reduce
       (fn [{:keys [a] :as acc} i]
         (let [{:keys [xup s parent]} (nth ls i)
               {:keys [c u d uu]} (nth pass2 i)
               ap (if (neg? parent) (base-acceleration gravity) (nth a parent))
               a' (mapv + (lin/mat-vec xup ap) c)
               qddi (if (< (abs (double d)) 1e-12)
                      0.0
                      (/ (- (double uu) (reduce + (map * u a'))) (double d)))]
           (-> acc
               (update :a conj (mapv + a' (scaled s qddi)))
               (update :qdd conj qddi))))
       {:a [] :qdd []}
       (range n))))))

;; ---------------------------------------------------------------------------
;; Moving it

(defn step
  "One semi-implicit Euler step of `{:q :qd}` under `tau`.

  Velocity first and then position, which is the same choice
  `allgo.physics.rigid` makes and for the same reason: it costs nothing
  and it does not pump energy into an oscillation the way explicit Euler
  does."
  ([state model dt] (step state model dt nil))
  ([{:keys [q qd]} model ^double dt {:keys [tau] :as opts}]
   (let [tau (or tau (vec (repeat (dof model) 0.0)))
         qdd (forward-dynamics model q qd tau opts)
         qd' (mapv (fn [a b] (+ (double a) (* dt (double b)))) qd qdd)]
     {:q (mapv (fn [a b] (+ (double a) (* dt (double b)))) q qd')
      :qd qd'})))

;; ---------------------------------------------------------------------------
;; Where the links are

(defn poses
  "Each link's frame in world coordinates, as `{:rot :pos}`.

  Reduced coordinates describe an arm by its angles, and something has to
  turn those back into places to draw."
  [model q-vec]
  (:out
   (reduce
    (fn [{:keys [frames] :as acc} i]
      (let [link (nth model i)
            parent (long (:parent link))
            {prot :rot ppos :pos} (if (neg? parent)
                                    {:rot q/identity-q :pos v/zero}
                                    (nth frames parent))
            {orot :rot opos :pos} (:origin link)
            {jrot :rot jpos :pos} (joint-transform link (nth q-vec i))
            rot (q/mul (q/mul prot (or orot q/identity-q)) jrot)
            pos (v/add (v/add ppos (q/rotate prot (or opos v/zero)))
                       (q/rotate (q/mul prot (or orot q/identity-q)) jpos))
            frame {:rot rot :pos pos}]
        (-> acc (update :frames conj frame) (update :out conj frame))))
    {:frames [] :out []}
    (range (count model)))))

(defn energy
  "Kinetic plus potential, for checking that nothing is being invented.

  A chain under no torque and no damping must keep this constant, and an
  error in the spatial algebra almost always shows up here before it
  shows up anywhere a person would notice."
  ([model q-vec qd] (energy model q-vec qd nil))
  ([model q-vec qd {:keys [gravity]}]
   (let [g (or gravity [0.0 -9.81 0.0])
         ls (links model q-vec)
         fs (poses model q-vec)
         vels (reduce (fn [acc i]
                        (let [{:keys [xup s parent]} (nth ls i)
                              vp (if (neg? parent) (vec (repeat 6 0.0)) (nth acc parent))]
                          (conj acc (mapv + (lin/mat-vec xup vp) (scaled s (nth qd i))))))
                      []
                      (range (count ls)))]
     (reduce
      +
      (map (fn [i]
             (let [vi (nth vels i)
                   inertia (:i (nth ls i))
                   ;; Half v^T I v, with both in the link's own frame.
                   ke (* 0.5 (reduce + (map * vi (lin/mat-vec inertia vi))))
                   {:keys [rot pos]} (nth fs i)
                   com (v/add pos (q/rotate rot (:com (nth model i))))
                   pe (- (* (double (:mass (nth model i))) (v/dot g com)))]
               (+ ke pe)))
           (range (count ls)))))))
