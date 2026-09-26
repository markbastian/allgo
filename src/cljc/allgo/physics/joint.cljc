(ns allgo.physics.joint
  "Joints between rigid bodies. After Ten Minute Physics 25.

  A joint is two frames -- a position and an orientation, one fixed to
  each body -- and a statement about how they are allowed to differ. Every
  kind here is some combination of three primitives:

    pin      hold the two origins together, which is a distance
             constraint with a target of zero
    align    drive one frame's axis onto the other's, by correcting
             about the axis perpendicular to both
    limit    hold an angle between two bounds, doing nothing while it is
             inside them

  From those: a ball joint pins and limits swing and twist; a hinge also
  aligns one axis, so only rotation about it survives; a servo is a hinge
  with a target angle; a motor is a servo whose target advances every
  step; prismatic and cylinder relax the pin along one axis so the bodies
  can slide; fixed aligns the orientations outright.

  These are `allgo.physics.rigid` constraints and nothing more -- each
  kind is a method on that namespace's `solve` multimethod, so a world can
  hold joints and plain distance constraints together and `rigid/step`
  does not need to know the difference.

  One departure from the reference, measured rather than assumed. It
  halves every position-level angular correction, under a comment
  reading stabilize rotation. That destabilizes these joints badly: on a ball joint
  pendulum at ten substeps the anchors, which should stay together, come
  18 units apart, against 0.005 at full strength. Halving is not a
  relaxation here because `allgo.physics.rigid/inverse-mass` sizes the
  correction assuming the *whole* angular response -- scale only the
  angular half and what is delivered no longer opposes the error. The
  giveaway is that a quarter behaves better than a half, which no genuine
  damping factor does. `:angular-relaxation` is a per-joint setting for
  anyone wanting the reference's behavior; it defaults to applying the
  correction in full."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]
            [clojure.math :as math]))

(def ^:private axis-x [1.0 0.0 0.0])
(def ^:private axis-y [0.0 1.0 0.0])

(def default-relaxation
  "How much of an angular correction to apply.

  One, meaning all of it, which is what keeps it consistent with the
  generalized inverse mass the correction was sized against. See the
  namespace docstring for why the reference's half is worse."
  1.0)

(defn- relaxation-of [j] (double (get j :angular-relaxation default-relaxation)))

;; ---------------------------------------------------------------------------
;; Frames

(defn frames
  "Where the joint's two frames are now: `[pos0 rot0 pos1 rot1]`.

  A frame is carried by its body, so it moves and turns with it. With no
  second body the far frame stays where it was put, which is how a joint
  anchors to the world."
  [bodies {:keys [a b local-pos-a local-rot-a local-pos-b local-rot-b]}]
  (let [ba (bodies a)]
    [(rigid/local->world ba local-pos-a)
     (q/mul (:rot ba) local-rot-a)
     (if b (rigid/local->world (bodies b) local-pos-b) local-pos-b)
     (if b (q/mul (:rot (bodies b)) local-rot-b) local-rot-b)]))

(defn joint
  "A joint holding two bodies at a frame given in world terms.

  `at` is where the joint sits and `rot` which way it faces; the joint's
  own x axis is the one a hinge turns about and a prismatic slides along,
  so `rot` is how you aim it. `b` may be nil to anchor to the world.

  Everything beyond that is per kind: `:target-distance` and
  `:distance-min`/`:distance-max` for the sliding kinds, `:swing-min`,
  `:swing-max`, `:twist-min`, `:twist-max` for the limits, `:target-angle`
  for a servo and `:velocity` for a motor."
  [bodies {:keys [kind a b at rot] :as opts}]
  (let [at (vec (or at v/zero))
        rot (vec (or rot q/identity-q))
        ba (bodies a)]
    (merge {:swing-min ##-Inf :swing-max ##Inf
            :twist-min ##-Inf :twist-max ##Inf
            :distance-min ##-Inf :distance-max ##Inf
            :target-distance 0.0
            :distance-compliance 0.0
            :target-angle-compliance 0.0
            :velocity 0.0
            :compliance 0.0}
           opts
           {:kind kind :a a :b b
            :local-pos-a (rigid/world->local ba at)
            :local-rot-a (q/mul (:inv-rot ba) rot)
            :local-pos-b (if b (rigid/world->local (bodies b) at) at)
            :local-rot-b (if b (q/mul (:inv-rot (bodies b)) rot) rot)})))

;; ---------------------------------------------------------------------------
;; Primitives

(defn- pin
  "Holds the two frame origins `target` apart, zero being coincident."
  [bodies {:keys [a b compliance] :as j} target dt]
  (let [[p0 _ p1 _] (frames bodies j)
        delta (v/sub p1 p0)
        d (v/length delta)
        target (double target)]
    (if (and (zero? d) (zero? target))
      {:bodies bodies :force 0.0}
      (let [n (if (zero? d)
                ;; Coincident but wanting to be apart: any direction will
                ;; do, so take the frame's own.
                (q/rotate (second (frames bodies j)) [0.0 0.0 1.0])
                (v/scale delta (/ 1.0 d)))]
        (rigid/correct bodies {:a a :b b
                               ;; Toward each other when too far apart:
                               ;; `delta` runs from the first frame to the
                               ;; second, so a positive multiple of it
                               ;; closes the gap. The opposite sign drives
                               ;; them apart and diverges within a second.
                               :corr (v/scale n (- d target))
                               :at p0 :other-at p1
                               :compliance (or compliance 0.0)
                               :dt dt
                               :angular-relaxation (relaxation-of j)})))))

(defn- align
  "Drives `a0` onto `a1` by correcting about the axis perpendicular to
  both -- which is exactly what their cross product is."
  [bodies {:keys [a b] :as j} a0 a1 compliance dt]
  (:bodies (rigid/correct bodies {:a a :b b
                                  :corr (v/cross a0 a1)
                                  :at nil :other-at nil
                                  :compliance compliance
                                  :dt dt
                                  :angular-relaxation (relaxation-of j)})))

(defn signed-angle
  "The angle from `a` to `b` about `n`, in radians, within a half turn.

  The sine alone only distinguishes a quarter turn either way, so the
  cosine settles which half of the circle it is in."
  ^double [n a b]
  (let [s (v/dot (v/cross a b) n)
        phi (math/asin (min 1.0 (max -1.0 s)))
        phi (if (neg? (v/dot a b)) (- math/PI phi) phi)]
    (cond (> phi math/PI) (- phi (* 2.0 math/PI))
          (< phi (- math/PI)) (+ phi (* 2.0 math/PI))
          :else phi)))

(defn- limit
  "Holds the angle from `a0` to `a1` about `n` within bounds.

  Does nothing while the angle is inside them, which is what makes it a
  limit rather than a spring: a hinge with limits swings freely until it
  hits one."
  [bodies j n a0 a1 lo hi compliance dt]
  (let [phi (signed-angle n a0 a1)]
    (if (<= (double lo) phi (double hi))
      bodies
      (let [phi (min (double hi) (max (double lo) phi))
            ;; Where a0 would be at the limit; correcting a1 onto that is
            ;; the same align as before.
            target (q/rotate (q/from-axis-angle n phi) a0)]
        (align bodies j target a1 compliance dt)))))

(defn- limited?
  [lo hi]
  (or (> (double lo) ##-Inf) (< (double hi) ##Inf)))

;; ---------------------------------------------------------------------------
;; Position

(defn- solve-slide
  "Lets the bodies separate along the joint's x axis, and holds them
  together across it.

  The offset is taken into the joint's own frame, the allowed component
  cleared, and the rest put back -- so what is left to correct is only the
  part the joint does not permit."
  [bodies {:keys [a b kind target-distance distance-min distance-max] :as j} dt]
  (let [[p0 r0 p1 _] (frames bodies j)
        local (q/rotate (q/conjugate r0) (v/sub p1 p0))
        [lx ly lz] local
        target (min (double distance-max) (max (double distance-min) (double target-distance)))
        lx' (cond
              (= kind :cylinder) (- lx target)
              (> lx (double distance-max)) (- lx (double distance-max))
              (< lx (double distance-min)) (- lx (double distance-min))
              :else 0.0)]
    (rigid/correct bodies {:a a :b b
                           :corr (q/rotate r0 [lx' ly lz])
                           :at p0 :other-at p1
                           :compliance 0.0
                           :dt dt
                           :angular-relaxation (relaxation-of j)})))

(defn- solve-position
  [bodies {:keys [kind target-distance] :as j} dt]
  (case kind
    (:prismatic :cylinder) (solve-slide bodies j dt)
    (pin bodies j (or target-distance 0.0) dt)))

;; ---------------------------------------------------------------------------
;; Orientation

(defn- hinge-orientation
  [bodies {:keys [target-angle has-target-angle? target-angle-compliance
                  swing-min swing-max] :as j} dt]
  (let [bodies (let [[_ r0 _ r1] (frames bodies j)]
                 ;; Only rotation about the hinge axis survives, so the
                 ;; two frames' x axes are driven together.
                 (align bodies j (q/rotate r0 axis-x) (q/rotate r1 axis-x) 0.0 dt))
        bodies (if has-target-angle?
                 (let [[_ r0 _ r1] (frames bodies j)
                       n (q/rotate r0 axis-x)]
                   (limit bodies j n (q/rotate r0 axis-y) (q/rotate r1 axis-y)
                          target-angle target-angle (or target-angle-compliance 0.0) dt))
                 bodies)]
    (if (limited? swing-min swing-max)
      (let [[_ r0 _ r1] (frames bodies j)
            n (q/rotate r0 axis-x)]
        (limit bodies j n (q/rotate r0 axis-y) (q/rotate r1 axis-y)
               swing-min swing-max 0.0 dt))
      bodies)))

(defn- ball-orientation
  "Swing and twist, separated.

  Swing is how far the two x axes have parted; twist is rotation about
  them once that is taken out. Separating them is what lets a shoulder be
  described -- free to swing widely, barely free to twist."
  [bodies {:keys [swing-min swing-max twist-min twist-max] :as j} dt]
  (let [bodies (if (limited? swing-min swing-max)
                 (let [[_ r0 _ r1] (frames bodies j)
                       a0 (q/rotate r0 axis-x)
                       a1 (q/rotate r1 axis-x)
                       n (v/normalize (v/cross a0 a1))]
                   (if (v/zero-vector? n)
                     bodies
                     (limit bodies j n a0 a1 swing-min swing-max 0.0 dt)))
                 bodies)]
    (if (limited? twist-min twist-max)
      (let [[_ r0 _ r1] (frames bodies j)
            a0 (q/rotate r0 axis-x)
            a1 (q/rotate r1 axis-x)
            ;; The average of the two axes, which is the one to measure
            ;; twist about when they do not quite line up.
            n (v/normalize (v/add a0 a1))]
        (if (v/zero-vector? n)
          bodies
          (let [;; Project the reference axes onto the plane across the
                ;; twist axis, or what is measured is partly swing.
                flatten (fn [u] (v/normalize (v/add-scaled u n (- (v/dot n u)))))
                b0 (flatten (q/rotate r0 axis-y))
                b1 (flatten (q/rotate r1 axis-y))]
            (limit bodies j n b0 b1 twist-min twist-max 0.0 dt))))
      bodies)))

(defn- fixed-orientation
  [bodies {:keys [a b] :as j} dt]
  (let [[_ r0 _ r1] (frames bodies j)
        [dx dy dz dw] (q/mul r0 (q/conjugate r1))
        corr (v/scale [dx dy dz] 2.0)
        ;; A quaternion and its negation are the same orientation; the
        ;; sign says which way round the shorter correction lies.
        corr (if (pos? dw) (v/negate corr) corr)]
    (:bodies (rigid/correct bodies {:a a :b b :corr corr
                                    :at nil :other-at nil
                                    :compliance 0.0 :dt dt
                                    :angular-relaxation (relaxation-of j)}))))

(defn- solve-orientation
  [bodies {:keys [kind] :as j} dt]
  (case kind
    ;; A cylinder is a hinge that may also slide, so its rotation is
    ;; constrained the same way -- `:swing-min`/`:swing-max` bound the turn
    ;; about the axis, as they do for a hinge.
    (:hinge :servo :motor :cylinder) (hinge-orientation bodies j dt)
    :ball (ball-orientation bodies j dt)
    ;; A prismatic joint slides and does not turn at all, so its
    ;; orientation is locked outright. The reference sends it to the ball
    ;; branch, which applies limits and nothing else -- so a prismatic
    ;; joint with no limits set is free to rotate however it likes, which
    ;; is not a prismatic joint. Its scenes come from files that set the
    ;; limits to zero; a constructor should not need that of its caller.
    (:prismatic :fixed) (fixed-orientation bodies j dt)
    bodies))

;; ---------------------------------------------------------------------------
;; As rigid body constraints

(defn- solve-joint [bodies j dt]
  (let [{:keys [bodies force]} (solve-position bodies j dt)]
    {:bodies (solve-orientation bodies j dt) :force force}))

(doseq [kind [:ball :hinge :servo :fixed :prismatic :cylinder]]
  (defmethod rigid/solve kind [bodies j dt] (solve-joint bodies j dt)))

(defmethod rigid/solve :motor [bodies j dt]
  ;; A motor is a servo whose target keeps moving. The step is capped so
  ;; that a fast motor cannot ask for more than half a turn in one
  ;; substep, which the angle measure could not tell from turning the
  ;; other way.
  (solve-joint bodies j dt))

(defn advance-motors
  "Moves every motor's target on by one step of `dt`.

  Kept out of `solve` deliberately: solving runs once per substep and per
  iteration, and a target that advanced each time would turn at a rate
  depending on how the solver was configured."
  [constraints dt]
  (mapv (fn [{:keys [kind velocity] :as c}]
          (if (= :motor kind)
            (let [step (min 1.0 (max -1.0 (* (double velocity) (double dt))))]
              (-> c
                  (update :target-angle (fnil + 0.0) step)
                  (assoc :has-target-angle? true)))
            c))
        constraints))

(defn step
  "A frame of a world holding joints.

  The same as `allgo.physics.rigid/step`, with motor targets advanced
  first."
  ([world] (step world {}))
  ([world opts]
   (let [dt (double (:dt (merge rigid/default-world world opts)))]
     (rigid/step (update world :constraints advance-motors dt) opts))))

;; ---------------------------------------------------------------------------
;; Constructors

(defn ball
  "Two bodies pinned at a point, free to turn any way within its limits."
  [bodies opts] (joint bodies (assoc opts :kind :ball)))

(defn hinge
  "Free to turn about the joint frame's x axis and nothing else."
  [bodies opts] (joint bodies (assoc opts :kind :hinge)))

(defn servo
  "A hinge driven to `:target-angle`."
  [bodies opts]
  (joint bodies (assoc opts :kind :servo :has-target-angle? true)))

(defn motor
  "A hinge whose target angle advances at `:velocity` radians a second."
  [bodies opts]
  (joint bodies (assoc opts :kind :motor :has-target-angle? true)))

(defn fixed
  "No relative motion at all -- useful for welding two bodies into one."
  [bodies opts] (joint bodies (assoc opts :kind :fixed)))

(defn prismatic
  "Slides along the joint frame's x axis, between `:distance-min` and
  `:distance-max`, and does not turn."
  [bodies opts] (joint bodies (assoc opts :kind :prismatic)))

(defn cylinder
  "Slides along the joint frame's x axis and turns about it."
  [bodies opts] (joint bodies (assoc opts :kind :cylinder)))
