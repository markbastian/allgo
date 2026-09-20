(ns allgo.physics.rigid
  "Rigid bodies by position based dynamics. After Ten Minute Physics 22.

  `allgo.physics.xpbd` moves particles, which have a position and nothing
  else. A rigid body also has an orientation, and that is the whole
  difference: a correction applied away from the centre of mass has to
  turn the body as well as move it, by an amount that depends on which way
  the body is facing. Everything here follows from that.

  The loop is the one from 10, unchanged in shape: guess where everything
  goes, fix up the constraints, and read the velocity back off the
  positions.

    integrate          gravity, then move and turn by the current velocity
    solve              push positions and orientations to satisfy the
                       constraints
    update-velocities  velocity is wherever the body ended up, over dt

  Reading velocity back rather than tracking it is what makes this stable.
  A constraint cannot inject energy, because it only ever moves things;
  whatever velocity results is by definition consistent with where the
  body actually is.

  `inverse-mass` is the piece worth understanding. For a correction along
  `n` applied at world point `p`, it answers how much the body will give
  -- and for a body that is easy to push but hard to spin, that depends
  entirely on where `p` is relative to the centre. A correction through
  the centre sees only the mass; one at arm's length sees mostly the
  inertia. Passing no point at all asks the purely angular question, which
  is what an orientation constraint needs and what `allgo.physics.joint`
  is built on.

  Bodies are immutable maps held in a vector, and constraints name them by
  index. That is the opposite of the flat arrays elsewhere in this
  library, and deliberately: a scene has tens of rigid bodies where a
  fluid has thousands of particles, and at tens the clarity is worth more
  than the allocation."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]))

;; ---------------------------------------------------------------------------
;; Bodies

(defn- with-rotation
  "Keeps the inverse orientation alongside the orientation.

  Every correction needs it -- to take a world direction into the body's
  own frame, where the inertia is diagonal -- so it is computed once per
  change rather than at each use."
  [b rot]
  (let [rot (q/normalize rot)]
    (assoc b :rot rot :inv-rot (q/inverse rot))))

(defn body
  "A rigid body from its mass and principal inertia.

  `inertia` is the diagonal of the inertia tensor in the body's own
  frame. That it is only a diagonal is not a simplification: every body
  has axes about which the tensor is diagonal, and a shape built
  symmetrically is already aligned with its own."
  [{:keys [pos rot vel omega mass inertia damping angular-damping]
    :or   {pos v/zero rot q/identity-q vel v/zero omega v/zero
           damping 0.0 angular-damping 0.0}}]
  (let [mass (double (or mass 0.0))
        [ix iy iz] (or inertia v/zero)]
    (-> {:pos (vec pos) :vel (vec vel) :omega (vec omega)
         :prev-pos (vec pos) :prev-rot (vec rot)
         :mass mass
         :inv-mass (if (pos? mass) (/ 1.0 mass) 0.0)
         :inertia (vec (or inertia v/zero))
         :inv-inertia [(if (pos? (double ix)) (/ 1.0 (double ix)) 0.0)
                       (if (pos? (double iy)) (/ 1.0 (double iy)) 0.0)
                       (if (pos? (double iz)) (/ 1.0 (double iz)) 0.0)]
         :damping (double damping)
         :angular-damping (double angular-damping)}
        (with-rotation rot))))

(defn box
  "A solid box of `size` `[x y z]`.

  A density of zero, or none, makes it static -- infinitely heavy, which
  is how a floor or an anchor is spelled."
  [{:keys [size density] :as opts}]
  (let [[sx sy sz] (map double size)
        d (double (or density 0.0))
        mass (* d sx sy sz)]
    (assoc (body (assoc opts
                        :mass mass
                        :inertia (when (pos? mass)
                                   [(* (/ 1.0 12.0) mass (+ (* sy sy) (* sz sz)))
                                    (* (/ 1.0 12.0) mass (+ (* sx sx) (* sz sz)))
                                    (* (/ 1.0 12.0) mass (+ (* sx sx) (* sy sy)))])))
           :shape :box :size [sx sy sz])))

(defn ball
  "A solid sphere of `radius`."
  [{:keys [radius density] :as opts}]
  (let [r (double radius)
        d (double (or density 0.0))
        mass (* d (/ 4.0 3.0) Math/PI r r r)
        i (* 0.4 mass r r)]
    (assoc (body (assoc opts :mass mass :inertia (when (pos? mass) [i i i])))
           :shape :ball :radius r :size [(* 2 r) (* 2 r) (* 2 r)])))

(defn static?
  "A body nothing can move."
  [{:keys [inv-mass]}]
  (zero? (double inv-mass)))

(defn sleeping?
  "A body that has settled and been taken out of the simulation.

  `allgo.physics.solver` puts bodies to sleep and wakes them; the flag
  lives here because it is a fact about the body, and because collision
  detection needs to read it too -- a pair of sleeping bricks is a pair
  nothing need be asked about."
  [b]
  (true? (:sleeping? b)))

(defn inert?
  "A body that will not move this step, whether it cannot or need not."
  [b]
  (or (static? b) (sleeping? b)))

(defn wake
  "Puts a body back in the simulation, with its sleep clock reset.

  Anything moving a body from outside the step -- placing it, giving it a
  velocity, firing it at a wall -- has to do this, or it will be ignored
  until something bumps into it."
  [b]
  (assoc b :sleeping? false :sleep-time 0.0))

(defn local->world [{:keys [pos rot]} p] (v/add (q/rotate rot p) pos))

(defn world->local [{:keys [pos inv-rot]} p] (q/rotate inv-rot (v/sub p pos)))

(defn set-pose
  "Places a body, keeping its inverse orientation in step."
  [b pos rot]
  (-> b (assoc :pos (vec pos)) (with-rotation rot)))

;; ---------------------------------------------------------------------------
;; Integration

(defn integrate
  "Gravity, then move and turn by the current velocity.

  The orientation is stepped by `rot += dt/2 * (omega, 0) * rot`, which is
  the quaternion derivative for a body spinning at `omega`, and then
  renormalized -- a first-order step walks off the unit sphere, and an
  orientation that is not unit length is not a rotation."
  [{:keys [pos rot vel omega] :as b} dt gravity]
  (if (static? b)
    b
    (let [dt (double dt)
          vel (v/add-scaled vel gravity dt)
          pos' (v/add-scaled pos vel dt)
          [dx dy dz dw] (q/mul [(omega 0) (omega 1) (omega 2) 0.0] rot)
          [rx ry rz rw] rot
          h (* 0.5 dt)]
      (-> b
          (assoc :prev-pos pos :prev-rot rot :vel vel :pos pos')
          (with-rotation [(+ rx (* h dx)) (+ ry (* h dy))
                          (+ rz (* h dz)) (+ rw (* h dw))])))))

(defn update-velocities
  "Velocity is wherever the body ended up, over the step.

  The angular part comes from the rotation that took the body from where
  it was to where it is: its vector part is the axis scaled by the sine of
  half the angle, so twice it over `dt` is the angular velocity for small
  steps. A negative scalar part means the quaternion took the long way
  round to the same orientation, and the axis has to be flipped."
  [{:keys [pos prev-pos rot prev-rot damping angular-damping] :as b} dt]
  (if (static? b)
    b
    (let [dt (double dt)
          vel (v/scale (v/sub pos prev-pos) (/ 1.0 dt))
          [dx dy dz dw] (q/mul rot (q/inverse prev-rot))
          s (/ 2.0 dt)
          omega (let [o [(* dx s) (* dy s) (* dz s)]]
                  (if (neg? dw) (v/negate o) o))
          lin (max 0.0 (- 1.0 (* (double damping) dt)))
          ang (max 0.0 (- 1.0 (* (double angular-damping) dt)))]
      (assoc b :vel (v/scale vel lin) :omega (v/scale omega ang)))))

;; ---------------------------------------------------------------------------
;; Corrections

(defn inverse-mass
  "How much this body gives to a correction along `n` applied at `at`.

  With a point, this is the linear and angular response together: the
  further `at` is from the centre, the more of the correction goes into
  spinning the body and the less into moving it. With `at` nil it is the
  purely angular response to a correction about the axis `n`, which is
  what an orientation constraint asks."
  ^double [{:keys [pos inv-rot inv-mass inv-inertia] :as b} n at]
  (if (static? b)
    0.0
    (let [rn (q/rotate inv-rot (if at (v/cross (v/sub at pos) n) n))
          [ix iy iz] inv-inertia
          [rx ry rz] rn
          w (+ (* rx rx (double ix)) (* ry ry (double iy)) (* rz rz (double iz)))]
      (if at (+ w (double inv-mass)) w))))

(defn- apply-impulse
  "Moves and turns one body by `p`, applied at world point `at`.

  With `at` nil the impulse is a pure couple: it turns the body and does
  not move it, which is how an orientation is corrected without dragging
  the body sideways."
  [{:keys [inv-mass inv-rot inv-inertia] :as b} p at velocity-level? angular-relaxation]
  (if (static? b)
    b
    (let [b (if (and at (not velocity-level?))
              (update b :pos v/add-scaled p inv-mass)
              b)
          b (if (and at velocity-level?)
              (update b :vel v/add-scaled p inv-mass)
              b)
          ;; The lever arm is measured from the moved centre, as the
          ;; reference measures it. It makes no difference: the centre
          ;; moved by `p * inv-mass`, which is parallel to `p`, and a
          ;; vector crossed with something parallel to itself is zero.
          torque (if at (v/cross (v/sub at (:pos b)) p) p)
          d-omega (q/rotate (:rot b) (v/mul inv-inertia (q/rotate inv-rot torque)))]
      (if velocity-level?
        (update b :omega v/add d-omega)
        (let [[ox oy oz] (v/scale d-omega (double angular-relaxation))
              [dx dy dz dw] (q/mul [ox oy oz 0.0] (:rot b))
              [rx ry rz rw] (:rot b)]
          (with-rotation b [(+ rx (* 0.5 dx)) (+ ry (* 0.5 dy))
                            (+ rz (* 0.5 dz)) (+ rw (* 0.5 dw))]))))))

(def default-correction
  {:compliance 0.0
   :velocity-level? false
   ;; The reference halves every position-level angular correction in its
   ;; joint solver and not in its rigid body one. Halving is a relaxation:
   ;; the constraint is not satisfied in one pass but converges over the
   ;; substeps, and converges more calmly. 1.0 is the rigid body
   ;; behaviour; `allgo.physics.joint` asks for 0.5.
   :angular-relaxation 1.0})

(defn correct
  "Applies a correction between two bodies, returning `{:bodies :force}`.

  `corr` is the whole error as a vector: its direction is which way to
  push and its length is how far. `a` and `b` are indices into `bodies`;
  `b` may be nil, which corrects `a` against the world. `at` and `other-at`
  are the world points the correction acts through, and nil for a
  correction that is purely about orientation.

  `compliance` is the inverse stiffness, in metres per newton: zero is
  rigid, and larger is springier in a way that does not depend on the step
  size -- which is the point of XPBD over pushing a fraction of the error
  each iteration."
  [bodies {:keys [a b corr at other-at dt] :as opts}]
  (let [{:keys [compliance velocity-level? angular-relaxation]}
        (merge default-correction opts)
        c (v/length corr)]
    (if (zero? c)
      {:bodies bodies :force 0.0}
      (let [n (v/scale corr (/ 1.0 c))
            dt (double dt)
            w (+ (inverse-mass (bodies a) n at)
                 (if b (inverse-mass (bodies b) n other-at) 0.0))]
        (if (zero? w)
          {:bodies bodies :force 0.0}
          (let [alpha (if velocity-level? 0.0 (/ (double compliance) (* dt dt)))
                lambda (/ (- c) (+ w alpha))
                p (v/scale n (- lambda))]
            {:bodies (cond-> bodies
                       true (update a apply-impulse p at velocity-level? angular-relaxation)
                       b    (update b apply-impulse (v/negate p) other-at
                                    velocity-level? angular-relaxation))
             :force (/ lambda (* dt dt))}))))))

;; ---------------------------------------------------------------------------
;; Constraints

(defmulti solve
  "Applies one constraint, returning `{:bodies :force}`.

  Dispatches on `:kind`, so a new constraint is a new method rather than a
  change to the step -- the same shape `allgo.physics.fire` uses for its
  sources."
  (fn [_bodies c _dt] (:kind c)))

(defmethod solve :default [_ c _]
  (throw (ex-info "unknown constraint kind" {:constraint c})))

(defn attach
  "Records where a constraint's endpoints sit on their bodies.

  Stored in each body's own frame, so that as the bodies move and turn the
  attachment points go with them."
  [bodies {:keys [a b at other-at] :as c}]
  (let [at (or at (:at c))
        other-at (or other-at at)]
    (assoc c
           :local-a (world->local (bodies a) at)
           :local-b (if b (world->local (bodies b) other-at) (vec other-at)))))

(defn endpoints
  "Where a constraint's two attachment points are now, in world terms."
  [bodies {:keys [a b local-a local-b]}]
  [(local->world (bodies a) local-a)
   (if b (local->world (bodies b) local-b) (vec local-b))])

(defmethod solve :distance
  [bodies {:keys [a b distance compliance unilateral?] :as c} dt]
  (let [[p0 p1] (endpoints bodies c)
        delta (v/sub p1 p0)
        d (v/length delta)]
    (if (or (zero? d)
            ;; A rope pulls and does not push.
            (and unilateral? (< d (double distance))))
      {:bodies bodies :force 0.0}
      (correct bodies {:a a :b b
                       :corr (v/scale delta (/ (- d (double distance)) d))
                       :at p0 :other-at p1
                       :compliance (or compliance 0.0)
                       :dt dt}))))

(defn distance-constraint
  "Holds two attachment points a fixed distance apart.

  `b` may be nil to pin a point in space. `:unilateral?` makes it a rope
  rather than a rod: it pulls when stretched and does nothing when slack."
  [bodies {:keys [a b at other-at distance compliance unilateral?]}]
  (let [at (vec at)
        other-at (vec (or other-at at))]
    (attach bodies {:kind :distance :a a :b b :at at :other-at other-at
                    :distance (double (or distance
                                          (v/distance at other-at)))
                    :compliance (double (or compliance 0.0))
                    :unilateral? (boolean unilateral?)})))

;; ---------------------------------------------------------------------------
;; The loop

(def default-world
  {:dt (/ 1.0 60.0)
   :gravity [0.0 -10.0 0.0]
   ;; Substeps beat iterations: ten small steps each solved once converges
   ;; better than one step solved ten times, because every substep also
   ;; re-integrates. This is the central result of the XPBD papers.
   :substeps 10
   :iterations 1})

(defn step
  "One frame. Returns the world, with `:bodies` and `:forces` replaced.

  `:forces` is the force each constraint carried on the final substep,
  which is the thing a scene wants when it is asking whether a rope would
  have snapped."
  ([world] (step world {}))
  ([{:keys [bodies constraints] :as world} opts]
   (let [{:keys [dt gravity substeps iterations]} (merge default-world world opts)
         sdt (/ (double dt) (long substeps))]
     (loop [s 0 bodies bodies forces (vec (repeat (count constraints) 0.0))]
       (if (= s (long substeps))
         (assoc world :bodies bodies :forces forces)
         (let [bodies (mapv #(integrate % sdt gravity) bodies)
               [bodies forces]
               (loop [i 0 bodies bodies forces forces]
                 (if (= i (long iterations))
                   [bodies forces]
                   (let [[bodies forces]
                         (reduce (fn [[bs fs] k]
                                   (let [{:keys [bodies force]}
                                         (solve bs (nth constraints k) sdt)]
                                     [bodies (assoc fs k force)]))
                                 [bodies forces]
                                 (range (count constraints)))]
                     (recur (inc i) bodies forces))))]
           (recur (inc s)
                  (mapv #(update-velocities % sdt) bodies)
                  forces)))))))
