(ns allgo.simulation.motorcycle
  "A motorcycle, ridden.

  ## Which engine

  Three were on the table: `allgo.physics.rigid` and `allgo.physics.joint`
  (XPBD: free bodies held together by constraints), `allgo.physics.solver`
  (sequential impulse, TGS or XPBD, for contacts between loose bodies),
  and `allgo.physics.articulated` (Featherstone, reduced coordinates). A
  motorcycle is a tree of joints -- frame, swingarm and rear wheel; frame,
  steering head, fork and front wheel -- with nothing that closes a loop,
  and that is exactly what reduced coordinates describe.

  It also has the mass ratios that constraint formulations lose joints
  over. A 180 kilogram frame and rider on a 12 kilogram wheel is fifteen
  to one across the axle, and a pinned axle that drifts a few millimetres
  is a wheel that wobbles. In reduced coordinates the axle is a single
  angle and cannot drift at all. So the bike is one articulated model and
  the joints are exact.

  The ground is not part of the tree, and that half is contacts:
  `allgo.physics.world`'s sequential impulse sweep, which asks the model
  the same three questions it asks a brick. That is the hybrid -- exact
  joints, iterated contacts -- and it is run at several substeps a frame,
  which is the thing TGS would have added: every sweep sees the bike
  where it now is rather than where it was at the top of the frame.

  ## Wheels

  Tori, not balls. A ball touches the ground straight under its centre
  however far it leans, which is a tyre as wide as the wheel is tall and
  a bike that corners like a beach ball. A torus touches at the bottom of
  its own rim, out in the plane of the wheel, which is where a tyre does.
  See `allgo.physics.contact/torus-box`.

  ## What is a force and what is a joint

  The joints are the structure: a hinge for the swingarm, one for each
  axle, one for the steering head raked back along the steering axis, and
  a slider for the fork tubes. Everything that *acts* on the bike is a
  generalised force on those joints, recomputed every substep from where
  the bike has got to and handed to the world as `:tau`:

    shocks      a spring and damper on the swingarm angle and on the
                fork's travel, preloaded so the bike sits at its sag
    bump stops  much stiffer springs past the ends of the travel
    motor       torque on the rear axle, the lesser of the engine's peak
                torque and what its power can make at that wheel speed.
                Its reaction lands on the swingarm, which is what squats
                the rear under power
    brakes      torque against each wheel's spin
    rider       torque on the steering head, which is the only way a
                rider keeps a bike up

  The one thing that is not a joint force is what slows the bike down --
  air and rolling loss. It acts on the machine as a whole, so it goes in
  as an impulse through the frame's centre each substep.

  The springs are explicit, so they want the substeps: at a substep of
  1/480s the stiffest of them, a bump stop, turns through a third of a
  radian of its own oscillation per step, which is inside what a
  semi-implicit step holds.

  ## Riding it

  A bike stays up by being steered under itself, and a rider does it
  without thinking: lean left, steer left, the wheels run back under the
  mass. `rider-steer` does that. It is given a lean to hold and steers
  toward whichever way the bike is falling relative to it -- so asking
  for a lean to the left first steers *right*, pushing the wheels out
  from under the bike and tipping it left, then follows it round. That
  is countersteering, and nobody wrote it down as a rule: it falls out of
  steering into the fall.

  The gains fall with the square of speed. The steering that rights a
  bike does it through the sideways acceleration of the wheels, which
  goes as speed squared over turning radius; at walking pace it takes a
  lot of handlebar and at motorway speed very little.

  ## Axes

  The bike faces +x with +y up, so its left is -z. Wheels turn about -z,
  which makes a positive wheel speed forward. A positive steering angle
  turns the front wheel toward -z, and a positive lean is toward -z too,
  so that steering *into* a lean is steering with the same sign."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.world :as pw]
            [allgo.simulation.rider :as rider]))

;; ---------------------------------------------------------------------------
;; The bike

(def defaults
  "Everything about the machine, in SI units. Positions are in the
  frame's coordinates, whose origin is the sprung mass's centre."
  {:wheel-radius 0.32            ; outside of the tyre
   :tyre-radius 0.06             ; the tube of the torus
   :rake (/ (* 25.0 Math/PI) 180.0)
   :head [0.45 0.25 0.0]         ; where the steering axis meets the frame
   :fork-length 0.65             ; head to axle, along the steering axis
   :fork-offset 0.03             ; axle ahead of the steering axis
   :pivot [-0.20 -0.12 0.0]      ; swingarm pivot
   :swingarm [-0.55 -0.21 0.0]   ; pivot to rear axle
   :frame-mass 180.0             ; frame, engine and a rider lumped in
   :carried-mass 0.0             ; a rider riding on it as a body of their own
   :swingarm-mass 8.0
   :steering-mass 6.0
   :fork-mass 5.0
   :wheel-mass 11.0
   ;; Suspension, as rates at the wheel.
   :rear-rate 25000.0            ; N/m
   :rear-damping 1600.0          ; N s/m
   :front-rate 20000.0
   :front-damping 1400.0
   :rear-travel [-0.22 0.10]     ; swingarm angle, radians
   :front-travel [-0.12 0.05]    ; fork slide, metres
   :bump-stiffness 15.0          ; times the spring rate
   ;; Drive.
   :peak-torque 320.0            ; N m at the rear wheel
   :peak-power 45000.0           ; W
   :rear-brake 350.0             ; N m at full lever
   :front-brake 700.0
   ;; What slows it with the throttle shut.
   :drag-area 0.6                ; Cd times frontal area, m^2
   :rolling-resistance 0.015     ; of the weight
   ;; The rider's hands on the bars.
   :steer-stiffness 500.0        ; N m per radian short of where they want it
   :steer-damping 25.0
   :max-steer 0.55})

(defn- box-inertia [[sx sy sz] m]
  (let [f (/ (double m) 12.0)
        sx (double sx) sy (double sy) sz (double sz)]
    [[(* f (+ (* sy sy) (* sz sz))) 0.0 0.0]
     [0.0 (* f (+ (* sx sx) (* sz sz))) 0.0]
     [0.0 0.0 (* f (+ (* sx sx) (* sy sy)))]]))

(defn- wheel-inertia
  "Most of a wheel's mass is out at the rim, so it is nearer a hoop than
  a disc about its axle."
  [m r]
  (let [m (double m) r (double r)
        axial (* 0.8 m r r)
        diametral (* 0.45 m r r)]
    [[diametral 0.0 0.0] [0.0 diametral 0.0] [0.0 0.0 axial]]))

(defn steering-axis
  "Up the steering axis, raked back from vertical."
  [{:keys [rake]}]
  [(- (Math/sin rake)) (Math/cos rake) 0.0])

(defn- fork-axis
  "Down the fork tubes: the way the slider extends."
  [cfg]
  (v/negate (steering-axis cfg)))

(defn- front-axle
  "The front axle in the steering head's own frame."
  [{:keys [rake fork-length fork-offset] :as cfg}]
  (v/add (v/scale (fork-axis cfg) fork-length)
         (v/scale [(Math/cos rake) (Math/sin rake) 0.0] fork-offset)))

(def ^:private idq [0.0 0.0 0.0 1.0])

(def swingarm 0)
(def rear-wheel 1)
(def steering 2)
(def fork 3)
(def front-wheel 4)

(defn model
  "The bike as an articulated model: the frame is the free root, and five
  joints hang off it -- every one of them a single number, so the
  joint's index in `:links` is also its index in `q`, `qd` and `tau`."
  ([] (model defaults))
  ([cfg]
   (let [{:keys [wheel-radius tyre-radius head pivot swingarm fork-length
                 frame-mass swingarm-mass steering-mass fork-mass wheel-mass]} cfg
         tyre {:shape :torus :major (- (double wheel-radius) (double tyre-radius))
               :minor tyre-radius :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
         spin [0.0 0.0 -1.0]]
     {:base {:mass frame-mass :com [0.0 0.0 0.0]
             :inertia (box-inertia [1.3 0.6 0.45] frame-mass)
             ;; What hits the ground when it goes over: the engine and the
             ;; tank, not the rider's knee.
             :shape :box :size [0.95 0.32 0.30]
             :shape-pose {:rot idq :pos [0.0 -0.02 0.0]}}
      :links
      [{:parent -1 :joint :revolute :axis [0.0 0.0 1.0]
        :origin {:rot nil :pos pivot}
        :mass swingarm-mass :com (v/scale swingarm 0.5)
        :inertia (box-inertia [0.6 0.08 0.22] swingarm-mass)}
       (merge {:parent 0 :joint :revolute :axis spin
               :origin {:rot nil :pos swingarm}
               :mass wheel-mass :com [0.0 0.0 0.0]
               :inertia (wheel-inertia wheel-mass wheel-radius)}
              tyre)
       {:parent -1 :joint :revolute :axis (steering-axis cfg)
        :origin {:rot nil :pos head}
        :mass steering-mass :com [0.0 0.08 0.0]
        :inertia (box-inertia [0.2 0.12 0.7] steering-mass)}
       {:parent 2 :joint :prismatic :axis (fork-axis cfg)
        :origin {:rot nil :pos [0.0 0.0 0.0]}
        :mass fork-mass :com (v/scale (fork-axis cfg) (* 0.5 (double fork-length)))
        :inertia (box-inertia [0.1 0.6 0.25] fork-mass)}
       (merge {:parent 3 :joint :revolute :axis spin
               :origin {:rot nil :pos (front-axle cfg)}
               :mass wheel-mass :com [0.0 0.0 0.0]
               :inertia (wheel-inertia wheel-mass wheel-radius)}
              tyre)]})))

;; ---------------------------------------------------------------------------
;; Where it starts

(defn ride-height
  "How high the frame's origin sits with both tyres on the ground and the
  suspension at its sag, which the springs are preloaded to make zero."
  [{:keys [wheel-radius pivot swingarm]}]
  (- (double wheel-radius) (+ (double (nth pivot 1)) (double (nth swingarm 1)))))

(defn start-pose
  "Upright, rolling forward at `speed`, facing along `heading` radians
  about the vertical from +x, at `pos` on the ground."
  ([cfg speed] (start-pose cfg speed 0.0 [0.0 0.0]))
  ([cfg speed heading [x z]]
   (let [spin (/ (double speed) (double (:wheel-radius cfg)))]
     {:q [0.0 0.0 0.0 0.0 0.0]
      :qd [0.0 spin 0.0 0.0 spin]
      :base {:rot (q/from-axis-angle [0.0 1.0 0.0] (double heading))
             :pos [(double x) (ride-height cfg) (double z)]
             ;; Angular velocity, then the velocity of the frame's origin,
             ;; both in the frame's own coordinates.
             :vel [0.0 0.0 0.0 (double speed) 0.0 0.0]}})))

;; ---------------------------------------------------------------------------
;; Reading it

(defn telemetry
  "What a rider would want to know, from a pose."
  [cfg {:keys [q qd base]}]
  (let [{:keys [rot vel]} base
        lateral (q/rotate rot [0.0 0.0 1.0])
        forward (q/rotate rot [1.0 0.0 0.0])
        [wx] vel]
    {:speed (double (nth vel 3))
     :lean (Math/asin (max -1.0 (min 1.0 (double (nth lateral 1)))))
     :lean-rate (- (double wx))
     :heading (Math/atan2 (- (double (nth forward 2))) (double (nth forward 0)))
     :steer (double (nth q steering))
     :steer-rate (double (nth qd steering))
     :rear-travel (double (nth q swingarm))
     :front-travel (double (nth q fork))
     :rear-wheel-speed (* (double (nth qd rear-wheel)) (double (:wheel-radius cfg)))
     :front-wheel-speed (* (double (nth qd front-wheel)) (double (:wheel-radius cfg)))
     :upright? (> (double (nth (q/rotate rot [0.0 1.0 0.0]) 1)) 0.3)}))

;; ---------------------------------------------------------------------------
;; What acts on it

(defn- lever
  "How far a vertical push at the rear axle is from the swingarm pivot:
  what turns a rate at the wheel into a rate at the pivot."
  ^double [{:keys [swingarm]}]
  (abs (double (nth swingarm 0))))

(defn- rear-preload
  "The torque the rear spring has to find at zero travel to hold its
  share of the sprung weight."
  ^double [{:keys [frame-mass swingarm-mass pivot swingarm head] :as cfg}]
  (let [rear-x (+ (double (nth pivot 0)) (double (nth swingarm 0)))
        front-x (double (nth (v/add head (front-axle cfg)) 0))
        share (/ front-x (- front-x rear-x))
        sprung (* (+ (double frame-mass) (double (:carried-mass cfg 0.0))
                     (* 0.5 (double swingarm-mass)))
                  9.81 share)]
    (* sprung (lever cfg))))

(defn- front-preload
  "The force along the fork the front spring has to find at zero travel."
  ^double [{:keys [frame-mass steering-mass fork-mass rake pivot swingarm head] :as cfg}]
  (let [rear-x (+ (double (nth pivot 0)) (double (nth swingarm 0)))
        front-x (double (nth (v/add head (front-axle cfg)) 0))
        share (/ (- rear-x) (- front-x rear-x))
        sprung (* (+ (double frame-mass) (double (:carried-mass cfg 0.0))
                     (double steering-mass) (* 0.5 (double fork-mass)))
                  9.81 share)]
    (* sprung (Math/cos (double rake)))))

(defn- spring
  "A preloaded spring and damper at `x`, with bump stops past `[lo hi]`.
  Positive is the way the spring pushes: extension."
  [x xd rate damping preload [lo hi] stop]
  (let [x (double x) xd (double xd) rate (double rate) damping (double damping)
        preload (double preload) stop (double stop)
        lo (double lo) hi (double hi)
        bump (cond (< x lo) (* stop rate (- lo x))
                   (> x hi) (* stop rate (- hi x))
                   :else 0.0)
        ;; Past a stop the damping goes up with the stiffness, or the
        ;; stop is a trampoline.
        damp (if (or (< x lo) (> x hi)) (* damping (Math/sqrt stop)) damping)]
    (+ preload (- (* rate x)) bump (- (* damp xd)))))

(defn- drive
  "The motor's torque at the rear wheel: its peak, until the wheel is
  turning fast enough that its power runs out first."
  ^double [{:keys [peak-torque peak-power]} ^double throttle ^double wheel-omega]
  (* throttle
     (min (double peak-torque)
          (/ (double peak-power) (max 1.0 (abs wheel-omega))))))

(defn- brake
  "Torque against a wheel's spin, faded in over the last few radians a
  second so that a stopped wheel is held rather than chattered."
  ^double [^double strength ^double lever ^double omega]
  (- (* strength lever (Math/tanh (/ omega 2.0)))))

(defn rider-steer
  "The steering angle a rider wants: into the fall, measured from the lean
  they are trying to hold, plus the angle that lean needs to hold a
  steady turn at this speed.

  `lean` is what the bike is doing, `target` what the rider wants, and
  `held` how long it has been off, integrated: the feed-forward is only
  the steer a rigid, point-mass bike would need, and a real one -- trail,
  a rider sat up on it -- needs a little more or less, which a rider
  finds by feel and `held` finds by accumulating. All the gains are
  scaled by `g L / v^2` -- the steer it takes, at this speed,
  to make a sideways acceleration that holds a radian of lean -- which
  is why the same rider is heavy-handed at a walk and barely moves the
  bars at speed."
  (^double [cfg t ^double target] (rider-steer cfg t target 0.0))
  (^double [cfg {:keys [lean lean-rate speed]} ^double target ^double held]
   (let [{:keys [max-steer]} cfg
         wheelbase 1.5
         v (max 2.5 (abs (double speed)))
         scale (/ (* 9.81 wheelbase) (* v v))
         feed-forward (* scale (Math/tan target))
         k-lean (* 5.0 scale)
         k-rate (* 1.3 scale)
         k-held (* 2.0 scale)
         want (+ feed-forward
                 (* k-lean (- (double lean) target))
                 (* k-rate (double lean-rate))
                 (* k-held held))]
     (max (- (double max-steer)) (min (double max-steer) want)))))

(defn forces
  "The generalised forces on every joint, given where the bike is and
  what the rider is asking for.

  `controls` is `{:throttle 0..1 :brake 0..1 :lean radians :hands? bool}`,
  and `:lean-held` is the rider's integrated lean error, which `step`
  keeps.
  With `:hands?` false nobody is steering: the bars are free, and the
  bike is left to whatever its own geometry does with them."
  [cfg {:keys [q qd] :as pose} {lever-pull :brake :keys [throttle lean hands?]
                                :or {throttle 0.0 lever-pull 0.0 lean 0.0 hands? true}
                                :as controls}]
  (let [{:keys [rear-rate rear-damping front-rate front-damping rear-travel
                front-travel bump-stiffness rear-brake front-brake
                steer-stiffness steer-damping]} cfg
        t (telemetry cfg pose)
        l (lever cfg)
        l2 (* l l)
        ;; A rate at the wheel becomes a rate at the pivot through the
        ;; square of the lever: once for the force, once for the travel.
        rear (spring (double (nth q swingarm)) (double (nth qd swingarm))
                     (* (double rear-rate) l2) (* (double rear-damping) l2)
                     (rear-preload cfg) rear-travel (double bump-stiffness))
        front (spring (double (nth q fork)) (double (nth qd fork))
                      (double front-rate) (double front-damping)
                      (front-preload cfg) front-travel (double bump-stiffness))
        rear-omega (double (nth qd rear-wheel))
        front-omega (double (nth qd front-wheel))
        motor (drive cfg (double throttle) rear-omega)
        steer (if hands?
                (let [want (rider-steer cfg t (double lean)
                                        (double (:lean-held controls 0.0)))]
                  (- (* (double steer-stiffness) (- want (double (nth q steering))))
                     (* (double steer-damping) (double (nth qd steering)))))
                ;; Free bars still have bearings.
                (- (* 0.5 (double (nth qd steering)))))]
    [rear
     (+ motor (brake (double lever-pull) (double rear-brake) rear-omega))
     steer
     front
     (brake (double lever-pull) (double front-brake) front-omega)]))

;; ---------------------------------------------------------------------------
;; A world to ride it in

(def world-defaults
  "Tyres grip better than bricks, and a bike wants more iterations than a
  pile of them: two contacts carry the whole machine."
  {:gravity [0.0 -9.81 0.0]
   :friction 1.1
   :restitution 0.0
   :iterations 12
   :self-collide? false})

(defn ground
  "A floor to ride on: a slab `size` metres square with its top at y = 0."
  [size]
  (rigid/box {:pos [0.0 -0.5 0.0] :size [size 1.0 size]}))

(def bike-frame-mass
  "The frame and engine without anyone on them."
  110.0)

(defn with-rider
  "`cfg` for a bike that carries a rider as a body of their own: the frame
  loses the rider it had lumped into it, and the suspension is preloaded
  for the weight that now arrives through the seat and pegs instead."
  [cfg]
  (assoc cfg :frame-mass bike-frame-mass
         :carried-mass (rider/mass (rider/model cfg))))

(def ^:private slack
  "How much of a thrown rider's muscle is left: none. What comes off the
  bike is a plain ragdoll, held together by its joint limits."
  0.0)

(defn scene
  "A world with the bike in it, on `obstacles` -- static boxes, the first
  of which is usually the floor. With `:rider? true` a rider is seated on
  it and pinned there; `cfg` should then come from `with-rider`."
  ([cfg obstacles] (scene cfg obstacles (start-pose cfg 6.0)))
  ([cfg obstacles pose] (scene cfg obstacles pose nil))
  ([cfg obstacles pose {:keys [rider?]}]
   (let [bike {:model (model cfg) :pose pose}]
     (if-not rider?
       (pw/world obstacles [bike] world-defaults)
       (let [r (rider/model cfg)]
         (assoc (pw/world obstacles
                          [bike {:model r :pose (rider/start-pose cfg r (:base pose))}]
                          ;; While seated the rider's knees are tucked in
                          ;; against the frame; the two only collide once
                          ;; they have come apart.
                          (assoc world-defaults
                                 :pins (rider/pins cfg 0 1 steering)
                                 :collide-models? false))
                :rider {:tone 1.0 :attached? true}))))))

(defn bike-pose [scene] (:pose (first (:models scene))))

(defn rider-pose [scene] (:pose (second (:models scene))))

(def bail-impulse
  "How hard a rider throws themselves off, in newton seconds: about what a
  jump sideways out of the saddle is."
  110.0)

(defn throw-rider
  "The rider comes off, now: every pin lets go, the muscles switch off,
  and from here on the rider collides with the bike like anything else.

  With `push` -- newton seconds, toward the bike's left for positive --
  the rider also launches themselves off the side, which is what bailing
  out is: an impulse at the chest, sideways and up, and the same impulse
  the other way into the frame, which is left kicked and wobbling. A
  bike that is already going over needs no push; the crash throws them."
  ([scene] (throw-rider scene 0.0))
  ([scene push]
   (if-not (get-in scene [:rider :attached?])
     scene
     (let [released (assoc scene :pins [] :collide-models? true
                           :rider {:tone slack :attached? false})]
       (if (zero? (double push))
         released
         (let [{bm :model bp :pose} (first (:models scene))
               {rm :model rp :pose} (second (:models scene))
               rot (:rot (:base bp))
               ;; Toward the bike's left, which is -z in its frame, and up.
               dir (v/normalize (v/add (q/rotate rot [0.0 0.0 (if (pos? (double push)) -1.0 1.0)])
                                       [0.0 0.7 0.0]))
               chest (let [f (ab/frame-of rm (:q rp) (:base rp) rider/torso)]
                       (v/add (:pos f) (q/rotate (:rot f) [0.0 0.25 0.0])))
               seat (:pos (:base bp))
               j (abs (double push))]
           (-> released
               (update-in [:models 1 :pose] #(ab/apply-impulse rm % rider/torso chest dir j))
               (update-in [:models 0 :pose] #(ab/apply-impulse bm % -1 seat (v/negate dir) j)))))))))

(defn- let-go
  "A rider still on a bike that has gone over, or whose grip on it has
  broken anywhere, is thrown."
  [cfg scene]
  (if (and (get-in scene [:rider :attached?])
           (or (seq (:broken-pins scene))
               (not (:upright? (telemetry cfg (bike-pose scene))))))
    (throw-rider scene)
    scene))

(defn- total-mass ^double [{:keys [frame-mass swingarm-mass steering-mass fork-mass wheel-mass
                                   carried-mass]}]
  (+ (double frame-mass) (double (or carried-mass 0.0)) (double swingarm-mass)
     (double steering-mass) (double fork-mass) (* 2.0 (double wheel-mass))))

(defn resistance
  "The force holding the bike back at `speed`: air, which goes as the
  square of it, and the tyres' rolling loss, which does not. Without
  these a steady throttle accelerates for ever."
  ^double [{:keys [drag-area rolling-resistance] :as cfg} ^double speed]
  (+ (* 0.5 1.2 (double drag-area) speed speed)
     (* (double rolling-resistance) (total-mass cfg) 9.81 (Math/tanh (/ speed 0.5)))))

(defn- slow
  "`pose` after one substep of resistance, as an impulse through the
  frame's centre against the way it is going over the ground."
  [cfg model pose dt]
  (let [{:keys [rot vel pos]} (:base pose)
        ground-vel (let [w (q/rotate rot (subvec (vec vel) 3 6))] [(nth w 0) 0.0 (nth w 2)])
        speed (v/length ground-vel)]
    (if (< speed 1e-3)
      pose
      (ab/apply-impulse model pose -1 pos (v/scale ground-vel (/ -1.0 speed))
                        (* (resistance cfg speed) (double dt))))))

(defn step
  "One substep: slow the bike by what the air and the tyres cost it, work
  out the forces from where it is now, hand them to the world, and step
  it."
  [cfg scene controls dt]
  (let [;; A thrown rider is no longer riding: nobody on the throttle,
        ;; the brakes or the bars. The bike is on its own, and this one
        ;; cannot stay up on its own.
        controls (if (and (:rider scene) (not (get-in scene [:rider :attached?])))
                   {:throttle 0.0 :brake 0.0 :lean 0.0 :hands? false}
                   controls)
        {:keys [model pose]} (first (:models scene))
        pose (slow cfg model pose dt)
        ;; The rider's feel for how far off the lean has been: integrated,
        ;; leaking away over twenty seconds or so, so an old error does not
        ;; outlive the turn it came from, and capped so it cannot wind up.
        err (- (double (:lean (telemetry cfg pose))) (double (:lean controls 0.0)))
        held (let [h (+ (* (double (:lean-held scene 0.0)) (- 1.0 (* 0.05 (double dt))))
                        (* err (double dt)))]
               (max -0.3 (min 0.3 h)))
        controls (assoc controls :lean-held held)
        rider (second (:models scene))
        scene (cond-> (-> scene
                          (assoc-in [:models 0 :pose] pose)
                          (assoc-in [:models 0 :tau] (forces cfg pose controls)))
                true (assoc :lean-held held)
                rider (assoc-in [:models 1 :tau]
                                (rider/muscles (:model rider) (:pose rider)
                                               (get-in scene [:rider :tone] 1.0)
                                               (when (get-in scene [:rider :attached?])
                                                 (rider/weight-shift (:lean controls 0.0))))))]
    (let-go cfg (pw/step scene dt))))

(defn run
  "`n` substeps of `dt` under fixed `controls`. For tests and tuning."
  [cfg scene controls dt n]
  (reduce (fn [s _] (step cfg s controls dt)) scene (range n)))
