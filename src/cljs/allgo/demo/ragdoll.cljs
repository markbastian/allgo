(ns allgo.demo.ragdoll
  "A figure with joints, dropped on the floor.

  It is drawn as a crash test dummy, and only drawn so: each part the
  physics collides is still a box, and the dummy's capsules, rounded
  blocks and targets sit inside them, but for the balls at the joints
  and the hands and feet at the ends.

  Everything under this demo is `allgo.physics.articulated`, which
  describes a jointed body by its joint angles rather than by a body per
  limb and constraints to hold them together. The difference is not
  subtle and it is worth knowing what to look for.

  ## What the other way would look like

  `allgo.physics.joint` can build the same figure out of eleven free
  bodies and ten constraints, and it is the right shape for a pile of
  bricks. For a chain of limbs it gives thirty-six degrees of freedom to
  a thing that has six, and then spends its iterations taking thirty
  away again. They never quite go, and the error grows with the mass
  ratio across the joint -- which is exactly what a hand on a forearm on
  an upper arm is. Two meter-long links, the lower a thousand times the
  heavier, come nearly three meters apart at five substeps.

  Here there is no such number to report. A joint angle cannot be
  violated because there is nowhere for the violation to live. Watch the
  shoulders as it falls from a height: the arms swing, and the point
  where the arm meets the body does not move at all.

  ## What to look at

  **The limits are doing most of the work.** Turn them off and the
  figure settles into a pile of sticks: the head folds back on itself,
  the knees bend the wrong way, the hips reach a hundred and fifty
  degrees. Every joint here has a range -- a cone and a twist for the
  ball joints, a pair of angles for the elbows and knees -- and those
  are solved as one-sided constraints in the same sweep as the floor
  contacts. A limit is a push between a link and its parent; a contact
  is a push between a link and the floor; the solver does not need to be
  told which.

  **The bricks are the point.** They are loose rigid bodies, simulated
  by `allgo.physics.solver`'s kind of arithmetic, and the figure is a
  jointed model simulated by nobody's idea of the same thing. Until
  `allgo.physics.world` existed the two could not touch: a ragdoll
  could fall down stairs, because stairs do not move, and could not
  knock a brick off a wall. Dropped on four courses, it knocks all
  sixteen askew.

  What makes that work is that a contact does not care what is on
  either side of it. It asks both sides the same three questions -- how
  fast is your surface here, what does a unit impulse buy, take this
  impulse -- and a brick answers from its inverse inertia while a
  forearm answers with a walk up its own tree to the pelvis. One
  Gauss-Seidel sweep covers every contact in the scene together, which
  is what lets a brick being pushed by a hand and by the brick under it
  see both pushes in the same iteration.

  **`shove` pushes the chest away from wherever you are standing**, by
  an impulse rather than a force -- `apply-impulse`, the same entry
  point a contact uses. It is worth trying on a figure that has already
  settled: the whole body answers, because an impulse at the chest
  reaches the feet through the joints in the same step rather than a
  link at a time.

  **The elbows and knees are hinges and the rest are ball joints.** An
  elbow that could swing sideways would not read as an elbow. The
  shoulders and hips are spherical, which is three numbers of velocity
  and a quaternion of configuration -- three stacked hinges would gimbal
  lock, and a tumbling figure finds the orientation where they do.

  **Limbs collide with each other**, and `self collision` turns that
  off. The difference is easier to count than to see: left to settle on
  the floor, the figure comes to rest with thirty-two contact points'
  worth of limb inside other limb when it is off, and two when it is
  on. An arm inside the ribcage is the one that shows.

  A link and its parent are never tested against each other -- they meet
  at the joint and overlap there by construction, so a contact between
  them would be permanent and would push the figure apart from the
  inside. Everything else is fair game, and the shoulders sit two
  centimeters wider than the ribcage for the same reason: a forearm
  resting five millimeters inside the chest at the start is a contact
  that never goes away.

  ## What it costs

  Eleven bodies, twenty-eight degrees of freedom, and as many loose
  bricks as `brick rows` asks for. On the JVM a settled step is 1.4ms
  with no bricks, 2.4ms with eight and 3.0ms with sixteen; in flight it
  is 0.17ms. Self collision adds about half of that last one and
  nothing at all to the others, because a figure whose limbs are not
  inside each other is a figure with fewer contacts to solve.

  In a browser it is a good deal slower, and this docstring used to
  name a figure for that which should not have been trusted -- it was
  measured in an automated tab, which Chrome backgrounds and gives less
  of a processor, and the same measurement taken twice varied by half.
  What is worth saying without a number is *why* it is slower. The
  arithmetic is on flat arrays but the recursion walking over it is
  not: `(nth ls j)`, a map destructured by keyword, a `conj` onto a
  vector, once per link per direction per contact. On a JVM that is
  nearly free beside the floating point. In JavaScript it is most of
  the cost. The frame counter is in the corner."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.world :as pw]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]
            ["three/examples/jsm/geometries/RoundedBoxGeometry.js" :refer [RoundedBoxGeometry]]))

(def ^:private identity-q [0.0 0.0 0.0 1.0])

(defn- deg [d] (* (double d) (/ js/Math.PI 180.0)))

(defn- box-inertia
  "A solid box's inertia about its own center."
  [[sx sy sz] m]
  (let [f (/ (double m) 12.0)
        sx (double sx) sy (double sy) sz (double sz)]
    [[(* f (+ (* sy sy) (* sz sz))) 0.0 0.0]
     [0.0 (* f (+ (* sx sx) (* sz sz))) 0.0]
     [0.0 0.0 (* f (+ (* sx sx) (* sy sy)))]]))

(defn- bone
  "A limb reaching `down` or up from its joint, as a box.

  The joint is at one end and the mass in the middle, which is what a
  limb is: a link's frame sits where it attaches to its parent, so the
  center of mass -- and with it the shape -- is half a length away."
  [dir parent joint axis origin size m]
  (let [[_ sy _] size]
    (cond-> {:parent parent :joint joint :origin {:rot nil :pos origin}
             :mass m :com [0.0 (* 0.5 (double dir) (double sy)) 0.0]
             :inertia (box-inertia size m)
             :shape :box :size size}
      axis (assoc :axis axis))))

(def ^:private down -1.0)
(def ^:private up 1.0)

(defn- ball [link cone twist] (assoc link :cone (deg cone) :twist (deg twist)))
(defn- hinge [link lo hi] (assoc link :limit [(deg lo) (deg hi)]))

(def ^:private joint-damping
  "Each joint's friction, in newton meters a radian a second."
  0.2)

(defn- figure
  "Eleven boxes and ten joints, roughly a person.

  The pelvis is the root and is free to move, which is what makes this a
  ragdoll rather than a mobile bolted to the ceiling. Masses are near
  enough a seventy kilogram adult; the exact numbers matter less than
  their ratios, and the ratios are what a constraint solver would
  struggle with.

  Every joint is damped a little, its own friction: undamped, a shin
  hanging over the slab's edge swings for ever and the figure never
  sleeps. The armature beside the damper, `d dt` at the demo's longest
  step, is what lets a limb as light as a forearm turning about its own
  length be damped without the damper overshooting."
  [limits?]
  (let [b (fn [l] (-> (if limits? l (dissoc l :cone :twist :limit))
                      (assoc :damping joint-damping :armature (/ joint-damping 60.0))))]
    {:base {:mass 12.0 :com [0.0 0.0 0.0]
            :inertia (box-inertia [0.28 0.20 0.18] 12.0)
            :shape :box :size [0.28 0.20 0.18]
            :shape-pose {:rot identity-q :pos [0.0 0.0 0.0]}}
     :links
     [(b (ball (bone up -1 :spherical nil [0.0 0.10 0.0] [0.34 0.42 0.20] 25.0) 25 35))
      (b (ball (bone up 0 :spherical nil [0.0 0.42 0.0] [0.20 0.22 0.20] 5.0) 40 60))
      (b (ball (bone down 0 :spherical nil [0.23 0.38 0.0] [0.10 0.28 0.10] 2.5) 85 70))
      (b (hinge (bone down 2 :revolute [1.0 0.0 0.0] [0.0 -0.28 0.0]
                      [0.09 0.26 0.09] 1.8) 0 150))
      (b (ball (bone down 0 :spherical nil [-0.23 0.38 0.0] [0.10 0.28 0.10] 2.5) 85 70))
      (b (hinge (bone down 4 :revolute [1.0 0.0 0.0] [0.0 -0.28 0.0]
                      [0.09 0.26 0.09] 1.8) 0 150))
      (b (ball (bone down -1 :spherical nil [0.09 -0.10 0.0] [0.14 0.40 0.14] 8.0) 60 35))
      (b (hinge (bone down 6 :revolute [1.0 0.0 0.0] [0.0 -0.40 0.0]
                      [0.12 0.38 0.12] 3.5) -150 0))
      (b (ball (bone down -1 :spherical nil [-0.09 -0.10 0.0] [0.14 0.40 0.14] 8.0) 60 35))
      (b (hinge (bone down 8 :revolute [1.0 0.0 0.0] [0.0 -0.40 0.0]
                      [0.12 0.38 0.12] 3.5) -150 0))]}))

(def ^:private floor
  (rigid/box {:pos [0.0 -0.5 0.0] :size [40.0 1.0 40.0]}))

(defn- slab
  "Something to fall off, tilted by `tilt` degrees."
  [tilt]
  (rigid/box {:pos [0.0 0.25 0.0] :size [1.6 0.5 1.6]
              :rot (q/from-axis-angle [0.0 0.0 1.0] (deg tilt))}))

(def ^:private brick-size [0.22 0.11 0.14])

(defn- wall
  "A running-bond wall of loose bricks, `rows` high, for the figure to
  land in. Small and light: the point is that they move when a ragdoll
  hits them, not that they stand up to one.

  Beside the slab's low edge, clear of it: it used to stand where the
  slab is, its bricks started up to 31 cm inside it, and the convex
  solve threw them out at 16 m/s -- a brick wedged between the slab and
  the floor is pushed from both sides at once."
  [rows]
  (let [[bw bh _] brick-size]
    (vec (for [r (range (long rows))
               c (range 4)
               :let [offset (if (odd? r) (* 0.5 (double bw)) 0.0)]]
           (rigid/box {:pos [(+ (* (- c 1.5) (double bw)) offset -1.45)
                             (+ (* 0.5 (double bh)) (* r (double bh)))
                             0.0]
                       :size brick-size
                       :density 500.0})))))

(def ^:private ^js controls
  #js {:height 1.8
       :tilt 18
       :obstacle true
       :limits true
       :selfCollide true
       :friction 0.7
       :restitution 0.0
       ;; The sequential impulse: a heap of twenty dummies is 15 ms a
       ;; step on the JVM under it and 109 under the convex solve, which
       ;; earns that on tall stacks and buys nothing visible in a pile of
       ;; ragdolls. One dummy is quick under either.
       :solver "sequential-impulse"
       :iterations 8
       :substeps 1
       :spin 2.0
       :shoveForce 120
       :brickRows 4
       :dummies 1
       :start (fn [])
       :reset (fn [])
       :shove (fn [])})

(defn- start-state
  "The figure standing upright, `height` above the floor at `x`, `z`,
  facing `yaw` about the vertical, given a turn."
  ([model height spin] (start-state model height spin 0.0 0.0 0.0))
  ([model height spin x z yaw]
   {:q (mapv (fn [l] (if (= :spherical (:joint l)) identity-q 0.0)) (ab/chain model))
    :qd (vec (repeat (ab/dof model) 0.0))
    :base {:rot (q/mul (q/from-axis-angle [0.0 1.0 0.0] yaw) (q/from-axis-angle [0.0 0.0 1.0] 0.15))
           :pos [x height z]
           :vel [0.0 spin 0.0 0.0 0.0 0.0]}}))

(defn- drop-points
  "Where `n` dummies start: one where one always has, or more in rows of
  five, three quarters of a meter apart -- clear of each other at the
  start -- each row a little higher than the last and each dummy turned
  a little from its neighbor. All of them in the rows: the first kept
  in the middle once stood inside the eighth, which started a third of
  a meter higher in the same place."
  [n height]
  (for [i (range n)
        :let [col (mod i 5) row (quot i 5)]]
    (if (== 1 (long n))
      [0.0 height 0.0 0.0]
      [(* 0.75 (- col 2)) (+ height (* 0.35 row)) (* 0.75 (- row 1)) (* 0.7 i)])))

;; ---------------------------------------------------------------------------
;; Looks: a crash test dummy
;;
;; Only what is drawn. Each part is still the box the physics collides,
;; and what is drawn for it stays inside that box, but for the balls at
;; the joints and the hands and feet at the ends.

(defn- target-texture
  "The quartered disc a crash test dummy wears to be tracked on film."
  []
  (let [c (js/document.createElement "canvas")
        ctx (.getContext c "2d")
        h 64.0 r 58.0]
    (set! (.-width c) 128) (set! (.-height c) 128)
    (dotimes [i 4]
      (.beginPath ctx)
      (.moveTo ctx h h)
      (.arc ctx h h r (* i 0.5 js/Math.PI) (* (inc i) 0.5 js/Math.PI))
      (.closePath ctx)
      (set! (.-fillStyle ctx) (if (even? i) "#17171a" "#f2c018"))
      (.fill ctx))
    (.beginPath ctx)
    (.arc ctx h h r 0 (* 2 js/Math.PI))
    (set! (.-lineWidth ctx) 6)
    (set! (.-strokeStyle ctx) "#17171a")
    (.stroke ctx)
    (THREE/CanvasTexture. c)))

(def ^:private materials
  (delay
    {:yellow (THREE/MeshStandardMaterial. #js {:color 0xf2c018 :roughness 0.5})
     :black (THREE/MeshStandardMaterial. #js {:color 0x17171a :roughness 0.6})
     :hand (THREE/MeshStandardMaterial. #js {:color 0x3b3f48 :roughness 0.7})
     :target (THREE/MeshStandardMaterial. #js {:map (target-texture) :transparent true :roughness 0.5
                                               :polygonOffset true :polygonOffsetFactor -2})
     :brick (THREE/MeshStandardMaterial. #js {:color 0xb98d5f :roughness 0.85})}))

(def ^:private reaching-up
  "The links whose joint is at the bottom of their box: the spine and the
  neck. The rest hang from theirs."
  #{0 1})

(defn- dummy-part
  "What is drawn for `link`, whose collision box is `body`: a group in the
  box's own frame, placed where the box is."
  [link body]
  (let [{:keys [yellow black hand target]} @materials
        [sx sy sz] (:size body)
        g (THREE/Group.)
        add (fn [geo mat x y z]
              (let [m (THREE/Mesh. geo mat)]
                (.set (.-position m) x y z)
                (.add g m)
                m))
        ;; the joint's end of the box, and the far end
        jy (if (reaching-up link) (* -0.5 sy) (* 0.5 sy))
        joint! (fn [r] (add (THREE/SphereGeometry. r 16 12) black 0 jy 0))
        ;; a target on the face at `x` (either side), at height `y`
        target! (fn [size x y]
                  (let [m (add (THREE/PlaneGeometry. size size) target x y 0)]
                    (set! (.. m -rotation -y) (if (pos? x) (* 0.5 js/Math.PI) (* -0.5 js/Math.PI)))))]
    (case (long link)
      ;; the pelvis: a belt, and a target on each hip
      -1 (do (add (RoundedBoxGeometry. sx sy sz 3 0.04) yellow 0 0 0)
             (add (THREE/BoxGeometry. (* 1.02 sx) 0.035 (* 1.02 sz)) black 0 (- (* 0.5 sy) 0.04) 0)
             (doseq [x [(* 0.5 sx) (* -0.5 sx)]] (target! 0.09 (+ x (* (js/Math.sign x) 0.002)) -0.01)))
      ;; the chest: ribbed, a target on each side, the waist a ball
      0 (do (add (RoundedBoxGeometry. sx sy sz 3 0.05) yellow 0 0 0)
            (doseq [y [0.02 0.07 0.12]]
              (add (THREE/BoxGeometry. (* 0.7 sx) 0.012 0.004) black 0 y (- (+ (* 0.5 sz) 0.002))))
            (doseq [x [(* 0.5 sx) (* -0.5 sx)]] (target! 0.09 (+ x (* (js/Math.sign x) 0.002)) (* 0.22 sy)))
            (joint! 0.075))
      ;; the head: round-edged, a target on each side, the neck a ball
      1 (do (add (RoundedBoxGeometry. (* 0.95 sx) sy (* 0.95 sz) 4 0.07) yellow 0 0 0)
            (doseq [x [(* 0.475 sx) (* -0.475 sx)]] (target! 0.11 (+ x (* (js/Math.sign x) 0.002)) 0.01))
            (joint! 0.05))
      ;; a limb: a capsule with a ball at its joint, and a hand or a foot
      (let [w (min sx sz)
            r (* 0.45 w)]
        (add (THREE/CapsuleGeometry. r (max 0.0 (- sy (* 2 r))) 6 12) yellow 0 0 0)
        (joint! (* 0.55 w))
        (case (long link)
          (3 5) (add (THREE/SphereGeometry. (* 0.6 sx) 14 10) hand 0 (- jy) 0)
          (7 9) (add (RoundedBoxGeometry. (* 1.15 sx) 0.07 (* 1.9 sz) 2 0.02) black
                     0 (+ (- jy) 0.035) (* -0.45 sz))
          nil)))
    g))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 50 1 0.1 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        ;; Held where it starts, in the air, until `start`; `reset` puts
        ;; it back there.
        paused? (atom true)
        tick-fps! (fps/meter! container)
        ;; `:pose` rather than `:state`, so that destructuring it cannot
        ;; shadow the atom holding it. It did, and `swap!` was called on
        ;; a map for as long as the shove button existed.
        world (atom {:scene nil :meshes {}})]
    (set! (.-background scene) (THREE/Color. 0x0b0d15))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.add scene (THREE/AmbientLight. 0xffffff 0.6))
    (let [sun (THREE/DirectionalLight. 0xfff1d8 1.2)]
      (.set (.-position sun) 4.0 8.0 5.0)
      (.add scene sun))
    (.set (.-position camera) 3.4 1.9 3.8)
    (let [orbit (OrbitControls. camera (.-domElement renderer))
          ground (THREE/Mesh. (THREE/BoxGeometry. 40 1 40)
                              (THREE/MeshStandardMaterial.
                               #js {:color 0x333a4d :roughness 0.95}))
          prop (THREE/Mesh. (THREE/BoxGeometry. 1.6 0.5 1.6)
                            (THREE/MeshStandardMaterial.
                             #js {:color 0x555f7a :roughness 0.9}))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0.0 0.7 0.0)
      (.set (.-position ground) 0.0 -0.5 0.0)
      (.add scene ground)
      (.add scene prop)
      (letfn [(obstacles []
                (cond-> [floor]
                  (.-obstacle controls) (conj (slab (.-tilt controls)))))
              (place! [^js m body]
                (let [[x y z] (:pos body)
                      [qx qy qz qw] (:rot body)]
                  (.set (.-position m) x y z)
                  (.set (.-quaternion m) qx qy qz qw)))
              (sync-meshes! []
                (let [{:keys [scene meshes]} @world]
                  (doseq [[mi {:keys [model pose]}] (map-indexed vector (:models scene))
                          {:keys [link body]} (ab/collision-bodies model (:q pose) (:base pose))]
                    (when-let [^js m (get meshes [:link mi link])]
                      (place! m body)))
                  (doseq [[i body] (map-indexed vector (:bodies scene))]
                    (when-let [^js m (get meshes [:rigid i])]
                      (place! m body)))))
              (sync-prop! []
                (let [ob (slab (.-tilt controls))
                      [x y z] (:pos ob)
                      [qx qy qz qw] (:rot ob)]
                  (set! (.-visible prop) (boolean (.-obstacle controls)))
                  (.set (.-position prop) x y z)
                  (.set (.-quaternion prop) qx qy qz qw)))
              (brick-mesh [body]
                (let [[sx sy sz] (:size body)]
                  (THREE/Mesh. (THREE/BoxGeometry. sx sy sz) (:brick @materials))))
              (rebuild! []
                (doseq [[_ ^js m] (:meshes @world)]
                  (.remove scene m)
                  ;; the materials are shared and stay
                  (.traverse m (fn [^js o] (when-let [geo (.-geometry o)] (.dispose geo)))))
                (let [model (figure (boolean (.-limits controls)))
                      spin (double (.-spin controls))
                      models (vec (for [[i [x h z yaw]] (map-indexed vector (drop-points (long (.-dummies controls))
                                                                                         (double (.-height controls))))]
                                    {:model model
                                     :pose (start-state model h (if (zero? i) spin (* spin (- 1.0 (* 0.15 (mod i 4)))))
                                                        x z yaw)}))
                      bodies (into (obstacles) (wall (.-brickRows controls)))
                      sc (pw/world bodies models (opts))
                      meshes (into (into {}
                                         (for [[mi {:keys [pose]}] (map-indexed vector models)
                                               {:keys [link body]}
                                               (ab/collision-bodies model (:q pose) (:base pose))]
                                           [[:link mi link] (dummy-part link body)]))
                                   (for [[i body] (map-indexed vector bodies)
                                         :when (not (rigid/static? body))]
                                     [[:rigid i] (brick-mesh body)]))]
                  (doseq [[_ ^js m] meshes] (.add scene m))
                  (reset! world {:scene sc :meshes meshes})
                  (reset! paused? true)
                  (sync-prop!)
                  (sync-meshes!)))
              (shove! []
                ;; A push on the torso, from wherever the camera is, so
                ;; it always shoves the figure away from the viewer.
                (let [dir (let [p (.-position camera)]
                            [(- (.-x p)) 0.15 (- (.-z p))])
                      len (js/Math.hypot (nth dir 0) (nth dir 1) (nth dir 2))
                      unit (mapv #(/ (double %) len) dir)]
                  (reset! paused? false)
                  ;; Every dummy, each woken first: a figure that has come
                  ;; to rest is asleep, and would ignore the push.
                  (swap! world update :scene
                         (fn [sc]
                           (reduce (fn [sc mi]
                                     (let [{:keys [model pose]} (nth (:models sc) mi)
                                           ;; At the chest rather than at the shoulder
                                           ;; joint, so it topples the figure instead
                                           ;; of spinning it about its own neck.
                                           chest (ab/frame-of model (:q pose) (:base pose) 0)
                                           where (v/add (:pos chest) (q/rotate (:rot chest) [0.0 0.25 0.0]))]
                                       (update-in (pw/wake sc [:link mi -1]) [:models mi :pose]
                                                  #(ab/apply-impulse model % 0 where unit
                                                                     (.-shoveForce controls)))))
                                   sc (range (count (:models sc))))))))
              (opts []
                {:gravity [0.0 -9.81 0.0]
                 :friction (.-friction controls)
                 :restitution (.-restitution controls)
                 :solver (keyword (.-solver controls))
                 :iterations (long (.-iterations controls))
                 :self-collide? (boolean (.-selfCollide controls))})
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)
                        n (max 1 (long (.-substeps controls)))
                        h (/ 1.0 (* 60.0 n))
                        o (opts)]
                    (when-not @paused?
                      (dotimes [_ n]
                        (swap! world update :scene #(pw/step (merge % o) h))))
                    (sync-meshes!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (set! (.-start controls) #(reset! paused? false))
        (set! (.-reset controls) rebuild!)
        (set! (.-shove controls) shove!)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "height" 0.4 4.0 0.1) (.onFinishChange rebuild!))
            (-> (.add controls "spin" -6.0 6.0 0.5) (.onFinishChange rebuild!))
            (.add controls "start")
            (.add controls "reset")
            (.add controls "shove")
            (-> (.add controls "shoveForce" 20 400 10) (.name "shove force"))
            (-> (.add controls "limits") (.name "joint limits") (.onChange rebuild!))
            (-> (.add controls "selfCollide") (.name "self collision"))
            ;; The slab is part of the simulation now rather than a
            ;; prop, so moving it means building the scene again.
            (-> (.add controls "obstacle") (.onChange rebuild!))
            (-> (.add controls "tilt" 0 40 1) (.onFinishChange rebuild!))
            (-> (.add controls "brickRows" 0 8 1) (.name "brick rows")
                (.onFinishChange rebuild!))
            (-> (.add controls "dummies" 1 40 1)
                (.onFinishChange (fn [_]
                                   ;; framed for the crowd, or for one as it was,
                                   ;; only when the count changes: a reset keeps
                                   ;; wherever the camera has been turned to
                                   (if (> (long (.-dummies controls)) 1)
                                     (do (.set (.-position camera) 5.5 4.2 7.0) (.set (.-target orbit) 0.0 1.2 0.0))
                                     (do (.set (.-position camera) 3.4 1.9 3.8) (.set (.-target orbit) 0.0 0.7 0.0)))
                                   (rebuild!))))
            (-> (.add controls "friction" 0.0 1.5 0.05))
            (-> (.add controls "restitution" 0.0 0.6 0.05))
            (-> (.add controls "solver" #js ["sequential-impulse" "convex"]))
            ;; Sweeps of the sequential impulse; the convex solve runs
            ;; Newton to convergence and ignores it.
            (-> (.add controls "iterations" 2 20 1))
            (-> (.add controls "substeps" 1 4 1))))
        {:start (fn [] (when-not @running?
                         (reset! running? true)
                         (on-resize)
                         (animate)))
         :stop (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "ragdoll")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
