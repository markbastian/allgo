(ns allgo.demo.ragdoll
  "A figure with joints, dropped on the floor.

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
  an upper arm is. Two metre-long links, the lower a thousand times the
  heavier, come nearly three metres apart at five substeps.

  Here there is no such number to report. A joint angle cannot be
  violated because there is nowhere for the violation to live. Watch the
  shoulders under `drop` from a height: the arms swing, and the point
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
  centimetres wider than the ribcage for the same reason: a forearm
  resting five millimetres inside the chest at the start is a contact
  that never goes away.

  ## What it costs

  Eleven bodies, twenty-eight degrees of freedom. On the JVM a step is
  0.17ms in flight and 1.3ms once it has landed on eighteen contacts;
  self collision adds about half of the first and nothing at all to the
  second, because a figure whose limbs are not inside each other is a
  figure with fewer contacts to solve.

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
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private identity-q [0.0 0.0 0.0 1.0])

(defn- deg [d] (* (double d) (/ js/Math.PI 180.0)))

(defn- box-inertia
  "A solid box's inertia about its own centre."
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
  centre of mass -- and with it the shape -- is half a length away."
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

(defn- figure
  "Eleven boxes and ten joints, roughly a person.

  The pelvis is the root and is free to move, which is what makes this a
  ragdoll rather than a mobile bolted to the ceiling. Masses are near
  enough a seventy kilogram adult; the exact numbers matter less than
  their ratios, and the ratios are what a constraint solver would
  struggle with."
  [limits?]
  (let [b (fn [l] (if limits? l (dissoc l :cone :twist :limit)))]
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

(def ^:private ^js controls
  #js {:height 1.8
       :tilt 18
       :obstacle true
       :limits true
       :selfCollide true
       :friction 0.7
       :restitution 0.0
       :iterations 8
       :substeps 1
       :spin 2.0
       :drop (fn [])
       :shove (fn [])})

(defn- start-state
  "The figure standing upright, `height` above the floor, given a turn."
  [model ^double height ^double spin]
  {:q (mapv (fn [l] (if (= :spherical (:joint l)) identity-q 0.0)) (ab/chain model))
   :qd (vec (repeat (ab/dof model) 0.0))
   :base {:rot (q/from-axis-angle [0.0 0.0 1.0] 0.15)
          :pos [0.0 height 0.0]
          :vel [0.0 spin 0.0 0.0 0.0 0.0]}})

(defn- material [kind]
  (case kind
    :head (THREE/MeshStandardMaterial. #js {:color 0xd8b08a :roughness 0.7})
    :torso (THREE/MeshStandardMaterial. #js {:color 0x4f6b9c :roughness 0.75})
    :limb (THREE/MeshStandardMaterial. #js {:color 0x8a9bbd :roughness 0.75})
    (THREE/MeshStandardMaterial. #js {:color 0x2b3040 :roughness 0.9})))

(defn- part-kind [link]
  (case (long link)
    -1 :torso
    0 :torso
    1 :head
    :limb))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 50 1 0.1 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom {:model nil :state nil :meshes {} :obstacles []})]
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
              (sync-meshes! []
                (let [{:keys [model state meshes]} @state]
                  (doseq [{:keys [link body]} (ab/collision-bodies
                                               model (:q state) (:base state))]
                    (when-let [^js m (get meshes link)]
                      (let [[x y z] (:pos body)
                            [qx qy qz qw] (:rot body)]
                        (.set (.-position m) x y z)
                        (.set (.-quaternion m) qx qy qz qw))))))
              (sync-prop! []
                (let [ob (slab (.-tilt controls))
                      [x y z] (:pos ob)
                      [qx qy qz qw] (:rot ob)]
                  (set! (.-visible prop) (boolean (.-obstacle controls)))
                  (.set (.-position prop) x y z)
                  (.set (.-quaternion prop) qx qy qz qw)))
              (rebuild! []
                (doseq [[_ ^js m] (:meshes @state)]
                  (.remove scene m)
                  (.dispose (.-geometry m)))
                (let [model (figure (boolean (.-limits controls)))
                      st (start-state model (.-height controls) (.-spin controls))
                      meshes (into {}
                                   (for [{:keys [link body]}
                                         (ab/collision-bodies model (:q st) (:base st))
                                         :let [[sx sy sz] (:size body)]]
                                     [link (doto (THREE/Mesh.
                                                  (THREE/BoxGeometry. sx sy sz)
                                                  (material (part-kind link)))
                                             (-> .-castShadow (set! true)))]))]
                  (doseq [[_ ^js m] meshes] (.add scene m))
                  (reset! state {:model model :state st :meshes meshes})
                  (sync-prop!)
                  (sync-meshes!)))
              (shove! []
                ;; A push on the torso, from wherever the camera is, so
                ;; it always shoves the figure away from the viewer.
                (let [{:keys [model state]} @state
                      dir (let [p (.-position camera)]
                            [(- (.-x p)) 0.15 (- (.-z p))])
                      len (js/Math.hypot (nth dir 0) (nth dir 1) (nth dir 2))
                      unit (mapv #(/ (double %) len) dir)
                      where (:pos (ab/frame-of model (:q state) (:base state) 0))]
                  (swap! state update :state
                         #(ab/apply-impulse model % 0 where unit 30.0))))
              (opts []
                {:gravity [0.0 -9.81 0.0]
                 :obstacles (obstacles)
                 :friction (.-friction controls)
                 :restitution (.-restitution controls)
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
                        {:keys [model]} @state
                        n (max 1 (long (.-substeps controls)))
                        h (/ 1.0 (* 60.0 n))
                        o (opts)]
                    (dotimes [_ n]
                      (swap! state update :state #(ab/step % model h o)))
                    (sync-meshes!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (set! (.-drop controls) rebuild!)
        (set! (.-shove controls) shove!)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "height" 0.4 4.0 0.1) (.onFinishChange rebuild!))
            (-> (.add controls "spin" -6.0 6.0 0.5) (.onFinishChange rebuild!))
            (.add controls "drop")
            (.add controls "shove")
            (-> (.add controls "limits") (.name "joint limits") (.onChange rebuild!))
            (-> (.add controls "selfCollide") (.name "self collision"))
            (-> (.add controls "obstacle") (.onChange sync-prop!))
            (-> (.add controls "tilt" 0 40 1) (.onChange sync-prop!))
            (-> (.add controls "friction" 0.0 1.5 0.05))
            (-> (.add controls "restitution" 0.0 0.6 0.05))
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
