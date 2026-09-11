(ns allgo.demo.arm
  "Inverse kinematics for a six-axis arm: click a target and watch it
  reach, after John J. Craig's Introduction to Robotics.

  The number in the corner is the point of this demo. A six-axis arm whose
  last three axes meet at a point -- which is how almost every industrial
  robot is built, because it is the condition that makes the algebra
  solvable -- can reach any pose it can reach in exactly **eight** ways.
  The shoulder can pass either side of the target, the elbow can be up or
  down, and the wrist can arrive flipped. Click a ball and all eight are
  found at once, exactly, with no guessing and no iteration.

  `showAll` draws them together as ghosts. They are not variations on a
  theme: two of them can have the arm folded completely differently while
  the gripper sits in precisely the same place, which is why a robot
  program has to name a configuration and not just a position.

  `useLimits` is where it stops being a curiosity. Real joints do not turn
  forever, and of the eight solutions typically only two or four are
  actually available -- that is the number an engineer cares about, and
  the overlay shows both.

  Targets outside the workspace report no solutions at all. That is not a
  failure to converge; the closed form knows there is nothing to find, and
  says so.

  `solver` switches to damped least squares, which is what you must use
  when a chain has no closed form. Watch the count collapse to one: a
  descent finds whichever configuration it happens to walk to, and can
  never tell you how many others there were.

  `drive` decides how the arm gets there. `trajectory` interpolates the
  joint angles, which is what a robot controller actually commands.
  `servo physics` hands the angles to the XPBD servo joints of
  `allgo.physics.joint` and lets the arm be dragged there by its motors --
  the same joints as the ragdoll demo, doing the job they are named for."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as quat]
            [allgo.geometry.vec3 :as v]
            [allgo.kinematics.analytic :as ik]
            [allgo.kinematics.chain :as k]
            [allgo.kinematics.numeric :as nik]
            [allgo.physics.joint :as joint]
            [allgo.physics.rigid :as rigid]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:solver     "closed form"
       :drive      "trajectory"
       :showAll    true
       :useLimits  true
       :showAxes   false
       :speed      1.6
       :targets    7})

(def ^:private half-pi (/ js/Math.PI 2))

(def ^:private puma
  "The PUMA 560, in Craig's modified DH parameters and metres. The arm his
  chapter 4 solves, and the one every later text compares against."
  (k/chain [{:alpha 0.0         :a 0.0    :d 0.0    :limits [-2.79 2.79]}
            {:alpha (- half-pi) :a 0.0    :d 0.0    :limits [-3.92 0.78]}
            {:alpha 0.0         :a 0.4318 :d 0.1245 :limits [-0.78 3.92]}
            {:alpha (- half-pi) :a 0.0203 :d 0.4318 :limits [-1.92 2.96]}
            {:alpha half-pi     :a 0.0    :d 0.0    :limits [-1.75 1.75]}
            {:alpha (- half-pi) :a 0.0    :d 0.0    :limits [-4.64 4.64]}]
           ;; A gripper sticking out of the wrist. The arm's own last three
           ;; links have no length -- their axes all meet at a point, which
           ;; is what makes the closed form possible -- so without a tool
           ;; there is nothing between the elbow and the end of the world
           ;; to draw, and nothing for touching a ball to mean.
           {:rot quat/identity-q :pos [0.0 0.0 0.11]}))

;; A hinge in `allgo.physics.joint` turns about its frame's x axis where a
;; Denavit-Hartenberg joint turns about z. This is the rotation between the
;; two conventions, and the only thing needed to hand one to the other.
(def ^:private x->z (quat/from-axis-angle [0.0 1.0 0.0] (- half-pi)))

(defn- reach
  "How far the arm can get, for scattering targets it can mostly reach."
  []
  (+ 0.4318 0.4318 0.0203))

;; ---------------------------------------------------------------------------
;; Solving

(defn- solve-all
  "Every configuration that reaches the target, by the chosen method."
  [target]
  (if (= "closed form" (.-solver controls))
    (ik/solutions puma target {:limits? (.-useLimits controls)})
    ;; Descent returns whatever it walked to. Several seeds find several,
    ;; and never say whether that is all of them.
    (nik/solve-from-many puma target 8 {:limits? (.-useLimits controls)
                                        :iterations 200})))

(defn- choose
  "The configuration that moves the joints least, which is what a
  controller picks and for the reason it picks it."
  [solutions from]
  (when (seq solutions)
    (apply min-key
           (fn [s] (reduce + (map (fn [a b] (js/Math.abs (k/wrap-angle (- a b)))) s from)))
           solutions)))

;; ---------------------------------------------------------------------------
;; Rendering the arm

(defn- segment!
  "Points a cylinder from `a` to `b`. Cylinders are built along y, so this
  is the rotation taking y onto the link."
  [^js mesh a b]
  (let [d (v/sub b a)
        len (v/length d)]
    (if (< len 1e-6)
      (set! (.-visible mesh) false)
      (let [mid (v/scale (v/add a b) 0.5)
            [qx qy qz qw] (quat/from-vectors [0.0 1.0 0.0] d)]
        (set! (.-visible mesh) true)
        (.set (.-position mesh) (nth mid 0) (nth mid 1) (nth mid 2))
        (.set (.-quaternion mesh) qx qy qz qw)
        (.set (.-scale mesh) 1.0 len 1.0)))))

(defn- make-arm
  "One cylinder per link and a disc per joint, reused every frame.

  Ghost arms are drawn thinner as well as fainter: eight of them over the
  same point otherwise read as one solid mass rather than as eight ways of
  folding the same linkage."
  [^js scene colour opacity]
  (let [mat (THREE/MeshPhongMaterial.
             #js {:color colour :transparent (< opacity 1.0) :opacity opacity})
        joint-mat (THREE/MeshPhongMaterial.
                   #js {:color 0x2f3646 :transparent (< opacity 1.0) :opacity opacity})
        links (vec (for [_ (range 7)]
                     (let [m (THREE/Mesh. (THREE/CylinderGeometry. (if (< opacity 1.0) 0.016 0.035)
                                                                   (if (< opacity 1.0) 0.016 0.035) 1.0 16) mat)]
                       (.add scene m) m)))
        discs (vec (for [_ (range 6)]
                     (let [m (THREE/Mesh. (THREE/CylinderGeometry. 0.055 0.055 0.07 20) joint-mat)]
                       (.add scene m) m)))]
    {:links links :discs discs :material mat :joint-material joint-mat}))

(defn- draw-tool!
  "The gripper, from the wrist out to the tip that does the touching."
  [^js mesh values]
  (let [fs (k/frames puma values)]
    (segment! mesh (:pos (peek fs)) (:pos (k/pose puma values)))))

(defn- draw-arm!
  "Cylinders along the links, discs across the joint axes.

  The discs are not decoration. Several of a PUMA's frames sit on top of
  one another -- the first three are all at the shoulder, and the last
  three all at the wrist, which is the condition that makes the closed
  form exist -- so the links between them have no length and nothing to
  draw. Without the discs the arm reads as a bare stick with a shoulder
  and a wrist that are not there."
  [{:keys [links discs]} values show-axes?]
  (let [fs (k/frames puma values)]
    (draw-tool! (nth links 6) values)
    (dotimes [i 6]
      (segment! (nth links i) (:pos (nth fs i)) (:pos (nth fs (inc i))))
      (let [^js disc (nth discs i)
            f (nth fs (inc i))
            ;; The joint's own axis, which a disc lying across it shows.
            [qx qy qz qw] (quat/mul (:rot f) (quat/from-axis-angle [1.0 0.0 0.0] half-pi))
            [px py pz] (:pos f)]
        (set! (.-visible disc) true)
        (.set (.-position disc) px py pz)
        (.set (.-quaternion disc) qx qy qz qw)
        ;; Stretched into a spindle, a disc becomes a visible joint axis.
        (.set (.-scale disc) 1.0 (if show-axes? 4.5 1.0) 1.0)))))

(defn- set-arm-visible! [{:keys [links discs]} on?]
  (doseq [^js m links] (set! (.-visible m) on?))
  ;; Ghosts show their links only: eight sets of joint discs over one point
  ;; bury the arm that is actually moving.
  (doseq [^js m discs] (set! (.-visible m) false)))

;; ---------------------------------------------------------------------------
;; The two ways of getting there

(defn- build-physics
  "A rigid link per joint, hinged and driven by servos.

  Every link is given the same mass whatever its length. The wrist links
  are geometrically zero -- their joints all meet at a point, which is the
  whole reason a closed form exists -- so sizing mass by length hands the
  solver a ratio of hundreds to one across a joint, and it comes apart."
  []
  (let [fs (k/frames puma (k/home puma))
        centres (vec (for [i (range 6)]
                       (v/scale (v/add (:pos (nth fs i)) (:pos (nth fs (inc i)))) 0.5)))
        base (rigid/box {:size [0.2 0.2 0.2] :pos (:pos (first fs))})
        links (vec (for [i (range 6)]
                     (rigid/box {:size [0.18 0.12 0.12] :density 700.0
                                 :pos (nth centres i)})))
        bodies (into [base] links)
        joints (vec (for [i (range 6)]
                      (joint/servo bodies {:a i :b (inc i)
                                           :at (:pos (nth fs (inc i)))
                                           :rot (quat/mul (:rot (nth fs (inc i))) x->z)
                                           :target-angle 0.0})))]
    {:bodies bodies :constraints joints :gravity [0.0 0.0 0.0]
     :dt (/ 1.0 60.0) :substeps 20}))

(defn- physics-angles
  "What the servos have actually achieved.

  The arm is drawn from these through forward kinematics rather than from
  where the bodies ended up, so the links stay exactly the length they
  are. The solver's positional drift is a property of the solver, not of
  a robot, and the angles are what it tracks well -- to a thousandth of a
  radian, where the bodies themselves wander centimetres."
  [world]
  (mapv (fn [jt]
          (let [[_ r0 _ r1] (joint/frames (:bodies world) jt)]
            (joint/signed-angle (quat/rotate r0 [1.0 0.0 0.0])
                                (quat/rotate r0 [0.0 1.0 0.0])
                                (quat/rotate r1 [0.0 1.0 0.0]))))
        (:constraints world)))

(defn- step-trajectory
  "A joint-space move: every joint starts and stops together, which is
  what a controller commands when it is told to go to a configuration
  rather than to follow a path."
  [current target dt speed]
  (let [deltas (map (fn [a b] (k/wrap-angle (- b a))) current target)
        biggest (reduce max 0.0 (map js/Math.abs deltas))]
    (if (< biggest 1e-5)
      (vec target)
      (let [step (* speed dt)
            t (min 1.0 (/ step biggest))]
        (mapv (fn [a d] (+ a (* t d))) current deltas)))))

;; ---------------------------------------------------------------------------

(defn- scatter-targets
  "Balls to reach for: most inside the workspace, one or two beyond it so
  that no solution is a visible answer too."
  [n]
  (let [r (reach)]
    (vec (for [i (range n)]
           (let [out? (and (pos? i) (zero? (mod i 5)))
                 radius (if out? (* r (+ 1.25 (* 0.3 (js/Math.random))))
                            (* r (+ 0.35 (* 0.5 (js/Math.random)))))
                 phi (* 2 js/Math.PI (/ i n))
                 z (* radius (- (* 1.1 (js/Math.random)) 0.35))]
             {:pos [(* radius (js/Math.cos phi))
                    (* radius (js/Math.sin phi) 0.55)
                    z]
              :beyond? out?})))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        raycaster (THREE/Raycaster.)
        running? (atom false)
        tick-fps! (fps/meter! container)
        overlay  (js/document.createElement "div")
        state    (atom {})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.8))
    (let [sun (THREE/DirectionalLight. 0xffffff 0.9)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 12 12)
                             (THREE/MeshPhongMaterial. #js {:color 0x121828}))]
      (set! (.-x (.-rotation floor)) (- half-pi))
      (set! (.-y (.-position floor)) -0.9)
      (.add scene floor))
    (.add scene (doto (THREE/GridHelper. 3 12 0x2b3450 0x1b2236)
                  (-> .-position (.setY -0.899))))
    ;; A column from the floor to the shoulder. The table puts the shoulder
    ;; at the origin, so without this the arm floats.
    (let [column (THREE/Mesh. (THREE/CylinderGeometry. 0.055 0.085 0.9 20)
                              (THREE/MeshPhongMaterial. #js {:color 0x3a4356}))]
      (.set (.-position column) 0 -0.45 0)
      (.add scene column))
    (set! (.-className overlay) "demo-readout")
    (set! (.-cssText (.-style overlay))
          (str "position:absolute;left:14px;top:14px;z-index:5;pointer-events:none;"
               "font:12px ui-monospace,monospace;color:#dfe7f5;background:rgba(0,0,0,0.55);"
               "padding:8px 10px;border-radius:8px;white-space:pre;line-height:1.5"))
    (.appendChild container overlay)
    (set! (.-position (.-style container)) "relative")
    (.set (.-position camera) 1.35 0.95 1.35)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 -0.05 0)
      (letfn [(rebuild! []
                (doseq [^js m (:target-meshes @state)] (.remove scene m))
                (let [targets (scatter-targets (.-targets controls))
                      meshes (mapv (fn [t]
                                     (let [m (THREE/Mesh.
                                              (THREE/SphereGeometry. 0.045 18 14)
                                              (THREE/MeshPhongMaterial.
                                               #js {:color (if (:beyond? t) 0x8b4a52 0x4f7ac0)}))
                                           [x y z] (:pos t)]
                                       (.set (.-position m) x y z)
                                       (.add scene m)
                                       m))
                                   targets)]
                  (swap! state assoc :targets targets :target-meshes meshes
                         :chosen nil :solutions [] :all-count 0 :picked nil
                         :goal (k/home puma))))
              (resize! []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (.setSize renderer w h)
                    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera))))
              (report! []
                (let [{:keys [picked all-count solutions]} @state]
                  (set! (.-textContent overlay)
                        (cond
                          (nil? picked) "click a ball for the arm to touch"
                          (zero? all-count)
                          (str "target " picked "\nout of reach\n0 solutions")
                          :else
                          (str "target " picked "\n"
                               all-count " solution" (if (= 1 all-count) "" "s")
                               (if (= "closed form" (.-solver controls))
                                 " (exact, all of them)"
                                 " (descent found these; there may be more)")
                               "\n"
                               (count solutions)
                               (if (.-useLimits controls)
                                 " within joint limits" " shown"))))))
              (pick! [^js e]
                (let [rect (.getBoundingClientRect (.-domElement renderer))
                      nx (- (* 2 (/ (- (.-clientX e) (.-left rect)) (.-width rect))) 1)
                      ny (- 1 (* 2 (/ (- (.-clientY e) (.-top rect)) (.-height rect))))]
                  (.setFromCamera raycaster (THREE/Vector2. nx ny) camera)
                  (let [hits (.intersectObjects raycaster (clj->js (:target-meshes @state)) false)]
                    (when (pos? (.-length hits))
                      (let [^js hit (aget hits 0)
                            idx (.indexOf (to-array (:target-meshes @state)) (.-object hit))
                            t (nth (:targets @state) idx)
                            ;; Approach the ball from wherever the arm can:
                            ;; a position with a free orientation would be a
                            ;; different problem, so aim the tool along the
                            ;; line from the base outward.
                            ;; Touching a ball fixes where the tool goes
                            ;; and not which way it faces, so the
                            ;; orientation is ours to choose. Pointing the
                            ;; tool outward along the radius is the natural
                            ;; choice and the one the joint limits like:
                            ;; it leaves a configuration reachable for
                            ;; about ninety-six targets in a hundred,
                            ;; where an arbitrary fixed orientation leaves
                            ;; most of them with none.
                            target {:pos (:pos t)
                                    :rot (quat/from-vectors [0.0 0.0 1.0] (:pos t))}
                            from (:goal @state)
                            all (if (= "closed form" (.-solver controls))
                                  (ik/solutions puma target)
                                  (nik/solve-from-many puma target 8 {:limits? false
                                                                      :iterations 200}))
                            shown (solve-all target)
                            best (choose (if (seq shown) shown all) from)]
                        (swap! state assoc :picked idx :all-count (count all)
                               :solutions shown :goal (or best (:goal @state)))
                        (report!))))))
              (tick []
                (when @running?
                  (js/requestAnimationFrame tick)
                  (let [t0 (js/performance.now)
                        {:keys [current goal ghosts arm physics]} @state
                        dt (/ 1.0 60.0)
                        servo? (= "servo physics" (.-drive controls))
                        [current physics]
                        (if servo?
                          (let [w (joint/step
                                   (update physics :constraints
                                           (fn [cs] (mapv #(assoc %1 :target-angle %2) cs goal))))]
                            [(physics-angles w) w])
                          [(step-trajectory current goal dt (.-speed controls)) physics])]
                    (swap! state assoc :current current :physics physics)
                    (draw-arm! arm current (.-showAxes controls))
                    ;; Ghosts of every other configuration that reaches the
                    ;; same point.
                    (let [shown (if (.-showAll controls) (:solutions @state) [])]
                      (doseq [[i g] (map-indexed vector ghosts)]
                        (if (< i (count shown))
                          (do (set-arm-visible! g true)
                              (draw-arm! g (nth shown i) false))
                          (set-arm-visible! g false))))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.addEventListener (.-domElement renderer) "pointerdown" pick!)
        (.observe (js/ResizeObserver. resize!) container)
        (resize!)
        (swap! state assoc
               :arm (make-arm scene 0xd9a15c 1.0)
               :ghosts (vec (repeatedly 8 #(make-arm scene 0x7fc4ff 0.32)))
               :current (k/home puma)
               :goal (k/home puma)
               :physics (build-physics))
        (rebuild!)
        (report!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "solver" #js ["closed form" "damped least squares"])
                (.onChange report!))
            (.add controls "drive" #js ["trajectory" "servo physics"])
            (-> (.add controls "useLimits") (.onChange report!))
            (.add controls "showAll")
            (.add controls "showAxes")
            (.add controls "speed" 0.2 6 0.1)
            (-> (.add controls "targets" 3 12 1) (.onFinishChange rebuild!))
            (.add #js {:reset rebuild!} "reset")))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "arm")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
