(ns allgo.demo.motorcycle
  "A motorcycle on a plane, ridden with the keyboard.

  Everything that moves is `allgo.simulation.motorcycle`: one articulated
  model -- frame, swingarm, two wheels, steering head and fork, joined by
  exact joints -- on a floor, with the tyres' grip solved by
  `allgo.physics.world`'s sequential impulse sweep at four substeps a
  frame. The docstring there says why that and not the other two
  solvers; this namespace only draws it and reads the keys.

  ## What to try

  **Lean it.** A and D ask the rider for a lean, and the rider gets there
  the way a real one does: the bars go the *wrong* way first, the wheels
  run out from under the bike, it tips, and then the bars follow it into
  the turn. Watch the handlebar at the moment you press a key.

  **Take your hands off.** `hands on bars` off and nobody is steering.
  This bike's geometry does not hold it up on its own at these speeds, so
  over it goes -- which is the evidence that the rider is doing the work
  when it stays up.

  **Brake hard.** The fork dives against its springs and the rear goes
  light and locks: a rear wheel with most of the weight gone from it
  skids, and you can see it stop turning while the bike is still moving.

  **Ride over the kerbs.** The front drops the fork a few centimetres,
  the rear shock -- the rod between the frame and the swingarm -- takes
  it a moment later, and both settle in a couple of bounces.

  **Watch the rider.** They are a jointed figure of their own, pinned to
  the bike at the seat, the grips and the pegs, holding the riding pose
  with a spring and damper at every joint -- so they sway as the bike
  pitches, their arms follow the bars, and in a turn they tip their upper
  body in toward the inside.

  **Bail out.** B, or the button, and the rider jumps off the side. The
  pins let go, the muscles switch off, and what leaves the saddle is a
  plain ragdoll held together by its joint limits, colliding with the
  bike on its way past it. Nobody is riding any more, so the bike -- kicked
  by the jump and with nobody on the bars -- goes over too. The same
  happens without the button whenever the bike goes down with the rider
  on it. `ragdoll rider` off puts the
  old lumped rider back, whose weight is simply part of the frame's.

  **Open it up.** The rear rises under power rather than squatting: the
  motor's torque reacts on the swingarm, and at this geometry that lifts
  the frame. The front goes light as weight comes off it.

  Until a key is pressed the rider laps a circle on their own, over a
  ring of kerbs, holding about 12 m/s -- so there is something to watch,
  and it is the same rider doing it: the autopilot only chooses the lean
  to ask for. Take over and the ramp is off to the right of the start."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            [allgo.simulation.motorcycle :as m]
            ["lil-gui" :default GUI]
            ["three" :as THREE]))

(def ^:private ^js controls
  #js {:riderBody true
       :handsOn true
       :autopilot true
       :maxLean 40
       :kerbs true
       :ramp true
       :startSpeed 10
       :substeps 4
       :rearRate 25
       :frontRate 20
       :damping 1.0
       :reset (fn [])
       :bail (fn [])})

;; ---------------------------------------------------------------------------
;; The course

(def ^:private loop-radius
  "The circle the autopilot laps, which starts where the bike does and
  turns left."
  35.0)

(def ^:private loop-centre [0.0 0.0 (- loop-radius)])

(defn- on-loop
  "The point `theta` radians round the loop from the start, and the angle
  of its tangent about the vertical."
  [theta]
  (let [theta (double theta)]
    [[(* loop-radius (Math/sin theta)) 0.0 (- (* loop-radius (Math/cos theta)) loop-radius)]
     theta]))

(defn- kerbs
  "Low kerbs laid across the loop, so the autopilot rides over them every
  lap and the suspension has something to do."
  []
  (vec (for [theta [0.55 1.5 2.5 3.4 4.5 5.5]
             :let [[[x _ z] across] (on-loop theta)]]
         (rigid/box {:pos [x 0.04 z] :size [0.45 0.08 4.0]
                     :rot (q/from-axis-angle [0.0 1.0 0.0] across)}))))

(defn- ramp
  "A kicker off to the side of the loop, for when you take over."
  []
  (rigid/box {:pos [0.0 0.35 22.0] :size [5.0 0.3 4.0]
              :rot (q/from-axis-angle [0.0 0.0 1.0] (/ (* 9.0 Math/PI) 180.0))}))

(defn- obstacles []
  (cond-> [(m/ground 4000.0)]
    (.-kerbs controls) (into (kerbs))
    (.-ramp controls) (conj (ramp))))

(defn- config []
  (let [d (.-damping controls)]
    (assoc (cond-> m/defaults (.-riderBody controls) m/with-rider)
           :rear-rate (* 1000.0 (.-rearRate controls))
           :front-rate (* 1000.0 (.-frontRate controls))
           :rear-damping (* d (:rear-damping m/defaults))
           :front-damping (* d (:front-damping m/defaults)))))

;; ---------------------------------------------------------------------------
;; Meshes

(defn- mat [color & [rough metal]]
  (THREE/MeshStandardMaterial. #js {:color color :roughness (or rough 0.6)
                                    :metalness (or metal 0.1)}))

(defn- box-mesh
  "A box of `size` centred at `pos` in its parent, turned by `rot`."
  [size pos material & [rot]]
  (let [[sx sy sz] size
        [x y z] pos
        m (THREE/Mesh. (THREE/BoxGeometry. sx sy sz) material)]
    (.set (.-position m) x y z)
    (when rot (let [[qx qy qz qw] rot] (.set (.-quaternion m) qx qy qz qw)))
    m))

(defn- rod
  "A cylinder from `a` to `b` in its parent's frame."
  [a b radius material]
  (let [d (v/sub b a)
        len (v/length d)
        m (THREE/Mesh. (THREE/CylinderGeometry. radius radius 1 10) material)
        mid (v/add a (v/scale d 0.5))
        up (THREE/Vector3. 0 1 0)
        dir (.normalize (THREE/Vector3. (nth d 0) (nth d 1) (nth d 2)))]
    (.set (.-scale m) 1 len 1)
    (.set (.-position m) (nth mid 0) (nth mid 1) (nth mid 2))
    (.setFromUnitVectors (.-quaternion m) up dir)
    m))

(defn- set-rod!
  "Stretch an existing `rod` mesh to run from `a` to `b`, in world
  coordinates."
  [^js m a b]
  (let [d (v/sub b a)
        len (max 1e-6 (v/length d))
        mid (v/add a (v/scale d 0.5))]
    (.set (.-scale m) 1 len 1)
    (.set (.-position m) (nth mid 0) (nth mid 1) (nth mid 2))
    (.setFromUnitVectors (.-quaternion m) (THREE/Vector3. 0 1 0)
                         (.normalize (THREE/Vector3. (nth d 0) (nth d 1) (nth d 2))))))

(defn- wheel-group
  "A tyre as the torus the physics uses, a rim, and spokes -- the spokes
  are what show it turning."
  [cfg]
  (let [{:keys [wheel-radius tyre-radius]} cfg
        major (- wheel-radius tyre-radius)
        g (THREE/Group.)
        tyre (THREE/Mesh. (THREE/TorusGeometry. major tyre-radius 12 40) (mat 0x15161a 0.9))
        rim (THREE/Mesh. (THREE/TorusGeometry. (- major tyre-radius 0.005) 0.012 6 40)
                         (mat 0xb8bcc6 0.3 0.8))
        hub (THREE/Mesh. (THREE/CylinderGeometry. 0.05 0.05 0.12 12) (mat 0x8d93a1 0.4 0.7))]
    (.set (.-rotation hub) (/ Math/PI 2) 0 0)
    (.add g tyre)
    (.add g rim)
    (.add g hub)
    (dotimes [i 6]
      (let [a (* i (/ Math/PI 3))
            spoke (rod [0.0 0.0 0.0]
                       [(* (- major tyre-radius) (Math/cos a)) (* (- major tyre-radius) (Math/sin a)) 0.0]
                       0.01 (mat 0xd0d4dc 0.3 0.8))]
        (.add g spoke)))
    g))

(defn- bike-meshes
  "One group per body, each built in that body's own frame so placing it
  is copying a pose. `lumped?` draws a rider sat on the frame, for when
  the rider's weight is part of the frame's rather than a body of their
  own."
  [cfg lumped?]
  (let [{:keys [pivot swingarm head fork-length]} cfg
        paint (mat 0xc2302b 0.35 0.3)
        dark (mat 0x23262e 0.7 0.2)
        metal (mat 0xb8bcc6 0.3 0.8)
        leathers (mat 0x2a2f3d 0.8)
        helmet (mat 0xe8e8ee 0.3 0.2)
        frame (THREE/Group.)
        arm (THREE/Group.)
        bars (THREE/Group.)
        slider (THREE/Group.)
        fork-dir (v/negate (m/steering-axis cfg))]
    ;; The frame: engine, tank, seat and tail, and the pegs.
    (doto frame
      (.add (box-mesh [0.46 0.30 0.28] [0.05 -0.08 0.0] dark))
      (.add (box-mesh [0.42 0.16 0.30] [0.18 0.16 0.0] paint))
      (.add (box-mesh [0.46 0.08 0.26] [-0.24 0.14 0.0] (mat 0x111217 0.9)))
      (.add (box-mesh [0.34 0.10 0.20] [-0.55 0.18 0.0] paint))
      (.add (rod [0.40 0.20 0.0] head 0.03 paint))
      (.add (rod [-0.05 -0.20 0.1] [-0.75 0.02 0.14] 0.035 metal))
      (.add (rod [-0.10 -0.14 -0.26] [-0.10 -0.14 0.26] 0.015 metal)))
    ;; A rider sat on the frame, when their mass is the frame's: only a
    ;; shape, moving with it.
    (when lumped?
      (doto frame
        (.add (box-mesh [0.24 0.46 0.34] [-0.18 0.44 0.0] leathers
                        (q/from-axis-angle [0.0 0.0 1.0] -0.5)))
        (.add (doto (THREE/Mesh. (THREE/SphereGeometry. 0.13 16 12) helmet)
                (-> .-position (.set 0.02 0.78 0.0))))
        (.add (box-mesh [0.42 0.13 0.13] [-0.05 0.20 0.17] leathers))
        (.add (box-mesh [0.42 0.13 0.13] [-0.05 0.20 -0.17] leathers))))
    ;; Swingarm: two arms from the pivot to the axle.
    (doseq [z [0.11 -0.11]]
      (.add arm (rod [0.0 0.0 z] [(nth swingarm 0) (nth swingarm 1) z] 0.025 metal)))
    ;; Steering head: the stanchions stay with the bars; the lower legs
    ;; ride on the slider, so the fork is seen to telescope.
    (doto bars
      (.add (rod [0.0 0.12 -0.36] [0.0 0.12 0.36] 0.016 metal))
      (.add (box-mesh [0.08 0.04 0.26] [0.0 0.0 0.0] dark)))
    (doseq [z [0.09 -0.09]]
      (.add bars (rod [0.0 0.0 z] (v/add (v/scale fork-dir 0.34) [0.0 0.0 z]) 0.022 metal))
      (.add slider (rod (v/add (v/scale fork-dir 0.30) [0.0 0.0 z])
                        (v/add (v/scale fork-dir fork-length) [0.0 0.0 z]) 0.03 dark)))
    {:frame frame :swingarm arm :bars bars :slider slider
     :rear (wheel-group cfg) :front (wheel-group cfg)
     :shock (rod [0.0 0.0 0.0] [0.0 1.0 0.0] 0.028 (mat 0xf0b43c 0.4 0.4))
     :pivot pivot}))

(defn- rider-meshes
  "The ragdoll rider, a mesh per part, keyed by link: the same shapes it
  collides with, so what is drawn is what hits the floor."
  [model pose]
  (let [leathers (mat 0x2a2f3d 0.8)
        helmet (mat 0xe8e8ee 0.3 0.2)]
    (into {}
          (for [{:keys [link body]} (ab/collision-bodies model (:q pose) (:base pose))]
            [link (if (= :ball (:shape body))
                    (THREE/Mesh. (THREE/SphereGeometry. (:radius body) 16 12) helmet)
                    (let [[sx sy sz] (:size body)]
                      (THREE/Mesh. (THREE/BoxGeometry. sx sy sz) leathers)))]))))

(defn- place! [^js o {:keys [pos rot]}]
  (let [[x y z] pos [qx qy qz qw] rot]
    (.set (.-position o) x y z)
    (.set (.-quaternion o) qx qy qz qw)))

;; ---------------------------------------------------------------------------
;; The rider's hands and feet

(def ^:private keymap
  {"ArrowUp" :throttle "KeyW" :throttle
   "ArrowDown" :brake "KeyS" :brake
   "ArrowLeft" :left "KeyA" :left
   "ArrowRight" :right "KeyD" :right})

(defn- approach
  "Move `x` toward `target` by at most `rate * dt`: keys are on or off,
  and a throttle or a lean that snapped would jolt the bike."
  [x target rate dt]
  (let [x (double x) step (* (double rate) (double dt))]
    (+ x (max (- step) (min step (- (double target) x))))))

(def ^:private cruise 12.0)

(defn- autopilot
  "Lap the loop: the lean a steady turn of its radius needs at this speed,
  plus whatever turns the bike toward where it should be heading -- along
  the circle, drawn back onto it if it has drifted off."
  [pose speed max-lean]
  (let [{:keys [pos rot]} (:base pose)
        fwd (let [f (q/rotate rot [1.0 0.0 0.0])] (v/normalize [(nth f 0) 0.0 (nth f 2)]))
        r (let [d (v/sub pos loop-centre)] [(nth d 0) 0.0 (nth d 2)])
        dist (max 1e-6 (v/length r))
        out (v/scale r (/ 1.0 dist))
        along [(nth out 2) 0.0 (- (double (nth out 0)))]
        want (v/normalize (v/sub along (v/scale out (* 0.08 (- dist loop-radius)))))
        err (Math/atan2 (nth (v/cross fwd want) 1) (v/dot fwd want))
        steady (Math/atan (/ (* speed speed) (* 9.81 loop-radius)))]
    (max (- max-lean) (min max-lean (+ steady (* 1.4 err))))))

(defn- rider-inputs
  "What the rider is asking for this frame, smoothed toward what the keys
  (or the autopilot) want."
  [inputs held pose speed dt]
  (let [max-lean (* (/ Math/PI 180.0) (.-maxLean controls))
        auto? (and (.-autopilot controls) (empty? held))
        want-throttle (cond (held :throttle) 1.0
                            auto? (max 0.0 (min 1.0 (+ 0.1 (* 0.3 (- cruise (double speed))))))
                            :else 0.0)
        want-brake (if (held :brake) 1.0 0.0)
        want-lean (cond
                    auto? (autopilot pose (double speed) max-lean)
                    (and (held :left) (not (held :right))) max-lean
                    (and (held :right) (not (held :left))) (- max-lean)
                    :else 0.0)]
    {:throttle (approach (:throttle inputs) want-throttle 2.5 dt)
     :brake (approach (:brake inputs) want-brake 4.0 dt)
     :lean (approach (:lean inputs) want-lean 1.2 dt)
     :hands? (boolean (.-handsOn controls))}))

;; ---------------------------------------------------------------------------
;; Readout

(defn- hud! [^js el cfg pose inputs rider]
  (let [{:keys [speed lean steer rear-travel front-travel upright?
                rear-wheel-speed]} (m/telemetry cfg pose)
        deg #(Math/round (* (/ 180.0 Math/PI) (double %)))
        bar (fn [x lo hi]
              (let [f (/ (- (double x) (double lo)) (- (double hi) (double lo)))]
                (str "<span class='moto-bar'><i style='width:"
                     (Math/round (* 100 (max 0.0 (min 1.0 f)))) "%'></i></span>")))
        [rlo rhi] (:rear-travel cfg)
        [flo fhi] (:front-travel cfg)
        skid? (and (> (abs (double speed)) 1.0)
                   (> (abs (- (double rear-wheel-speed) (double speed))) 1.0))]
    (set! (.-innerHTML el)
          (str "<b>" (Math/round (* 3.6 (double speed))) "</b> km/h"
               (when-not upright? " &nbsp;<em>down &mdash; R to reset</em>")
               (when (and upright? (false? (:attached? rider)))
                 " &nbsp;<em>rider thrown &mdash; R to reset</em>")
               (when skid? " &nbsp;<em>rear skid</em>")
               "<br>lean " (deg lean) "&deg; &nbsp; steer " (deg steer) "&deg;"
               "<br>rear " (bar rear-travel rlo rhi) " front " (bar front-travel flo fhi)
               "<br>throttle " (bar (:throttle inputs) 0 1) " brake " (bar (:brake inputs) 0 1)
               "<br><span class='moto-keys'>W/&uarr; throttle &middot; S/&darr; brake &middot; "
               "A D/&larr; &rarr; lean &middot; B bail &middot; R reset</span>"))))

;; ---------------------------------------------------------------------------

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 55 1 0.1 600)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        hud (doto (js/document.createElement "div")
              (-> .-className (set! "numeric-readout moto-hud")))
        held (atom #{})
        state (atom nil)
        course (THREE/Group.)]
    (set! (.-background scene) (THREE/Color. 0x0b0d15))
    (set! (.-fog scene) (THREE/Fog. 0x0b0d15 40 160))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.appendChild container hud)
    (.add scene (THREE/HemisphereLight. 0xcfd8ff 0x20242e 0.7))
    (let [sun (THREE/DirectionalLight. 0xfff1d8 1.1)]
      (.set (.-position sun) 30 60 20)
      (.add scene sun))
    (let [ground (THREE/Mesh. (THREE/PlaneGeometry. 4000 4000) (mat 0x2c3242 0.95))
          grid (THREE/GridHelper. 200 100 0x4a5570 0x363e52)]
      (.set (.-rotation ground) (- (/ Math/PI 2)) 0 0)
      (.add scene ground)
      (set! (.. grid -position -y) 0.002)
      (.add scene grid)
      (.add scene course)
      (letfn [(build-course! []
                (doseq [c (vec (.-children course))] (.remove course c))
                (doseq [b (rest (obstacles))]
                  (let [[sx sy sz] (:size b)]
                    (.add course (doto (box-mesh [sx sy sz] (:pos b) (mat 0x7d869c 0.8) (:rot b)))))))
              (reset-bike! []
                (build-course!)
                (let [cfg (config)
                      rider? (boolean (.-riderBody controls))
                      sim (m/scene cfg (obstacles) (m/start-pose cfg (.-startSpeed controls))
                                   {:rider? rider?})
                      meshes (bike-meshes cfg (not rider?))
                      people (when rider?
                               (let [{:keys [model pose]} (second (:models sim))]
                                 (rider-meshes model pose)))]
                  (when-let [old (:meshes @state)]
                    (doseq [[k ^js o] old :when (not= k :pivot)] (.remove scene o)))
                  (doseq [[_ ^js o] (:people @state)] (.remove scene o))
                  (doseq [[k ^js o] meshes :when (not= k :pivot)] (.add scene o))
                  (doseq [[_ ^js o] people] (.add scene o))
                  (reset! state {:cfg cfg
                                 :meshes meshes
                                 :people people
                                 :scene sim
                                 :inputs {:throttle 0.0 :brake 0.0 :lean 0.0 :hands? true}})))
              (sync! []
                (let [{:keys [cfg meshes scene people]} @state
                      model (:model (first (:models scene)))
                      pose (m/bike-pose scene)
                      {:keys [q base]} pose
                      frames (ab/poses model q base)
                      root (ab/frame-of model q base -1)]
                  (place! (:frame meshes) root)
                  (place! (:swingarm meshes) (nth frames m/swingarm))
                  (place! (:rear meshes) (nth frames m/rear-wheel))
                  (place! (:bars meshes) (nth frames m/steering))
                  (place! (:slider meshes) (nth frames m/fork))
                  (place! (:front meshes) (nth frames m/front-wheel))
                  (when-let [{rm :model rp :pose} (second (:models scene))]
                    (doseq [{:keys [link body]} (ab/collision-bodies rm (:q rp) (:base rp))]
                      (when-let [^js mesh (get people link)]
                        (place! mesh body))))
                  ;; The shock runs from a mount on the frame to one
                  ;; partway down the swingarm, so its length is the
                  ;; suspension's travel.
                  (let [top (v/add (:pos root) (q/rotate (:rot root) [-0.12 0.12 0.0]))
                        arm (nth frames m/swingarm)
                        bottom (v/add (:pos arm) (q/rotate (:rot arm)
                                                           (v/scale (:swingarm cfg) 0.55)))]
                    (set-rod! (:shock meshes) top bottom))
                  ;; Keep the grid under the bike, snapped so it does not
                  ;; appear to slide with it.
                  (let [[x _ z] (:pos root)]
                    (.set (.-position grid) (* 2.0 (Math/round (/ x 2.0))) 0.002
                          (* 2.0 (Math/round (/ z 2.0))))
                    pose)))
              (chase! [pose dt]
                ;; Behind and above, along the bike's heading but not its
                ;; lean, eased so a bump does not shake the picture.
                (let [{:keys [pos rot]} (:base pose)
                      fwd (q/rotate rot [1.0 0.0 0.0])
                      flat (v/normalize [(nth fwd 0) 0.0 (nth fwd 2)])
                      want (v/add pos (v/add (v/scale flat -4.2) [0.0 1.7 0.0]))
                      p (.-position camera)
                      k (- 1.0 (Math/exp (* -3.0 dt)))]
                  (.set p (+ (.-x p) (* k (- (nth want 0) (.-x p))))
                        (+ (.-y p) (* k (- (nth want 1) (.-y p))))
                        (+ (.-z p) (* k (- (nth want 2) (.-z p)))))
                  (.lookAt camera (nth pos 0) (+ 0.6 (double (nth pos 1))) (nth pos 2))))
              (on-key [down? ^js e]
                (let [^js target (.-target e)
                      tag (some-> target .-tagName)]
                  ;; Keys typed into the controls, or used to walk the
                  ;; demo picker, are not the rider's.
                  (when-not (or (#{"INPUT" "SELECT" "TEXTAREA"} tag)
                                (and target (.-closest target)
                                     (.closest target ".demo-picker, .lil-gui")))
                    (case (.-code e)
                      "KeyR" (when down? (reset-bike!))
                      "KeyB" (when down? (bail!))
                      (when-let [k (keymap (.-code e))]
                        (.preventDefault e)
                        (swap! held (if down? conj disj) k))))))
              (bail! []
                ;; Off the side the bike is leaning toward, or the left if
                ;; it is upright.
                (let [{:keys [cfg scene]} @state
                      lean (double (:lean (m/telemetry cfg (m/bike-pose scene))))]
                  (swap! state update :scene
                         #(m/throw-rider % (* m/bail-impulse (if (neg? lean) -1.0 1.0))))))
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
                        frame-dt (/ 1.0 60.0)
                        n (max 1 (long (.-substeps controls)))
                        h (/ frame-dt n)
                        ;; Not `scene`: that name is the one being drawn.
                        {:keys [cfg inputs] sim :scene} @state
                        pose (m/bike-pose sim)
                        inputs (rider-inputs inputs @held pose
                                             (:speed (m/telemetry cfg pose)) frame-dt)]
                    (swap! state assoc :inputs inputs)
                    (dotimes [_ n]
                      (swap! state update :scene #(m/step cfg % inputs h)))
                    (let [pose (sync!)]
                      (chase! pose frame-dt)
                      (hud! hud cfg pose inputs (:rider (:scene @state))))
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (let [kd (partial on-key true)
              ku (partial on-key false)]
          (set! (.-reset controls) reset-bike!)
          (set! (.-bail controls) bail!)
          (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
          (reset-bike!)
          (.set (.-position camera) -6 2.4 0)
          (let [gui (GUI. #js {:container container})]
            (doto gui
              (.add controls "reset")
              (-> (.add controls "bail") (.name "bail out (B)"))
              (-> (.add controls "riderBody") (.name "ragdoll rider") (.onChange reset-bike!))
              (-> (.add controls "handsOn") (.name "hands on bars"))
              (-> (.add controls "autopilot") (.name "lap until a key"))
              (-> (.add controls "maxLean" 10 55 1) (.name "max lean°"))
              (-> (.add controls "startSpeed" 2 25 1) (.name "start speed m/s")
                  (.onFinishChange reset-bike!))
              (-> (.add controls "kerbs") (.onChange reset-bike!))
              (-> (.add controls "ramp") (.onChange reset-bike!))
              (-> (.add controls "rearRate" 8 60 1) (.name "rear spring kN/m")
                  (.onFinishChange reset-bike!))
              (-> (.add controls "frontRate" 8 60 1) (.name "front spring kN/m")
                  (.onFinishChange reset-bike!))
              (-> (.add controls "damping" 0.2 3.0 0.1) (.name "damping ×")
                  (.onFinishChange reset-bike!))
              ;; No fewer than four: at 1/180s the pins holding the rider
              ;; on push the bike along by about 30 newtons that are not
              ;; there, measured against 1/240 and 1/480, which agree.
              (-> (.add controls "substeps" 4 8 1))))
          {:start (fn []
                    (when-not @running?
                      (reset! running? true)
                      (js/window.addEventListener "keydown" kd)
                      (js/window.addEventListener "keyup" ku)
                      (on-resize)
                      (animate)))
           :stop (fn []
                   (reset! running? false)
                   (reset! held #{})
                   (js/window.removeEventListener "keydown" kd)
                   (js/window.removeEventListener "keyup" ku))})))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "motorcycle")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
