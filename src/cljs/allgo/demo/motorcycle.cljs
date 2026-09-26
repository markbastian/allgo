(ns allgo.demo.motorcycle
  "A motorcycle on a plane, ridden with the keyboard.

  Everything that moves is `allgo.simulation.motorcycle`: one articulated
  model -- frame, swingarm, two wheels, steering head and fork, joined by
  exact joints -- on a floor, with the tires' grip solved by
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

  **Ride over the curbs.** The front drops the fork a few centimeters,
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
  ring of curbs, holding about 12 m/s -- so there is something to watch,
  and it is the same rider doing it: the autopilot only chooses the lean
  to ask for. Take over and the ramp is off to the right of the start."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.simulation.courses :as c]
            [allgo.simulation.motorcycle :as m]
            ["lil-gui" :default GUI]
            ["three" :as THREE]))

(def ^:private ^js controls
  #js {:course "Excitebike"
       :riderBody true
       :handsOn true
       :autopilot true
       :maxLean 40
       :substeps 4
       :rearRate 25
       :frontRate 20
       :damping 1.0
       :reset (fn [])
       :bail (fn [])})

;; ---------------------------------------------------------------------------
;; The course

(defn- current-course []
  (or (some #(when (= (.-course controls) (:name %)) %) c/courses)
      (first c/courses)))

(def ^:private kind-colors
  {:ramp 0x8a6a45 :whoop 0x7d5f3e :barrier 0xc9a74a :curb 0xc0c4cc
   :log 0x6b4a2b :rock 0x7d7f86})

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
  "A box of `size` centered at `pos` in its parent, turned by `rot`."
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
  "A tire as the torus the physics uses, a rim, and spokes -- the spokes
  are what show it turning."
  [cfg]
  (let [{:keys [wheel-radius tire-radius]} cfg
        major (- wheel-radius tire-radius)
        g (THREE/Group.)
        tire (THREE/Mesh. (THREE/TorusGeometry. major tire-radius 12 40) (mat 0x15161a 0.9))
        rim (THREE/Mesh. (THREE/TorusGeometry. (- major tire-radius 0.005) 0.012 6 40)
                         (mat 0xb8bcc6 0.3 0.8))
        hub (THREE/Mesh. (THREE/CylinderGeometry. 0.05 0.05 0.12 12) (mat 0x8d93a1 0.4 0.7))]
    (.set (.-rotation hub) (/ Math/PI 2) 0 0)
    (.add g tire)
    (.add g rim)
    (.add g hub)
    (dotimes [i 6]
      (let [a (* i (/ Math/PI 3))
            spoke (rod [0.0 0.0 0.0]
                       [(* (- major tire-radius) (Math/cos a)) (* (- major tire-radius) (Math/sin a)) 0.0]
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

(defn- drag-controls
  "Lean, throttle and brake from a drag `dx` `dy` pixels from where it
  began, each at full once the drag reaches `reach`: left and right lean,
  up -- forward, away from the rider -- opens the throttle, and down
  brakes. Proportional, not on and off, which is the thing a finger can
  do that a key cannot. A few pixels either way count as nothing, so a
  tap or a wobble of the thumb does not twitch the bike."
  [dx dy reach]
  (let [dead 6.0
        scale (fn [d] (let [m (- (abs (double d)) dead)]
                        (if (pos? m) (* (Math/sign (double d)) (min 1.0 (/ m (- (double reach) dead)))) 0.0)))
        x (scale dx)
        y (scale dy)]
    ;; Screen y runs down, so up the screen is negative.
    {:lean (- x) :throttle (max 0.0 (- y)) :brake (max 0.0 y)}))

(defn- rider-inputs
  "What the rider is asking for this frame, smoothed toward what the keys
  or a drag -- or, with neither, the course's own autopilot -- want."
  [inputs held drag pose cfg dt]
  (let [max-lean (* (/ Math/PI 180.0) (.-maxLean controls))
        auto (when (and (.-autopilot controls) (empty? held) (nil? drag))
               ((:autopilot (current-course)) pose (m/telemetry cfg pose)
                                              {:cfg cfg :max-lean max-lean}))
        want-throttle (cond (held :throttle) 1.0 drag (:throttle drag) auto (:throttle auto 0.0) :else 0.0)
        want-brake (cond (held :brake) 1.0 drag (:brake drag) auto (:brake auto 0.0) :else 0.0)
        want-lean (cond
                    auto (:lean auto 0.0)
                    (and (held :left) (not (held :right))) max-lean
                    (and (held :right) (not (held :left))) (- max-lean)
                    drag (* max-lean (:lean drag))
                    :else 0.0)]
    (if auto
      ;; The autopilot's asks go straight through. Smoothed like a key
      ;; press, the throttle and brake it flies the bike with in the air
      ;; arrive a quarter of a second late, and the landing is missed.
      {:throttle want-throttle :brake want-brake :lean want-lean
       :hands? (boolean (.-handsOn controls))}
      {:throttle (approach (:throttle inputs) want-throttle 2.5 dt)
       :brake (approach (:brake inputs) want-brake 4.0 dt)
       :lean (approach (:lean inputs) want-lean 1.2 dt)
       :hands? (boolean (.-handsOn controls))})))

;; ---------------------------------------------------------------------------
;; Readout

(defn- hud! [^js el cfg pose inputs rider {:keys [t finished]}]
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
               (when (:finish (current-course))
                 (if finished
                   (str "<br><em>finished in " (.toFixed finished 2) " s</em>")
                   (str "<br>time " (.toFixed (double t) 1) " s")))
               "<br>lean " (deg lean) "&deg; &nbsp; steer " (deg steer) "&deg;"
               "<br>rear " (bar rear-travel rlo rhi) " front " (bar front-travel flo fhi)
               "<br>throttle " (bar (:throttle inputs) 0 1) " brake " (bar (:brake inputs) 0 1)
               "<br><span class='moto-keys'>W/&uarr; throttle &middot; S/&darr; brake &middot; "
               "A D/&larr; &rarr; lean &middot; B bail &middot; R reset"
               "<br>or drag: &larr; &rarr; lean &middot; &uarr; throttle &middot; &darr; brake</span>"))))

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
        ;; Lean, throttle and brake from a mouse or finger held down on
        ;; the scene, or nil when nothing is.
        drag (atom nil)
        ;; Where a drag began and where it has got to: a ring and a dot.
        ring (doto (js/document.createElement "div")
               (-> .-className (set! "moto-drag")))
        state (atom nil)
        course (THREE/Group.)]
    (set! (.-background scene) (THREE/Color. 0x0b0d15))
    (set! (.-fog scene) (THREE/Fog. 0x0b0d15 40 160))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.appendChild container hud)
    (.appendChild container ring)
    ;; Mouse and touch alike, through pointer events: press on the scene
    ;; and drag. Left and right lean, up opens the throttle, down brakes,
    ;; each in proportion to how far; let go and the bike is left to
    ;; the keys, or to itself. The ring shows how far is full.
    (let [^js canvas (.-domElement renderer)
          ^js dot (doto (js/document.createElement "div") (-> .-className (set! "moto-drag__dot")))
          origin (atom nil)
          reach (fn [] (max 60.0 (* 0.18 (min (.-clientWidth container) (.-clientHeight container)))))
          local (fn [^js e]
                  (let [^js r (.getBoundingClientRect container)]
                    [(- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r))]))
          show! (fn [[x0 y0] [x y]]
                  (let [r (reach)
                        dx (- x x0) dy (- y y0)
                        d (Math/hypot dx dy)
                        k (if (> d r) (/ r d) 1.0)
                        ^js s (.-style ring)]
                    (set! (.-display s) "block")
                    (set! (.-left s) (str x0 "px"))
                    (set! (.-top s) (str y0 "px"))
                    (set! (.-width s) (str (* 2 r) "px"))
                    (set! (.-height s) (str (* 2 r) "px"))
                    (set! (.. dot -style -transform)
                          (str "translate(" (* k dx) "px," (* k dy) "px)"))))
          end! (fn [^js e]
                 (when (= (:id @origin) (.-pointerId e))
                   (reset! origin nil)
                   (reset! drag nil)
                   (set! (.. ring -style -display) "none")))]
      (.appendChild ring dot)
      ;; Without this a finger dragged on the scene scrolls the page.
      (set! (.. canvas -style -touchAction) "none")
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (when (and (nil? @origin) (zero? (.-button e)))
                             (.preventDefault e)
                             ;; So the drag goes on when the finger slides off
                             ;; the scene. Some browsers refuse it; the drag
                             ;; still works, just not past the edge.
                             (try (.setPointerCapture canvas (.-pointerId e)) (catch :default _ nil))
                             (let [p (local e)]
                               (reset! origin {:id (.-pointerId e) :at p})
                               (reset! drag (drag-controls 0.0 0.0 (reach)))
                               (show! p p)))))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (when-let [{:keys [id at]} @origin]
                             (when (= id (.-pointerId e))
                               (let [[x0 y0] at
                                     [x y :as p] (local e)]
                                 (reset! drag (drag-controls (- x x0) (- y y0) (reach)))
                                 (show! at p))))))
      (.addEventListener canvas "pointerup" end!)
      (.addEventListener canvas "pointercancel" end!))
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
                (doseq [^js c (vec (.-children course))]
                  (.remove course c)
                  (.dispose (.-geometry c)))
                (let [crs (current-course)]
                  (doseq [b ((:obstacles crs))]
                    (let [[sx sy sz] (:size b)]
                      (.add course (box-mesh [sx sy sz] (:pos b)
                                             (mat (kind-colors (:kind b) 0x7d869c) 0.85)
                                             (:rot b)))))
                  ;; A course with a finish is a lane: dirt down it, and
                  ;; lines across at the start and the end.
                  (when-let [finish (:finish crs)]
                    (let [width (if (= :excitebike (:id crs)) 7.0 4.0)
                          len (+ (double finish) 30.0)]
                      (.add course (box-mesh [len 0.01 width] [(- (* 0.5 len) 10.0) 0.004 0.0]
                                             (mat 0x5a4632 1.0)))
                      (doseq [x [0.0 finish]]
                        (.add course (box-mesh [0.4 0.012 width] [x 0.006 0.0]
                                               (mat 0xe8e8ee 0.6))))))))
              (reset-bike! []
                (build-course!)
                (let [cfg (config)
                      rider? (boolean (.-riderBody controls))
                      sim (c/scene cfg (current-course) {:rider? rider?})
                      meshes (bike-meshes cfg (not rider?))
                      people (when rider?
                               (let [{:keys [model pose]} (second (:models sim))]
                                 (rider-meshes model pose)))]
                  (when-let [old (:meshes @state)]
                    (doseq [[k ^js o] old :when (not= k :pivot)] (.remove scene o)))
                  (doseq [[_ ^js o] (:people @state)] (.remove scene o))
                  (doseq [^js o (:balls @state)] (.remove scene o))
                  (doseq [[k ^js o] meshes :when (not= k :pivot)] (.add scene o))
                  (doseq [[_ ^js o] people] (.add scene o))
                  (reset! state {:cfg cfg
                                 :meshes meshes
                                 :people people
                                 :scene sim
                                 :inputs {:throttle 0.0 :brake 0.0 :lean 0.0 :hands? true}
                                 ;; The clock, when the course finished,
                                 ;; and since when the rider has been off.
                                 :t 0.0 :finished nil :off-since nil})))
              (sync! []
                (let [{:keys [cfg meshes people] sim :scene} @state
                      model (:model (first (:models sim)))
                      pose (m/bike-pose sim)
                      {:keys [q base]} pose
                      frames (ab/poses model q base)
                      root (ab/frame-of model q base -1)]
                  (place! (:frame meshes) root)
                  (place! (:swingarm meshes) (nth frames m/swingarm))
                  (place! (:rear meshes) (nth frames m/rear-wheel))
                  (place! (:bars meshes) (nth frames m/steering))
                  (place! (:slider meshes) (nth frames m/fork))
                  (place! (:front meshes) (nth frames m/front-wheel))
                  (when-let [{rm :model rp :pose} (second (:models sim))]
                    (doseq [{:keys [link body]} (ab/collision-bodies rm (:q rp) (:base rp))]
                      (when-let [^js mesh (get people link)]
                        (place! mesh body))))
                  ;; Cannonballs, a mesh each, made the first time one is
                  ;; seen. They are the only loose bodies in the scene.
                  (let [balls (filterv #(= :cannonball (:kind %)) (:bodies sim))]
                    (while (< (count (:balls @state)) (count balls))
                      (let [^js m (THREE/Mesh. (THREE/SphereGeometry. (:radius (first balls)) 20 14)
                                               (mat 0x2a2c33 0.35))]
                        (.add scene m)
                        (swap! state update :balls (fnil conj []) m)))
                    (doseq [[^js m b] (map vector (:balls @state) balls)]
                      (place! m b)))
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
                      ;; Not in the key list on purpose.
                      "KeyC" (when down?
                               (swap! state (fn [{:keys [cfg] :as st}]
                                              (update st :scene #(m/fire-at cfg %)))))
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
                        {:keys [cfg inputs t finished off-since] sim :scene} @state
                        pose (m/bike-pose sim)
                        inputs (rider-inputs inputs @held @drag pose cfg frame-dt)
                        t (+ (double t) frame-dt)
                        crs (current-course)
                        finished (or finished (when (c/finished? crs pose) t))
                        off? (or (false? (get-in sim [:rider :attached?]))
                                 (not (:upright? (m/telemetry cfg pose))))
                        off-since (when off? (or off-since t))]
                    (swap! state assoc :inputs inputs :t t :finished finished :off-since off-since)
                    (dotimes [_ n]
                      (swap! state update :scene #(m/step cfg % inputs h)))
                    (let [pose (sync!)]
                      (chase! pose frame-dt)
                      (hud! hud cfg pose inputs (:rider (:scene @state)) @state))
                    ;; Riding itself, it goes round again: a few seconds
                    ;; after the finish, or after coming off.
                    (when (and (.-autopilot controls) (empty? @held) (nil? @drag)
                               (or (and finished (> (- t finished) 3.0))
                                   (and off-since (> (- t off-since) 4.0))))
                      (reset-bike!))
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
              (-> (.add controls "course" (clj->js (mapv :name c/courses)))
                  (.onChange reset-bike!))
              (-> (.add controls "autopilot") (.name "ride itself until a key"))
              (-> (.add controls "maxLean" 10 55 1) (.name "max lean°"))
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
