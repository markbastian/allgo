(ns allgo.demo.human-arm
  "A seven-jointed arm, and the one number that says which way the elbow
  went.

  The six-axis robot next door reaches a pose in exactly eight ways. This
  one reaches it in infinitely many, and the difference is a joint: a ball
  at the shoulder, a hinge at the elbow, a ball at the wrist is seven
  degrees of freedom for a hand pose that needs only six. The spare one
  has to go somewhere, and where it goes is the elbow.

  Drag `swivel` and watch. The hand does not move. It cannot -- every
  swivel angle is an exact solution, and the ring drawn through the elbow
  is the set of places the elbow is allowed to be. Try it on yourself:
  put your palm flat on the desk and swing your elbow out. Same hand, same
  table, different arm.

  Two things worth noticing while you drag.

  The elbow's *bend* never changes. Only how far apart your shoulder and
  wrist are decides that, and swinging the elbow round does not move
  either of them. The readout shows it holding still while everything
  else moves.

  And the shoulder is not free, despite having three degrees of freedom.
  Two of them aim the upper arm at the elbow. The third is spent keeping
  the elbow hinge across the plane of the arm, because a hinge cannot bend
  sideways out of its own plane. Pick the swivel and the entire arm
  follows exactly -- there is nothing left over to choose.

  `sweep` runs the swivel round on its own, which is the clearest way to
  see a hand held still by an arm that will not stop moving."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as quat]
            [allgo.geometry.vec3 :as v]
            [allgo.kinematics.arm :as arm]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:swivel     0.0
       :sweep      true
       :speed      0.5
       :ghosts     9
       :showCircle true
       :reach      0.62
       :height     0.30})

(def ^:private the-arm (arm/arm {:upper 0.30 :forearm 0.27}))

(defn- target-pose
  "Where the hand is being held. Reach and height place it; the hand keeps
  a fixed orientation so that the only thing free is the elbow."
  []
  (let [[lo hi] (arm/span the-arm)
        d (max (+ lo 0.02) (min (- hi 0.02) (* (.-reach controls) hi)))
        up (.-height controls)
        dir (v/normalize [1.0 (* 2.0 (- up 0.5)) 0.25])]
    {:pos (v/scale dir d)
     :rot (quat/from-vectors [1.0 0.0 0.0] dir)}))

;; ---------------------------------------------------------------------------
;; Drawing

(defn- bone!
  "Points a capsule from `a` to `b`."
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

(defn- make-limb [^js scene colour opacity radius]
  (let [mat (THREE/MeshPhongMaterial.
             #js {:color colour :transparent (< opacity 1.0) :opacity opacity})
        bone (fn [] (let [m (THREE/Mesh. (THREE/CylinderGeometry. radius radius 1.0 14) mat)]
                      (.add scene m) m))
        blob (fn [r] (let [m (THREE/Mesh. (THREE/SphereGeometry. r 16 12) mat)]
                       (.add scene m) m))]
    {:upper (bone) :forearm (bone)
     :shoulder (blob (* 1.5 radius)) :elbow (blob (* 1.3 radius)) :hand (blob (* 1.2 radius))}))

(defn- draw-limb! [{:keys [upper forearm shoulder elbow hand]} p]
  (let [{s :shoulder e :elbow-position w :wrist} p]
    (bone! upper s e)
    (bone! forearm e w)
    (doseq [[^js m at] [[shoulder s] [elbow e] [hand w]]]
      (set! (.-visible m) true)
      (.set (.-position m) (nth at 0) (nth at 1) (nth at 2)))))

(defn- hide-limb! [{:keys [upper forearm shoulder elbow hand]}]
  (doseq [^js m [upper forearm shoulder elbow hand]] (set! (.-visible m) false)))

(defn- draw-circle!
  "The ring of places the elbow may be."
  [^js line circle]
  (if-not circle
    (set! (.-visible line) false)
    (let [{:keys [centre radius zero ninety]} circle
          pts (for [i (range 65)]
                (let [t (* 2 js/Math.PI (/ i 64))]
                  (v/add centre (v/scale (v/add (v/scale zero (js/Math.cos t))
                                                (v/scale ninety (js/Math.sin t)))
                                         radius))))
          arr (js/Float32Array. (* 3 65))]
      (doseq [[i p] (map-indexed vector pts)]
        (aset arr (* 3 i) (nth p 0))
        (aset arr (+ (* 3 i) 1) (nth p 1))
        (aset arr (+ (* 3 i) 2) (nth p 2)))
      (set! (.-visible line) true)
      (.setAttribute (.-geometry line) "position" (THREE/BufferAttribute. arr 3))
      (set! (.. line -geometry -attributes -position -needsUpdate) true))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        overlay  (js/document.createElement "div")
        state    (atom {:phase 0.0})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.85))
    (let [sun (THREE/DirectionalLight. 0xffffff 0.9)]
      (.set (.-position sun) 2 3 2)
      (.add scene sun))
    (.add scene (doto (THREE/GridHelper. 1.6 16 0x2b3450 0x1b2236)
                  (-> .-position (.setY -0.5))))
    (set! (.-cssText (.-style overlay))
          (str "position:absolute;left:14px;top:14px;z-index:5;pointer-events:none;"
               "font:12px ui-monospace,monospace;color:#dfe7f5;background:rgba(0,0,0,0.55);"
               "padding:8px 10px;border-radius:8px;white-space:pre;line-height:1.5"))
    (.appendChild container overlay)
    (set! (.-position (.-style container)) "relative")
    (.set (.-position camera) 0.62 0.42 0.72)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))
          ;; The hand, drawn where it is being held so that it is visibly
          ;; the thing not moving.
          held (THREE/Mesh. (THREE/SphereGeometry. 0.028 20 14)
                            (THREE/MeshPhongMaterial. #js {:color 0xf2b134}))
          ring (THREE/Line. (THREE/BufferGeometry.)
                            (THREE/LineBasicMaterial. #js {:color 0x6fd0ff}))]
      (.add scene held)
      (.add scene ring)
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0.18 0.0 0.0)
      (letfn [(resize! []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (.setSize renderer w h)
                    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera))))
              (tick []
                (when @running?
                  (js/requestAnimationFrame tick)
                  (let [t0 (js/performance.now)
                        {:keys [limb ghosts phase]} @state
                        phase (if (.-sweep controls)
                                (+ phase (* 0.016 (.-speed controls)))
                                phase)
                        psi (if (.-sweep controls)
                              (* 2 js/Math.PI (mod phase 1.0))
                              (.-swivel controls))
                        target (target-pose)
                        circle (arm/elbow-circle the-arm (:pos target))
                        solution (arm/solve the-arm target psi)]
                    (swap! state assoc :phase phase)
                    (when (.-sweep controls)
                      (set! (.-swivel controls) psi))
                    (let [[hx hy hz] (:pos target)]
                      (.set (.-position held) hx hy hz))
                    (set! (.-visible ring) (boolean (and circle (.-showCircle controls))))
                    (when (.-showCircle controls) (draw-circle! ring circle))
                    ;; Ghosts at evenly spaced swivels: the family, sampled.
                    (let [samples (if (pos? (.-ghosts controls))
                                    (arm/swivel-samples the-arm target (.-ghosts controls))
                                    [])]
                      (doseq [[i g] (map-indexed vector ghosts)]
                        (if (< i (count samples))
                          (draw-limb! g (arm/pose the-arm (nth samples i)))
                          (hide-limb! g))))
                    (if solution
                      (let [p (arm/pose the-arm solution)]
                        (draw-limb! limb p)
                        (set! (.-textContent overlay)
                              (str "swivel " (.toFixed psi 2) " rad\n"
                                   "elbow bend " (.toFixed (:elbow solution) 3)
                                   " rad -- unchanged by the swivel\n"
                                   "hand held at "
                                   (.toFixed (v/distance (:pos target) [0 0 0]) 3)
                                   " m from the shoulder")))
                      (do (hide-limb! limb)
                          (set! (.-textContent overlay) "out of reach")))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. resize!) container)
        (resize!)
        (swap! state assoc
               :limb (make-limb scene 0xd9a15c 1.0 0.026)
               :ghosts (vec (repeatedly 16 #(make-limb scene 0x7fc4ff 0.16 0.012))))
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (.add controls "sweep")
            (.add controls "swivel" (- js/Math.PI) js/Math.PI 0.01)
            (.add controls "speed" 0.05 2 0.05)
            (.add controls "ghosts" 0 16 1)
            (.add controls "showCircle")
            (.add controls "reach" 0.25 0.98 0.01)
            (.add controls "height" 0.0 1.0 0.01)))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "human-arm")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
