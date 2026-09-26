(ns allgo.demo.cloth
  "WebGL demo for cloth, after Ten Minute Physics 14.

  Cloth is `allgo.physics.xpbd` with no new machinery at all: a triangle
  sheet held by two sets of distance constraints. One set holds the edges,
  which is what stops it stretching; the other holds the pairs of vertices
  opposite each shared edge, which is what stops it folding. `bending` is
  the compliance of the second set, and it is the whole difference between
  card and silk.

  `pinning` decides what the sheet is attached to, and is worth trying
  against `bending`: a hanging sheet is shaped mostly by gravity whatever
  its stiffness, while a sheet held along one edge shows the difference
  plainly.

  `selfCollision` is Ten Minute Physics 15, and shows best under `squeeze`,
  which drags the two pinned edges together until the sheet has nowhere to
  go but into itself. Without it the folds pass straight through one
  another; with it they stack."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.tri-mesh :as tri]
            [allgo.physics.xpbd :as xpbd]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:resolution        24
       :bending           0.05
       :stretch           0.0
       :pinning           "none"
       :gravity           9.8
       :damping           0.4
       :substeps          12
       :obstacle          true
       :friction          0.35
       :selfCollision     true
       :squeeze           false
       :showWireframe     false})

(def ^:private sheet-size 1.6)
(def ^:private ball-radius 0.32)
(def ^:private ball-center [0.0 0.55 0.0])

(defn- pinned-ids
  "Which particles are held, for the current pinning mode."
  [mesh n]
  (let [row (vec (range (inc n)))]                       ; the j = 0 edge
    (case (.-pinning controls)
      "none"        []
      "one corner"  [(first (tri/corner-vertices mesh))]
      "two corners" (take 2 (tri/corner-vertices mesh))
      "four corners" (tri/corner-vertices mesh)
      "one edge"    row
      "two edges"   (concat row (map #(+ % (* (inc n) n)) (range (inc n))))
      [])))

(defn- build! [^js scene state]
  (doseq [k [:mesh-obj :wire]]
    (when-let [^js old (k @state)]
      (.remove scene old)
      (.dispose (.-geometry old))
      (.dispose (.-material old))))
  (let [n     (.-resolution controls)
        mesh  (tri/grid n n (/ sheet-size n) 1.25)
        body  (xpbd/cloth mesh {:stretch-compliance (.-stretch controls)
                                :bend-compliance    (.-bending controls)})
        _     (xpbd/pin! body (pinned-ids mesh n))
        ;; A constraint, so contact is enforced on every substep. Applied
        ;; once a frame instead, the sheet tunnels straight through.
        thickness (* 0.7 (/ sheet-size n))
        body  (cond-> body
                (.-obstacle controls)
                (xpbd/add-constraint
                 (xpbd/sphere-constraint ball-center ball-radius (.-friction controls)))
                (.-selfCollision controls)
                ;; Detected once a frame, resolved every substep. Needs the
                ;; speed limit too, or a fast fold crosses its own thickness
                ;; between detections.
                (xpbd/add-constraint
                 (xpbd/self-collision-constraint (xpbd/cloth mesh {}) thickness 0.15)))
        ;; WebGL has no double-precision attribute, so the solver's array
        ;; is mirrored into a single-precision one each frame.
        render (js/Float32Array. (alength (:pos body)))
        geom  (doto (THREE/BufferGeometry.)
                (.setAttribute "position" (THREE/BufferAttribute. render 3))
                (.setIndex (clj->js (:tri-ids mesh))))
        mat   (THREE/MeshPhongMaterial. #js {:color 0xe0a85c :side THREE/DoubleSide
                                             :flatShading false})
        obj   (THREE/Mesh. geom mat)
        wire  (THREE/LineSegments.
               (doto (THREE/BufferGeometry.)
                 (.setAttribute "position" (THREE/BufferAttribute. render 3))
                 (.setIndex (clj->js (:edge-ids mesh))))
               (THREE/LineBasicMaterial. #js {:color 0x4a7ab0}))]
    (.add scene obj)
    (.add scene wire)
    (swap! state assoc :mesh mesh :body body :render render
           :mesh-obj obj :wire wire :geom geom :thickness thickness
           :pinned (vec (pinned-ids mesh n)) :frame 0
           :tris (quot (count (:tri-ids mesh)) 3))))

(defn- refresh! [{:keys [^js geom ^js wire ^js render body]}]
  (.set render (:pos body))
  (set! (.-visible wire) (.-showWireframe controls))
  (.computeVertexNormals geom)
  (set! (.. geom -attributes -position -needsUpdate) true)
  (set! (.. wire -geometry -attributes -position -needsUpdate) true)
  (.computeBoundingSphere geom))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.75))
    (let [sun (THREE/DirectionalLight. 0xffffff 1.0)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 20 20)
                             (THREE/MeshPhongMaterial. #js {:color 0x141a2a}))]
      (set! (.-x (.-rotation floor)) (- (/ js/Math.PI 2)))
      (.add scene floor))
    (.add scene (doto (THREE/GridHelper. 6 12 0x2b3450 0x1b2236)
                  (-> .-position (.setY 0.002))))
    (let [ball (THREE/Mesh. (THREE/SphereGeometry. ball-radius 32 24)
                            (THREE/MeshPhongMaterial. #js {:color 0x4f7ac0}))]
      (.set (.-position ball) (nth ball-center 0) (nth ball-center 1) (nth ball-center 2))
      (.add scene ball)
      (.set (.-position camera) 1.9 1.5 2.3)
      (.appendChild container (.-domElement renderer))
      (let [orbit (OrbitControls. camera (.-domElement renderer))]
        (set! (.-enableDamping orbit) true)
        (.set (.-target orbit) 0 0.7 0)
        (letfn [(rebuild! []
                  (set! (.-visible ball) (.-obstacle controls))
                  (build! scene state)
                  (refresh! @state))
                (resize! []
                  (let [w (.-clientWidth container) h (.-clientHeight container)]
                    (when (and (pos? w) (pos? h))
                      (.setSize renderer w h)
                      (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                      (set! (.-aspect camera) (/ w h))
                      (.updateProjectionMatrix camera))))
                (squeeze! []
                  ;; Drag the two pinned edges together, which is the only
                  ;; way to make a sheet meet itself hard enough to see the
                  ;; difference self-collision makes.
                  (let [{:keys [body frame]} @state
                        n (.-resolution controls)
                        row-a (range (inc n))
                        row-b (map #(+ % (* (inc n) n)) (range (inc n)))
                        t (min 1.0 (/ frame 180.0))
                        z (* 0.5 sheet-size (- 1.0 (* 0.95 t)))]
                    (doseq [i row-a] (let [[px py _] (xpbd/particle body i)]
                                       (xpbd/set-particle! body i [px py (- z)])))
                    (doseq [i row-b] (let [[px py _] (xpbd/particle body i)]
                                       (xpbd/set-particle! body i [px py z])))
                    (swap! state update :frame inc)))
                (tick []
                  (when @running?
                    (js/requestAnimationFrame tick)
                    (let [t0   (js/performance.now)
                          body (:body @state)]
                      (when (.-squeeze controls) (squeeze!))
                      (xpbd/step! body {:gravity  [0.0 (- (.-gravity controls)) 0.0]
                                        :substeps (.-substeps controls)
                                        :damping  (.-damping controls)
                                        :max-velocity (when (.-selfCollision controls)
                                                        (xpbd/speed-limit (:thickness @state)
                                                                          (/ 1.0 60.0)))
                                        :floor    0.0})
                      (refresh! @state)
                      (.update orbit)
                      (.render renderer scene camera)
                      (tick-fps! (- (js/performance.now) t0)))))]
          (.observe (js/ResizeObserver. resize!) container)
          (resize!)
          (rebuild!)
          (let [gui (GUI. #js {:container container})]
            (doto gui
              (-> (.add controls "resolution" 8 48 2) (.onFinishChange rebuild!))
              (-> (.add controls "pinning" #js ["none" "two corners" "one edge" "two edges"
                                                "four corners" "one corner"])
                  (.onChange rebuild!))
              (-> (.add controls "bending" 0 0.5 0.005) (.onChange rebuild!))
              (-> (.add controls "stretch" 0 0.01 0.0005) (.onChange rebuild!))
              (.add controls "gravity" 0 20 0.5)
              (.add controls "damping" 0 3 0.1)
              (.add controls "substeps" 1 25 1)
              (-> (.add controls "obstacle") (.onChange rebuild!))
              (-> (.add controls "friction" 0 1 0.05) (.onChange rebuild!))
              (-> (.add controls "selfCollision") (.onChange rebuild!))
              (-> (.add controls "squeeze")
                  (.onChange (fn [on]
                               ;; Squeezing needs the two edges held.
                               (when on (set! (.-pinning controls) "two edges"))
                               (rebuild!))))
              (.add controls "showWireframe")
              (.add #js {:reset rebuild!} "reset")))
          {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
           :stop  (fn [] (reset! running? false))})))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "cloth")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
