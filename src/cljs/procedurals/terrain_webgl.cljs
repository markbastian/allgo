(ns procedurals.terrain-webgl
  (:require [procedurals.fps :as fps]
            [procedurals.mesh :as mesh]
            [procedurals.terrain :as terrain]
            [procedurals.tin :as tin]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def config
  {:generator    :diamond-square
   :iterations   8
   :width        1.0
   :tin-points   500
   :dim          129
   :height-scale 200
   :cell-scale   8
   :noise-scale  180.0})

(defn- set-vec3! [^js arr idx x y z]
  (let [base (* idx 3)]
    (aset arr base x)
    (aset arr (+ base 1) y)
    (aset arr (+ base 2) z)))

(defn- diamond-square-grid [{:keys [iterations width]}]
  (let [grid-data (terrain/generate {:width width :iterations iterations
                                     :corners [0.0 (rand) (rand) (rand)]})]
    {:grid (terrain/cells->grid grid-data) :dim (:dim grid-data)}))

(defn- tin-grid [{:keys [tin-points dim smooth-passes] :or {tin-points 500 dim 129 smooth-passes 3}}]
  {:grid (-> (tin/generate {:n tin-points :size 1.0}) (tin/sample-grid dim) (tin/smooth-grid smooth-passes))
   :dim  dim})

(defn build-terrain-data
  "Builds a heightmap via `:generator` (`:diamond-square`, the default
  regular-grid midpoint-displacement fractal, or `:tin`, a Delaunay-
  triangulated irregular network per `procedurals.tin`)."
  [{:keys [generator] :or {generator :diamond-square} :as config}]
  (let [{:keys [grid dim]} (case generator
                             :diamond-square (diamond-square-grid config)
                             :tin (tin-grid config))
        hs (flatten grid)]
    {:grid grid
     :dim  dim
     :lo   (apply min hs)
     :hi   (apply max hs)}))

(defn build-geometry
  "Bakes the heightmap into an indexed THREE.BufferGeometry: one vertex per
  grid cell (shared across its adjacent triangles) so normals interpolate
  smoothly, matching `procedurals.terrain-shape/build-mesh-shape` but using
  the standard y-up world (+height, not the Quil renderer's -height)."
  [{:keys [grid dim lo hi]} {:keys [cell-scale height-scale noise-scale]}]
  (let [n          (* dim dim)
        offset     (* 0.5 (dec dim) cell-scale)
        positions  (js/Float32Array. (* n 3))
        normals    (js/Float32Array. (* n 3))
        colors     (js/Float32Array. (* n 3))
        slopes     (for [i (range dim) j (range dim)] (mesh/slope-at grid dim i j cell-scale height-scale))
        max-slope  (max 1e-6 (apply max slopes))]
    (doseq [i (range dim) j (range dim)]
      (let [idx        (+ (* i dim) j)
            h          (get-in grid [i j])
            height-t   (if (= lo hi) 0.5 (/ (double (- h lo)) (- hi lo)))
            slope-t    (/ (mesh/slope-at grid dim i j cell-scale height-scale) max-slope)
            wx         (- (* j cell-scale) offset)
            wz         (- (* i cell-scale) offset)
            noise-t    (mesh/noise-at wx wz noise-scale)
            [nx ny nz] (mesh/normal-at grid dim i j cell-scale height-scale 1.0)
            [r g b]    (mesh/terrain-rgb height-t slope-t noise-t)]
        (set-vec3! positions idx wx (* h height-scale) wz)
        (set-vec3! normals idx nx ny nz)
        (set-vec3! colors idx (/ r 255) (/ g 255) (/ b 255))))
    (let [quads (* (dec dim) (dec dim))
          index (js/Uint32Array. (* quads 6))]
      (doseq [i (range (dec dim)) j (range (dec dim))]
        (let [a    (+ (* i dim) j)
              b    (+ (* i dim) (inc j))
              c    (+ (* (inc i) dim) j)
              d    (+ (* (inc i) dim) (inc j))
              base (* 6 (+ (* i (dec dim)) j))]
          (aset index base a) (aset index (+ base 1) c) (aset index (+ base 2) b)
          (aset index (+ base 3) b) (aset index (+ base 4) c) (aset index (+ base 5) d)))
      (doto (THREE/BufferGeometry.)
        (.setAttribute "position" (THREE/BufferAttribute. positions 3))
        (.setAttribute "normal" (THREE/BufferAttribute. normals 3))
        (.setAttribute "color" (THREE/BufferAttribute. colors 3))
        (.setIndex (THREE/BufferAttribute. index 1))))))

(defn build-mesh [config]
  (let [terrain-data (build-terrain-data config)
        geometry     (build-geometry terrain-data config)
        material     (THREE/MeshStandardMaterial. #js {:vertexColors true
                                                       :side         THREE/DoubleSide
                                                       :roughness    0.9
                                                       :metalness    0.0})]
    {:span (* (:cell-scale config) (dec (:dim terrain-data)))
     :mesh (THREE/Mesh. geometry material)}))

(defn- dispose-mesh! [^js terrain-mesh]
  (.dispose (.-geometry terrain-mesh))
  (.dispose (.-material terrain-mesh)))

(defn init! [^js container config]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 1 100000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        {:keys [span mesh]} (build-mesh config)
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {:terrain-mesh mesh :wireframe? false :generator (:generator config) :hovering? false})]
    (set! (.-background scene) (THREE/Color. 0x0f0f19))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (.appendChild container (.-domElement renderer))
    (.add scene mesh)
    (.add scene (THREE/AmbientLight. 0xffffff 0.5))
    (let [sun (THREE/DirectionalLight. 0xfff4e0 1.2)]
      (.set (.-position sun) (* span 0.4) (* span 1.0) (* span 0.3))
      (.add scene sun))
    (.set (.-position camera) (* span 0.9) (* span 0.7) (* span 0.9))
    (let [controls (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping controls) true)
      (.set (.-target controls) 0 0 0)
      (letfn [(apply-wireframe! [^js m]
                (set! (.-wireframe (.-material m)) (:wireframe? @state)))
              (regenerate! []
                (dispose-mesh! (:terrain-mesh @state))
                (.remove scene (:terrain-mesh @state))
                (let [{:keys [mesh]} (build-mesh (assoc config :generator (:generator @state)))]
                  (apply-wireframe! mesh)
                  (.add scene mesh)
                  (swap! state assoc :terrain-mesh mesh)))
              (on-key-down [^js e]
                (when (and (:hovering? @state) (= (.-key e) "r")) (regenerate!)))
              (on-resize []
                ;; A hidden card measures 0x0, which would make the aspect NaN.
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (.update controls)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.addEventListener container "mouseenter" #(swap! state assoc :hovering? true))
        (.addEventListener container "mouseleave" #(swap! state assoc :hovering? false))
        (.addEventListener js/window "keydown" on-key-down)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui       (GUI. #js {:container container})
              gui-state #js {:generator (name (:generator @state)) :wireframe false :regenerate regenerate!}]
          (-> (.add gui gui-state "generator" #js ["tin" "diamond-square"])
              (.onChange (fn [v]
                           (swap! state assoc :generator (keyword v))
                           (regenerate!))))
          (-> (.add gui gui-state "wireframe")
              (.onChange (fn [v]
                           (swap! state assoc :wireframe? v)
                           (apply-wireframe! (:terrain-mesh @state)))))
          (.add gui gui-state "regenerate"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start!
  "Build the scene on first selection, then resume the render loop."
  []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "terrain")]
                     (reset! controller (init! container config))))]
    ((:start c))))

(defn stop!
  "Halt the render loop. The WebGL context is kept -- it costs nothing idle,
  and rebuilding it would throw away the camera the user has framed."
  []
  (when-let [c @controller] ((:stop c))))
