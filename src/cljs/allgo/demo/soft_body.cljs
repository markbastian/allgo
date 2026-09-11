(ns allgo.demo.soft-body
  "WebGL demo for `allgo.physics.xpbd`, after Ten Minute Physics 10.

  Tetrahedral bodies dropped on a floor and solved by XPBD. The two
  compliance sliders are the interesting controls: edge compliance is how
  far the material will stretch, volume compliance how far it can be
  squashed, and zero on both is as rigid as the solver gets. `squash`
  flattens every particle onto the floor plane, which is the demonstration
  the tutorial is named for -- a body destroyed that thoroughly is back in
  shape within a second, because nothing in the method can diverge.

  The solver keeps positions in a flat array laid out exactly as a WebGL
  position buffer, but at double precision, which WebGL has no attribute
  format for -- handing three.js the solver's own array fails with
  `Unsupported buffer data format`. So the geometry owns a single-precision
  mirror and each frame is one native `TypedArray.set` over it, which
  converts as it copies."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.tet-mesh :as tet]
            [allgo.physics.xpbd :as xpbd]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:bodies            3
       :resolution        3
       :shape             "ball"
       :substeps          10
       :edgeCompliance    0.0
       :volumeCompliance  0.0
       :gravity           10.0
       :damping           0.0
       :wireframe         false})

(def ^:private colors [0xe0a85c 0x7ed4a0 0x6f9fe0 0xd47ea0 0xa0d47e])

(defn- make-mesh
  "A tetrahedral body, optionally rounded off into a ball."
  [n]
  (let [box (tet/lattice-box n n n (/ 1.0 n))]
    (if (= "ball" (.-shape controls))
      (tet/deform box (tet/sphere-deformation [0.0 0.5 0.0] 0.5))
      box)))

(defn- refresh-geometry!
  "Copy the solved positions into the render buffer and re-light the skin."
  [{:keys [body ^js geometry ^js render-pos]}]
  (.set render-pos (:pos body))
  (.computeVertexNormals geometry)
  (set! (.. geometry -attributes -position -needsUpdate) true)
  (.computeBoundingSphere geometry))

(defn- build-body
  "One soft body plus the three.js mesh that draws its skin."
  [i scene]
  (let [mesh  (tet/translate (make-mesh (.-resolution controls))
                             ;; Spread across the floor and staggered only
                             ;; slightly in height, so they all land in
                             ;; shot and at different moments.
                             [(* 1.1 (- i (/ (dec (.-bodies controls)) 2.0)))
                              (+ 0.8 (* 0.45 i))
                              (* 0.35 (- i (/ (dec (.-bodies controls)) 2.0)))])
        body  (xpbd/soft-body mesh {:edge-compliance   (.-edgeCompliance controls)
                                    :volume-compliance (.-volumeCompliance controls)})
        ;; Same layout as the solver's array, one precision down.
        render-pos (js/Float32Array. (alength (:pos body)))
        geometry (doto (THREE/BufferGeometry.)
                   (.setAttribute "position" (THREE/BufferAttribute. render-pos 3))
                   (.setIndex (clj->js (:surface-tri-ids mesh))))
        material (THREE/MeshPhongMaterial.
                  #js {:color (nth colors (mod i (count colors)))
                       :flatShading true
                       :wireframe (.-wireframe controls)})
        surface  (THREE/Mesh. geometry material)]
    (.add scene surface)
    (let [entry {:body body :mesh mesh :surface surface
                 :geometry geometry :render-pos render-pos}]
      (refresh-geometry! entry)
      entry)))

(defn- world []
  {:gravity  [0.0 (- (.-gravity controls)) 0.0]
   :substeps (.-substeps controls)
   :damping  (.-damping controls)
   :floor    0.0})

;; ---------------------------------------------------------------------------

(defn- build-scene [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 500)
        renderer (THREE/WebGLRenderer. #js {:antialias true})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.7))
    (let [sun (THREE/DirectionalLight. 0xffffff 1.1)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    ;; A floor to land on, matching the solver's y = 0 plane.
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 20 20)
                             (THREE/MeshPhongMaterial. #js {:color 0x141a2a}))]
      (set! (.-x (.-rotation floor)) (- (/ js/Math.PI 2)))
      (.add scene floor))
    (let [grid (THREE/GridHelper. 8 16 0x2b3450 0x1b2236)]
      (set! (.-y (.-position grid)) 0.002)
      (.add scene grid))
    (.set (.-position camera) 2.6 2.2 4.2)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0.0 0.45 0.0)
      (.update orbit)
      {:scene scene :camera camera :renderer renderer :orbit orbit})))

(defn init! [^js container]
  (let [{:keys [^js scene ^js camera ^js renderer ^js orbit]} (build-scene container)
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom {:bodies []})]
    (letfn [(clear! []
              (doseq [{:keys [^js surface ^js geometry]} (:bodies @state)]
                (.remove scene surface)
                (.dispose geometry)
                (.dispose (.-material surface)))
              (swap! state assoc :bodies []))
            (rebuild! []
              (clear!)
              (swap! state assoc :bodies
                     (mapv #(build-body % scene) (range (.-bodies controls)))))
            (squash! []
              ;; The tutorial's party trick: destroy the body completely and
              ;; watch it recover. Nothing here is integrated, so there is
              ;; no state left to be invalid.
              (doseq [{:keys [body]} (:bodies @state)]
                (dotimes [i (:n body)]
                  (let [[px _ pz] (xpbd/particle body i)]
                    (xpbd/set-particle! body i [px 0.01 pz])))))
            (resize! []
              (let [w (.-clientWidth container)
                    h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (.setSize renderer w h)
                  (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                  (set! (.-aspect camera) (/ w h))
                  (.updateProjectionMatrix camera))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      w  (world)]
                  (doseq [{:keys [body] :as entry} (:bodies @state)]
                    (xpbd/step! body w)
                    (refresh-geometry! entry))
                  (.update orbit)
                  (.render renderer scene camera)
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (rebuild!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "bodies" 1 5 1) (.onFinishChange rebuild!))
          (-> (.add controls "resolution" 2 6 1) (.onFinishChange rebuild!))
          (-> (.add controls "shape" #js ["ball" "box"]) (.onChange rebuild!))
          (.add controls "substeps" 1 25 1)
          (-> (.add controls "edgeCompliance" 0 0.005 0.0001) (.onChange rebuild!))
          (-> (.add controls "volumeCompliance" 0 0.005 0.0001) (.onChange rebuild!))
          (.add controls "gravity" 0 30 1)
          (.add controls "damping" 0 5 0.1)
          (-> (.add controls "wireframe") (.onChange rebuild!))
          (.add #js {:squash squash!} "squash")
          (.add #js {:drop rebuild!} "drop")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "soft-body")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
