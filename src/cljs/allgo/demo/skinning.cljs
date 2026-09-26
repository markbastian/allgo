(ns allgo.demo.skinning
  "WebGL demo for `allgo.physics.skinning`, after Ten Minute Physics 12.

  A detailed surface carried by a coarse tetrahedral cage. The cage is
  what XPBD simulates; the surface is bound to it once and thereafter
  follows for a weighted sum of four points per vertex. Turn on `showCage`
  to see how little is actually being simulated.

  The caption gives the ratio the technique trades on. A cage of a couple
  of hundred tetrahedra costs the same whether the surface it carries has
  two thousand vertices or thirty thousand, because the surface is not
  simulated at all -- so the detail is very nearly free, and raising
  `detail` changes the picture without changing the physics."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.tet-mesh :as tet]
            [allgo.physics.skinning :as skinning]
            [allgo.physics.xpbd :as xpbd]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:shape            "knot"
       :detail           6
       :cage             3
       :substeps         10
       :edgeCompliance   0.0005
       :volumeCompliance 0.0
       :gravity          10.0
       :showCage         false
       :showSurface      true})

(defn- visual-geometry
  "The detailed surface. `detail` scales its tessellation, not its size."
  []
  (let [d (.-detail controls)]
    (case (.-shape controls)
      "sphere" (THREE/SphereGeometry. 0.42 (* 8 d) (* 6 d))
      "torus"  (THREE/TorusGeometry. 0.34 0.15 (* 6 d) (* 12 d))
      (THREE/TorusKnotGeometry. 0.30 0.10 (* 24 d) (* 6 d)))))

(defn- cage-for
  "A lattice cage enclosing the geometry's bounding box, with a margin so
  the surface sits inside rather than on the boundary."
  [^js geometry n]
  (.computeBoundingBox geometry)
  (let [bb   (.-boundingBox geometry)
        lo   (.-min bb) hi (.-max bb)
        span (max (- (.-x hi) (.-x lo)) (- (.-y hi) (.-y lo)) (- (.-z hi) (.-z lo)))
        size (/ (* 1.25 span) n)]
    (-> (tet/lattice-box n n n size)
        ;; lattice-box is centered in x and z with its base at y=0; move it
        ;; onto the geometry's own center.
        (tet/translate [(/ (+ (.-x lo) (.-x hi)) 2)
                        (- (/ (+ (.-y lo) (.-y hi)) 2) (/ (* n size) 2))
                        (/ (+ (.-z lo) (.-z hi)) 2)]))))

;; ---------------------------------------------------------------------------

(defn- build! [^js scene state]
  (doseq [k [:surface :cage-lines]]
    (when-let [^js old (k @state)]
      (.remove scene old)
      (.dispose (.-geometry old))
      (.dispose (.-material old))))
  (let [geometry (visual-geometry)
        vis      (vec (.. geometry -attributes -position -array))
        cage     (cage-for geometry (.-cage controls))
        ;; Lift the whole thing clear of the floor so it has somewhere to
        ;; fall from.
        drop-to  0.85
        [lo]     (tet/bounds cage)
        cage     (tet/translate cage [0.0 (- drop-to (nth lo 1)) 0.0])
        shift    (- drop-to (nth lo 1))
        vis      (vec (map-indexed (fn [i v] (if (= 1 (mod i 3)) (+ v shift) v)) vis))
        body     (xpbd/soft-body cage {:edge-compliance   (.-edgeCompliance controls)
                                       :volume-compliance (.-volumeCompliance controls)})
        skin     (skinning/bind cage vis)
        surface  (THREE/Mesh. geometry
                              (THREE/MeshPhongMaterial.
                               #js {:color 0xe0a85c :flatShading false
                                    :side THREE/DoubleSide}))
        cage-lines (THREE/LineSegments.
                    (doto (THREE/BufferGeometry.)
                      (.setAttribute "position"
                                     (THREE/BufferAttribute.
                                      (js/Float32Array. (* 3 (quot (count (:verts cage)) 3))) 3))
                      (.setIndex (clj->js (:edge-ids cage))))
                    (THREE/LineBasicMaterial. #js {:color 0x4a7ab0}))]
    (.add scene surface)
    (.add scene cage-lines)
    (swap! state assoc
           :geometry geometry :surface surface :cage-lines cage-lines
           :cage cage :body body :skin skin
           :vis-count (quot (count vis) 3)
           :tet-count (quot (count (:tet-ids cage)) 4))))

(defn- refresh! [{:keys [^js geometry ^js surface ^js cage-lines body cage skin]}]
  (set! (.-visible surface) (.-showSurface controls))
  (set! (.-visible cage-lines) (.-showCage controls))
  (when (.-showSurface controls)
    ;; The whole per-frame cost of the surface: one weighted sum of four
    ;; cage corners per vertex.
    (skinning/skin! skin (:tet-ids cage) (:pos body)
                    (.. geometry -attributes -position -array))
    (.computeVertexNormals geometry)
    (set! (.. geometry -attributes -position -needsUpdate) true)
    (.computeBoundingSphere geometry))
  (when (.-showCage controls)
    (let [^js arr (.. cage-lines -geometry -attributes -position -array)
          pos     (:pos body)]
      (.set arr pos)
      (set! (.. cage-lines -geometry -attributes -position -needsUpdate) true))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.7))
    (let [sun (THREE/DirectionalLight. 0xffffff 1.1)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 20 20)
                             (THREE/MeshPhongMaterial. #js {:color 0x141a2a}))]
      (set! (.-x (.-rotation floor)) (- (/ js/Math.PI 2)))
      (.add scene floor))
    (.add scene (doto (THREE/GridHelper. 6 12 0x2b3450 0x1b2236)
                  (-> .-position (.setY 0.002))))
    (.set (.-position camera) 2.0 1.5 2.8)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0.55 0)
      (letfn [(rebuild! [] (build! scene state) (refresh! @state))
              (squash! []
                (let [body (:body @state)]
                  (dotimes [i (:n body)]
                    (let [[px _ pz] (xpbd/particle body i)]
                      (xpbd/set-particle! body i [px 0.01 pz])))))
              (resize! []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (.setSize renderer w h)
                    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera))))
              (caption []
                (let [{:keys [vis-count tet-count]} @state]
                  (str tet-count " tets carrying " vis-count " vertices — "
                       (.toFixed (/ vis-count (max 1 tet-count)) 0) ":1")))
              (tick []
                (when @running?
                  (js/requestAnimationFrame tick)
                  (let [t0 (js/performance.now)]
                    (xpbd/step! (:body @state)
                                {:gravity  [0.0 (- (.-gravity controls)) 0.0]
                                 :substeps (.-substeps controls)
                                 :floor    0.0})
                    (refresh! @state)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. resize!) container)
        (resize!)
        (rebuild!)
        (let [gui (GUI. #js {:container container})
              readout (.add gui #js {:ratio (caption)} "ratio")]
          (.disable readout)
          (letfn [(regenerate! [] (rebuild!) (.setValue readout (caption)))]
            (doto gui
              (-> (.add controls "shape" #js ["knot" "torus" "sphere"]) (.onChange regenerate!))
              (-> (.add controls "detail" 1 8 1) (.onFinishChange regenerate!))
              (-> (.add controls "cage" 2 6 1) (.onFinishChange regenerate!))
              (.add controls "substeps" 1 25 1)
              (-> (.add controls "edgeCompliance" 0 0.002 0.0001) (.onChange regenerate!))
              (-> (.add controls "volumeCompliance" 0 0.002 0.0001) (.onChange regenerate!))
              (.add controls "gravity" 0 30 1)
              (.add controls "showCage")
              (.add controls "showSurface")
              (.add #js {:squash squash!} "squash")
              (.add #js {:drop regenerate!} "drop"))))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "skinning")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
