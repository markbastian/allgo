(ns allgo.demo.geometry
  "Three.js demo for `allgo.procedural.instancing` and
  `allgo.procedural.lsystem`: trees, a forest, two fractals and four
  L-system plants, each grown from a few rules.

  Both namespaces produce the same thing -- a list of primitives, each a
  shape, a transform and a color -- and every shape is one
  `InstancedMesh`, so a forest of a hundred thousand limbs is a handful
  of draw calls.

  For the instanced models, `detail` is the smallest bounding sphere, in
  pixels, that is worth expanding: anything smaller is drawn as its
  proxy. With `follow camera` on, that is measured from wherever the
  camera is, and the model is expanded again whenever the camera stops
  moving -- walk into the forest and the trees around you sprout twigs
  and leaves while the far ones stay crowns. `depth` caps the recursion
  outright. For the plants, `depth` is how many rounds the string is
  rewritten. The caption counts what was generated and how long it took."
  (:require [allgo.demo.fps :as fps]
            [allgo.procedural.instancing :as in]
            [allgo.procedural.lsystem :as ls]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]
            [clojure.math :as math]))

(def ^:private models
  (into (mapv #(assoc % :kind :instanced) in/catalog)
        (map #(assoc % :kind :plant) ls/plants)))

(def ^:private ^js controls
  #js {:model        "tree"
       :seed         1
       :detail       3
       :depth        12
       :followCamera true})

(def ^:private fov 50)

;; A cap on what one expansion may generate: it runs on the page's own
;; thread, and past this a fine detail in the forest takes seconds.
(def ^:private limit 120000)

(defn- model [] (first (filter #(= (.-model controls) (:name %)) models)))

(defn- default-depth
  "A plant's own number of rounds; for an instanced model, deep enough
  that detail is what stops it."
  [{:keys [kind rounds]}]
  (if (= kind :plant) rounds 12))

(defn- geometries
  "One unit geometry per shape, in the frame its transform expects."
  []
  {:cylinder    (doto (THREE/CylinderGeometry. 0.85 1 1 7 1 true) (.translate 0 0.5 0))
   :leaf        (let [s (doto (THREE/Shape.)
                          (.moveTo 0 0)
                          (.quadraticCurveTo 0.32 0.45 0 1)
                          (.quadraticCurveTo -0.32 0.45 0 0))]
                  (THREE/ShapeGeometry. s 3))
   :sphere      (THREE/IcosahedronGeometry. 1 1)
   :tetrahedron (THREE/TetrahedronGeometry. 1)
   :box         (THREE/BoxGeometry. 2 2 2)})

(defn- materials []
  {:leaf   (THREE/MeshLambertMaterial. #js {:side THREE/DoubleSide})
   :sphere (THREE/MeshLambertMaterial. #js {:flatShading true})
   :other  (THREE/MeshLambertMaterial.)})

(defn- write-mesh!
  "Copies `prims` into `mesh`'s instance buffers: our transforms are
  row-major 3x4, three's are column-major 4x4."
  [^js mesh prims]
  (let [^js m (.. mesh -instanceMatrix -array)
        ^js c (.. mesh -instanceColor -array)]
    (loop [i 0 prims (seq prims)]
      (when prims
        (let [{:keys [transform color]} (first prims)
              [a b cc x d e f y g h k z] transform
              o (* 16 i)
              [r gr bl] (or color [0.8 0.8 0.8])
              q (* 3 i)]
          (aset m o a) (aset m (+ o 1) d) (aset m (+ o 2) g) (aset m (+ o 3) 0)
          (aset m (+ o 4) b) (aset m (+ o 5) e) (aset m (+ o 6) h) (aset m (+ o 7) 0)
          (aset m (+ o 8) cc) (aset m (+ o 9) f) (aset m (+ o 10) k) (aset m (+ o 11) 0)
          (aset m (+ o 12) x) (aset m (+ o 13) y) (aset m (+ o 14) z) (aset m (+ o 15) 1)
          (aset c q r) (aset c (+ q 1) gr) (aset c (+ q 2) bl)
          (recur (inc i) (next prims)))))
    (set! (.-count mesh) (count prims))
    (set! (.. mesh -instanceMatrix -needsUpdate) true)
    (set! (.. mesh -instanceColor -needsUpdate) true)))

(defn- caption-text [{:keys [kind primitives instances proxies string-length ms]}]
  (str (count primitives) " primitives"
       (if (= kind :plant)
         (str " · string of " string-length " symbols")
         (str " · " instances " instances, " proxies " of them drawn as proxies"))
       " · " (.toFixed ms 0) " ms"))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. fov 1 0.05 500)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        geos     (geometries)
        mats     (materials)
        meshes   (atom {})
        ground   (doto (THREE/Mesh. (THREE/CircleGeometry. 40 48)
                                    (THREE/MeshLambertMaterial. #js {:color 0x4a5a34}))
                   (.rotateX (- (/ math/PI 2))))
        caption  (js/document.createElement "div")
        running? (atom false)
        dirty?   (atom true)
        tick-fps! (fps/meter! container)
        settle   (atom nil)]
    (set! (.-background scene) (THREE/Color. 0x9fb8cf))
    (set! (.-fog scene) (THREE/Fog. 0x9fb8cf 25 70))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (set! (.. renderer -domElement -style -borderRadius) "8px")
    (.appendChild container (.-domElement renderer))
    (set! (.-style caption)
          "position:absolute;left:10px;bottom:8px;font:11px ui-monospace,Menlo,monospace;color:#1d2735;pointer-events:none")
    (.appendChild container caption)
    (.add scene ground)
    (.add scene (THREE/HemisphereLight. 0xdfeaf5 0x3a3020 0.75))
    (let [sun (THREE/DirectionalLight. 0xfff2dc 0.8)]
      (.set (.-position sun) 30 50 20)
      (.add scene sun))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (letfn [(mesh-for [shape n]
                ;; Grown, never shrunk: a buffer big enough for the largest
                ;; expansion so far, drawn to `count`.
                (let [^js have (get @meshes shape)]
                  (if (and have (>= (.. have -instanceMatrix -count) n))
                    have
                    (let [cap (max 64 (* 2 n))
                          mat (get mats shape (:other mats))
                          ^js mesh (THREE/InstancedMesh. (get geos shape) mat cap)]
                      (.setColorAt mesh 0 (THREE/Color.))
                      ;; Three 0.150 culls an instanced mesh by its one
                      ;; geometry, not by where the instances are.
                      (set! (.-frustumCulled mesh) false)
                      (when have (.remove scene have) (.dispose have))
                      (.add scene mesh)
                      (swap! meshes assoc shape mesh)
                      mesh))))
              (pixel-angle []
                ;; The angle one pixel subtends at the middle of the view.
                (/ (math/to-radians fov) (max 1 (.-clientHeight container))))
              (generate []
                (let [{:keys [kind make] :as m} (model)
                      seed (.-seed controls)]
                  (if (= kind :plant)
                    (assoc (ls/grow m {:seed seed :rounds (.-depth controls)}) :kind kind)
                    (let [p (.-position camera)]
                      (assoc (in/expand (make {:seed seed})
                                        (cond-> {:max-depth (.-depth controls)
                                                 :limit limit}
                                          (.-followCamera controls)
                                          (assoc :eye [(.-x p) (.-y p) (.-z p)]
                                                 :min-size (* (.-detail controls) (pixel-angle)))
                                          (not (.-followCamera controls))
                                          (assoc :min-size (* 0.01 (.-detail controls)))))
                             :kind kind)))))
              (rebuild! []
                (let [t0 (js/performance.now)
                      result (generate)
                      by-shape (group-by :shape (:primitives result))]
                  (doseq [[shape ^js mesh] @meshes]
                    (when-not (contains? by-shape shape) (set! (.-count mesh) 0)))
                  (doseq [[shape prims] by-shape]
                    (write-mesh! (mesh-for shape (count prims)) prims))
                  (set! (.-visible ground) (or (= :plant (:kind result))
                                               (contains? #{"tree" "forest"} (.-model controls))))
                  (set! (.-textContent caption)
                        (caption-text (assoc result :ms (- (js/performance.now) t0))))
                  (reset! dirty? true)))
              (frame! []
                (let [{:keys [distance target]} (model)
                      [tx ty tz] target]
                  (.set (.-target orbit) tx ty tz)
                  (.set (.-position camera) (+ tx (* 0.75 distance)) (+ ty (* 0.12 distance)) (+ tz (* 0.6 distance)))
                  (.update orbit)))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h)
                    (reset! dirty? true))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)
                        moved? (.update orbit)]
                    (when (or moved? @dirty?)
                      (reset! dirty? false)
                      (.render renderer scene camera)
                      (tick-fps! (- (js/performance.now) t0))))))]
        ;; Expand again once the camera has come to rest, not on every
        ;; frame of a drag.
        (.addEventListener orbit "change"
                           (fn []
                             (reset! dirty? true)
                             (when (and (.-followCamera controls) (= :instanced (:kind (model))))
                               (js/clearTimeout @settle)
                               (reset! settle (js/setTimeout rebuild! 250)))))
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})
              ^js detail (.add gui controls "detail" 0.5 20 0.5)
              ^js follow (.add gui controls "followCamera")]
          (letfn [(show-for-kind! []
                    (let [inst? (= :instanced (:kind (model)))]
                      (.show detail inst?)
                      (.show follow inst?)))]
            (-> (.add gui controls "model" (clj->js (mapv :name models)))
                (.onChange (fn [_]
                             (set! (.-depth controls) (default-depth (model)))
                             (show-for-kind!) (frame!) (rebuild!))))
            (-> (.add gui controls "seed" 1 99 1) (.onFinishChange rebuild!))
            (.onFinishChange detail rebuild!)
            (-> (.add gui controls "depth" 1 14 1) (.listen) (.onFinishChange rebuild!))
            (.onChange follow rebuild!)
            (show-for-kind!)))
        (on-resize)
        (frame!)
        (rebuild!)
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "geometry")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
