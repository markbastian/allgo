(ns allgo.demo.terrain-webgl
  (:require [allgo.demo.fps :as fps]
            [allgo.procedural.mesh :as mesh]
            [allgo.procedural.terrain :as terrain]
            [allgo.procedural.tin :as tin]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def grid-sizes
  "The sizes offered, per generator.

  Diamond-square only produces squares of `2^n + 1` on a side, since every
  round doubles the resolution and keeps the shared edge. The triangulated
  network can be sampled at any size and is given the same ladder, cut
  short: it resamples by asking which triangle covers each cell and it
  asks by looking through all of them, so its cost grows with the grid
  *and* with the mesh. Measured in this demo, a 1025 grid takes eighty-five
  seconds that way against five for diamond-square -- long enough to look
  like the page has hung, so it is not offered. Giving that resampling a
  spatial index would lift the limit."
  {:diamond-square [33 65 129 257 513 1025]
   :tin            [33 65 129 257]})

(defn sizes-for [generator]
  (get grid-sizes generator (:diamond-square grid-sizes)))

(def config
  {:generator    :diamond-square
   :size         257
   :width        1.0
   :tin-points   500
   :height-scale 200
   ;; How wide the terrain is in world units, held fixed as the grid size
   ;; changes. Resolution should decide how much detail a landscape has,
   ;; not how large it is -- and a landscape that grew with its grid would
   ;; leave the camera in the wrong place every time the size changed.
   :world-span   2048
   :noise-scale  180.0})

(defn size->iterations
  "How many rounds of diamond-square give a grid this size."
  [size]
  (long (js/Math.round (/ (js/Math.log (dec size)) (js/Math.log 2)))))

(defn- set-vec3! [^js arr idx x y z]
  (let [base (* idx 3)]
    (aset arr base x)
    (aset arr (+ base 1) y)
    (aset arr (+ base 2) z)))

(defn- diamond-square-grid [{:keys [size width]}]
  (let [grid-data (terrain/generate {:width width :iterations (size->iterations size)
                                     :corners [0.0 (rand) (rand) (rand)]})]
    {:grid (terrain/cells->grid grid-data) :dim (:dim grid-data)}))

(defn- tin-grid [{:keys [tin-points size smooth-passes] :or {tin-points 500 smooth-passes 3}}]
  {:grid (-> (tin/generate {:n tin-points :size 1.0})
             (tin/sample-grid size)
             (tin/smooth-grid smooth-passes))
   :dim  size})

(defn build-terrain-data
  "Builds a heightmap via `:generator` (`:diamond-square`, the default
  regular-grid midpoint-displacement fractal, or `:tin`, a Delaunay-
  triangulated irregular network per `allgo.procedural.tin`)."
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
  smoothly, matching `allgo.desktop.terrain-shape/build-mesh-shape` but using
  the standard y-up world (+height, not the Quil renderer's -height)."
  [{:keys [grid dim lo hi]} {:keys [cell-scale height-scale noise-scale]}]
  (let [n          (* dim dim)
        offset     (* 0.5 (dec dim) cell-scale)
        positions  (js/Float32Array. (* n 3))
        normals    (js/Float32Array. (* n 3))
        colors     (js/Float32Array. (* n 3))
        ;; Once each, not twice: the slope was being computed to find the
        ;; steepest and then again for every vertex, which at a thousand
        ;; cells a side is two million gradient evaluations for one
        ;; million vertices.
        slopes     (let [a (js/Float32Array. n)]
                     (dotimes [i dim]
                       (dotimes [j dim]
                         (aset a (+ (* i dim) j)
                               (mesh/slope-at grid dim i j cell-scale height-scale))))
                     a)
        max-slope  (let [m (areduce slopes k acc 1e-6 (max acc (aget slopes k)))] m)]
    (doseq [i (range dim) j (range dim)]
      (let [idx        (+ (* i dim) j)
            h          (get-in grid [i j])
            height-t   (if (= lo hi) 0.5 (/ (double (- h lo)) (- hi lo)))
            slope-t    (/ (aget slopes idx) max-slope)
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
        ;; The world stays the same size; the cells get smaller.
        config       (assoc config :cell-scale
                            (/ (:world-span config)
                               (max 1 (dec (long (:dim terrain-data))))))
        geometry     (build-geometry terrain-data config)
        material     (THREE/MeshStandardMaterial. #js {:vertexColors true
                                                       :side         THREE/DoubleSide
                                                       :roughness    0.9
                                                       :metalness    0.0})]
    {:span (:world-span config)
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
        state    (atom {:terrain-mesh mesh :wireframe? false :generator (:generator config)
                        :size (:size config) :hovering? false})]
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
                (let [{:keys [mesh]} (build-mesh (assoc config
                                                        :generator (:generator @state)
                                                        :size (:size @state)))]
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
              gui-state #js {:generator (name (:generator @state))
                             :size (:size @state)
                             :wireframe false :regenerate regenerate!}]
          ;; The size control is rebuilt when the generator changes, because
          ;; the two do not offer the same sizes.
          (let [size-ctl (atom nil)
                generator-ctl (atom nil)
                on-size (fn [v] (swap! state assoc :size (long v)) (regenerate!))
                make-size! (fn [generator]
                             (let [allowed (sizes-for generator)
                                   chosen (if (some #{(:size @state)} allowed)
                                            (:size @state)
                                            (last allowed))]
                               (swap! state assoc :size chosen)
                               (set! (.-size gui-state) chosen)
                               (let [^js c (-> (.add gui gui-state "size" (clj->js allowed))
                                               (.onChange on-size))
                                     ^js after (.-domElement ^js @generator-ctl)]
                                 ;; A rebuilt controller is appended at the
                                 ;; end; put it back beside the generator it
                                 ;; belongs to.
                                 (.insertBefore (.-parentNode after)
                                                (.-domElement c)
                                                (.-nextSibling after))
                                 (reset! size-ctl c))))]
            (reset! generator-ctl
                    (-> (.add gui gui-state "generator" #js ["tin" "diamond-square"])
                        (.onChange (fn [v]
                                     (let [generator (keyword v)
                                           allowed (sizes-for generator)
                                           ;; Keep the size if the new
                                           ;; generator offers it, and come
                                           ;; down to its largest if not.
                                           chosen (if (some #{(:size @state)} allowed)
                                                    (:size @state)
                                                    (last allowed))]
                                       (swap! state assoc :generator generator :size chosen)
                                       (when-let [^js c @size-ctl] (.destroy c))
                                       (make-size! generator)
                                       (regenerate!))))))
            (make-size! (:generator @state)))
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
