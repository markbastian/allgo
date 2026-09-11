(ns allgo.demo.water
  "WebGL demo for `allgo.physics.height-field`, after Ten Minute Physics 20.

  A tank of water with things floating in it. The water is one height per
  column, not a volume, which is why a whole pond runs at full rate and
  also why nothing here can ever break or overturn -- a height field has
  no way to represent water above water.

  `density` is the parameter to play with, and it is worth watching
  against the waterline: a body settles with exactly that fraction of
  itself under, because the water it pushes aside weighs what the body
  weighs. That is Archimedes, and nothing in the code puts it there -- it
  falls out of summing the displaced depth column by column.

  Click the water to drop a splash into it. `waveSpeed` is capped by what
  the grid can carry, so winding it up stops making waves faster once the
  cap is reached rather than blowing the surface apart; `velDamping` is
  how fast the pond goes quiet, and `posDamping` smooths the grid-scale
  chop the stencil cannot resolve.

  `drag` decides only how long a body takes to settle, not where it ends
  up. A body much lighter than water bobs for a long time at the default,
  because once it leaves the water there is nothing left to damp it."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.height-field :as hf]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:spacing     0.04
       :waveSpeed   2.0
       :posDamping  1.0
       :velDamping  0.3
       :alpha       0.5
       :drag        12.0
       :bodies      3
       :density     0.7
       :radius      0.22
       :showWireframe false})

(def ^:private tank-size [2.5 1.0 3.0])
(def ^:private tank-border 0.03)
(def ^:private water-depth 0.8)

(defn- make-bodies
  "`n` balls spread across the tank, alternating either side of the
  chosen density so the difference is visible side by side."
  [n]
  (let [d (.-density controls)
        r (.-radius controls)]
    (vec (for [i (range n)]
           (let [t (if (= 1 n) 0.5 (/ i (dec n)))
                 ;; Spread the densities around the slider, staying inside
                 ;; a range that still floats.
                 dens (max 0.1 (min 3.0 (* d (+ 0.5 (* 1.2 t)))))]
             (hf/sphere {:pos [(+ -0.7 (* 1.4 t)) 1.1 (* 0.6 (- t 0.5))]
                         :radius r :density dens}))))))

(defn- build! [^js scene state]
  (doseq [k [:water-obj :wire]]
    (when-let [^js old (k @state)]
      (.remove scene old)
      (.dispose (.-geometry old))
      (.dispose (.-material old))))
  (doseq [^js old (:ball-objs @state)]
    (.remove scene old)
    (.dispose (.-geometry old))
    (.dispose (.-material old)))
  (let [[sx _ sz] tank-size
        surface (hf/surface {:size-x sx :size-z sz
                             :spacing (.-spacing controls)
                             :depth water-depth})
        {:keys [nx nz]} surface
        verts (js/Float32Array. (* 3 nx nz))
        uvs   (js/Float32Array. (* 2 nx nz))
        idx   (js/Uint32Array. (* 6 (dec nx) (dec nz)))]
    (dotimes [i nx]
      (dotimes [j nz]
        (let [id (+ (* i nz) j)]
          (aset verts (* 3 id) (hf/column-x surface i))
          (aset verts (+ (* 3 id) 1) water-depth)
          (aset verts (+ (* 3 id) 2) (hf/column-z surface j))
          (aset uvs (* 2 id) (/ i nx))
          (aset uvs (+ (* 2 id) 1) (/ j nz)))))
    (let [p (atom 0)]
      (dotimes [i (dec nx)]
        (dotimes [j (dec nz)]
          (let [id0 (+ (* i nz) j)
                id1 (+ (* i nz) j 1)
                id2 (+ (* (inc i) nz) j 1)
                id3 (+ (* (inc i) nz) j)]
            (doseq [v [id0 id1 id2 id0 id2 id3]]
              (aset idx @p v)
              (swap! p inc))))))
    (let [geom (doto (THREE/BufferGeometry.)
                 (.setAttribute "position" (THREE/BufferAttribute. verts 3))
                 (.setAttribute "uv" (THREE/BufferAttribute. uvs 2))
                 (.setIndex (THREE/BufferAttribute. idx 1)))
          mat  (THREE/MeshPhongMaterial.
                #js {:color 0x3a7ebf :transparent true :opacity 0.85
                     :shininess 90 :side THREE/DoubleSide})
          obj  (THREE/Mesh. geom mat)
          wire (THREE/Mesh. geom (THREE/MeshBasicMaterial.
                                  #js {:color 0x8fd0ff :wireframe true}))
          bodies (make-bodies (.-bodies controls))
          ball-objs
          (mapv (fn [b]
                  (let [m (THREE/Mesh.
                           (THREE/SphereGeometry. (:radius b) 24 18)
                           (THREE/MeshPhongMaterial.
                            #js {:color (if (< (/ (:mass b) (hf/volume b)) 1.0)
                                          0xf2b134 0xd94f3d)}))]
                    (.add scene m)
                    m))
                bodies)]
      (.add scene obj)
      (.add scene wire)
      (swap! state assoc
             :surface surface :bodies bodies :verts verts
             :geom geom :water-obj obj :wire wire :ball-objs ball-objs))))

(defn- refresh! [{:keys [^js geom ^js wire surface bodies ^js verts ball-objs]}]
  (let [{:keys [nx nz]} surface]
    (dotimes [i nx]
      (dotimes [j nz]
        (let [id (+ (* i nz) j)]
          (aset verts (+ (* 3 id) 1) (hf/height surface i j))))))
  (set! (.-visible wire) (.-showWireframe controls))
  (.computeVertexNormals geom)
  (set! (.. geom -attributes -position -needsUpdate) true)
  (.computeBoundingSphere geom)
  (doseq [[^js m b] (map vector ball-objs bodies)]
    (let [[x y z] (:pos b)]
      (.set (.-position m) x y z))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 200)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        raycaster (THREE/Raycaster.)
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.7))
    (let [sun (THREE/DirectionalLight. 0xffffff 1.0)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    ;; The tank: four walls and a floor, drawn only so the water has
    ;; something to sit in.
    (let [[sx sy sz] tank-size
          b tank-border
          mat (THREE/MeshPhongMaterial. #js {:color 0x6b7280
                                             :transparent true :opacity 0.35
                                             :side THREE/DoubleSide})
          add! (fn [gx gy gz px py pz]
                 (let [m (THREE/Mesh. (THREE/BoxGeometry. gx gy gz) mat)]
                   (.set (.-position m) px py pz)
                   (.add scene m)))]
      (add! b sy sz (* -0.5 sx) (* 0.5 sy) 0)
      (add! b sy sz (* 0.5 sx) (* 0.5 sy) 0)
      (add! sx sy b 0 (* 0.5 sy) (* -0.5 sz))
      (add! sx sy b 0 (* 0.5 sy) (* 0.5 sz))
      (let [floor (THREE/Mesh. (THREE/BoxGeometry. sx 0.02 sz)
                               (THREE/MeshPhongMaterial. #js {:color 0x1b2236}))]
        (.set (.-position floor) 0 -0.01 0)
        (.add scene floor)))
    (.set (.-position camera) 0 2.4 3.6)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0.4 0)
      (letfn [(rebuild! []
                (build! scene state)
                (refresh! @state))
              (resize! []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (.setSize renderer w h)
                    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera))))
              (splash-at! [^js e]
                ;; Raycast onto the resting water plane; close enough to
                ;; the real surface to aim at, and it cannot miss.
                (let [rect (.getBoundingClientRect (.-domElement renderer))
                      nx (- (* 2 (/ (- (.-clientX e) (.-left rect)) (.-width rect))) 1)
                      ny (- 1 (* 2 (/ (- (.-clientY e) (.-top rect)) (.-height rect))))
                      plane (THREE/Plane. (THREE/Vector3. 0 1 0) (- water-depth))
                      hit (THREE/Vector3.)]
                  (.setFromCamera raycaster (THREE/Vector2. nx ny) camera)
                  (when (.intersectPlane (.-ray raycaster) plane hit)
                    (hf/splash! (:surface @state) (.-x hit) (.-z hit) 0.18 0.25))))
              (tick []
                (when @running?
                  (js/requestAnimationFrame tick)
                  (let [t0 (js/performance.now)
                        {:keys [surface bodies]} @state
                        w (hf/step! {:surface surface :bodies bodies}
                                    {:tank {:size tank-size :border tank-border}
                                     :wave-speed  (.-waveSpeed controls)
                                     :pos-damping (.-posDamping controls)
                                     :vel-damping (.-velDamping controls)
                                     :alpha       (.-alpha controls)
                                     :drag        (.-drag controls)})]
                    (swap! state assoc :bodies (:bodies w))
                    (refresh! @state)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.addEventListener (.-domElement renderer) "dblclick" splash-at!)
        (.observe (js/ResizeObserver. resize!) container)
        (resize!)
        (rebuild!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "spacing" 0.02 0.08 0.005) (.onFinishChange rebuild!))
            (.add controls "waveSpeed" 0.5 8 0.1)
            (.add controls "posDamping" 0 4 0.1)
            (.add controls "velDamping" 0 3 0.05)
            (.add controls "alpha" 0 1 0.05)
            (.add controls "drag" 0 40 1)
            (-> (.add controls "bodies" 0 6 1) (.onChange rebuild!))
            (-> (.add controls "density" 0.1 2.0 0.05) (.onChange rebuild!))
            (-> (.add controls "radius" 0.1 0.35 0.01) (.onChange rebuild!))
            (.add controls "showWireframe")
            (.add #js {:splash (fn []
                                 (hf/splash! (:surface @state) 0.0 0.0 0.25 0.35))}
                  "splash")
            (.add #js {:reset rebuild!} "reset")))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "water")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
