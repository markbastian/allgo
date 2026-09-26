(ns allgo.demo.spatial-hash
  "WebGL demo for `allgo.spatial.hash`, after Ten Minute Physics 11.

  Thousands of balls bouncing in a box, every pair of overlaps found and
  resolved each frame. The `broadPhase` control switches between the hash
  and the obvious double loop over all pairs, which is the whole point:
  the simulation is identical either way, and only one of them survives a
  few thousand balls.

  Balls are drawn with one instanced mesh -- a single draw call for the
  lot, since issuing thousands of them would cost more than the physics
  and would misattribute the difference the broad phase makes."
  (:require [allgo.demo.fps :as fps]
            [allgo.spatial.hash :as hash]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:balls        1500
       :radius       0.025
       :gravity      3.0
       :broadPhase   "spatial hash"
       :showContacts true})

(def ^:private bounds 1.0)

;; ---------------------------------------------------------------------------
;; The simulation, on flat arrays

(defn- spawn
  "`n` balls on a loose grid so none start overlapping, with a random
  drift."
  [n radius]
  (let [pos     (js/Float64Array. (* 3 n))
        vel     (js/Float64Array. (* 3 n))
        spacing (* 3.0 radius)
        per-row (js/Math.floor (/ (* 2 bounds) spacing))]
    (dotimes [i n]
      (let [b  (* 3 i)
            ;; Fill the floor of the box first and stack upward, so a
            ;; large count starts as a slab across the base rather than a
            ;; wall down one side.
            ix (mod i per-row)
            iz (mod (js/Math.floor (/ i per-row)) per-row)
            iy (js/Math.floor (/ i (* per-row per-row)))]
        (aset pos b (+ (- bounds) radius (* spacing ix)))
        (aset pos (+ b 1) (+ radius (* spacing iy)))
        (aset pos (+ b 2) (+ (- bounds) radius (* spacing iz)))
        (aset vel b (- (rand 2.0) 1.0))
        (aset vel (+ b 1) (- (rand 2.0) 1.0))
        (aset vel (+ b 2) (- (rand 2.0) 1.0))))
    {:pos pos :vel vel :n n}))

(defn- integrate!
  "Move under gravity and bounce off the walls of the box."
  [^js pos ^js vel n dt gravity radius]
  (dotimes [i n]
    (let [b (* 3 i)]
      (aset vel (+ b 1) (- (aget vel (+ b 1)) (* gravity dt)))
      (dotimes [d 3]
        (let [k (+ b d)]
          (aset pos k (+ (aget pos k) (* (aget vel k) dt)))
          (let [lo (if (= d 1) radius (+ (- bounds) radius))
                hi (- bounds radius)]
            (cond
              (< (aget pos k) lo) (do (aset pos k lo) (aset vel k (- (aget vel k))))
              (> (aget pos k) hi) (do (aset pos k hi) (aset vel k (- (aget vel k)))))))))))

(defn- resolve-pair!
  "Separate two overlapping balls and exchange the velocity along the line
  between them -- an equal-mass elastic collision, which is just a swap of
  the normal components."
  [^js pos ^js vel i j min-dist]
  (let [a  (* 3 i) b (* 3 j)
        dx (- (aget pos a) (aget pos b))
        dy (- (aget pos (+ a 1)) (aget pos (+ b 1)))
        dz (- (aget pos (+ a 2)) (aget pos (+ b 2)))
        d2 (+ (* dx dx) (* dy dy) (* dz dz))]
    (when (and (pos? d2) (< d2 (* min-dist min-dist)))
      (let [d  (js/Math.sqrt d2)
            nx (/ dx d) ny (/ dy d) nz (/ dz d)
            corr (* 0.5 (- min-dist d))]
        (aset pos a (+ (aget pos a) (* nx corr)))
        (aset pos (+ a 1) (+ (aget pos (+ a 1)) (* ny corr)))
        (aset pos (+ a 2) (+ (aget pos (+ a 2)) (* nz corr)))
        (aset pos b (- (aget pos b) (* nx corr)))
        (aset pos (+ b 1) (- (aget pos (+ b 1)) (* ny corr)))
        (aset pos (+ b 2) (- (aget pos (+ b 2)) (* nz corr)))
        (let [vi (+ (* (aget vel a) nx) (* (aget vel (+ a 1)) ny) (* (aget vel (+ a 2)) nz))
              vj (+ (* (aget vel b) nx) (* (aget vel (+ b 1)) ny) (* (aget vel (+ b 2)) nz))
              dv (- vj vi)]
          (aset vel a (+ (aget vel a) (* nx dv)))
          (aset vel (+ a 1) (+ (aget vel (+ a 1)) (* ny dv)))
          (aset vel (+ a 2) (+ (aget vel (+ a 2)) (* nz dv)))
          (aset vel b (- (aget vel b) (* nx dv)))
          (aset vel (+ b 1) (- (aget vel (+ b 1)) (* ny dv)))
          (aset vel (+ b 2) (- (aget vel (+ b 2)) (* nz dv))))
        true))))

(defn- collide-hashed!
  "Overlaps via the spatial hash. Returns `[contacts checks]`."
  [h ^js pos ^js vel n min-dist ^js touched]
  (hash/rebuild! h pos n)
  (let [contacts (js/Int32Array. 1)
        checks   (js/Int32Array. 1)]
    (dotimes [i n]
      (let [found (hash/query! h pos i min-dist)]
        (dotimes [k found]
          (let [j (hash/neighbor h k)]
            (when (< i j)
              (aset checks 0 (inc (aget checks 0)))
              (when (resolve-pair! pos vel i j min-dist)
                (aset touched i 1)
                (aset touched j 1)
                (aset contacts 0 (inc (aget contacts 0)))))))))
    [(aget contacts 0) (aget checks 0)]))

(defn- collide-brute!
  "The same, comparing every pair. Correct, and quadratic."
  [^js pos ^js vel n min-dist ^js touched]
  (let [contacts (js/Int32Array. 1)
        checks   (js/Int32Array. 1)]
    (dotimes [i n]
      (loop [j (inc i)]
        (when (< j n)
          (aset checks 0 (inc (aget checks 0)))
          (when (resolve-pair! pos vel i j min-dist)
            (aset touched i 1)
            (aset touched j 1)
            (aset contacts 0 (inc (aget contacts 0))))
          (recur (inc j)))))
    [(aget contacts 0) (aget checks 0)]))

;; ---------------------------------------------------------------------------

(defn- build-scene [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 1 0.01 100)
        renderer (THREE/WebGLRenderer. #js {:antialias true})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.add scene (THREE/AmbientLight. 0x8899bb 0.75))
    (let [sun (THREE/DirectionalLight. 0xffffff 1.0)]
      (.set (.-position sun) 2 4 3)
      (.add scene sun))
    (let [box (THREE/Box3. (THREE/Vector3. (- bounds) 0 (- bounds))
                           (THREE/Vector3. bounds (* 2 bounds) bounds))]
      (.add scene (THREE/Box3Helper. box (THREE/Color. 0x2b3450))))
    (.set (.-position camera) 2.4 2.0 3.4)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 bounds 0)
      (.update orbit)
      {:scene scene :camera camera :renderer renderer :orbit orbit})))

(defn init! [^js container]
  (let [{:keys [^js scene ^js camera ^js renderer ^js orbit]} (build-scene container)
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom {})
        scratch   (THREE/Matrix4.)
        plain     (THREE/Color. 0x6f9fe0)
        hit       (THREE/Color. 0xe0a85c)]
    (letfn [(rebuild! []
              (when-let [^js old (:instanced @state)]
                (.remove scene old)
                (.dispose (.-geometry old))
                (.dispose (.-material old)))
              (let [n      (.-balls controls)
                    radius (.-radius controls)
                    {:keys [pos vel]} (spawn n radius)
                    geom   (THREE/SphereGeometry. radius 6 5)
                    mat    (THREE/MeshPhongMaterial. #js {:color 0xffffff})
                    inst   (THREE/InstancedMesh. geom mat n)]
                (.add scene inst)
                (dotimes [i n] (.setColorAt inst i plain))
                (swap! state assoc
                       :pos pos :vel vel :n n :radius radius
                       :instanced inst
                       :touched (js/Int32Array. n)
                       ;; Cells one diameter across: the query radius, so a
                       ;; ball's neighbors are never more than a cell away.
                       :hash (hash/spatial-hash (* 2 radius) n))))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (.setSize renderer w h)
                  (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
                  (set! (.-aspect camera) (/ w h))
                  (.updateProjectionMatrix camera))))
            (draw! []
              (let [{:keys [^js instanced ^js pos n ^js touched]} @state]
                (dotimes [i n]
                  (let [b (* 3 i)]
                    (.setPosition scratch (aget pos b) (aget pos (+ b 1)) (aget pos (+ b 2)))
                    (.setMatrixAt instanced i scratch)
                    (when (.-showContacts controls)
                      (.setColorAt instanced i (if (pos? (aget touched i)) hit plain)))))
                (set! (.. instanced -instanceMatrix -needsUpdate) true)
                (when (and (.-showContacts controls) (.-instanceColor instanced))
                  (set! (.. instanced -instanceColor -needsUpdate) true))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      {:keys [pos vel n radius hash ^js touched]} @state
                      min-dist (* 2 radius)]
                  (.fill touched 0)
                  (integrate! pos vel n (/ 1.0 60.0) (.-gravity controls) radius)
                  (let [[contacts checks]
                        (if (= "brute force" (.-broadPhase controls))
                          (collide-brute! pos vel n min-dist touched)
                          (collide-hashed! hash pos vel n min-dist touched))]
                    (swap! state assoc :contacts contacts :checks checks))
                  (draw!)
                  (.update orbit)
                  (.render renderer scene camera)
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (rebuild!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "balls" 100 6000 100) (.onFinishChange rebuild!))
          (-> (.add controls "radius" 0.01 0.06 0.005) (.onFinishChange rebuild!))
          (.add controls "gravity" 0 20 0.5)
          (.add controls "broadPhase" #js ["spatial hash" "brute force"])
          (.add controls "showContacts")
          (.add #js {:restart rebuild!} "restart")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "spatial-hash")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
