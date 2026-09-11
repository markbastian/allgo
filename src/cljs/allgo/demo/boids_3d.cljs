(ns allgo.demo.boids-3d
  "Three.js demo for `allgo.simulation.boids`, driven by the same rules as the 2D
  canvas viewer -- the flocking namespace is arity-generic, so going to three
  dimensions is purely a matter of passing 3-element bounds.

  Boids are instanced cones oriented along their heading and tinted by it
  (hue from compass bearing, brightness from climb), so a flock settling into
  alignment reads as the swarm converging on a single colour."
  (:require [allgo.demo.fps :as fps]
            [allgo.simulation.boids :as boids]
            [allgo.simulation.boids-flat :as flat]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private world [100 100 100])

;; Obstacle geometry declared once: it builds both the convex bodies GJK is
;; queried against and the meshes drawn for them.

(def ^:private shapes
  [{:kind :sphere :c [50 50 50] :r 17}
   {:kind :box    :lo [8 8 62] :hi [38 38 92]}])

(def ^:private obstacles
  (mapv (fn [{:keys [kind c r lo hi]}]
          (case kind
            :sphere (boids/sphere-obstacle c r)
            :box    (boids/box-obstacle lo hi)))
        shapes))
;; The instanced mesh is allocated once at this size and its draw count
;; set per frame, so this is the ceiling the slider can reach.
(def ^:private capacity 4000)

(def ^:private ^js controls
  #js {:boids         60
       :separation    1.6
       :alignment     1.0
       :cohesion      0.9
       :perception    22
       :personalSpace 9
       :speed         0.8
       :avoidance     3.0
       :obstacles     true
       ;; The same rules on flat arrays, obstacles included. Several times
       ;; faster and identical to within floating-point noise, which is
       ;; what `allgo.simulation.flock/divergence` checks.
       :flatArrays    true})

(defn- flat? [] (.-flatArrays controls))

(defn- params []
  {:separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (.-personalSpace controls)
   :max-speed         (.-speed controls)
   :max-force         0.02
   :edges             :bounce
   :boid-radius       1.4
   :avoid-radius      15.0
   :avoid-weight      (.-avoidance controls)
   :obstacles         (if (.-obstacles controls) obstacles [])})

(defn- boid-geometry
  "A cone nose-up rotated onto +Z, the axis `Object3D.lookAt` aligns."
  []
  (doto (THREE/ConeGeometry. 0.9 2.8 6)
    (.rotateX (/ js/Math.PI 2))))

(defn- obstacle-meshes
  "Meshes for the obstacles, shifted from boid space into the centred world."
  [[w h d]]
  (let [material (THREE/MeshLambertMaterial. #js {:color 0x39456b})
        offset   [(/ w 2.0) (/ h 2.0) (/ d 2.0)]]
    (for [{:keys [kind c r lo hi]} shapes
          :let [geo (case kind
                      :sphere (THREE/SphereGeometry. r 24 16)
                      :box    (let [[x0 y0 z0] lo [x1 y1 z1] hi]
                                (THREE/BoxGeometry. (- x1 x0) (- y1 y0) (- z1 z0))))
                [cx cy cz] (boids/v- (case kind
                                       :sphere c
                                       :box    (boids/v* (boids/v+ (vec lo) (vec hi)) 0.5))
                                     offset)]]
      (doto (THREE/Mesh. geo material)
        (-> .-position (.set cx cy cz))))))

(defn- bounds-box [[w h d]]
  (THREE/LineSegments.
   (THREE/EdgesGeometry. (THREE/BoxGeometry. w h d))
   (THREE/LineBasicMaterial. #js {:color 0x2c3a5c})))

(defn- resize-to-flock [flock n max-speed]
  (let [have (count flock)]
    (cond
      (= have n) flock
      (< n have) (subvec flock 0 n)
      :else      (into flock (repeatedly (- n have) #(boids/random-boid world max-speed))))))

(defn- write-instances!
  "Push flock state into the InstancedMesh: one matrix and colour per boid.
  Boid space is [0,size); the mesh is centred on the origin, hence the offset."
  [^js mesh ^js scratch ^js color flock [w h d]]
  (let [offset [(/ w 2) (/ h 2) (/ d 2)]]
    (doseq [[i {:keys [pos vel]}] (map-indexed vector flock)]
      (let [[x y z]    (boids/v- pos offset)
            [nx ny nz] (boids/normalize vel)]
        (.set (.-position scratch) x y z)
        (.lookAt scratch (+ x nx) (+ y ny) (+ z nz))
        (.updateMatrix scratch)
        (.setMatrixAt mesh i (.-matrix scratch))
        (.setHSL color
                 (+ 0.5 (/ (js/Math.atan2 nz nx) (* 2 js/Math.PI)))
                 0.8
                 (+ 0.55 (* 0.18 ny)))
        (.setColorAt mesh i color)))
    (set! (.-count mesh) (count flock))
    (set! (.-needsUpdate (.-instanceMatrix mesh)) true)
    (when-let [ic (.-instanceColor mesh)] (set! (.-needsUpdate ic) true))))

(defn- write-instances-flat!
  "The same, straight from the position and velocity arrays.

  Going through `to-boids` first would allocate a map and two vectors per
  boid per frame, which is most of what the flat representation is for."
  [^js mesh ^js scratch ^js color {:keys [n ^js pos ^js vel]} [w h d]]
  (let [ox (/ w 2) oy (/ h 2) oz (/ d 2)]
    (dotimes [i n]
      (let [b  (* 3 i)
            x  (- (aget pos b) ox)
            y  (- (aget pos (+ b 1)) oy)
            z  (- (aget pos (+ b 2)) oz)
            vx (aget vel b) vy (aget vel (+ b 1)) vz (aget vel (+ b 2))
            m  (js/Math.sqrt (+ (* vx vx) (* vy vy) (* vz vz)))
            nx (if (pos? m) (/ vx m) 0.0)
            ny (if (pos? m) (/ vy m) 0.0)
            nz (if (pos? m) (/ vz m) 0.0)]
        (.set (.-position scratch) x y z)
        (.lookAt scratch (+ x nx) (+ y ny) (+ z nz))
        (.updateMatrix scratch)
        (.setMatrixAt mesh i (.-matrix scratch))
        (.setHSL color
                 (+ 0.5 (/ (js/Math.atan2 nz nx) (* 2 js/Math.PI)))
                 0.8
                 (+ 0.55 (* 0.18 ny)))
        (.setColorAt mesh i color)))
    (set! (.-count mesh) n)
    (set! (.-needsUpdate (.-instanceMatrix mesh)) true)
    (when-let [ic (.-instanceColor mesh)] (set! (.-needsUpdate ic) true))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 0.5 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        mesh     (THREE/InstancedMesh. (boid-geometry)
                                       (THREE/MeshLambertMaterial.)
                                       capacity)
        scratch  (THREE/Object3D.)
        obstacle-group (reduce (fn [^js g m] (.add g m) g)
                               (THREE/Group.) (obstacle-meshes world))
        color    (THREE/Color.)
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {:flock (let [reference (boids/flock (.-boids controls) world (params))]
                                 (if (flat?) (flat/from-boids reference) reference))})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (.add scene mesh)
    (.add scene (bounds-box world))
    (set! (.-visible obstacle-group) (.-obstacles controls))
    (.add scene obstacle-group)
    (.add scene (THREE/AmbientLight. 0xffffff 0.65))
    (let [sun (THREE/DirectionalLight. 0xfff4e0 0.9)]
      (.set (.-position sun) 80 120 60)
      (.add scene sun))
    (.set (.-position camera) 115 75 115)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(reset-flock! []
                (let [reference (boids/flock (.-boids controls) world (params))]
                  (swap! state assoc :flock
                         (if (flat?) (flat/from-boids reference) reference))))
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
                  (let [t0    (js/performance.now)
                        p     (params)
                        flock (if (flat?)
                                ;; A flat flock is a fixed size, so changing
                                ;; the count is a rebuild rather than a resize.
                                (let [f (:flock @state)]
                                  (if (= (:n f) (.-boids controls))
                                    (flat/step f world p)
                                    (flat/step (flat/from-boids
                                                (boids/flock (.-boids controls) world p))
                                               world p)))
                                (-> (:flock @state)
                                    (resize-to-flock (.-boids controls) (:max-speed p))
                                    (boids/step world p)))]
                    (swap! state assoc :flock flock)
                    (if (flat?)
                      (write-instances-flat! mesh scratch color flock world)
                      (write-instances! mesh scratch color flock world))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "flatArrays") (.onChange reset-flock!))
          (.add gui controls "boids" 10 capacity 10)
          (.add gui controls "separation" 0 3 0.1)
          (.add gui controls "alignment" 0 3 0.1)
          (.add gui controls "cohesion" 0 3 0.1)
          (.add gui controls "perception" 8 50 1)
          (.add gui controls "personalSpace" 2 30 1)
          (.add gui controls "speed" 0.2 2.5 0.1)
          (.add gui controls "avoidance" 0 6 0.1)
          (-> (.add gui controls "obstacles")
              (.onChange (fn [v] (set! (.-visible obstacle-group) v))))
          (.add gui #js {:reset reset-flock!} "reset"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "boids-3d")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
