(ns allgo.demo.rigid
  "WebGL demo for `allgo.physics.rigid`, after Ten Minute Physics 22.

  Boxes linked by distance constraints. The solver is the one from the
  soft body demo with one thing added -- an orientation -- and the whole
  consequence of that addition is visible in the `chain`: pull on a link
  off-centre and it turns, because a correction applied away from the
  centre of mass is partly a torque. A particle solver cannot do that at
  all.

  `substeps` is the parameter to play with, and the lesson is that ten
  small steps solved once each beat one step solved ten times. Drop it to
  1 and the chain goes slack and stretchy; the constraints are being
  solved exactly as often, but each solve is against a position that has
  not been re-integrated.

  `compliance` is inverse stiffness in metres per newton, and it means the
  same thing whatever the step size -- that is the difference between XPBD
  and pushing a fraction of the error each iteration, where stiffness
  quietly depended on how often you pushed.

  `show forces` colours each link by what it is carrying. On the bridge,
  walk the load along and watch where it goes."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.rigid :as rigid]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:scene      "chain"
       :links      12
       :substeps   10
       :compliance 0.0
       :gravity    10.0
       :damping    0.2
       :showForces true})

(defn- chain-scene
  "A row of links, the first one pinned in the air."
  [n]
  (let [w 0.34 h 0.1 gap 0.36
        bodies (vec (for [i (range n)]
                      (rigid/box {:size [w h h] :density 800.0
                                  :pos [(* gap (+ i 0.5)) 2.2 0.0]})))
        cs (vec (for [i (range n)]
                  (if (zero? i)
                    (rigid/distance-constraint
                     bodies {:a 0 :b nil
                             :at [0.0 2.2 0.0] :other-at [0.0 2.2 0.0]
                             :distance 0.0})
                    (rigid/distance-constraint
                     bodies {:a (dec i) :b i
                             :at [(* gap i) 2.2 0.0] :other-at [(* gap i) 2.2 0.0]
                             :distance 0.0}))))]
    {:bodies bodies :constraints cs}))

(defn- bridge-scene
  "Planks linked end to end, pinned at both banks."
  [n]
  (let [w 0.4 h 0.08 span (* w n)
        x0 (* -0.5 span)
        bodies (vec (for [i (range n)]
                      (rigid/box {:size [w h 1.2] :density 500.0
                                  :pos [(+ x0 (* w (+ i 0.5))) 1.6 0.0]})))
        joint (fn [i] (+ x0 (* w i)))
        cs (vec (concat
                 [(rigid/distance-constraint
                   bodies {:a 0 :b nil :at [(joint 0) 1.6 0.0]
                           :other-at [(joint 0) 1.6 0.0] :distance 0.0})]
                 (for [i (range 1 n)]
                   (rigid/distance-constraint
                    bodies {:a (dec i) :b i :at [(joint i) 1.6 0.0]
                            :other-at [(joint i) 1.6 0.0] :distance 0.0}))
                 [(rigid/distance-constraint
                   bodies {:a (dec n) :b nil :at [(joint n) 1.6 0.0]
                           :other-at [(joint n) 1.6 0.0] :distance 0.0})]))]
    {:bodies bodies :constraints cs}))

(defn- pendulum-scene
  "Three arms, pinned at the top. Chaotic, and a good test of whether the
  solver is quietly adding energy -- it should swing forever without ever
  swinging higher."
  [n]
  (let [len 0.7 h 0.08
        bodies (vec (for [i (range n)]
                      (rigid/box {:size [len h h] :density 600.0
                                  :pos [(* len (+ i 0.5)) 2.4 0.0]})))
        cs (vec (for [i (range n)]
                  (let [x (* len i)]
                    (if (zero? i)
                      (rigid/distance-constraint
                       bodies {:a 0 :b nil :at [0.0 2.4 0.0]
                               :other-at [0.0 2.4 0.0] :distance 0.0})
                      (rigid/distance-constraint
                       bodies {:a (dec i) :b i :at [x 2.4 0.0]
                               :other-at [x 2.4 0.0] :distance 0.0})))))]
    {:bodies bodies :constraints cs}))

(defn- make-scene []
  (let [n (.-links controls)]
    (case (.-scene controls)
      "chain"    (chain-scene n)
      "bridge"   (bridge-scene (max 4 n))
      "pendulum" (pendulum-scene (max 2 (min 4 n)))
      (chain-scene n))))

(defn- build! [^js three-scene state]
  (doseq [^js m (:meshes @state)]
    (.remove three-scene m) (.dispose (.-geometry m)) (.dispose (.-material m)))
  (doseq [^js l (:links @state)]
    (.remove three-scene l) (.dispose (.-geometry l)) (.dispose (.-material l)))
  (let [{:keys [bodies constraints]} (make-scene)
        meshes (mapv (fn [b]
                       (let [[sx sy sz] (:size b)
                             m (THREE/Mesh.
                                (THREE/BoxGeometry. sx sy sz)
                                (THREE/MeshPhongMaterial. #js {:color 0xd9b26f}))]
                         (.add three-scene m)
                         m))
                     bodies)
        links (mapv (fn [_]
                      (let [l (THREE/Mesh.
                               (THREE/SphereGeometry. 0.045 12 8)
                               (THREE/MeshBasicMaterial. #js {:color 0xff5544}))]
                        (.add three-scene l)
                        l))
                    constraints)]
    (swap! state assoc
           :world {:bodies bodies :constraints constraints}
           :meshes meshes :links links)))

(defn- force-colour
  "Red for heavily loaded, blue for slack."
  [f peak]
  (let [t (min 1.0 (/ (js/Math.abs f) (max 1e-6 peak)))]
    (THREE/Color. t (* 0.35 (- 1.0 t)) (- 1.0 t))))

(defn- refresh! [{:keys [world meshes links]}]
  (let [{:keys [bodies constraints forces]} world
        peak (reduce (fn [m f] (max m (js/Math.abs f))) 1e-6 (or forces []))]
    (doseq [[^js m b] (map vector meshes bodies)]
      (let [[x y z] (:pos b)
            [qx qy qz qw] (:rot b)]
        (.set (.-position m) x y z)
        (.set (.-quaternion m) qx qy qz qw)))
    (doseq [[k ^js l] (map-indexed vector links)]
      (let [[p0 _] (rigid/endpoints bodies (nth constraints k))
            [x y z] p0]
        (.set (.-position l) x y z)
        (set! (.-visible l) (.-showForces controls))
        (when (and (.-showForces controls) forces)
          (.copy (.. l -material -color)
                 (force-colour (nth forces k 0.0) peak)))))))

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
      (.set (.-position sun) 3 5 4)
      (.add scene sun))
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 40 40)
                             (THREE/MeshPhongMaterial. #js {:color 0x141a2a}))]
      (set! (.-x (.-rotation floor)) (- (/ js/Math.PI 2)))
      (.add scene floor))
    (.add scene (doto (THREE/GridHelper. 12 24 0x2b3450 0x1b2236)
                  (-> .-position (.setY 0.002))))
    (.set (.-position camera) 2.6 2.4 4.4)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0.8 1.4 0)
      (letfn [(rebuild! [] (build! scene state) (refresh! @state))
              (resize! []
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
                        w (:world @state)
                        damping (.-damping controls)
                        w (update w :bodies
                                  (fn [bs] (mapv #(assoc % :damping damping
                                                         :angular-damping damping) bs)))
                        w (update w :constraints
                                  (fn [cs] (mapv #(assoc % :compliance (.-compliance controls)) cs)))
                        w (rigid/step (assoc w
                                             :gravity [0.0 (- (.-gravity controls)) 0.0]
                                             :dt (/ 1.0 60.0)
                                             :substeps (.-substeps controls)))]
                    (swap! state assoc :world w)
                    (refresh! @state)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. resize!) container)
        (resize!)
        (rebuild!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "scene" #js ["chain" "bridge" "pendulum"])
                (.onChange rebuild!))
            (-> (.add controls "links" 2 24 1) (.onFinishChange rebuild!))
            (.add controls "substeps" 1 25 1)
            (.add controls "compliance" 0 0.02 0.0005)
            (.add controls "gravity" 0 20 0.5)
            (.add controls "damping" 0 3 0.1)
            (.add controls "showForces")
            (.add #js {:reset rebuild!} "reset")))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "rigid")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
