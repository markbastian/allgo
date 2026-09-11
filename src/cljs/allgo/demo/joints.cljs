(ns allgo.demo.joints
  "WebGL demo for `allgo.physics.joint`, after Ten Minute Physics 25.

  A joint is two frames -- one fixed to each body -- and a statement about
  how they may differ. Every kind here is built from three primitives:
  pin the origins together, drive one frame's axis onto the other's, and
  hold an angle between bounds. A ball joint pins; a hinge also aligns,
  so only one rotation survives; a servo is a hinge with a target angle; a
  motor is a servo whose target keeps moving.

  `scene` walks through them. The `ragdoll` is the one that shows why
  limits matter: turn `limits` off and the joints still hold it together,
  but it folds through itself in ways a body does not.

  `substeps` is worth moving on the `chain` scene. A joint is a hard
  constraint, and how well it holds is how many substeps it gets: the
  anchor gap falls roughly by half every time you double them.

  `show frames` draws each joint's two frames. Where they are drawn apart
  is the constraint error you are watching the solver work off."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as q]
            [allgo.physics.joint :as joint]
            [allgo.physics.rigid :as rigid]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private ^js controls
  #js {:scene      "ragdoll"
       :substeps   15
       :gravity    10.0
       :motorSpeed 3.0
       :limits     true
       :showFrames false})

(def ^:private up-axis
  "A joint frame whose x axis points along world z, so a hinge made with
  it swings in the xy plane."
  (q/from-axis-angle [0 1 0] (- (/ js/Math.PI 2))))

(defn- ragdoll-scene []
  (let [limits? (.-limits controls)
        anchor (rigid/box {:size [0.3 0.3 0.3] :pos [0.0 3.0 0.0]})
        torso  (rigid/box {:size [0.45 0.7 0.28] :density 900.0 :pos [0.0 2.5 0.0]})
        head   (rigid/ball {:radius 0.19 :density 900.0 :pos [0.0 3.05 0.0]})
        limb   (fn [x y w h] (rigid/box {:size [w h 0.16] :density 900.0 :pos [x y 0.0]}))
        upper-l (limb -0.42 2.55 0.5 0.16)
        lower-l (limb -0.92 2.55 0.5 0.16)
        upper-r (limb 0.42 2.55 0.5 0.16)
        lower-r (limb 0.92 2.55 0.5 0.16)
        thigh-l (limb -0.14 1.85 0.18 0.6)
        shin-l  (limb -0.14 1.25 0.18 0.6)
        thigh-r (limb 0.14 1.85 0.18 0.6)
        shin-r  (limb 0.14 1.25 0.18 0.6)
        bodies [anchor torso head upper-l lower-l upper-r lower-r
                thigh-l shin-l thigh-r shin-r]
        swing (fn [opts] (if limits? opts (dissoc opts :swing-min :swing-max)))
        cs [(joint/ball bodies (swing {:a 0 :b 1 :at [0.0 2.85 0.0]
                                       :swing-min -0.35 :swing-max 0.35}))
            (joint/ball bodies (swing {:a 1 :b 2 :at [0.0 2.88 0.0]
                                       :swing-min -0.5 :swing-max 0.5}))
            (joint/ball bodies (swing {:a 1 :b 3 :at [-0.2 2.55 0.0]
                                       :swing-min -1.2 :swing-max 1.2}))
            (joint/hinge bodies (swing {:a 3 :b 4 :at [-0.67 2.55 0.0] :rot up-axis
                                        :swing-min -2.2 :swing-max 0.0}))
            (joint/ball bodies (swing {:a 1 :b 5 :at [0.2 2.55 0.0]
                                       :swing-min -1.2 :swing-max 1.2}))
            (joint/hinge bodies (swing {:a 5 :b 6 :at [0.67 2.55 0.0] :rot up-axis
                                        :swing-min 0.0 :swing-max 2.2}))
            (joint/ball bodies (swing {:a 1 :b 7 :at [-0.14 2.15 0.0]
                                       :swing-min -0.9 :swing-max 0.9}))
            (joint/hinge bodies (swing {:a 7 :b 8 :at [-0.14 1.55 0.0] :rot up-axis
                                        :swing-min -2.0 :swing-max 0.0}))
            (joint/ball bodies (swing {:a 1 :b 9 :at [0.14 2.15 0.0]
                                       :swing-min -0.9 :swing-max 0.9}))
            (joint/hinge bodies (swing {:a 9 :b 10 :at [0.14 1.55 0.0] :rot up-axis
                                        :swing-min -2.0 :swing-max 0.0}))]]
    {:bodies bodies :constraints cs}))

(defn- chain-scene []
  (let [n 10
        anchor (rigid/box {:size [0.2 0.2 0.2] :pos [0.0 3.0 0.0]})
        links (vec (for [i (range n)]
                     (rigid/box {:size [0.18 0.34 0.18] :density 900.0
                                 :pos [0.0 (- 2.8 (* 0.36 i)) 0.0]})))
        bodies (into [anchor] links)
        cs (vec (for [i (range n)]
                  (joint/ball bodies {:a (if (zero? i) 0 i)
                                      :b (inc i)
                                      :at [0.0 (- 2.98 (* 0.36 i)) 0.0]})))]
    {:bodies bodies :constraints cs}))

(defn- motor-scene []
  (let [base (rigid/box {:size [0.5 0.3 0.5] :pos [0.0 1.2 0.0]})
        arm  (rigid/box {:size [1.2 0.14 0.14] :density 800.0 :pos [0.6 1.2 0.0]})
        fore (rigid/box {:size [0.9 0.12 0.12] :density 800.0 :pos [1.65 1.2 0.0]})
        bodies [base arm fore]
        cs [(joint/motor bodies {:a 0 :b 1 :at [0.0 1.2 0.0] :rot up-axis
                                 :velocity (.-motorSpeed controls)})
            (joint/hinge bodies {:a 1 :b 2 :at [1.2 1.2 0.0] :rot up-axis
                                 :swing-min -1.4 :swing-max 1.4})]]
    {:bodies bodies :constraints cs}))

(defn- slider-scene []
  (let [base (rigid/box {:size [0.4 0.4 0.4] :pos [0.0 2.6 0.0]})
        car  (rigid/box {:size [0.5 0.3 0.3] :density 900.0 :pos [0.0 2.0 0.0]})
        bob  (rigid/box {:size [0.3 0.3 0.3] :density 900.0 :pos [0.7 2.0 0.0]})
        bodies [base car bob]
        cs [(joint/prismatic bodies {:a 0 :b 1 :at [0.0 2.0 0.0]
                                     :rot (q/from-axis-angle [0 0 1] (/ js/Math.PI 2))
                                     :distance-min -0.9 :distance-max 0.9})
            (joint/ball bodies {:a 1 :b 2 :at [0.0 2.0 0.0] :target-distance 0.7})]]
    {:bodies bodies :constraints cs}))

(defn- make-scene []
  (case (.-scene controls)
    "ragdoll" (ragdoll-scene)
    "chain"   (chain-scene)
    "motor"   (motor-scene)
    "slider"  (slider-scene)
    (ragdoll-scene)))

(defn- body-mesh [b]
  (let [mat (THREE/MeshPhongMaterial.
             #js {:color (if (rigid/static? b) 0x6b7280 0xd9a15c)})]
    (if (= :ball (:shape b))
      (THREE/Mesh. (THREE/SphereGeometry. (:radius b) 20 16) mat)
      (let [[sx sy sz] (:size b)]
        (THREE/Mesh. (THREE/BoxGeometry. sx sy sz) mat)))))

(defn- build! [^js three-scene state]
  (doseq [^js m (concat (:meshes @state) (:frames @state))]
    (.remove three-scene m)
    (when (.-geometry m) (.dispose (.-geometry m)))
    (when (.-material m) (.dispose (.-material m))))
  (let [{:keys [bodies constraints]} (make-scene)
        meshes (mapv (fn [b] (let [m (body-mesh b)] (.add three-scene m) m)) bodies)
        frames (vec (mapcat (fn [_]
                              (for [colour [0x55ff88 0xff5588]]
                                (let [m (THREE/Mesh.
                                         (THREE/SphereGeometry. 0.05 10 8)
                                         (THREE/MeshBasicMaterial. #js {:color colour}))]
                                  (.add three-scene m)
                                  m)))
                            constraints))]
    (swap! state assoc
           :world {:bodies bodies :constraints constraints}
           :meshes meshes :frames frames)))

(defn- refresh! [{:keys [world meshes frames]}]
  (let [{:keys [bodies constraints]} world]
    (doseq [[^js m b] (map vector meshes bodies)]
      (let [[x y z] (:pos b)
            [qx qy qz qw] (:rot b)]
        (.set (.-position m) x y z)
        (.set (.-quaternion m) qx qy qz qw)))
    (doseq [[k c] (map-indexed vector constraints)]
      (let [[p0 _ p1 _] (joint/frames bodies c)
            ^js m0 (nth frames (* 2 k) nil)
            ^js m1 (nth frames (inc (* 2 k)) nil)]
        (when (and m0 m1)
          (set! (.-visible m0) (.-showFrames controls))
          (set! (.-visible m1) (.-showFrames controls))
          (.set (.-position m0) (nth p0 0) (nth p0 1) (nth p0 2))
          (.set (.-position m1) (nth p1 0) (nth p1 1) (nth p1 2)))))))

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
      (.set (.-position sun) 3 6 4)
      (.add scene sun))
    (let [floor (THREE/Mesh. (THREE/PlaneGeometry. 40 40)
                             (THREE/MeshPhongMaterial. #js {:color 0x141a2a}))]
      (set! (.-x (.-rotation floor)) (- (/ js/Math.PI 2)))
      (.add scene floor))
    (.add scene (doto (THREE/GridHelper. 12 24 0x2b3450 0x1b2236)
                  (-> .-position (.setY 0.002))))
    (.set (.-position camera) 0.2 2.4 5.0)
    (.appendChild container (.-domElement renderer))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 2.0 0)
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
                        w (update w :constraints
                                  (fn [cs]
                                    (mapv #(if (= :motor (:kind %))
                                             (assoc % :velocity (.-motorSpeed controls))
                                             %)
                                          cs)))
                        w (joint/step (assoc w
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
            (-> (.add controls "scene" #js ["ragdoll" "chain" "motor" "slider"])
                (.onChange rebuild!))
            (.add controls "substeps" 2 40 1)
            (.add controls "gravity" 0 20 0.5)
            (.add controls "motorSpeed" -8 8 0.5)
            (-> (.add controls "limits") (.onChange rebuild!))
            (.add controls "showFrames")
            (.add #js {:reset rebuild!} "reset")))
        {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "joints")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
