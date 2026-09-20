(ns allgo.demo.bricks
  "Shoot a ball at a wall, with the solver switchable underneath.

  The same scene, the same contacts, three different answers to what to
  do about them -- which is the only honest way to compare
  `:sequential-impulse`, `:tgs` and `:xpbd`, because each one's character
  shows in how a *pile* behaves rather than in any single collision.

  ## What to look for

  Build a tall wall, drop `iterations` to three or four, and watch the
  bottom courses. Sequential impulse sags: it linearised at the top of
  the step and its last iteration is still solving the geometry the first
  one saw, so weight leaks downwards through the stack. TGS holds, for
  one reason -- it moves the bodies between iterations and re-measures,
  so the later passes are solving the wall as it now is.

  Turn `iterations` back up and the difference mostly disappears, which
  is the honest summary: TGS buys accuracy per iteration, not accuracy
  you cannot otherwise have.

  XPBD is the calmest of the three at rest and the least like the others
  on impact. It never solves a velocity -- it moves bodies until they
  stop overlapping and reads the velocity back off the movement -- so it
  cannot inject energy, and a wall under it settles dead rather than
  shivering. What it gives up is control over the bounce.

  ## Cost, and where it actually goes

  The wall starts small because this is not yet fast, and it is worth
  being exact about why rather than leaving it to be rediscovered.
  Measured in a browser on a settled nine-by-eight wall -- 72 bricks, 581
  contacts -- a step costs about 142ms, and it splits:

      collision detection   78ms
      everything else       64ms

  So it is `allgo.physics.contact` rather than any of the solvers. The
  separating axis test asks fifteen questions of every touching pair and
  each is a handful of dot and cross products, which on the JVM the
  escape analysis makes nearly free and in JavaScript allocates a
  three-element vector apiece. The solve itself is already written in
  primitive doubles for that reason; the collision detection is not.

  Turning the iteration count down does not help -- it makes things
  *slower*, because a wall that is not held up spreads out and touches
  more. That is the clearest evidence that the solver is not the
  bottleneck.

  About twenty bricks runs at fifty frames a second, and the sliders go
  up from there if you want to watch it struggle. What would fix it is
  the same treatment the solve already had, plus a broad phase that is
  not every pair against every other -- `allgo.spatial.hash` and
  `allgo.spatial.sweep` are both sitting there unused."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.solver :as solver]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private brick-size [0.9 0.45 0.45])
(def ^:private ball-radius 0.35)

(def ^:private ^js controls
  #js {:solver "tgs"
       :rows 5
       :cols 4
       :iterations 8
       :substeps 4
       :friction 0.55
       :restitution 0.0
       :speed 26
       :fire (fn [])
       :reset (fn [])})

(defn- wall
  "A floor and a running-bond wall of bricks standing on it."
  [rows cols]
  (let [[bw bh bd] brick-size]
    (vec (cons (rigid/box {:pos [0.0 -0.5 0.0] :size [60.0 1.0 30.0]})
               (for [r (range rows)
                     c (range cols)
                     ;; Every other course offset by half a brick, which
                     ;; is what stops the wall being a set of independent
                     ;; columns with nothing tying them together.
                     :let [offset (if (odd? r) (* 0.5 bw) 0.0)
                           x (+ (* (- c (* 0.5 (dec cols))) bw) offset)]]
                 (rigid/box {:pos [x (+ (* 0.5 bh) (* r bh)) 0.0]
                             :size [bw bh bd]
                             :density 1.6}))))))

(defn- projectile [speed]
  (rigid/ball {:pos [0.0 2.2 14.0]
               :radius ball-radius
               :density 7.8
               :vel [0.0 1.0 (- (double speed))]}))

(defn- body-mesh [b]
  (let [mat (THREE/MeshStandardMaterial.
             #js {:color (if (rigid/static? b) 0x2b3040 0xb98d5f)
                  :roughness 0.82 :metalness 0.05})]
    (if (= :ball (:shape b))
      (THREE/Mesh. (THREE/SphereGeometry. (:radius b) 20 14)
                   (THREE/MeshStandardMaterial.
                    #js {:color 0x9fb6d4 :roughness 0.35 :metalness 0.5}))
      (let [[sx sy sz] (:size b)]
        (THREE/Mesh. (THREE/BoxGeometry. sx sy sz) mat)))))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 55 1 0.1 500)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom {:world nil :meshes [] :cost 0.0})]
    (set! (.-background scene) (THREE/Color. 0x0b0d15))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.add scene (THREE/AmbientLight. 0xffffff 0.55))
    (let [sun (THREE/DirectionalLight. 0xfff1d8 1.25)]
      (.set (.-position sun) 8.0 14.0 10.0)
      (.add scene sun))
    (.set (.-position camera) 9.0 5.5 15.0)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0.0 1.8 0.0)
      (letfn [(sync-meshes! []
                (let [{:keys [world meshes]} @state]
                  (doseq [[^js m b] (map vector meshes (:bodies world))]
                    (let [[x y z] (:pos b)
                          [qx qy qz qw] (:rot b)]
                      (.set (.-position m) x y z)
                      (.set (.-quaternion m) qx qy qz qw)))))
              (rebuild! []
                (doseq [^js m (:meshes @state)]
                  (.remove scene m)
                  (.dispose (.-geometry m))
                  (.dispose (.-material m)))
                (let [bodies (wall (long (.-rows controls)) (long (.-cols controls)))
                      meshes (mapv body-mesh bodies)]
                  (doseq [^js m meshes] (.add scene m))
                  (swap! state assoc
                         :world (solver/world bodies (solver-opts))
                         :meshes meshes
                         :cost 0.0)
                  (sync-meshes!)))
              (solver-opts []
                {:solver (keyword (.-solver controls))
                 :iterations (long (.-iterations controls))
                 :substeps (long (.-substeps controls))
                 :friction (double (.-friction controls))
                 :restitution (double (.-restitution controls))})
              (fire! []
                (let [b (projectile (.-speed controls))
                      m (body-mesh b)]
                  (.add scene m)
                  (swap! state (fn [s]
                                 (-> s
                                     (update-in [:world :bodies] conj b)
                                     (update :meshes conj m))))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (swap! state update :world
                           (fn [w] (solver/step (merge w (solver-opts)) (/ 1.0 60.0))))
                    (let [t1 (js/performance.now)]
                      (swap! state update :cost #(+ (* 0.9 %) (* 0.1 (- t1 t0)))))
                    (sync-meshes!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (set! (.-fire controls) fire!)
        (set! (.-reset controls) rebuild!)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (let [gui (GUI. #js {:container container})]
          (doto gui
            (-> (.add controls "solver" #js ["sequential-impulse" "tgs" "xpbd"]))
            (-> (.add controls "iterations" 1 20 1))
            (-> (.add controls "substeps" 1 8 1))
            (-> (.add controls "friction" 0.0 1.2 0.05))
            (-> (.add controls "restitution" 0.0 0.8 0.05))
            (-> (.add controls "speed" 8 60 1))
            (.add controls "fire")
            (-> (.add controls "rows" 3 16 1) (.onFinishChange rebuild!))
            (-> (.add controls "cols" 3 16 1) (.onFinishChange rebuild!))
            (.add controls "reset")))
        {:start (fn [] (when-not @running?
                         (reset! running? true)
                         (on-resize)
                         (animate)))
         :stop (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "bricks")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
