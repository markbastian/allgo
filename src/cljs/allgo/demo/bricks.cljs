(ns allgo.demo.bricks
  "Shoot a ball at a stack of bricks, with the solver switchable underneath.

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

  ## The scenes

  `scene` picks what is standing there to be knocked down. Each one is a
  different question to ask a solver, and each brings its own controls
  and its own camera:

  - **wall** -- a running-bond wall, `rows` by `cols`. The stack test:
    how much weight a solver can hold up without letting the bottom sag.

  ## How tall a wall each of them will hold

  Turn the rows up and they part company. Four bricks wide, gripping at
  0.55, left alone for twenty seconds, counting the bricks still where
  they were laid -- one short is full marks, because a running bond ends
  every course with a brick half over the edge and that one falls off by
  itself:

      courses            4      6      8     12     16
      sequential impulse 15/16  23/24  31/32  10/48   0/64
      tgs                15/16  23/24  31/32  47/48  63/64
      xpbd               15/16  23/24  31/32   2/48   0/64

  Up to eight courses it makes no difference which is selected. Past
  that, TGS is the only one that holds, because it is the only one that
  substeps: support reaches the top of the wall within the step instead
  of a course per step. Sequential impulse cannot be iterated into it --
  sixteen iterations and sixty-four give the same nineteen bricks -- and
  XPBD gives out at the same height for its own reason. Neither is
  broken; the difference is the thing the demo is for.

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

(defn- floor []
  (rigid/box {:pos [0.0 -0.5 0.0] :size [60.0 1.0 30.0]}))

(defn- brick
  "One brick of `size`, standing at `pos`, optionally turned by `rot`."
  ([pos size] (brick pos size nil))
  ([pos size rot]
   (rigid/box (cond-> {:pos pos :size size :density 1.6}
                rot (assoc :rot rot)))))

(defmulti build-scene
  "The named scene to knock down, read off the GUI `controls`.

  Returns `{:bodies :span :height}`: the bodies with the floor first, and
  how far the stack reaches sideways and upwards. Those two numbers are
  all the camera and the projectile need -- a scene says how big it is
  and everything else is framed from that, so adding one does not mean
  also picking a viewpoint and a firing line by hand."
  (fn [name _controls] name))

(defmethod build-scene "wall" [_ ^js c]
  (let [rows (long (.-rows c))
        cols (long (.-cols c))
        [bw bh bd] brick-size]
    {:bodies (into [(floor)]
                   (for [r (range rows)
                         c (range cols)
                         ;; Every other course offset by half a brick,
                         ;; which is what stops the wall being a set of
                         ;; independent columns with nothing tying them
                         ;; together.
                         :let [offset (if (odd? r) (* 0.5 bw) 0.0)
                               x (+ (* (- c (* 0.5 (dec cols))) bw) offset)]]
                     (brick [x (+ (* 0.5 bh) (* r bh)) 0.0] [bw bh bd])))
     :span (* 0.5 cols bw)
     :height (* rows bh)}))

(def ^:private scene-names
  "In the order they are offered, simplest first."
  ["wall"])

(def ^:private ^js controls
  #js {:scene "wall"
       :solver "tgs"
       :rows 5
       :cols 4
       :iterations 8
       :substeps 4
       :friction 0.55
       :restitution 0.0
       :speed 26
       :fire (fn [])
       :reset (fn [])})

(def ^:private scene-sliders
  "The sliders each scene brings with it, as `[property min max step]`.

  Everything not here -- the solver, the material, the shot -- is common
  to all of them and always shown. A property belongs to exactly one
  scene, so that each can pick its own range without having to agree
  with anybody else about what, say, a sensible height is."
  {"wall" [["rows" 3 16 1] ["cols" 3 16 1]]})

(defn- projectile
  "A shot high on the stack, from far enough out to see it coming.

  High rather than square in the middle: a stack that is about to fail
  fails from the top, so that is where the interesting answer is."
  [speed {:keys [span height]}]
  (rigid/ball {:pos [0.0 (max 1.0 (* 0.72 height)) (+ 8.0 (* 1.6 span))]
               :radius ball-radius
               :density 7.8
               :vel [0.0 (* 0.04 speed) (- (double speed))]}))

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
        state (atom {:world nil :meshes [] :cost 0.0 :framing nil})]
    (set! (.-background scene) (THREE/Color. 0x0b0d15))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.add scene (THREE/AmbientLight. 0xffffff 0.55))
    (let [sun (THREE/DirectionalLight. 0xfff1d8 1.25)]
      (.set (.-position sun) 8.0 14.0 10.0)
      (.add scene sun))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (letfn [(sync-meshes! []
                (let [{:keys [world meshes]} @state]
                  (doseq [[^js m b] (map vector meshes (:bodies world))]
                    (let [[x y z] (:pos b)
                          [qx qy qz qw] (:rot b)]
                      (.set (.-position m) x y z)
                      (.set (.-quaternion m) qx qy qz qw)))))
              (look-at-scene! []
                ;; Stand back far enough to see all of it, and look at
                ;; the middle rather than the floor.
                (let [{:keys [span height]} (:framing @state)
                      d (+ 9.0 (* 1.7 span) (* 0.8 height))]
                  (.set (.-position camera) (* 0.6 d) (+ 2.0 (* 1.3 height)) d)
                  (.set (.-target orbit) 0.0 (* 0.55 height) 0.0)
                  (.update orbit)))
              (rebuild! []
                (doseq [^js m (:meshes @state)]
                  (.remove scene m)
                  (.dispose (.-geometry m))
                  (.dispose (.-material m)))
                (let [{:keys [bodies] :as built} (build-scene (.-scene controls) controls)
                      meshes (mapv body-mesh bodies)]
                  (doseq [^js m meshes] (.add scene m))
                  (swap! state assoc
                         :world (solver/world bodies (solver-opts))
                         :meshes meshes
                         :cost 0.0
                         :framing (select-keys built [:span :height]))
                  (sync-meshes!)))
              (solver-opts []
                {:solver (keyword (.-solver controls))
                 :iterations (long (.-iterations controls))
                 :substeps (long (.-substeps controls))
                 :friction (double (.-friction controls))
                 :restitution (double (.-restitution controls))})
              (fire! []
                (let [b (projectile (.-speed controls) (:framing @state))
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
        (look-at-scene!)
        (let [gui (GUI. #js {:container container})
              ^js scene-ctl (.add gui controls "scene" (clj->js scene-names))]
          (doto gui
            (-> (.add controls "solver" #js ["sequential-impulse" "tgs" "xpbd"]))
            (-> (.add controls "iterations" 1 20 1))
            (-> (.add controls "substeps" 1 8 1))
            (-> (.add controls "friction" 0.0 1.2 0.05))
            (-> (.add controls "restitution" 0.0 0.8 0.05))
            (-> (.add controls "speed" 8 60 1))
            (.add controls "fire")
            (.add controls "reset"))
          ;; Every scene's sliders are built once and all but the
          ;; current scene's are hidden, which keeps a rebuild out of the
          ;; business of tearing down and rebuilding the GUI.
          (let [sliders (into {}
                              (for [[_ specs] scene-sliders
                                    [prop lo hi step] specs]
                                [prop (-> (.add gui controls prop lo hi step)
                                          (.onFinishChange rebuild!))]))
                show-sliders! (fn []
                                (let [mine (into #{} (map first)
                                                 (scene-sliders (.-scene controls)))]
                                  (doseq [[prop ^js ctl] sliders]
                                    (.show ctl (contains? mine prop)))))]
            (.onChange scene-ctl (fn [& _]
                                   (rebuild!)
                                   (look-at-scene!)
                                   (show-sliders!)))
            (show-sliders!)))
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
