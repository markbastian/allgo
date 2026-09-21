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

  - **column** -- one brick on another, `courses` high. Nothing is asked
    of it but to hold still, which is why it is the sharpest test in
    here: see below.
  - **jenga** -- three blocks to a level, `levels` of them, each level
    turned ninety degrees from the one below. Wide enough not to lean
    like the column does, so it asks a different question: whether a
    solver can rest a block across two others and leave the third alone.
  - **keep** -- a round tower inside a half circle of rampart, which is
    the board of Crossbows and Catapults laid out in bricks. Both are
    rings of boxes rather than anything curved, so every contact in it
    is between two faces that are not parallel -- the one scene here
    that is not built out of right angles.
  - **wall** -- a running-bond wall, `rows` by `cols`. The stack test:
    how much weight a solver can hold up without letting the bottom sag.

  ## The ring, which reverses the answer

  The keep is the only scene here not built out of right angles: a ring
  of boxes meets its neighbours face to face only at the inner edge and
  fans apart outwards, so almost every contact in it is between two
  planes at an angle.

  At the size it starts -- five courses of keep inside two of rampart,
  59 blocks -- that costs nothing. Every block is still where it was
  laid after a minute and the whole scene is asleep, under all three
  solvers, for about 0.4ms a step.

  One course more of each is where they part company. Six and three, 77
  blocks, after a minute:

      sequential impulse   23 of 48 keep blocks still laid; the top
                           three courses have unwound off the tower
      tgs                  all 48 laid, but the rampart never goes
                           quiet, so nothing sleeps -- 5.3ms a step
      xpbd                 all 48 laid, everything asleep, 0.5ms

  Which is the opposite of the Jenga answer below, where XPBD is the
  one that cannot settle and TGS is the one that can. Neither result is
  the general one: they are two different scenes and the character that
  wins in each is different.

  Two of 29 rampart blocks come off under all three, and that is not a
  fault: an arc bonded course over course leaves a block at each end of
  the top course with a joint under it, and it falls off the same way
  the wall's end brick does.

  ## Which of them can leave a Jenga tower alone

  The tower is wide enough that none of them has trouble holding it up
  -- twelve levels stands under all three, and only sequential impulse
  gives out at eighteen. What it separates them on instead is quiet.
  Twelve levels, thirty-six blocks, left alone for a minute:

      sequential impulse   never sleeps, all 36 awake at 60s
      tgs                  asleep at 2.2s
      xpbd                 never sleeps, all 36 awake at 60s

  Turn sleeping off and the reason is plain -- mean block speed once it
  has stopped going anywhere, against the 0.01 m/s the sleep test wants:

      sequential impulse   0.028
      tgs                  0.007
      xpbd                 0.009

  XPBD misses by a hair and pays the full price, because an island
  sleeps only when everything in it is still: eleven of thirty-six
  blocks over the line is as bad as all of them. Sequential impulse is
  not close. And TGS's two seconds are partly self-fulfilling -- it is
  quiet enough early to freeze, and frozen is free, so the tremble that
  would have built up over the next minute never does.

  ## How tall a column each of them will hold

  A column is the demo's hardest question, and it is worth being clear
  about why. Its bricks are laid square, so in exact arithmetic it would
  stand for ever -- every contact is symmetric and there is nothing to
  tip it. What it actually does is lean, slowly, along its narrow axis,
  and past a certain height the lean runs away with it. The seconds
  below are how long it stood before the top brick left the footprint,
  watched for twenty, gripping at 0.55:

      courses             5    6    7    8    9   10   11   12
      sequential impulse  --   --   --  4.8  2.8  2.7  2.5  2.9
      tgs                 --   --   --   --   --   --   --  7.3
      xpbd                --   --   --   --   --  8.4  5.3  4.2

  So seven courses for sequential impulse, nine for XPBD and eleven for
  TGS, and the demo starts at eight -- tall enough that switching to
  sequential impulse knocks it down while you watch.

  It is drift, not a kick. Trace the top brick and it walks one way in
  millimetres a second and never comes back; the column is standing on
  its own rounding error, and once the centre of mass is outside the
  footprint gravity does the rest. Iterations do not buy it back --
  sixteen and four fall at the same height -- because the error is not
  in how well each step is solved but in the fact that there are sixty
  of them a second, each leaving a little behind.

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
  Measured on the JVM on a settled sixteen-by-sixteen wall -- 257 bodies,
  around 2,400 contacts -- a step and where it goes:

                            sequential impulse    tgs    xpbd
      total (ms)                         14.7   19.6    32.9
      collision detection                 6.8    6.8     7.5
      the constraint solve                2.2    3.0    16.4
      everything else                     5.7    9.8     9.0

  So it is `allgo.physics.contact` rather than any of the solvers -- the
  solve is a sixth of a sequential impulse step. The separating axis test
  asks fifteen questions of every touching pair and each is a handful of
  dot and cross products, which on the JVM the escape analysis makes
  nearly free and in JavaScript allocates a three-element vector apiece.
  The solve is already written in primitive doubles for that reason; the
  collision detection is not, and that is the next thing.

  It was not always the top of the list. Until recently most of a step
  was neither detection nor solve but the cost of carrying bodies and
  contacts through persistent maps and vectors -- the XPBD position pass
  rebuilding a 257 element vector of body maps twice per contact, most of
  all. `allgo.physics.solver` keeps the poses in flat arrays now and the
  three solvers went from 19.9 / 32.9 / 91.0 ms to the table above.

  Turning the iteration count down does not help -- it makes things
  *slower*, because a wall that is not held up spreads out and touches
  more. That is the clearest evidence that the solver is not the
  bottleneck.

  These are JVM numbers, from `make test`'s JDK rather than from a
  browser. The browser is slower and has not been re-measured since; a
  tab that is not in front is throttled about sixteenfold, so anything
  timed there has to be taken as a ratio against a calibration loop
  rather than read off the clock."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.quaternion :as q]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.solver :as solver]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private brick-size [0.9 0.45 0.45])

(def ^:private jenga-size
  "A Jenga block, at the real 5 : 1.67 : 1 -- three of them laid side by
  side are exactly as wide as one is long, which is what lets the tower
  turn ninety degrees every level and still sit square."
  [1.2 0.24 0.4])
(def ^:private ball-radius 0.35)

(defn- floor []
  (rigid/box {:pos [0.0 -0.5 0.0] :size [60.0 1.0 30.0]}))

(defn- brick
  "One brick of `size`, standing at `pos`, optionally turned by `rot`."
  ([pos size] (brick pos size nil))
  ([pos size rot]
   (rigid/box (cond-> {:pos pos :size size :density 1.6}
                rot (assoc :rot rot)))))

(def ^:private keep-radius 1.3)
(def ^:private rampart-radius 3.2)

(defn- chord
  "How long a block may be to sit in a ring and still touch its
  neighbours rather than start inside them.

  A ring of boxes is wedges the wrong way round: the blocks meet at
  their inner faces and fan apart towards the outside. So the length to
  take is the chord at the inner radius, `step` radians apart. Taking it
  at the centres instead buries each block a few millimetres in the one
  beside it, and a ring laid like that does not settle -- it springs."
  ^double [^double radius ^double depth ^double step]
  (* 2.0 (- radius (* 0.5 depth)) (Math/sin (* 0.5 step))))

(defn- ring-brick
  "A block `theta` radians round a ring of `radius`, lying along it.

  Turned about the upright so that its length runs along the tangent and
  its depth points out from the centre, which is a quarter turn less
  `theta` -- at `theta` of zero the block faces down +x and has to be
  turned ninety degrees to get there."
  [^double radius ^double theta ^double y [len height depth]]
  (brick [(* radius (Math/cos theta)) y (* radius (Math/sin theta))]
         [len height depth]
         (q/from-axis-angle [0.0 1.0 0.0] (- (* 0.5 Math/PI) theta))))

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

(defmethod build-scene "keep" [_ ^js c]
  (let [courses (long (.-keepCourses c))
        rampart (long (.-rampartCourses c))
        [_ bh bd] brick-size
        ;; Eight blocks to a course round the keep, ten to a course
        ;; round the half circle of the rampart -- which at these two
        ;; radii comes out at very nearly the same block for both.
        keep-step (/ (* 2.0 Math/PI) 8)
        wall-step (/ Math/PI 10)
        keep-len (chord keep-radius bd keep-step)
        wall-len (chord rampart-radius bd wall-step)]
    {:bodies (-> [(floor)]
                 (into (for [r (range courses)
                             i (range 8)
                             ;; Half a block round on every other
                             ;; course, for the same reason the wall
                             ;; does it: without the bond the keep is
                             ;; eight independent columns.
                             :let [theta (* keep-step (+ i (if (odd? r) 0.5 0.0)))]]
                         (ring-brick keep-radius theta (+ (* 0.5 bh) (* r bh))
                                     [keep-len bh bd])))
                 (into (for [r (range rampart)
                             ;; The rampart is an arc, not a ring, so
                             ;; the bonded course cannot just be shifted
                             ;; -- that would hang a block off each end
                             ;; over nothing. It drops one block instead
                             ;; and sits half a block in at both ends.
                             :let [odd? (odd? r)
                                   n (if odd? 9 10)]
                             i (range n)
                             :let [theta (* wall-step (+ i (if odd? 1.0 0.5)))]]
                         (ring-brick rampart-radius theta (+ (* 0.5 bh) (* r bh))
                                     [wall-len bh bd]))))
     :span (+ rampart-radius (* 0.5 bd))
     :height (* (max courses rampart) bh)}))

(defmethod build-scene "jenga" [_ ^js c]
  (let [levels (long (.-jengaLevels c))
        [bl bh bw] jenga-size]
    {:bodies (into [(floor)]
                   (for [l (range levels)
                         i [-1 0 1]
                         ;; Turned ninety degrees from the level below.
                         ;; A box's size carries its orientation here --
                         ;; swapping the two horizontal extents is the
                         ;; same body as turning it, and spares the
                         ;; contacts a rotation they cannot use.
                         :let [across? (odd? l)
                               o (* i bw)
                               y (+ (* 0.5 bh) (* l bh))]]
                     (if across?
                       (brick [o y 0.0] [bw bh bl])
                       (brick [0.0 y o] [bl bh bw]))))
     :span (* 0.5 bl)
     :height (* levels bh)}))

(defmethod build-scene "column" [_ ^js c]
  (let [courses (long (.-columnCourses c))
        [bw bh bd] brick-size]
    {:bodies (into [(floor)]
                   (for [r (range courses)]
                     (brick [0.0 (+ (* 0.5 bh) (* r bh)) 0.0] [bw bh bd])))
     :span (* 0.5 bw)
     :height (* courses bh)}))

(def ^:private scene-names
  "In the order they are offered, simplest first."
  ["column" "jenga" "keep" "wall"])

(def ^:private ^js controls
  #js {:scene "wall"
       :solver "tgs"
       :rows 5
       :cols 4
       :columnCourses 8
       :jengaLevels 12
       :keepCourses 5
       :rampartCourses 2
       :iterations 8
       :substeps 4
       :friction 0.55
       :restitution 0.0
       :speed 26
       :fire (fn [])
       :reset (fn [])})

(def ^:private scene-sliders
  "The sliders each scene brings with it, `[property min max step label]`.

  Everything not here -- the solver, the material, the shot -- is common
  to all of them and always shown. A property belongs to exactly one
  scene, so that each can pick its own range without having to agree
  with anybody else about what, say, a sensible height is; `label` is
  how it gets a plain name back in the GUI regardless."
  {"column" [["columnCourses" 2 40 1 "courses"]]
   "jenga"  [["jengaLevels" 3 24 1 "levels"]]
   "keep"   [["keepCourses" 2 14 1 "keep courses"]
             ["rampartCourses" 1 8 1 "rampart courses"]]
   "wall"   [["rows" 3 16 1] ["cols" 3 16 1]]})

(defn- projectile
  "A shot high on the stack, from far enough out to see it coming.

  High rather than square in the middle: a stack that is about to fail
  fails from the top, so that is where the interesting answer is."
  [speed {:keys [span height]}]
  (rigid/ball {:pos [0.0 (max 1.0 (* 0.72 height)) (+ 8.0 (* 1.6 span))]
               :radius ball-radius
               :density 7.8
               :vel [0.0 (* 0.04 speed) (- (double speed))]}))

(defn- clay
  "The colour of brick `i`, jittered a little about the same clay.

  Not decoration. A Jenga tower is blocks that fit exactly, so in one
  flat colour it renders as a single brown slab and the thing the scene
  is about -- that each level lies across the one below -- is invisible.
  The step is the plastic constant's reciprocal, which is the
  one-dimensional golden ratio's better-behaved cousin: consecutive
  indices land about as far apart as they can, so no two bricks laid
  next to each other come out the same shade, and it is arithmetic
  rather than a random number, so a scene looks the same every time it
  is built."
  [^long i]
  (let [jitter (- (mod (* (inc i) 0.7548776662) 1.0) 0.5)]
    ;; Brightness only. Turning the hue as well looked like three kinds
    ;; of brick rather than one kind badly fired, and the seams read
    ;; just as clearly without it.
    (.multiplyScalar (THREE/Color. 0xb98d5f) (+ 1.0 (* 0.45 jitter)))))

(defn- body-mesh [i b]
  (if (= :ball (:shape b))
    (THREE/Mesh. (THREE/SphereGeometry. (:radius b) 20 14)
                 (THREE/MeshStandardMaterial.
                  #js {:color 0x9fb6d4 :roughness 0.35 :metalness 0.5}))
    (let [[sx sy sz] (:size b)]
      (THREE/Mesh. (THREE/BoxGeometry. sx sy sz)
                   (THREE/MeshStandardMaterial.
                    #js {:color (if (rigid/static? b) (THREE/Color. 0x2b3040) (clay i))
                         :roughness 0.82 :metalness 0.05})))))

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
                      meshes (vec (map-indexed body-mesh bodies))]
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
                      m (body-mesh 0 b)]
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
                                    [prop lo hi step label] specs]
                                [prop (-> (.add gui controls prop lo hi step)
                                          (.name (or label prop))
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
