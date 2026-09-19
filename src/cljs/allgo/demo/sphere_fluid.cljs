(ns allgo.demo.sphere-fluid
  "Incompressible flow on a sphere, after `allgo.physics.sphere-fluid`.

  The simulation grid *is* latitude and longitude, so the texture wrapped
  round the globe is the state itself rather than a projection of it --
  one texel per cell, no resampling anywhere. That is the opposite of the
  case `allgo.demo.planet` makes, and the difference is worth noticing:
  there the field was a continuous function and the grid was a lossy way
  of looking at it, while here the grid is what exists.

  What to look at:

  * **Rotation.** Turn it up and the vortices stop being vortices and
    become bands. Nothing in the code says \"make bands\"; it falls out of
    advecting `omega + f` instead of `omega`.
  * **Advection.** Switch from corrected to plain semi-Lagrangian and
    watch the whole field dissolve into a smudge over a few seconds. Both
    are stable; only one of them is worth looking at.
  * **Vorticity against tracer.** The vorticity is what is solved for; the
    tracer is a passive dye that shows what the velocity does to material,
    and stretches into filaments far finer than the vortices that make
    them.

  ## Which frame you are watching from

  The simulation is solved in the frame that turns with the sphere -- that
  is what the Coriolis term is for -- so `co-rotating` is not a convenience
  view, it is the frame the numbers are already in, and holding the sphere
  still is what lets the weather be visible at all.

  `inertial` turns the sphere at its real rate, and at any setting that
  produces bands it will strobe. That is the honest answer rather than a
  defect: the planet turns roughly thirty times while material crosses it
  once, which is about Jupiter's ratio, and there is no way to show both
  motions at one speed. `drifting` makes no physical claim -- it is a slow
  turn so the far side comes round without dragging.

  The `speed` control is the other half of this. The jets reach about 1.1
  radians per unit time, which at full speed laps the planet in under five
  seconds, so the default runs at a quarter of that. Nothing about the
  physics changes with it; it only sets how much simulated time a frame
  covers.

  The polar axis is drawn because banded flow looks like stripes from
  everywhere, and stripes say nothing about which way is up."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.sphere-fluid :as sf]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def config
  {:radius 1.0
   ;; About two thirds of a cell per step at the default resolution.
   ;; Semi-Lagrangian advection has no stability limit, so this is chosen
   ;; for how fast the flow reads on screen and for how much interpolation
   ;; error a step is worth -- not by a CFL condition, which is the whole
   ;; point of the method.
   :dt 0.02
   :steps-per-frame 1
   ;; The jets reach about 1.1 radians per unit time, so at full speed
   ;; they lap the planet in under five seconds of wall clock -- which
   ;; reads as a blur rather than as weather. Slowed to a quarter by
   ;; default; the control goes back up.
   :speed 0.25
   ;; Cosmetic drift, for seeing the far side without dragging: one turn
   ;; every seventy seconds.
   :drift 0.0015
   ;; nlon must be a power of two -- the longitude transform is radix two.
   :resolutions {"64 x 32" [32 64] "128 x 64" [64 128] "256 x 128" [128 256]}
   :default-resolution "128 x 64"
   ;; Jupiter, near enough. The jets are put in rather than waited for --
   ;; on a real planet they are maintained by forcing from below over
   ;; geological time, and a decaying simulation started from scattered
   ;; vortices just decays. What the rotation then does to them is the
   ;; thing worth watching.
   :jets 9
   :jet-amplitude 8.0
   ;; Rossby's answer to "how fast should it spin": the jets roll up only
   ;; while their own vorticity gradient beats the planetary one, which
   ;; puts the threshold at amplitude * jets / 2. Just under it and the
   ;; bands survive while their shear lines curl; over it and beta wins,
   ;; the bands go glassy and the whorls disappear entirely.
   :rotation-fraction 0.85
   :seed-vortices 3})

(defn- kuo-threshold
  "The rotation at which the jets stop being unstable, and so stop
  producing vortices."
  [amplitude jets]
  (* 0.5 (double amplitude) (double jets)))

;; ---------------------------------------------------------------------------
;; Colour

(defn- clamp01 [x] (max 0.0 (min 1.0 (double x))))

(defn- ramp
  "Piecewise-linear colour ramp through `stops`, each `[t r g b]`."
  [stops t]
  (let [t (clamp01 t)]
    (loop [[[t0 r0 g0 b0] & more] stops]
      (let [[t1 r1 g1 b1] (first more)]
        (cond
          (nil? t1) [r0 g0 b0]
          (<= t t1) (let [f (if (== t1 t0) 0.0 (/ (- t t0) (- t1 t0)))]
                      [(+ r0 (* f (- r1 r0)))
                       (+ g0 (* f (- g1 g0)))
                       (+ b0 (* f (- b1 b0)))])
          :else (recur more))))))

(def ^:private vorticity-stops
  "Diverging, and dark in the middle: the interesting thing about a
  vorticity field is where it is *not* zero, and a ramp that runs through
  black rather than through white makes the two signs read as two
  separate things rather than as one continuum."
  [[0.0 0.60 0.92 1.00]
   [0.25 0.10 0.45 0.85]
   [0.5 0.04 0.05 0.09]
   [0.75 0.90 0.35 0.10]
   [1.0 1.00 0.93 0.75]])

(def ^:private tracer-stops
  [[0.0 0.03 0.05 0.10]
   [0.4 0.10 0.32 0.45]
   [0.75 0.55 0.78 0.80]
   [1.0 0.98 0.99 0.96]])

;; ---------------------------------------------------------------------------
;; Simulation

(defn- seed!
  "Alternating zonal jets, a few stray vortices to spoil the symmetry, and
  a banded tracer to show what the flow does to material."
  [sim n]
  (let [g (:grid sim)]
    (sf/zonal-jets! g (:vorticity sim)
                    {:n-jets (:jets config)
                     :amplitude (:jet-amplitude config)
                     :perturbation 0.04
                     :wave 6})
    (dotimes [k n]
      ;; Evenly over the sphere rather than uniformly in latitude and
      ;; longitude, which would crowd them at the poles.
      (let [y (- 1.0 (/ (* 2.0 (+ k 0.5)) n))
            r (js/Math.sqrt (max 0.0 (- 1.0 (* y y))))
            th (* 2.39996 k)
            sign (if (even? k) 1.0 -1.0)]
        (sf/add-vortex! g (:vorticity sim)
                        [(* r (js/Math.cos th)) y (* r (js/Math.sin th))]
                        (* sign 1.5)
                        0.2)))
    (when (:tracer sim) (sf/seed-tracer! g (:tracer sim) (* 2 (:jets config))))
    (sf/update-velocity! sim)
    sim))

(defn- build-sim [{:keys [resolution rotation viscosity advection]}]
  (let [[nlat nlon] (get (:resolutions config) resolution
                         (get (:resolutions config) (:default-resolution config)))]
    (seed! (sf/simulation {:grid {:nlat nlat :nlon nlon :radius (:radius config)}
                           :rotation rotation
                           :viscosity viscosity
                           :advection (keyword advection)
                           :tracer? true})
           (:seed-vortices config))))

;; ---------------------------------------------------------------------------
;; Texture

(defn- make-texture [nlon nlat]
  (let [data (js/Uint8Array. (* 4 nlon nlat))
        tex (THREE/DataTexture. data nlon nlat THREE/RGBAFormat THREE/UnsignedByteType)]
    (set! (.-wrapS tex) THREE/RepeatWrapping)
    (set! (.-wrapT tex) THREE/ClampToEdgeWrapping)
    (set! (.-minFilter tex) THREE/LinearFilter)
    (set! (.-magFilter tex) THREE/LinearFilter)
    (set! (.-needsUpdate tex) true)
    {:texture tex :data data}))

(defn- paint!
  "Writes a field into the texture, scaled by its own extreme.

  Rescaling every frame rather than fixing the range: the vorticity decays
  and the tracer mixes, and a fixed scale would have the picture fade to
  nothing while the interesting structure was still there."
  [^js data ^js field nlon nlat stops signed?]
  (let [n (* nlon nlat)
        [mn mx] (loop [i 0 mn 1e30 mx -1e30]
                  (if (= i n)
                    [mn mx]
                    (let [v (aget field i)]
                      (recur (inc i) (min mn v) (max mx v)))))
        ;; A signed field is centred so that zero lands in the middle of
        ;; the ramp, where it is dark; an unsigned one just fills it.
        peak (max 1e-12 (max (js/Math.abs mn) (js/Math.abs mx)))
        lo (if signed? (- peak) mn)
        hi (if signed? peak mx)
        span (max 1e-12 (- hi lo))]
    (dotimes [i n]
      (let [t (/ (- (aget field i) lo) span)
            [r g b] (ramp stops t)
            base (* 4 i)]
        (aset data base (js/Math.round (* 255 (clamp01 r))))
        (aset data (+ base 1) (js/Math.round (* 255 (clamp01 g))))
        (aset data (+ base 2) (js/Math.round (* 255 (clamp01 b))))
        (aset data (+ base 3) 255)))
    data))

(defn- axis-object
  "The rotation axis, drawn so that north is obvious.

  Worth having even though the axis is fixed: once the flow is banded
  every view looks like stripes, and stripes give no clue which way is up
  or which pole you are over."
  []
  (let [group (THREE/Group.)
        pts (js/Float32Array. #js [0.0 -1.55 0.0 0.0 1.55 0.0])
        geo (doto (THREE/BufferGeometry.)
              (.setAttribute "position" (THREE/BufferAttribute. pts 3)))
        line (THREE/Line. geo (THREE/LineBasicMaterial.
                               #js {:color 0x86b5ff :transparent true :opacity 0.7}))
        north (THREE/Mesh. (THREE/ConeGeometry. 0.05 0.14 18)
                           (THREE/MeshBasicMaterial. #js {:color 0xa9d0ff}))
        south (THREE/Mesh. (THREE/SphereGeometry. 0.035 14 10)
                           (THREE/MeshBasicMaterial.
                            #js {:color 0x3f6ea8 :transparent true :opacity 0.8}))]
    (.set (.-position north) 0.0 1.62 0.0)
    (.set (.-position south) 0.0 -1.58 0.0)
    (.add group line)
    (.add group north)
    (.add group south)
    group))

;; ---------------------------------------------------------------------------
;; Scene

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 45 (/ (.-clientWidth container)
                                               (max 1 (.-clientHeight container)))
                                         0.01 100)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        running? (atom false)
        tick-fps! (fps/meter! container)
        group (THREE/Group.)
        state (atom {})
        ^js controls #js {:resolution (:default-resolution config)
                          :rotation (js/Math.round
                                     (* (:rotation-fraction config)
                                        (kuo-threshold (:jet-amplitude config)
                                                       (:jets config))))
                          :viscosity 0.0
                          :advection "corrected"
                          :show "vorticity"
                          :speed (:speed config)
                          :view "co-rotating"
                          :axis true
                          :stir (fn [])
                          :reset (fn [])}]
    (set! (.-background scene) (THREE/Color. 0x05060d))
    (.setSize renderer (.-clientWidth container) (max 1 (.-clientHeight container)))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.add scene group)
    (let [axis (axis-object)]
      (.add scene axis)
      (swap! state assoc :axis axis))
    (.set (.-position camera) 0 0.9 3.0)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (set! (.-minDistance orbit) 1.3)
      (set! (.-maxDistance orbit) 10)
      (letfn [(dispose! []
                (when-let [^js m (:mesh @state)]
                  (.remove group m)
                  (.dispose (.-geometry m))
                  (.dispose (.-map (.-material m)))
                  (.dispose (.-material m))))
              (rebuild! [& _]
                (dispose!)
                (let [sim (build-sim {:resolution (.-resolution controls)
                                      :rotation (double (.-rotation controls))
                                      :viscosity (double (.-viscosity controls))
                                      :advection (.-advection controls)})
                      {:keys [nlat nlon]} (:grid sim)
                      {:keys [texture data]} (make-texture nlon nlat)
                      mesh (THREE/Mesh.
                            (THREE/SphereGeometry. (:radius config) 96 64)
                            (THREE/MeshBasicMaterial. #js {:map texture}))]
                  (.add group mesh)
                  (swap! state merge {:sim sim :mesh mesh :texture texture :data data
                                      :nlat nlat :nlon nlon})))
              (repaint! []
                (let [{:keys [sim ^js data ^js texture nlat nlon]} @state
                      vorticity? (= (.-show controls) "vorticity")
                      field (if vorticity? (:vorticity sim) (:tracer sim))]
                  (paint! data field nlon nlat
                          (if vorticity? vorticity-stops tracer-stops)
                          vorticity?)
                  (set! (.-needsUpdate texture) true)))
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
                  (let [t0 (js/performance.now)
                        dt (* (:dt config) (double (.-speed controls)))]
                    (dotimes [_ (:steps-per-frame config)]
                      (swap! state update :sim sf/step! dt))
                    (repaint!)
                    ;; Three frames to watch this from, and only one of
                    ;; them is a physical statement.
                    ;;
                    ;; The simulation is solved in the *co-rotating*
                    ;; frame -- that is what the Coriolis term is for --
                    ;; so holding the sphere still is not a trick, it is
                    ;; the frame the numbers are already in, and it is the
                    ;; one that lets you see the weather.
                    ;;
                    ;; "inertial" turns the sphere at its actual rate,
                    ;; which at these settings is thirty-four degrees a
                    ;; frame. That it strobes is the honest answer rather
                    ;; than a defect: the planet really does turn about
                    ;; thirty times in the while it takes material to
                    ;; cross it once, which is roughly Jupiter's ratio.
                    ;;
                    ;; "drifting" is neither, and makes no claim to be --
                    ;; a slow turn so the far side comes round.
                    (let [v (.-view controls)]
                      (cond
                        (= v "drifting")
                        (set! (.-y (.-rotation group))
                              (+ (.-y (.-rotation group)) (:drift config)))
                        (= v "inertial")
                        (set! (.-y (.-rotation group))
                              (+ (.-y (.-rotation group))
                                 (* (double (.-rotation controls)) dt)))
                        :else nil))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (set! (.-stir controls)
              (fn []
                (let [{:keys [sim]} @state
                      g (:grid sim)]
                  (dotimes [_ 3]
                    (let [y (- (* 2.0 (js/Math.random)) 1.0)
                          r (js/Math.sqrt (max 0.0 (- 1.0 (* y y))))
                          th (* 2.0 js/Math.PI (js/Math.random))]
                      (sf/add-vortex! g (:vorticity sim)
                                      [(* r (js/Math.cos th)) y (* r (js/Math.sin th))]
                                      (* (if (< (js/Math.random) 0.5) 1.0 -1.0)
                                         (+ 5.0 (* 5.0 (js/Math.random))))
                                      (+ 0.15 (* 0.1 (js/Math.random))))))
                  (sf/update-velocity! sim))))
        (set! (.-reset controls) (fn [] (rebuild!)))
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (repaint!)
        (let [gui (GUI. #js {:container container})
              flow (.addFolder gui "flow")
              view (.addFolder gui "view")]
          (doto flow
            (-> (.add controls "rotation" 0.0 90.0 1.0)
                (.onChange (fn [v] (swap! state assoc-in [:sim :rotation] (double v)))))
            (-> (.add controls "viscosity" 0.0 0.005 0.0001)
                (.onChange (fn [v] (swap! state assoc-in [:sim :viscosity] (double v)))))
            (-> (.add controls "advection" #js ["corrected" "semi-lagrangian"])
                (.onChange (fn [v] (swap! state assoc-in [:sim :advection] (keyword v)))))
            (-> (.add controls "speed" 0.05 2.0 0.05))
            (.add controls "stir")
            (.add controls "reset"))
          (doto view
            (-> (.add controls "show" #js ["vorticity" "tracer"]))
            (-> (.add controls "view" #js ["co-rotating" "drifting" "inertial"]))
            (-> (.add controls "resolution" (clj->js (vec (keys (:resolutions config)))))
                (.onChange rebuild!))
            (-> (.add controls "axis")
                (.onChange (fn [v] (set! (.-visible (:axis @state)) v))))))
        {:start (fn [] (when-not @running?
                         (reset! running? true)
                         (on-resize)
                         (animate)))
         :stop (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start!
  "Build the scene on first selection, then resume the render loop."
  []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "sphere-fluid")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop!
  "Halt the render loop, keeping the context and the framing."
  []
  (when-let [c @controller] ((:stop c))))
