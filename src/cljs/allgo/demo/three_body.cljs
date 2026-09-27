(ns allgo.demo.three-body
  "Three bodies under their mutual gravity, integrated in real time.

  The special solutions -- the figure eight, Lagrange's triangle, the
  butterfly -- repeat forever on paper and only as long as the integrator
  keeps them there. The Pythagorean problem and the random free falls are
  the general case: close encounters that need steps a thousand times
  shorter than the rest, and chaos, which the ghost shows. It is the same
  system started a hair to one side, and the distance between the two
  grows until it is as large as the orbits themselves."
  (:require [allgo.astro.three-body :as tb]
            [allgo.demo.error-chart :as chart]
            [allgo.demo.fps :as fps]
            [allgo.geometry.vec3 :as v3]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private max-points 4000)
(def ^:private frame-view 30.0)
(def ^:private step-budget 4000)
(def ^:private colors [0xff6b6b 0x6c8cff 0xffd479])
(def ^:private by-name (into {} (map (juxt :name identity)) tb/scenarios))

(def ^:private ^js controls
  #js {:scenario  "Figure eight"
       :seed      1
       :method    "DOPRI5(4)"
       :tolerance -10
       :speed     1.0
       :trail     1500
       :ghost     false
       :nudge     -6
       :running   true})

(defn- scenario []
  (let [sc (by-name (.-scenario controls))]
    (if (= :free-fall (:id sc))
      (merge sc (tb/free-fall (.-seed controls)))
      sc)))

(defn- build [sc]
  (tb/integrator sc (.-method controls) (js/Math.pow 10 (.-tolerance controls)) 0.02))

(defn- fresh []
  (let [{:keys [masses pos vel] :as sc} (scenario)]
    {:scenario sc
     :integ    (build sc)
     :ghost    (build (tb/nudge sc (js/Math.pow 10 (.-nudge controls))))
     :energy   (tb/energy masses pos vel)
     :steps    0
     :chart    []}))

;; ------------------------------------------------------------------ trails

(defn- trail [color]
  (let [geo (doto (THREE/BufferGeometry.)
              (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* max-points 3)) 3)))
        line (THREE/Line. geo (THREE/LineBasicMaterial. #js {:color color :transparent true :opacity 0.8}))]
    (set! (.-frustumCulled line) false)
    (.setDrawRange geo 0 0)
    line))

(defn- extend-trail!
  "Adds `points` (scene coordinates, flat) to the end of a trail, keeping at
  most `keep` of the newest. The oldest are dropped by sliding the buffer
  down once a frame rather than once a point, so a frame costs one copy of
  the buffer however many steps it took."
  [^js line ^js points keep]
  (let [geo   (.-geometry line)
        attr  (.getAttribute geo "position")
        ^js arr (.-array attr)
        n     (.. geo -drawRange -count)
        k     (/ (.-length points) 3)
        keep  (min max-points keep)
        k     (min k keep)
        points (.subarray points (- (.-length points) (* 3 k)))
        drop  (max 0 (- (+ n k) keep))
        n'    (- n drop)]
    (when (pos? drop) (.copyWithin arr 0 (* 3 drop) (* 3 n)))
    (.set arr points (* 3 n'))
    (.setDrawRange geo 0 (+ n' k))
    (set! (.-needsUpdate attr) true)))

(defn- clear-trail! [^js line] (.setDrawRange (.-geometry line) 0 0))

;; ------------------------------------------------------------------- scene

(defn- to-scene
  "The orbital plane is the scene's floor: x across, y into the screen, so
  seen from above the motion turns the way it does on paper."
  [scale [x y z]]
  [(* scale x) (* scale z) (* scale (- y))])

(defn init! [^js container]
  (let [scene     (THREE/Scene.)
        camera    (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.1 5000)
        renderer  (THREE/WebGLRenderer. #js {:antialias true})
        trails    (mapv trail colors)
        spheres   (mapv #(THREE/Mesh. (THREE/SphereGeometry. 1.0 24 16)
                                      (THREE/MeshBasicMaterial. #js {:color %}))
                        colors)
        ghosts    (mapv #(THREE/Mesh. (THREE/SphereGeometry. 1.0 16 10)
                                      (THREE/MeshBasicMaterial. #js {:color % :wireframe true
                                                                     :transparent true :opacity 0.55}))
                        colors)
        grid      (THREE/GridHelper. 120 24 0x1c2233 0x121726)
        readout   (js/document.createElement "div")
        ctx       (chart/mount! container)
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom (fresh))]
    (set! (.-className readout) "numeric-readout")
    (.appendChild container readout)
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (.appendChild container (.-domElement renderer))
    (doseq [o (concat [grid] trails spheres ghosts)] (.add scene o))
    (.set (.-position camera) 0 62 42)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(scale [] (/ frame-view (:view (:scenario @state))))
              (restart! []
                (reset! state (fresh))
                (run! clear-trail! trails)
                (doseq [[^js s m] (map vector spheres (:masses (:scenario @state)))]
                  ;; volume in proportion to mass
                  (.setScalar (.-scale s) (* 1.1 (js/Math.cbrt m))))
                (doseq [[^js g ^js s] (map vector ghosts spheres)]
                  (.copy (.-scale g) (.-scale s))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (integrate! []
                (let [{:keys [integ ghost scenario energy]} @state
                      t-end   (+ (:t integ) (/ (.-speed controls) 60.0))
                      [integ path] (tb/advance integ t-end step-budget)
                      ;; the ghost is held to the main run's clock, so the
                      ;; two are compared at the same moment even when the
                      ;; step budget ran out before `t-end`
                      [ghost _] (tb/advance ghost (:t integ) step-budget)
                      s       (scale)
                      {:keys [masses]} scenario
                      rel     (abs (/ (- (tb/energy masses (tb/positions integ) (tb/velocities integ))
                                         energy)
                                      energy))]
                  (doseq [b (range 3)]
                    (extend-trail! (trails b)
                                   (js/Float32Array.from
                                    (clj->js (mapcat #(to-scene s (tb/body % b)) path)))
                                   (.-trail controls)))
                  (swap! state assoc :integ integ :ghost ghost :rel rel
                         :steps (+ (:steps @state) (count path))
                         :chart (chart/push (:chart @state) rel))))
              (publish! []
                (let [{:keys [integ ghost rel steps chart scenario]} @state
                      s     (scale)
                      ps    (tb/bodies (tb/positions integ))
                      gs    (tb/bodies (tb/positions ghost))
                      apart (reduce max (map v3/distance ps gs))
                      {:keys [period]} scenario]
                  (doseq [[^js o p] (map vector spheres ps)]
                    (let [[x y z] (to-scene s p)] (.set (.-position o) x y z)))
                  (doseq [[^js o p] (map vector ghosts gs)]
                    (set! (.-visible o) (.-ghost controls))
                    (let [[x y z] (to-scene s p)] (.set (.-position o) x y z)))
                  (doseq [^js t trails] (set! (.-visible t) (pos? (.-trail controls))))
                  (chart/draw! ctx chart)
                  (set! (.-textContent readout)
                        (str (.-method controls)
                             "   t=" (.toFixed (:t integ) 2)
                             (when period (str " (" (.toFixed (/ (:t integ) period) 2) " periods)"))
                             "   |dE/E|=" (.toExponential (or rel 0) 1)
                             "\n" (.toLocaleString steps) " steps, h=" (.toExponential (:h integ) 1)
                             (when (.-ghost controls)
                               (str "   ghost apart " (.toExponential apart 1)))))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (when (.-running controls) (integrate!))
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (restart!)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "scenario" (clj->js (mapv :name tb/scenarios))) (.onChange restart!))
          (-> (.add gui controls "seed" 1 100 1) (.name "seed (free fall)") (.onChange restart!))
          (-> (.add gui controls "method" (clj->js (mapv :name tb/integrators))) (.onChange restart!))
          (-> (.add gui controls "tolerance" -13 -5 1) (.name "tolerance (log10)") (.onChange restart!))
          (.add gui controls "speed" 0.1 20.0 0.1)
          (.add gui controls "trail" 0 max-points 50)
          (.add gui controls "ghost")
          (-> (.add gui controls "nudge" -12 -2 1) (.name "nudge (log10)") (.onChange restart!))
          (.add gui controls "running")
          (.add gui #js {:restart restart!} "restart"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "three-body")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
