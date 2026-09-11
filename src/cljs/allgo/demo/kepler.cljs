(ns allgo.demo.kepler
  "Two-body motion integrated by the Nystrom family in real time.

  This is the problem Montenbruck & Gill's chapter is written for, and it is
  the one where the choice of method shows up as something other than
  accuracy. The strip chart tracks relative energy error, because that is
  where the interesting distinction lives: Verlet is only second order and
  starts far behind RKN4, but its error oscillates within a bound while
  RKN4's walks steadily upward. Run either for a few hundred orbits and the
  more accurate method is the one that has drifted further."
  (:require [allgo.demo.fps :as fps]
            [allgo.numerics :as num]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private mu 1.0)
(def ^:private max-points 60000)
(def ^:private view-scale 34.0)
(def ^:private by-name (into {} (map (juxt :name identity)) num/second-order))

(def ^:private ^js controls
  #js {:method         "RKN4"
       :eccentricity   0.6
       :stepsPerOrbit  300
       :stepsPerFrame  8
       :trail          true
       :running        true})

(defn- period [] (num/orbital-period mu 1.0))

(defn- fresh []
  (let [[r0 v0] (num/periapsis-state mu 1.0 (.-eccentricity controls))]
    {:integ  (num/integrator-2 (by-name (.-method controls)) (num/kepler mu)
                               0.0 r0 v0 (/ (period) (.-stepsPerOrbit controls)))
     :energy (num/specific-energy mu r0 v0)
     :n      0
     :steps  0
     :chart  []}))

;; ------------------------------------------------------------- strip chart

(def ^:private chart-w 210)
(def ^:private chart-h 74)
(def ^:private chart-lo -14.0)
(def ^:private chart-hi -1.0)

(defn- draw-chart! [^js ctx samples]
  (.clearRect ctx 0 0 chart-w chart-h)
  (set! (.-fillStyle ctx) "rgba(5,7,13,0.72)")
  (.fillRect ctx 0 0 chart-w chart-h)
  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.10)")
  (set! (.-lineWidth ctx) 1)
  (doseq [decade (range (Math/ceil chart-lo) chart-hi 3)]
    (let [y (* chart-h (- 1.0 (/ (- decade chart-lo) (- chart-hi chart-lo))))]
      (.beginPath ctx) (.moveTo ctx 0 y) (.lineTo ctx chart-w y) (.stroke ctx)))
  (when (> (count samples) 1)
    (set! (.-strokeStyle ctx) "#6c8cff")
    (set! (.-lineWidth ctx) 1.5)
    (.beginPath ctx)
    (doseq [[i v] (map-indexed vector samples)]
      (let [x (* chart-w (/ i (max 1 (dec (count samples)))))
            y (* chart-h (- 1.0 (/ (- v chart-lo) (- chart-hi chart-lo))))]
        (if (zero? i) (.moveTo ctx x y) (.lineTo ctx x y))))
    (.stroke ctx))
  (set! (.-fillStyle ctx) "rgba(148,148,171,0.9)")
  (set! (.-font ctx) "9px ui-monospace, Menlo, monospace")
  (.fillText ctx "log |dE/E|" 6 11))

(defn- trail-geometry []
  (doto (THREE/BufferGeometry.)
    (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* max-points 3)) 3))))

(defn init! [^js container]
  (let [scene     (THREE/Scene.)
        camera    (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.1 5000)
        renderer  (THREE/WebGLRenderer. #js {:antialias true})
        geo       (trail-geometry)
        trail     (THREE/Line. geo (THREE/LineBasicMaterial. #js {:color 0x6c8cff :transparent true :opacity 0.75}))
        body      (THREE/Mesh. (THREE/SphereGeometry. 1.4 20 14)
                               (THREE/MeshBasicMaterial. #js {:color 0xffd479}))
        primary   (THREE/Mesh. (THREE/SphereGeometry. 3.4 28 20)
                               (THREE/MeshBasicMaterial. #js {:color 0xffb03a}))
        readout   (js/document.createElement "div")
        chart     (js/document.createElement "canvas")
        ctx       (.getContext chart "2d")
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom (fresh))]
    (set! (.-className readout) "numeric-readout")
    (set! (.-className chart) "numeric-chart")
    (set! (.-width chart) chart-w)
    (set! (.-height chart) chart-h)
    (.appendChild container readout)
    (.appendChild container chart)
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (set! (.-frustumCulled trail) false)
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (doseq [o [trail body primary]] (.add scene o))
    (.set (.-position camera) 0 70 105)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(restart! []
                (reset! state (fresh))
                (.setDrawRange geo 0 0))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (record! [{:keys [integ energy n chart] :as st}]
                (let [rel (abs (/ (- (num/specific-energy mu (:y integ) (:dy integ)) energy) energy))
                      arr (.-array (.getAttribute geo "position"))
                      ;; the trail overlaps itself every orbit, so precession
                      ;; from integration error is visible as the ellipse
                      ;; failing to retrace its own path
                      i   (if (>= n max-points) 0 n)
                      [x y z] (:y integ)]
                  (aset arr (* i 3) (* view-scale x))
                  (aset arr (+ (* i 3) 1) (* view-scale z))
                  (aset arr (+ (* i 3) 2) (* view-scale y))
                  (assoc st :n (inc i) :steps (inc (:steps st 0))
                         :chart (let [v (Math/log10 (max 1e-16 rel))]
                                  (if (>= (count chart) chart-w)
                                    (conj (subvec chart 1) v)
                                    (conj chart v)))
                         :rel rel)))
              (integrate! []
                (swap! state
                       (fn [st]
                         (reduce (fn [s _] (record! (update s :integ num/step)))
                                 st
                                 (range (.-stepsPerFrame controls))))))
              (publish! []
                (let [{:keys [integ n rel chart steps]} @state
                      [x y z] (:y integ)]
                  (set! (.-visible trail) (.-trail controls))
                  (.setDrawRange geo 0 n)
                  (set! (.-needsUpdate (.getAttribute geo "position")) true)
                  (.set (.-position body) (* view-scale x) (* view-scale z) (* view-scale y))
                  (draw-chart! ctx chart)
                  (set! (.-textContent readout)
                        ;; Evaluations, not steps: GBS8-2 holds its energy far
                        ;; better than Verlet but spends ten times the force
                        ;; calls doing it, and the chart alone hides that.
                        (let [stages (:stages (:method integ))]
                          (str (.-method controls)
                               "   orbits=" (.toFixed (/ (:t integ) (period)) 2)
                               "   |dE/E|=" (.toExponential (or rel 0) 2)
                               "\n" stages " f-evals/step"
                               "   " (.toLocaleString (* stages steps)) " total")))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (when (.-running controls) (integrate!))
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "method" (clj->js (mapv :name num/second-order)))
              (.onChange restart!))
          (-> (.add gui controls "eccentricity" 0.0 0.9 0.01) (.onChange restart!))
          (-> (.add gui controls "stepsPerOrbit" 30 2000 10) (.onChange restart!))
          (.add gui controls "stepsPerFrame" 1 60 1)
          (.add gui controls "trail")
          (.add gui controls "running")
          (.add gui #js {:restart restart!} "restart"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "kepler")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
