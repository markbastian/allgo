(ns allgo.demo.satellite
  "A satellite under the force model of Montenbruck & Gill chapter 3,
  propagated by the integrators of chapter 4.

  The thing to watch is the orbit plane. Under a spherical Earth an orbit is
  a closed ellipse that never moves; the equatorial bulge torques it, and the
  plane wheels steadily about the pole -- about five degrees a day for an
  ISS-like orbit. Turn the harmonics down to degree zero and the trail
  retraces itself exactly; turn them up and it sweeps out a band."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.forces :as forces]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.time :as atime]
            [allgo.demo.fps :as fps]
            [allgo.numerics :as num]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private scale (/ 20.0 c/R-earth))
(def ^:private max-points 90000)
(def ^:private by-name (into {} (map (juxt :name identity)) num/first-order))

(def ^:private ^js controls
  #js {:altitude        400.0
       :inclination     51.6
       :degree          4
       :sun             false
       :moon            false
       :srp             false
       :drag            false
       :tides           false
       :relativity      false
       :method          "RK4"
       :stepSeconds     20.0
       :secondsPerFrame 300.0
       :running         true})

(defn- config []
  {:degree      (int (.-degree controls))
   :sun?        (.-sun controls)
   :moon?       (.-moon controls)
   :tides?      (.-tides controls)
   :relativity? (.-relativity controls)
   :field       geo/earth
   :srp         (when (.-srp controls) {:area-to-mass 0.02 :cr 1.3})
   :drag        (when (.-drag controls) {:area-to-mass 0.02 :cd 2.2})})

(defn- fresh []
  (let [[r0 v0] (forces/circular-state (+ c/R-earth (.-altitude controls))
                                       (* (.-inclination controls) c/degrees))]
    {:integ (num/integrator (by-name (.-method controls))
                            (forces/first-order (config) c/mjd-J2000)
                            0.0 (into r0 v0) (.-stepSeconds controls)
                            {:adaptive? false})
     :n 0 :raan0 nil}))

(defn- raan
  "Right ascension of the ascending node -- the angle the orbit plane makes
  with the vernal equinox, and the element J2 drives."
  [[x y z] [vx vy vz]]
  (let [hx (- (* y vz) (* z vy))
        hy (- (* z vx) (* x vz))]
    (js/Math.atan2 hx (- hy))))

(defn- graticule
  "A wireframe Earth. A plain sphere would hide its rotation, and the field
  is fixed to the planet, not the sky."
  []
  (THREE/LineSegments.
   (THREE/WireframeGeometry. (THREE/SphereGeometry. (* scale c/R-earth) 24 16))
   (THREE/LineBasicMaterial. #js {:color 0x2a4a7a :transparent true :opacity 0.35})))

(defn init! [^js container]
  (let [scene     (THREE/Scene.)
        camera    (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.1 5000)
        renderer  (THREE/WebGLRenderer. #js {:antialias true})
        geo-buf   (doto (THREE/BufferGeometry.)
                    (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* max-points 3)) 3)))
        trail     (THREE/Line. geo-buf (THREE/LineBasicMaterial. #js {:color 0x6c8cff :transparent true :opacity 0.8}))
        sat       (THREE/Mesh. (THREE/SphereGeometry. 0.35 12 10)
                               (THREE/MeshBasicMaterial. #js {:color 0xffd479}))
        earth     (graticule)
        equator   (THREE/Line. (doto (THREE/BufferGeometry.)
                                 (.setFromPoints (clj->js (for [k (range 65)]
                                                            (let [a (* 2 js/Math.PI (/ k 64))]
                                                              (THREE/Vector3. (* scale c/R-earth (js/Math.cos a))
                                                                              0
                                                                              (* scale c/R-earth (js/Math.sin a))))))))
                               (THREE/LineBasicMaterial. #js {:color 0xff7a5c :transparent true :opacity 0.6}))
        readout   (js/document.createElement "div")
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom (fresh))]
    (set! (.-className readout) "numeric-readout")
    (.appendChild container readout)
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (set! (.-frustumCulled trail) false)
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (doseq [o [trail sat earth equator]] (.add scene o))
    (.set (.-position camera) 46 30 46)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(restart! []
                (reset! state (fresh))
                (.setDrawRange geo-buf 0 0))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (advance! []
                (swap! state
                       (fn [{:keys [integ n raan0] :as st}]
                         (let [target (+ (:t integ) (.-secondsPerFrame controls))
                               s'     (num/step-until integ target 4000)
                               y      (:y s')
                               r      (subvec y 0 3)
                               ;; three.js is y-up; the orbital frame is z-up
                               [x yy z] r
                               i      (if (>= n max-points) 0 n)
                               arr    (.-array (.getAttribute geo-buf "position"))]
                           (aset arr (* i 3) (* scale x))
                           (aset arr (+ (* i 3) 1) (* scale z))
                           (aset arr (+ (* i 3) 2) (* scale yy))
                           (assoc st :integ s' :n (inc i)
                                  :raan0 (or raan0 (raan r (subvec y 3 6))))))))
              (publish! []
                (let [{:keys [integ n raan0]} @state
                      y   (:y integ)
                      r   (subvec y 0 3)
                      v   (subvec y 3 6)
                      [x yy z] r
                      days (/ (:t integ) 86400.0)
                      alt (- (js/Math.sqrt (reduce + (map * r r))) c/R-earth)
                      dr  (let [d (- (raan r v) raan0)]
                            ;; unwrap: the node can pass through +/- pi
                            (js/Math.atan2 (js/Math.sin d) (js/Math.cos d)))]
                  (.setDrawRange geo-buf 0 n)
                  (set! (.-needsUpdate (.getAttribute geo-buf "position")) true)
                  (.set (.-position sat) (* scale x) (* scale z) (* scale yy))
                  ;; the field is fixed to the Earth, so show it turning
                  (set! (.-y (.-rotation earth)) (- (atime/gmst (+ c/mjd-J2000 days))))
                  (set! (.-textContent readout)
                        (str (.-method controls) "  deg " (int (.-degree controls))
                             "   t=" (.toFixed days 3) " d"
                             "   alt=" (.toFixed alt 1) " km"
                             "\\nnode drift " (.toFixed (/ (* dr 180.0 (if (> days 0.05) (/ 1.0 days) 0.0)) js/Math.PI) 3)
                             " deg/day   (J2 predicts "
                             (.toFixed (let [a (+ c/R-earth (.-altitude controls))
                                             q (/ c/R-earth a)
                                             nm (js/Math.sqrt (/ c/GM-earth (* a a a)))]
                                         (* -1.5 geo/J2 q q nm (js/Math.cos (* (.-inclination controls) c/degrees))
                                            86400.0 (/ 180.0 js/Math.PI)))
                                       3)
                             ")"))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (when (.-running controls) (advance!))
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "altitude" 300 2000 10) (.onChange restart!))
          (-> (.add gui controls "inclination" 0 180 0.1) (.onChange restart!))
          (-> (.add gui controls "degree" 0 8 1) (.onChange restart!))
          (.add gui controls "secondsPerFrame" 30 3000 10)
          (.add gui controls "running")
          (let [f (.addFolder gui "forces")]
            (doseq [k ["sun" "moon" "srp" "drag" "tides" "relativity"]]
              (-> (.add f controls k) (.onChange restart!))))
          (let [f (.addFolder gui "integrator")]
            (-> (.add f controls "method" (clj->js (mapv :name num/first-order))) (.onChange restart!))
            (-> (.add f controls "stepSeconds" 1 120 1) (.onChange restart!))
            (.close f))
          (.add gui #js {:restart restart!} "restart"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "satellite")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
