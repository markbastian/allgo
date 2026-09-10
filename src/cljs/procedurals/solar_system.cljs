(ns procedurals.solar-system
  "The solar system, the Moon and a satellite drawn in one frame.

  Everything here comes from `procedurals.astro` already referred to the
  mean equator and equinox of J2000, so no body needs rotating to sit
  beside another. That is the point of the demo as much as the picture is:
  the planets come from Vallado's Keplerian elements, the Moon from an
  entirely separate analytical series, and they can simply be added.

  Scale is the real difficulty. Neptune is 30 AU out and the Moon is 0.0026
  AU from Earth -- a ratio of twelve thousand to one, so any single linear
  scale showing both puts one of them under a pixel. Two honest answers are
  offered rather than one dishonest one: a logarithmic radial scale, which
  keeps every orbit on screen and every angle true while making distances
  unreadable, and a linear scale to zoom around in. The Moon carries its own
  magnification, shown as a number, because at true scale it is invisible
  and pretending otherwise would be the one thing worth avoiding.

  The ecliptic sits tilted 23.4 degrees here, which looks wrong until you
  remember the frame is equatorial: this is the sky as the Earth's equator
  sees it, not as the planets would draw it."
  (:require [procedurals.astro.constants :as c]
            [procedurals.astro.ephemeris :as eph]
            [procedurals.astro.kepler :as kep]
            [procedurals.astro.planets :as pl]
            [procedurals.astro.time :as atime]
            [procedurals.fps :as fps]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private view-radius 26.0)
(def ^:private outer-au 31.0)

(def ^:private palette
  {:mercury 0xb0a08c :venus 0xe8c47a :earth 0x5b9bd5 :mars 0xc1553a
   :jupiter 0xd8a06a :saturn 0xe3d08a :uranus 0x8fd0d8 :neptune 0x5570d0})

(def ^:private dot-size
  ;; Not to scale -- at true scale every planet is far under a pixel. Ordered
  ;; by actual radius so the giants still read as giants.
  {:mercury 0.16 :venus 0.24 :earth 0.26 :mars 0.19
   :jupiter 0.62 :saturn 0.54 :uranus 0.38 :neptune 0.37})

(def ^:private controls
  #js {:scale         "logarithmic"
       :daysPerSecond 12.0
       :moonZoom      400.0
       :showOrbits    true
       :showEcliptic  true
       :running       true})

(defn- radial
  "Map a distance in AU onto the view. Logarithmic keeps all eight orbits on
  screen at once; linear is true and needs zooming."
  [au]
  (* view-radius
     (if (= "logarithmic" (.-scale controls))
       (/ (js/Math.log (+ 1.0 (* 3.0 au))) (js/Math.log (+ 1.0 (* 3.0 outer-au))))
       (/ au outer-au))))

(defn- place
  "A heliocentric position in km, mapped into view coordinates. Direction is
  preserved exactly; only the radius is rescaled, so angles stay true
  whichever mode is chosen."
  [[x y z]]
  (let [r  (js/Math.sqrt (+ (* x x) (* y y) (* z z)))]
    (if (< r 1e-9)
      [0.0 0.0 0.0]
      (let [s (/ (radial (/ r c/AU)) r)]
        ;; three.js is y-up; the frame is z-up
        [(* s x) (* s z) (* s y)]))))

(defn- orbit-points
  "One full revolution of a planet, sampled in mean anomaly rather than in
  time, so a slow outer planet costs no more samples than a fast inner one."
  [planet mjd n]
  (let [el (pl/elements-at planet mjd)
        ce (js/Math.cos eph/obliquity-J2000)
        se (js/Math.sin eph/obliquity-J2000)]
    (for [k (range (inc n))]
      (let [M  (* 2.0 js/Math.PI (/ k n))
            nu (kep/mean->true M (:e el))
            [r _] (kep/elements->state 1.0 (assoc el :nu nu :M nil))
            [x y z] r]
        (place [x (- (* ce y) (* se z)) (+ (* se y) (* ce z))])))))

(defn- line-of [points colour opacity]
  (THREE/Line.
   (doto (THREE/BufferGeometry.)
     (.setFromPoints (clj->js (map (fn [[x y z]] (THREE/Vector3. x y z)) points))))
   (THREE/LineBasicMaterial. #js {:color colour :transparent true :opacity opacity})))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.05 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        sun      (THREE/Mesh. (THREE/SphereGeometry. 0.85 24 16)
                              (THREE/MeshBasicMaterial. #js {:color 0xffcf5c}))
        bodies   (into {} (for [k pl/order]
                            [k (THREE/Mesh. (THREE/SphereGeometry. (dot-size k) 16 12)
                                            (THREE/MeshBasicMaterial. #js {:color (palette k)}))]))
        moon     (THREE/Mesh. (THREE/SphereGeometry. 0.12 12 10)
                              (THREE/MeshBasicMaterial. #js {:color 0xcfd6e0}))
        moon-line (THREE/Line. (doto (THREE/BufferGeometry.)
                                 (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* 130 3)) 3)))
                               (THREE/LineBasicMaterial. #js {:color 0x8fa0b8 :transparent true :opacity 0.6}))
        orbits   (THREE/Group.)
        ecliptic (THREE/Group.)
        readout  (js/document.createElement "div")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {:mjd c/mjd-J2000 :last nil})]
    (set! (.-className readout) "numeric-readout")
    (.appendChild container readout)
    (set! (.-background scene) (THREE/Color. 0x04050a))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (.add scene sun)
    (doseq [[_ m] bodies] (.add scene m))
    (.add scene moon)
    (.add scene moon-line)
    (.add scene orbits)
    (.add scene ecliptic)
    (.set (.-position camera) 0 34 34)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(rebuild-orbits! []
                (doseq [^js child (vec (.-children orbits))] (.remove orbits child))
                (doseq [k pl/order]
                  (.add orbits (line-of (orbit-points k (:mjd @state) 128) (palette k) 0.45)))
                (doseq [^js child (vec (.-children ecliptic))] (.remove ecliptic child))
                ;; A ring in the ecliptic plane, tilted into the equatorial
                ;; frame -- the visible reminder of which frame this is.
                (let [ce (js/Math.cos eph/obliquity-J2000)
                      se (js/Math.sin eph/obliquity-J2000)
                      r  (radial 31.0)]
                  (.add ecliptic
                        (line-of (for [k (range 129)]
                                   (let [a (* 2.0 js/Math.PI (/ k 128.0))
                                         x (* r (js/Math.cos a))
                                         y (* r (js/Math.sin a))]
                                     [x (* se y) (* ce y)]))
                                 0x445070 0.5))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (update-bodies! []
                (let [mjd (:mjd @state)
                      earth (pl/heliocentric :earth mjd)]
                  (doseq [k pl/order]
                    (let [[x y z] (place (pl/heliocentric k mjd))
                          ^js m (bodies k)]
                      (.set (.-position m) x y z)))
                  ;; The Moon is geocentric already; adding it to the Earth's
                  ;; heliocentric place is legitimate only because both are in
                  ;; the same frame, which is the whole point.
                  (let [[ex ey ez] (place earth)
                        [mx my mz] (eph/moon mjd)
                        k (* (.-moonZoom controls) (/ view-radius outer-au) (/ 1.0 c/AU))]
                    (.set (.-position moon) (+ ex (* k mx)) (+ ey (* k mz)) (+ ez (* k my)))
                    ;; and its orbit, drawn at the same magnification
                    (let [arr (.-array (.getAttribute (.-geometry moon-line) "position"))]
                      (doseq [i (range 130)]
                        (let [[ox oy oz] (eph/moon (+ mjd (* 27.32 (/ i 129.0))))]
                          (aset arr (* i 3) (+ ex (* k ox)))
                          (aset arr (+ (* i 3) 1) (+ ey (* k oz)))
                          (aset arr (+ (* i 3) 2) (+ ez (* k oy)))))
                      (.setDrawRange (.-geometry moon-line) 0 130)
                      (set! (.-needsUpdate (.getAttribute (.-geometry moon-line) "position")) true)))))
              (publish! []
                (let [mjd (:mjd @state)
                      [yr mo dy hr] (atime/mjd->calendar mjd)
                      moon-km (let [[a b cc] (eph/moon mjd)]
                                (js/Math.sqrt (+ (* a a) (* b b) (* cc cc))))]
                  (set! (.-visible orbits) (.-showOrbits controls))
                  (set! (.-visible ecliptic) (.-showEcliptic controls))
                  (set! (.-textContent readout)
                        (str (.toFixed yr 0) "-" (.padStart (str mo) 2 "0") "-" (.padStart (str dy) 2 "0")
                             "  " (.toFixed hr 0) "h   MJD " (.toFixed mjd 1)
                             "\\n" (.-scale controls) " radial scale"
                             "   Moon at " (.toFixed moon-km 0) " km, drawn "
                             (.toFixed (.-moonZoom controls) 0) "x"))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)
                        {:keys [last]} @state]
                    (when (and (.-running controls) last)
                      (swap! state update :mjd + (* (.-daysPerSecond controls) (/ (- t0 last) 1000.0))))
                    (swap! state assoc :last t0)
                    (update-bodies!)
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild-orbits!)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "scale" #js ["logarithmic" "linear"]) (.onChange rebuild-orbits!))
          (.add gui controls "daysPerSecond" 0 200 1)
          (.add gui controls "moonZoom" 1 2000 1)
          (.add gui controls "showOrbits")
          (.add gui controls "showEcliptic")
          (.add gui controls "running")
          (.add gui #js {:toJ2000 (fn [] (swap! state assoc :mjd c/mjd-J2000) (rebuild-orbits!))} "toJ2000"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "solar-system")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
