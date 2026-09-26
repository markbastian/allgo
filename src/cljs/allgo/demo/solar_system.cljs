(ns allgo.demo.solar-system
  "The solar system, the Moon and a satellite drawn in one frame.

  Everything here comes from `allgo.astro` already referred to the
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
  sees it, not as the planets would draw it.

  Behind it all is the naked-eye sky, the nine thousand stars of the
  Bright Star Catalogue (`allgo.astro.stars`), in the same frame and so
  needing no turning either: the planets pass in front of the zodiac's
  constellations because that is where they are. The stars sit on a
  sphere centered on the camera, which is what being infinitely far away
  looks like -- orbit the view and they turn with it, zoom and they do
  not move. The catalogue is fetched when the demo starts rather than
  built into the page, being half a megabyte the other demos have no use
  for.

  The constellations' figures and names, and the IAU's names for the
  brightest stars, can be drawn over them. The figures come from
  d3-celestial rather than from the catalogue, so a figure landing on its
  stars is two independent sources agreeing about where the stars are --
  and a figure the right way round, Orion with Betelgeuse at his upper
  left, is the check that the sky is not drawn in a mirror. North is up,
  as a star chart has it, until the view is turned."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.frames :as frames]
            [allgo.astro.kepler :as kep]
            [allgo.astro.planets :as pl]
            [allgo.astro.stars :as stars]
            [allgo.astro.time :as atime]
            [allgo.demo.fps :as fps]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]
            ["three/examples/jsm/renderers/CSS2DRenderer.js" :refer [CSS2DObject CSS2DRenderer]]))

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

(def ^:private ^js controls
  #js {:scale         "logarithmic"
       :daysPerSecond 12.0
       :moonZoom      400.0
       :showOrbits    true
       :showEcliptic  true
       :showStars     true
       :starLimit     6.5
       :showLines     false
       :showConNames  false
       :showStarNames false
       :running       true})

(def ^:private to-view
  "EME2000 into three.js's y-up axes. See `allgo.astro.frames/y-up`, and
  its test: this was once a swap of y and z, which is a mirror."
  frames/y-up)

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
        (to-view [(* s x) (* s y) (* s z)])))))

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

(defn- line-of [points color opacity]
  (THREE/Line.
   (doto (THREE/BufferGeometry.)
     (.setFromPoints (clj->js (map (fn [[x y z]] (THREE/Vector3. x y z)) points))))
   (THREE/LineBasicMaterial. #js {:color color :transparent true :opacity opacity})))

(def ^:private sky-radius
  "How far out the stars are drawn: well beyond the planets, well inside
  the camera's far plane. It does not matter which, since they move with
  the camera."
  1000.0)

(def ^:private ^js faintest
  "The shader's magnitude limit, kept so the slider can reach it."
  #js {:value 6.5})

(defn- sky
  "The catalogue's stars as one cloud of points: each placed at its J2000
  direction, sized and dimmed by its magnitude and tinted by its color.

  Every star is drawn at least a pixel and a bit across, so size alone
  cannot say how faint the faint ones are; brightness does the rest,
  falling with the square root of the flux so that a sixth-magnitude star
  is dim rather than gone. Stars fainter than the slider's limit are
  discarded in the shader, so moving it costs nothing."
  [catalogue]
  (let [n (count catalogue)
        pos (js/Float32Array. (* 3 n))
        col (js/Float32Array. (* 3 n))
        size (js/Float32Array. n)
        mag (js/Float32Array. n)]
    (doseq [[i s] (map-indexed vector catalogue)]
      (let [[x y z] (to-view (stars/direction s))
            v (double (:vmag s))
            [r g b] (stars/color (:bv s))
            glow (min 1.0 (+ 0.4 (* 0.6 (js/Math.sqrt (/ (stars/flux v) (stars/flux 1.0))))))]
        (aset pos (* 3 i) (* sky-radius x))
        (aset pos (+ 1 (* 3 i)) (* sky-radius y))
        (aset pos (+ 2 (* 3 i)) (* sky-radius z))
        (aset col (* 3 i) (* glow r))
        (aset col (+ 1 (* 3 i)) (* glow g))
        (aset col (+ 2 (* 3 i)) (* glow b))
        (aset size i (max 2.4 (* 9.0 (js/Math.pow 10.0 (* -0.2 (+ v 1.46))))))
        (aset mag i v)))
    (doto (THREE/Points.
           (doto (THREE/BufferGeometry.)
             (.setAttribute "position" (THREE/BufferAttribute. pos 3))
             (.setAttribute "color" (THREE/BufferAttribute. col 3))
             (.setAttribute "size" (THREE/BufferAttribute. size 1))
             (.setAttribute "vmag" (THREE/BufferAttribute. mag 1)))
           (THREE/ShaderMaterial.
            #js {:uniforms #js {:limit faintest
                                :scale #js {:value (or js/window.devicePixelRatio 1)}}
                 :vertexColors true
                 :transparent true
                 :depthWrite false
                 :blending THREE/AdditiveBlending
                 :vertexShader "
                   uniform float limit;
                   uniform float scale;
                   attribute float size;
                   attribute float vmag;
                   varying vec3 vColor;
                   varying float vHidden;
                   void main() {
                     vColor = color;
                     vHidden = vmag > limit ? 1.0 : 0.0;
                     gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
                     gl_PointSize = size * scale;
                   }"
                 :fragmentShader "
                   varying vec3 vColor;
                   varying float vHidden;
                   void main() {
                     if (vHidden > 0.5) discard;
                     float d = length(gl_PointCoord - 0.5);
                     if (d > 0.5) discard;
                     float a = 1.0 - smoothstep(0.25, 0.5, d);
                     gl_FragColor = vec4(vColor * a, a);
                   }"}))
      ;; Always drawn, and first: it is everything's background.
      (-> .-frustumCulled (set! false))
      (-> .-renderOrder (set! -1)))))

(defn- figures
  "The constellations' figures as one set of faint line segments on the
  sky sphere, a little inside the stars so a line never covers one."
  [constellations]
  (let [segs (for [{:keys [lines]} constellations
                   pl lines
                   [a b] (partition 2 1 pl)
                   p [a b]]
               (let [[x y z] (to-view (apply stars/unit p))
                     r (* 0.995 sky-radius)]
                 (THREE/Vector3. (* r x) (* r y) (* r z))))]
    (doto (THREE/LineSegments.
           (doto (THREE/BufferGeometry.) (.setFromPoints (clj->js segs)))
           (THREE/LineBasicMaterial. #js {:color 0x5d7bb0 :transparent true :opacity 0.45
                                          :depthWrite false}))
      (-> .-frustumCulled (set! false))
      (-> .-renderOrder (set! -1)))))

(defn- label
  "An HTML label pinned to the sky at `[ra dec]`: three.js places it every
  frame, and hides it when it is behind the camera."
  [text cls [ra dec]]
  (let [el (doto (js/document.createElement "div")
             (-> .-className (set! cls))
             ;; The text in a span of its own: three.js owns the div's
             ;; transform, to place it, so any offset has to go inside.
             (.appendChild (doto (js/document.createElement "span")
                             (-> .-textContent (set! text)))))
        [x y z] (to-view (stars/unit ra dec))
        ^js o (CSS2DObject. el)]
    (.set (.-position o) (* sky-radius x) (* sky-radius y) (* sky-radius z))
    o))

(def ^:private named-brighter-than
  "The faintest star given its name: the IAU has named 333 of the
  catalogue's stars, and all of them at once would bury the sky in
  text. The ninety or so brighter than this are the ones a person knows."
  2.5)

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
        state    (atom {:mjd c/mjd-J2000 :last nil})
        stars-at (atom nil)
        ;; Everything infinitely far away, which moves with the camera.
        sky-group (THREE/Group.)
        ;; The figures, and the two kinds of label, once loaded.
        lines-at (atom nil)
        con-labels (atom [])
        star-labels (atom [])
        ^js labels (doto (CSS2DRenderer.)
                     (-> .-domElement .-className (set! "sky-labels")))]
    (.add scene sky-group)
    (.appendChild container (.-domElement labels))
    ;; The catalogue, the figures and the names, fetched once; the demo
    ;; runs without them until they arrive, and without them at all if
    ;; they cannot be had.
    (let [text (fn [url] (-> (js/fetch url)
                             (.then (fn [^js r] (if (.-ok r) (.text r)
                                                    (throw (js/Error. (str url ": HTTP " (.-status r)))))))))]
      (-> (js/Promise.all #js [(text "data/bsc5.tsv") (text "data/constellations.tsv")
                               (text "data/star-names.tsv")])
          (.then (fn [^js texts]
                   (let [catalogue (stars/parse (aget texts 0))
                         constellations (stars/parse-constellations (aget texts 1))
                         names (stars/parse-names (aget texts 2))
                         ^js points (sky catalogue)
                         ^js figs (figures constellations)]
                     (reset! stars-at points)
                     (reset! lines-at figs)
                     (.add sky-group points)
                     (.add sky-group figs)
                     (reset! con-labels
                             (vec (for [{nm :name at :label} constellations]
                                    (label nm "sky-label sky-label--constellation" at))))
                     (reset! star-labels
                             (vec (for [s catalogue
                                        :let [nm (names (:hr s))]
                                        :when (and nm (<= (:vmag s) named-brighter-than))]
                                    (label nm "sky-label sky-label--star" [(:ra s) (:dec s)]))))
                     (doseq [^js o (concat @con-labels @star-labels)] (.add sky-group o)))))
          (.catch (fn [e] (js/console.warn "No star catalogue:" e)))))
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
                                     (to-view [x (* ce y) (* se y)])))
                                 0x445070 0.5))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h)
                    (.setSize labels w h))))
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
                        [mx my mz] (to-view (eph/moon mjd))
                        k (* (.-moonZoom controls) (/ view-radius outer-au) (/ 1.0 c/AU))]
                    (.set (.-position moon) (+ ex (* k mx)) (+ ey (* k my)) (+ ez (* k mz)))
                    ;; and its orbit, drawn at the same magnification
                    (let [arr (.-array (.getAttribute (.-geometry moon-line) "position"))]
                      (doseq [i (range 130)]
                        (let [[ox oy oz] (to-view (eph/moon (+ mjd (* 27.32 (/ i 129.0)))))]
                          (aset arr (* i 3) (+ ex (* k ox)))
                          (aset arr (+ (* i 3) 1) (+ ey (* k oy)))
                          (aset arr (+ (* i 3) 2) (+ ez (* k oz)))))
                      (.setDrawRange (.-geometry moon-line) 0 130)
                      (set! (.-needsUpdate (.getAttribute (.-geometry moon-line) "position")) true)))))
              (publish! []
                (let [mjd (:mjd @state)
                      [yr mo dy hr] (atime/mjd->calendar mjd)
                      moon-km (let [[a b cc] (eph/moon mjd)]
                                (js/Math.sqrt (+ (* a a) (* b b) (* cc cc))))]
                  (set! (.-visible orbits) (.-showOrbits controls))
                  (set! (.-visible ecliptic) (.-showEcliptic controls))
                  (when-let [^js s @stars-at]
                    (set! (.-visible s) (.-showStars controls)))
                  (when-let [^js l @lines-at]
                    (set! (.-visible l) (.-showLines controls)))
                  ;; Each label on its own: the label renderer asks the
                  ;; label whether it is visible, not its group.
                  (doseq [^js o @con-labels] (set! (.-visible o) (.-showConNames controls)))
                  (doseq [^js o @star-labels] (set! (.-visible o) (.-showStarNames controls)))
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
                    ;; The sky goes where the camera goes: infinitely far.
                    (.copy (.-position sky-group) (.-position camera))
                    (.render labels scene camera)
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
          (-> (.add gui controls "showStars") (.name "stars"))
          (-> (.add gui controls "starLimit" 1.0 6.5 0.1) (.name "faintest star")
              (.onChange (fn [v] (set! (.-value faintest) v))))
          (-> (.add gui controls "showLines") (.name "constellation lines"))
          (-> (.add gui controls "showConNames") (.name "constellation names"))
          (-> (.add gui controls "showStarNames") (.name "star names"))
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
