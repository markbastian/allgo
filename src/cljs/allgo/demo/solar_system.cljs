(ns allgo.demo.solar-system
  "The solar system, the Moon and a satellite drawn in one frame.

  Everything here comes from `allgo.astro` already referred to the
  mean equator and equinox of J2000, so no body needs rotating to sit
  beside another. That is the point of the demo as much as the picture is:
  the planets come from VSOP87 and the Moon from ELP, both as Meeus
  abridges them, and they can simply be added. The panel switches both to
  a coarser pair -- Standish's Keplerian elements and Montenbruck and
  Gill's short lunar series -- and the readout says, for the body in
  focus, how far apart the two put it: arcminutes, which is what the
  coarse ones are worth.

  Pluto is there too, from Meeus's chapter 37 -- a fit to a numerical
  ephemeris over 1885 to 2099 and to nothing else, so outside those years
  it is not drawn, and its orbit is drawn only along the arc the fit
  covers.

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
  as a star chart has it, until the view is turned.

  Click a body, or pick it from the menu, and the view flies to it and
  centers on it: the same model, the orbits and the rest of the planets
  still in it, but turning about that body instead of the Sun. Up close
  each wears Solar System Scope's map of it, is flattened at the poles as
  it really is, and is turned the way it faces at that moment by the
  IAU's rotation models (`allgo.astro.rotation`). Escape, or the button,
  goes back out. The maps, five megabytes in all, are fetched when the
  page opens."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.frames :as frames]
            [allgo.astro.kepler :as kep]
            [allgo.astro.moon :as lunar]
            [allgo.astro.planet-orbits :as orbits]
            [allgo.astro.planets :as pl]
            [allgo.astro.rotation :as rot]
            [allgo.astro.stars :as stars]
            [allgo.astro.time :as atime]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.demo.fps :as fps]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]
            ["three/examples/jsm/renderers/CSS2DRenderer.js" :refer [CSS2DObject CSS2DRenderer]]
            [clojure.math :as math]))

(def ^:private view-radius 26.0)
(def ^:private outer-au 31.0)

(def ^:private palette
  {:mercury 0xb0a08c :venus 0xe8c47a :earth 0x5b9bd5 :mars 0xc1553a
   :jupiter 0xd8a06a :saturn 0xe3d08a :uranus 0x8fd0d8 :neptune 0x5570d0
   :pluto 0xc9b8a4})

(def ^:private dot-size
  ;; Not to scale -- at true scale every planet is far under a pixel. Ordered
  ;; by actual radius so the giants still read as giants.
  {:sun 0.85 :moon 0.12
   :mercury 0.16 :venus 0.24 :earth 0.26 :mars 0.19
   :jupiter 0.62 :saturn 0.54 :uranus 0.38 :neptune 0.37 :pluto 0.1})

(def ^:private all-bodies (into [:sun] (concat pl/order [:moon :pluto])))

(def ^:private body-names
  {:sun "Sun" :mercury "Mercury" :venus "Venus" :earth "Earth" :moon "Moon" :mars "Mars"
   :jupiter "Jupiter" :saturn "Saturn" :uranus "Uranus" :neptune "Neptune" :pluto "Pluto"})

(def ^:private focus-order
  [:sun :mercury :venus :earth :moon :mars :jupiter :saturn :uranus :neptune :pluto])

(def ^:private texture-file
  ;; Solar System Scope's 2K maps (CC BY 4.0; see the README). Venus wears
  ;; its clouds, which are what anyone has ever seen of it.
  {:sun "2k_sun.jpg" :mercury "2k_mercury.jpg" :venus "2k_venus_atmosphere.jpg"
   :earth "2k_earth_daymap.jpg" :moon "2k_moon.jpg" :mars "2k_mars.jpg"
   :jupiter "2k_jupiter.jpg" :saturn "2k_saturn.jpg" :uranus "2k_uranus.jpg"
   :neptune "2k_neptune.jpg"})

(def ^:private saturn-rings
  "The main rings' extent, km from Saturn's center: the C ring's inner edge
  to the A ring's outer one."
  [74658.0 136775.0])

(def ^:private flight-ms 1400.0)

(def ^:private ^js controls
  #js {:focus         "Overview"
       :ephemeris     "VSOP87 + ELP"
       :scale         "logarithmic"
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

(def ^:private precise-source "VSOP87 + ELP")

(defn- precise? [] (= precise-source (.-ephemeris controls)))

(defn- ephemeris
  "Where the bodies are at `mjd`, from whichever source is chosen: the
  heliocentric positions of the planets and the geocentric position of the
  Moon, km, in EME2000. VSOP87 and ELP are Meeus's, good to an arcsecond
  or so; the elements and the short lunar series are good to arcminutes,
  and the readout says how far apart the two are for the body in focus."
  ([mjd] (ephemeris mjd (precise?)))
  ([mjd precise]
   (if precise
     (let [m (vsop87/ecliptic-of-date->J2000 mjd)]
       {:planet (fn [k] (mapv #(* c/AU %) (vsop87/equatorial-J2000 k mjd m)))
        :moon   (fn [] (lunar/geocentric-J2000 mjd m))})
     {:planet (fn [k] (pl/heliocentric k mjd))
      :moon   (fn [] (eph/moon mjd))})))

(defn- pluto-position
  "Pluto's heliocentric position, km, in EME2000, or nil outside 1885-2099:
  Meeus's Pluto is a fit to a numerical ephemeris over those years and
  means nothing beyond them. Both sources share it; neither has another."
  [mjd]
  (when-let [[l b r] (orbits/pluto mjd)]
    (let [ce (math/cos eph/obliquity-J2000) se (math/sin eph/obliquity-J2000)
          x (* r c/AU (math/cos b) (math/cos l))
          y (* r c/AU (math/cos b) (math/sin l))
          z (* r c/AU (math/sin b))]
      [x (- (* ce y) (* se z)) (+ (* se y) (* ce z))])))

(def ^:private pluto-span
  "The years Meeus's Pluto covers, as MJD."
  [(atime/calendar->mjd 1885 1 1) (atime/calendar->mjd 2099 12 31)])

(defn- radial
  "Map a distance in AU onto the view. Logarithmic keeps all eight orbits on
  screen at once; linear is true and needs zooming."
  [au]
  (* view-radius
     (if (= "logarithmic" (.-scale controls))
       (/ (math/log (+ 1.0 (* 3.0 au))) (math/log (+ 1.0 (* 3.0 outer-au))))
       (/ au outer-au))))

(defn- place
  "A heliocentric position in km, mapped into view coordinates. Direction is
  preserved exactly; only the radius is rescaled, so angles stay true
  whichever mode is chosen."
  [[x y z]]
  (let [r  (math/sqrt (+ (* x x) (* y y) (* z z)))]
    (if (< r 1e-9)
      [0.0 0.0 0.0]
      (let [s (/ (radial (/ r c/AU)) r)]
        (to-view [(* s x) (* s y) (* s z)])))))

(defn- orbit-points
  "One full revolution of a planet, sampled in mean anomaly rather than in
  time, so a slow outer planet costs no more samples than a fast inner one.
  The orbit is the mean one of the chosen source: Meeus's mean elements
  referred to J2000 alongside VSOP87, Standish's alongside his."
  [planet mjd n]
  (let [el (if (precise?)
             (update (orbits/mean-elements-J2000 planet mjd) :a * c/AU)
             (pl/elements-at planet mjd))
        ce (math/cos eph/obliquity-J2000)
        se (math/sin eph/obliquity-J2000)]
    (for [k (range (inc n))]
      (let [M  (* 2.0 math/PI (/ k n))
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
            glow (min 1.0 (+ 0.4 (* 0.6 (math/sqrt (/ (stars/flux v) (stars/flux 1.0))))))]
        (aset pos (* 3 i) (* sky-radius x))
        (aset pos (+ 1 (* 3 i)) (* sky-radius y))
        (aset pos (+ 2 (* 3 i)) (* sky-radius z))
        (aset col (* 3 i) (* glow r))
        (aset col (+ 1 (* 3 i)) (* glow g))
        (aset col (+ 2 (* 3 i)) (* glow b))
        (aset size i (max 2.4 (* 9.0 (math/pow 10.0 (* -0.2 (+ v 1.46))))))
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

(defn- asset
  "The URL of `path` in the site, wherever the page asking for it is.

  Relative to the bundle, not to the page: the demo page is at the top of
  the site and this demo's own full-window page is one folder down, and
  a plain `data/bsc5.tsv` would be looked for next to each. The bundle is
  in one place, `js/compiled/allgo.js` under the site's root, so that is
  where the root is found."
  [path]
  (let [src (some (fn [^js s] (let [u (.-src s)] (when (re-find #"js/compiled/allgo\.js" u) u)))
                  (array-seq (js/document.getElementsByTagName "script")))]
    (if src
      (str (subs src 0 (.indexOf src "js/compiled/allgo.js")) path)
      path)))

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

(defn- orient!
  "Turns `mesh` to face the way `body` does at `mjd`. The mesh is built
  with its pole along y and its prime meridian along x -- the body-fixed
  frame carried into three.js's axes by the same `y-up` as everything
  else -- so the rotation it needs is the model's, conjugated by `y-up`."
  [^js mesh body mjd]
  (let [r (rot/body->icrf body mjd)
        cols (mapv (fn [j] (to-view (mapv #(nth % j) r))) (range 3))
        ;; body x, y, z axes in view coordinates; local x = X, y = Z, z = -Y
        [bx by bz] cols
        lx bx ly bz lz (mapv - by)
        m (THREE/Matrix4.)]
    (.set m (lx 0) (ly 0) (lz 0) 0
          (lx 1) (ly 1) (lz 1) 0
          (lx 2) (ly 2) (lz 2) 0
          0 0 0 1)
    (.setFromRotationMatrix (.-quaternion mesh) m)))

(defn- shape!
  "Sizes `mesh` to an equatorial radius of `r` view units, flattened at the
  poles as `body` is."
  [^js mesh body r]
  (let [[a _ c] (:radii (rot/models body))]
    (.set (.-scale mesh) r (* r (/ c a)) r)))

(defn- ring
  "Saturn's rings, as a flat annulus in its equatorial plane, in units of
  Saturn's equatorial radius. The texture is a strip running outward, so
  each vertex's u is its fraction of the way across."
  []
  (let [a (rot/radius :saturn)
        [r0 r1] (map #(/ % a) saturn-rings)
        geo (THREE/RingGeometry. r0 r1 160 1)
        pos (.getAttribute geo "position")
        uv  (.getAttribute geo "uv")]
    (dotimes [i (.-count pos)]
      (let [x (.getX pos i) y (.getY pos i)]
        (.setXY uv i (/ (- (math/sqrt (+ (* x x) (* y y))) r0) (- r1 r0)) 0.5)))
    (.rotateX geo (- (/ math/PI 2.0)))
    (THREE/Mesh. geo (THREE/MeshStandardMaterial.
                      #js {:color 0xd8c8a0 :transparent true :opacity 0.5
                           :side THREE/DoubleSide :depthWrite false
                           :roughness 1.0 :metalness 0.0}))))

(def ^:private named-brighter-than
  "The faintest star given its name: the IAU has named 333 of the
  catalogue's stars, and all of them at once would bury the sky in
  text. The ninety or so brighter than this are the ones a person knows."
  2.5)

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.05 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        sphere   (THREE/SphereGeometry. 1.0 64 40)
        ;; Every body a unit sphere, sized and turned each frame. The Sun
        ;; makes its own light; the rest are lit by it.
        meshes   (into {} (for [k all-bodies]
                            [k (THREE/Mesh. sphere
                                            (if (= k :sun)
                                              (THREE/MeshBasicMaterial. #js {:color 0xffcf5c})
                                              (THREE/MeshStandardMaterial.
                                               #js {:color (get palette k 0xcfd6e0)
                                                    :roughness 1.0 :metalness 0.0})))]))
        bodies   (select-keys meshes pl/order)
        rings    (ring)
        ;; The Sun is a point of light at the center, and since the radial
        ;; scale keeps every direction from the Sun true, it lights each
        ;; planet from the side it really does. (The Moon, pushed out from
        ;; the Earth by its magnification, is lit a few degrees wrong.)
        sunlight (THREE/PointLight. 0xffffff 1.3 0 0)
        ambient  (THREE/AmbientLight. 0xffffff 0.07)
        loaded   (atom #{})
        flight   (atom nil)
        moon-path (atom nil)
        focus-ctl (atom nil)
        back     (doto (js/document.createElement "button")
                   (-> .-className (set! "focus-back"))
                   (-> .-type (set! "button"))
                   (-> .-textContent (set! "\u2190 Solar system"))
                   (-> .-hidden (set! true)))
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
      (-> (js/Promise.all #js [(text (asset "data/bsc5.tsv")) (text (asset "data/constellations.tsv"))
                               (text (asset "data/star-names.tsv"))])
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
    (.appendChild container back)
    (doseq [[_ m] meshes] (.add scene m))
    (.add (bodies :saturn) rings)
    (.add scene sunlight)
    (.add scene ambient)
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
                ;; Pluto's is the arc its theory covers -- 214 of its 248
                ;; years -- and no more: the gap is honest.
                (let [[t0 t1] pluto-span]
                  (.add orbits (line-of (keep #(some-> (pluto-position %) place)
                                              (range t0 t1 60.0))
                                        (palette :pluto) 0.35)))
                (doseq [^js child (vec (.-children ecliptic))] (.remove ecliptic child))
                ;; A ring in the ecliptic plane, tilted into the equatorial
                ;; frame -- the visible reminder of which frame this is.
                (let [ce (math/cos eph/obliquity-J2000)
                      se (math/sin eph/obliquity-J2000)
                      r  (radial 31.0)]
                  (.add ecliptic
                        (line-of (for [k (range 129)]
                                   (let [a (* 2.0 math/PI (/ k 128.0))
                                         x (* r (math/cos a))
                                         y (* r (math/sin a))]
                                     (to-view [x (* ce y) (* se y)])))
                                 0x445070 0.5))))
              (focused []
                (some (fn [[k n]] (when (= n (.-focus controls)) k)) body-names))
              (load-textures! [k]
                (let [pair (if (#{:earth :moon} k) [:earth :moon] [k])]
                  (doseq [b pair :when (and (texture-file b) (not (@loaded b)))]
                    (swap! loaded conj b)
                    (.load (THREE/TextureLoader.) (asset (str "textures/" (texture-file b)))
                           (fn [^js tex]
                             (set! (.-anisotropy tex) (.. renderer -capabilities (getMaxAnisotropy)))
                             (let [^js mat (.-material (meshes b))]
                               (set! (.-map mat) tex)
                               (.setHex (.-color mat) 0xffffff)
                               (set! (.-needsUpdate mat) true)))
                           nil
                           (fn [e] (js/console.warn "No texture for" (name b) e))))
                  (when (and (= k :saturn) (not (@loaded :rings)))
                    (swap! loaded conj :rings)
                    (.load (THREE/TextureLoader.) (asset "textures/2k_saturn_ring_alpha.png")
                           (fn [^js tex]
                             (let [^js mat (.-material rings)]
                               (set! (.-map mat) tex)
                               (.setHex (.-color mat) 0xffffff)
                               (set! (.-opacity mat) 1.0)
                               (set! (.-needsUpdate mat) true)))))))
              (view-direction
                ;; Which side to look at a body from: the side the Earth
                ;; sees, so the Moon shows its near side at tonight's phase
                ;; and Saturn its rings at the tilt a telescope would. The
                ;; Earth itself is seen from sunward and a little north,
                ;; and the overview from wherever the view already is.
                [k]
                (let [{:keys [planet moon]} (ephemeris (:mjd @state))
                      v (fn [[x y z]] (.normalize (THREE/Vector3. x y z)))
                      earth (to-view (planet :earth))]
                  (case k
                    nil nil
                    :pluto (v (mapv - earth (to-view (or (pluto-position (:mjd @state)) [1.0 0.0 0.0]))))
                    :sun (v earth)
                    :earth (.normalize (.add (.negate (v earth)) (THREE/Vector3. 0 0.6 0)))
                    :moon (.negate (v (to-view (moon))))
                    (v (mapv - earth (to-view (planet k)))))))
              (set-focus!
                ;; Through the panel's own control, so the menu shows it
                ;; and its change handler does the flying.
                [label]
                (when-let [^js ctl @focus-ctl] (.setValue ctl label)))
              (pick
                ;; The body under a point on the screen, if any: the nearest
                ;; whose disc -- or a fingertip around its center, for the
                ;; ones a few pixels across -- contains it.
                [cx cy]
                (let [^js el (.-domElement renderer)
                      rect (.getBoundingClientRect el)
                      w (.-width rect) h (.-height rect)
                      x (- cx (.-left rect)) y (- cy (.-top rect))
                      per-unit (/ (* 0.5 h) (math/tan (* 0.5 (/ (* math/PI (.-fov camera)) 180.0))))]
                  (->> meshes
                       (keep (fn [[k ^js m]]
                               (when (.-visible m)
                                 (let [p (.clone (.-position m))
                                       dist (.distanceTo p (.-position camera))
                                       v (.project p camera)]
                                   (when (< (.-z v) 1.0)
                                     (let [sx (* 0.5 (+ 1.0 (.-x v)) w)
                                           sy (* 0.5 (- 1.0 (.-y v)) h)
                                           r (max 18.0 (/ (* (.. m -scale -x) per-unit) dist))
                                           d (math/hypot (- sx x) (- sy y))]
                                       (when (<= d r) [k d])))))))
                       (sort-by second)
                       ffirst)))
              (refocus! []
                (set! (.-hidden back) (nil? (focused)))
                (when-let [k (focused)] (load-textures! k))
                (let [target (.clone (.-target orbit))
                      offset (.sub (.clone (.-position camera)) target)]
                  (reset! flight {:t0 (js/performance.now) :from target :offset offset
                                  :toward (view-direction (focused))})))
              (focus-point
                ;; Where the view should be centered: the focused body, or
                ;; the Sun for the overview.
                []
                (if-let [k (focused)]
                  (.clone (.-position (meshes k)))
                  (THREE/Vector3. 0 0 0)))
              (follow! []
                (let [p (focus-point)]
                  (if-let [{:keys [t0 ^js from ^js offset ^js toward]} @flight]
                    ;; Flying: the target slides to the body and the camera
                    ;; closes to viewing distance along the line it looked
                    ;; down, eased at both ends.
                    (let [u (min 1.0 (/ (- (js/performance.now) t0) flight-ms))
                          e (- (* 3.0 u u) (* 2.0 u u u))
                          dist (if-let [k (focused)] (* 4.0 (dot-size k)) 48.0)
                          len (+ (.length offset) (* e (- dist (.length offset))))
                          target (.lerp (.clone from) p e)
                          dir0 (.normalize (.clone offset))
                          dir (if toward
                                (let [d (.lerp (.clone dir0) toward e)]
                                  (if (< (.length d) 1e-6) toward (.normalize d)))
                                dir0)]
                      (.copy (.-target orbit) target)
                      (.copy (.-position camera) (.add (.multiplyScalar dir len) target))
                      (when (>= u 1.0) (reset! flight nil)))
                    ;; Following: carried along with the body, however the
                    ;; view has been turned since.
                    (let [delta (.sub (.clone p) (.-target orbit))]
                      (.add (.-position camera) delta)
                      (.copy (.-target orbit) p))))
                (set! (.-minDistance orbit) (if-let [k (focused)] (* 1.2 (dot-size k)) 0.0)))
              (layout!
                ;; Every body turned to face the way it does now, and sized.
                []
                (let [mjd (:mjd @state)]
                  (doseq [[b ^js m] meshes]
                    (orient! m b mjd)
                    (shape! m b (dot-size b)))
                  (set! (.-visible orbits) (.-showOrbits controls))
                  (set! (.-visible ecliptic) (.-showEcliptic controls))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h)
                    (.setSize labels w h))))
              (update-bodies! []
                (let [mjd (:mjd @state)
                      {:keys [planet moon]} (ephemeris mjd)
                      earth (planet :earth)]
                  (doseq [k pl/order]
                    (let [[x y z] (place (if (= k :earth) earth (planet k)))
                          ^js m (bodies k)]
                      (.set (.-position m) x y z)))
                  (let [^js m (meshes :pluto)]
                    (if-let [p (pluto-position mjd)]
                      (let [[x y z] (place p)]
                        (set! (.-visible m) true)
                        (.set (.-position m) x y z))
                      (set! (.-visible m) false)))
                  ;; The Moon is geocentric already; adding it to the Earth's
                  ;; heliocentric place is legitimate only because both are in
                  ;; the same frame, which is the whole point.
                  (let [[ex ey ez] (place earth)
                        [mx my mz] (to-view (moon))
                        k (* (.-moonZoom controls) (/ view-radius outer-au) (/ 1.0 c/AU))
                        precise (precise?)]
                    (.set (.-position (meshes :moon)) (+ ex (* k mx)) (+ ey (* k my)) (+ ez (* k mz)))
                    ;; and its orbit, drawn at the same magnification: a
                    ;; month ahead, sampled afresh only when the Moon has
                    ;; moved on by a sample or the source has changed --
                    ;; the ELP series is 120 terms a point
                    (let [key [precise (math/round (/ mjd 0.2))]]
                      (when (not= key (:key @moon-path))
                        (reset! moon-path
                                {:key key
                                 :points (vec (for [i (range 130)]
                                                (let [t (+ mjd (* 27.32 (/ i 129.0)))]
                                                  (to-view ((:moon (ephemeris t precise)))))))})))
                    (let [arr (.-array (.getAttribute (.-geometry moon-line) "position"))]
                      (doseq [[i [ox oy oz]] (map-indexed vector (:points @moon-path))]
                        (aset arr (* i 3) (+ ex (* k ox)))
                        (aset arr (+ (* i 3) 1) (+ ey (* k oy)))
                        (aset arr (+ (* i 3) 2) (+ ez (* k oz))))
                      (.setDrawRange (.-geometry moon-line) 0 130)
                      (set! (.-needsUpdate (.getAttribute (.-geometry moon-line) "position")) true)))))
              (sources-apart
                ;; How far apart the two ephemerides put `k` as seen from
                ;; the Earth -- the Sun and the Earth itself as seen from
                ;; the other -- which is the honest measure of the coarse
                ;; one's error.
                [k mjd]
                (let [a (ephemeris mjd true) b (ephemeris mjd false)
                      seen (fn [{:keys [planet moon]}]
                             (case k
                               :moon (moon)
                               (:sun :earth) (planet :earth)
                               (mapv - (planet k) (planet :earth))))
                      u (seen a) v (seen b)
                      cosang (/ (reduce + (map * u v))
                                (math/sqrt (* (reduce + (map * u u)) (reduce + (map * v v)))))
                      arcmin (/ (math/acos (min 1.0 cosang)) (/ math/PI 10800.0))]
                  (str "\n" (if (= k :moon) "The short lunar series is " "Standish's elements are ")
                       (.toFixed arcmin 1) "' from " (if (= k :moon) "ELP" "VSOP87"))))
              (publish! []
                (let [mjd (:mjd @state)
                      [yr mo dy hr] (atime/mjd->calendar mjd)
                      moon-km (let [[a b cc] ((:moon (ephemeris mjd)))]
                                (math/sqrt (+ (* a a) (* b b) (* cc cc))))]
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
                             (if-let [k (focused)]
                               (let [h (rot/sidereal-day k)]
                                 (str "\n" (body-names k) ": radius " (.toLocaleString (math/round (rot/radius k))) " km,"
                                      " turns in " (if (< (abs h) 48) (str (.toFixed (abs h) 2) " h") (str (.toFixed (/ (abs h) 24) 2) " d"))
                                      (when (neg? h) " (backward)")
                                      (if (= k :pluto)
                                        (when-not (pluto-position mjd) "\nnot shown: its theory covers 1885-2099 only")
                                        (sources-apart k mjd))))
                               (str "\n" (.-scale controls) " radial scale"
                                    "   Moon at " (.toFixed moon-km 0) " km, drawn "
                                    (.toFixed (.-moonZoom controls) 0) "x"))))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)
                        {:keys [last]} @state]
                    (when (and (.-running controls) last)
                      (swap! state update :mjd + (* (.-daysPerSecond controls) (/ (- t0 last) 1000.0))))
                    (swap! state assoc :last t0)
                    (update-bodies!)
                    (layout!)
                    (follow!)
                    (publish!)
                    (.update orbit)
                    ;; The sky goes where the camera goes: infinitely far.
                    (.copy (.-position sky-group) (.-position camera))
                    (.render labels scene camera)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        ;; A click or tap on a body focuses it; a press that moves is the
        ;; view being turned, and is left to the orbit controls.
        (let [^js el (.-domElement renderer)
              down (atom nil)]
          (.addEventListener el "pointerdown"
                             (fn [^js e] (reset! down [(.-clientX e) (.-clientY e) (js/performance.now)])))
          (.addEventListener el "pointerup"
                             (fn [^js e]
                               (when-let [[x0 y0 t0] @down]
                                 (reset! down nil)
                                 (when (and (< (math/hypot (- (.-clientX e) x0) (- (.-clientY e) y0)) 6)
                                            (< (- (js/performance.now) t0) 500))
                                   (when-let [k (pick (.-clientX e) (.-clientY e))]
                                     (when (not= k (focused)) (set-focus! (body-names k))))))))
          (.addEventListener el "pointermove"
                             (fn [^js e]
                               (when (and (= "mouse" (.-pointerType e)) (zero? (.-buttons e)))
                                 (set! (.. el -style -cursor)
                                       (if (pick (.-clientX e) (.-clientY e)) "pointer" ""))))))
        (.addEventListener back "click" (fn [] (set-focus! "Overview")))
        (.addEventListener js/window "keydown"
                           (fn [^js e]
                             (when (and (= "Escape" (.-key e)) (focused) @running?
                                        (not (#{"INPUT" "SELECT" "TEXTAREA"} (.. e -target -tagName))))
                               (set-focus! "Overview"))))
        (rebuild-orbits!)
        ;; Every map up front, so a body is dressed however it is reached --
        ;; by the focus menu or by zooming in on it by hand.
        (doseq [k focus-order] (load-textures! k))
        (let [gui (GUI. #js {:container container})]
          (reset! focus-ctl
                  (-> (.add gui controls "focus" (clj->js (into ["Overview"] (map body-names focus-order))))
                      (.onChange refocus!)))
          (-> (.add gui controls "ephemeris" #js [precise-source "Keplerian elements"])
              (.onChange (fn [_] (reset! moon-path nil) (rebuild-orbits!))))
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
          (.add gui #js {:toJ2000 (fn [] (swap! state assoc :mjd c/mjd-J2000) (rebuild-orbits!))} "toJ2000")
          ;; On a phone the panel would cover half the sky; it starts
          ;; closed, a tap on its title away.
          (when (.-matches (js/window.matchMedia "(max-width: 640px)"))
            (.close gui)))
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
