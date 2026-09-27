(ns allgo.demo.solar-sky
  "The sky as seen from one place on the Earth, for the solar system demo.

  Everything is turned into the observer's horizon frame: the stars from
  J2000 by precession and nutation to the true equator of date, then by
  the local apparent sidereal time and the latitude; the Sun, the Moon and
  the planets from their apparent places of date, the Moon topocentric --
  seen from the observer, not from the Earth's center, which moves it by
  up to a degree. Near the horizon everything is lifted by refraction, the
  stars in their shader and the bodies by the same formula (Saemundsson's,
  Meeus 16.4), so a setting Sun is still seen whole when it has
  geometrically gone.

  The Moon is a lit sphere at its true angular size, so zooming in -- the
  wheel narrows the field of view -- shows its phase and its librations as
  they are. The planets are points, brightest largest, with their names.
  Below the horizon is ground; above it the sky brightens through
  twilight as the Sun rises, drowning the stars.

  Three.js axes here: x east, y up, z south."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.frames :as frames]
            [allgo.astro.geodesy :as geodesy]
            [allgo.astro.illumination :as illum]
            [allgo.astro.moon :as moon]
            [allgo.astro.rise :as rise]
            [allgo.astro.rotation :as rot]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as atime]
            [allgo.numerics.linear :as lin]
            ["three" :as THREE]
            ["three/examples/jsm/renderers/CSS2DRenderer.js" :refer [CSS2DObject CSS2DRenderer]]
            [clojure.math :as math]))

(def ^:private planets [:mercury :venus :mars :jupiter :saturn :uranus :neptune])

(def ^:private names
  {:sun "Sun" :moon "Moon" :mercury "Mercury" :venus "Venus" :mars "Mars" :jupiter "Jupiter"
   :saturn "Saturn" :uranus "Uranus" :neptune "Neptune"})

(def ^:private body-radius 600.0)
(def ^:private ground-radius 400.0)

;; ------------------------------------------------------------ the frames

(defn- horizon-matrix
  "Rows taking true-of-date equatorial coordinates into three.js's horizon
  axes, x east, y up, z south, for latitude `lat` and local apparent
  sidereal time `lst`."
  [lat lst]
  (let [sp (math/sin lat) cp (math/cos lat)
        ;; into the frame of the local meridian, then east-north-up
        enu (lin/mat-mul [[0.0 1.0 0.0] [(- sp) 0.0 cp] [cp 0.0 sp]] (frames/rz lst))
        [e n u] enu]
    [e u (mapv - n)]))

(defn- refracted
  "A three.js horizon direction lifted by refraction toward the zenith."
  [[x y z]]
  (let [alt (math/asin (max -1.0 (min 1.0 y)))]
    (if (< alt (* -1.0 c/degrees))
      [x y z]
      (let [alt' (+ alt (coord/refraction alt))
            h (math/hypot x z)
            [ax az] (if (< h 1e-9) [1.0 0.0] [(/ x h) (/ z h)])]
        [(* ax (math/cos alt')) (math/sin alt') (* az (math/cos alt'))]))))

(defn- unit [[ra dec]]
  [(* (math/cos dec) (math/cos ra)) (* (math/cos dec) (math/sin ra)) (math/sin dec)])

(def ^:private refraction-glsl
  "Saemundsson's refraction in GLSL, the same formula `coordinates/
  refraction` evaluates for the bodies."
  "vec3 refractDir(vec3 d) {
     float alt = asin(clamp(d.y, -1.0, 1.0));
     float h = degrees(alt);
     if (h < -1.0) return d;
     float R = 1.02 / tan(radians(h + 10.3 / (h + 5.11))) / 60.0;
     float alt2 = alt + radians(R);
     float horiz = length(d.xz);
     vec2 az = horiz > 1e-6 ? d.xz / horiz : vec2(1.0, 0.0);
     return vec3(az.x * cos(alt2), sin(alt2), az.y * cos(alt2));
   }")

(defn- star-material
  "The main demo's star shader, with each star refracted after the sky is
  turned onto the horizon."
  [limit]
  (THREE/ShaderMaterial.
   #js {:uniforms #js {:limit limit :scale #js {:value (or js/window.devicePixelRatio 1)}}
        :vertexColors true :transparent true :depthWrite false :blending THREE/AdditiveBlending
        :vertexShader (str "uniform float limit; uniform float scale;
                            attribute float size; attribute float vmag;
                            varying vec3 vColor; varying float vHidden;\n"
                           refraction-glsl
                           "void main() {
                              vColor = color;
                              vHidden = vmag > limit ? 1.0 : 0.0;
                              vec3 w = (modelMatrix * vec4(position, 1.0)).xyz;
                              float r = length(w);
                              vec3 p = refractDir(w / r) * r;
                              gl_Position = projectionMatrix * viewMatrix * vec4(p, 1.0);
                              gl_PointSize = size * scale;
                            }")
        :fragmentShader "varying vec3 vColor; varying float vHidden;
                         void main() {
                           if (vHidden > 0.5) discard;
                           float d = length(gl_PointCoord - 0.5);
                           if (d > 0.5) discard;
                           float a = 1.0 - smoothstep(0.25, 0.5, d);
                           gl_FragColor = vec4(vColor * a, a);
                         }"}))

(defn- line-material []
  (THREE/ShaderMaterial.
   #js {:transparent true :depthWrite false
        :vertexShader (str refraction-glsl
                           "void main() {
                              vec3 w = (modelMatrix * vec4(position, 1.0)).xyz;
                              float r = length(w);
                              gl_Position = projectionMatrix * viewMatrix * vec4(refractDir(w / r) * r, 1.0);
                            }")
        :fragmentShader "void main() { gl_FragColor = vec4(0.36, 0.48, 0.69, 0.45); }"}))

(defn- label [text cls]
  (let [el (doto (js/document.createElement "div")
             (-> .-className (set! cls))
             (.appendChild (doto (js/document.createElement "span") (-> .-textContent (set! text)))))]
    (CSS2DObject. el)))

(defn- clock
  "Hours and minutes of a day fraction, or a dash when there is no event."
  [f]
  (if (or (nil? f) (js/isNaN f))
    "  --  "
    (let [m (math/round (* 1440.0 (- f (math/floor f))))]
      (str (.padStart (str (mod (quot m 60) 24)) 2 "0") ":" (.padStart (str (mod m 60)) 2 "0")))))

;; ------------------------------------------------------------ the view

(defn create
  "The sky view, drawing into `renderer` in `container` alongside the main
  demo. Returns a map of functions: `:set-stars!` hands it the catalogue's
  geometry once loaded, `:frame!` draws it for an MJD (TT) and place
  `[lat lon]` (radians) and returns the readout text, `:show!` turns it
  on and off, and `:resize!`."
  [^js container ^js renderer limit]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 60 1 0.1 3000)
        labels (doto (CSS2DRenderer.)
                 (-> .-domElement .-className (set! "sky-labels")))
        _ (set! (.. labels -domElement -style -display) "none")
        _ (.appendChild container (.-domElement labels))
        sky (doto (THREE/Group.) (-> .-matrixAutoUpdate (set! false)))
        ground (THREE/Mesh. (THREE/SphereGeometry. ground-radius 64 16 0 (* 2 math/PI) (/ math/PI 2) (/ math/PI 2))
                            (THREE/MeshBasicMaterial. #js {:color 0x10140f :side THREE/BackSide}))
        horizon-line (THREE/Line. (doto (THREE/BufferGeometry.)
                                    (.setFromPoints (clj->js (for [k (range 129)]
                                                               (let [a (* 2 math/PI (/ k 128))]
                                                                 (THREE/Vector3. (* 0.999 ground-radius (math/cos a)) 0.0
                                                                                 (* 0.999 ground-radius (math/sin a))))))))
                                  (THREE/LineBasicMaterial. #js {:color 0x3a4a3a}))
        sunlight (THREE/DirectionalLight. 0xffffff 2.2)
        sun (THREE/Mesh. (THREE/SphereGeometry. 1.0 32 16) (THREE/MeshBasicMaterial. #js {:color 0xfff1c0}))
        moon-mesh (THREE/Mesh. (THREE/SphereGeometry. 1.0 64 32)
                               (THREE/MeshStandardMaterial. #js {:color 0xcfd6e0 :roughness 1.0 :metalness 0.0}))
        planet-points (THREE/Points. (doto (THREE/BufferGeometry.)
                                       (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* 3 (count planets))) 3))
                                       (.setAttribute "size" (THREE/BufferAttribute. (js/Float32Array. (count planets)) 1)))
                                     (THREE/ShaderMaterial.
                                      #js {:transparent true :depthWrite false
                                           :vertexShader "attribute float size;
                                                          void main() {
                                                            gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
                                                            gl_PointSize = size;
                                                          }"
                                           :fragmentShader "void main() {
                                                              float d = length(gl_PointCoord - 0.5);
                                                              if (d > 0.5) discard;
                                                              gl_FragColor = vec4(1.0, 0.95, 0.85, 1.0 - smoothstep(0.3, 0.5, d));
                                                            }"}))
        body-labels (into {} (for [k (cons :moon (cons :sun planets))]
                               [k (label (names k) "sky-label sky-label--planet")]))
        compass (for [[t x z] [["N" 0 -1] ["E" 1 0] ["S" 0 1] ["W" -1 0]]]
                  (doto (label t "sky-label sky-label--compass")
                    (-> .-position (.set (* 0.98 ground-radius x) 4.0 (* 0.98 ground-radius z)))))
        stars (atom nil)
        rises (atom nil)
        look (atom {:az math/PI :alt (* 25 c/degrees)})
        dirs (atom {})
        shown (atom false)]
    (.add scene sky)
    (.add scene ground)
    (.add scene horizon-line)
    (.add scene sunlight)
    (.add scene (.-target sunlight))
    (.add scene sun)
    (.add scene moon-mesh)
    (.add scene planet-points)
    (set! (.-frustumCulled planet-points) false)
    (set! (.-renderOrder ground) 2)
    (doseq [[_ o] body-labels] (.add scene o))
    (doseq [o compass] (.add scene o))
    ;; drag to look about; the wheel zooms by narrowing the field of view
    (let [^js el (.-domElement renderer)
          drag (atom nil)]
      (.addEventListener el "pointerdown"
                         (fn [^js e] (when @shown (reset! drag [(.-clientX e) (.-clientY e) @look]))))
      (.addEventListener js/window "pointerup" (fn [_] (reset! drag nil)))
      (.addEventListener el "pointermove"
                         (fn [^js e]
                           (when-let [[x0 y0 {:keys [az alt]}] (and @shown @drag)]
                             (let [k (/ (* (.-fov camera) c/degrees) (.-clientHeight el))]
                               (reset! look {:az (- az (* k (- (.-clientX e) x0)))
                                             :alt (max (* -10 c/degrees)
                                                       (min (* 90 c/degrees) (+ alt (* k (- (.-clientY e) y0)))))})))))
      (.addEventListener el "wheel"
                         (fn [^js e]
                           (when @shown
                             (.preventDefault e)
                             (set! (.-fov camera) (max 0.5 (min 100.0 (* (.-fov camera) (math/exp (* 0.001 (.-deltaY e)))))))
                             (.updateProjectionMatrix camera)))
                         #js {:passive false}))
    {:set-stars!
     (fn [^js points ^js lines]
       (let [p (THREE/Points. (.-geometry points) (star-material limit))
             l (THREE/LineSegments. (.-geometry lines) (line-material))]
         (set! (.-frustumCulled p) false)
         (set! (.-frustumCulled l) false)
         (reset! stars {:points p :lines l})
         (.add sky p)
         (.add sky l)))

     :look-at!
     ;; turn to face a body, where it was last drawn -- and narrow the
     ;; field for the Sun and Moon, whose disks are half a degree
     (fn [k]
       (if (= k :south)
         (do (reset! look {:az math/PI :alt (* 25 c/degrees)})
             (set! (.-fov camera) 60.0)
             (.updateProjectionMatrix camera))
         (when-let [[x y z] (some-> (@dirs k) refracted)]
           (reset! look {:az (math/atan2 x (- z)) :alt (math/asin y)})
           (set! (.-fov camera) (if (#{:sun :moon} k) 2.0 40.0))
           (.updateProjectionMatrix camera))))

     :bodies (vec (cons :sun (cons :moon planets)))

     :show!
     (fn [on?]
       (reset! shown on?)
       (set! (.. labels -domElement -style -display) (if on? "" "none")))

     :resize!
     (fn [w h]
       (set! (.-aspect camera) (/ w h))
       (.updateProjectionMatrix camera)
       (.setSize labels w h))

     :frame!
     (fn [mjd [lat lon] {:keys [lines? moon-map]}]
       (let [dt (atime/delta-t (atime/decimal-year mjd))
             ut (- mjd (/ dt 86400.0))
             lst (+ (frames/gast ut mjd) lon)
             H (horizon-matrix lat lst)
             ;; J2000 -> true of date -> horizon, and the star geometry is
             ;; in the main demo's y-up axes: [x z -y] undone first
             M (lin/mat-mul H (frames/precession-nutation mjd))
             yup [[1.0 0.0 0.0] [0.0 0.0 -1.0] [0.0 1.0 0.0]]
             Ms (lin/mat-mul M yup)
             ^js m4 (.-matrix sky)
             to-sky (fn [radec] (lin/mat-vec H (unit radec)))
             rho (geodesy/parallax-constants lat 0.0)
             ;; the bodies, apparent and topocentric
             [sra sdec sR] (solar/apparent-equatorial mjd)
             sun-dir (to-sky [sra sdec])
             sun-alt (math/asin (second sun-dir))
             [mra mdec mdist] (moon/apparent-equatorial mjd)
             [mra' mdec'] (coord/topocentric [mra mdec] (/ mdist c/AU) rho lon (- lst lon))
             moon-dir (to-sky [mra' mdec'])
             planet-dirs (into {} (for [k planets]
                                    (let [[ra dec] (elliptic/apparent k mjd)] [k (to-sky [ra dec])])))
             _ (reset! dirs (assoc planet-dirs :sun sun-dir :moon moon-dir))
             ;; twilight: the sky brightens from astronomical dusk
             day (max 0.0 (min 1.0 (/ (+ (/ sun-alt c/degrees) 18.0) 23.0)))]
         (.set m4 (get-in Ms [0 0]) (get-in Ms [0 1]) (get-in Ms [0 2]) 0
               (get-in Ms [1 0]) (get-in Ms [1 1]) (get-in Ms [1 2]) 0
               (get-in Ms [2 0]) (get-in Ms [2 1]) (get-in Ms [2 2]) 0
               0 0 0 1)
         (set! (.-matrixWorldNeedsUpdate sky) true)
         (set! (.-background scene) (.lerpColors (THREE/Color.) (THREE/Color. 0x04050a) (THREE/Color. 0x5f8fc8)
                                                 (* day day)))
         (when-let [{:keys [^js points ^js lines]} @stars]
           (set! (.. points -material -uniforms -limit -value)
                 (min (.-value limit) (- 6.5 (* 9.0 day))))
           (set! (.-visible lines) (boolean lines?)))
         ;; the Sun: a disk at its true size, and the light the Moon is lit by
         (let [[x y z] (refracted sun-dir)
               r (* body-radius (math/tan (/ (* 959.63 c/arcsec) sR)))]
           (.set (.-position sun) (* body-radius x) (* body-radius y) (* body-radius z))
           (.setScalar (.-scale sun) r)
           (.set (.-position sunlight) (* body-radius (sun-dir 0)) (* body-radius (sun-dir 1))
                 (* body-radius (sun-dir 2))))
         ;; the Moon: a sphere at its true size, turned as it faces us
         (let [[x y z] (refracted moon-dir)
               ^js mat (.-material moon-mesh)
               r (* body-radius (math/tan (illum/moon-semidiameter mdist)))
               B (rot/body->icrf :moon mjd)
               ;; mesh-local axes to the body's: x, z, -y
               C [[1.0 0.0 0.0] [0.0 0.0 -1.0] [0.0 1.0 0.0]]
               R (lin/mat-mul (lin/mat-mul M B) C)
               ^js m (THREE/Matrix4.)]
           (when (and moon-map (not (.-map mat)))
             (set! (.-map mat) moon-map)
             (.setHex (.-color mat) 0xffffff)
             (set! (.-needsUpdate mat) true))
           (.set (.-position moon-mesh) (* body-radius x) (* body-radius y) (* body-radius z))
           (.setScalar (.-scale moon-mesh) r)
           (.set m (get-in R [0 0]) (get-in R [0 1]) (get-in R [0 2]) 0
                 (get-in R [1 0]) (get-in R [1 1]) (get-in R [1 2]) 0
                 (get-in R [2 0]) (get-in R [2 1]) (get-in R [2 2]) 0 0 0 0 1)
           (.setFromRotationMatrix (.-quaternion moon-mesh) m))
         ;; the planets, sized by their brightness
         (let [^js pos (.getAttribute (.-geometry planet-points) "position")
               ^js size (.getAttribute (.-geometry planet-points) "size")]
           (doseq [[i k] (map-indexed vector planets)]
             (let [[x y z] (refracted (planet-dirs k))]
               (.setXYZ pos i (* body-radius x) (* body-radius y) (* body-radius z))
               (.setX size i (if (#{:uranus :neptune} k) 2.5 (if (#{:venus :jupiter} k) 7.0 5.0)))))
           (set! (.-needsUpdate pos) true)
           (set! (.-needsUpdate size) true))
         (doseq [[k ^js o] body-labels]
           (let [[x y z] (refracted (case k :sun sun-dir :moon moon-dir (planet-dirs k)))
                 up? (> y -0.01)]
             (set! (.-visible o) up?)
             (.set (.-position o) (* body-radius x) (* body-radius y) (* body-radius z))))
         ;; where the camera looks
         (let [{:keys [az alt]} @look]
           (.set (.-position camera) 0 0 0)
           (.lookAt camera (* (math/sin az) (math/cos alt)) (math/sin alt) (* -1 (math/cos az) (math/cos alt))))
         (.render renderer scene camera)
         (.render labels scene camera)
         ;; the readout: rising and setting, worked out once per day and place
         (let [day0 (math/floor ut)
               key [day0 lat lon]]
           (when (not= key (:key @rises))
             (let [place [lat lon]]
               (reset! rises {:key key
                              :table (concat [[:sun (rise/sun place day0)] [:moon (rise/moon place day0)]]
                                             (for [k planets] [k (rise/planet k place day0)]))})))
           (let [alt-az (fn [[x y z]] (str (.toFixed (/ (math/asin y) c/degrees) 1) "° up, "
                                           (.toFixed (mod (/ (math/atan2 x (- z)) c/degrees) 360.0) 0) "° azimuth"))
                 lst-h (mod (/ lst (/ math/PI 12.0)) 24.0)]
             (str (.toFixed (/ lat c/degrees) 2) "° " (if (neg? lon) "W " "E ")
                  (.toFixed (abs (/ lon c/degrees)) 2) "°   local sidereal time "
                  (clock (/ lst-h 24.0))
                  "\nSun " (alt-az sun-dir) "   Moon " (alt-az moon-dir)
                  "\nUT          rise   transit  set"
                  (apply str (for [[k t] (:table @rises)]
                               (str "\n" (.padEnd (names k) 10) (if t (str (clock (:rise t)) "  " (clock (:transit t)) "  " (clock (:set t)))
                                                                    "   never crosses the horizon")))))))))}))
