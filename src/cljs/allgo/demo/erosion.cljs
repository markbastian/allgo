(ns allgo.demo.erosion
  "Rain falling on a fractal heightmap, one batch of droplets a frame.

  The point of watching it rather than looking at a finished picture is
  that erosion is not a filter. Nothing here is smoothing the terrain or
  drawing valleys onto it: the valleys are where the water went, and the
  water went there because of where the terrain already was. Each frame's
  droplets are steered by what the last frame's droplets cut, so the
  channels sharpen themselves. Leave it running and the fractal noise
  sorts itself into ridges and drainage with a direction.

  ## Three ways to look at it

  `terrain` shades by height and steepness, which is what the landscape
  would look like.

  `drainage` colors by `:flux` -- how much water has crossed each cell
  over the whole run, on a log scale, since a main channel carries orders
  of magnitude more than a hillside. This is the view worth watching from
  the start: the network appears within the first few frames, long before
  the terrain visibly changes, because the water finds the valleys before
  it has finished cutting them.

  `change` is height against the heightmap we began with -- red where
  material has been taken, blue where it has been put down. It separates
  the two halves of the model that the shaded view mixes together, and it
  is the one that shows sediment fanning out where the grade eases.

  ## The control that matters

  `min-slope` is the floor under the slope in the capacity formula, and
  it is scale dependent -- see `allgo.procedural.erosion`. Take it up to
  the 0.01 that is the usual published default and the demo stops being
  erosion: that floor sits above almost every slope in this heightmap, so
  every droplet everywhere carries the same load, and what you get is the
  terrain blurred rather than carved. It is worth doing once, because the
  failure looks like an effect rather than like a bug."
  (:require [allgo.demo.fps :as fps]
            [allgo.procedural.erosion :as erosion]
            [allgo.procedural.terrain :as terrain]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def config
  {:sizes [65 129 257]
   :size 129
   :width 1.0
   ;; How wide the terrain is in world units, held fixed as the grid size
   ;; changes, so resolution decides how much detail there is rather than
   ;; how large the thing is.
   :world-span 2048
   ;; Diamond-square at :width 1.0 spans about half a unit of height, so
   ;; this puts the relief at roughly a quarter of the span -- mountains
   ;; rather than a plain. Erosion does not know about it either way: the
   ;; model works in heightmap units, and this is only the camera's idea
   ;; of how tall they are.
   :height-scale 800
   ;; Per frame. The cost is about 45ms of droplets and 15ms of rewriting
   ;; the mesh at 129 on a side, which leaves room to spare at sixty
   ;; frames a second; the slider goes up for anyone who would rather
   ;; have the result than watch it arrive.
   :droplets-per-frame 800})

(defn- size->iterations [size]
  (long (js/Math.round (/ (js/Math.log (dec size)) (js/Math.log 2)))))

(defn- fresh-field
  "A new diamond-square heightmap, and a copy of it to measure against."
  [{:keys [width size]}]
  (let [g (terrain/generate {:width width
                             :iterations (size->iterations size)
                             :corners [0.0 (rand) (rand) (rand)]})
        ^js hs (:heights g)
        original (js/Float64Array. hs)]
    (assoc g :original original)))

;; ---------------------------------------------------------------------------
;; Color

(defn- lerp3 [[r1 g1 b1] [r2 g2 b2] t]
  (let [t (max 0.0 (min 1.0 t))]
    [(+ r1 (* t (- r2 r1))) (+ g1 (* t (- g2 g1))) (+ b1 (* t (- b2 b1)))]))

(def ^:private rock [0.42 0.40 0.38])
(def ^:private grass [0.30 0.40 0.22])
(def ^:private dirt [0.46 0.36 0.24])
(def ^:private snow [0.92 0.93 0.95])
(def ^:private silt [0.55 0.48 0.33])

(defn- terrain-color
  "Height picks the band, steepness pulls it toward bare rock.

  Without the steepness term the snow line is a contour and the cliffs
  wear it too, which is the one thing that stops eroded terrain from
  reading as eroded -- the whole point is that the steep faces are steep."
  [height-t slope-t]
  (let [base (cond
               (< height-t 0.30) (lerp3 grass dirt (/ height-t 0.30))
               (< height-t 0.65) (lerp3 dirt rock (/ (- height-t 0.30) 0.35))
               :else (lerp3 rock snow (min 1.0 (/ (- height-t 0.65) 0.25))))]
    (lerp3 base rock (min 1.0 (* 1.6 slope-t)))))

;; ---------------------------------------------------------------------------
;; Geometry

(defn- build-geometry
  "One vertex per cell, shared between the triangles that meet there so
  the normals interpolate instead of faceting."
  [dim]
  (let [n (* dim dim)
        quads (* (dec dim) (dec dim))
        index (js/Uint32Array. (* quads 6))]
    (dotimes [i (dec dim)]
      (dotimes [j (dec dim)]
        (let [a (+ (* i dim) j)
              b (+ (* i dim) (inc j))
              c (+ (* (inc i) dim) j)
              d (+ (* (inc i) dim) (inc j))
              base (* 6 (+ (* i (dec dim)) j))]
          (aset index base a) (aset index (+ base 1) c) (aset index (+ base 2) b)
          (aset index (+ base 3) b) (aset index (+ base 4) c) (aset index (+ base 5) d))))
    (doto (THREE/BufferGeometry.)
      (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* n 3)) 3))
      (.setAttribute "normal" (THREE/BufferAttribute. (js/Float32Array. (* n 3)) 3))
      (.setAttribute "color" (THREE/BufferAttribute. (js/Float32Array. (* n 3)) 3))
      (.setIndex (THREE/BufferAttribute. index 1)))))

(defn- refresh-geometry!
  "Rewrites every vertex from the heightmap as it stands.

  Reads the flat array directly rather than going through
  `allgo.procedural.mesh`, whose helpers take a nested vector and reach
  into it with `get-in`. That is the right shape to bake a mesh once; this
  runs every frame that rains, and at 257 on a side it is a quarter of a
  million lookups a frame."
  [^js geometry {:keys [^js heights ^js original ^js flux dim]} opts]
  (let [{:keys [show exaggeration cell-scale height-scale]} opts
        dim (long dim)
        n (* dim dim)
        ^js pos (.-array (.getAttribute geometry "position"))
        ^js nrm (.-array (.getAttribute geometry "normal"))
        ^js col (.-array (.getAttribute geometry "color"))
        offset (* 0.5 (dec dim) cell-scale)
        y-scale (* height-scale exaggeration)
        at (fn [i j] (aget heights (+ (* (min (dec dim) (max 0 i)) dim)
                                      (min (dec dim) (max 0 j)))))
        ;; One pass for the ranges the coloring needs, so the second pass
        ;; can be a straight write.
        [lo hi] (loop [k 0 lo ##Inf hi ##-Inf]
                  (if (= k n)
                    [lo hi]
                    (let [v (aget heights k)]
                      (recur (inc k) (min lo v) (max hi v)))))
        span (max 1e-9 (- hi lo))
        max-flux (if flux
                   (loop [k 0 m 0.0] (if (= k n) m (recur (inc k) (max m (aget flux k)))))
                   0.0)
        log-max (js/Math.log (+ 1.0 max-flux))
        max-change (loop [k 0 m 1e-9]
                     (if (= k n)
                       m
                       (recur (inc k) (max m (js/Math.abs (- (aget heights k) (aget original k)))))))
        drainage? (= show "drainage")
        change? (= show "change")]
    (dotimes [i dim]
      (dotimes [j dim]
        (let [idx (+ (* i dim) j)
              base (* idx 3)
              h (aget heights idx)
              ;; Central differences, in world units, so the normal and
              ;; the steepness agree with the mesh actually drawn.
              dzdx (/ (* y-scale (- (at i (inc j)) (at i (dec j)))) (* 2.0 cell-scale))
              dzdz (/ (* y-scale (- (at (inc i) j) (at (dec i) j))) (* 2.0 cell-scale))
              len (js/Math.sqrt (+ (* dzdx dzdx) 1.0 (* dzdz dzdz)))
              slope-t (min 1.0 (/ (js/Math.sqrt (+ (* dzdx dzdx) (* dzdz dzdz))) 2.0))
              [r g b] (cond
                        drainage?
                        ;; Log, because a trunk channel carries orders of
                        ;; magnitude more than the slope beside it and a
                        ;; linear scale shows one bright thread and
                        ;; nothing else.
                        (let [t (if (pos? log-max)
                                  (/ (js/Math.log (+ 1.0 (aget flux idx))) log-max)
                                  0.0)]
                          (lerp3 [0.13 0.15 0.20] [0.45 0.80 1.0] (js/Math.pow t 1.5)))

                        change?
                        (let [d (/ (- h (aget original idx)) max-change)]
                          (if (neg? d)
                            (lerp3 [0.55 0.55 0.55] [0.85 0.22 0.18] (min 1.0 (* -3.0 d)))
                            (lerp3 [0.55 0.55 0.55] [0.25 0.55 0.95] (min 1.0 (* 3.0 d)))))

                        :else
                        (let [t (/ (- h lo) span)
                              c (terrain-color t slope-t)]
                          ;; A wash of silt in the channels, so the shaded
                          ;; view carries some of what the drainage view
                          ;; says outright.
                          (if (and flux (pos? log-max))
                            (lerp3 c silt
                                   (max 0.0 (- (/ (js/Math.log (+ 1.0 (aget flux idx))) log-max)
                                               0.55)))
                            c)))]
          (aset pos base (- (* j cell-scale) offset))
          (aset pos (+ base 1) (* h y-scale))
          (aset pos (+ base 2) (- (* i cell-scale) offset))
          (aset nrm base (/ (- dzdx) len))
          (aset nrm (+ base 1) (/ 1.0 len))
          (aset nrm (+ base 2) (/ (- dzdz) len))
          (aset col base r)
          (aset col (+ base 1) g)
          (aset col (+ base 2) b))))
    (set! (.-needsUpdate (.getAttribute geometry "position")) true)
    (set! (.-needsUpdate (.getAttribute geometry "normal")) true)
    (set! (.-needsUpdate (.getAttribute geometry "color")) true)
    (.computeBoundingSphere geometry)))

;; ---------------------------------------------------------------------------
;; Scene

(defn init! [^js container config]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 55 1 1 100000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        span (:world-span config)
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom {:size (:size config)})
        controls #js {:rain true
                      :droplets (:droplets-per-frame config)
                      :show "drainage"
                      :size (:size config)
                      :exaggeration 1.0
                      ;; No hyphen in this name, and that is not a style
                      ;; choice. `#js {}` keeps the key it is given, but
                      ;; `(.-min-slope o)` compiles to `o.min_slope`, so a
                      ;; hyphenated property writes one name and reads
                      ;; another -- lil-gui sets "min-slope" and the read
                      ;; comes back undefined. It reached the solver as a
                      ;; NaN, which spread to every height and left the
                      ;; geometry with nothing finite to draw. The label
                      ;; below puts the hyphen back for the reader.
                      :minSlope (:min-slope erosion/defaults)
                      :capacity (:capacity erosion/defaults)
                      :erosion (:erosion erosion/defaults)
                      :deposition (:deposition erosion/defaults)
                      :inertia (:inertia erosion/defaults)
                      :radius (:radius erosion/defaults)
                      :border (:border erosion/defaults)
                      :rained "0.0 per cell"
                      :reset (fn [])}]
    (set! (.-background scene) (THREE/Color. 0x0b0d14))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (.add scene (THREE/AmbientLight. 0xffffff 0.55))
    (let [sun (THREE/DirectionalLight. 0xfff2dc 1.15)]
      (.set (.-position sun) (* span 0.5) (* span 0.9) (* span 0.35))
      (.add scene sun))
    (.set (.-position camera) (* span 0.75) (* span 0.55) (* span 0.75))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(view-opts []
                {:show (.-show controls)
                 :exaggeration (double (.-exaggeration controls))
                 :cell-scale (/ span (max 1 (dec (long (:size @state)))))
                 :height-scale (:height-scale config)})
              (rain-opts []
                {:droplets (long (.-droplets controls))
                 :min-slope (double (.-minSlope controls))
                 :capacity (double (.-capacity controls))
                 :erosion (double (.-erosion controls))
                 :deposition (double (.-deposition controls))
                 :inertia (double (.-inertia controls))
                 :radius (long (.-radius controls))
                 :border (double (.-border controls))})
              (repaint! []
                (refresh-geometry! (:geometry @state) (:field @state) (view-opts)))
              (report! []
                (let [{:keys [field droplets]} @state
                      cells (* (long (:dim field)) (long (:dim field)))]
                  (set! (.-rained controls)
                        (str (.toFixed (/ (double droplets) cells) 2) " per cell"))))
              (rebuild! []
                (when-let [^js m (:mesh @state)]
                  (.remove scene m)
                  (.dispose (.-geometry m))
                  (.dispose (.-material m)))
                (let [size (long (:size @state))
                      field (fresh-field (assoc config :size size))
                      geometry (build-geometry (:dim field))
                      mesh (THREE/Mesh. geometry
                                        (THREE/MeshStandardMaterial.
                                         #js {:vertexColors true
                                              :roughness 0.95
                                              :metalness 0.0
                                              :side THREE/DoubleSide}))]
                  (.add scene mesh)
                  (swap! state assoc :field field :geometry geometry :mesh mesh :droplets 0)
                  (repaint!)
                  (report!)))
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
                  (let [t0 (js/performance.now)]
                    (when (.-rain controls)
                      (let [n (long (.-droplets controls))]
                        (swap! state update :field
                               #(erosion/erode! % (assoc (rain-opts) :droplets n)))
                        (swap! state update :droplets + n))
                      (repaint!)
                      (report!))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (set! (.-reset controls) rebuild!)
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (let [gui (GUI. #js {:container container})
              rain-folder (.addFolder gui "rain")
              model (.addFolder gui "model")
              view (.addFolder gui "view")]
          (doto rain-folder
            (.add controls "rain")
            (-> (.add controls "droplets" 100 6000 100))
            (-> (.add controls "rained") (.listen) (.disable))
            (.add controls "reset"))
          (doto model
            ;; Logarithmic, because the interesting range spans three
            ;; decades and the top of it is the blur described in the
            ;; namespace docstring -- worth being able to reach, and not
            ;; worth half the slider.
            (-> (.add controls "minSlope" 0.00005 0.02)
                (.step 0.00005)
                (.name "min-slope"))
            (-> (.add controls "capacity" 0.5 12.0 0.5))
            (-> (.add controls "erosion" 0.0 1.0 0.05))
            (-> (.add controls "deposition" 0.0 1.0 0.05))
            (-> (.add controls "inertia" 0.0 0.6 0.01))
            (-> (.add controls "radius" 1 6 1))
            (-> (.add controls "border" 0 32 1)))
          (doto view
            (-> (.add controls "show" #js ["drainage" "terrain" "change"])
                (.onChange (fn [_] (repaint!))))
            (-> (.add controls "exaggeration" 0.25 2.5 0.05)
                (.onChange (fn [_] (repaint!))))
            (-> (.add controls "size" (clj->js (:sizes config)))
                (.onChange (fn [v]
                             (swap! state assoc :size (long v))
                             (rebuild!)))))
          (.close model))
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
                   (when-let [container (js/document.getElementById "erosion")]
                     (reset! controller (init! container config))))]
    ((:start c))))

(defn stop!
  "Halt the render loop, keeping the terrain as the rain has left it."
  []
  (when-let [c @controller] ((:stop c))))
