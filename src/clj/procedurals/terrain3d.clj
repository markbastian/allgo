(ns procedurals.terrain3d
  (:require [procedurals.perlin :as perlin]
            [procedurals.terrain :as terrain]
            [quil.applet :as applet]
            [quil.core :as q]
            [quil.middleware :as m])
  (:import [processing.core PApplet PConstants PShape]))

(defn smoothstep [edge0 edge1 x]
  (let [t (-> (/ (double (- x edge0)) (- edge1 edge0)) (max 0.0) (min 1.0))]
    (* t t (- 3.0 (* 2.0 t)))))

(defn mix3 [[r1 g1 b1] [r2 g2 b2] t]
  [(+ r1 (* t (- r2 r1))) (+ g1 (* t (- g2 g1))) (+ b1 (* t (- b2 b1)))])

(def water-color [35 60 55])
(def grass-color-a [70 105 48])
(def grass-color-b [92 120 58])
(def dirt-color [118 92 58])
(def rock-color-a [108 102 96])
(def rock-color-b [90 86 82])
(def snow-color [232 234 238])

(defn terrain-color
  "Blends a natural-looking color from three independent signals: `height-t`
  (normalized elevation) picks the base biome band, `slope-t` (normalized
  steepness) pulls steep faces toward bare rock regardless of elevation so
  cliffs and peaks don't get an unbroken snowcap, and `noise-t` (fractal
  Perlin noise sampled in world space) warps the band boundaries and mixes
  in a second palette color so regions read as patchy rather than banded."
  [height-t slope-t noise-t]
  (let [t     (-> height-t (+ (* 0.12 (- noise-t 0.5))) (max 0.0) (min 1.0))
        grass (mix3 grass-color-a grass-color-b noise-t)
        rock  (mix3 rock-color-a rock-color-b noise-t)
        base  (cond
                (< t 0.12) (mix3 water-color grass (smoothstep 0.02 0.12 t))
                (< t 0.55) (mix3 grass dirt-color (smoothstep 0.35 0.55 t))
                (< t 0.78) (mix3 dirt-color rock (smoothstep 0.55 0.78 t))
                :else (mix3 rock snow-color (smoothstep 0.78 0.9 t)))
        rockiness (smoothstep 0.18 0.55 slope-t)
        [r g b] (mix3 base rock rockiness)]
    (q/color (max 0 (min 255 r)) (max 0 (min 255 g)) (max 0 (min 255 b)))))

(defn slope-at
  "Gradient magnitude (world-space rise/run) of the heightmap at grid cell
  [i j], via central differences. Clamped to the grid at the edges."
  [grid dim i j cell-scale height-scale]
  (let [clamp (fn [v] (max 0 (min (dec dim) v)))
        h     (fn [ii jj] (get-in grid [(clamp ii) (clamp jj)]))
        dzdx  (/ (* height-scale (- (h i (inc j)) (h i (dec j)))) (* 2.0 cell-scale))
        dzdz  (/ (* height-scale (- (h (inc i) j) (h (dec i) j))) (* 2.0 cell-scale))]
    (Math/sqrt (+ (* dzdx dzdx) (* dzdz dzdz)))))

(defn noise-at
  "Low-frequency fractal (multi-octave) Perlin noise sampled at a world
  position, normalized to 0..1, used to break up the height/slope bands."
  [wx wz noise-scale]
  (perlin/operlin (/ wx noise-scale) (/ wz noise-scale) 0.0 0.5 4))

(defn build-mesh-shape
  "Bakes the heightmap into a single retained PShape (GPU-buffered) so `draw`
  can render the whole terrain with one `q/shape` call instead of re-issuing
  tens of thousands of vertex/fill calls every frame."
  [grid dim lo hi cell-scale height-scale noise-scale]
  (let [offset     (* 0.5 (dec dim) cell-scale)
        slopes     (for [i (range dim) j (range dim)] (slope-at grid dim i j cell-scale height-scale))
        max-slope  (max 1e-6 (apply max slopes))
        ^PShape sh (.createShape ^PApplet (applet/current-applet))]
    (.beginShape sh PConstants/QUADS)
    (.noStroke sh)
    (doseq [i (range (dec dim))
            j (range (dec dim))]
      (doseq [[ii jj] [[i j] [i (inc j)] [(inc i) (inc j)] [(inc i) j]]]
        (let [h        (get-in grid [ii jj])
              height-t (if (= lo hi) 0.5 (/ (double (- h lo)) (- hi lo)))
              slope-t  (/ (slope-at grid dim ii jj cell-scale height-scale) max-slope)
              wx       (- (* jj cell-scale) offset)
              wz       (- (* ii cell-scale) offset)
              noise-t  (noise-at wx wz noise-scale)]
          (.fill sh (unchecked-int (terrain-color height-t slope-t noise-t)))
          (.vertex sh (float wx) (float (- (* h height-scale))) (float wz)))))
    (.endShape sh)
    sh))

(defn build-terrain [{:keys [iterations width height-scale cell-scale noise-scale] :or {noise-scale 180.0}}]
  (let [grid-data (terrain/generate {:width width :iterations iterations
                                     :corners [0.0 (rand) (rand) (rand)]})
        grid       (terrain/cells->grid grid-data)
        dim        (:dim grid-data)
        hs         (terrain/heights grid-data)
        lo         (apply min hs)
        hi         (apply max hs)]
    {:dim dim
     :terrain-shape (build-mesh-shape grid dim lo hi cell-scale height-scale noise-scale)}))

(defn setup [config]
  (q/frame-rate 30)
  (merge (build-terrain config)
         {:config config :rot-x 0.55 :rot-y 0.8 :zoom 1.0}))

(defn regenerate [state]
  (merge state (build-terrain (:config state))))

(defn key-pressed [state {:keys [key]}]
  (if (= key :r) (regenerate state) state))

(def max-elevation (- (/ Math/PI 2) 0.1))

(defn mouse-dragged [state {:keys [x y p-x p-y]}]
  (-> state
      (update :rot-y - (* 0.01 (- x p-x)))
      (update :rot-x (fn [rx] (-> rx (+ (* 0.01 (- y p-y)))
                                  (max (- max-elevation))
                                  (min max-elevation))))))

(defn mouse-wheel [state rotation]
  (update state :zoom #(-> % (- (* rotation 0.05)) (max 0.3) (min 4.0))))

(defn draw [{:keys [dim terrain-shape rot-x rot-y zoom]
             {:keys [cell-scale]} :config}]
  (q/background 15 15 25)
  (let [span   (* cell-scale (dec dim))
        radius (/ (* span 1.4) zoom)
        ex     (* radius (Math/cos rot-x) (Math/sin rot-y))
        ey     (* (- radius) (Math/sin rot-x))
        ez     (* radius (Math/cos rot-x) (Math/cos rot-y))]
    (q/camera ex ey ez 0 0 0 0 1 0)
    (q/perspective (/ Math/PI 3) (/ (q/width) (q/height)) 1 (* radius 10))
    (q/directional-light 220 220 200 -0.4 -1 -0.3)
    (q/ambient-light 70 70 80)
    (q/point-light 255 255 240 ex (- radius) ez)
    (q/shape terrain-shape)))

(defn launch-sketch [config]
  (q/sketch
   :title "Fractal Terrain"
   :setup #(setup config)
   :draw #'draw
   :key-pressed #'key-pressed
   :mouse-dragged #'mouse-dragged
   :mouse-wheel #'mouse-wheel
   :size [(:width-px config 900) (:height-px config 700)]
   :renderer :p3d
   :middleware [m/fun-mode]))

(comment
  (launch-sketch
   {:iterations   7
    :width        1.0
    :height-scale 200
    :cell-scale   8
    :width-px     900
    :height-px    700}))
