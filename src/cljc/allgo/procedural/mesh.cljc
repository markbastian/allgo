(ns allgo.procedural.mesh
  (:require [allgo.procedural.perlin :as perlin]
            #?(:clj [clojure.math :as math] :cljs [cljs.math :as math])))

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

(defn terrain-rgb
  "Blends a natural-looking [r g b] color from three independent signals:
  `height-t` (normalized elevation) picks the base biome band, `slope-t`
  (normalized steepness) pulls steep faces toward bare rock regardless of
  elevation so cliffs and peaks don't get an unbroken snowcap, and `noise-t`
  (fractal Perlin noise sampled in world space) warps the band boundaries
  and mixes in a second palette color so regions read as patchy rather than
  banded."
  [height-t slope-t noise-t]
  (let [t         (-> height-t (+ (* 0.12 (- noise-t 0.5))) (max 0.0) (min 1.0))
        grass     (mix3 grass-color-a grass-color-b noise-t)
        rock      (mix3 rock-color-a rock-color-b noise-t)
        base      (cond
                    (< t 0.12) (mix3 water-color grass (smoothstep 0.02 0.12 t))
                    (< t 0.55) (mix3 grass dirt-color (smoothstep 0.35 0.55 t))
                    (< t 0.78) (mix3 dirt-color rock (smoothstep 0.55 0.78 t))
                    :else (mix3 rock snow-color (smoothstep 0.78 0.9 t)))
        rockiness (smoothstep 0.18 0.55 slope-t)
        [r g b]   (mix3 base rock rockiness)]
    [(max 0 (min 255 r)) (max 0 (min 255 g)) (max 0 (min 255 b))]))

(defn noise-at
  "Low-frequency fractal (multi-octave) Perlin noise sampled at a world
  position, normalized to 0..1, used to break up the height/slope bands."
  [wx wz noise-scale]
  (perlin/operlin (/ wx noise-scale) (/ wz noise-scale) 0.0 0.5 4))

(defn gradient-at
  "Central-difference gradient [dzdx dzdz] (world-space rise/run) of the
  heightmap at grid cell [i j]. Clamped to the grid at the edges."
  [grid dim i j cell-scale height-scale]
  (let [clamp (fn [v] (max 0 (min (dec dim) v)))
        h     (fn [ii jj] (get-in grid [(clamp ii) (clamp jj)]))]
    [(/ (* height-scale (- (h i (inc j)) (h i (dec j)))) (* 2.0 cell-scale))
     (/ (* height-scale (- (h (inc i) j) (h (dec i) j))) (* 2.0 cell-scale))]))

(defn slope-at
  "Gradient magnitude (world-space rise/run) of the heightmap at grid cell
  [i j], via central differences. Clamped to the grid at the edges."
  [grid dim i j cell-scale height-scale]
  (let [[dzdx dzdz] (gradient-at grid dim i j cell-scale height-scale)]
    (math/sqrt (+ (* dzdx dzdx) (* dzdz dzdz)))))

(defn normal-at
  "Unit surface normal of the heightmap at grid cell [i j], derived from the
  same central-difference gradient as `slope-at`. `y-sign` is the sign a
  renderer applies to height when building vertex y (+1.0 for a standard
  y-up world, -1.0 if a renderer negates height like `terrain-shape` does);
  a flat cell's normal is [0 y-sign 0]."
  [grid dim i j cell-scale height-scale y-sign]
  (let [[dzdx dzdz] (gradient-at grid dim i j cell-scale height-scale)
        n           [(- dzdx) y-sign (- dzdz)]
        len         (math/sqrt (reduce + (map * n n)))]
    (mapv #(/ % len) n)))
