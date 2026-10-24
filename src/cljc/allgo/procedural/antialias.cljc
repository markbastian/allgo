(ns allgo.procedural.antialias
  "Antialiasing a procedural texture: answering what a pixel *covers*
  rather than what is at its middle.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, Worley's
  \"Advanced Antialiasing\", with Peachey's filtered patterns from
  chapter 2. A texture is a function of a point, and a pixel is not a
  point: it is a patch of the surface, which on a plane running off to
  the horizon grows from a fraction of a brick near the camera to whole
  fields of them far away. Sampling the function once at the patch's
  middle reports whatever detail happens to be there, and as the camera
  moves that answer flickers. The right answer is the function's average
  over the patch, and there are two ways to get it.

  **Integrate it.** Some patterns can be averaged exactly, in closed
  form. A pulse train's integral is a staircase with ramps, and the
  average over any interval is a difference of two of them, so a stripe
  is box-filtered at the cost of two evaluations (`filtered-pulse-train`).
  A checkerboard is the exclusive or of two of those, and since a box
  filter over a square is the product of two over lines, it is filtered
  exactly too (`filtered-checker`). For a sum of noise octaves the
  integral is not to be had, but the next best thing is: an octave too
  fine for the pixel averages to its mean, zero, so it is faded out as
  its wavelength nears the pixel's width (`clamped-fbm`). That is
  `fractal/octaves-for` done per pixel instead of per surface.

  **Sample it many times.** Anything can be averaged by evaluating it at
  several points across the patch. Where those points go matters more
  than how many there are: on a regular grid (`:grid`) the error is
  aliasing again, only finer -- a pattern near the grid's own frequency
  turns into moiré. Jittering each sample within its cell of the grid
  (`:jittered`) trades that structured error for noise, which the eye
  forgives. And most pixels do not need many samples at all: `adaptive`
  takes five and subdivides only where they disagree, which puts the
  work on the edges and spends next to nothing on flat color.

  Integration is exact and costs one sample; it works only for what can
  be integrated. Sampling works for everything and costs in proportion
  to the quality. A real shader uses the first wherever it can and the
  second for whatever is left.

  ## Footprints

  Every function here takes the pixel's footprint as a width `w` in the
  texture's own units: how far the texture coordinate moves from one
  pixel to the next. A footprint is really a parallelogram, long and thin
  at a grazing angle; taking its longer side, as the demo does, blurs
  across the short one -- the classic overblur that anisotropic filtering
  exists to fix."
  (:require [allgo.math :as am]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as shaping]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Filtering by integration

(defn- pulse-integral
  "The integral from 0 to `x` of a period-1 pulse train that is 1 for the
  first `duty` of each period: `duty` for every whole period, and the
  part of the current one that is on."
  ^double [^double duty ^double x]
  (let [i (math/floor x)
        f (- x i)]
    (+ (* i duty) (min f duty))))

(defn filtered-pulse-train
  "`shaping/pulse-train` with period 1, box-filtered over width `w`: the
  fraction of `[x - w/2, x + w/2]` that is on. Exact; and as `w` grows
  past a period it settles to `duty`, the pattern's average, which is
  exactly the gray a distant field of stripes should be."
  ^double [^double duty ^double x ^double w]
  (if (< w 1e-9)
    (shaping/pulse-train duty 1.0 x)
    (let [h (* 0.5 w)]
      (/ (- (pulse-integral duty (+ x h)) (pulse-integral duty (- x h))) w))))

(defn checker
  "A unit checkerboard: 1 on the squares whose corner coordinates sum
  even, 0 on the rest."
  ^double [^double x ^double y]
  (let [a (shaping/pulse-train 0.5 1.0 x)
        b (shaping/pulse-train 0.5 1.0 y)]
    (+ a b (* -2.0 a b))))

(defn filtered-checker
  "`checker` box-filtered over a `w` square. The board is `a xor b` of
  two stripe patterns, `a + b - 2ab`; a box filter over a square is the
  product of one along x and one along y, so the average of `ab` is the
  product of the averages and the whole is exact."
  ^double [^double x ^double y ^double w]
  (let [a (filtered-pulse-train 0.5 x w)
        b (filtered-pulse-train 0.5 y w)]
    (+ a b (* -2.0 a b))))

(defn clamped-fbm
  "fBm over `basis` that leaves out what a footprint of width `w` cannot
  hold: `(fn [x y z w])`. Octave `i` has frequency `lacunarity^i`, and is
  faded out as `frequency * w` runs from 1/4 to 1/2 -- half a cycle a
  pixel being where it would start to alias. At `w` = 0 it is plain fBm,
  `octaves` deep. Options as for `allgo.procedural.fractal/fbm`: `:H`,
  `:lacunarity`, `:octaves`."
  ([basis] (clamped-fbm basis {}))
  ([basis {:keys [H lacunarity octaves] :or {H 1.0 lacunarity 2.0 octaves 8}}]
   (let [lacunarity (double lacunarity)
         n (long octaves)
         weights (double-array (for [i (range n)] (math/pow lacunarity (- (* (double H) i)))))]
     (fn ^double [^double x ^double y ^double z ^double w]
       (loop [i 0 f 1.0 acc 0.0]
         (if (>= i n)
           acc
           (let [fade (- 1.0 (shaping/smoothstep 0.25 0.5 (* f w)))]
             (if (<= fade 0.0)
               acc
               (recur (inc i) (* f lacunarity)
                      (+ acc (* fade (aget weights i)
                                (noise/sample basis (* x f) (* y f) (* z f)))))))))))))

;; ---------------------------------------------------------------------------
;; Filtering by sampling

(defn- hash01
  "A repeatable number in [0, 1) for integers `i j k`: the fractional
  part of a sine, enough to scatter jitter without a pattern."
  ^double [^long i ^long j ^long k]
  (let [s (* 43758.5453 (math/sin (+ (* 12.9898 i) (* 78.233 j) (* 37.719 k))))]
    (- s (math/floor s))))

(defn- add-into [acc v]
  (if (number? v) (+ (double acc) (double v)) (mapv + acc v)))

(defn- scale-by [v s]
  (if (number? v) (* (double v) s) (mapv #(* % s) v)))

(defn supersample
  "The average of `f` -- `(fn [x y])` to a number or a color vector --
  over the `w` square centered on `x y`, from `n` by `n` samples.
  `pattern` is `:grid`, a sample in the middle of each cell, or
  `:jittered`, a sample anywhere in each, placed by a hash of `seed`
  (pass the pixel's coordinates, so the jitter differs from pixel to
  pixel and is the same every frame)."
  [f x y w n pattern seed]
  (let [n (long n)
        cell (/ (double w) n)
        x0 (- (double x) (* 0.5 (double w)))
        y0 (- (double y) (* 0.5 (double w)))
        [si sj] seed]
    (loop [i 0 j 0 acc nil]
      (cond
        (>= j n) (scale-by acc (/ 1.0 (* n n)))
        (>= i n) (recur 0 (inc j) acc)
        :else
        (let [[dx dy] (if (= pattern :jittered)
                        (let [k (+ i (* n j))]
                          [(hash01 si sj (* 2 k)) (hash01 si sj (inc (* 2 k)))])
                        [0.5 0.5])
              v (f (+ x0 (* cell (+ i dx))) (+ y0 (* cell (+ j dy))))]
          (recur (inc i) j (if acc (add-into acc v) v)))))))

(defn- disagreement
  "How far apart the largest and smallest of `vs` are: a number's range,
  or a color's largest range over its channels."
  ^double [vs]
  (if (number? (first vs))
    (- (double (apply max vs)) (double (apply min vs)))
    (apply max (for [c (range (count (first vs)))]
                 (let [cs (map #(nth % c) vs)] (- (double (apply max cs)) (double (apply min cs))))))))

(defn adaptive
  "The average of `f` over the `w` square at `x y`, sampled where it is
  needed: `[value samples]`.

  Samples the four corners and the middle; if they agree to within
  `threshold`, their average is the answer. If not, the square is split
  in four and each quarter asked the same, down to `depth` levels. The
  corners are shared with the neighbors in a real renderer; here they
  are not, which overstates the count but not the idea: a flat region
  costs five samples, and an edge costs what it needs."
  [f x y w threshold depth]
  (let [h (* 0.5 (double w))
        x (double x) y (double y)
        vs [(f (- x h) (- y h)) (f (+ x h) (- y h)) (f (- x h) (+ y h)) (f (+ x h) (+ y h)) (f x y)]]
    (if (or (<= depth 0) (<= (disagreement vs) threshold))
      [(scale-by (reduce add-into vs) 0.2) 5]
      (let [q (* 0.5 h)
            parts (for [dx [(- q) q] dy [(- q) q]]
                    (adaptive f (+ x dx) (+ y dy) h threshold (dec depth)))]
        [(scale-by (reduce add-into (map first parts)) 0.25)
         (reduce + 5 (map second parts))]))))

;; ---------------------------------------------------------------------------
;; Patterns, filtered both ways

(def patterns
  "What the demo draws, each as `{:point (fn [x y]) :filtered (fn [x y w])}`:
  the texture at a point, and its exact or approximate box filter over a
  footprint `w`. Values in [0, 1]."
  (let [fb (clamped-fbm (noise/gradient-basis {:seed 21}) {:octaves 9 :H 0.8})]
    {"checker" {:point checker :filtered filtered-checker}
     "stripes" {:point (fn [x _] (shaping/pulse-train 0.5 1.0 x))
                :filtered (fn [x _ w] (filtered-pulse-train 0.5 x w))}
     "fbm" {:point (fn [x y] (am/clamp (+ 0.5 (* 0.6 (fb (* 0.5 x) (* 0.5 y) 0.5 0.0))) 0.0 1.0))
            :filtered (fn [x y w] (am/clamp (+ 0.5 (* 0.6 (fb (* 0.5 x) (* 0.5 y) 0.5 (* 0.5 w)))) 0.0 1.0))}
     "tiles" {:point (fn [x y]
                       (* (checker x y)
                          (am/clamp (+ 0.6 (* 0.5 (fb (* 2.0 x) (* 2.0 y) 0.1 0.0))) 0.0 1.0)))
              ;; A product is not filtered by filtering its factors -- but
              ;; the noise is smooth at the scale of a square, so nearly.
              :filtered (fn [x y w]
                          (* (filtered-checker x y w)
                             (am/clamp (+ 0.6 (* 0.5 (fb (* 2.0 x) (* 2.0 y) 0.1 (* 2.0 w)))) 0.0 1.0)))}}))
