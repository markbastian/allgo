(ns allgo.procedural.fractal
  "Fractal constructions: how to build a world out of one basis function.

  Musgrave, in Ebert et al., *Texturing & Modeling: A Procedural
  Approach* -- chapters 14 (\"A Brief Introduction to Fractals\"), 16
  (\"Procedural Fractal Terrains\") and 20 (\"MojoWorld: Building
  Procedural Planets\"). Every function here takes a basis from
  `allgo.procedural.noise` -- or another fractal from here, which is the
  same thing -- and returns another basis. That closure is the design: a
  fractal *is* a basis, so they nest without any special provision for it.

  ## The one idea

  Add up copies of a basis, each smaller and finer than the last:

      value = sum over i of  basis(point * lacunarity^i) * lacunarity^(-H*i)

  `lacunarity` is the gap between successive frequencies, conventionally
  2 -- each octave is twice as fine as the one before. `H` is the Hurst
  exponent, which decides how fast the amplitudes fall off and is the only
  parameter that changes what *kind* of surface you get. H = 1 gives a
  smooth, strongly correlated surface; H = 0 gives amplitudes that do not
  fall off at all and a surface that is rough at every scale, which is
  white-ish noise. The fractal dimension is 3 - H for a surface in space.

  ## Monofractal against multifractal

  `fbm` and `turbulence` are *monofractal*: one H, so the roughness is the
  same everywhere. Real terrain is not like that -- plains are smooth,
  mountains are rough, and the roughness varies with where you are and
  especially with altitude. The rest of this namespace is Musgrave's
  answer: make each octave's contribution depend on what the sum so far
  came to. `hetero-terrain` scales the next octave by the current
  altitude, so low ground stays flat. `hybrid-multifractal` and
  `ridged-multifractal` weight it by the previous octave's local value, so
  detail appears only where there is already something to put it on. That
  one change is most of the difference between \"fractal noise\" and
  \"landscape\".

  ## Ranges

  The monofractals are roughly zero-centered and roughly as wide as their
  basis. The multifractals are not: they multiply, so their range depends
  on every parameter at once and cannot usefully be predicted. That is not
  a flaw to be fixed, it is what makes them do what they do -- but it does
  mean you have to measure before you can map one onto a color or a
  height. `value-range` is here for that."
  (:require [allgo.procedural.noise :as noise]
            [clojure.math :as math]))

(def defaults
  "The book's parameters, under the book's names.

  `:lacunarity` 2.0 is the conventional choice and a slightly unfortunate
  one: an exact doubling lands every octave on the same lattice planes, so
  their features line up and leave a faint regularity in the result.
  Musgrave suggests nudging it off the integer -- 1.93, 2.137 -- when that
  shows. It is left at 2.0 here because that is what the book's listings
  use and what every published picture was made with."
  {:H 1.0 :lacunarity 2.0 :octaves 8.0 :offset 0.7 :gain 2.0})

(defn spectral-weights
  "The amplitude allowed to each octave: `lacunarity^(-H*i)`.

  The book calls this the exponent array and builds it once per set of
  parameters rather than calling `pow` inside the octave loop, which is
  worth doing when the loop runs once per pixel."
  ^doubles [{:keys [H lacunarity octaves] :or {H (:H defaults)
                                               lacunarity (:lacunarity defaults)
                                               octaves (:octaves defaults)}}]
  (let [n (inc (long (math/ceil (double octaves))))
        a (double-array n)
        H (double H)
        lacunarity (double lacunarity)]
    (dotimes [i n]
      (aset a i (math/pow lacunarity (* (- H) (double i)))))
    a))

(defn octaves-for
  "The most octaves worth summing when the result is sampled every
  `spacing` apart, in the basis's own domain units.

  This is the antialiasing story of chapters 2 and 17, and it is the
  difference between a fractal that looks detailed and one that looks
  noisy. Octave `i` has features about `lacunarity^-i` across. Once those
  are smaller than two samples apart they cannot be reconstructed, and
  what arrives instead is not detail but a different, lower-frequency
  pattern -- speckle that moves when the camera does. Summing them costs
  time to make the picture worse.

      lacunarity^-i >= 2 * spacing   =>   i <= -log(2 * spacing) / log(lacunarity)

  The answer is deliberately fractional. Every construction here takes
  non-integer `:octaves` and fades the last one in by its remainder, so a
  surface sampled more finely gains detail continuously rather than
  gaining a whole octave at once -- which would pop.

  Never less than one octave, however coarse the sampling: a fractal with
  no octaves at all is zero, and a flat answer is a worse lie than a
  coarse one."
  ^double [spacing lacunarity]
  (let [spacing (double spacing)
        lacunarity (double lacunarity)]
    (if (or (<= spacing 0.0) (<= lacunarity 1.0))
      ##Inf
      (max 1.0 (/ (- (math/log (* 2.0 spacing))) (math/log lacunarity))))))

(defn- params
  "Splits an options map into what every octave loop below needs."
  [opts]
  (let [{:keys [lacunarity octaves offset gain]} (merge defaults opts)
        octaves (double octaves)
        whole (long octaves)]
    {:weights (spectral-weights (merge defaults opts))
     :lacunarity (double lacunarity)
     :octaves octaves
     :whole whole
     :remainder (- octaves whole)
     :offset (double offset)
     :gain (double gain)}))

(defn fbm
  "Fractional Brownian motion: the plain sum of scaled copies of `basis`.

  The monofractal, and the thing every other construction here is a
  deviation from. Statistically identical everywhere, which is exactly
  what makes it read as texture rather than as terrain -- but it is the
  right answer for clouds, for water, and for any surface that has no
  reason to prefer one altitude to another.

  Non-integer `:octaves` are allowed and are how you fade detail in
  smoothly as a camera approaches: the last octave is added at partial
  strength rather than appearing all at once."
  ([basis] (fbm basis {}))
  ([basis opts]
   (let [{:keys [^doubles weights ^double lacunarity ^long whole ^double remainder]} (params opts)]
     (fn ^double [^double x ^double y ^double z]
       (loop [i 0 x x y y z z acc 0.0]
         (if (< i whole)
           (recur (inc i)
                  (* x lacunarity) (* y lacunarity) (* z lacunarity)
                  (+ acc (* (noise/sample basis x y z) (aget weights i))))
           (if (pos? remainder)
             (+ acc (* remainder (noise/sample basis x y z) (aget weights whole)))
             acc)))))))

(defn turbulence
  "Perlin's turbulence: `fbm` over `|basis|`.

  Folding the basis at zero before summing leaves a crease along every
  zero crossing of every octave, and those creases are what make
  turbulence look like smoke and flame where `fbm` looks like haze. Note
  that it is no longer zero-centered -- every term is positive.

  Written as the composition it is, rather than as its own loop, because
  that is the point of `allgo.procedural.noise/absolute` existing."
  ([basis] (turbulence basis {}))
  ([basis opts] (fbm (noise/absolute basis) opts)))

(defn multifractal
  "Musgrave's first multifractal: the octaves multiplied rather than added.

  `value *= offset + basis(point) * weight`. Because each octave scales
  everything below it, the amplitude of the detail follows the value of
  the sum so far -- high ground gets rough and low ground stays smooth,
  with no explicit rule saying so. `:offset` keeps the running product
  away from zero, where it would collapse and never recover."
  ([basis] (multifractal basis {}))
  ([basis opts]
   (let [{:keys [^doubles weights ^double lacunarity ^long whole
                 ^double remainder ^double offset]} (params opts)]
     (fn ^double [^double x ^double y ^double z]
       (loop [i 0 x x y y z z acc 1.0]
         (if (< i whole)
           (recur (inc i)
                  (* x lacunarity) (* y lacunarity) (* z lacunarity)
                  (* acc (+ offset (* (noise/sample basis x y z) (aget weights i)))))
           (if (pos? remainder)
             (* acc (* remainder (+ offset (* (noise/sample basis x y z) (aget weights whole)))))
             acc)))))))

(defn hetero-terrain
  "Heterogeneous terrain: octaves scaled by the altitude reached so far.

  The first octave is taken unscaled and becomes a rough altitude; every
  octave after it is multiplied by that running altitude before being
  added. Valleys therefore get almost no detail and peaks get all of it,
  which is the single most recognizable property of real landscape and
  the reason this is the one to reach for first when you want ground.

  `:offset` sets the altitude of the sea floor, in effect: it is added to
  every octave, so raising it lifts the whole function away from zero and
  lets detail survive in more places."
  ([basis] (hetero-terrain basis {}))
  ([basis opts]
   (let [{:keys [^doubles weights ^double lacunarity ^long whole
                 ^double remainder ^double offset]} (params opts)]
     (fn ^double [^double x ^double y ^double z]
       (let [value (+ offset (noise/sample basis x y z))]
         (loop [i 1
                x (* x lacunarity) y (* y lacunarity) z (* z lacunarity)
                acc value]
           (if (< i whole)
             (let [increment (* (+ (noise/sample basis x y z) offset)
                                (aget weights i)
                                acc)]
               (recur (inc i)
                      (* x lacunarity) (* y lacunarity) (* z lacunarity)
                      (+ acc increment)))
             (if (pos? remainder)
               (+ acc (* remainder
                         (* (+ (noise/sample basis x y z) offset) (aget weights whole))
                         acc))
               acc))))))))

(defn hybrid-multifractal
  "The hybrid: octaves weighted by the *previous octave's* local value.

  Where `hetero-terrain` weights by the running altitude, this weights by
  the last octave alone, and lets that weight decay multiplicatively. The
  effect is terrain whose roughness varies over the ground rather than
  strictly with height -- smooth plains beside rough hills at the same
  altitude, which `hetero-terrain` cannot produce. Musgrave's own
  favorite for continents with mountain ranges on them.

  The weight is clamped at 1 each octave because without it the product
  diverges: a run of large values would feed on itself."
  ([basis] (hybrid-multifractal basis {}))
  ([basis opts]
   (let [{:keys [^doubles weights ^double lacunarity ^long whole
                 ^double remainder ^double offset]} (params opts)]
     (fn ^double [^double x ^double y ^double z]
       (let [first-signal (* (+ (noise/sample basis x y z) offset) (aget weights 0))]
         (loop [i 1
                x (* x lacunarity) y (* y lacunarity) z (* z lacunarity)
                acc first-signal
                weight first-signal]
           (if (< i whole)
             (let [w (min 1.0 weight)
                   signal (* (+ (noise/sample basis x y z) offset) (aget weights i))]
               (recur (inc i)
                      (* x lacunarity) (* y lacunarity) (* z lacunarity)
                      (+ acc (* w signal))
                      (* weight signal)))
             (if (pos? remainder)
               (+ acc (* remainder (noise/sample basis x y z) (aget weights whole)))
               acc))))))))

(defn ridged-multifractal
  "Ridged multifractal: the one that makes mountain ranges.

  Two things happen per octave. The basis is folded at zero, turned upside
  down and squared -- `(offset - |basis|)^2` -- which converts its zero
  crossings into sharp ridge lines, since a crease pointing up reads as a
  watershed. Then each octave is weighted by the previous octave's signal
  scaled by `:gain`, so that fine ridges only appear along the flanks of
  coarse ones. Together those give branching ridge networks that look
  eroded without any erosion having been simulated.

  `:offset` should be near 1 -- it is what the fold is subtracted from,
  so it sets where the ridges sit. `:gain` near 2 is the book's value."
  ([basis] (ridged-multifractal basis {}))
  ([basis opts]
   ;; The book's listings use 0.7 for the offset of the two terrain
   ;; multifractals and 1.0 for this one, where the offset is what the
   ;; folded basis is subtracted from and so decides where the ridges sit.
   (let [{:keys [^doubles weights ^double lacunarity ^long whole
                 ^double offset ^double gain]} (params (merge {:offset 1.0} opts))
         fold (fn ^double [^double v]
                (let [s (- offset (abs v))] (* s s)))]
     (fn ^double [^double x ^double y ^double z]
       (let [signal (fold (noise/sample basis x y z))]
         (loop [i 1
                x (* x lacunarity) y (* y lacunarity) z (* z lacunarity)
                acc signal
                prev signal]
           (if (>= i whole)
             acc
             (let [w (min 1.0 (max 0.0 (* prev gain)))
                   s (* (fold (noise/sample basis x y z)) w)]
               (recur (inc i)
                      (* x lacunarity) (* y lacunarity) (* z lacunarity)
                      (+ acc (* s (aget weights i)))
                      s)))))))))

(def constructions
  "The fractal constructions by name, so a demo or a parameter file can
  name one without knowing which var it is."
  {:fbm fbm
   :turbulence turbulence
   :multifractal multifractal
   :hetero-terrain hetero-terrain
   :hybrid-multifractal hybrid-multifractal
   :ridged-multifractal ridged-multifractal})

(defn build
  "Builds `[construction basis opts]` from data. `(build :fbm b {})`."
  [construction basis opts]
  ((get constructions construction fbm) basis opts))

;; ---------------------------------------------------------------------------
;; Working with the result

(defn value-range
  "`[lowest highest]` of `f` over `n` samples of the box `[lo hi]^3`.

  The multifractals have no range you can derive, so measuring is the
  only honest way to map one onto heights or colors. Samples are taken
  on a deterministic low-discrepancy walk rather than at random, so that
  the answer does not wobble between calls and a parameter change that
  moves the range shows up as a change in the range."
  ([f n] (value-range f n -4.0 4.0))
  ([f n lo hi]
   (let [n (long n)
         lo (double lo) hi (double hi) span (- hi lo)
         ;; Three mutually irrational strides: the additive recurrence
         ;; covers the box evenly without ever repeating.
         wrap (fn ^double [^double t] (+ lo (* span (- t (math/floor t)))))]
     (loop [i 0 mn ##Inf mx ##-Inf]
       (if (= i n)
         [mn mx]
         (let [t (double i)
               v (noise/sample f
                               (wrap (* t 0.7548776662))
                               (wrap (* t 0.5698402909))
                               (wrap (* t 0.3819660112)))]
           (recur (inc i) (min mn v) (max mx v))))))))

(defn normalized
  "`f` rescaled so that `[lo hi]` becomes [0, 1], and clamped there.

  The other half of `value-range`, and the last step before a fractal
  becomes a color."
  [f lo hi]
  (let [lo (double lo)
        span (let [s (- (double hi) lo)] (if (zero? s) 1.0 s))]
    (fn ^double [^double x ^double y ^double z]
      (let [t (/ (- (noise/sample f x y z) lo) span)]
        (cond (< t 0.0) 0.0 (> t 1.0) 1.0 :else t)))))

(defn gradient
  "The numerical gradient of `f` at a point, by central differences.

  `[dx dy dz]`. Six evaluations, which for an eight-octave multifractal is
  not cheap -- but a procedural function has no analytic derivative to
  offer, and a surface normal has to come from somewhere."
  [f [x y z] epsilon]
  (let [x (double x) y (double y) z (double z) e (double epsilon)
        d (/ 1.0 (* 2.0 e))]
    [(* d (- (noise/sample f (+ x e) y z) (noise/sample f (- x e) y z)))
     (* d (- (noise/sample f x (+ y e) z) (noise/sample f x (- y e) z)))
     (* d (- (noise/sample f x y (+ z e)) (noise/sample f x y (- z e))))]))
