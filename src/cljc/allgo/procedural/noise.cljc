(ns allgo.procedural.noise
  "Basis functions: the raw material every procedural fractal is built from.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*. Chapter 2
  and chapter 12 for lattice noise, chapter 20 (Musgrave, \"MojoWorld:
  Building Procedural Planets\") for the idea that organises this
  namespace. MojoWorld's central move is to separate *what is repeated*
  from *how it is repeated*: a basis function is a band-limited random
  function of space, and a fractal is a recipe for summing scaled copies of
  one. Swap the basis and the same recipe gives you a different world --
  the same ridged multifractal over gradient noise is mountains, and over
  a cellular basis is shattered rock.

  So everything here has one shape:

      basis : (fn [x y z] -> double)      roughly in [-1, 1]
      vector basis : (fn [x y z] -> [x y z])

  and `allgo.procedural.fractal` consumes them without caring which one it
  was handed. The constructors take an options map and *return* a basis,
  because a basis is a seed plus a table and you want to build it once and
  call it a million times.

  ## What a lattice basis is

  Give every point of the integer lattice a repeatable random something --
  a value for value noise, a gradient direction for gradient noise -- and
  interpolate between the eight (or sixteen, in 4D) corners of the cell a
  point falls in. Repeatable is the whole trick, and it is done by hashing
  the integer coordinates rather than storing anything: the table here is
  256 long and the world is unbounded.

  Gradient noise is the one to reach for. Value noise puts its extremes at
  the lattice points, so the grid is visible in the output; gradient noise
  is zero at every lattice point and takes its extremes between them, which
  hides the grid far better. The cost is a dot product per corner.

  ## The fourth dimension

  `gradient-basis-4d` is here for the reason MojoWorld is built on 4D
  noise: a 3D world carved out of a 4D function can be *moved* through the
  fourth coordinate, so clouds evolve and a planet's terrain can be
  dialled continuously from one variant to another, without any of the
  popping a reseed would cause."
  (:require [allgo.procedural.shaping :as shaping]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Repeatable randomness
;;
;; Park-Miller, done in doubles rather than in bit twiddling, so that the
;; JVM's 64-bit longs and JavaScript's 32-bit bitwise operators cannot
;; disagree about what a seed means. Every product below stays under 2^53
;; and so is exact on both.

(def ^:private prng-modulus 2147483647.0)
(def ^:private prng-multiplier 16807.0)

(defn advance
  "The next state of the random stream. Park-Miller's minimal standard.

  The modulus is taken by hand rather than with `rem`, and that is worth a
  paragraph because it is the innermost loop of every cellular and
  sparse-convolution sample. A remainder of two large doubles misses the
  integer fast path in a JavaScript engine and costs on the order of a
  hundred nanoseconds; subtracting `m * floor(x / m)` costs a divide and a
  floor, both of which are single instructions, and measured eight times
  faster here.

  It is also exact, which is the part that has to be argued rather than
  measured. `16807 * s` never exceeds 2^53, so the product is an exact
  integer. Its true quotient by `m` is either an integer -- in which case
  the division is exact -- or at least `1/m` away from one, which is four
  hundred times the largest error the division can introduce. So the floor
  is never off by one, and neither platform can disagree with the other
  about what a seed means."
  ^double [^double s]
  (let [x (* prng-multiplier s)
        v (- x (* prng-modulus (math/floor (/ x prng-modulus))))]
    (if (pos? v) v 1.0)))

(defn unit
  "A random stream state as a number in [0, 1)."
  ^double [^double s]
  (/ s prng-modulus))

(defn signed
  "A random stream state as a number in [-1, 1)."
  ^double [^double s]
  (- (* 2.0 (unit s)) 1.0))

(defn- fold
  "Mixes one integer into a stream state.

  Coordinates are taken modulo 65536, so the lattice repeats every 65536
  cells along each axis. That is no worse than the 256-entry permutation
  table underneath the classic noise functions and far past anything a
  scene will reach."
  ^double [^double s ^long v]
  ;; `bit-and` with a mask is `mod` by a power of two, on both platforms
  ;; and for negative coordinates too, and it is a single instruction
  ;; where `mod` is a call.
  (advance (advance (+ s (double (bit-and v 0xFFFF)) 1.0))))

(defn cell-seed
  "A repeatable random stream for the lattice cell `[i j k]` of world `seed`.

  This is what lets a procedural texture be unbounded and still be a
  function: no cell is stored anywhere, and asking twice gives the same
  answer because the answer is derived from the coordinates."
  ^double [^long seed ^long i ^long j ^long k]
  (-> (+ 1.0 (double (mod seed 2147483646)))
      (fold i) (fold j) (fold k)))

(defn poisson-count
  "How many feature points a cell gets, drawn from a Poisson distribution
  with the given `mean` using the single uniform `u`.

  Worley's cellular basis (chapter 4) wants a Poisson process -- points
  scattered with no structure at all, since any regularity in the count
  per cell would show up as regularity in the texture. Inverting the CDF
  costs one uniform and a short loop, and is capped at 32 because the tail
  beyond that has no measurable probability for the means anyone uses."
  ^long [^double mean ^double u]
  (let [p0 (math/exp (- mean))]
    (loop [k 0 p p0 cdf p0]
      (if (or (<= u cdf) (>= k 32))
        k
        (let [p' (* p (/ mean (double (inc k))))]
          (recur (inc k) p' (+ cdf p')))))))

;; ---------------------------------------------------------------------------
;; Lattice tables

(defn permutation
  "A shuffled 0..255, doubled to 512 so lookups never need a wrap.

  Perlin's table, but built by a seeded Fisher-Yates rather than written
  out as a constant, because MojoWorld gives every fractal its own seed
  and a constant table would give every world the same mountains."
  ^ints [seed]
  (let [n 256
        a (int-array (range n))]
    (loop [i (dec n) s (double (cell-seed (long seed) 0 0 0))]
      (when (pos? i)
        (let [j (long (* (unit s) (inc i)))
              t (aget a i)]
          (aset a i (aget a j))
          (aset a j t)
          (recur (dec i) (advance s)))))
    (int-array (concat (seq a) (seq a)))))

(defn value-table
  "256 random values in [-1, 1], the payload of value noise."
  ^doubles [seed]
  (let [a (double-array 256)]
    (loop [i 0 s (advance (double (cell-seed (long seed) 1 1 1)))]
      (when (< i 256)
        (aset a i (signed s))
        (recur (inc i) (advance s))))
    a))

(defn- wrap
  "A lattice coordinate folded into `[0, period)`, or left alone if
  `period` is zero.

  What makes a noise field tileable: fold the lattice and the field
  repeats exactly, seams included, because the cell at `period - 1`
  interpolates towards the cell at 0 -- which is the same cell the tile
  next door starts from."
  ^long [^long v ^long period]
  (if (zero? period) v (mod v period)))

(defn- hash3
  "The permutation table folded over three lattice coordinates."
  ^long [^ints p ^long i ^long j ^long k]
  (aget p (+ (aget p (+ (aget p (bit-and i 0xFF)) (bit-and j 0xFF)))
             (bit-and k 0xFF))))

(defn- hash4
  ;; Five arguments, so the primitive hints have to go: Clojure's
  ;; primitive function interfaces stop at four.
  [^ints p i j k l]
  (long (aget p (+ (aget p (+ (aget p (+ (aget p (bit-and (long i) 0xFF))
                                         (bit-and (long j) 0xFF)))
                              (bit-and (long k) 0xFF)))
                   (bit-and (long l) 0xFF)))))

(defn- lerp
  "Hoisted out of the bases below, which would otherwise allocate one of
  these per sample."
  ^double [^double a ^double b ^double t]
  (+ a (* t (- b a))))

(defn- fade
  "Perlin's quintic fade curve.

  `allgo.procedural.shaping/smootherstep` on [0, 1] with the clamp and the
  division dropped -- the argument here is a cell-local fraction and is
  already in range, and this runs once per axis per sample per octave."
  ^double [^double t]
  (* t t t (+ (* t (- (* t 6.0) 15.0)) 10.0)))

(defn- grad3
  "Perlin's improved gradient set: the twelve midpoints of a cube's edges.

  Twelve rather than the whole sphere, because a dot product with a vector
  whose components are 0 and +/-1 is an add and a negate, and because a
  set this small and this even shows no directional bias."
  ^double [^long h ^double x ^double y ^double z]
  (let [h (bit-and h 15)
        u (if (< h 8) x y)
        v (cond (< h 4) y (or (== h 12) (== h 14)) x :else z)]
    (+ (if (zero? (bit-and h 1)) u (- u))
       (if (zero? (bit-and h 2)) v (- v)))))

(def ^:private grad4-table
  "The 32 gradients of 4D noise: every vector with one zero component and
  +/-1 in the other three, flattened four doubles to a gradient."
  (double-array
   (for [axis (range 4)
         signs (range 8)
         c (range 4)]
     (if (== c axis)
       0.0
       (let [;; The three non-zero slots take the three bits of `signs`.
             slot (if (< c axis) c (dec c))]
         (if (zero? (bit-and signs (bit-shift-left 1 slot))) 1.0 -1.0))))))

(defn- grad4
  [h x y z w]
  (let [b (* 4 (bit-and (long h) 31))
        x (double x) y (double y) z (double z) w (double w)
        ^doubles g grad4-table]
    (+ (* (aget g b) x)
       (* (aget g (+ b 1)) y)
       (* (aget g (+ b 2)) z)
       (* (aget g (+ b 3)) w))))

;; ---------------------------------------------------------------------------
;; The bases

(defn sample
  "Evaluates a basis at a point.

  `(basis x y z)` would do, and on the JVM it would box all three
  arguments and the result on every call -- four allocations per octave
  per sample, which a fractal of eight octaves turns into thirty-two. A
  basis built here is a `^double` function of three `^double`s, so it
  already implements the primitive interface; this reaches for it when it
  is there and falls back when it is not. ClojureScript has no such
  distinction and no such cost, so there it is the plain call."
  #?(:clj ^double [f ^double x ^double y ^double z]
     :cljs [f x y z])
  #?(:clj (if (instance? clojure.lang.IFn$DDDD f)
            (.invokePrim ^clojure.lang.IFn$DDDD f x y z)
            (double (f x y z)))
     :cljs (f x y z)))

(defn value-basis
  "Value noise: a random value at every lattice point, interpolated.

  Cheaper than gradient noise and worse: the extremes sit on the lattice,
  so the texture carries a visible axis-aligned grid. Worth having because
  the difference between the two is the clearest illustration of why
  gradient noise is built the way it is.

  `:period`, if given, makes the field repeat every that many lattice
  cells along each axis. Useful when the result has to be baked into a
  texture that will be tiled: a non-repeating field baked into a repeating
  texture shows a crease at every seam."
  ([] (value-basis {}))
  ([{:keys [seed period] :or {seed 0}}]
   (let [^ints p (permutation seed)
         ^doubles vals (value-table seed)
         period (long (or period 0))
         at (fn ^double [^long i ^long j ^long k]
              (aget vals (hash3 p (wrap i period) (wrap j period) (wrap k period))))]
     (fn ^double [^double x ^double y ^double z]
       (let [i (long (math/floor x)) j (long (math/floor y)) k (long (math/floor z))
             fx (- x i) fy (- y j) fz (- z k)
             u (fade fx) v (fade fy) w (fade fz)]
         (lerp (lerp (lerp (at i j k) (at (inc i) j k) u)
                     (lerp (at i (inc j) k) (at (inc i) (inc j) k) u) v)
               (lerp (lerp (at i j (inc k)) (at (inc i) j (inc k)) u)
                     (lerp (at i (inc j) (inc k)) (at (inc i) (inc j) (inc k)) u) v)
               w))))))

(defn gradient-basis
  "Perlin's improved gradient noise. The default basis everywhere here.

  Zero at every lattice point, with a random gradient direction there; the
  value between points comes from those gradients rather than from stored
  values, which is what keeps the lattice out of the picture. Measured
  extremes land near +/-0.8 rather than exactly +/-1, which matters only
  when choosing the `:offset` of a multifractal.

  `:period`, if given, makes the field repeat every that many lattice
  cells along each axis. Useful when the result has to be baked into a
  texture that will be tiled: a non-repeating field baked into a repeating
  texture shows a crease at every seam."
  ([] (gradient-basis {}))
  ([{:keys [seed period] :or {seed 0}}]
   (let [^ints p (permutation seed)
         period (long (or period 0))]
     (fn ^double [^double x ^double y ^double z]
       (let [i (long (math/floor x)) j (long (math/floor y)) k (long (math/floor z))
             fx (- x i) fy (- y j) fz (- z k)
             u (fade fx) v (fade fy) w (fade fz)
             g (fn ^double [^long di ^long dj ^long dk]
                 (grad3 (hash3 p
                               (wrap (+ i di) period)
                               (wrap (+ j dj) period)
                               (wrap (+ k dk) period))
                        (- fx di) (- fy dj) (- fz dk)))]
         (lerp (lerp (lerp (g 0 0 0) (g 1 0 0) u)
                     (lerp (g 0 1 0) (g 1 1 0) u) v)
               (lerp (lerp (g 0 0 1) (g 1 0 1) u)
                     (lerp (g 0 1 1) (g 1 1 1) u) v)
               w))))))

(defn gradient-basis-4d
  "Gradient noise in four dimensions: `(fn [x y z w] -> double)`.

  Sixteen corners instead of eight. The fourth coordinate is what lets a
  world be animated or morphed continuously -- hold `w` fixed and this is
  an ordinary 3D basis, sweep it and the whole 3D function deforms
  smoothly into an unrelated one. `slice` turns one back into a 3D basis."
  ([] (gradient-basis-4d {}))
  ([{:keys [seed] :or {seed 0}}]
   (let [^ints p (permutation seed)]
     (fn ^double [^double x ^double y ^double z ^double w]
       (let [i (long (math/floor x)) j (long (math/floor y))
             k (long (math/floor z)) l (long (math/floor w))
             fx (- x i) fy (- y j) fz (- z k) fw (- w l)
             su (fade fx) sv (fade fy) sw (fade fz) st (fade fw)
             g (fn ^double [^long di ^long dj ^long dk ^long dl]
                 (grad4 (hash4 p (+ i di) (+ j dj) (+ k dk) (+ l dl))
                        (- fx di) (- fy dj) (- fz dk) (- fw dl)))
             face (fn ^double [^long dl]
                    (lerp (lerp (lerp (g 0 0 0 dl) (g 1 0 0 dl) su)
                                (lerp (g 0 1 0 dl) (g 1 1 0 dl) su) sv)
                          (lerp (lerp (g 0 0 1 dl) (g 1 0 1 dl) su)
                                (lerp (g 0 1 1 dl) (g 1 1 1 dl) su) sv)
                          sw))]
         (lerp (face 0) (face 1) st))))))

(defn slice
  "A 3D basis that is the 4D basis `basis4` held at `w`."
  [basis4 w]
  (let [w (double w)]
    (fn ^double [^double x ^double y ^double z] (basis4 x y z w))))

(defn vector-basis
  "A random *direction* at every point: `(fn [x y z] -> [x y z])`.

  The book's `VecNoise3`, and the input to every kind of domain
  distortion. Three gradient bases would do it; three widely separated
  reads of one basis do it as well and build one table instead of three.
  The offsets are large and irrational-looking on purpose, so that the
  three components never line up into a visible correlation."
  ([] (vector-basis {}))
  ([{:keys [seed] :or {seed 0}}]
   (let [n (gradient-basis {:seed seed})]
     (fn [^double x ^double y ^double z]
       [(n x y z)
        (n (+ y 31.416) (+ z 47.853) (+ x 19.271))
        (n (+ z 73.129) (+ x 11.937) (+ y 59.642))]))))

(defn sparse-convolution-basis
  "Lewis's sparse convolution noise: random impulses convolved with a kernel.

  Not a lattice basis at all. Scatter a Poisson process of impulses through
  space, each with a random weight, and sum a smooth radial kernel around
  each one. The lattice survives only as bookkeeping -- impulses are drawn
  per cell so that the sum stays finite and repeatable -- and none of it
  shows in the output, so this is the basis with no grid whatsoever. It
  costs 27 cells' worth of impulses per sample, which is why it is not the
  default.

  `:density` is the mean impulses per unit cell, `:radius` the kernel's
  support measured in cells."
  ([] (sparse-convolution-basis {}))
  ([{:keys [seed density radius] :or {seed 0 density 3.0 radius 1.0}}]
   (let [seed (long seed)
         density (double density)
         radius (double radius)
         r2 (* radius radius)
         ;; The neighbourhood searched, flattened to one index so that the
         ;; sample loop is a loop and not three of them.
         reach (long (math/ceil radius))
         side (inc (* 2 reach))
         cells (* side side side)
         ;; The sum of `n` independent impulses has variance proportional
         ;; to `n`, and `n` is proportional to the density and to the
         ;; volume of the kernel's support, so dividing by the root of both
         ;; holds the spread steady however they are set. The constant is
         ;; measured, not derived: it puts the extremes near +/-1, which is
         ;; where the other bases sit and where a multifractal's `:offset`
         ;; expects to find them.
         norm (/ 0.55 (math/sqrt (* density radius radius radius)))]
     (fn ^double [^double x ^double y ^double z]
       (let [ci (long (math/floor x)) cj (long (math/floor y)) ck (long (math/floor z))]
         (loop [c 0 sum 0.0]
           (if (= c cells)
             (* norm sum)
             (let [i (+ ci (- (rem c side) reach))
                   j (+ cj (- (rem (quot c side) side) reach))
                   k (+ ck (- (quot c (* side side)) reach))
                   s0 (cell-seed seed i j k)
                   m (poisson-count density (unit s0))]
               (recur
                (inc c)
                ;; Boxing the running sum would cost more than the impulses
                ;; do, so the inner loop is pinned back to a double.
                (double
                 (loop [n 0 s (advance s0) acc sum]
                   (if (>= n m)
                     acc
                     (let [px (+ i (unit s))
                           s (advance s)
                           py (+ j (unit s))
                           s (advance s)
                           pz (+ k (unit s))
                           s (advance s)
                           w (signed s)
                           dx (- x px) dy (- y py) dz (- z pz)
                           d2 (+ (* dx dx) (* dy dy) (* dz dz))]
                       (recur (inc n)
                              (advance s)
                              (if (< d2 r2)
                               ;; 1 - 3t^2 + 2t^3: one at the impulse, and
                               ;; zero with zero slope at the edge of its
                               ;; support, so neighbouring cells join with
                               ;; no seam.
                                (let [t (math/sqrt (/ d2 r2))]
                                  (+ acc (* w (+ 1.0 (* t t (- (* 2.0 t) 3.0))))))
                                acc)))))))))))))))

;; ---------------------------------------------------------------------------
;; Transforming a basis into another basis
;;
;; Everything below takes a basis and returns a basis, so they compose with
;; each other and with the fractals in `allgo.procedural.fractal` in any
;; order. This is the whole of MojoWorld's "function space": a small set of
;; bases, a small set of fractals, and free composition between them.

(defn absolute
  "`|basis|`. The fold that turns noise into Perlin's turbulence.

  Reflecting the negative half upwards leaves a crease wherever the
  original crossed zero, and those creases are the thin dark filaments
  that read as smoke, flame and marble veins."
  [basis]
  (fn ^double [^double x ^double y ^double z] (abs (sample basis x y z))))

(defn ridge
  "`(offset - |basis|)^2`. The fold that turns noise into ridge lines.

  Same crease as `absolute`, then turned upside down so it becomes a
  ridge rather than a valley, then squared to sharpen it. Musgrave's
  ridged multifractal is this applied octave by octave; applied once it is
  a basis in its own right."
  ([basis] (ridge basis 1.0))
  ([basis offset]
   (let [offset (double offset)]
     (fn ^double [^double x ^double y ^double z]
       (let [s (- offset (abs (sample basis x y z)))]
         (* s s))))))

(defn scaled
  "`basis` with its domain multiplied by `f` and its value by `a`."
  ([basis f] (scaled basis f 1.0))
  ([basis f a]
   (let [f (double f) a (double a)]
     (fn ^double [^double x ^double y ^double z]
       (* a (sample basis (* f x) (* f y) (* f z)))))))

(defn translated
  "`basis` with its domain shifted by `[dx dy dz]`.

  The cheapest way to decorrelate two uses of one basis."
  [basis [dx dy dz]]
  (let [dx (double dx) dy (double dy) dz (double dz)]
    (fn ^double [^double x ^double y ^double z]
      (basis (+ x dx) (+ y dy) (+ z dz)))))

(defn distorted
  "Musgrave's variable-lacunarity noise: `basis` read at a distorted point.

  Chapter 15. Look up a random vector at the point, walk `amount` of the
  way along it, and evaluate the basis *there*. The result is the same
  basis with its frequency varying from place to place -- stretched here,
  compressed there -- which is the single cheapest way to make a texture
  stop looking machine-made. It is also the most expensive knob in the
  book: every sample now costs four noise evaluations instead of one.

  `warp` is a vector basis; the domain is misregistered by half a cell
  first so that the distortion and the basis it distorts do not share
  their zeroes."
  ([basis warp] (distorted basis warp 1.0))
  ([basis warp amount]
   (let [amount (double amount)]
     (fn ^double [^double x ^double y ^double z]
       (let [[dx dy dz] (warp (+ x 0.5) (+ y 0.5) (+ z 0.5))]
         (basis (+ x (* amount (double dx)))
                (+ y (* amount (double dy)))
                (+ z (* amount (double dz)))))))))

(defn stepped
  "`basis` quantised into `n` terraces, with `smooth` of each riser blended.

  MojoWorld's stepped bases, and the direct route to sedimentary strata
  and rice-terrace landscapes: the terrain function is unchanged, but its
  range is stepped, so contour lines become cliffs."
  ([basis n] (stepped basis n 0.0))
  ([basis n smooth]
   (let [n (double n)
         smooth (double smooth)]
     (fn ^double [^double x ^double y ^double z]
       (let [v (* 0.5 (+ 1.0 (sample basis x y z)))
             t (* v n)
             i (math/floor t)
             f (- t i)
             blend (if (pos? smooth) (shaping/smoothstep 0.0 smooth f) 0.0)]
         (- (* 2.0 (/ (+ i blend) n)) 1.0))))))

(defn shaped
  "`basis` with `f` applied to its value, mapped through [0, 1] and back.

  The bridge to `allgo.procedural.shaping`: `(shaped b #(shaping/gain 0.7 %))`
  is a basis with more contrast, `(shaped b #(shaping/bias 0.3 %))` one
  that spends more of its time low. Bases run [-1, 1] and the shaping
  functions want [0, 1], and doing that conversion here is the reason this
  exists rather than an inline `comp`."
  [basis f]
  (fn ^double [^double x ^double y ^double z]
    (let [v (* 0.5 (+ 1.0 (sample basis x y z)))]
      (- (* 2.0 (double (f v))) 1.0))))
