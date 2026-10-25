(ns allgo.procedural.surflet
  "Noise taken apart: a sum of surflets, one to a lattice point.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, Perlin's
  chapter 12. Gradient noise is usually written as an interpolation --
  dot each corner's gradient with the offset to it, then blend the eight
  results across the cell. Perlin's other reading of the same arithmetic
  is that every lattice point owns a small wave, a *surflet*: its
  gradient's ramp, `g . (p - corner)`, multiplied by a hump that is 1 at
  the point and falls to 0 one cell away in every direction. Noise is
  the sum of all of them. Nothing changes in the numbers -- the blend
  weight a corner gets along an axis is exactly the hump's height there
  -- and `separable` reproduces `noise/gradient-basis` to rounding.

  What changes is what can be done with it:

  - **It is local.** A surflet is zero more than a cell from its point,
    so the noise near a point is the business of the few surflets around
    it and no others. Change one lattice point's gradient and only the
    cells touching it change: `:overrides` does exactly that, which is
    noise you can edit by hand.
  - **The hump can be any shape.** The separable hump is a product of
    one-dimensional fades, square in plan, and the cube lattice shows
    through faintly along the axes. A radial hump, a function of
    distance alone, is round (`radial`) -- but a round hump that reaches
    one cell leaves the middle of each cell nearly unreached, so it must
    reach further and sum more surflets: 64 in 3D at a radius of 1.5,
    against 8. That cost, growing as 4^n with the dimension, is what
    sent Perlin to the simplex lattice (`noise/simplex-basis`), whose
    cells have n + 1 corners rather than 2^n.

  Both kinds are bases in the sense of `allgo.procedural.noise`:
  `(fn [x y z])` to a number near [-1, 1], and so fractals take them."
  (:require [allgo.procedural.noise :as noise]
            [clojure.math :as math]))

(defn hump
  "A surflet's falloff along one axis, or with distance: 1 at 0, falling
  to 0 at `|t|` = 1 and staying there, flat at both ends. One minus
  Perlin's quintic fade, which is why the separable sum is his noise."
  ^double [^double t]
  (let [t (abs t)]
    (if (>= t 1.0)
      0.0
      (- 1.0 (* t t t (+ (* t (- (* t 6.0) 15.0)) 10.0))))))

(defn surflet
  "One surflet, as a function of a point: the gradient `g`'s ramp through
  `center`, times the hump -- separable, the product of a hump along
  each axis, or `:radial`, one hump of distance over `radius`."
  ([g center] (surflet g center {}))
  ([[gx gy gz] [cx cy cz] {:keys [kernel radius] :or {kernel :separable radius 1.5}}]
   (let [gx (double gx) gy (double gy) gz (double gz)
         cx (double cx) cy (double cy) cz (double cz)
         radius (double radius)]
     (fn ^double [^double x ^double y ^double z]
       (let [dx (- x cx) dy (- y cy) dz (- z cz)
             w (if (= kernel :radial)
                 (hump (/ (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz))) radius))
                 (* (hump dx) (hump dy) (hump dz)))]
         (if (zero? w) 0.0 (* w (+ (* gx dx) (* gy dy) (* gz dz)))))))))

(defn gradient-at
  "The gradient of the lattice point `[i j k]` under `opts`: its override
  if it has one, or the one `noise/gradient-basis` would give it."
  [^ints p overrides ijk]
  (or (get overrides ijk) (let [[i j k] ijk] (noise/lattice-gradient p i j k))))

(defn reach
  "How many cells either side of a point's own a surflet can reach from:
  1 for the separable kind, whose cell's eight corners are all there is,
  and more for a radial one wider than a cell."
  ^long [{:keys [kernel radius] :or {kernel :separable radius 1.5}}]
  (if (= kernel :radial) (long (math/ceil (double radius))) 1))

(defn basis
  "Noise as a sum of surflets: a basis `(fn [x y z])`.

  Options: `:seed` for the lattice's gradients, as `noise/gradient-basis`
  takes it; `:kernel`, `:separable` (the default, and equal to gradient
  noise) or `:radial`; `:radius` of a radial hump, in cells; and
  `:overrides`, a map of lattice points `[i j k]` to gradients to use
  there instead; and `:dimensions` 2 for noise in the plane z = 0,
  summing only the lattice points on it (evaluate it at z 0). A radial basis is scaled by `:gain` (default 0.45 at the
  default radius, where it sums 2.2 times larger), so that its spread
  matches the separable one's."
  ([] (basis {}))
  ([{:keys [seed kernel radius overrides gain dimensions] :or {seed 0 kernel :separable radius 1.5 dimensions 3} :as opts}]
   (let [^ints p (noise/permutation seed)
         radial? (= kernel :radial)
         radius (double radius)
         n (reach opts)
         gain (double (or gain (if radial? 0.45 1.0)))
         overrides (or overrides {})]
     (fn ^double [^double x ^double y ^double z]
       (let [i0 (long (math/floor x)) j0 (long (math/floor y)) k0 (long (math/floor z))
             lo (- 1 n)
             [klo khi] (if (= 2 dimensions) [(- k0) (- k0)] [lo n])]
         (loop [di lo dj lo dk (long klo) acc 0.0]
           (cond
             (> dk (long khi)) (* gain acc)
             (> dj n) (recur lo lo (inc dk) acc)
             (> di n) (recur lo (inc dj) dk acc)
             :else
             (let [i (+ i0 di) j (+ j0 dj) k (+ k0 dk)
                   dx (- x i) dy (- y j) dz (- z k)
                   w (if radial?
                       (hump (/ (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz))) radius))
                       (* (hump dx) (hump dy) (hump dz)))]
               (recur (inc di) dj dk
                      (if (zero? w)
                        acc
                        (let [[gx gy gz] (gradient-at p overrides [i j k])]
                          (+ acc (* w (+ (* (double gx) dx) (* (double gy) dy) (* (double gz) dz)))))))))))))))

(defn surflets-at
  "The lattice points whose surflets reach `[x y z]` under `opts`, with
  each one's share of the value there: `[[[i j k] value] ...]`, nonzero
  ones only. What the demo draws when it takes a point apart."
  [{:keys [seed overrides dimensions] :or {seed 0 dimensions 3} :as opts} [x y z]]
  (let [^ints p (noise/permutation seed)
        n (reach opts)
        x (double x) y (double y) z (double z)
        i0 (long (math/floor x)) j0 (long (math/floor y)) k0 (long (math/floor z))
        gain (double (or (:gain opts) (if (= :radial (:kernel opts)) 0.45 1.0)))]
    (vec (for [di (range (- 1 n) (inc n)) dj (range (- 1 n) (inc n))
               dk (if (= 2 dimensions) [(- k0)] (range (- 1 n) (inc n)))
               :let [ijk [(+ i0 di) (+ j0 dj) (+ k0 dk)]
                     v (* gain (double ((surflet (gradient-at p (or overrides {}) ijk) ijk opts) x y z)))]
               :when (not (zero? v))]
           [ijk v]))))
