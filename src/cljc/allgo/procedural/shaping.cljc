(ns allgo.procedural.shaping
  "The small functions every procedural texture is assembled from.

  Ebert, Musgrave, Peachey, Perlin and Worley, *Texturing & Modeling: A
  Procedural Approach*, chapter 2 (Peachey, \"Building Procedural
  Textures\"). Nothing here is deep. What makes them worth a namespace is
  that the rest of this library is written in terms of them: a procedural
  texture is a pile of arithmetic, and these are the shapes the arithmetic
  is made of.

  Three families, which is the whole vocabulary:

  * **Selectors** -- `step`, `pulse`, `pulse-train`. Hard edges, and the
    only things here that are discontinuous.
  * **Blends** -- `boxstep`, `smoothstep`, `smootherstep`, `mix`, `spline`.
    The same edges made continuous, which is what antialiasing wants and
    what makes a texture look like a material rather than a decal.
  * **Warps** -- `bias`, `gain`, `gamma`. These take [0, 1] to [0, 1] and
    move the values around inside it without reordering them, which is how
    you say \"more of this, but the same thing\" -- clouds with a harder
    edge, terrain with flatter plains.

  The argument order follows the book, which follows the RenderMan shading
  language: the value being shaped comes last. That reads backward for
  Clojure's threading macros and is worth keeping anyway, because these are
  the names and the signatures a reader of the book already has."
  (:require [allgo.math :as am]
            [clojure.math :as math]))

(defn mix
  "Linear blend: `a` at `t` = 0, `b` at `t` = 1.

  The book calls this `lerp`; RenderMan calls it `mix`. It is the join in
  every layered texture -- two materials and a mask."
  ^double [^double a ^double b ^double t]
  (am/lerp a b t))

(defn remap
  "`x`, which lies in [`a0`, `a1`], moved to where it falls in [`b0`, `b1`].

  Not in the book by name; it is the line that appears at the seam of every
  two functions here that disagree about their range -- noise coming out of
  [-1, 1] and into the [0, 1] the blends want."
  ;; Five arguments, so no primitive hints: Clojure's primitive function
  ;; interfaces stop at four.
  [a0 a1 b0 b1 x]
  (let [a0 (double a0) a1 (double a1) b0 (double b0) b1 (double b1)]
    (if (== a0 a1)
      b0
      (+ b0 (* (- b1 b0) (/ (- (double x) a0) (- a1 a0)))))))

(defn step
  "0 below `a`, 1 at or above it. The hard edge."
  ^double [^double a ^double x]
  (if (< x a) 0.0 1.0))

(defn pulse
  "1 between `a` and `b`, 0 outside. `step(a) - step(b)`."
  ^double [^double a ^double b ^double x]
  (- (step a x) (step b x)))

(defn boxstep
  "The linear ramp from 0 at `a` to 1 at `b`, clamped outside.

  The cheapest antialiased `step` there is: `a` and `b` are placed a pixel
  apart and the edge is exactly one pixel wide. `smoothstep` is the same
  idea with the corners taken off."
  ^double [^double a ^double b ^double x]
  (am/clamp (/ (- x a) (if (== a b) 1e-30 (- b a))) 0.0 1.0))

(defn smoothstep
  "The Hermite ramp from 0 at `a` to 1 at `b`, flat at both ends.

  `3t^2 - 2t^3`. First derivative zero at either end, so two of these back
  to back do not show a crease where they meet."
  ^double [^double a ^double b ^double x]
  (let [t (boxstep a b x)]
    (* t t (- 3.0 (* 2.0 t)))))

(defn smootherstep
  "`smoothstep` with the second derivative zero at the ends as well.

  `6t^5 - 15t^4 + 10t^3`, Perlin's improved-noise fade curve (chapter 12).
  Gradient noise interpolated with the Hermite curve shows a faint grid of
  creases along the lattice planes, because the second derivative jumps
  there; this is the fix, and `allgo.procedural.noise` uses it throughout."
  ^double [^double a ^double b ^double x]
  (let [t (boxstep a b x)]
    (* t t t (+ (* t (- (* t 6.0) 15.0)) 10.0))))

(defn bias
  "Perlin's bias: bends [0, 1] onto itself, pushing values up or down.

  `b` = 0.5 is the identity. Above it, values rise; below, they fall.
  `pow(t, ln(b)/ln(0.5))`, which is the power that maps 0.5 to `b`."
  ^double [^double b ^double t]
  (cond
    (<= t 0.0) 0.0
    (>= t 1.0) 1.0
    (<= b 0.0) 0.0
    (>= b 1.0) 1.0
    :else (math/pow t (/ (math/log b) (math/log 0.5)))))

(defn gain
  "Perlin's gain: pushes values away from the middle, or toward it.

  `g` = 0.5 is the identity. Above it the curve steepens around 0.5 and
  contrast rises; below it the curve flattens and everything drifts to
  gray. Built from two halves of `bias` so that it fixes 0, 0.5 and 1."
  ^double [^double g ^double t]
  (if (< t 0.5)
    (* 0.5 (bias (- 1.0 g) (* 2.0 t)))
    (- 1.0 (* 0.5 (bias (- 1.0 g) (- 2.0 (* 2.0 t)))))))

(defn gamma
  "`pow(x, 1/g)`. Brightens the midtones without moving the ends."
  ^double [^double g ^double x]
  (if (<= x 0.0) 0.0 (math/pow x (/ 1.0 g))))

(defn sawtooth
  "The fractional part of `x / period`, in [0, 1). Repeats the domain."
  ^double [^double period ^double x]
  (let [p (if (zero? period) 1.0 period)
        t (/ x p)]
    (- t (math/floor t))))

(defn triangle
  "A triangle wave over `period`, in [0, 1]. `sawtooth` folded in half."
  ^double [^double period ^double x]
  (let [t (sawtooth period x)]
    (if (< t 0.5) (* 2.0 t) (- 2.0 (* 2.0 t)))))

(defn pulse-train
  "`pulse` repeated every `period`: 1 for the first `edge` of each period."
  ^double [^double edge ^double period ^double x]
  (pulse 0.0 (/ edge (if (zero? period) 1.0 period)) (sawtooth period x)))

(defn spline
  "Catmull-Rom through `knots`, with `t` in [0, 1] spanning the whole run.

  The book's color-ramp primitive: four or more control values, and a
  curve that passes through each of them. The first and last knots are
  phantoms -- the curve starts at the second and ends at the second to
  last -- which is what makes the interior tangents well defined, so give
  it at least four and expect to duplicate the ends."
  ^double [^double t coll]
  (let [knots (vec coll)
        n (count knots)]
    (cond
      (zero? n) 0.0
      (< n 4) (double (nth knots (long (am/clamp (math/floor (* t (double n))) 0.0 (double (dec n))))))
      :else
      (let [spans (double (- n 3))
            x (* (am/clamp t 0.0 1.0) spans)
            i (long (min (double (dec spans)) (math/floor x)))
            u (- x i)
            k (fn ^double [^long j] (double (nth knots (+ i j))))
            p0 (k 0) p1 (k 1) p2 (k 2) p3 (k 3)]
        ;; The Catmull-Rom basis, written out rather than as a matrix
        ;; product: one fused polynomial in u, evaluated by Horner.
        (* 0.5
           (+ (* 2.0 p1)
              (* u (+ (- p2 p0)
                      (* u (+ (- (+ (* 2.0 p0) (* 4.0 p2))
                                 (+ (* 5.0 p1) p3))
                              (* u (+ (- (* 3.0 p1) p0)
                                      (- (* 3.0 p2)) p3))))))))))))
