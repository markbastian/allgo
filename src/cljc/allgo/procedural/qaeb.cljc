(ns allgo.procedural.qaeb
  "Ray tracing a surface that exists only as a function.

  Musgrave, in Ebert et al., *Texturing & Modeling: A Procedural
  Approach*, chapter 17, \"QAEB Rendering for Procedural Models\" --
  quasi-analytic, error-bounded. The problem it solves is the one chapter
  20 leaves you with: a procedural planet has no geometry to intersect. It
  has a function that says how far above the surface a point in space is,
  and nothing else. There is no polygon to hit, no bounding volume that is
  not the whole planet, and no way to know the function's Lipschitz bound,
  so sphere tracing is out too.

  ## What QAEB does instead

  March along the ray in steps, and take the first step that ends up below
  the surface as a bracket around the crossing. The two ideas that make
  that more than brute force:

  * **The step grows with distance.** `dt = epsilon * t`. A step at twice
    the distance projects to the same size on screen, so the error is
    bounded *in the image* rather than in the world -- which is the only
    place it matters. A fixed step fine enough for the foreground would
    waste almost all of its samples in the distance, and one coarse enough
    for the distance would miss everything near the camera. This is the
    \"error-bounded\" half of the name, and `epsilon` is, near enough, the
    pixel size in radians.

  * **The crossing is interpolated, not searched.** Once a step straddles
    the surface, the altitude at either end is known, and one linear
    interpolation lands within the error bound already accepted. That is
    the \"quasi-analytic\" half. Bisection afterwards is available and
    usually unnecessary.

  What you give up is the guarantee. A ridge thinner than a step can be
  stepped over, and that shows as a dropout along a silhouette rather than
  as noise. Smaller `epsilon`, or `:max-slope` to bound the step by how
  fast the surface can rise, buys it back at a price.

  ## Coherence

  Neighbouring rays hit at nearly the same distance, so starting a ray at
  the previous one's hit distance less a margin removes most of the march.
  Chapter 17 gets a large constant factor that way. Nothing here does it
  for you -- it is a property of the traversal order, which is the
  caller's -- but `:t-min` is the hook: pass the last hit, less enough to
  be safe."
  (:require [allgo.geometry.vec3 :as v3]))

(defn intersect
  "Marches `altitude` along a ray and returns where it first goes below zero.

  `altitude` is a function of a point to its signed height above the
  surface -- positive outside, negative inside.
  `allgo.procedural.planet/altitude` is one; a height field's is
  `(- (:y p) (h (:x p) (:z p)))`.

  Returns `{:t :point :steps}`, or nil if the ray reaches `:t-max` without
  crossing.

  * `:t-min` / `:t-max` -- the segment to search. This is not a detail: a
    step proportional to distance has no scale of its own at the origin,
    so a march that starts at the eye spends a step per `epsilon` of
    *relative* distance -- from 0.001 to 4 at an epsilon of 0.001 is
    eight thousand of them, almost all through empty space. Start at the
    bounding volume (`sphere-entry`) or at the previous ray's hit less a
    margin, and the same picture costs a few hundred. The default of
    0.001 is a near plane, not a recommendation.
  * `:epsilon` -- step as a fraction of distance; the pixel's angular
    size. 1/1000 is a reasonable place to start.
  * `:min-step` -- a floor under the step, in world units. Zero by
    default, which leaves the march purely proportional.
  * `:max-slope` -- if given, the largest the surface may rise per unit
    along the ray. The step is then also capped at `altitude / max-slope`,
    which cannot skip over anything and turns the march into a true sphere
    trace where the bound is known.
  * `:max-steps` -- the safety net.
  * `:refine` -- bisections after the interpolation. Zero is normally
    right."
  ([altitude origin direction] (intersect altitude origin direction {}))
  ([altitude origin direction
    {:keys [t-min t-max epsilon min-step max-slope max-steps refine]
     :or {t-min 1e-3 t-max 1e9 epsilon 1e-3 min-step 0.0 max-steps 20000 refine 0}}]
   (let [dir (v3/normalize direction)
         t-max (double t-max)
         epsilon (double epsilon)
         min-step (double min-step)
         max-slope (when max-slope (double max-slope))
         max-steps (long max-steps)
         at (fn ^double [^double t] (double (altitude (v3/add-scaled origin dir t))))]
     (loop [t (double t-min)
            a (at (double t-min))
            steps 0]
       (cond
         (neg? a) {:t t :point (v3/add-scaled origin dir t) :steps steps}
         (or (> t t-max) (>= steps max-steps)) nil
         :else
         (let [dt (max min-step
                       (let [d (* epsilon t)]
                         (if max-slope (min d (/ a max-slope)) d))
                       ;; A purely proportional step at t = 0 is no step at
                       ;; all; this is the floor that gets the march moving
                       ;; when a caller has not given it a near plane.
                       1e-12)
               t' (+ t dt)
               a' (at t')]
           (if (neg? a')
             ;; Bracketed. One interpolation lands inside the error already
             ;; accepted; bisect afterwards only if asked.
             (let [t0 (loop [lo t hi t' alo a ahi a' n (long refine)]
                        (let [mid (if (== alo ahi)
                                    (* 0.5 (+ lo hi))
                                    (+ lo (* (- hi lo) (/ alo (- alo ahi)))))]
                          (if (zero? n)
                            mid
                            (let [am (at mid)]
                              (if (neg? am)
                                (recur lo mid alo am (dec n))
                                (recur mid hi am ahi (dec n)))))))]
               {:t t0 :point (v3/add-scaled origin dir t0) :steps (inc steps)})
             (recur t' a' (inc steps)))))))))

(defn shadowed?
  "Whether anything blocks the ray from `point` towards `light`.

  The same march with the answer thrown away, which is all a shadow ray
  is. `:t-min` has to be far enough off the surface to clear the point the
  ray started from, or every point shadows itself."
  ([altitude point light] (shadowed? altitude point light {}))
  ([altitude point light opts]
   (some? (intersect altitude point light
                     (merge {:t-min 1e-4 :epsilon 4e-3} opts)))))

(defn sphere-entry
  "Where a ray enters and leaves the sphere of radius `r` about the origin.

  `[near far]`, both possibly negative if the ray starts inside, or nil if
  it misses. Not part of QAEB, but every march against a planet wants it:
  the terrain lives in a shell, and starting the march at the shell rather
  than at the eye removes the empty space between them."
  [origin direction r]
  (let [d (v3/normalize direction)
        r (double r)
        b (v3/dot origin d)
        c (- (v3/dot origin origin) (* r r))
        disc (- (* b b) c)]
    (when-not (neg? disc)
      (let [s (Math/sqrt disc)]
        [(- (- b) s) (+ (- b) s)]))))
