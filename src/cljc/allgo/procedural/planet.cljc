(ns allgo.procedural.planet
  "Procedural planets: the whole of a world as a function of a direction.

  Musgrave, in Ebert et al., *Texturing & Modeling: A Procedural
  Approach*, chapter 20, \"MojoWorld: Building Procedural Planets\".
  Everything under this namespace rests on one decision, and it is worth
  stating before any of the code:

  **The terrain is a function of a point in space, not of a point on a
  map.** A height field indexed by latitude and longitude is the obvious
  representation and it cannot be made to work. The poles are singular --
  every longitude meets there -- so features are stretched into streaks
  approaching them and pinched to nothing at them, and the seam at 180
  degrees never quite closes. Sample the terrain at the 3D unit vector
  instead and none of that exists: there is no seam, because there is
  nowhere in the domain that two different coordinates name the same
  place, and there are no poles, because the sphere has no preferred axis.
  That is what \"solid texture\" means and it is the reason chapter 20
  spends its effort on 3D bases.

  What a planet is here, then:

      radius R, relief A, and a terrain basis f
      surface(d) = R + A * f(d)        for a unit direction d

  Add a sea level and the land/water division falls out for free -- the
  coastline is a contour of `f`, which is why procedural coastlines have
  the right fractal wiggle without anyone drawing one. Add a second
  fractal for cloud cover and a third for the surface color and you have
  a planet you can fly to, at any magnification, having stored nothing.

  ## Composing one

  The parts come from elsewhere: `allgo.procedural.noise` for bases,
  `allgo.procedural.fractal` for the constructions, `allgo.procedural.shaping`
  for the blends. `default-terrain` is a worked example of putting them
  together rather than a privileged implementation -- read it, then write
  your own, which is the entire argument of the book.

  ## Scale

  `:radius` and `:relief` are in whatever units you like, as long as they
  agree. Earth's numbers are a useful sanity check and a good reminder of
  how flat planets really are: 6378 km of radius against about 9 km of
  relief, so the terrain is a little over a thousandth of the sphere.
  Every rendering of a planet you have seen exaggerates it."
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as shaping]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Getting about on a sphere

(defn direction
  "The unit vector at `latitude` and `longitude`, both in radians.

  Only for talking to the outside world -- nothing inside a planet is
  parameterized this way, for the reason in the namespace docstring."
  [latitude longitude]
  (let [lat (double latitude) lon (double longitude)
        c (math/cos lat)]
    [(* c (math/cos lon)) (math/sin lat) (* c (math/sin lon))]))

(defn latitude
  "Latitude of a direction, in radians, from -pi/2 to pi/2."
  ^double [[x y z]]
  (math/atan2 (double y)
              (math/sqrt (+ (* (double x) (double x)) (* (double z) (double z))))))

(defn longitude
  "Longitude of a direction, in radians, from -pi to pi."
  ^double [[x _ z]]
  (math/atan2 (double z) (double x)))

(defn fibonacci-directions
  "`n` unit vectors spread as evenly over the sphere as a simple rule can
  manage.

  Successive points are separated by the golden angle, which is the
  irrational rotation that packs most evenly; the result has no clusters
  and no bands, unlike sampling latitude and longitude independently. Used
  here to measure a terrain function fairly, and useful for anything that
  needs to ask a planet a question at every point equally."
  [n]
  (let [n (long n)
        golden (* math/PI (- 3.0 (math/sqrt 5.0)))]
    (mapv (fn [i]
            (let [i (double i)
                  y (- 1.0 (/ (* 2.0 i) (max 1.0 (double (dec n)))))
                  r (math/sqrt (max 0.0 (- 1.0 (* y y))))
                  theta (* golden i)]
              [(* r (math/cos theta)) y (* r (math/sin theta))]))
          (range n))))

(defn tangent-basis
  "Two unit vectors perpendicular to `d` and to each other.

  Picked by crossing with whichever axis `d` is least aligned to, which is
  the standard dodge for the fact that you cannot comb a sphere: any fixed
  choice of reference axis fails somewhere, and this fails nowhere by not
  fixing one."
  [d]
  (let [[_ y _] d
        axis (if (< (abs (double y)) 0.9) [0.0 1.0 0.0] [1.0 0.0 0.0])
        u (v3/normalize (v3/cross d axis))
        v (v3/cross d u)]
    [u v]))

;; ---------------------------------------------------------------------------
;; Terrain

(defn terrain-range
  "`[lowest highest]` of `f` over the unit sphere, from `n` directions.

  A multifractal's range depends on every one of its parameters at once
  and cannot be predicted, so a planet built on one has to measure before
  it can say what its own sea level means. This is that measurement, taken
  on the sphere rather than in a box because the sphere is the only part
  of the domain a planet ever visits."
  ([f] (terrain-range f 4096))
  ([f n]
   (reduce (fn [[lo hi] [x y z]]
             (let [v (noise/sample f (double x) (double y) (double z))]
               [(min (double lo) v) (max (double hi) v)]))
           [##Inf ##-Inf]
           (fibonacci-directions n))))

(defn default-terrain
  "A continents-and-mountains terrain, as an example of composing one.

  Three layers, which between them are most of what chapter 20 says about
  making a world look like a world:

  * **Continents.** A low-frequency `fbm`, biased so that it spends more
    of its time low than high. That is what puts a minority of the surface
    above sea level and gives oceans that are mostly one connected body.
  * **Mountains.** A second fractal -- a ridged multifractal by default,
    for its branching ridge lines -- normalized and then multiplied by the
    continental mask, so that ranges rise from land and the sea floor
    stays comparatively smooth. Real ocean floor is smoother than real
    land for a reason procedural terrain does not have, namely that it is
    buried in sediment, but the result reads correctly.
  * **Distortion.** The whole domain warped by a vector basis before any
    of it is evaluated, which breaks up the residual isotropy of summed
    noise. The most expensive knob in the chapter -- four noise
    evaluations per sample instead of one -- and the one whose absence is
    most obvious once you have seen it switched on.

  Everything is swappable, which is the point of it being here rather than
  inlined into `planet`:

  * `:basis-fn` -- a function of a seed to a basis. Anything from
    `allgo.procedural.noise` or `allgo.procedural.cellular` will do, and
    the two non-lattice bases cost about thirty times what the lattice
    ones do, so mind the octave count.
  * `:construction` -- which fractal builds the mountains, by the names
    in `allgo.procedural.fractal/constructions`.
  * `:seed` -- makes a different planet.
  * `:relief-mix` -- how much of the range the mountains are allowed.
  * `:sample-spacing` -- how far apart, in radians on the unit sphere,
    this terrain is going to be sampled. Given it, each layer's octave
    count is cut to what that sampling can actually carry, by
    `allgo.procedural.fractal/octaves-for`. This is not an optimization
    with a quality cost, it is the opposite: octaves finer than the
    samples do not arrive as detail, they arrive as speckle that moves
    when the camera does, and they are paid for at full price. Leave it
    out only if the caller does not know its own sampling rate."
  ([] (default-terrain {}))
  ([{:keys [seed basis-fn construction continent-frequency mountain-frequency
            distortion continent-bias relief-mix octaves lacunarity H offset gain
            sample-spacing]
     :or {seed 0 construction :ridged-multifractal
          continent-frequency 0.8 mountain-frequency 2.2
          distortion 0.35 continent-bias 0.42 relief-mix 0.55
          octaves 9.0 lacunarity 2.0 H 0.9 offset 1.0 gain 2.0}}]
   (let [seed (long seed)
         ;; A layer's sampling spacing in its own domain is the spacing on
         ;; the sphere times whatever frequency it was scaled by. The
         ;; margin is for the domain distortion, which compresses the
         ;; domain as much as it stretches it, and so puts some of the
         ;; terrain at a higher frequency than its nominal one.
         limit (fn ^double [^double asked ^double frequency]
                 (if sample-spacing
                   (min asked
                        (fractal/octaves-for (* (double sample-spacing) frequency
                                                (+ 1.0 (double distortion)))
                                             lacunarity))
                   asked))
         basis-fn (or basis-fn #(noise/gradient-basis {:seed %}))
         warp (noise/vector-basis {:seed (+ seed 101)})
         base (fn [s]
                (let [b (basis-fn (+ seed (long s)))]
                  (if (pos? (double distortion))
                    (noise/distorted b warp distortion)
                    b)))
         continents (-> (fractal/fbm (base 0) {:octaves (limit 5.0 continent-frequency)
                                               :H 1.0 :lacunarity lacunarity})
                        (noise/shaped #(shaping/bias continent-bias %))
                        (noise/scaled continent-frequency))
         raw (-> (fractal/build construction (base 7)
                                {:octaves (limit octaves mountain-frequency)
                                 :H H :lacunarity lacunarity
                                 :offset offset :gain gain})
                 (noise/scaled mountain-frequency))
         ;; Every construction has a different natural range -- the ridged
         ;; one never goes below zero, fBm straddles it -- so the mountains
         ;; are measured and recentered rather than assumed. Without this,
         ;; changing the construction would also move the sea level.
         [mlo mhi] (terrain-range raw 1024)
         mspan (let [d (- (double mhi) (double mlo))] (if (zero? d) 1.0 d))
         relief-mix (double relief-mix)]
     (fn ^double [^double x ^double y ^double z]
       (let [c (noise/sample continents x y z)
             ;; Mountains only where there is land to put them on, faded in
             ;; across the coastal shelf so the shore is not a cliff
             ;; everywhere.
             mask (shaping/smoothstep -0.15 0.45 c)
             m (- (* 2.0 (/ (- (noise/sample raw x y z) (double mlo)) mspan)) 1.0)]
         (+ c (* relief-mix mask m)))))))

;; ---------------------------------------------------------------------------
;; The planet

(defn planet
  "Assembles a planet from its parts.

  * `:radius` -- of the datum, the sphere the terrain displaces.
  * `:relief` -- the displacement at full scale, in the same units. The
    terrain function is normalized onto [-1, 1] first, so this is the
    real peak-to-datum height whatever the fractal underneath does.
  * `:sea-level` -- where the water sits, in those same units, so 0 is
    halfway up the range and negative values give a drier world.
  * `:terrain` -- a basis; `default-terrain` if not given.
  * `:clouds` -- a basis for cloud density, or nil for a clear sky.
  * `:cloud-cover` -- how much of the sky the clouds take, in [0, 1].
  * `:samples` -- how many directions to measure the terrain over.
  * `:sample-spacing` -- how far apart, in radians, the planet is going to
    be asked about. Passed on to the default terrain and the default
    clouds so that neither generates detail finer than will survive.

  The measurement happens here, once, so that everything downstream can
  treat elevation as a number in known units. It is a measurement and not
  a proof: `:samples` directions are tried, and a direction that was not
  among them can come out a percent or two past the relief. Raise
  `:samples` if that matters, and do not rely on the bound being exact."
  ([] (planet {}))
  ([{:keys [radius relief sea-level terrain clouds cloud-cover samples seed
            sample-spacing]
     :or {radius 1.0 relief 0.035 sea-level 0.0 cloud-cover 0.45
          samples 4096 seed 0}
     :as opts}]
   (let [terrain (or terrain (default-terrain {:seed seed
                                               :sample-spacing sample-spacing}))
         clouds (if (contains? opts :clouds)
                  clouds
                  (fractal/fbm (noise/gradient-basis {:seed (+ (long seed) 991)})
                               {:octaves (if sample-spacing
                                           (min 7.0 (fractal/octaves-for sample-spacing 2.0))
                                           7.0)
                                :H 0.85}))
         [lo hi] (terrain-range terrain samples)
         span (let [s (- (double hi) (double lo))] (if (zero? s) 1.0 s))]
     {:radius (double radius)
      :relief (double relief)
      :sea-level (double sea-level)
      :terrain terrain
      :clouds clouds
      :cloud-cover (double cloud-cover)
      :terrain-lo (double lo)
      :terrain-hi (double hi)
      ;; Terrain value -> [-1, 1], so that :relief means what it says.
      :terrain-scale (/ 2.0 span)
      :terrain-shift (- (+ (/ (* 2.0 (double lo)) span) 1.0))})))

(defn elevation
  "Height of the ground above the datum at `d`, in world units.

  Signed: negative is below the datum, which is usually under water. This
  is the ground, not the surface -- the sea does not raise it."
  ^double [{:keys [terrain ^double relief ^double terrain-scale ^double terrain-shift]}
           [x y z]]
  (* relief (+ (* terrain-scale (noise/sample terrain (double x) (double y) (double z)))
               terrain-shift)))

(defn ocean?
  "Whether `d` is under water."
  [{:keys [^double sea-level] :as planet} d]
  (< (elevation planet d) sea-level))

(defn depth
  "How deep the water is at `d`, or 0 on land."
  ^double [{:keys [^double sea-level] :as planet} d]
  (max 0.0 (- sea-level (elevation planet d))))

(defn surface-radius
  "Distance from the center to what you would see at `d`: the ground, or
  the sea surface where the ground is below it."
  ^double [{:keys [^double radius ^double sea-level] :as planet} d]
  (+ radius (max sea-level (elevation planet d))))

(defn surface-point
  "The visible surface at `d`, as a point in space."
  [planet d]
  (v3/scale (v3/normalize d) (surface-radius planet d)))

(defn altitude
  "Height of a point in space above the visible surface below it.

  Negative inside the planet. This is the function a ray marcher asks
  about, and the reason the planet is a value rather than a mesh:
  `allgo.procedural.qaeb` needs exactly this and nothing else."
  ^double [planet p]
  (let [r (v3/length p)]
    (if (zero? r)
      (- (surface-radius planet [0.0 1.0 0.0]))
      (- r (surface-radius planet (v3/scale p (/ 1.0 r)))))))

(defn surface-normal
  "The unit normal of the displaced surface at `d`.

  Two finite differences in the tangent plane, crossed. A procedural
  surface has no analytic derivative on offer, so this costs four more
  terrain evaluations -- which for an eight-octave multifractal is the
  dominant cost of rendering one, and the reason `epsilon` is worth
  tuning: too small and the difference is lost in rounding, too large and
  the fine detail is smoothed away before it is ever shaded."
  ([planet d] (surface-normal planet d 1e-3))
  ([planet d epsilon]
   (let [d (v3/normalize d)
         e (double epsilon)
         [u v] (tangent-basis d)
         at (fn [off]
              (let [p (v3/normalize (v3/add d off))]
                (v3/scale p (surface-radius planet p))))
         pu+ (at (v3/scale u e)) pu- (at (v3/scale u (- e)))
         pv+ (at (v3/scale v e)) pv- (at (v3/scale v (- e)))
         n (v3/cross (v3/sub pu+ pu-) (v3/sub pv+ pv-))]
     ;; The cross product comes out along the outward radial for a sphere
     ;; traversed this way; flip if the terrain has folded it over.
     (if (neg? (v3/dot n d))
       (v3/normalize (v3/negate n))
       (v3/normalize n)))))

(defn slope
  "Angle between the surface normal and straight up, in radians.

  Zero on the flat, pi/2 on a vertical face. The other half of Musgrave's
  texturing rule: altitude says what a place is, slope says whether
  anything can stay on it."
  (^double [planet d] (slope planet d 1e-3))
  (^double [planet d epsilon]
   (let [d (v3/normalize d)
         n (surface-normal planet d epsilon)]
     (math/acos (shaping/clamp -1.0 1.0 (v3/dot n d))))))

(defn cloud-cover
  "Cloud opacity at `d`, in [0, 1].

  A fractal thresholded with a soft edge rather than a hard one: the
  threshold is what makes a sky mostly clear with distinct clouds in it
  instead of uniform haze, and the softness is what keeps the clouds from
  looking cut out with scissors. `:cloud-cover` moves the threshold."
  ^double [{:keys [clouds ^double cloud-cover]} [x y z]]
  (if (nil? clouds)
    0.0
    (let [v (noise/sample clouds (double x) (double y) (double z))
          ;; Cover 0 puts the threshold above everything, cover 1 below it.
          edge (- 1.0 (* 2.0 cloud-cover))]
      (shaping/smoothstep edge (+ edge 0.45) v))))

;; ---------------------------------------------------------------------------
;; Color

(def terran
  "Musgrave's Terran palette: what a temperate planet is made of.

  Seven colors and the rules below are the whole of it. There is no
  climate model here and there does not need to be one -- altitude, slope
  and latitude between them predict what covers a piece of ground well
  enough to fool the eye, which is the standard the chapter sets."
  {:abyss   [0.01 0.05 0.16]
   :ocean   [0.04 0.20 0.38]
   :shallow [0.10 0.42 0.52]
   :beach   [0.76 0.70 0.50]
   :lowland [0.18 0.34 0.15]
   :upland  [0.35 0.33 0.20]
   :rock    [0.40 0.38 0.36]
   :snow    [0.96 0.97 0.99]})

(defn surface-color
  "The color of the ground or sea at `d`, as `[r g b]` in [0, 1].

  Layered exactly as the chapter describes, each layer a `smoothstep`
  between two of the palette entries:

  1. Under water, blend from shallow to abyssal with depth.
  2. On land, a ramp from beach through lowland and upland to bare rock,
     keyed on height above sea level as a fraction of the relief.
  3. Snow above a line that comes *down* toward the poles, because the
     snow line is a temperature contour and temperature falls with both
     altitude and latitude.
  4. Steep ground is bare rock whatever else it would have been, and too
     steep for snow to lie on, which is what puts dark faces on white
     mountains and stops the whole range looking like icing.

  Pass `:slope` and `:elevation` if you already know them. Between them
  they are the entire cost of this function: the slope is four more
  terrain evaluations, and anything drawing a mesh has the vertex normals
  in hand already -- a better answer anyway, since it is the slope of the
  surface actually being drawn -- and the elevation was needed to place
  the vertex in the first place. Supplying both makes coloring free."
  ([planet d] (surface-color planet d {}))
  ([{:keys [^double relief ^double sea-level] :as planet} d
    {:keys [snow-line steep pole-effect palette epsilon]
     :or {snow-line 0.55 steep 0.6 pole-effect 0.45 epsilon 1e-3}
     given-slope :slope given-elevation :elevation}]
   (let [{:keys [abyss ocean shallow beach lowland upland rock snow]}
         (merge terran palette)
         d (v3/normalize d)
         h (double (or given-elevation (elevation planet d)))]
     (if (< h sea-level)
       ;; Water: two blends rather than one, because the eye reads the
       ;; shelf break. Depth is measured against the deepest the terrain
       ;; can go, so the shelf stays a shelf as the sea level moves.
       (let [deepest (max 1e-9 (+ relief sea-level))
             t (shaping/clamp 0.0 1.0 (/ (- sea-level h) deepest))]
         (if (< t 0.25)
           (v3/lerp shallow ocean (shaping/smoothstep 0.0 0.25 t))
           (v3/lerp ocean abyss (shaping/smoothstep 0.25 1.0 t))))
       (let [;; Height above the waterline as a fraction of what is left
             ;; of the relief, so raising the sea level does not also
             ;; repaint the mountains.
             span (max 1e-9 (- relief sea-level))
             t (shaping/clamp 0.0 1.0 (/ (- h sea-level) span))
             steepness (/ (double (or given-slope (slope planet d epsilon)))
                          (* 0.5 math/PI))
             ground (cond
                      (< t 0.02) (v3/lerp beach lowland (shaping/smoothstep 0.0 0.02 t))
                      (< t 0.35) (v3/lerp lowland upland (shaping/smoothstep 0.02 0.35 t))
                      :else (v3/lerp upland rock (shaping/smoothstep 0.35 0.7 t)))
             bare (shaping/smoothstep (* 0.7 steep) steep steepness)
             ground (v3/lerp ground rock bare)
             ;; The snow line falls toward the poles.
             lat (abs (latitude d))
             line (max 0.0 (- snow-line (* pole-effect (/ lat (* 0.5 math/PI)))))
             snowy (* (shaping/smoothstep line (+ line 0.12) t)
                      (- 1.0 (shaping/smoothstep (* 0.8 steep) steep steepness)))]
         (v3/lerp ground snow snowy))))))

;; ---------------------------------------------------------------------------
;; Baking

(defn equirectangular
  "Samples `f` over a `width` by `height` latitude/longitude grid.

  Row 0 is the north pole, column 0 is longitude -pi. Returns whatever `f`
  returns, row by row, in a vector -- a map to hand a renderer, once, in
  exchange for not evaluating the fractal per frame.

  This is the one place the projection is allowed, and only because
  something outside wants a rectangle. Note what it costs: the polar rows
  sample the same few square kilometers hundreds of times over while the
  equatorial rows undersample, which is exactly the waste the solid
  formulation exists to avoid."
  [f width height]
  (let [w (long width) h (long height)]
    (into []
          (for [j (range h)
                i (range w)]
            (let [lat (- (* 0.5 math/PI) (* math/PI (/ (+ 0.5 (double j)) h)))
                  lon (- (* 2.0 math/PI (/ (+ 0.5 (double i)) w)) math/PI)]
              (f (direction lat lon)))))))
