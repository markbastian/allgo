(ns allgo.procedural.atmosphere
  "Why the sky is blue, the sunset is red, and distant hills are pale.

  Musgrave, in Ebert et al., *Texturing & Modeling: A Procedural
  Approach*, chapter 18, \"Atmospheric Models\", and the last thing a
  procedural planet needs before it stops looking like a painted ball.
  Terrain rendered against black is unconvincing at any level of detail,
  because most of what tells the eye how big and how far away a landscape
  is comes from the air in front of it, not from the landscape.

  Two ingredients, and they are both in the chapter's title:

  ## Density

  The gaseous atmospheric density distribution, which for any planet that
  has settled down is an exponential: pressure at a height supports the
  weight of everything above it, so density falls by a constant factor
  every `:scale-height`. Earth's is about 8 km against 6378 km of radius,
  which is why the atmosphere is a film and not a shell -- and why the
  sunset, which looks through hundreds of kilometers of it, is so much
  redder than noon, which looks through eight.

  ## Scattering

  Rayleigh scattering, off molecules much smaller than the wavelength,
  goes as 1 / lambda^4. Blue is scattered about five and a half times as
  strongly as red, so blue is what arrives from directions other than the
  sun -- the sky -- and red is what is left of the sun after a long path
  through air has taken the blue out of it. One constant, one exponent,
  and both of the sky's famous colors.

  Mie scattering, off droplets and dust comparable to the wavelength, is
  near enough wavelength-independent and strongly forward-biased. It is
  the white haze around the sun and the gray of a humid horizon.

  ## What this computes

  Single scattering: light from the sun, scattered once toward the eye,
  attenuated on both legs of that journey. Multiple scattering is what
  fills in the deep blue overhead and the light in shadows, and is not
  here -- it needs a solver, not a quadrature. Single scattering is what
  the chapter describes and is enough to be convincing.

  Units are the planet's: pass radii and scale heights in whatever
  `allgo.procedural.planet` is using. The scattering coefficients are
  per unit length in those units, so they have to be rescaled along with
  everything else."
  (:require [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(def wavelengths-rgb
  "Representative wavelengths for red, green and blue, in nanometers."
  [680.0 550.0 440.0])

(defn rayleigh-coefficients
  "Scattering coefficients for RGB, in proportion to 1 / lambda^4.

  `strength` is the coefficient at 550 nm, in reciprocal length units of
  whatever the planet is measured in. Everything characteristic about the
  sky's color is in the ratios, which come out near 1 : 2.3 : 5.6 --
  that is the whole of why it is blue."
  [strength]
  (let [strength (double strength)
        green (nth wavelengths-rgb 1)]
    (mapv (fn [l] (* strength (math/pow (/ green (double l)) 4.0)))
          wavelengths-rgb)))

(defn atmosphere
  "The parameters of one atmosphere, with Earth-like defaults scaled to a
  planet of radius 1.

  * `:radius` -- of the planet's datum.
  * `:thickness` -- how far up the model goes. Density is negligible
    above a handful of scale heights; there is no point integrating
    through vacuum.
  * `:scale-height` -- the e-folding height of the density.
  * `:rayleigh` -- coefficients per channel, from `rayleigh-coefficients`.
  * `:mie` / `:mie-g` -- haze strength and how forward-biased it is.
    `:mie-g` near 0.76 is the usual fit for atmospheric aerosol.
  * `:sun` -- the sun's color and intensity, per channel."
  ([] (atmosphere {}))
  ([{:keys [radius thickness scale-height rayleigh mie mie-g sun]
     :or {radius 1.0 thickness 0.025 scale-height 0.0013
          mie 0.0021 mie-g 0.76 sun [22.0 22.0 22.0]}}]
   {:radius (double radius)
    :thickness (double thickness)
    :top (+ (double radius) (double thickness))
    :scale-height (double scale-height)
    :rayleigh (or rayleigh (rayleigh-coefficients 21.0))
    :mie (double mie)
    :mie-g (double mie-g)
    :sun (mapv double sun)}))

(defn density
  "Relative density at a point: 1 at the surface, falling by `e` every
  scale height. The GADD, and the only thing about the air's structure
  this model knows."
  ^double [{:keys [^double radius ^double scale-height]} p]
  (let [h (- (v3/length p) radius)]
    (math/exp (- (/ (max 0.0 h) scale-height)))))

(defn- ray-sphere
  "`[near far]` where a unit-direction ray meets the sphere of radius `r`
  about the origin, or nil."
  [origin dir ^double r]
  (let [b (v3/dot origin dir)
        c (- (v3/dot origin origin) (* r r))
        disc (- (* b b) c)]
    (when-not (neg? disc)
      (let [s (math/sqrt disc)]
        [(- (- b) s) (+ (- b) s)]))))

(defn optical-depth
  "Integrated density from `origin` along `dir` for `length`, in `steps`.

  The quantity every transmittance is an exponential of: how much air the
  light went through, in units of how much it would go through in one
  length unit at the surface. Trapezoidal, because the integrand is
  smooth and the sample count is what costs."
  ;; Five arguments, so no primitive hints on them: Clojure's primitive
  ;; function interfaces stop at four.
  [atm origin dir length steps]
  (let [dir (v3/normalize dir)
        n (max 1 (long steps))
        dt (/ (double length) (double n))]
    (loop [i 0 acc 0.0]
      (if (> i n)
        (* acc dt)
        (let [w (if (or (= i 0) (= i n)) 0.5 1.0)
              p (v3/add-scaled origin dir (* dt (double i)))]
          (recur (inc i) (+ acc (* w (density atm p)))))))))

(defn transmittance
  "What fraction of each channel survives an optical depth.

  Beer's law, per channel, with Mie extinction folded in -- it is gray, so
  it darkens without coloring."
  [{:keys [rayleigh ^double mie]} ^double depth]
  (mapv (fn [b] (math/exp (- (* (+ (double b) mie) depth)))) rayleigh))

(defn rayleigh-phase
  "How much Rayleigh scattering goes at an angle: `3/(16pi) (1 + cos^2)`.

  Symmetric forward and backward, and only twice as strong along the
  axis as across it -- which is why the whole sky glows rather than just
  the part near the sun."
  ^double [^double cos-theta]
  (* (/ 3.0 (* 16.0 math/PI)) (+ 1.0 (* cos-theta cos-theta))))

(defn mie-phase
  "Henyey-Greenstein, the standard stand-in for the Mie phase function.

  `g` is how forward-biased it is: 0 is isotropic, and 0.76 puts most of
  the scattering within a few degrees of the sun's direction, which is
  the bright halo you cannot look at."
  ^double [^double g ^double cos-theta]
  (let [g2 (* g g)]
    (/ (* (/ 3.0 (* 8.0 math/PI)) (- 1.0 g2) (+ 1.0 (* cos-theta cos-theta)))
       (* (+ 2.0 g2)
          (math/pow (max 1e-9 (- (+ 1.0 g2) (* 2.0 g cos-theta))) 1.5)))))

(defn in-scatter
  "Light scattered into the view ray along `[t-near t-far]`, per channel.

  The single-scattering integral, as directly as it can be written: walk
  the segment, and at each sample ask how much sunlight reaches that
  point, how much of it turns toward the eye, and how much of that
  survives the rest of the way. Everything expensive about an atmosphere
  is in the inner sun-ray optical depth, which is why `:sun-steps` is
  separate from `:steps` and much smaller.

  Returns `[r g b]`."
  ([atm origin dir sun-dir t-near t-far] (in-scatter atm origin dir sun-dir t-near t-far {}))
  ([{:keys [rayleigh ^double mie ^double mie-g ^double top sun] :as atm}
    origin dir sun-dir t-near t-far
    {:keys [steps sun-steps] :or {steps 16 sun-steps 6}}]
   (let [t-near (double t-near)
         t-far (double t-far)
         dir (v3/normalize dir)
         sun-dir (v3/normalize sun-dir)
         cos-theta (v3/dot dir sun-dir)
         pr (rayleigh-phase cos-theta)
         pm (mie-phase mie-g cos-theta)
         n (max 1 (long steps))
         dt (/ (- t-far t-near) (double n))
         [sr sg sb] sun]
     (loop [i 0 view-depth 0.0 acc [0.0 0.0 0.0]]
       (if (= i n)
         (mapv * acc [sr sg sb])
         (let [p (v3/add-scaled origin dir (+ t-near (* dt (+ 0.5 (double i)))))
               d (density atm p)
               view-depth (+ view-depth (* d dt))
               ;; How far the sunlight traveled through air to get here.
               sun-depth (if-let [[_ far] (ray-sphere p sun-dir top)]
                           (optical-depth atm p sun-dir (max 0.0 (double far)) (long sun-steps))
                           0.0)
               total (+ view-depth sun-depth)
               [tr tg tb] (transmittance atm total)
               ;; Rayleigh is colored and Mie is not, so they cannot share
               ;; a coefficient even though they share everything else.
               contrib (fn [^double beta ^double t]
                         (* d dt t (+ (* beta pr) (* mie pm))))]
           (recur (inc i) view-depth
                  [(+ (nth acc 0) (contrib (double (nth rayleigh 0)) tr))
                   (+ (nth acc 1) (contrib (double (nth rayleigh 1)) tg))
                   (+ (nth acc 2) (contrib (double (nth rayleigh 2)) tb))])))))))

(defn sky-color
  "The color of the sky looking from `origin` along `dir`, or nil if that
  ray never enters the atmosphere.

  `:t-far` cuts the integral short where the ground is, which is what
  makes this do aerial perspective as well as sky: pass the distance to
  the terrain and you get the haze in front of it.

  The result is radiance, not a color: it is unbounded above, and the sun
  seen directly is orders of magnitude brighter than the sky beside it.
  Run it through `tone-map` before it reaches a pixel."
  ([atm origin dir sun-dir] (sky-color atm origin dir sun-dir {}))
  ([{:keys [^double top] :as atm} origin dir sun-dir {:keys [t-far] :as opts}]
   (let [dir (v3/normalize dir)]
     (when-let [[near far] (ray-sphere origin dir top)]
       (let [t0 (max 0.0 (double near))
             t1 (if t-far (min (double t-far) (double far)) (double far))]
         (when (> t1 t0)
           (in-scatter atm origin dir sun-dir t0 t1 opts)))))))

(defn aerial-perspective
  "`surface-color` as seen through the air between `origin` and `point`.

  The two halves of what distance does to a color, and they pull in
  opposite directions: the air takes light out of the view ray, and it
  puts scattered sunlight back in. Near enough, nothing happens; far
  enough, everything you see is the second term and the surface has
  vanished into the haze. Getting this right is most of the difference
  between a landscape that has a sense of scale and one that does not."
  ([atm origin point sun-dir surface-color]
   (aerial-perspective atm origin point sun-dir surface-color {}))
  ([atm origin point sun-dir surface-color {:keys [steps] :or {steps 12} :as opts}]
   (let [seg (v3/sub point origin)
         len (v3/length seg)]
     (if (zero? len)
       surface-color
       (let [dir (v3/scale seg (/ 1.0 len))
             depth (optical-depth atm origin dir len (long steps))
             t (transmittance atm depth)
             haze (or (sky-color atm origin dir sun-dir (assoc opts :t-far len))
                      [0.0 0.0 0.0])]
         (mapv (fn [c tc h] (+ (* (double c) (double tc)) (double h)))
               surface-color t haze))))))

(defn tone-map
  "Radiance to something a screen can show: exposure, then gamma.

  `1 - exp(-exposure * c)` rolls the highlights off smoothly instead of
  clipping them, which is what keeps the sun from becoming a flat white
  disk with a hard edge; the 1/2.2 afterward is the sRGB transfer curve,
  near enough. Not from chapter 18 -- the chapter computes radiance and
  stops, quite properly -- but every use of this namespace needs it, and
  leaving it out is how atmospheric renderings end up looking burned."
  ([color] (tone-map 1.0 color))
  ([exposure color]
   (let [e (double exposure)]
     (mapv (fn [c] (math/pow (- 1.0 (math/exp (* (- e) (max 0.0 (double c))))) (/ 1.0 2.2)))
           color))))
