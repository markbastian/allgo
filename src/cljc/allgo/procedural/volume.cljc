(ns allgo.procedural.volume
  "Things with no surface: steam, cloud, fire, and solids gone soft at the
  edges, as densities in space, and a ray marcher to see them by.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*. Perlin's
  chapter 12 for hypertexture -- a shape given as a density that falls
  from 1 inside to 0 outside, with noise added before it is sharpened, so
  that the boundary itself becomes the texture. Ebert's chapters 6 to 8
  for gases: steam and smoke as turbulence shaped by where the gas can
  be, and clouds -- cumulus, cirrus and stratus -- as his two levels, a
  coarse shape (blobs, streaks or a layer) and turbulence for everything
  smaller.

  A volume here is a map:

      {:density  (fn [x y z] -> d)     how much is there, 0 for nothing
       :extinction k                   how opaque a unit of d is, per unit length
       :albedo   [r g b]               what fraction of light it scatters
       :emission (fn [d] -> [r g b])   light it gives off itself, or nil
       :surface  s}                    0 for a gas, 1 for something solid

  and everything stays inside the unit ball, which is what the marcher
  bounds its rays by.

  `:surface` decides how a sample is lit. A gas is lit only by how much
  light reaches it through the gas between it and the lamp -- that is all
  that shapes a cloud. A hypertexture is a solid with a soft edge, and
  reads as one only if it is also lit like a surface, by the angle
  between the lamp and the way its density falls off, which is Perlin's
  normal: the gradient."
  (:require [allgo.math :as am]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as s]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn- len ^double [^double x ^double y ^double z] (math/sqrt (+ (* x x) (* y y) (* z z))))

(defn- turbulence
  "Perlin's turbulence over gradient noise, scaled to land near [0, 1]."
  [seed octaves]
  (let [f (fractal/turbulence (noise/gradient-basis {:seed seed}) {:octaves octaves})]
    (fn ^double [x y z] (* 1.4 (double (f (double x) (double y) (double z)))))))

;; ---------------------------------------------------------------------------
;; Chapter 12: hypertexture

(defn noisy-sphere
  "A ball whose surface is noise.

  Perlin's first hypertexture. The ball is a density falling linearly
  from 1 at the center to 0 at the edge; fBm is added to it, and `gain`
  sharpens the sum, so what was a soft edge becomes a surface again --
  but one wherever the noisy density crosses a half, which is lumpy at
  every scale. `:roughness` is how much noise; `:sharpness` how hard the
  result is, from fog at 0.5 to rock near 1."
  ([] (noisy-sphere {}))
  ([{:keys [seed roughness sharpness frequency]
     :or {seed 0 roughness 0.35 sharpness 0.95 frequency 2.5}}]
   (let [f (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves 5})
         roughness (double roughness) sharpness (double sharpness) k (double frequency)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       base (- 1.0 (/ (len x y z) 0.66))]
                   (s/gain sharpness (am/clamp (+ 0.5 base (* roughness (double (f (* k x) (* k y) (* k z)))))
                                               0.0 1.0))))
      :extinction 15.0
      :albedo [0.70 0.62 0.52]
      :surface 1.0})))

(defn eroded-cube
  "A cube worn away: hypertexture on a box rather than a ball.

  The same construction as `noisy-sphere`, on the density of a cube --
  1 in the middle, falling to 0 at its faces -- with turbulence in place
  of fBm, whose creases bite pits and channels into the faces and round
  off the corners first, as weather does."
  ([] (eroded-cube {}))
  ([{:keys [seed erosion] :or {seed 0 erosion 0.3}}]
   (let [t (turbulence seed 5) erosion (double erosion)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       base (- 1.0 (/ (max (abs x) (abs y) (abs z)) 0.5))]
                   (s/gain 0.97 (am/clamp (+ 0.6 base (- (* erosion (t (* 2.5 x) (* 2.5 y) (* 2.5 z)))))
                                          0.0 1.0))))
      :extinction 15.0
      :albedo [0.62 0.60 0.56]
      :surface 1.0})))

(defn- fire-color
  "Black-body-ish: red where the fire is thin, through orange and yellow
  to white where it is dense."
  [d]
  (let [d (am/clamp (double d) 0.0 1.0)]
    [(* 6.0 (s/smoothstep 0.0 0.35 d))
     (* 4.0 (s/smoothstep 0.2 0.75 d))
     (* 3.0 (s/smoothstep 0.6 1.0 d))]))

(defn fireball
  "A ball of flame, pushed out ragged by turbulence and burning brighter
  toward its core.

  Perlin's flame, made solid: the distance from the center has
  turbulence added to it before it is read, so the edge is torn, and the
  turbulence drifts upward with `:time`, which is how the flames climb.
  It gives off light rather than scattering it: the color is a ramp on
  the density, from dull red at the fringe to white at the heart."
  ([] (fireball {}))
  ([{:keys [seed time] :or {seed 0 time 0.0}}]
   (let [t (turbulence seed 5) time (double time)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       r (len x y z)
                       tu (t (* 2.5 x) (- (* 2.5 y) (* 1.5 time)) (* 2.5 z))
                       ;; The flames lean upward: the top of the ball is
                       ;; torn further out than the bottom.
                       q (- r (* tu (+ 0.25 (* 0.2 (max 0.0 y)))))]
                   (- 1.0 (s/smoothstep 0.05 0.35 q))))
      :extinction 6.0
      :albedo [0.0 0.0 0.0]
      :emission fire-color
      :surface 0.0})))

;; ---------------------------------------------------------------------------
;; Chapters 6 to 8: gases

(defn steam
  "Steam rising from a cup: turbulence, held to a column that widens,
  thins and twists as it rises.

  Ebert's steam. The gas itself is turbulence; what makes it steam is
  everything multiplied onto it -- falling off with distance from the
  column, which spreads as it rises, and with height, as the steam
  disperses -- and the helix it is turned through, more the higher it
  goes. `:time` carries the turbulence upward through the column, so the
  steam is animated by moving the noise, not the gas: Ebert's chapter 6."
  ([] (steam {}))
  ([{:keys [seed time twist] :or {seed 0 time 0.0 twist 2.5}}]
   (let [t (turbulence seed 5) time (double time) twist (double twist)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       h (/ (+ y 0.85) 1.7)]
                   (if (or (neg? h) (> h 1.0))
                     0.0
                     (let [a (* twist h)
                           c (math/cos a) sn (math/sin a)
                           u (- (* c x) (* sn z)) w (+ (* sn x) (* c z))
                           spread (+ 0.14 (* 0.4 h))
                           r (/ (len u 0.0 w) spread)
                           column (- 1.0 (s/smoothstep 0.4 1.0 r))
                           gas (t (* 4.0 u) (- (* 3.0 y) (* 2.0 time)) (* 4.0 w))]
                       (* column (math/pow (- 1.0 h) 1.5) (s/smoothstep 0.25 0.9 gas)
                          (s/smoothstep 0.0 0.08 h))))))
      :extinction 22.0
      :albedo [0.92 0.93 0.95]
      :surface 0.0})))

(defn- wyvill
  "Wyvill's soft-object falloff, of the squared distance over the squared
  radius: 1 at the center, 0 from the radius out, smooth at both ends."
  ^double [^double r2]
  (if (>= r2 1.0)
    0.0
    (let [r4 (* r2 r2) r6 (* r4 r2)]
      (+ 1.0 (* (/ -4.0 9.0) r6) (* (/ 17.0 9.0) r4) (* (/ -22.0 9.0) r2)))))

(defn cumulus
  "A fair-weather cloud: Ebert's two levels.

  The shape is implicit: a cluster of soft blobs, piled up from a flat
  base, their fields summed so they merge where they overlap. The detail
  is procedural: the point is moved by turbulence before the blobs are
  asked about it, which tears the smooth outline into billows, and more
  turbulence thins the inside. Neither level alone looks like a cloud --
  blobs are cotton wool, turbulence is fog.

  `:time` boils it: the noise of both levels is read further down each
  instant, so billows roll up through a cloud that stays where it is."
  ([] (cumulus {}))
  ([{:keys [seed blobs billow time] :or {seed 0 blobs 14 billow 0.22 time 0.0}}]
   (let [rng (random/rng seed 41 0 0)
         blobs (vec (for [i (range blobs)]
                      (let [a (* 2.0 math/PI (random/uniform rng 1.0))
                            d (* (random/uniform rng 0.0 0.6) (- 1.0 (/ i (double blobs) 2.0)))
                            r (random/uniform rng 0.3 0.45)]
                        [(* d (math/cos a))
                         (+ -0.25 r (* 0.25 (random/uniform rng 1.0)))
                         (* 0.7 d (math/sin a))
                         r])))
         warp (noise/vector-basis {:seed (+ seed 1)})
         t (turbulence (+ seed 2) 4)
         billow (double billow) time (double time)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       ;; Up through the noise: the billows rise.
                       ny (- y (* 0.3 time))
                       ;; Two octaves of warp: the billows, and the billows on them.
                       [wx wy wz] (warp (* 2.5 x) (* 2.5 ny) (* 2.5 z))
                       [vx vy vz] (warp (* 7.0 x) (* 7.0 ny) (* 7.0 z))
                       px (+ x (* billow (+ (double wx) (* 0.4 (double vx)))))
                       py (+ y (* billow (+ (double wy) (* 0.4 (double vy)))))
                       pz (+ z (* billow (+ (double wz) (* 0.4 (double vz)))))
                       field (reduce (fn [acc [cx cy cz r]]
                                       (let [dx (- px cx) dy (- py cy) dz (- pz cz)
                                             r (double r)]
                                         (+ acc (wyvill (/ (+ (* dx dx) (* dy dy) (* dz dz)) (* r r))))))
                                     0.0 blobs)]
                   (if (< field 0.15)
                     0.0
                     ;; Turbulence eats into the field before it is
                     ;; thresholded, so the edge breaks into puffs, and
                     ;; thins what is left.
                     (let [e (t (* 6.0 x) (* 6.0 ny) (* 6.0 z))]
                       (* (s/smoothstep 0.3 0.55 (- field (* 0.3 e)))
                          ;; The flat base, where the air stops rising.
                          (s/smoothstep -0.32 -0.22 y)
                          (+ 0.6 (* 0.4 e)))))))
      :extinction 20.0
      :albedo [0.96 0.96 0.98]
      :surface 0.0})))

(defn cirrus
  "Mares' tails: a few thin streaks high up, combed out by the wind.

  Ebert's two levels again, drawn out flat. The shape is a handful of
  long, thin implicit streaks lying along the wind, each bending as it
  goes and hooked down at its end, where its ice falls out of it. The
  detail is turbulence stretched along the wind -- slow along it, fast
  across it and up and down -- so that what is inside a streak is
  fibers rather than billows. And it is thin: little density and little
  extinction, so the sky shows through and the sun lights all of it
  nearly alike, which is why cirrus is white and flat where cumulus is
  modeled in light and shade. `:wind` is the direction the streaks run,
  in radians about the vertical; `:time` blows the fibers along them."
  ([] (cirrus {}))
  ([{:keys [seed streaks wind time] :or {seed 0 streaks 6 wind -0.5 time 0.0}}]
   (let [rng (random/rng seed 43 0 0)
         streaks (vec (for [_ (range streaks)]
                        {:angle (random/uniform rng -0.25 0.25)
                         :offset (random/uniform rng -0.7 0.7)
                         :along (random/uniform rng -0.2 0.2)
                         :height (random/uniform rng -0.05 0.25)
                         :length (random/uniform rng 0.5 0.8)
                         :bend (random/uniform rng -0.25 0.25)
                         :hook (random/uniform rng 0.08 0.18)
                         :width (random/uniform rng 0.07 0.12)}))
         t (turbulence (+ seed 3) 5)
         warp (noise/vector-basis {:seed (+ seed 4)})
         wind (double wind) time (double time)
         cw (math/cos wind) sw (math/sin wind)]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       ;; Along the wind, and across it.
                       a (+ (* cw x) (* sw z))
                       c (- (* cw z) (* sw x))
                       [wx wy _] (warp (* 1.5 x) (* 1.5 y) (* 1.5 z))
                       c (+ c (* 0.06 (double wx)))
                       y (+ y (* 0.02 (double wy)))
                       field (reduce
                              (fn [^double acc {:keys [angle offset along height length bend hook width]}]
                                (let [angle (double angle)
                                      ca (math/cos angle) sa (math/sin angle)
                                      ;; s along the streak, -1 to 1; n across it.
                                      s (/ (- (+ (* ca a) (* sa c)) (double along)) (double length))
                                      n (- (- (* ca c) (* sa a)) (double offset) (* (double bend) s s))
                                      ;; The fall streak: the downwind end sags.
                                      dy (+ (- y (double height)) (* (double hook) (s/smoothstep 0.1 1.0 s)))
                                      w (double width)
                                      th (* 0.35 w)]
                                  (+ acc (wyvill (+ (* s s) (/ (* n n) (* w w)) (/ (* dy dy) (* th th)))))))
                              0.0 streaks)]
                   (if (< field 0.02)
                     0.0
                     (let [fibers (t (* 1.2 (- a time)) (* 14.0 y) (* 14.0 c))]
                       (* (s/smoothstep 0.02 0.4 field)
                          (s/smoothstep 0.3 0.8 fibers)
                          (- 1.0 (s/smoothstep 0.8 0.95 (len x y z))))))))
      :extinction 5.0
      :albedo [0.97 0.97 1.0]
      :surface 0.0})))

(defn stratus
  "An overcast: a sheet of cloud, flat underneath, rolling on top, and
  broken here and there into patches.

  Ebert's stratus is the cloud with the least shape -- a layer, so the
  macrostructure is a plane and no blobs, and everything else is noise.
  Low-frequency fBm across the sheet decides where there is cloud
  (`:cover`, roughly the fraction of the sky it covers) and how thick it
  is there, so that it thins to nothing at a gap rather than ending in a
  wall; turbulence gives the top its rolls and the whole its mottling.
  Lit from above through itself, the sheet is pale on top and gray
  underneath, as an overcast is. `:time` drifts it."
  ([] (stratus {}))
  ([{:keys [seed cover time] :or {seed 0 cover 0.8 time 0.0}}]
   (let [patch (fractal/fbm (noise/gradient-basis {:seed (+ seed 5)}) {:octaves 5})
         t (turbulence (+ seed 7) 4)
         ;; fBm here spreads about zero with a deviation near a quarter,
         ;; so this is the level about `cover` of the sheet lies above.
         edge (- (* -0.8 (- (double cover) 0.5)) 0.13)
         time (double time)
         base -0.1
         thickness 0.11]
     {:density (fn [x y z]
                 (let [x (double x) y (double y) z (double z)
                       r (len x 0.0 z)]
                   (if (or (> r 0.97) (< y (- base 0.03)) (> y (+ base thickness 0.05)))
                     0.0
                     (let [px (- x (* 0.15 time))
                           p (double (patch (* 1.4 px) 0.0 (* 1.4 z)))
                           ;; Ragged at the rim instead of a hard circle.
                           p (- p (* 0.8 (s/smoothstep 0.7 0.97 r)))
                           cover (s/smoothstep edge (+ edge 0.2) p)]
                       (if (<= cover 0.0)
                         0.0
                         (let [e (t (* 9.0 px) (* 9.0 y) (* 9.0 z))
                               top (+ base (* thickness cover (+ 0.8 (* 0.3 e))))]
                           (* (s/smoothstep (- base 0.03) (+ base 0.02) y)
                              (- 1.0 (s/smoothstep (- top 0.06) top y))
                              (+ 0.65 (* 0.35 e)))))))))
      :extinction 18.0
      :albedo [0.9 0.91 0.93]
      :surface 0.0})))

(def catalog
  "Every volume here, by name. `:moves?` marks the ones that `:time`
  animates; the two hypertextures are solids, and keep still."
  [{:name "noisy sphere" :make noisy-sphere :source "chapter 12"}
   {:name "eroded cube" :make eroded-cube :source "chapter 12"}
   {:name "fireball" :make fireball :source "chapter 12" :moves? true}
   {:name "steam" :make steam :source "chapters 5, 6" :moves? true}
   {:name "cumulus" :make cumulus :source "chapter 8" :moves? true}
   {:name "cirrus" :make cirrus :source "chapter 8" :moves? true}
   {:name "stratus" :make stratus :source "chapter 8" :moves? true}])

;; ---------------------------------------------------------------------------
;; Rendering

(defn camera-ray
  "The ray `[origin direction]` through `(u, v)` in [-1, 1] on the
  picture, from a camera `distance` from the origin at `azimuth` and
  `elevation` radians, looking at it, `fov` radians across."
  [azimuth elevation distance fov u v]
  (let [az (double azimuth) el (double elevation) d (double distance)
        ce (math/cos el)
        ox (* d ce (math/sin az)) oy (* d (math/sin el)) oz (* d ce (math/cos az))
        ;; Forward, right, up.
        fx (/ (- ox) d) fy (/ (- oy) d) fz (/ (- oz) d)
        rl (len fz 0.0 (- fx))
        rx (/ (- fz) rl) rz (/ fx rl)
        ux (- (* fy rz)) uy (- (* rz fx) (* rx fz)) uz (* fy rx)
        k (math/tan (* 0.5 (double fov)))
        dx (+ fx (* k (+ (* (double u) rx) (* (double v) ux))))
        dy (+ fy (* k (double v) uy))
        dz (+ fz (* k (+ (* (double u) rz) (* (double v) uz))))
        dl (len dx dy dz)]
    [[ox oy oz] [(/ dx dl) (/ dy dl) (/ dz dl)]]))

(defn- ball-span
  "Where a ray enters and leaves the unit ball, or nil if it misses."
  [[ox oy oz] [dx dy dz]]
  (let [ox (double ox) oy (double oy) oz (double oz)
        b (+ (* ox (double dx)) (* oy (double dy)) (* oz (double dz)))
        c (- (+ (* ox ox) (* oy oy) (* oz oz)) 1.0)
        disc (- (* b b) c)]
    (when (pos? disc)
      (let [q (math/sqrt disc)]
        [(max 0.0 (- (- b) q)) (- q b)]))))

(defn pixel-jitter
  "A number in [0, 1) for the pixel `(px, py)`, unrelated to its
  neighbors': a hash of the two, so white noise. Sequences that are more
  even than that -- the golden ratio along a diagonal, Roberts' R2 --
  leave their own lattice in a volume's bands."
  ^double [px py]
  (let [h (* 43758.5453 (math/sin (+ (* 12.9898 (double px)) (* 78.233 (double py)))))]
    (- h (math/floor h))))

(def default-light
  "Up, to the right and toward the viewer, so that what faces us is lit
  and the undersides fall into shadow."
  (let [[x y z] [0.6 0.7 0.4] l (len x y z)] [(/ x l) (/ y l) (/ z l)]))

(defn march
  "The color `[r g b]` seen along the ray `[origin direction]` through
  `volume`, over `:background`.

  Front to back: at each of `:steps` samples across the unit ball, the
  slab of medium there absorbs a share of what is behind it and adds its
  own -- scattered light, and emitted if it burns -- and the march stops
  once what is in front is opaque. `:jitter` in [0, 1) shifts where the
  samples fall, which turns the bands that regular steps leave into
  noise the eye forgives -- so long as it differs from pixel to pixel
  without a pattern: see `pixel-jitter`.

  Light reaching a sample is the sun through `:shadow-steps` samples of
  the medium toward it, over `:shadow-length`, or all of it if that is 0,
  plus `:ambient` from the sky. A dense medium is dark a short way in, so
  a short shadow march finely sampled sees more than a long coarse one,
  whose few samples land in the billows or between them by chance and
  speckle the picture.

  A volume with `:surface` is solid where its density passes a half.
  Short of that it is marched as a gas, which is its soft edge; at the
  first sample past it, the crossing is found by bisection and lit as a
  surface, by Lambert's law on the density's gradient there. Lit at
  whichever sample happened to land past the half instead, its surface
  moved by up to a step from pixel to pixel, and the gradients taken
  there turned it to fur."
  [{:keys [density extinction albedo emission surface]} [origin direction]
   {:keys [steps shadow-steps shadow-length light sun ambient background jitter]
    :or {steps 64 shadow-steps 6 shadow-length 0.5 light default-light sun 1.0 ambient 0.35
         background [0.0 0.0 0.0] jitter 0.5}}]
  (let [[br bg bb] background]
    (if-let [[t0 t1] (ball-span origin direction)]
      (let [[ox oy oz] origin [dx dy dz] direction [lx ly lz] light
            ox (double ox) oy (double oy) oz (double oz)
            dx (double dx) dy (double dy) dz (double dz)
            lx (double lx) ly (double ly) lz (double lz)
            ext (double extinction) surface (double (or surface 0.0))
            [ar ag ab] albedo
            ar (double ar) ag (double ag) ab (double ab)
            dt (/ 2.0 (double steps))
            sdt (/ (double shadow-length) (max 1.0 (double shadow-steps)))
            sun (double sun) ambient (double ambient)
            density (fn ^double [x y z] (double (density x y z)))
            shadowed (fn ^double [^double px ^double py ^double pz]
                       (if (pos? (long shadow-steps))
                         (loop [i 1 acc 0.0]
                           (if (> i (long shadow-steps))
                             (math/exp (- (* ext acc sdt)))
                             (let [s (* (- i 0.5) sdt)
                                   qx (+ px (* s lx)) qy (+ py (* s ly)) qz (+ pz (* s lz))]
                               (recur (inc i)
                                      (if (< (+ (* qx qx) (* qy qy) (* qz qz)) 1.0)
                                        (+ acc (density qx qy qz))
                                        acc)))))
                         1.0))
            solid? (pos? surface)
            hit (fn [^double ta ^double tb]
                  ;; The crossing of a half between ta, short of it, and tb,
                  ;; past it; its color as a lit surface.
                  (let [^double t (loop [i 0 ta ta tb tb]
                                    (let [tm (* 0.5 (+ ta tb))]
                                      (if (= i 10)
                                        tm
                                        (if (< (density (+ ox (* tm dx)) (+ oy (* tm dy)) (+ oz (* tm dz))) 0.5)
                                          (recur (inc i) tm tb)
                                          (recur (inc i) ta tm)))))
                        px (+ ox (* t dx)) py (+ oy (* t dy)) pz (+ oz (* t dz))
                        e 0.004
                        gx (- (density (+ px e) py pz) (density (- px e) py pz))
                        gy (- (density px (+ py e) pz) (density px (- py e) pz))
                        gz (- (density px py (+ pz e)) (density px py (- pz e)))
                        gl (len gx gy gz)
                        [nx ny nz] (if (pos? gl) [(/ (- gx) gl) (/ (- gy) gl) (/ (- gz) gl)] [0.0 0.0 0.0])
                        nx (double nx) ny (double ny) nz (double nz)
                        lambert (max 0.0 (+ (* nx lx) (* ny ly) (* nz lz)))
                        lambert (s/mix 1.0 lambert surface)
                        ;; Off the surface a little before looking for the
                        ;; sun, or it shadows itself.
                        sh (if (pos? lambert) (shadowed (+ px (* 0.02 nx)) (+ py (* 0.02 ny)) (+ pz (* 0.02 nz))) 0.0)
                        l (+ ambient (* sun sh lambert))]
                    [(* ar l) (* ag l) (* ab l)]))]
        (loop [t (+ (double t0) (* dt (double jitter))) tr 1.0 r 0.0 g 0.0 b 0.0]
          (if (or (> t (double t1)) (< tr 0.01))
            [(+ r (* tr (double br))) (+ g (* tr (double bg))) (+ b (* tr (double bb)))]
            (let [px (+ ox (* t dx)) py (+ oy (* t dy)) pz (+ oz (* t dz))
                  d (density px py pz)]
              (cond
                (<= d 1e-4)
                (recur (+ t dt) tr r g b)

                (and solid? (>= d 0.5))
                (let [[sr sg sb] (hit (max (double t0) (- t dt)) t)]
                  [(+ r (* tr (double sr))) (+ g (* tr (double sg))) (+ b (* tr (double sb)))])

                :else
                (let [a (- 1.0 (math/exp (- (* ext d dt))))
                      l (if (or (pos? ar) (pos? ag) (pos? ab))
                          (+ ambient (* sun (shadowed px py pz)))
                          0.0)
                      [er eg eb] (if emission (emission d) [0.0 0.0 0.0])
                      w (* tr a)]
                  (recur (+ t dt) (* tr (- 1.0 a))
                         (+ r (* w (+ (* ar l) (double er))))
                         (+ g (* w (+ (* ag l) (double eg))))
                         (+ b (* w (+ (* ab l) (double eb)))))))))))
      [(double br) (double bg) (double bb)])))

(defn tone
  "A color brought into [0, 1]: 1 - e^-c a channel, so that fire, which
  emits well past 1, saturates to white the way film does instead of
  clipping to flat yellow."
  [[r g b]]
  [(- 1.0 (math/exp (- (* 1.6 (double r)))))
   (- 1.0 (math/exp (- (* 1.6 (double g)))))
   (- 1.0 (math/exp (- (* 1.6 (double b)))))])
