(ns allgo.procedural.texture
  "Procedural materials: brick, marble, wood, stone, metal, as functions of
  a point in space.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*. Peachey's
  chapter 2 for the method and for brick, the book's first example, bumps
  and all; Worley's chapters 3 and 4 for building materials out of layers
  and out of cells; Musgrave's chapter 15, \"Fractal Solid Textures\", for
  marble, clouds and cirrus, water, an earth, a gas giant's bands, random
  fBm coloring and the impressionist filter that paints any of the others.
  The rest -- the metals, wood, granite -- are the same parts put together
  the same way.

  What every one of them is:

      (fn [x y z] -> [r g b height shine])

  a color in [0, 1], a height for the bumps, and how glossy the surface is
  there, 0 for brick and 1 for polished steel. It is a *solid* texture: it
  is defined everywhere in space, so an object takes it on by being
  evaluated where its surface is, and carves its pattern out of the solid
  the way a statue is carved out of a block of marble, with no seams and
  nothing stretched. A plane at z = 0 shows the face of the block.

  Height is what sells a material. A brick wall drawn flat is a picture of
  bricks; with the mortar sunk and the brick edges rounded it is a wall,
  and the light across it says so. `shade` lights a texture from its own
  height, by differences, which is Peachey's bump-mapped brick: the
  geometry stays flat and only the normal is bent.

  Units are the texture's own. A brick is about a unit long and a stone
  about a unit across, and a caller scales the point to size them.

  Every constructor takes an options map, with `:seed` among the options,
  and builds its noise once: call the constructor once, the texture a
  million times."
  (:require [allgo.math :as am]
            [allgo.procedural.cellular :as cellular]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as s]
            [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Parts

(defn color-map
  "`(fn [t] -> [r g b])` through `knots`, a seq of colors, with `t` in
  [0, 1] spanning them.

  Peachey's color spline: Catmull-Rom per channel, so the ramp passes
  through every color given and turns smoothly between them. As with
  `shaping/spline`, the first and last knots shape the ends and are not
  reached, so repeat an end color to land on it."
  [knots]
  (let [cs (mapv (fn [c] (mapv #(double (nth % c)) knots)) [0 1 2])
        [rs gs bs] cs]
    (fn [t]
      (let [t (am/clamp (double t) 0.0 1.0)]
        [(am/clamp (s/spline t rs) 0.0 1.0)
         (am/clamp (s/spline t gs) 0.0 1.0)
         (am/clamp (s/spline t bs) 0.0 1.0)]))))

(defn- mixc
  "Color `a` toward color `b` by `t`."
  [[ar ag ab] [br bg bb] t]
  (let [t (double t)]
    [(s/mix ar br t) (s/mix ag bg t) (s/mix ab bb t)]))

(defn- scalec [[r g b] k]
  (let [k (double k)] [(* r k) (* g k) (* b k)]))

(defn- texel
  "A color, a height and a shine, as a texture returns them."
  [[r g b] h shine]
  [(am/clamp r 0.0 1.0) (am/clamp g 0.0 1.0) (am/clamp b 0.0 1.0)
   (double h) (am/clamp shine 0.0 1.0)])

(defn- cell-random
  "A repeatable random number in [0, 1) for the integer cell `[i j]`: one
  per brick, so each can be its own shade."
  ^double [seed i j]
  (random/unit (random/advance (random/cell-seed (long seed) (long i) (long j) 7))))

(defn- fbm
  "fBm over gradient noise, from `seed`, mapped onto roughly [0, 1]."
  [seed octaves]
  (let [f (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves octaves})]
    (fn ^double [x y z] (+ 0.5 (* 0.7 (double (f (double x) (double y) (double z))))))))

(defn- cells
  "The nearest two distances and the nearest point's value, as a
  function, at `density` points a unit cell."
  [seed density]
  (cellular/nearest {:seed seed :density density :n 2}))

;; ---------------------------------------------------------------------------
;; Chapter 2: brick

(defn brick
  "Running bond: courses of bricks, each offset half a brick from the
  last, in mortar.

  Peachey's brick, which is where the book begins. Which brick a point is
  in is two floors; how far it is from that brick's edge says brick or
  mortar, smoothed over a short width rather than stepped, so the edge
  can be lit. Each brick takes its own shade from its cell, and a coarse
  noise over the whole wall is the grime that ties them together.

  `:width` and `:height` are a brick's, mortar included, and `:mortar` is
  the joint's width."
  ([] (brick {}))
  ([{:keys [seed width height mortar color dark mortar-color]
     :or {seed 0 width 1.0 height 0.36 mortar 0.05
          color [0.62 0.25 0.16] dark [0.36 0.14 0.10] mortar-color [0.72 0.70 0.65]}}]
   (let [width (double width) height (double height) half (* 0.5 (double mortar))
         grime (fbm seed 5)
         grain (noise/gradient-basis {:seed (+ seed 1)})]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             row (math/floor (/ y height))
             xs (if (odd? (long row)) (+ x (* 0.5 width)) x)
             col (math/floor (/ xs width))
             u (- (/ xs width) col)
             v (- (/ y height) row)
             edge (min (* width (min u (- 1.0 u))) (* height (min v (- 1.0 v))))
             ;; 0 in the joint, 1 on the face, and a rounded shoulder
             ;; between, which is what catches the light.
             face (s/smoothstep half (+ half 0.03) edge)
             r (cell-random seed col row)
             g (grime x y z)
             sand (* 0.5 (+ 1.0 (double (grain (* 40.0 x) (* 40.0 y) (* 40.0 z)))))
             brick-color (scalec (mixc color dark (* r r)) (+ 0.75 (* 0.45 g) (* 0.12 sand)))
             mortar (scalec mortar-color (+ 0.8 (* 0.25 sand)))]
         (texel (mixc mortar brick-color face)
                (+ (* 0.6 face) (* 0.05 sand) (* 0.15 face g))
                (* 0.08 face)))))))

;; ---------------------------------------------------------------------------
;; Chapters 2 and 15: fractal solid textures

(defn marble
  "Veins through stone: a sine across the block, pushed about by
  turbulence, read through a color spline.

  The bands are where the sine crosses zero, and turbulence is what makes
  them wander and branch the way veins do; with no turbulence it is
  stripes. `:veining` is how far they wander."
  ([] (marble {}))
  ([{:keys [seed veining frequency]
     :or {seed 0 veining 0.9 frequency 4.0}}]
   (let [turb (fractal/turbulence (noise/gradient-basis {:seed seed}) {:octaves 6})
         veining (double veining) frequency (double frequency)
         ramp (color-map [[0.20 0.20 0.24] [0.20 0.20 0.24] [0.45 0.45 0.50]
                          [0.80 0.79 0.78] [0.93 0.92 0.90] [0.95 0.94 0.92]
                          [0.95 0.94 0.92]])]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             t (+ x (* 0.4 y) (* veining (double (turb x y z))))
             v (abs (math/sin (* frequency t)))]
         ;; `abs` folds the sine so both of its zero crossings are veins,
         ;; and `gain` keeps them thin.
         (texel (ramp (math/pow v 0.5)) 0.0 0.7))))))

(defn wood
  "Growth rings around a trunk running along x, and the grain along it.

  The rings are the distance from the trunk's axis, cut into periods;
  noise moves them about, which is the difference between wood and a
  target. The plane z = 0 lies `:offset` from the axis, so it cuts the
  rings the way a board is sawn from a log -- in long arches -- rather
  than through the pith."
  ([] (wood {}))
  ([{:keys [seed rings offset light dark]
     :or {seed 0 rings 6.0 offset 2.2
          light [0.80 0.58 0.35] dark [0.50 0.30 0.15]}}]
   (let [warp (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves 3})
         grain (noise/gradient-basis {:seed (+ seed 1)})
         rings (double rings) offset (double offset)]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             zz (+ z offset)
             r (math/sqrt (+ (* y y) (* zz zz)))
             r (+ r (* 0.12 (double (warp (* 0.3 x) (* 2.0 y) (* 2.0 zz)))))
             ring (s/sawtooth 1.0 (* rings r))
             ;; Early wood is pale and wide, late wood dark and narrow.
             late (s/smoothstep 0.65 0.95 ring)
             fine (* 0.5 (+ 1.0 (double (grain (* 0.6 x) (* 60.0 y) (* 60.0 zz)))))
             color (scalec (mixc light dark late) (+ 0.85 (* 0.2 fine)))]
         (texel color (* 0.08 (+ late fine)) 0.25))))))

(defn granite
  "Speckled stone: a field of grains, each crystal its own color.

  Three noises at grain size pick, in turn, dark mica, pale quartz and a
  pink feldspar ground. Thresholded noise is the book's way of making
  spots: noise is smooth, so cutting it at a level gives blobs of every
  size with smooth edges."
  ([] (granite {}))
  ([{:keys [seed grain] :or {seed 0 grain 30.0}}]
   (let [k (double grain)
         dark (fbm seed 3) pale (fbm (+ seed 1) 3) ground (fbm (+ seed 2) 4)]
     (fn [x y z]
       (let [x (* k (double x)) y (* k (double y)) z (* k (double z))
             g (ground (* 0.3 x) (* 0.3 y) (* 0.3 z))
             base (mixc [0.62 0.50 0.48] [0.80 0.70 0.66] g)
             q (s/smoothstep 0.66 0.70 (pale x y z))
             m (s/smoothstep 0.70 0.74 (dark (+ x 31.0) y z))
             color (-> base (mixc [0.92 0.92 0.92] q) (mixc [0.10 0.10 0.12] m))]
         (texel color 0.0 0.6))))))

(defn clouds
  "fBm as a sky: blue where the sum is low, white where it is high.

  Musgrave's first fractal solid texture and the plainest: the sum itself,
  with `:cover` sliding the threshold so the same noise gives a clear
  day or an overcast one."
  ([] (clouds {}))
  ([{:keys [seed cover] :or {seed 0 cover 0.5}}]
   (let [f (fbm seed 7) cover (double cover)]
     (fn [x y z]
       (let [d (f x y z)
             c (s/smoothstep (- 1.0 cover) (- 1.25 cover) d)]
         (texel (mixc [0.26 0.48 0.82] [0.97 0.97 0.98] c) 0.0 0.0))))))

(defn cirrus
  "High wisps: clouds whose noise is read through a distortion, then
  pulled out along x.

  Musgrave's distortion: look up a random vector at the point, move the
  point along it, and evaluate there. The same fBm becomes uneven --
  stretched in one place, crowded in the next -- which is how cirrus
  looks, combed out by the wind aloft."
  ([] (cirrus {}))
  ([{:keys [seed cover distortion] :or {seed 0 cover 0.5 distortion 1.2}}]
   (let [f (fractal/fbm (noise/distorted (noise/gradient-basis {:seed seed})
                                         (noise/vector-basis {:seed (+ seed 1)})
                                         distortion)
                        {:octaves 7})
         cover (double cover)]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             d (+ 0.5 (* 0.8 (double (f (* 0.25 x) (* 2.5 y) (* 2.5 z)))))
             c (s/smoothstep (- 1.1 cover) (- 1.45 cover) d)]
         (texel (mixc [0.22 0.45 0.80] [0.98 0.98 1.0] (* c c)) 0.0 0.0))))))

(defn gas-giant
  "Belts and zones: color by latitude, stirred.

  Musgrave's global circulation, as on Jupiter: the color is a function
  of height up the axis alone, so the planet is banded, and turbulence
  added to that height before it is read is the weather -- the bands
  wave, shear and eddy where they meet. `:bands` is how many there are a
  unit up the axis."
  ([] (gas-giant {}))
  ([{:keys [seed bands stir] :or {seed 0 bands 2.2 stir 0.12}}]
   (let [turb (fractal/fbm (noise/distorted (noise/gradient-basis {:seed seed})
                                            (noise/vector-basis {:seed (+ seed 1)})
                                            0.6)
                           {:octaves 6})
         bands (double bands) stir (double stir)
         ramp (color-map [[0.80 0.72 0.60] [0.80 0.72 0.60] [0.93 0.88 0.78]
                          [0.70 0.48 0.32] [0.88 0.80 0.66] [0.55 0.38 0.28]
                          [0.92 0.86 0.74] [0.75 0.58 0.42] [0.75 0.58 0.42]])]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             h (+ y (* stir (double (turb (* 1.5 x) (* 5.0 y) (* 1.5 z)))))]
         (texel (ramp (s/sawtooth 1.0 (* bands h))) 0.0 0.1))))))

(defn water
  "Open water: deep color, and waves that only the light reveals.

  The color hardly varies; all of it is in the height, an fBm of waves
  longer one way than the other, and the shine. Under `shade` the
  highlight breaks up across the slopes as sunlight does on a sea."
  ([] (water {}))
  ([{:keys [seed] :or {seed 0}}]
   (let [waves (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves 6 :H 0.6})]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             h (double (waves (* 1.5 x) (* 4.0 y) (* 4.0 z)))]
         (texel (mixc [0.04 0.20 0.34] [0.10 0.40 0.50] (+ 0.5 h)) (* 0.8 h) 0.8))))))

(defn earth
  "A world: fBm as altitude, the sea at a level, land colored by height,
  and ice toward the poles.

  The texture Musgrave grows into whole planets. Below `:sea-level` it is
  ocean, darker the deeper, flat and glossy; above it, a ramp from beach
  through green and brown to rock and snow. Ice comes in where `|y|`
  passes `:ice-line`, raggedly, and lower on high ground -- the poles are
  along y. `allgo.procedural.planet` is the same idea at full size."
  ([] (earth {}))
  ([{:keys [seed sea-level ice-line] :or {seed 0 sea-level 0.52 ice-line 1.25}}]
   (let [alt (fbm seed 8)
         sea (double sea-level) ice-line (double ice-line)
         land (color-map [[0.76 0.70 0.50] [0.76 0.70 0.50] [0.30 0.52 0.22]
                          [0.20 0.40 0.16] [0.45 0.40 0.25] [0.50 0.45 0.40]
                          [0.95 0.95 0.97] [0.95 0.95 0.97]])]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             a (alt (* 0.8 x) (* 0.8 y) (* 0.8 z))
             up (/ (- a sea) (- 1.0 sea))
             ice (s/smoothstep 0.0 0.08 (- (abs y) (- ice-line (* 0.5 (max 0.0 up)) (* 0.3 (- a 0.5)))))]
         (if (neg? up)
           (texel (mixc (mixc [0.02 0.08 0.28] [0.10 0.35 0.55] (s/smoothstep -0.5 0.0 up))
                        [0.92 0.95 0.98] ice)
                  0.0 (s/mix 0.6 0.1 ice))
           (texel (mixc (land (* 1.15 up)) [0.94 0.95 0.98] ice) (* 2.0 up) 0.05)))))))

(defn random-coloring
  "Three fBms, one a channel: color itself made fractal.

  Musgrave's random fBm coloring. Each channel wanders independently, so
  hues drift into one another at every scale, the way they do in a
  stained cloth or an oil slick, with nothing about it that suggests a
  palette was chosen."
  ([] (random-coloring {}))
  ([{:keys [seed] :or {seed 0}}]
   (let [[r g b] (map #(fbm (+ seed %) 6) [0 11 23])]
     (fn [x y z]
       (let [c (fn [f] (s/smoothstep 0.15 0.85 (f x y z)))]
         (texel [(c r) (c g) (c b)] 0.0 0.2))))))

(defn impressionist
  "Any texture, painted: each dab of paint takes its color from somewhere
  a little way off.

  Musgrave's impressionistic filter, which perturbs where an image is
  read by a random vector. Here the vector is held constant over each of
  Worley's cells, so the picture breaks into dabs, each a patch of the
  original moved by its own stroke, and the dabs keep the original's
  colors while losing its outlines. `:dab` is a dab's size and `:smear`
  how far it moves."
  ([texture] (impressionist texture {}))
  ([texture {:keys [seed dab smear] :or {seed 0 dab 0.18 smear 0.8}}]
   (let [k (/ 1.0 (double dab))
         near (cells seed 1.0)
         wobble (noise/vector-basis {:seed (+ seed 1)})
         smear (* (double smear) (double dab))]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             {[d1 d2] :distances [px py pz] :point} (near (* k x) (* k y) (* k z))
             [dx dy dz] (wobble px py pz)
             [r g b _ shine] (texture (+ x (* smear (double dx)))
                                      (+ y (* smear (double dy)))
                                      (+ z (* smear (double dz))))]
         ;; Dabs are a little raised, as thick paint is, toward their
         ;; middles and less at their edges.
         (texel [r g b] (* 0.15 (s/smoothstep 0.0 0.3 (- (double d2) (double d1)))) (* 0.5 (double shine))))))))

;; ---------------------------------------------------------------------------
;; Chapter 4: cells

(defn flagstone
  "Paving: irregular flat stones in a bed of mortar.

  Worley's F2 - F1 is zero on the boundary between two cells and grows
  away from it, so a threshold on it is mortar of an even width all
  round every stone. Each stone takes its color from the value its
  feature point carries, and a little fBm over that keeps it from looking
  painted."
  ([] (flagstone {}))
  ([{:keys [seed density gap] :or {seed 0 density 1.0 gap 0.07}}]
   (let [near (cells seed density) gap (double gap)
         grain (fbm (+ seed 1) 5)
         ramp (color-map [[0.42 0.40 0.38] [0.42 0.40 0.38] [0.55 0.50 0.42]
                          [0.48 0.47 0.48] [0.60 0.56 0.50] [0.60 0.56 0.50]])]
     (fn [x y z]
       (let [{[d1 d2] :distances v :value} (near x y z)
             face (s/smoothstep gap (+ gap 0.04) (- (double d2) (double d1)))
             g (grain (* 3.0 (double x)) (* 3.0 (double y)) (* 3.0 (double z)))
             stone (scalec (ramp (* 0.5 (+ 1.0 (double v)))) (+ 0.8 (* 0.35 g)))]
         (texel (mixc [0.30 0.28 0.25] stone face)
                (+ (* 0.6 face) (* 0.2 face g))
                (* 0.1 face)))))))

(defn cobblestone
  "Rounded stones set close: the same cells, but each one domed.

  F1, the distance to a cell's own point, is small in the middle of the
  cell and largest at its edge, which is exactly the shape of a worn
  cobble seen from above."
  ([] (cobblestone {}))
  ([{:keys [seed density] :or {seed 0 density 1.5}}]
   (let [near (cells seed density)
         grain (fbm (+ seed 1) 4)]
     (fn [x y z]
       (let [{[d1 d2] :distances v :value} (near x y z)
             d1 (double d1) d2 (double d2)
             face (s/smoothstep 0.02 0.12 (- d2 d1))
             dome (math/sqrt (max 0.0 (- 1.0 (math/pow (/ d1 (+ d1 (* 0.5 (- d2 d1)) 1e-9)) 2.0))))
             g (grain (* 4.0 (double x)) (* 4.0 (double y)) (* 4.0 (double z)))
             stone (scalec (mixc [0.45 0.43 0.40] [0.62 0.58 0.52] (* 0.5 (+ 1.0 (double v))))
                           (+ 0.75 (* 0.4 g)))]
         (texel (mixc [0.20 0.18 0.16] stone face) (* face (+ (* 0.7 dome) (* 0.1 g))) (* 0.15 face)))))))

(defn cracked-mud
  "Dried mud: plates split along their cell boundaries, curling up at the
  edges as they shrink.

  The cracks are a thin threshold on F2 - F1. The curl is the same
  quantity read near the crack instead of at it: a plate rises toward its
  rim before it drops into the gap."
  ([] (cracked-mud {}))
  ([{:keys [seed density] :or {seed 0 density 2.0}}]
   (let [near (cells seed density)
         grain (fbm (+ seed 1) 5)]
     (fn [x y z]
       (let [{[d1 d2] :distances v :value} (near x y z)
             e (- (double d2) (double d1))
             crack (s/smoothstep 0.015 0.05 e)
             curl (- 1.0 (s/smoothstep 0.05 0.4 e))
             g (grain (* 5.0 (double x)) (* 5.0 (double y)) (* 5.0 (double z)))
             mud (scalec (mixc [0.62 0.50 0.36] [0.72 0.60 0.45] (* 0.5 (+ 1.0 (double v))))
                         (+ 0.8 (* 0.3 g)))]
         (texel (mixc [0.18 0.13 0.09] mud crack)
                (* crack (+ 0.5 (* 0.25 curl) (* 0.1 g)))
                0.0))))))

(defn stained-glass
  "Leaded glass: a flat color a cell, and a dark lead line between.

  The same F2 - F1 as flagstone, with the stones made panes and given
  saturated colors from a ring of hues."
  ([] (stained-glass {}))
  ([{:keys [seed density] :or {seed 0 density 1.2}}]
   (let [near (cells seed density)
         ramp (color-map [[0.10 0.25 0.70] [0.10 0.25 0.70] [0.75 0.10 0.15]
                          [0.95 0.70 0.10] [0.10 0.55 0.25] [0.45 0.15 0.60]
                          [0.10 0.25 0.70] [0.10 0.25 0.70]])
         ripple (fbm (+ seed 1) 3)]
     (fn [x y z]
       (let [{[d1 d2] :distances v :value} (near x y z)
             e (- (double d2) (double d1))
             glass (s/smoothstep 0.06 0.09 e)
             r (ripple (* 6.0 (double x)) (* 6.0 (double y)) (* 6.0 (double z)))
             pane (scalec (ramp (* 0.5 (+ 1.0 (double v)))) (+ 0.85 (* 0.3 r)))]
         (texel (mixc [0.08 0.08 0.09] pane glass)
                (+ (* 0.3 (- 1.0 glass)) (* 0.05 glass r))
                (s/mix 0.5 0.8 glass)))))))

;; ---------------------------------------------------------------------------
;; Metals

(defn brushed-metal
  "Steel drawn under an abrasive: fine scratches all running one way.

  Noise squeezed along one axis and stretched along the other is a
  streak, and a few octaves of streaks are a brushed finish. What makes it
  metal rather than gray is the shine, which `shade` gives the color of
  the metal itself."
  ([] (brushed-metal {}))
  ([{:keys [seed color] :or {seed 0 color [0.62 0.64 0.67]}}]
   (let [streak (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves 5})]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             st (double (streak (* 0.4 x) (* 60.0 y) (* 60.0 z)))]
         (texel (scalec color (+ 0.92 (* 0.15 st))) (* 0.03 st) 0.9))))))

(defn hammered-metal
  "Copper beaten with a round hammer: overlapping shallow dimples.

  Each dimple is a cell, and its depth is the distance to the cell's
  point, squared -- a dish, deepest in the middle, meeting its neighbors
  in the ridges left between blows."
  ([] (hammered-metal {}))
  ([{:keys [seed density color] :or {seed 0 density 2.5 color [0.78 0.46 0.30]}}]
   (let [near (cells seed density)
         spacing (cellular/expected-spacing (double density))
         tarnish (fbm (+ seed 1) 4)]
     (fn [x y z]
       (let [{[d1] :distances v :value} (near x y z)
             d1 (/ (double d1) spacing)
             t (tarnish x y z)]
         (texel (scalec (mixc color [0.55 0.30 0.20] (* 0.6 t)) (+ 0.95 (* 0.08 (double v))))
                (* 2.5 (min 2.0 (* d1 d1)))
                0.85))))))

(defn rust
  "Steel going to rust: rough orange-brown patches eating into a bright
  surface.

  One fBm decides where the rust is and a finer one what it looks like,
  which is Worley's way of building a material in layers -- each layer a
  mask and a look, composed in order. The rust is raised, rough and dull;
  the steel between it smooth and bright."
  ([] (rust {}))
  ([{:keys [seed amount] :or {seed 0 amount 0.5}}]
   (let [where (fbm seed 6)
         look (fbm (+ seed 1) 6)
         pits (noise/gradient-basis {:seed (+ seed 2)})
         amount (double amount)
         ramp (color-map [[0.30 0.12 0.06] [0.30 0.12 0.06] [0.55 0.24 0.08]
                          [0.72 0.38 0.14] [0.45 0.20 0.10] [0.45 0.20 0.10]])]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             m (s/smoothstep (- 1.05 amount) (- 1.15 amount) (+ (where x y z) (* 0.15 (look (* 4.0 x) (* 4.0 y) (* 4.0 z)))))
             l (look (* 6.0 x) (* 6.0 y) (* 6.0 z))
             p (double (pits (* 30.0 x) (* 30.0 y) (* 30.0 z)))
             steel (scalec [0.60 0.62 0.64] (+ 0.92 (* 0.1 p)))]
         (texel (mixc steel (ramp l) m)
                (* m (+ 0.25 (* 0.3 l) (* 0.1 p)))
                (s/mix 0.85 0.05 m)))))))

(defn diamond-plate
  "Tread plate: raised lozenges on steel, alternate ones turned a quarter
  turn.

  Pure pattern, no noise in the shape at all: a grid, a checkerboard of
  orientations, and in each cell an elongated lozenge, its height falling
  off from a ridge down its middle."
  ([] (diamond-plate {}))
  ([{:keys [seed size] :or {seed 0 size 0.5}}]
   (let [size (double size)
         scratch (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves 4})]
     (fn [x y z]
       (let [x (double x) y (double y) z (double z)
             gx (/ x size) gy (/ y size)
             i (math/floor gx) j (math/floor gy)
             a (- gx i 0.5) b (- gy j 0.5)
             ;; Rotated 45 degrees one way or the other, by checkerboard.
             [a b] (if (even? (long (+ i j)))
                     [(* 0.7071 (+ a b)) (* 0.7071 (- b a))]
                     [(* 0.7071 (- a b)) (* 0.7071 (+ a b))])
             d (+ (abs (/ a 0.5)) (abs (/ b 0.12)))
             lug (- 1.0 (s/smoothstep 0.75 1.0 d))
             st (double (scratch (* 0.5 x) (* 20.0 y) z))]
         (texel (scalec [0.66 0.68 0.70] (+ 0.92 (* 0.1 st) (* 0.06 lug)))
                (* lug (- 1.0 (* 0.5 d)))
                0.9))))))

(def catalog
  "Every texture here, by name, with what it is built from. `scale` is a
  size that suits it in a view a unit or two across."
  [{:name "brick" :make brick :scale 3.0 :source "chapter 2"}
   {:name "marble" :make marble :scale 1.5 :source "chapters 2, 15"}
   {:name "wood" :make wood :scale 1.5 :source "chapter 2"}
   {:name "granite" :make granite :scale 1.0 :source "chapter 2"}
   {:name "clouds" :make clouds :scale 2.0 :source "chapter 15"}
   {:name "cirrus" :make cirrus :scale 2.0 :source "chapter 15"}
   {:name "gas giant" :make gas-giant :scale 3.0 :source "chapter 15"}
   {:name "water" :make water :scale 2.0 :source "chapter 15"}
   {:name "earth" :make earth :scale 3.0 :source "chapter 15"}
   {:name "random coloring" :make random-coloring :scale 2.0 :source "chapter 15"}
   {:name "impressionist" :make #(impressionist (earth %) %) :scale 3.0
    :source "chapter 15"}
   {:name "flagstone" :make flagstone :scale 3.0 :source "chapter 4"}
   {:name "cobblestone" :make cobblestone :scale 3.0 :source "chapter 4"}
   {:name "cracked mud" :make cracked-mud :scale 2.0 :source "chapter 4"}
   {:name "stained glass" :make stained-glass :scale 3.0 :source "chapter 4"}
   {:name "brushed metal" :make brushed-metal :scale 1.0 :source "chapters 2, 3"}
   {:name "hammered metal" :make hammered-metal :scale 3.0 :source "chapter 4"}
   {:name "rust" :make rust :scale 1.5 :source "chapter 3"}
   {:name "diamond plate" :make diamond-plate :scale 2.0 :source "chapter 2"}])

;; ---------------------------------------------------------------------------
;; Lighting

(def default-light
  "Up and to the left, toward the viewer, which is where a picture of a
  surface wants its light so that what is raised reads as raised."
  (let [[x y z] [-0.5 0.6 0.62]
        l (math/sqrt (+ (* x x) (* y y) (* z z)))]
    [(/ x l) (/ y l) (/ z l)]))

(defn shade
  "The lit color `[r g b]` of `texture` at point `p`, on a surface with
  normal `n` and tangents `tu`, `tv`, seen from straight ahead.

  The bumps come from the texture's own height, by forward differences
  `eps` along each tangent -- best about the size of a pixel, so finer
  bumps than can be shown are not lit into noise. The normal is tilted by
  the slope times `:bump`, and lit by a diffuse term and a specular one
  with `shine` as its strength, tinted by the color where the surface is
  glossy, as a metal's highlight is. Sunk places get less of the ambient
  light, which is cheap and most of what makes mortar look deep."
  [texture [px py pz] [nx ny nz] [ux uy uz] [vx vy vz] eps
   {:keys [bump light ambient] :or {bump 0.04 light default-light ambient 0.3}}]
  (let [px (double px) py (double py) pz (double pz) eps (double eps)
        [r g b h shine] (texture px py pz)
        h (double h)
        hu (double (nth (texture (+ px (* eps ux)) (+ py (* eps uy)) (+ pz (* eps uz))) 3))
        hv (double (nth (texture (+ px (* eps vx)) (+ py (* eps vy)) (+ pz (* eps vz))) 3))
        k (/ (double bump) eps)
        su (* k (- hu h)) sv (* k (- hv h))
        mx (- nx (* su ux) (* sv vx))
        my (- ny (* su uy) (* sv vy))
        mz (- nz (* su uz) (* sv vz))
        ml (math/sqrt (+ (* mx mx) (* my my) (* mz mz)))
        mx (/ mx ml) my (/ my ml) mz (/ mz ml)
        [lx ly lz] light
        diffuse (max 0.0 (+ (* mx lx) (* my ly) (* mz lz)))
        ;; Blinn's half vector between the light and a viewer on +z.
        hx lx hy ly hz (+ 1.0 (double lz))
        hl (math/sqrt (+ (* hx hx) (* hy hy) (* hz hz)))
        spec (* (double shine)
                (math/pow (max 0.0 (/ (+ (* mx hx) (* my hy) (* mz hz)) hl))
                          (s/mix 8.0 80.0 shine)))
        amb (* (double ambient) (+ 0.55 (* 0.45 (am/clamp h 0.0 1.0))))
        lit (+ amb (* (- 1.0 (* 0.5 (double shine))) diffuse))
        tint (fn [c] (+ (* (double c) lit) (* spec (s/mix 1.0 (double c) shine))))]
    [(am/clamp (tint r) 0.0 1.0) (am/clamp (tint g) 0.0 1.0) (am/clamp (tint b) 0.0 1.0)]))

(defn plane-pixel
  "The lit color of `texture` on the plane z = 0, at `(x, y)`, for a view
  where a pixel is `pixel` units across."
  [texture x y pixel opts]
  (shade texture [x y 0.0] [0.0 0.0 1.0] [1.0 0.0 0.0] [0.0 1.0 0.0] pixel opts))

(defn sphere-pixel
  "The lit color of `texture` on a sphere of radius `scale` at the
  origin, seen from +z, where `(sx, sy)` in the unit disk is where on
  the picture of it; nil off the sphere.

  The point is on the sphere, so the solid texture is what is carved
  there: no seam, no pinching at the poles."
  [texture sx sy scale pixel opts]
  (let [sx (double sx) sy (double sy) r2 (+ (* sx sx) (* sy sy))]
    (when (< r2 1.0)
      (let [sz (math/sqrt (- 1.0 r2))
            ;; East and north on the sphere; at the poles, any.
            ul (math/sqrt (+ (* sz sz) (* sx sx)))
            [ux uy uz] (if (< ul 1e-9) [1.0 0.0 0.0] [(/ sz ul) 0.0 (/ (- sx) ul)])
            vx (* sy uz) vy (- (* sz ux) (* sx uz)) vz (- (* sy ux))
            scale (double scale)]
        (shade texture [(* scale sx) (* scale sy) (* scale sz)] [sx sy sz]
               [ux uy uz] [vx vy vz] (* scale (double pixel)) opts)))))
