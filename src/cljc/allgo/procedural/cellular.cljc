(ns allgo.procedural.cellular
  "Worley's cellular basis: texture built from the distances to scattered
  feature points.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, chapter 4
  (Worley, \"Cellular Texturing\"), and one of the bases MojoWorld offers
  in chapter 20. Every other basis in this library interpolates something
  stored on a lattice. This one does not interpolate at all: scatter
  points through space, and at any given point report how far away the
  nearest ones are.

  That gives a completely different vocabulary of shapes. `f1` alone is a
  field of rounded bumps -- flagstones, pebbles, cell walls seen from
  inside. `f2-f1` is zero exactly on the boundaries between points and
  positive away from them, which draws the Voronoi diagram of the feature
  points as a network of sharp creases: cracked mud, crazed glaze, dragon
  scales, the veins in a leaf. Neither is achievable by summing noise, and
  both are recognizable at a glance, which is the argument for the basis.

  ## Scattering the points

  A Poisson process, and it has to be: any regularity in how many points a
  cell gets, or where in the cell they sit, shows up in the texture as a
  grid. Space is cut into unit cells only so that the search can be
  finite; a cell's points are derived from its integer coordinates by
  `allgo.random/cell-seed`, so nothing is stored and the field
  is unbounded.

  The search covers the 3x3x3 cells around the sample. That is an
  approximation -- a point in the next ring out could in principle be
  nearer than one found, if the near cells happened to come up empty --
  and it is the approximation Worley's paper makes too, with the bound
  that makes it exact available for the asking. At the default density of
  around three points per cell the chance of it mattering is not
  measurable in a picture."
  (:require [allgo.procedural.noise :as noise]
            [allgo.random :as random]
            [clojure.math :as math]))

(def metrics
  "Distance functions, keyed by name. The metric changes the *shape* of a
  cell as much as the point scatter does.

  * `:euclidean` -- round cells, the default and the one that looks
    physical.
  * `:manhattan` -- cells with corners on the axes; reads as something
    machined or crystalline.
  * `:chebyshev` -- square cells; the most obviously artificial, and the
    right answer for anything meant to look like tiling or circuitry."
  {:euclidean (fn ^double [^double dx ^double dy ^double dz]
                (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz))))
   :manhattan (fn ^double [^double dx ^double dy ^double dz]
                (+ (abs dx) (abs dy) (abs dz)))
   :chebyshev (fn ^double [^double dx ^double dy ^double dz]
                (max (abs dx) (max (abs dy) (abs dz))))})

(defn feature-points
  "The feature points of one lattice cell, as `[[x y z v] ...]`.

  `v` is a random value in [-1, 1] carried by the point -- Worley's
  per-feature attribute, which is what lets a cell be given a flat color
  or a height of its own rather than only a distance."
  [{:keys [seed density jitter] :or {seed 0 density 3.0 jitter 1.0}} i j k]
  (let [seed (long seed) i (long i) j (long j) k (long k)
        density (double density) jitter (double jitter)
        s0 (random/cell-seed seed i j k)
        m (noise/poisson-count density (random/unit s0))
        center (* 0.5 (- 1.0 jitter))]
    (loop [n 0 s (random/advance s0) acc []]
      (if (>= n m)
        acc
        (let [px (+ i center (* jitter (random/unit s)))
              s (random/advance s)
              py (+ j center (* jitter (random/unit s)))
              s (random/advance s)
              pz (+ k center (* jitter (random/unit s)))
              s (random/advance s)
              v (random/signed s)]
          (recur (inc n) (random/advance s) (conj acc [px py pz v])))))))

(defn- insert!
  "Slides `d` into `best`, an ascending array of the smallest distances so
  far. Insertion sort, which is the right sort for four elements."
  [^doubles best ^double d]
  (let [n (alength best)]
    (when (< d (aget best (dec n)))
      (loop [i (dec n)]
        (if (and (pos? i) (> (aget best (dec i)) d))
          (do (aset best i (aget best (dec i)))
              (recur (dec i)))
          (aset best i d))))))

(defn- searcher
  "The search itself, shared by everything below.

  `(search best bp x y z)` fills `best` with the smallest distances found,
  ascending, and `bp` with the nearest feature point and the value it
  carries -- four doubles, `[x y z value]`. Both arrays belong to the
  caller.

  Written to fill arrays rather than to return a value because of what it
  costs not to: this runs once per sample per octave, and a version that
  allocated a vector per cell examined -- twenty-seven of them per call --
  measured an order of magnitude slower than this one."
  [{:keys [seed density jitter metric]
    :or {seed 0 density 3.0 jitter 1.0 metric :euclidean}}]
  (let [seed (long seed) density (double density) jitter (double jitter)
        dist (get metrics metric (:euclidean metrics))
        center (* 0.5 (- 1.0 jitter))]
    (fn [^doubles best ^doubles bp x y z]
      (let [x (double x) y (double y) z (double z)
            ci (long (math/floor x)) cj (long (math/floor y)) ck (long (math/floor z))]
        (dotimes [i (alength best)] (aset best i ##Inf))
        (loop [c 0 bd ##Inf]
          (if (= c 27)
            bd
            (let [i (+ ci (- (rem c 3) 1))
                  j (+ cj (- (rem (quot c 3) 3) 1))
                  k (+ ck (- (quot c 9) 1))
                  s0 (random/cell-seed seed i j k)
                  m (noise/poisson-count density (random/unit s0))]
              (recur
               (inc c)
               (double
                (loop [p 0 s (random/advance s0) bd bd]
                  (if (>= p m)
                    bd
                    (let [px (+ i center (* jitter (random/unit s)))
                          s (random/advance s)
                          py (+ j center (* jitter (random/unit s)))
                          s (random/advance s)
                          pz (+ k center (* jitter (random/unit s)))
                          s (random/advance s)
                          v (random/signed s)
                          d (double (dist (- x px) (- y py) (- z pz)))]
                      (insert! best d)
                      (if (< d bd)
                        (do (aset bp 0 px) (aset bp 1 py) (aset bp 2 pz) (aset bp 3 v)
                            (recur (inc p) (random/advance s) d))
                        (recur (inc p) (random/advance s) bd))))))))))))))

(defn nearest
  "`(fn [x y z] -> {:distances [...] :value v :point [...]})`.

  `:distances` are the `:n` smallest, ascending, in lattice units.
  `:value` and `:point` belong to the single nearest feature point. Use it
  when you want more than one of those answers and would rather not pay
  for the search twice; the scalar bases below go through the same search
  without building the map."
  ([] (nearest {}))
  ([{:keys [n] :or {n 2} :as opts}]
   (let [search (searcher opts)
         n (long n)]
     (fn [^double x ^double y ^double z]
       (let [best (double-array n)
             bp (double-array 4)]
         (search best bp x y z)
         {:distances (vec best)
          :value (aget bp 3)
          :point [(aget bp 0) (aget bp 1) (aget bp 2)]})))))

(defn expected-spacing
  "The mean distance from a point in space to the nearest feature point,
  for a Poisson process of `density` points per unit cell.

  `Gamma(4/3) / (4*pi*density/3)^(1/3)`, which is the standard result and
  worth having in closed form: it is what the bases below divide by, so
  that changing the density changes how big the cells are without also
  changing how light or dark the texture is."
  ^double [^double density]
  (/ 0.8929795115692493                  ; Gamma(4/3)
     (math/cbrt (* (/ (* 4.0 math/PI) 3.0) density))))

(defn distance-basis
  "`(fn [x y z] -> [d1 d2 ...])`, the `:n` nearest distances in units of
  `expected-spacing`, so that the numbers mean the same thing at any
  density."
  ([] (distance-basis {}))
  ([{:keys [density n] :or {density 3.0 n 2} :as opts}]
   (let [search (searcher opts)
         n (long n)
         k (/ 1.0 (expected-spacing (double density)))]
     (fn [^double x ^double y ^double z]
       (let [best (double-array n)
             bp (double-array 4)]
         (search best bp x y z)
         (mapv #(* k (double %)) best))))))

(defn- scalar-basis
  "Wraps a reduction over the nearest distances as an ordinary basis,
  centered on zero and scaled to land near [-1, 1].

  The centring and scaling constants at the call sites are measured rather
  than derived -- the distributions of F1, F2 and F2-F1 have closed forms
  only for the Euclidean metric, and the point of them is that a cellular
  basis and a gradient basis can be swapped for one another inside a
  fractal without every other parameter having to move too. Under
  `:manhattan` and `:chebyshev`, distances run larger and the result sits
  off center."
  [opts n f k]
  (let [search (searcher opts)
        n (long n)
        e (expected-spacing (double (:density opts 3.0)))
        k (double k)]
    (fn ^double [^double x ^double y ^double z]
      (let [best (double-array n)
            bp (double-array 4)]
        (search best bp x y z)
        (dotimes [i n] (aset best i (/ (aget best i) e)))
        (* k (double (f best)))))))

(defn f1
  "Distance to the nearest feature point: a field of rounded bumps.

  Zero at the mean spacing, negative in a feature point's immediate
  neighborhood and positive out toward the cell boundaries. Inverted, it
  is the classic pebbles-and-bubbles texture."
  ([] (f1 {}))
  ([opts] (scalar-basis opts 1 (fn [^doubles d] (- (aget d 0) 1.0)) 1.1)))

(defn f2
  "Distance to the second-nearest feature point."
  ([] (f2 {}))
  ([opts] (scalar-basis opts 2 (fn [^doubles d] (- (aget d 1) 1.32)) 1.12)))

(defn f2-f1
  "`F2 - F1`: zero exactly on the Voronoi boundaries, positive elsewhere.

  The cracked-mud basis. Every crease is a line equidistant from two
  feature points, so this draws the Voronoi diagram of the scatter without
  ever computing the diagram."
  ([] (f2-f1 {}))
  ([opts] (scalar-basis opts 2 (fn [^doubles d] (- (- (aget d 1) (aget d 0)) 0.33)) 1.0)))

(defn combination
  "Worley's linear combination `c1*F1 + c2*F2 + ...` of the nearest
  distances.

  The whole family in one function: `[-1 1]` is `f2-f1`, `[1]` is `f1`,
  `[-1 0 1]` picks out a different set of creases again. Worley's chapter
  is largely a catalog of what the coefficients do, and this is the knob
  that catalog turns."
  ([coeffs] (combination {} coeffs))
  ([opts coeffs]
   (let [coeffs (mapv double coeffs)
         n (count coeffs)]
     (scalar-basis opts n
                   (fn [^doubles ds]
                     (loop [i 0 acc 0.0]
                       (if (= i n)
                         acc
                         (recur (inc i) (+ acc (* (double (nth coeffs i)) (aget ds i)))))))
                   1.0))))

(defn id-basis
  "The random value carried by the nearest feature point: flat cells.

  Constant within a cell and discontinuous at its boundaries, so this is
  the basis for anything made of discrete pieces that happen to tile
  space -- a mosaic, a field of crystal grains, islands each with their
  own elevation."
  ([] (id-basis {}))
  ([opts]
   (let [search (searcher opts)]
     (fn ^double [^double x ^double y ^double z]
       (let [best (double-array 1)
             bp (double-array 4)]
         (search best bp x y z)
         (aget bp 3))))))
