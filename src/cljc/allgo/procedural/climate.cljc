(ns allgo.procedural.climate
  "Wind, rain, rivers and biomes for a planet of plates.

  `allgo.procedural.tectonics` says where the ground is high. This says
  what that does to the weather: the ranges catch the rain on their
  windward side and leave deserts in their lee, the subtropics dry out
  under sinking air whatever the ground does, and the rain that falls
  runs down to the sea in rivers that join as they go.

  ## The steps, over the plates' cells

      temperature  warmest at the equator, falling toward the poles and
                   with height -- six and a half degrees a kilometer,
                   the standard atmosphere's lapse rate
      wind         the three cells of the general circulation: trade
                   winds from the east up to thirty degrees, drifting
                   toward the equator; westerlies from thirty to sixty,
                   drifting poleward; polar easterlies beyond
      moisture     the sea saturates the air above it, the more the
                   warmer it is, and each land cell takes its air from
                   the cells upwind of it. Some of what arrives falls:
                   more where the air rises -- up a mountainside, and in
                   the rising belts at the equator and around sixty
                   degrees -- and less in the sinking belt around thirty,
                   which is where the world's great deserts are. What does
                   not fall moves on, less what the ground breathes back
      rivers       every cell drains to its lowest neighbor, after the
                   pits are filled so that water always finds the sea --
                   a filled pit is a lake -- and the rain is summed down
                   those paths, so a river is wherever enough of it meets
      biomes       Whittaker's diagram: temperature on one axis,
                   precipitation on the other

  The moisture is the one step that is not a single pass. A cell's air
  comes from its upwind neighbors, which take theirs from further upwind,
  and on a sphere there is no order that visits every cell after all the
  cells upwind of it -- the westerlies go round and round. So the passes
  repeat until nothing changes, each carrying the air a few cells further
  inland.

  ## What it is not

  A climate rather than weather, and a caricature of one: no seasons, no
  ocean currents, no monsoon. Winds follow latitude alone, as they would
  on a planet with no continents. It is enough to put the deserts and
  rainforests where they belong relative to each other and to the
  mountains, which is what a map needs.

  The rivers are only as long as the terrain lets them be. The plates'
  ground rises steadily from the coast inland, with no basins in the
  interior to gather water, so it drains straight out to the nearest
  shore: on ten thousand cells the longest river is some half a dozen
  cells and most are two. Long rivers need a long, shallow basin, which
  is a thing erosion makes and plates alone do not."
  (:require [clojure.math :as math]))

(def defaults
  {:sea-level 0.0
   ;; how many kilometers one unit of the plates' elevation stands for
   :km-per-unit 6.0
   ;; how many times more readily air rains out per unit it is lifted up a
   ;; slope: rising air rains more, in proportion to how readily it rains
   ;; at that latitude at all, so a range in the trades is wet on its
   ;; windward side and a range under the subtropical high much less so
   :orographic 30.0
   ;; how much of the rain the ground gives back to the air
   :recycling 0.5
   :passes 120
   ;; the share of the land, by how much rain it carries, that is river
   :rivers 0.03})

;; ---------------------------------------------------------------------------
;; Vectors and tables

(defn- dot ^double [[ax ay az] [bx by bz]]
  (+ (* (double ax) (double bx)) (* (double ay) (double by)) (* (double az) (double bz))))

(defn- cross [[ax ay az] [bx by bz]]
  (let [ax (double ax) ay (double ay) az (double az)
        bx (double bx) by (double by) bz (double bz)]
    [(- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))]))

(defn- unit [v]
  (let [l (math/sqrt (dot v v))]
    (if (< l 1e-12) [0.0 0.0 0.0] (mapv #(/ (double %) l) v))))

(defn- interpolate
  "The piecewise-linear `table`, `[[x y] ...]` sorted by x, at `x`."
  ^double [table ^double x]
  (let [[x0 y0] (first table)]
    (if (<= x (double x0))
      (double y0)
      (loop [[[xa ya] [xb yb] & _ :as t] table]
        (cond
          (nil? xb) (double ya)
          (<= x (double xb)) (+ (double ya) (* (- (double yb) (double ya))
                                               (/ (- x (double xa)) (- (double xb) (double xa)))))
          :else (recur (rest t)))))))

(defn- quantile
  "The value a fraction `q` of the way up `xs`, sorted; 0 for none."
  ^double [xs ^double q]
  (let [v (vec (sort xs))]
    (if (empty? v)
      0.0
      (double (v (min (dec (count v)) (long (* q (count v)))))))))

(defn latitude
  "Latitude of a unit vector, radians, the pole along y as in
  `allgo.procedural.planet`."
  ^double [[_ y _]]
  (math/asin (max -1.0 (min 1.0 (double y)))))

;; ---------------------------------------------------------------------------
;; Temperature and wind

(defn temperature
  "Mean annual temperature, degrees Celsius, at latitude `lat` (radians)
  and `height` km above the sea: 28 at the equator, -30 at the poles,
  falling as the square of latitude, which puts the Earth's 13 degrees
  near forty-five -- and 6.5 degrees colder for every kilometer up."
  ^double [^double lat ^double height]
  (let [f (/ (abs lat) (* 0.5 math/PI))]
    (- 28.0 (* 58.0 f f) (* 6.5 (max 0.0 height)))))

(def ^:private zonal
  "The prevailing wind from west to east, by latitude in degrees:
  negative is an easterly. Trades, westerlies, polar easterlies, with the
  calms between them."
  [[0.0 -0.7] [15.0 -1.0] [30.0 0.0] [45.0 1.0] [60.0 0.0] [75.0 -0.6] [90.0 -0.3]])

(defn wind
  "The prevailing wind at unit vector `p`, a vector along the surface.

  The east-west part from the bands above; the north-south part is the
  drift that closes each cell of the circulation, toward the equator in
  the trades and the polar easterlies and toward the pole in the
  westerlies."
  [p]
  (let [lat (latitude p)
        east (unit (cross p [0.0 1.0 0.0]))
        north (cross east p)
        u (interpolate zonal (math/to-degrees (abs lat)))
        v (* -0.3 (math/sin (* 6.0 (abs lat))) (math/signum lat))]
    (mapv + (map #(* u (double %)) east) (map #(* v (double %)) north))))

(def ^:private rising
  "How readily air rains out at each latitude, degrees, before any
  mountain lifts it: the rising belts at the equator and near sixty,
  the sinking one near thirty that makes the subtropical deserts."
  [[0.0 0.11] [10.0 0.09] [22.0 0.025] [32.0 0.02] [45.0 0.05] [58.0 0.07] [75.0 0.03] [90.0 0.015]])

;; ---------------------------------------------------------------------------
;; Moisture

(defn- upwind-weights
  "For each cell, its neighbors weighted by how directly the wind blows
  from each of them toward it -- or, where the air is calm, all alike."
  [points neighbors winds]
  (mapv (fn [i]
          (let [p (points i)
                w (unit (winds i))
                ws (mapv (fn [j] (max 0.0 (dot (unit (mapv - p (points j))) w)))
                         (nth neighbors i))
                total (reduce + ws)]
            (if (< total 1e-9)
              (let [k (count ws)] (vec (repeat k (/ 1.0 k))))
              (mapv #(/ % total) ws))))
        (range (count points))))

(defn- moisture
  "Each cell's rainfall and the humidity of the air leaving it, as
  `[rain humidity]`, by repeated passes until the air has reached
  everywhere it is going to."
  [{:keys [points neighbors]} land? elevation temperatures weights
   {:keys [orographic recycling passes sea-level]}]
  (let [n (count points)
        humid (double-array n)
        rain (double-array n)
        orographic (double orographic)
        keep-rate (- 1.0 (double recycling))
        sea (double sea-level)
        ;; warm sea saturates the air more than cold
        sea-air (fn [t] (max 0.1 (min 1.0 (+ 0.35 (/ (double t) 40.0)))))
        base (mapv #(interpolate rising (math/to-degrees (abs (latitude %)))) points)]
    (dotimes [i n]
      (when-not (land? i) (aset humid i (double (sea-air (temperatures i))))))
    (dotimes [_ (long passes)]
      (dotimes [i n]
        (when (land? i)
          (let [nb (nth neighbors i)
                ws (nth weights i)
                [h-in e-up] (loop [k 0 h 0.0 e 0.0]
                              (if (= k (count nb))
                                [h e]
                                (let [j (nth nb k) w (double (nth ws k))]
                                  (recur (inc k)
                                         (+ h (* w (aget humid (long j))))
                                         (+ e (* w (max sea (double (elevation j)))))))))
                rise (max 0.0 (- (double (elevation i)) e-up))
                fraction (min 1.0 (* (double (base i)) (+ 1.0 (* orographic rise))))
                p (* h-in fraction)]
            (aset rain i p)
            (aset humid i (- h-in (* keep-rate p)))))))
    [(vec rain) (vec humid)]))

;; ---------------------------------------------------------------------------
;; Rivers

(defn- drainage
  "Each cell's downstream neighbor, or -1 for the sea, after the pits are
  filled: `[downstream filled]`.

  Barnes, Lehman and Mulla's priority flood. The sea is the outlet, and
  cells are taken lowest first from its edge inward; each cell reached
  drains to the one that reached it, and stands at least a hair above it.
  A cell lower than that has been filled -- it is in a pit, and the pit
  is a lake whose surface is the filled height."
  [neighbors elevation land?]
  (let [n (count neighbors)
        down (int-array n -1)
        filled (double-array n)
        seen (int-array n)
        start (for [i (range n) :when (not (land? i))] i)]
    (doseq [i start]
      (aset seen i 1)
      (aset filled i (double (elevation i))))
    (loop [queue (into (sorted-set) (for [i start
                                          :when (some land? (nth neighbors i))]
                                      [(double (elevation i)) i]))]
      (when-let [[h c :as top] (first queue)]
        (let [fresh (remove #(pos? (aget seen (long %))) (nth neighbors c))]
          (doseq [f fresh]
            (aset seen (long f) 1)
            (aset down (long f) (int c))
            (aset filled (long f) (max (double (elevation f)) (+ (double h) 1e-6))))
          (recur (into (disj queue top) (map (fn [f] [(aget filled (long f)) f]) fresh))))))
    [(vec down) (vec filled)]))

(defn- accumulate
  "How much rain each cell carries: its own and everything upstream of
  it, summed from the highest cell down."
  [down filled rain land?]
  (let [n (count down)
        flow (double-array n)
        order (sort-by #(- (double (filled %))) (filter land? (range n)))]
    (doseq [i order] (aset flow (long i) (+ (aget flow (long i)) (double (rain i)))))
    (doseq [i order]
      (let [d (long (down i))]
        (when (and (>= d 0) (land? d))
          (aset flow d (+ (aget flow d) (aget flow (long i)))))))
    (vec flow)))

;; ---------------------------------------------------------------------------
;; Biomes

(defn biome
  "Whittaker's biome for mean temperature `t` (Celsius) and precipitation
  `p`, as a fraction of the wettest the model produces."
  [^double t ^double p]
  (cond
    (< t -10.0) :ice
    (< t 0.0) (if (< p 0.15) :cold-desert :tundra)
    (< t 8.0) (cond (< p 0.12) :cold-desert (< p 0.3) :tundra :else :taiga)
    (< t 20.0) (cond (< p 0.1) :desert (< p 0.25) :grassland (< p 0.55) :temperate-forest
                     :else :temperate-rainforest)
    :else (cond (< p 0.1) :desert (< p 0.3) :savanna (< p 0.55) :seasonal-forest
                :else :tropical-rainforest)))

(def biome-colors
  "Each biome's color, `[r g b]` in [0, 1] -- what it looks like from
  orbit rather than on a chart, so a desert is sand and a taiga nearly
  black."
  {:ice [0.92 0.94 0.97]
   :tundra [0.55 0.55 0.45]
   :cold-desert [0.62 0.58 0.50]
   :taiga [0.13 0.24 0.16]
   :desert [0.80 0.70 0.48]
   :grassland [0.52 0.56 0.30]
   :temperate-forest [0.20 0.36 0.15]
   :temperate-rainforest [0.10 0.30 0.14]
   :savanna [0.62 0.58 0.32]
   :seasonal-forest [0.25 0.40 0.14]
   :tropical-rainforest [0.06 0.28 0.08]
   :lake [0.10 0.30 0.50]})

;; ---------------------------------------------------------------------------
;; Assembling

(defn climate
  "The climate of a `allgo.procedural.tectonics` world.

  `:sea-level` is in the plates' own units, zero being their sea;
  the rest of `defaults` tune it. Returns, per cell:

      :temperature    degrees Celsius
      :rain           how much falls, as a fraction of the wettest cell
      :downstream     the cell it drains to, -1 for the sea
      :flow           the rain it carries, as a fraction of the most any
                      cell does
      :river?         whether it carries enough to be a river
      :lake?          whether it is in a filled pit
      :biome          a key of `biome-colors`, nil at sea

  and `:wind`, the prevailing wind at each."
  ([world] (climate world {}))
  ([{:keys [points neighbors elevation] :as world} opts]
   (let [{:keys [sea-level km-per-unit rivers] :as opts} (merge defaults opts)
         sea (double sea-level)
         n (count points)
         land? (fn [i] (>= (double (elevation i)) sea))
         temperatures (mapv (fn [i] (temperature (latitude (points i))
                                                 (* (double km-per-unit) (- (double (elevation i)) sea))))
                            (range n))
         winds (mapv wind points)
         weights (upwind-weights points neighbors winds)
         [rain _] (moisture world land? elevation temperatures weights opts)
         land (filter land? (range n))
         ;; Rain is measured against the wettest twentieth of the land
         ;; rather than the wettest cell: rain is heavy-tailed, and one
         ;; windward slope taking the scale would leave every other cell
         ;; reading as desert.
         wet (max 1e-9 (quantile (map rain land) 0.95))
         rain (mapv #(min 1.0 (/ (double %) wet)) rain)
         [down filled] (drainage neighbors elevation land?)
         flow (accumulate down filled rain land?)
         most (max 1e-9 (reduce max 0.0 flow))
         flow (mapv #(/ (double %) most) flow)
         ;; The rivers are the land carrying the most water: a share of
         ;; the cells, so a wet world and a dry one both get some.
         river-flow (quantile (map flow land) (- 1.0 (double rivers)))
         lake? (mapv (fn [i] (and (land? i) (> (double (filled i)) (+ (double (elevation i)) 1e-6))))
                     (range n))]
     {:temperature temperatures
      :rain rain
      :wind winds
      :downstream down
      :flow flow
      :river? (mapv (fn [i] (and (land? i) (not (lake? i)) (>= (double (flow i)) (double river-flow))))
                    (range n))
      :lake? lake?
      :biome (mapv (fn [i]
                     (cond
                       (not (land? i)) nil
                       (lake? i) :lake
                       :else (biome (temperatures i) (rain i))))
                   (range n))})))
