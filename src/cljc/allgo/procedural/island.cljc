(ns allgo.procedural.island
  "An island as a polygon map: Voronoi cells with a coastline, elevation,
  rivers, moisture and biomes, after Amit Patel's *Polygonal Map
  Generation for Games*.

  ## Why not just sample a function

  `allgo.procedural.planet` samples a continuous fractal and asks it for
  a height at a point, which is the right shape for terrain and the wrong
  shape for geography. A noise field has no notion of *this* island,
  *that* lake, the river between them, or how far inland you are. Those
  are facts about a region and its neighbours, and to have them at all
  you need the map to be made of countable things with edges between
  them.

  So this builds a graph and then walks it, and every interesting
  property falls out of a traversal rather than a formula:

      ocean      flood fill inwards from the edge of the map, which is
                 what makes a lake a lake -- not its shape, but the fact
                 that the fill never reached it
      elevation  breadth-first from the coast, so height is distance from
                 the sea. Mountains end up in the middle because that is
                 what being far from the sea means
      watershed  follow the same pointer to the end and you have the
                 stretch of coast this corner drains to; corners that
                 agree are one basin, and the ridges between basins fall
                 out of the disagreements
      rivers     follow the downhill pointer from a corner to the water
      moisture   breadth-first from rivers and lakes, decaying with each
                 step, so rain shadows and dry interiors happen on their
                 own

  Only the last step is a formula, and it is a lookup table: biome is a
  function of elevation and moisture, which is Whittaker's diagram.

  ## Corners carry the water, centers carry the land

  Elevation, downhill and rivers live on the *corners* of cells, and
  biomes live on the cells. That is not an implementation detail. Water
  runs along the boundaries between regions rather than through the
  middle of them -- a river is a border, which is why so many real ones
  are -- and a river that ran centre-to-centre would cut its own cells in
  half. Cells take their elevation and moisture as the average over their
  corners, which is the only place the two representations meet.

  ## What it is not

  No plate tectonics, no climate simulation, no erosion. Elevation here
  is distance from the coast with a curve applied, which produces a
  plausible island and not a physical one -- `allgo.procedural.erosion`
  is the namespace that simulates a cause. The two are complementary and
  do not currently meet: erosion carves a grid and this builds a graph."
  (:require [allgo.geometry.dual-mesh :as dm]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [clojure.math :as math]))

(def defaults
  {:points 900
   :bounds [[0.0 0.0] [1000.0 1000.0]]
   ;; Two is the usual dose: enough to even the cell sizes out, few
   ;; enough that the map still looks organic rather than like a
   ;; honeycomb. Each round costs another triangulation.
   :relax 2
   :shape :noise
   :seed 1
   :rng rand
   ;; A cell is water when this fraction of its corners are. Below a half
   ;; because a cell with a little water in it is a shore, and shores
   ;; should read as wet rather than as land that happens to be damp.
   :lake-threshold 0.3
   ;; Attempts, as a fraction of the corner count. Most are rejected for
   ;; starting too low, too high, or in the sea.
   :river-attempts 0.6
   ;; Draw the cell boundaries as wandering paths rather than straight
   ;; segments. On by default: a polygon map that looks polygonal is the
   ;; thing everyone wants to fix first, and it changes nothing but where
   ;; the ink goes. See `allgo.geometry.dual-mesh/noisy-edges`.
   :noisy? true})

;; ---------------------------------------------------------------------------
;; Island shapes
;;
;; A shape is `(fn [[nx ny]] -> land?)` over normalised coordinates, where
;; the map is [-1, 1] on both axes. Everything downstream asks only this
;; one question, so the difference between an archipelago and a single
;; round island is entirely here.

(defn- normalise [[[x0 y0] [x1 y1]]]
  (fn [[x y]]
    [(- (* 2.0 (/ (- (double x) x0) (- x1 x0))) 1.0)
     (- (* 2.0 (/ (- (double y) y0) (- y1 y0))) 1.0)]))

(defn radial-shape
  "Overlapping sine waves in polar coordinates: one blobby island with
  bays and headlands, always centred, always a single landmass.

  Predictable in a way the noise shape is not, which is what it is for --
  a map that must have exactly one island with a coast all the way round."
  [{:keys [seed bumps]}]
  (let [rng-seed (double (noise/cell-seed (long seed) 7 7 7))
        phase (* 2.0 math/PI (noise/unit rng-seed))
        s2 (noise/advance rng-seed)
        phase2 (* 2.0 math/PI (noise/unit s2))
        ;; The lobe count comes from the seed unless it is asked for.
        ;; Fixing it and varying only the phase gives every seed the same
        ;; starfish rotated, which is worse than no variation at all --
        ;; it looks like a bug rather than like a family of islands.
        bumps (long (or bumps (+ 3 (long (* 4.0 (noise/unit (noise/advance s2)))))))]
    (fn [[nx ny]]
      (let [r (math/sqrt (+ (* nx nx) (* ny ny)))
            a (math/atan2 ny nx)
            ;; Two waves of different frequency so the coast does not
            ;; repeat with the period of a single sine.
            edge (+ 0.58
                    (* 0.18 (math/sin (+ phase (* bumps a))))
                    (* 0.10 (math/sin (+ phase2 (* (+ bumps 3) a)))))]
        (< r edge)))))

(defn noise-shape
  "Fractal noise pulled down towards the edges of the map.

  The falloff is what turns a noise field into an island: without it the
  land runs off every side, and with it the same field becomes a coast.
  Raise `:falloff` for a smaller island, `:threshold` for a wetter world,
  and expect archipelagos rather than one landmass -- which is the reason
  to choose this over `radial-shape`.

  Few octaves on purpose. The shape is asked one question per corner, so
  detail finer than a cell cannot be drawn -- it can only be sampled, and
  what it produces is speckle: a coast that frays into single-cell islets
  and an interior pocked with one-cell lakes. At six octaves and this
  frequency the finest band is some seventy cycles across a map about
  thirty cells wide, which is all noise and no shape. Three octaves puts
  the smallest feature at a few cells, where the polygons can represent
  it."
  [{:keys [seed frequency octaves falloff power threshold]
    :or {frequency 1.2 octaves 3 falloff 0.6 power 2.0 threshold 0.18}}]
  ;; Simplex rather than Perlin, for the reason the basis exists: Perlin
  ;; fades along each axis and leaves a faint squareness, which on a
  ;; coastline reads as shores that prefer the compass points.
  (let [f (fractal/fbm (noise/simplex-basis {:seed seed}) {:octaves octaves})]
    (fn [[nx ny]]
      (let [d (math/sqrt (+ (* nx nx) (* ny ny)))
            e (* 0.5 (+ 1.0 (f (* frequency nx) (* frequency ny) 0.5)))]
        (> (- e (* (double falloff) (math/pow d (double power)))) (double threshold))))))

(def shapes
  "The shapes by name, so a demo can offer them without knowing the vars."
  {:radial radial-shape
   :noise noise-shape})

;; ---------------------------------------------------------------------------
;; Water, ocean and coast

(defn- assign-water
  "Which corners are under water, and then which cells are.

  The shape decides the corners; a cell follows its corners. Border
  corners are water whatever the shape says, because the map has to end
  in sea or the flood fill below has nowhere to start."
  [{:keys [centers corners] :as island} shape bounds lake-threshold]
  (let [to-unit (normalise bounds)
        corners (mapv (fn [c]
                        (assoc c :water? (or (:border? c)
                                             (not (shape (to-unit (:point c)))))))
                      corners)
        centers (mapv (fn [c]
                        (let [cs (:corners c)
                              wet (count (filter #(:water? (corners %)) cs))]
                          (assoc c :water? (>= (/ (double wet) (max 1 (count cs)))
                                               (double lake-threshold)))))
                      centers)]
    (assoc island :centers centers :corners corners)))

(defn- assign-ocean
  "Flood the sea inwards from the edge, then read off what that made.

  A lake is water the flood did not reach, and nothing about its shape
  says so -- only its connectivity. This is the pass that a heightmap
  cannot do, and the reason the whole map is a graph."
  [{:keys [centers corners] :as island}]
  (let [ocean (loop [queue (into #?(:clj clojure.lang.PersistentQueue/EMPTY :cljs #queue [])
                                 (comp (filter #(and (:border? %) (:water? %))) (map :id))
                                 centers)
                     seen (into #{} (comp (filter #(and (:border? %) (:water? %))) (map :id)) centers)]
                (if-let [id (peek queue)]
                  (let [nbrs (remove #(or (contains? seen %) (not (:water? (centers %))))
                                     (:neighbors (centers id)))]
                    (recur (into (pop queue) nbrs) (into seen nbrs)))
                  seen))
        centers (mapv (fn [c]
                        (let [o (contains? ocean (:id c))]
                          (assoc c :ocean? o :lake? (and (:water? c) (not o)))))
                      centers)
        centers (mapv (fn [c]
                        (assoc c :coast? (and (not (:water? c))
                                              (some #(:ocean? (centers %)) (:neighbors c))
                                              true)))
                      centers)
        corners (mapv (fn [v]
                        (let [ts (:touches v)
                              n (count ts)
                              n-ocean (count (filter #(:ocean? (centers %)) ts))
                              n-land (count (remove #(:water? (centers %)) ts))]
                          (assoc v
                                 :ocean? (and (pos? n) (= n-ocean n))
                                 :coast? (and (pos? n-ocean) (pos? n-land)))))
                      corners)]
    (assoc island :centers centers :corners corners)))

;; ---------------------------------------------------------------------------
;; Elevation

(defn- assign-corner-elevation
  "Height as distance from the sea.

  Breadth-first from the edge of the map, where elevation is zero, and
  every step inland costs something: a lot when both ends are dry, almost
  nothing when crossing water. So lakes come out flat, and the land rises
  the further it is from a coast -- which puts the mountains in the
  middle without anything having to decide where the middle is."
  [{:keys [corners] :as island}]
  (let [start (into [] (comp (filter :border?) (map :id)) corners)
        init (reduce (fn [v id] (assoc v id 0.0))
                     (vec (repeat (count corners) ##Inf))
                     start)
        elev (loop [queue (into #?(:clj clojure.lang.PersistentQueue/EMPTY :cljs #queue []) start)
                    elev init]
               (if-let [id (peek queue)]
                 (let [q (corners id)
                       here (elev id)
                       [elev pushed]
                       (reduce (fn [[elev pushed] sid]
                                 (let [s (corners sid)
                                       step (if (and (not (:water? q)) (not (:water? s)))
                                              1.0
                                              0.01)
                                       next-e (+ here step)]
                                   (if (< next-e (elev sid))
                                     [(assoc elev sid next-e) (conj pushed sid)]
                                     [elev pushed])))
                               [elev []]
                               (:adjacent q))]
                   (recur (into (pop queue) pushed) elev))
                 elev))]
    (assoc island :corners (mapv #(assoc % :elevation (elev (:id %))) corners))))

(defn- redistribute-elevation
  "Reshapes the land's height distribution so that low ground is common.

  The breadth-first pass gives elevations that are essentially a count of
  steps, and those come out roughly uniform -- as much high ground as
  low, which no real landmass has. Sorting the land corners and laying
  them on the curve `x = sqrt(S) - sqrt(S(1 - y))` gives the opposite:
  lots of coastal plain, a little high ground, and the peaks rare. It is
  the histogram that is being set, so the *shape* of the terrain -- which
  corner is higher than which -- comes through untouched.

  Lake corners are in the redistribution, not excluded from it. They are
  not at sea level -- a tarn sits at the height of the shelf it is on --
  and dropping them to zero with the ocean puts every lake at the coast
  and drags the cells around it down to meet them. Only ocean and coast
  corners are truly at sea level, and those are set afterwards."
  [{:keys [corners] :as island}]
  (let [land (into [] (comp (remove :ocean?) (remove :coast?) (map :id)) corners)
        ordered (vec (sort-by #(:elevation (corners %)) land))
        n (count ordered)
        scale 1.1
        elev (reduce (fn [m [i id]]
                       (let [y (if (> n 1) (/ (double i) (dec n)) 0.0)
                             x (- (math/sqrt scale) (math/sqrt (* scale (- 1.0 y))))]
                         (assoc m id (min 1.0 x))))
                     {}
                     (map-indexed vector ordered))]
    (assoc island :corners
           (mapv (fn [v] (assoc v :elevation (get elev (:id v) 0.0))) corners))))

(defn- assign-center-elevation
  "A cell is as high as its corners average."
  [{:keys [centers corners] :as island}]
  (assoc island :centers
         (mapv (fn [c]
                 (let [cs (:corners c)]
                   (assoc c :elevation (if (seq cs)
                                         (/ (reduce + (map #(:elevation (corners %)) cs))
                                            (count cs))
                                         0.0))))
               centers)))

(defn- assign-downslope
  "The neighbour each corner drains to, which may be itself.

  One pointer per corner and every river in the map follows from it. A
  corner that is its own downslope is a pit, and the river that reaches
  one stops there."
  [{:keys [corners] :as island}]
  (assoc island :corners
         (mapv (fn [v]
                 (assoc v :downslope
                        (reduce (fn [best a]
                                  (if (< (:elevation (corners a)) (:elevation (corners best)))
                                    a best))
                                (:id v)
                                (:adjacent v))))
               corners)))

;; ---------------------------------------------------------------------------
;; Watersheds

(defn- assign-watersheds
  "Which stretch of coast each corner's water eventually reaches.

  Follow the downhill pointer far enough and you arrive at the sea. The
  corner where you arrive names the basin, so two corners share a
  watershed exactly when the rain landing on them ends up in the same
  place -- which is a fact about the whole path between them and not
  about how close together they are. Basins meet along ridges, and the
  ridge is wherever two neighbouring corners give different answers.

  mapgen2 relaxes this by re-reading every corner's downhill neighbour a
  hundred times over and stopping when nothing moves. Walking each chain
  once and remembering what the walk found gets the same answer without
  the iteration count -- and it has somewhere to put the case that
  iteration quietly rounds off, which is a pit with no outlet. A corner
  that drains into one belongs to no basin, and says so with `nil`."
  [{:keys [corners] :as island}]
  (let [terminal
        (fn [cache start]
          ;; Walk down, remembering the path, then credit the whole path
          ;; with whatever the walk ended in. Each corner is resolved
          ;; once however many chains run through it.
          (loop [v start path [] seen #{}]
            (let [c (corners v)]
              (cond
                (contains? @cache v)
                (let [t (@cache v)]
                  (swap! cache into (zipmap path (repeat t)))
                  t)

                (or (:coast? c) (:ocean? c))
                (do (swap! cache into (zipmap (conj path v) (repeat v))) v)

                ;; A pit, or a loop the pointers fell into: no outlet, so
                ;; no basin.
                (or (= (:downslope c) v) (contains? seen v))
                (do (swap! cache into (zipmap (conj path v) (repeat nil))) nil)

                :else
                (recur (:downslope c) (conj path v) (conj seen v))))))
        cache (atom {})
        sheds (mapv #(terminal cache (:id %)) corners)
        sizes (frequencies (remove nil? sheds))]
    (assoc island :corners
           (mapv (fn [v]
                   (let [w (sheds (:id v))]
                     (assoc v :watershed w :watershed-size (get sizes w 0))))
                 corners))))

;; ---------------------------------------------------------------------------
;; Rivers

(defn- edge-between [{:keys [corners edges]} a b]
  (some (fn [eid] (let [[v0 v1] (:corners (edges eid))]
                    (when (or (= v0 b) (= v1 b)) eid)))
        (:protrudes (corners a))))

(defn- assign-rivers
  "Drop springs on the high ground and let them run down the pointers.

  A river is not routed or planned; it is the downhill pointer followed
  until it reaches the sea. Where two of them land on the same corner the
  counts add, which is why the trunk is wider than the tributaries
  without anything having to model discharge."
  [{:keys [corners edges] :as island} rng attempts]
  (let [n (long (* (double attempts) (count corners)))
        [corner-flow edge-flow]
        (loop [i 0
               cf (vec (repeat (count corners) 0))
               ef (vec (repeat (count edges) 0))]
          (if (>= i n)
            [cf ef]
            (let [start ((:corners island) (long (* (rng) (count corners))))]
              (if (or (:ocean? start)
                      (< (:elevation start) 0.3)
                      (> (:elevation start) 0.9))
                (recur (inc i) cf ef)
                ;; Walk down until the coast, a pit, or a loop.
                (let [[cf ef]
                      (loop [q (:id start) cf cf ef ef steps 0]
                        (let [v (corners q)
                              d (:downslope v)]
                          (if (or (:coast? v) (:ocean? v) (= d q) (> steps 500))
                            [cf ef]
                            (let [eid (edge-between island q d)]
                              (recur d
                                     (-> cf (update q inc) (update d inc))
                                     (if eid (update ef eid inc) ef)
                                     (inc steps))))))]
                  (recur (inc i) cf ef))))))]
    (assoc island
           :corners (mapv #(assoc % :river (corner-flow (:id %))) corners)
           :edges (mapv #(assoc % :river (edge-flow (:id %))) edges))))

;; ---------------------------------------------------------------------------
;; Moisture

(defn- assign-moisture
  "Wet where the fresh water is, and drier with every step away from it.

  Breadth-first from rivers and lakes, losing a tenth each hop. Nothing
  models wind or rain, and yet interiors come out dry and river valleys
  green, because distance from fresh water is most of what moisture is."
  [{:keys [centers corners] :as island}]
  (let [seed-moisture (fn [v]
                        (cond
                          (pos? (long (:river v 0))) (min 3.0 (* 0.2 (:river v)))
                          (:water? v) 1.0
                          :else 0.0))
        init (mapv seed-moisture corners)
        start (into [] (comp (filter #(pos? (init (:id %)))) (map :id)) corners)
        moist (loop [queue (into #?(:clj clojure.lang.PersistentQueue/EMPTY :cljs #queue []) start)
                     m init]
                (if-let [id (peek queue)]
                  (let [here (* 0.9 (m id))
                        [m pushed] (reduce (fn [[m pushed] a]
                                             (if (> here (m a))
                                               [(assoc m a here) (conj pushed a)]
                                               [m pushed]))
                                           [m []]
                                           (:adjacent (corners id)))]
                    (recur (into (pop queue) pushed) m))
                  m))
        ;; The sea is not a source of fresh water, but a corner in it is
        ;; not dry either, and leaving it at zero puts desert on the beach.
        moist (reduce (fn [m v] (if (:ocean? v) (assoc m (:id v) 1.0) m)) moist corners)
        ;; Spread over the full range, for the same reason elevation is.
        land (vec (sort-by moist (into [] (comp (remove :ocean?) (remove :water?) (map :id)) corners)))
        n (count land)
        moist (reduce (fn [m [i id]]
                        (assoc m id (if (> n 1) (/ (double i) (dec n)) 0.0)))
                      moist
                      (map-indexed vector land))
        ;; Everything wet is simply wet. The breadth-first pass runs up to
        ;; 3.0 where rivers pile up, and the redistribution above only
        ;; covers land -- so without this a lake with a river through it
        ;; keeps a raw count, the cells around it average to more than
        ;; one, and the biome table, whose thresholds are all fractions,
        ;; reads off the end of itself.
        moist (reduce (fn [m v] (if (or (:ocean? v) (:water? v)) (assoc m (:id v) 1.0) m))
                      moist corners)
        corners (mapv #(assoc % :moisture (moist (:id %))) corners)]
    (assoc island
           :corners corners
           :centers (mapv (fn [c]
                            (let [cs (:corners c)]
                              (assoc c :moisture
                                     (if (seq cs)
                                       (/ (reduce + (map #(:moisture (corners %)) cs)) (count cs))
                                       0.0))))
                          centers))))

;; ---------------------------------------------------------------------------
;; Biomes

(defn biome
  "Whittaker's diagram, as a lookup on elevation and moisture.

  The one formula in the namespace, and it is a table. Everything it
  reads was arrived at by walking the graph."
  [{:keys [ocean? water? coast? elevation moisture]}]
  (let [e (double (or elevation 0.0))
        m (double (or moisture 0.0))]
    (cond
      ocean? :ocean
      water? (cond (< e 0.1) :marsh
                   (> e 0.8) :ice
                   :else :lake)
      coast? :beach
      (> e 0.8) (cond (> m 0.50) :snow
                      (> m 0.33) :tundra
                      (> m 0.16) :bare
                      :else :scorched)
      (> e 0.6) (cond (> m 0.66) :taiga
                      (> m 0.33) :shrubland
                      :else :temperate-desert)
      (> e 0.3) (cond (> m 0.83) :temperate-rain-forest
                      (> m 0.50) :temperate-deciduous-forest
                      (> m 0.16) :grassland
                      :else :temperate-desert)
      :else (cond (> m 0.66) :tropical-rain-forest
                  (> m 0.33) :tropical-seasonal-forest
                  (> m 0.16) :grassland
                  :else :subtropical-desert))))

(def biome-colors
  "Patel's palette, which is chosen so that the diagram reads as a
  diagram: the dry end warm, the wet end green, and altitude draining the
  colour out towards bare rock and snow."
  {:ocean [68 68 122] :lake [51 102 153] :marsh [47 102 102] :ice [153 255 255]
   :beach [160 144 119] :snow [255 255 255] :tundra [187 187 170]
   :bare [136 136 136] :scorched [85 85 85] :taiga [153 170 119]
   :shrubland [136 153 119] :temperate-desert [201 210 155]
   :temperate-rain-forest [68 136 85] :temperate-deciduous-forest [103 148 89]
   :grassland [136 170 85] :tropical-rain-forest [51 119 85]
   :tropical-seasonal-forest [85 153 68] :subtropical-desert [210 185 139]})

(defn- assign-biomes [{:keys [centers] :as island}]
  (assoc island :centers (mapv #(assoc % :biome (biome %)) centers)))

;; ---------------------------------------------------------------------------
;; The whole thing

(defn scatter
  "`n` points in `bounds`, drawn from `rng`."
  [n bounds rng]
  (let [[[x0 y0] [x1 y1]] bounds]
    (vec (repeatedly n #(vector (+ x0 (* (- x1 x0) (rng)))
                                (+ y0 (* (- y1 y0) (rng))))))))

(defn generate
  "A finished polygon map.

  Every pass in order, each one reading what the last one wrote:

      mesh -> water -> ocean -> elevation -> downslope -> watersheds
           -> rivers -> moisture -> biomes -> noisy edges

  Returns the `allgo.geometry.dual-mesh` value with the centers, corners
  and edges decorated -- so `dual-mesh/polygon` still draws a cell, and
  everything a renderer wants is on the cell it belongs to.

  Noisy edges come last, after every decision has been made, so that
  turning them off changes the picture and not the map: the same seed
  gives the same coastline either way, drawn straight or drawn
  wandering."
  ([] (generate {}))
  ([opts]
   (let [{:keys [points bounds relax shape seed rng lake-threshold river-attempts
                 noisy?]
          :as opts} (merge defaults opts)
         shape-fn (if (fn? shape)
                    shape
                    ((get shapes shape noise-shape) (assoc opts :seed seed)))]
     (-> (dm/relaxed (scatter points bounds rng) bounds relax)
         (assign-water shape-fn bounds lake-threshold)
         assign-ocean
         assign-corner-elevation
         redistribute-elevation
         assign-center-elevation
         assign-downslope
         assign-watersheds
         (assign-rivers rng river-attempts)
         assign-moisture
         assign-biomes
         (cond-> noisy? (dm/noisy-edges {:rng rng}))))))
