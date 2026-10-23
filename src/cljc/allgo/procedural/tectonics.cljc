(ns allgo.procedural.tectonics
  "Plate tectonics for a procedural planet: continents and mountain ranges
  that are where they are for a reason.

  `allgo.procedural.planet/default-terrain` puts land wherever a
  low-frequency noise happens to be high and mountains wherever a ridged
  one is, which makes a planet that is plausible everywhere and
  explicable nowhere. On a real one the mountains stand in lines along
  the edges of continents, the oceans have trenches beside island arcs
  and a ridge down their middles, and every one of those is a place where
  two plates meet. This builds the plates and lets their meetings decide
  the ground. The method is Amit Patel's, from his notes on procedural
  planet generation, with the plates moved as rigid bodies on the sphere.

  ## The steps

      cells       a geodesic sphere: every vertex a cell, its mesh edges
                  the adjacencies. Nearly uniform, with no poles and no
                  seam -- `allgo.geometry.tri-mesh/geodesic`.
      plates      seeds grown outward together, each step taking a random
                  cell off the frontier rather than the oldest, so that the
                  plates come out ragged instead of as Voronoi polygons
      continents  a share of the plates, by area, made continental -- thick,
                  light crust that floats high -- and the rest oceanic
      motion      each plate turning about its own Euler pole, which is
                  the only way a rigid plate can move on a sphere: a cell's
                  velocity is the pole crossed with its position, scaled
      boundaries  every edge between two plates, judged by how fast the
                  plates close across it: converging, pulling apart, or
                  sliding past
      ground      a base level for each kind of crust, falling from
                  continent to shelf to abyss away from the coast, and then
                  what each boundary builds, fading with distance from it

  What a converging boundary builds depends on what meets there:

      continent and continent    a range, the Himalaya
      continent and ocean        a coastal range on the continent and a
                                 trench where the ocean goes under, the
                                 Andes and the Peru-Chile Trench
      ocean and ocean            an arc of islands on one and a trench on
                                 the other, the Marianas
      and pulling apart:
      ocean                      a mid-ocean ridge
      continent                  a rift valley

  ## What it is not

  A picture of the present, not a history: nothing moves, nothing
  collides over time, nothing erodes. The motions are used once, to
  decide what each boundary is, and the ground is built from that. It
  gets the arrangement right -- mountains along margins, trenches beside
  arcs -- and leaves the shapes to the fractal detail laid over it."
  (:require [allgo.geometry.tri-mesh :as tri-mesh]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.random :as random]
            [clojure.math :as math]))

(def defaults
  {:resolution 5
   :plates 14
   :continental 0.35
   :seed 1})

;; ---------------------------------------------------------------------------
;; Vectors, as plain triples

(defn- dot ^double [[ax ay az] [bx by bz]]
  (+ (* (double ax) (double bx)) (* (double ay) (double by)) (* (double az) (double bz))))

(defn- cross [[ax ay az] [bx by bz]]
  (let [ax (double ax) ay (double ay) az (double az)
        bx (double bx) by (double by) bz (double bz)]
    [(- (* ay bz) (* az by)) (- (* az bx) (* ax bz)) (- (* ax by) (* ay bx))]))

(defn- sub [a b] (mapv - a b))

(defn- dot-any
  "`dot` for vectors of any length: the plates work on a flat map's points
  as well as the sphere's."
  ^double [a b]
  (reduce + (map * a b)))

(defn- unit [v]
  (let [l (math/sqrt (dot-any v v))]
    (mapv #(/ (double %) l) v)))

(defn- random-unit [rng]
  (let [z (random/uniform rng -1.0 1.0)
        a (random/uniform rng 0.0 (* 2.0 math/PI))
        r (math/sqrt (- 1.0 (* z z)))]
    [(* r (math/cos a)) (* r (math/sin a)) z]))

;; ---------------------------------------------------------------------------
;; The cells

(defn- cells
  "The geodesic sphere as a graph: `:points`, unit vectors; `:neighbors`,
  each cell's adjacent cells; `:tris`, the triangles as index triples,
  turned so that each runs counterclockwise seen from outside; and
  `:across`, for each triangle, the triangle beyond each of its three
  edges (edge `j` running from corner `j` to corner `j + 1`)."
  [resolution]
  (let [{:keys [verts tri-ids]} (tri-mesh/geodesic resolution)
        points (mapv vec (partition 3 verts))
        tris (mapv (fn [[a b c]]
                     (if (pos? (dot (cross (points a) (points b)) (points c)))
                       [a b c]
                       [a c b]))
                   (partition 3 tri-ids))
        n (count points)
        key-of (fn [a b] (+ (* (long (min a b)) n) (long (max a b))))
        by-edge (persistent!
                 (reduce (fn [m [t [a b c]]]
                           (reduce (fn [m [j u w]] (assoc! m (key-of u w) (conj (get m (key-of u w) []) [t j])))
                                   m [[0 a b] [1 b c] [2 c a]]))
                         (transient {})
                         (map-indexed vector tris)))
        across (mapv (fn [[t [a b c]]]
                       (mapv (fn [[u w]]
                               (some (fn [[t2 _]] (when (not= t2 t) t2)) (by-edge (key-of u w))))
                             [[a b] [b c] [c a]]))
                     (map-indexed vector tris))
        neighbors (persistent!
                   (reduce (fn [nb [a b c]]
                             (-> nb
                                 (assoc! a (conj (nth nb a) b c))
                                 (assoc! b (conj (nth nb b) c a))
                                 (assoc! c (conj (nth nb c) a b))))
                           (transient (vec (repeat n #{})))
                           tris))]
    {:points points
     :neighbors (mapv vec neighbors)
     :tris tris
     :across across
     ;; the angle across one edge, which turns a count of steps into a
     ;; distance on the unit sphere
     :spacing (math/acos (dot (points 0) (points (first (nth neighbors 0)))))}))

;; ---------------------------------------------------------------------------
;; The plates

(defn- grow-plates
  "Each cell's plate: `k` seeds grown out together, a random frontier cell
  at a time."
  [neighbors k rng]
  (let [n (count neighbors)
        plate (int-array n -1)
        seeds (take k (random/shuffle rng (range n)))]
    (doseq [[p s] (map-indexed vector seeds)] (aset plate (long s) (int p)))
    (loop [frontier (vec seeds)]
      (when (seq frontier)
        (let [i (random/below rng (count frontier))
              c (nth frontier i)
              frontier (pop (assoc frontier i (peek frontier)))
              p (aget plate (long c))
              fresh (filter #(neg? (aget plate (long %))) (nth neighbors c))]
          (doseq [f fresh] (aset plate (long f) p))
          (recur (into frontier fresh)))))
    (vec plate)))

(defn- plate-motions
  "Each plate's crust and motion: `{:continental? :pole :spin}`, the plate
  turning about the unit vector `:pole` at `:spin` radians per unit time;
  or on a flat map `{:continental? :velocity}`, the plate drifting, at up
  to a unit a unit time.

  Continents are taken in a random order until they cover `share` of the
  surface, so the land fraction is about what was asked whatever sizes
  the plates came out."
  [plate k share rng flat?]
  (let [area (frequencies plate)
        n (count plate)
        order (random/shuffle rng (range k))
        continental (loop [[p & more] order covered 0 chosen #{}]
                      (if (or (nil? p) (>= (/ (double covered) n) (double share)))
                        chosen
                        (recur more (+ covered (long (area p 0))) (conj chosen p))))]
    (mapv (fn [p]
            (if flat?
              (let [a (random/uniform rng 0.0 (* 2.0 math/PI))
                    s (random/uniform rng 0.0 1.0)]
                {:continental? (contains? continental p)
                 :velocity [(* s (math/cos a)) (* s (math/sin a))]})
              {:continental? (contains? continental p)
               :pole (random-unit rng)
               :spin (random/uniform rng -1.0 1.0)}))
          (range k))))

(defn- velocity
  "How fast the plate is carrying the point `x`: turning it about the
  plate's pole on the sphere, or carrying it along on a flat map."
  [{:keys [pole spin] :as motion} x]
  (if pole
    (mapv #(* (double spin) (double %)) (cross pole x))
    (:velocity motion)))

(defn- crust
  "Whether each cell's crust is continental.

  Not simply its plate's: a plate carries ocean floor and continent
  together, and most coasts on the Earth -- every one around the Atlantic
  -- are in the middle of a plate rather than at its edge. A continental
  plate makes its cells likely to be continent and an oceanic one
  unlikely, and a broad noise field decides between them, so continents
  spill past their plates and fall short of them, and now and then an
  oceanic plate carries a sliver of continent of its own. The threshold
  is set by rank, so the share of continent is `share` exactly."
  [{:keys [points]} plate motions share seed]
  (let [field (-> (fractal/fbm (noise/gradient-basis {:seed (+ 31 (long seed))})
                               {:octaves 4.0 :H 0.9})
                  (noise/scaled 1.6))
        ;; a flat map's points sample the field on a plane through it
        score (mapv (fn [[x y z] p]
                      (+ (if (:continental? (motions p)) 0.5 -0.5)
                         (* 1.2 (noise/sample field (double x) (double y) (double (or z 0.5))))))
                    points plate)
        n (count score)
        cut (nth (vec (sort score)) (min (dec n) (long (* n (- 1.0 (double share))))))]
    (mapv #(>= (double %) (double cut)) score)))

;; ---------------------------------------------------------------------------
;; Boundaries

(defn boundary-kind
  "What kind of boundary runs between cell `a` on a plate moving at `va`
  and cell `b` on one moving at `vb`, and how hard: `[kind strength]`.

  The plates' relative velocity resolved along the line from one cell to
  the other is how fast they close; across it, how fast they slide. A
  boundary is converging or diverging when closing is at least half of
  the motion, and a transform when sliding dominates -- then it builds
  nothing. Strength is the closing speed, in the same units as the
  plates' spins."
  [a b va vb]
  (let [d (unit (sub b a))
        r (sub va vb)
        closing (dot-any r d)
        speed (math/sqrt (dot-any r r))]
    (cond
      (< speed 1e-12) [:transform 0.0]
      (> closing (* 0.5 speed)) [:converging closing]
      (< closing (* -0.5 speed)) [:diverging (- closing)]
      :else [:transform 0.0])))

(def ^:private features
  "What each kind of boundary builds on the cell beside it: a height, in
  the same units as the base levels below, and a width in radians on the
  unit sphere -- on the Earth a radian is 6400 km, so a range here is a
  few hundred kilometers across and a trench a hundred."
  {:range  {:height 0.9 :width 0.05}
   :arc    {:height 0.7 :width 0.03}
   :trench {:height -0.5 :width 0.02}
   :ridge  {:height 0.3 :width 0.04}
   :rift   {:height -0.3 :width 0.03}})

(defn- feature-of
  "What a boundary builds on cell `a`'s side of it, given which of the two
  plates is continental: one of the keys of `features`, or nil.

  Where ocean meets ocean one has to go under the other, and which is
  not something this model knows -- real plates decide it by age -- so
  the plate with the higher number rides over."
  [kind a-cont? b-cont? a-plate b-plate]
  (case kind
    :converging (cond
                  (and a-cont? b-cont?) :range
                  a-cont? :range
                  b-cont? :trench
                  (> (long a-plate) (long b-plate)) :arc
                  :else :trench)
    :diverging (cond
                 (and a-cont? b-cont?) :rift
                 (or a-cont? b-cont?) nil
                 :else :ridge)
    nil))

;; ---------------------------------------------------------------------------
;; Distances over the graph

(defn- nearest
  "From each source cell outward, the number of steps to the nearest
  source and that source's value: `[steps values]`, unreached cells
  at `##Inf` and nil."
  [neighbors sources]
  (let [n (count neighbors)
        steps (double-array n ##Inf)
        values (object-array n)]
    (doseq [[c v] sources]
      (aset steps (long c) 0.0)
      (aset values (long c) v))
    (loop [ring (mapv first sources) d 1.0]
      (when (seq ring)
        (let [next-ring (persistent!
                         (reduce (fn [acc c]
                                   (reduce (fn [acc f]
                                             (if (= ##Inf (aget steps (long f)))
                                               (do (aset steps (long f) d)
                                                   (aset values (long f) (aget values (long c)))
                                                   (conj! acc f))
                                               acc))
                                           acc (nth neighbors c)))
                                 (transient []) ring))]
          (recur next-ring (inc d)))))
    [(vec steps) (vec values)]))

;; ---------------------------------------------------------------------------
;; The ground

(defn- find-boundaries
  "Every edge between two plates, once, as `{:a :b :kind :strength}` with
  `a` the lower cell.

  An edge is judged against the boundary's direction where it is, not
  its own. The cells are a geodesic grid, so a boundary between them
  zigzags, each step up to sixty degrees off the line the boundary
  actually takes -- and sixty degrees is exactly the angle at which
  closing falls to half the motion, the line between converging and
  sliding past. Judged edge by edge, one straight boundary came out red
  and yellow by turns. So the direction across is from the middle of the
  nearby cells on one plate to the middle of those on the other, two
  rings out from the edge."
  [{:keys [points neighbors]} plate motions]
  (let [ring (fn [c] (into #{c} (mapcat #(cons % (nth neighbors %))) (nth neighbors c)))
        middle (fn [cs] (mapv #(/ (double %) (count cs)) (reduce #(mapv + %1 %2) (map points cs))))]
    (vec (for [a (range (count points))
               b (nth neighbors a)
               :when (< (long a) (long b))
               :let [pa (plate a) pb (plate b)]
               :when (not= pa pb)
               :let [near (into (ring a) (ring b))
                     [kind strength] (boundary-kind (middle (filter #(= pa (plate %)) near))
                                                    (middle (filter #(= pb (plate %)) near))
                                                    (velocity (motions pa) (points a))
                                                    (velocity (motions pb) (points b)))]]
           {:a a :b b :kind kind :strength strength}))))

(defn- elevations
  "Each cell's height, sea level at zero: the crust's base level and every
  boundary feature near it."
  [{:keys [points neighbors spacing]} plate crust boundaries]
  (let [spacing (double spacing)
        cont? (fn [c] (nth crust c))
        ;; what each boundary builds, on each side of it -- the kind of a
        ;; boundary is the same from either side, what it builds is not
        found (for [{:keys [a b kind strength]} boundaries
                    [x y] [[a b] [b a]]
                    :let [f (feature-of kind (cont? x) (cont? y) (plate x) (plate y))]
                    :when f]
                [f x strength])
        by-feature (group-by first found)
        ;; how far each cell is from crust of the other kind: the coast,
        ;; for land, and the foot of the continental slope, for sea
        [to-other] (nearest neighbors (for [a (range (count points))
                                            :when (some #(not= (cont? a) (cont? %)) (nth neighbors a))]
                                        [a true]))
        ;; strengths are normalized by the fastest closing anywhere, so
        ;; that the boundary with most going on builds the full height
        top (max 1e-9 (reduce max 0.0 (map #(nth % 2) found)))
        fields (into {}
                     (for [[f entries] by-feature]
                       (let [best (reduce (fn [m [_ a s]] (update m a (fnil max 0.0) s)) {} entries)]
                         [f (nearest neighbors (map (fn [[a s]] [a (/ (double s) top)]) best))])))]
    (mapv (fn [c]
            (let [d (* spacing (double (nth to-other c)))
                  base (if (cont? c)
                         ;; shelf at the coast, rising to the interior
                         (- 0.25 (* 0.33 (math/exp (- (/ d 0.06)))))
                         ;; the continental slope down to the abyss
                         (+ -0.6 (* 0.45 (math/exp (- (/ d 0.05))))))]
              (reduce (fn [e [f [steps values]]]
                        (let [s (nth steps c)]
                          (if (= ##Inf s)
                            e
                            (let [{:keys [height width]} (features f)
                                  x (/ (* spacing (double s)) (double width))]
                              (+ e (* (double height) (double (nth values c)) (math/exp (- (* x x)))))))))
                      base fields)))
          (range (count points)))))

;; ---------------------------------------------------------------------------
;; Assembling

(defn on-graph
  "Plates over any graph of cells, and the ground they make.

  `graph` is `{:points :neighbors :spacing}`: each cell's position --
  unit vectors on the sphere, or points of a flat map, which should run
  about -1 to 1 -- its adjacent cells, and the distance one step between
  cells stands for, in radians of a planet the size of the Earth's, which
  is what the features' widths are measured in. `:flat?` moves the plates
  by drifting rather than turning, as a map that is not a sphere needs.

  Returns `:plate`, `:motions`, `:continental`, `:boundaries` and
  `:elevation`, as `plates` describes."
  [graph {:keys [seed continental flat?] k :plates}]
  (let [rng (random/rng seed)
        plate (grow-plates (:neighbors graph) k rng)
        motions (plate-motions plate k continental rng flat?)
        crust (crust graph plate motions continental seed)
        boundaries (find-boundaries graph plate motions)]
    {:plate plate
     :motions motions
     :continental crust
     :boundaries boundaries
     :elevation (elevations graph plate crust boundaries)}))

(defn plates
  "A planet's plates and the ground they make.

  * `:resolution` -- the geodesic sphere's subdivisions: 5 is ten
    thousand cells, 6 forty thousand.
  * `:plates` -- how many plates.
  * `:continental` -- the share of the surface that is continental crust.
  * `:seed` -- makes a different world.

  Returns the cells (`:points`, `:neighbors`, `:tris`, `:across`), each
  cell's `:plate` and whether its crust is `:continental`, each plate's
  `:motions`, the `:boundaries` between plates as `{:a :b :kind
  :strength}`, and each cell's `:elevation`, sea level at zero.

  `:continental` is a share of crust, not of land: the shelves are
  continental crust under water, so somewhat less than that stands dry."
  ([] (plates {}))
  ([opts]
   (let [opts (merge defaults opts)
         mesh (cells (:resolution opts))]
     (merge mesh (on-graph mesh opts)))))

;; ---------------------------------------------------------------------------
;; Asking about any direction

(defn- walk
  "The triangle holding unit vector `p`, walked to from triangle `t`
  across whichever edge `p` lies beyond."
  ^long [{:keys [points tris across]} p t]
  (loop [t (long t) steps 0]
    (let [[a b c] (nth tris t)
          pa (points a) pb (points b) pc (points c)
          s0 (dot p (cross pa pb))
          s1 (dot p (cross pb pc))
          s2 (dot p (cross pc pa))
          worst (min s0 s1 s2)]
      (if (or (>= worst -1e-12) (> steps 100000))
        t
        (recur (long (nth (nth across t) (cond (= worst s0) 0 (= worst s1) 1 :else 2)))
               (inc steps))))))

(defn locator
  "A function of a unit vector giving `[[i j k] [wi wj wk]]`: the corners
  of the triangle it falls in and its barycentric weights there.

  Walks from wherever the last question was answered, which for a mesh
  being swept vertex by vertex is a step or two away."
  [world]
  (let [last-t (volatile! 0)]
    (fn [p]
      (let [t (walk world p @last-t)
            [a b c :as corners] (nth (:tris world) t)
            {:keys [points]} world
            wa (dot p (cross (points b) (points c)))
            wb (dot p (cross (points c) (points a)))
            wc (dot p (cross (points a) (points b)))
            sum (+ wa wb wc)]
        (vreset! last-t t)
        [corners [(/ wa sum) (/ wb sum) (/ wc sum)]]))))

(defn plate-at
  "Which plate is under unit vector `p`: the plate of the nearest corner of
  its triangle."
  [world locate p]
  (let [[corners weights] (locate p)]
    ((:plate world) (nth corners (apply max-key #(nth weights %) [0 1 2])))))

(defn terrain
  "The plates' ground as a terrain basis for `allgo.procedural.planet`: a
  `^double` function of a direction, zero at sea level.

  Between cells the ground is interpolated across each triangle, and over
  it go two fractals. A broad one moves the shoreline: without it every
  coast would trace a plate boundary, where real continents overrun
  their plates' edges and fall short of them in bays, peninsulas and
  shelf seas. And a fine one, rougher the higher the ground stands --
  the plates say where the mountains are and the noise what they look
  like. `:detail` scales both; `:sample-spacing` is what
  `allgo.procedural.planet/default-terrain` takes it to mean."
  ([world] (terrain world {}))
  ([world {:keys [seed detail sample-spacing] :or {seed 0 detail 1.0}}]
   (let [locate (locator world)
         elevation (:elevation world)
         layer (fn [s frequency octaves]
                 (-> (fractal/fbm (noise/gradient-basis {:seed (+ s (long seed))})
                                  {:octaves (if sample-spacing
                                              (min octaves (fractal/octaves-for
                                                            (* (double sample-spacing) frequency) 2.0))
                                              octaves)
                                   :H 0.8})
                     (noise/scaled frequency)))
         coast (layer 23 2.0 4.0)
         fine (layer 17 6.0 7.0)
         detail (double detail)]
     (fn ^double [^double x ^double y ^double z]
       (let [l (math/sqrt (+ (* x x) (* y y) (* z z)))
             p [(/ x l) (/ y l) (/ z l)]
             [[i j k] [wi wj wk]] (locate p)
             e (+ (* (double wi) (double (elevation i)))
                  (* (double wj) (double (elevation j)))
                  (* (double wk) (double (elevation k))))]
         (+ e (* detail (+ (* 0.3 (noise/sample coast x y z))
                           (* (noise/sample fine x y z) (+ 0.12 (* 0.35 (max 0.0 e))))))))))))
