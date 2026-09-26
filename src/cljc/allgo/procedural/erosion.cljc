(ns allgo.procedural.erosion
  "Hydraulic erosion: rain on a heightmap, and what the water takes with it.

  Fractal terrain is self-similar everywhere, and that is exactly what is
  wrong with it. Real landscape is not: it has been sorted by something
  that only runs downhill. Ridges are sharp because water cannot climb
  them, valleys are flat-bottomed because that is where what the water
  carried had to stop, and the two are connected by a drainage network
  with a direction. `allgo.procedural.fractal`'s ridged multifractal fakes
  the look of that by folding the noise; this simulates the cause.

  ## Droplets, not a pipe model

  There are two ways to do this, and they answer different questions.

  The *pipe model* -- Mei, Decaudin and Hu (2007), following Musgrave,
  Kolb and Mace (1989) -- keeps a water depth, a sediment load and a
  velocity field on the grid, moves water between cells through virtual
  pipes, and advects sediment through the velocity field. It is shallow
  water with dirt in it. It gives lakes, deltas and flooding for free,
  because standing water is a state it can represent.

  This is the other one: Beyer (2015). A droplet is a particle with a
  position, a direction, a speed, a volume of water and a sediment load,
  and it walks downhill until it evaporates or leaves the map. Nothing is
  solved; the drainage network is not computed but *accumulated*, out of
  the statistics of many independent paths. Rivers appear where droplets
  agree, and they agree because the terrain they are cutting steers the
  ones that come after.

  Two reasons to take that road here. It is the cheaper one by a wide
  margin -- a droplet touches thirty cells over its whole life, where the
  grid solver touches every cell every step, and the interesting erosion
  takes thousands of steps. And this repository already has two grid water
  solvers in `allgo.physics.fluid` and `allgo.physics.height-field`; a
  third would mostly restate them. What is given up is real: no lakes, no
  standing water, no flood. A pit here is filled by the sediment of
  droplets that die trying to climb out of it, which is a plausible
  imitation of a lake silting up and is not a lake.

  ## What makes it carve rather than smooth

  One asymmetry does all the work. A droplet can hold

      capacity = max(slope, min-slope) * speed * water * k

  and slope, speed and water are not independent: water accelerates going
  down, so fast water on a steep grade can hold far more than it is
  carrying and takes the difference out of the bed. Where the grade eases
  the same droplet is suddenly over capacity and has to put the difference
  down. Erosion and deposition are the same rule read in two directions,
  and the `min-slope` floor is what stops a droplet on dead-flat ground
  from concluding it can carry nothing and dumping its whole load in one
  cell.

  Both halves are clamped by the height difference the droplet just
  crossed. It may not cut deeper than the drop it fell, nor fill higher
  than the step it is climbing -- without that the rule digs spikes and
  builds towers, since nothing else in it knows what a surface is.

  ## Why the two halves are not shaped alike

  Erosion is taken from a whole disk, `:radius` cells across, weighted by
  distance. A droplet's running deposits go into just the four cells
  beneath it, with the bilinear weights its height was read with.

  That asymmetry is deliberate and it is the difference between terrain
  and a bed of nails. Deposition is self-correcting -- fill a cell and the
  next droplet flows around it -- so it can afford to be sharp. Erosion is
  not: dig a cell and every later droplet is steered *into* it, digs it
  further, and the map ends up a field of one-cell wells with the drainage
  network lost between them. Spreading the cut over a disk is what makes a
  valley have sides.

  The exception is the last deposit a droplet makes. A droplet that runs
  out of water or out of lifetime still on the map sets down everything it
  is carrying at once, and that is far too much to put in four cells --
  it is a spike, and a quarter of a million of them read as grit strewn
  over the landscape. A drying puddle spreads, so that one goes over the
  disk as well.

  ## The edge of the map

  The border is not a feature of the landscape; it is where the map was
  cut. A droplet that runs off it is gone, and so is the sediment it was
  carrying -- which is right at a coastline and is a problem at an
  arbitrary crop, because the cells just inside the border are mined by
  every droplet that leaves and refilled by none.

  Left alone that is not a blemish but an instability. Capacity goes as
  slope times speed and speed goes as the square root of the drop, so
  capacity goes as the drop to the three halves, while the cut is limited
  only by the drop itself; every step therefore pays for a cut deeper than
  the step that prompted it. In the interior that feeds back into nothing,
  because what a droplet mines from a rim it sets down a few cells later.
  Against a border it compounds, and an edge pit is the one kind that
  never refills. Unchecked it dug a hole fifty times the depth of the
  terrain and still going, in a handful of droplets per cell.

  So erosion -- not deposition, which is free to build the coastal plain
  it would build anyway -- is faded out over the last `:border` cells,
  reaching zero at the edge. What that zeroes is a cut *centered* on a
  border cell, not one that reaches it: a droplet a few cells inside
  still swings its disk across the edge, so the outermost ring is worn
  lightly rather than not at all. Lightly is the point. Declining to cut
  what we cannot model the other side of costs a margin of barely
  weathered ground and buys a map without a trench around it.

  ## Units, and the one parameter to get right

  There is no length scale here. Slope is height per cell, so every
  parameter is read against the ratio of the heightmap's vertical units
  to its horizontal spacing, and none of them are physical constants.

  `:min-slope` is the one that will bite. It is the floor under the slope
  in the capacity formula, and it exists so that a droplet on dead-flat
  ground does not conclude it can carry nothing and drop its whole load in
  one cell. It has to sit *well below* the typical difference between
  neighboring cells, and the failure when it does not is quiet rather
  than loud: with the floor above the real slopes, every droplet
  everywhere has the same capacity, erosion stops depending on the terrain
  it is cutting, and what comes out is the heightmap smoothed rather than
  carved. It still looks like something happened. It is a blur.

  Worth measuring rather than guessing. A diamond-square grid from
  `allgo.procedural.terrain` at `:width 1.0` and 257 cells on a side has a
  median neighbor-to-neighbor difference near 0.0036; the 0.01 that is
  the usual published default puts the floor above 94% of the map, and
  gives exactly that blur. The default here is an order of magnitude
  below the median instead, and between one and two orders below it the
  result stops changing.

  ## What comes out

  A heightmap, and `:flux`: how much water crossed each cell, summed over
  every droplet that ever did. That is the drainage accumulation, and it
  is the thing to threshold if you want to know where the rivers are --
  the carved channel alone does not tell you, because a dry valley and a
  wet one have the same shape."
  (:require [allgo.array :as a]
            [allgo.procedural.shaping :as shaping]
            [clojure.math :as math]))

(def defaults
  "Parameters, with the meanings the docstring gives them.

  `:droplets` defaults to one per cell, which is enough to establish a
  drainage network and not enough to finish cutting it; several per cell
  is the usual range. `:rng` is a function of no arguments returning a
  number in [0, 1), and the hook that makes a run reproducible."
  {:droplets      nil
   ;; Steps, not seconds. A droplet moves one cell per step, so this is
   ;; how far it may travel -- and it is a cap on how far sediment can be
   ;; carried from where it was picked up.
   :lifetime      30
   ;; 0 follows the gradient exactly, which makes every droplet in a cell
   ;; take the same path and braids nothing; 1 ignores the terrain. Small
   ;; is right: enough to carry a droplet across a ridge it nearly met.
   :inertia       0.05
   :capacity      4.0
   ;; The floor under the slope in the capacity formula, in height units
   ;; per cell, and the parameter most worth getting right -- see the
   ;; note on units in the namespace docstring.
   :min-slope     0.0005
   :deposition    0.3
   :erosion       0.3
   :evaporation   0.01
   :gravity       4.0
   :initial-water 1.0
   :initial-speed 1.0
   ;; Cells. 1 is visibly pitted, 2 to 4 is the useful range, and large
   ;; radii cost the whole disk per erosion step for a softer cut.
   :radius        3
   ;; Cells of margin over which erosion fades to nothing at the map's
   ;; edge. Zero turns the taper off, and digs a trench around the map.
   :border        8
   :rng           rand})

;; ---------------------------------------------------------------------------
;; The erosion brush

(defn brush
  "The disk a droplet erodes from: offsets from the center cell, and the
  share of the cut each one takes.

  Weights fall linearly to zero at the rim and are normalized to sum to
  one, so the brush moves a fixed amount of material regardless of its
  size. Cells at the rim contribute nothing and are dropped rather than
  visited."
  [radius]
  (let [r  (long radius)
        rd (double radius)
        cells (vec (for [dr (range (- r) (inc r))
                         dc (range (- r) (inc r))
                         :let [d (math/sqrt (double (+ (* dr dr) (* dc dc))))
                               w (- 1.0 (/ d (+ rd 1e-9)))]
                         :when (> w 1e-12)]
                     [dr dc w]))
        total (reduce + 0.0 (map peek cells))]
    {:n (count cells)
     :dr (a/i32 (map first cells))
     :dc (a/i32 (map second cells))
     :w  (a/f64 (map #(/ (double (peek %)) total) cells))}))

(defn- spread-at!
  "Adds `amount` across the disk around cell `ri ci`, in the brush's
  proportions, and returns how much actually landed. Negative cuts.

  Less than asked for where the disk hangs off the grid, and the caller
  is told so -- for a cut that means the droplet is credited only with
  what it really took, so mass is conserved and erosion fades out toward
  the border instead of being concentrated there.

  The other reading, squeezing the whole cut into whatever part of the
  disk is in bounds, is what the first version of this did. It conserves
  mass too, and it is unstable: it cuts an edge cell at up to twice the
  rate of an interior one, and an edge pit is the one kind that never
  refills."
  ;; Primitives are coerced in the body rather than hinted on the
  ;; parameters: a fn taking primitives is limited to four arguments.
  [^doubles heights rows cols disk ri ci amount]
  (let [rows   (long rows)
        cols   (long cols)
        ri     (long ri)
        ci     (long ci)
        amount (double amount)
        n      (long (:n disk))
        ^ints bdr   (:dr disk)
        ^ints bdc   (:dc disk)
        ^doubles bw (:w disk)]
    (loop [k 0 applied 0.0]
      (if (= k n)
        applied
        (let [r (+ ri (aget bdr k))
              c (+ ci (aget bdc k))]
          (if (and (<= 0 r) (< r rows) (<= 0 c) (< c cols))
            (let [id    (+ (* r cols) c)
                  share (* amount (aget bw k))]
              (aset heights id (+ (aget heights id) share))
              (recur (inc k) (+ applied share)))
            (recur (inc k) applied)))))))

(defn- deposit-at!
  "Puts `amount` into the four cells around the droplet, with the bilinear
  weights its height was read with.

  No bounds check: a droplet is only ever inside the grid's last row and
  column, which `run-droplet!` maintains, so all four exist."
  [^doubles heights cols ri ci fr fc amount]
  (let [cols (long cols)
        fr (double fr)
        fc (double fc)
        amount (double amount)
        id (+ (* (long ri) cols) (long ci))
        add! (fn [^long i ^double w]
               (aset heights i (+ (aget heights i) (* amount w))))]
    (add! id             (* (- 1.0 fr) (- 1.0 fc)))
    (add! (+ id 1)       (* (- 1.0 fr) fc))
    (add! (+ id cols)    (* fr (- 1.0 fc)))
    (add! (+ id cols 1)  (* fr fc))))

;; ---------------------------------------------------------------------------
;; One droplet

(defn- run-droplet!
  "Walks one droplet downhill from `r0 c0` until it dries up, runs out of
  lifetime, or leaves the map.

  A droplet that leaves the map takes its sediment with it, which is what
  a river does at the coast. One that stops on the map puts its load down
  where it stopped. Between the two, nothing is created or destroyed."
  [^doubles heights ^doubles flux rows cols disk params r0 c0]
  (let [rows        (long rows)
        cols        (long cols)
        lifetime    (long (:lifetime params))
        inertia     (double (:inertia params))
        cap-k       (double (:capacity params))
        min-slope   (double (:min-slope params))
        deposition  (double (:deposition params))
        erosion     (double (:erosion params))
        evaporation (double (:evaporation params))
        gravity     (double (:gravity params))
        border      (double (:border params))
        rng         (:rng params)
        max-r       (- rows 1.0)
        max-c       (- cols 1.0)]
    (loop [step  0
           r     (double r0)
           c     (double c0)
           dr    0.0
           dc    0.0
           speed (double (:initial-speed params))
           water (double (:initial-water params))
           sed   0.0]
      (let [ri (long (math/floor r))
            ci (long (math/floor c))
            fr (- r ri)
            fc (- c ci)
            id (+ (* ri cols) ci)
            ;; Erosion is faded out toward the border. The border is not
            ;; a feature of the landscape, it is where the map was cut,
            ;; and a droplet that runs off it exports its load for good --
            ;; so the cells beside it are mined by every droplet that
            ;; leaves and refilled by none, and they sink. Declining to
            ;; cut what we cannot model the other side of costs a margin
            ;; of unweathered terrain and buys a map without a trench
            ;; around it.
            taper (if (pos? border)
                    (shaping/smoothstep
                     0.0 border
                     (double (min ri ci (- rows 1 ri) (- cols 1 ci))))
                    1.0)]
        (if (or (>= step lifetime) (< water 1e-6))
          ;; Out of time or out of water, and still on the map: what it
          ;; was carrying settles where it stopped. Water that evaporates
          ;; does not take its sediment with it -- that is how an alluvial
          ;; fan happens, and letting the load vanish instead is not a
          ;; rounding error. Every droplet picks up close to a full load
          ;; and most of them end this way, so dropping it on the floor
          ;; cost the terrain four fifths of its material within four
          ;; droplets per cell, planing the map toward flat.
          ;; Over the whole disk, not the four cells under the droplet.
          ;; An incremental deposit is small and self-correcting -- fill a
          ;; cell and the next droplet flows around it -- so it can afford
          ;; to be sharp. This one is a whole load at once, and put down
          ;; sharply it is a spike; a quarter of a million of them read as
          ;; grit strewn over the terrain rather than as sediment. A
          ;; drying puddle spreads, and so does this.
          (when (pos? sed) (spread-at! heights rows cols disk ri ci sed))
          (let [h00 (aget heights id)
                h01 (aget heights (+ id 1))
                h10 (aget heights (+ id cols))
                h11 (aget heights (+ id cols 1))
                h   (+ (* h00 (- 1.0 fr) (- 1.0 fc))
                       (* h01 (- 1.0 fr) fc)
                       (* h10 fr (- 1.0 fc))
                       (* h11 fr fc))
                ;; The bilinear surface's exact gradient, not a difference
                ;; of samples: the droplet is steered by the same surface
                ;; its height is read from, so it cannot be told it is
                ;; going downhill by one and uphill by the other.
                gr  (+ (* (- h10 h00) (- 1.0 fc)) (* (- h11 h01) fc))
                gc  (+ (* (- h01 h00) (- 1.0 fr)) (* (- h11 h10) fr))
                ;; Momentum against gradient, unnormalized first: the two
                ;; can cancel exactly -- a droplet that has run itself
                ;; onto a flat spot -- and then there is no downhill to
                ;; speak of and a direction has to be invented.
                ndr (- (* dr inertia) (* gr (- 1.0 inertia)))
                ndc (- (* dc inertia) (* gc (- 1.0 inertia)))
                len (math/sqrt (+ (* ndr ndr) (* ndc ndc)))
                flat? (< len 1e-12)
                ang (if flat? (* 2.0 math/PI (double (rng))) 0.0)
                udr (if flat? (math/sin ang) (/ ndr len))
                udc (if flat? (math/cos ang) (/ ndc len))
                nr  (+ r udr)
                nc  (+ c udc)]
            (aset flux id (+ (aget flux id) water))
            ;; Strictly inside the last row and column: a droplet's height
            ;; is read from the four cells starting at its own, so the
            ;; last of each has no cell beyond it to pair with.
            (when-not (or (< nr 0.0) (>= nr max-r) (< nc 0.0) (>= nc max-c))
              (let [nri (long (math/floor nr))
                    nci (long (math/floor nc))
                    nfr (- nr nri)
                    nfc (- nc nci)
                    nid (+ (* nri cols) nci)
                    nh  (+ (* (aget heights nid) (- 1.0 nfr) (- 1.0 nfc))
                           (* (aget heights (+ nid 1)) (- 1.0 nfr) nfc)
                           (* (aget heights (+ nid cols)) nfr (- 1.0 nfc))
                           (* (aget heights (+ nid cols 1)) nfr nfc))
                    dh  (- nh h)
                    cap (* (max (- dh) min-slope) speed water cap-k)
                    sed' (if (or (pos? dh) (< cap sed))
                           ;; Climbing, or over capacity. Climbing wins: a
                           ;; droplet heading uphill is in a pit, and what
                           ;; it drops is what fills the pit -- capped at
                           ;; the step it is climbing, so it settles level
                           ;; with the rim rather than over it.
                           (let [amt (if (pos? dh)
                                       (min dh sed)
                                       (* (- sed cap) deposition))]
                             (deposit-at! heights cols ri ci fr fc amt)
                             (- sed amt))
                           ;; Under capacity on a descent. `(- dh)` is the
                           ;; drop, and the cut may not exceed it. The
                           ;; droplet is credited with what the disk
                           ;; actually gave up, not what it asked for.
                           (let [amt (* taper (min (* (- cap sed) erosion) (- dh)))
                                 ;; Credited with what the disk actually
                                 ;; gave up, not with what it was asked
                                 ;; for. `spread-at!` signs its answer the
                                 ;; way it signs its argument, so a cut
                                 ;; comes back negative.
                                 cut (- (spread-at! heights rows cols disk ri ci (- amt)))]
                             (+ sed cut)))]
                (recur (inc step) nr nc udr udc
                       ;; Kinetic energy from the drop. Clamped because a
                       ;; droplet can be carried uphill by its momentum
                       ;; further than that momentum can pay for.
                       (math/sqrt (max 0.0 (+ (* speed speed) (* (- dh) gravity))))
                       (* water (- 1.0 evaporation))
                       sed')))))))))

;; ---------------------------------------------------------------------------
;; Rain

(defn- copy-doubles ^doubles [^doubles src]
  (let [n (alength src)
        ^doubles dst (a/f64 n)]
    (dotimes [i n] (aset dst i (aget src i)))
    dst))

(defn- dimensions
  "`[rows cols]` of a field, from `:rows`/`:cols` or a square `:dim`."
  [{:keys [dim rows cols]}]
  (let [rows (long (or rows dim))
        cols (long (or cols dim))]
    (when (or (< rows 2) (< cols 2))
      (throw (ex-info "erosion needs a grid at least 2x2"
                      {:rows rows :cols cols})))
    [rows cols]))

(defn erode!
  "Rains on the heightmap in place, and returns the field with `:flux`.

  The field is anything with `:heights` -- a flat `double` array indexed
  `row * cols + col` -- and either a square `:dim` or explicit `:rows` and
  `:cols`. `allgo.procedural.terrain`'s grids are already that shape, so

      (-> (terrain/generate {:width 1.0 :iterations 9}) erode!)

  works with nothing in between. `:flux` accumulates across calls if one
  is already there, so erosion can be run in rounds and still report the
  drainage of all of them."
  ([field] (erode! field {}))
  ([{:keys [heights] :as field} opts]
   (let [params    (merge defaults opts)
         [rows cols] (dimensions field)
         rng       (:rng params)
         droplets  (long (or (:droplets params) (* rows cols)))
         b         (brush (:radius params))
         ^doubles flux (or (:flux field) (a/f64 (* rows cols)))
         ;; Half-open: a droplet starting exactly on the last row or
         ;; column has no cell beyond it to interpolate against.
         span-r    (- rows 1.0)
         span-c    (- cols 1.0)]
     (dotimes [_ droplets]
       (run-droplet! heights flux rows cols b params
                     (* (double (rng)) span-r)
                     (* (double (rng)) span-c)))
     (assoc field :flux flux))))

(defn erode
  "`erode!` on a copy, leaving the field it was given alone.

  The same bargain `allgo.procedural.terrain/step` makes: the array inside
  is a concession to speed, and a caller holding a heightmap should not
  find it changed underneath."
  ([field] (erode field {}))
  ([field opts]
   (-> field
       (assoc :heights (copy-doubles (:heights field)))
       (dissoc :flux)
       (erode! opts))))

(defn drainage
  "How much water crossed cell `i j`, summed over every droplet.

  Zero where no droplet ever went. High along the channels, and rising
  downstream -- this is drainage accumulation, so thresholding it is how
  you decide what counts as a river."
  ^double [{:keys [flux dim cols]} i j]
  (aget ^doubles flux (+ (* (long i) (long (or cols dim))) (long j))))

(defn drainage-bounds
  "`[lowest highest]` flux over the whole field.

  Flux is unbounded above and heavily skewed -- a main channel carries
  orders of magnitude more than a hillside -- so this is mostly the
  denominator for a log scale rather than a linear one."
  [{:keys [flux]}]
  (let [^doubles f flux
        n (alength f)]
    (loop [i 0 lo ##Inf hi ##-Inf]
      (if (= i n)
        [lo hi]
        (let [v (aget f i)]
          (recur (inc i) (min lo v) (max hi v)))))))
