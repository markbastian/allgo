(ns allgo.astro.porkchop
  "Launch windows: when to leave one planet for another, and what it costs.

  Every pair of departure and arrival dates fixes one ballistic transfer,
  Lambert's orbit between where the planets are then, and with it a
  launch energy C3 and an arrival excess speed. Laid over a plane of
  departure against arrival dates, their contours make the porkchop plot
  of mission design (the name is for the contours' shape): two lobes of
  low energy -- type I, the short way round, under 180 degrees of
  heliocentric travel, and type II, the long way -- split by a ridge where
  the travel is near 180 degrees and the transfer plane has to stand up
  steeply to reach a planet a little out of the ecliptic. The pattern
  repeats every synodic period, as the planets come back to the same
  places relative to one another, but never quite the same, since the
  orbits are neither circular nor coplanar.

  This namespace draws the plane (`grid`), finds its minima (`optimum`)
  and lists the opportunities over years (`opportunities`), each against
  a choice of cost (`objectives`). It follows the way NASA's mission
  design handbooks build their tables -- Burke, Falck and McGuire,
  NASA/TM-2010-216764, for Earth to Mars -- and reproduces them.

  Dates are TT MJD; speeds km/s; C3 km^2/s^2."
  (:require [allgo.astro.interplanetary :as ip]
            [allgo.astro.planets :as planets]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.optimize :as opt]
            [clojure.math :as math]))

(defn synodic-period
  "Days for two planets to return to the same arrangement relative to the
  Sun: 1/|1/P1 - 1/P2| from their mean motions -- 780 days for the Earth
  and Mars, 584 for the Earth and Venus, 399 for the Earth and Jupiter."
  [a b]
  (let [rate #(second (get-in planets/elements [% :L]))]
    (/ (* 360.0 36525.0) (abs (- (rate a) (rate b))))))

(def objectives
  "Costs a transfer can be judged by, each a function of what
  `interplanetary/transfer` returns:

    :c3            the launch energy, which sets what a rocket can lift
    :v-inf-arrive  the arrival excess speed, which sets the braking at the
                   far end (or the heat shield's work)
    :v-inf-sum     the two excess speeds together, a first measure of the
                   whole trip"
  {:c3 :c3
   :v-inf-arrive #(v3/length (:v-inf-arrive %))
   :v-inf-sum #(+ (v3/length (:v-inf-depart %)) (v3/length (:v-inf-arrive %)))})

(defn- cost-fn [objective]
  (if (keyword? objective) (objectives objective) objective))

(defn- summary
  "What a porkchop cell keeps of a transfer."
  [t from to depart arrive type]
  (assoc t :from from :to to :depart depart :arrive arrive :tof (- arrive depart) :type type))

(defn- travel
  "Degrees of heliocentric travel from `r1` to `r2` the short way."
  [r1 r2]
  (math/to-degrees (math/acos (max -1.0 (min 1.0 (/ (v3/dot r1 r2) (* (v3/length r1) (v3/length r2))))))))

(defn transfer
  "The transfer from `from` at `depart` to `to` at `arrive` of `type`, 1
  (the short way) or 2 (the long way), or with `:best` whichever costs
  less by `objective`: `interplanetary/transfer`'s map with `:from :to
  :depart :arrive :tof :type` and `:travel`, the degrees of heliocentric
  travel, added; nil if none fits.

  `:ridge`, degrees (default 1), leaves out transfers within that much of
  180 degrees of travel. There the transfer plane is fixed only by how far
  the far planet sits out of the departure plane, and it stands up toward
  the pole at great cost -- except at the rare nodal transfers, where the
  far planet is crossing the ecliptic and the cost drops to almost nothing
  on a knife edge a day wide, too sensitive to fly (NASA/TM-2010-216764
  calls them undesirable). 0 keeps everything. Other options are
  `interplanetary/transfer`'s."
  ([from depart to arrive type] (transfer from depart to arrive type {}))
  ([from depart to arrive type {:keys [objective ephemeris ridge]
                                :or {objective :c3 ephemeris planets/heliocentric-state ridge 1.0}
                                :as opts}]
   (if (= type :best)
     (let [cost (cost-fn objective)
           ts (keep #(transfer from depart to arrive % opts) [1 2])]
       (when (seq ts) (apply min-key cost ts)))
     (let [angle (travel (first (ephemeris from depart)) (first (ephemeris to arrive)))
           angle (if (= type 2) (- 360.0 angle) angle)]
       (when (>= (abs (- angle 180.0)) ridge)
         (some-> (ip/transfer from depart to arrive (assoc opts :long? (= type 2)))
                 (summary from to depart arrive type)
                 (assoc :travel angle)))))))

(defn dates
  "Dates from `start` up to `end`, `step` days apart."
  [start end step]
  (vec (range start (+ end (* 0.5 step)) step)))

(defn grid
  "The porkchop plane: every transfer from `from` to `to` for each
  departure date in `departs` and arrival date in `arrives`, as
  `{:departs :arrives :cells}` with `:cells` a vector by arrival of
  vectors by departure -- rows up the plot, columns across it -- each a
  `transfer` or nil where none fits. Options: `:type` (default `:best`)
  and those of `transfer`."
  ([from to departs arrives] (grid from to departs arrives {}))
  ([from to departs arrives {:keys [type] :or {type :best} :as opts}]
   {:departs departs :arrives arrives
    :cells (mapv (fn [a] (mapv (fn [d] (transfer from d to a type opts)) departs)) arrives)}))

(defn contour-values
  "A grid's cells reduced to numbers by `objective` (a key of `objectives`
  or a function of a transfer), nil where there is no transfer -- the
  array a contour plotter wants."
  [{:keys [cells]} objective]
  (let [cost (cost-fn objective)]
    (mapv (fn [row] (mapv #(some-> % cost) row)) cells)))

(defn optimum
  "The least-cost transfer of `type` (1 or 2) near departure `depart` and
  arrival `arrive`, both dates free: the minimum of `objective` (default
  `:c3`) found by the simplex method from there. Options as `transfer`'s."
  ([from to depart arrive type] (optimum from to depart arrive type {}))
  ([from to depart arrive type {:keys [objective] :or {objective :c3} :as opts}]
   (let [cost (cost-fn objective)
         f (fn [[d a]] (some-> (transfer from d to a type opts) cost))
         {[d a] :x} (opt/nelder-mead f [depart arrive] {:step [2.0 2.0] :tol-x 1e-4 :tol-f 1e-9})]
     (transfer from d to a type opts))))

(defn- local-minima
  "Indices of the entries of `xs` (numbers, or nil for none) that are the
  least within `half` entries either side."
  [xs half]
  (let [n (count xs)]
    (filter (fn [i]
              (when-let [x (xs i)]
                (every? #(let [y (xs %)] (or (nil? y) (< x y) (and (= x y) (< i %))))
                        (remove #{i} (range (max 0 (- i half)) (min n (+ i half 1)))))))
            (range n))))

(defn opportunities
  "Every launch opportunity from `from` to `to` departing between `start`
  and `end`: for each synodic period and each type, the least-cost
  transfer, as `transfer`s sorted by departure. Found by scanning a grid
  and refining its minima with `optimum`. Options:

    :types      the transfer types to look for (default [1 2])
    :objective  what to minimize (default :c3)
    :tof        [shortest longest] flight in days (default 0.35 to 1.75
                times the Hohmann transfer's -- 90 to 450 days to Mars)
    :step       the scan's spacing in days, departure and flight time
                (default the synodic period / 200 and the Hohmann time /
                100)

  and those of `transfer`."
  ([from to start end] (opportunities from to start end {}))
  ([from to start end {:keys [types objective tof step] :or {types [1 2] objective :c3} :as opts}]
   (let [cost (cost-fn objective)
         syn (synodic-period from to)
         th (/ (:tof (ip/hohmann from to)) 86400.0)
         [tmin tmax] (or tof [(* 0.35 th) (* 1.75 th)])
         [ds dt] (or step [(max 1.0 (/ syn 200.0)) (max 1.0 (/ th 100.0))])
         ;; margin so a minimum near either end is still found whole
         departs (dates (- start (* 0.25 syn)) (+ end (* 0.25 syn)) ds)
         tofs (dates tmin tmax dt)]
     (->> (for [type types
                :let [best (mapv (fn [d]
                                   (let [ts (keep #(transfer from d to (+ d %) type opts) tofs)]
                                     (when (seq ts) (apply min-key cost ts))))
                                 departs)
                      costs (mapv #(some-> % cost) best)]
                i (local-minima costs (long (/ (* 0.5 syn) ds)))
                :let [{:keys [depart arrive]} (best i)
                      o (optimum from to depart arrive type opts)]
                :when (and o (<= start (:depart o) end))]
            o)
          (sort-by :depart)
          vec))))
