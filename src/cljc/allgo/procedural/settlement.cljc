(ns allgo.procedural.settlement
  "Where people would live on a polygon map, what they build between
  those places, and how much of the land each one ends up holding.

  Takes a finished `allgo.procedural.island` and adds the human layer
  that Azgaar's generator is mostly made of. It is a separate namespace
  because it is a separate question: the island decides what the land is
  like, and this decides what that implies about anybody living on it.
  Nothing here feeds back -- no town moves a river.

  Three passes, and the pleasing thing is that each is the same kind of
  traversal the physical map was built from:

      towns        score every cell for how good a place it is, then take
                   the best, refusing any that is too close to one
                   already taken
      roads        a minimum spanning tree over the towns, with each link
                   routed by A* across the cells -- so a road bends round
                   a mountain rather than over it
      territories  multi-source Dijkstra outward from every town at once.
                   A cell belongs to whichever town is cheapest to reach
                   from, which is a border drawn by travel time

  ## Borders come from the cost, not from the distance

  The last one is the reason to do this on a graph. Territory assigned by
  straight-line distance gives you a Voronoi diagram of the towns and
  nothing else -- a partition that knows nothing about the land. Assigned
  by travel cost, a ridge pushes the border away from itself, because
  crossing it is expensive from both sides; a valley pulls the border
  along it, because the valley is cheap. Frontiers end up on watersheds
  and mountain chains without anything ever looking for one, which is
  where real ones tend to be.

  ## What is not here

  Names come from `allgo.procedural.naming`, one invented language per
  *culture* -- and cultures are deliberately coarser than realms. A realm
  here has exactly one town in it, so a language per realm would mean no
  two places on the map ever share one, which is the bag of unrelated
  syllables that namespace exists to avoid. Grouped into a handful of
  cultures instead, neighboring realms sound alike and the places where
  that changes are audible.

  No population, no economy, no politics, no history. Azgaar's generator
  has all of those and they are what make it a *fantasy* map generator;
  this stops at the geography people impose on land."
  (:require [allgo.graph :as graph]
            [allgo.procedural.naming :as naming]
            [allgo.search :as search]
            [clojure.math :as math]))

(def defaults
  {:towns 14
   :seed 1
   ;; Cells of separation, measured in graph hops rather than distance,
   ;; so towns spread out by how connected the land is rather than by how
   ;; the map happens to be scaled.
   :spacing 3})

;; ---------------------------------------------------------------------------
;; Where the land is worth living on

(def ^:private biome-appeal
  "Roughly, how much of a living the land offers. Not a yield model --
  an ordering, so that towns go to the meadows before the tundra."
  {:grassland 1.0 :temperate-deciduous-forest 0.9 :tropical-seasonal-forest 0.8
   :shrubland 0.6 :temperate-rain-forest 0.6 :tropical-rain-forest 0.5
   :taiga 0.35 :beach 0.5 :temperate-desert 0.15 :subtropical-desert 0.1
   :tundra 0.1 :bare 0.0 :scorched 0.0 :snow 0.0 :ice 0.0 :marsh 0.3})

(defn- river-cells
  "The cells with a river along one of their borders.

  Rivers run on the boundaries between cells, so being *on* a river is a
  fact about a cell's edges rather than about the cell."
  [{:keys [edges]}]
  (into #{}
        (comp (filter #(pos? (long (:river % 0))))
              (mapcat :centers)
              (remove nil?))
        edges))

(defn appeal
  "How good a cell is to settle, as a number with no units.

  Fresh water first, because everything else is negotiable and that is
  not: a river or a lakeshore outweighs every other term here. Then flat
  and low, then a coast for a harbor, then whether the land grows
  anything."
  [island]
  (let [on-river (river-cells island)
        {:keys [centers]} island]
    (fn [c]
      (if (:water? c)
        ##-Inf
        (+ (if (contains? on-river (:id c)) 2.2 0.0)
           (* 1.6 (double (:moisture c)))
           (* -1.4 (double (:elevation c)))
           (if (:coast? c) 1.0 0.0)
           (get biome-appeal (:biome c) 0.0)
           ;; A cell hemmed in by water is a rock, not a site.
           (let [land-nbrs (count (remove #(:water? (centers %)) (:neighbors c)))]
             (if (< land-nbrs 2) -3.0 0.0)))))))

(defn- within
  "Cell ids reachable from `start` in at most `hops` steps over land."
  [{:keys [centers]} start hops]
  (loop [frontier #{start} seen #{start} n 0]
    (if (or (= n (long hops)) (empty? frontier))
      seen
      (let [next-frontier (into #{}
                                (comp (mapcat #(:neighbors (centers %)))
                                      (remove seen)
                                      (remove #(:water? (centers %))))
                                frontier)]
        (recur next-frontier (into seen next-frontier) (inc n))))))

(defn assign-towns
  "Picks the sites, best first, refusing any within `:spacing` hops of one
  already taken.

  Greedy rather than optimal, and the refusal is what matters: without it
  every town lands in the same river valley, because that valley really
  is the best place and the second best place in it is still better than
  anywhere else."
  [island {:keys [towns spacing]}]
  (let [{:keys [centers]} island
        score (appeal island)
        ranked (->> centers
                    (remove :water?)
                    (map (fn [c] [(:id c) (score c)]))
                    (remove #(= ##-Inf (second %)))
                    (sort-by second >))
        chosen (loop [[[id _] & more] ranked taken [] blocked #{}]
                 (cond
                   (or (nil? id) (= (count taken) (long towns))) taken
                   (contains? blocked id) (recur more taken blocked)
                   :else (recur more (conj taken id)
                                (into blocked (within island id spacing)))))]
    (assoc island
           :towns (vec chosen)
           :centers (reduce (fn [cs id] (assoc-in cs [id :town?] true))
                            centers chosen))))

;; ---------------------------------------------------------------------------
;; Getting about

(defn- distance [[ax ay] [bx by]]
  (let [dx (- (double bx) ax) dy (- (double by) ay)]
    (math/sqrt (+ (* dx dx) (* dy dy)))))

(defn travel
  "`(fn [from to])`, the cost of moving between two neighboring cells.

  Distance, multiplied up by how much climbing it involves and by how
  unfriendly the destination is. The climb term is steep on purpose: it
  is what bends a road along a valley instead of over the ridge beside
  it, and with a gentle one the roads come out straight and the
  territories come out as a Voronoi diagram of the towns."
  [{:keys [centers]}]
  (fn [from to]
    (let [a (centers from)
          b (centers to)
          climb (abs (- (double (:elevation b)) (double (:elevation a))))]
      (* (distance (:point a) (:point b))
         (+ 1.0
            (* 14.0 climb)
            (- 1.0 (get biome-appeal (:biome b) 0.5)))))))

(defn- land-neighbors [{:keys [centers]}]
  (fn [id]
    (into [] (remove #(:water? (centers %))) (:neighbors (centers id)))))

(defn nearest-source
  "Multi-source Dijkstra: for every land cell, which of `sources` is
  cheapest to reach it from, and at what cost.

  One sweep, not one per source. Seed the queue with all of them and the
  first to arrive claims the cell, which is exactly the cheapest-source
  rule at the cost of a single search. Returns `[owner cost]`, both maps
  keyed by cell."
  [island sources]
  (let [step (travel island)
        nbrs (land-neighbors island)
        [best owner]
        (loop [;; Ordered by cost then id, so ties break the same way on
               ;; every platform rather than on map iteration order.
               queue (into (sorted-set) (map (fn [t] [0.0 t t])) sources)
               best (into {} (map (fn [t] [t 0.0])) sources)
               owner (into {} (map (fn [t] [t t])) sources)]
          (if-let [[cost id town :as entry] (first queue)]
            (let [queue (disj queue entry)]
              (if (> cost (get best id ##Inf))
                (recur queue best owner)
                (let [[queue best owner]
                      (reduce (fn [[q b o] nb]
                                (let [c (+ (double cost) (double (step id nb)))]
                                  (if (< c (get b nb ##Inf))
                                    [(conj q [c nb town]) (assoc b nb c) (assoc o nb town)]
                                    [q b o])))
                              [queue best owner]
                              (nbrs id))]
                  (recur queue best owner))))
            [best owner]))]
    [owner best]))

(defn assign-territories
  "Who holds what: every land cell goes to the town it is cheapest to
  reach.

  `:territory` is the town id and `:reach` the cost of getting there,
  which is what a border is drawn from."
  [{:keys [centers towns] :as island}]
  (let [[owner best] (nearest-source island towns)]
    (assoc island :centers
           (mapv (fn [c]
                   (assoc c
                          :territory (get owner (:id c))
                          :reach (get best (:id c))))
                 centers))))

(defn- spread-picks
  "`n` of `ids`, chosen to be as far apart as they can be.

  Farthest-point sampling: take one, then repeatedly take whichever is
  furthest from everything taken so far. Taking the first `n` instead
  would cluster them, because the list is ordered by how good a site is
  and the good sites are not evenly distributed."
  [island ids n]
  (let [pts (fn [id] (:point ((:centers island) id)))
        ids (vec ids)]
    (if (<= (count ids) (long n))
      ids
      (loop [taken [(first ids)]]
        (if (= (count taken) (long n))
          taken
          (recur (conj taken
                       (apply max-key
                              (fn [id] (reduce min (map #(distance (pts id) (pts %)) taken)))
                              (remove (set taken) ids)))))))))

(defn assign-cultures
  "Groups the realms into a handful of cultures, spatially.

  Without this each realm is its own culture, and since a realm has one
  town in it, no two places on the map ever share a language -- which
  makes the naming a bag of unrelated syllables again, the exact thing
  `allgo.procedural.naming` exists to avoid. Cultures are coarser than
  realms on purpose: several neighboring realms sound alike, and the
  places where that changes are the interesting borders."
  [{:keys [centers towns] :as island} {:keys [cultures]}]
  (let [n (max 1 (long (or cultures (max 2 (quot (count towns) 3)))))
        seeds (spread-picks island towns n)
        [owner _] (nearest-source island seeds)]
    (assoc island
           :cultures (vec seeds)
           :centers (mapv (fn [c] (assoc c :culture (get owner (:id c)))) centers))))

(defn- centers->edge
  "Which mesh edge joins each pair of neighboring cells."
  [{:keys [edges]}]
  (reduce (fn [m {:keys [id centers]}]
            (let [[a b] centers]
              (if (and a b) (assoc m #{a b} id) m)))
          {}
          edges))

(defn assign-roads
  "A spanning tree over the towns, each link walked by A*.

  The tree decides *which* towns are joined -- the cheapest set of links
  that still connects everything, so no two roads run parallel down the
  same valley. A* decides where each one goes, and because it pays the
  same travel cost the territories do, roads and borders agree about
  what the difficult ground is.

  Towns on different landmasses simply are not joined. Kruskal gives a
  spanning forest rather than failing, and a route A* cannot find is
  dropped -- there is no road from one island to another, which is
  correct, and a bridge is not something this models."
  [{:keys [towns] :as island}]
  (let [step (travel island)
        nbrs (land-neighbors island)
        by-pair (centers->edge island)
        pts (fn [id] (:point ((:centers island) id)))
        ;; Straight-line distance to choose the links; real travel cost to
        ;; route them. Measuring the tree properly would mean a search per
        ;; pair before knowing which pairs matter.
        pairs (for [a towns b towns :when (< (long a) (long b))] #{a b})
        tree (graph/minimum-spanning-tree
              (set towns) (set pairs)
              (fn [e] (let [[a b] (vec e)] (distance (pts a) (pts b)))))
        traffic (reduce
                 (fn [acc e]
                   (let [[a b] (vec e)
                         route (search/a-star
                                {:start a
                                 :goal b
                                 :neighbors nbrs
                                 :cost step
                                 ;; Straight-line distance never exceeds
                                 ;; the cost, whose multiplier is at least
                                 ;; one -- so A* stays optimal.
                                 :heuristic (fn [x g] (distance (pts x) (pts g)))})]
                     (if (nil? route)
                       acc
                       (reduce (fn [acc [u v]]
                                 (if-let [eid (by-pair #{u v})]
                                   (update acc eid (fnil inc 0))
                                   acc))
                               acc
                               (partition 2 1 route)))))
                 {}
                 tree)]
    (assoc island
           :edges (mapv (fn [e] (assoc e :road (get traffic (:id e) 0))) (:edges island))
           :roads (vec tree))))

;; ---------------------------------------------------------------------------
;; Names

(defn assign-names
  "One invented language per territory, and a name for everything that
  has one.

  A town is named out of its own culture's language, a realm out of its
  capital's, and a river out of whichever culture holds the land at its
  mouth. The rivers themselves come from `allgo.procedural.island`, which
  traced them; this only names them. Names are keyed by the id of the thing named rather than drawn
  from a running stream, so adding a town does not rename the rest of the
  map."
  [{:keys [centers towns cultures] :as island} {:keys [seed]}]
  (let [seed (long (or seed 1))
        langs (into {} (map (fn [c] [c (naming/language (+ (* 7919 seed) (long c)))]))
                    cultures)
        fallback (naming/language (+ (* 7919 seed) 104729))
        ;; A place speaks the language of its culture, not of its realm.
        lang-at (fn [id] (get langs (:culture (centers id)) fallback))
        realms (into {} (map (fn [t] [t (naming/name-for (lang-at t) (+ 5000 (long t)))])) towns)
        centers (mapv (fn [c]
                        (cond-> c
                          (:town? c) (assoc :name (naming/name-for (lang-at (:id c)) (:id c)))
                          (:territory c) (assoc :realm (get realms (:territory c)))))
                      centers)
        ;; The rivers already have bodies; this only gives them names.
        ;; A river is named out of the culture holding the land at its
        ;; mouth, which is where anyone naming it would have been
        ;; standing.
        rivers (mapv (fn [{:keys [mouth] :as r}]
                       (let [owner (some #(when (:culture (centers %)) %)
                                         (:touches ((:corners island) mouth)))]
                         (assoc r :name (naming/name-for (if owner (lang-at owner) fallback)
                                                         (+ 20011 (long mouth))))))
                     (:rivers island))]
    (assoc island
           :centers centers
           :realms realms
           :rivers rivers
           :languages langs)))

(defn populate
  "Towns, the roads between them, who holds what, and what it is all
  called."
  ([island] (populate island {}))
  ([island opts]
   (let [opts (merge defaults opts)]
     (-> island
         (assign-towns opts)
         assign-roads
         assign-territories
         (assign-cultures opts)
         (assign-names opts)))))
