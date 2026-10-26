(ns allgo.procedural.city.lots
  "Blocks between the roads, and the lots they are divided into:
  Parish and Mueller's section 4.1.

  A block is a face of the road network, which is a plane graph -- every
  crossing was made by cutting both roads at one point -- so the faces
  are found by walking it: leave along an edge, and at every junction
  take the sharpest left, until back where you began. Each walk goes
  once around a face, counterclockwise for a block and clockwise for the
  outside of the whole network. Dead ends are pruned first; they bound
  nothing.

  The block is then shrunk by half the width of each road around it,
  which is the paper's \"blocks created by scaling from street
  crossings\", and cut into lots by their simple recursive rule: split
  across the longest side until every piece is under a threshold area.
  Each side of a lot remembers whether it faces a street, and then, as
  in the paper, lots that are too small, that are not convex, or that
  have no street to face are discarded -- the middle of a deep block is
  courtyard, not building land. The threshold area follows the zone, so
  downtown lots are large, for towers, and houses get narrow plots."
  (:require [allgo.geometry.polygon :as poly]
            [allgo.random :as random]))

(def defaults
  {:seed 1
   ;; Half widths of each kind of road, plus a pavement, meters.
   :half-width {:highway 11.0 :street 6.5}
   ;; Largest lot, by zone, square meters.
   :lot-area {:downtown 2600.0 :commercial 1300.0 :residential 520.0}
   ;; Blocks bigger than this are not blocks but country.
   :max-block-area 90000.0})

;; ---------------------------------------------------------------------------
;; Faces of the road graph

(defn- prune-dead-ends
  "The graph's edges once every road end has been cut back to the last
  junction, repeatedly: an adjacency map of node to neighbor to kind."
  [edges]
  (let [adj (reduce (fn [m [i j kind]]
                      (-> m (assoc-in [i j] kind) (assoc-in [j i] kind)))
                    {} edges)]
    (loop [adj adj]
      (let [leaves (keep (fn [[n nbrs]] (when (<= (count nbrs) 1) n)) adj)]
        (if (empty? leaves)
          adj
          (recur (reduce (fn [m n]
                           (let [m (reduce (fn [m o] (update m o dissoc n)) m (keys (get m n)))]
                             (dissoc m n)))
                         adj leaves)))))))

(defn faces
  "The faces of the plane graph `{:points :edges}`, each as `{:points
  :kinds}`: the corners counterclockwise, and the kind of road along
  each side. The unbounded outside is left out."
  [{:keys [points edges]}]
  (let [adj (prune-dead-ends edges)
        angle (fn [i j] (poly/angle-of (poly/sub (nth points j) (nth points i))))
        ;; Neighbors in counterclockwise order around each node.
        order (into {} (for [[i nbrs] adj]
                         [i (vec (sort-by #(angle i %) (keys nbrs)))]))
        next-edge (fn [u v]
                    ;; At v, coming from u: the neighbor just clockwise of
                    ;; u, which is the sharpest left.
                    (let [ring (order v)
                          k (first (keep-indexed #(when (= %2 u) %1) ring))]
                      [v (nth ring (mod (dec (+ k (count ring))) (count ring)))]))
        half-edges (for [[i nbrs] adj j (keys nbrs)] [i j])]
    (loop [todo half-edges seen #{} out []]
      (if-let [e (first todo)]
        (if (seen e)
          (recur (rest todo) seen out)
          (let [walk (loop [e e acc []]
                       (if (and (seq acc) (= e (first acc)))
                         acc
                         (if (> (count acc) 2000)
                           acc
                           (recur (apply next-edge e) (conj acc e)))))
                pts (mapv #(nth points (first %)) walk)
                kinds (mapv (fn [[a b]] (get-in adj [a b])) walk)]
            (recur (rest todo) (into seen walk)
                   (if (pos? (poly/signed-area pts))
                     (conj out {:points pts :kinds kinds})
                     out))))
        out))))

;; ---------------------------------------------------------------------------
;; Blocks

(defn blocks
  "The building land of each face: shrunk away from the roads around it
  by their half widths, and tidied. Faces too large to be blocks, and
  those that shrink to nothing, are dropped."
  ([graph] (blocks graph {}))
  ([graph opts]
   (let [{:keys [half-width max-block-area]} (merge defaults opts)]
     (into []
           (keep (fn [{:keys [points kinds]}]
                   (when (< (poly/area points) (double max-block-area))
                     (let [[pts kinds] (poly/simplify points kinds 2.0)]
                       (when (>= (count pts) 3)
                         (when-let [inner (poly/inset pts (mapv #(get half-width % 6.0) kinds))]
                           (when (> (poly/area inner) 60.0)
                             {:points inner :kinds kinds})))))))
           (faces graph)))))

;; ---------------------------------------------------------------------------
;; Lots

(defn- longest-edge [pts]
  (apply max-key #(let [[a b] %] (poly/distance a b)) (poly/edges pts)))

(defn subdivide
  "Cuts `pts` (with per-side `tags`) across its longest side, again and
  again, until every piece is under `max-area` -- `(fn [piece])` for an
  area that depends on where the piece is. The cut falls near the middle
  of that side, `jitter` of its length either way. Returns the pieces,
  `{:points :tags}`."
  [pts tags max-area rng {:keys [jitter] :or {jitter 0.12}}]
  (let [area-of (if (number? max-area) (constantly max-area) max-area)]
    (loop [todo [[pts tags 0]] out []]
      (if-let [[p t depth] (peek todo)]
        (let [todo (pop todo)
              a (poly/area p)]
          (if (or (<= a (double (area-of p))) (> depth 24))
            (recur todo (conj out {:points p :tags t}))
            (let [[e0 e1] (longest-edge p)
                  at (poly/lerp e0 e1 (+ 0.5 (random/uniform rng (- jitter) jitter)))
                  along (poly/unit (poly/sub e1 e0))
                  ;; The cut runs across that side.
                  d [(- (double (second along))) (double (first along))]
                  [lp lt] (poly/clip p t at d nil)
                  [rp rt] (poly/clip p t at (poly/scale d -1.0) nil)]
              (recur (cond-> todo
                       (seq lp) (conj [lp lt (inc depth)])
                       (seq rp) (conj [rp rt (inc depth)]))
                     out))))
        out))))

(defn- street-sides
  "Indices of the sides of a lot that face a street, at least `min-len`
  long."
  [{:keys [points tags]} min-len]
  (filterv (fn [i]
             (let [[a b] (nth (poly/edges points) i)]
               (and (nth tags i) (>= (poly/distance a b) min-len))))
           (range (count points))))

(defn lots
  "The lots of every block, `{:points :street :zone :area}`: `:street`
  the indices of the sides that face a road. `zone` is a function of
  `[x y]` giving the zone a lot is in, which sets how large it may be."
  ([blocks zone] (lots blocks zone {}))
  ([blocks zone opts]
   (let [{:keys [seed lot-area]} (merge defaults opts)
         rng (random/rng seed 29 0 0)
         size-at (fn [p] (let [[x y] (poly/centroid p)] (get lot-area (zone x y) 600.0)))]
     (into []
           (comp
            (mapcat (fn [{:keys [points kinds]}]
                      (subdivide points kinds size-at rng {})))
            (keep (fn [{:keys [points tags]}]
                    (let [[pts tgs] (poly/simplify points tags 0.5)
                          piece {:points pts :tags tgs}
                          a (poly/area pts)
                          [x y] (poly/centroid pts)
                          z (zone x y)
                          sides (street-sides piece 5.0)]
                      (when (and (>= (count pts) 3)
                                 (poly/convex? pts)
                                 (> a (* 0.22 (double (get lot-area z 600.0))))
                                 (seq sides))
                        {:points pts :street sides :zone z :area a})))))
           blocks))))
