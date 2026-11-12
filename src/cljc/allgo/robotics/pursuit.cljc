(ns allgo.robotics.pursuit
  "Visibility-based pursuit-evasion: LaValle, *Planning Algorithms*,
  Cambridge University Press, 2006, section 12.4 -- hide-and-seek, played
  so the seeker can prove that no one is left hiding.

  A pursuer moves through a region bounded by a simple polygon, seeing
  everything its straight lines of sight reach; an evader, as fast as it
  likes, may be anywhere. The pursuer's knowledge -- the book's
  nondeterministic information state -- is where it is and which parts of
  the shadow, the region it cannot see, might hide the evader. Only when
  the shadow's components change does that knowledge change, and in four
  ways (section 12.4.2):

    appear      a new component, just seen, so cleared
    disappear   a component comes into view and is gone
    split       a component parts in two, each keeping its label
    merge       two join, contaminated if either was -- the one event
                that undoes the pursuer's work

  The book's complete algorithm decomposes the region into cells where the
  components cannot change -- bounded by inflection and bitangent rays --
  and searches a directed graph whose vertices are a cell and a labeling
  of its components, from every component contaminated to none. Here the
  cells are replaced by a lattice of probe points: the shadow from a point
  is the probes it cannot see, its components those joined through
  neighboring probes, and a pursuer step from one probe to the next
  carries the labels across by overlap. A component is contaminated if it
  holds or touches a probe that was contaminated -- the lattice cannot
  tell a thin sliver of shadow from a gap, so contamination is allowed to
  spread to a neighbor, and a plan never claims a corner clear that it is
  not. The graph is searched breadth first, revealed as the search goes,
  as the book suggests; with it a plan has the fewest steps.

  An environment is built from a polygon and a spacing by `environment`;
  points are `[x y]`, probes their indices into its `:probes`. An
  information state is `{:at :shadow :contaminated}`: the pursuer's point,
  the shadow's components as sorted vectors of probes, and the set of
  indices of the components that might hide the evader."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- blocks?
  "Whether the edge `c`-`d` cuts the sight line `a`-`b`: crosses it
  between its ends, a corner touching it counting."
  [[ax ay] [bx by] [[cx cy] [dx dy]]]
  ;; polygon/segment-intersection, unrolled: this is the inner loop
  (let [ax (double ax) ay (double ay) cx (double cx) cy (double cy)
        rx (- (double bx) ax) ry (- (double by) ay)
        sx (- (double dx) cx) sy (- (double dy) cy)
        den (- (* rx sy) (* ry sx))]
    (when-not (< (abs den) 1e-12)
      (let [qx (- cx ax) qy (- cy ay)
            t (/ (- (* qx sy) (* qy sx)) den)
            u (/ (- (* qx ry) (* qy rx)) den)]
        (and (< 1e-9 t (- 1.0 1e-9)) (<= -1e-9 u (+ 1.0 1e-9)))))))

(defn sees?
  "Whether points `a` and `b`, inside the environment, see each other: no
  edge of its polygon crosses the segment between them."
  [{:keys [edges]} a b]
  (not-any? #(blocks? a b %) edges))

(def ^:private four [[1 0] [-1 0] [0 1] [0 -1]])

(defn environment
  "The region inside the simple polygon `poly`, with probes on a lattice
  `spacing` apart, offset half a spacing from the polygon's lowest
  corner so none falls on an edge through two lattice-aligned corners.

    :polygon :edges     the boundary
    :probes             the points inside, in order
    :cells              each probe's lattice `[i j]`
    :neighbors          each probe's lattice neighbors that it sees, four
                        ways -- the steps a pursuer takes, and the joins
                        between probes of one shadow component"
  [poly spacing]
  (let [x0 (reduce min (map first poly)) y0 (reduce min (map second poly))
        ni (long (math/ceil (/ (- (reduce max (map first poly)) x0) spacing)))
        nj (long (math/ceil (/ (- (reduce max (map second poly)) y0) spacing)))
        inside (vec (for [i (range ni) j (range nj)
                          :let [p [(+ x0 (* (+ i 0.5) spacing)) (+ y0 (* (+ j 0.5) spacing))]]
                          :when (polygon/contains-point? poly p)]
                      [[i j] p]))
        index (zipmap (map first inside) (range))
        env {:polygon poly :spacing spacing :edges (polygon/edges poly)
             :probes (mapv second inside) :cells (mapv first inside)}]
    (assoc env :neighbors
           (mapv (fn [[[i j] p]]
                   (vec (for [[di dj] four
                              :let [k (index [(+ i di) (+ j dj)])]
                              :when (and k (sees? env p (second (inside k))))]
                          k)))
                 inside))))

(defn nearest-probe
  "The probe nearest point `p`."
  [{:keys [probes]} p]
  (apply min-key #(polygon/distance p (probes %)) (range (count probes))))

(defn visible
  "The set of probes seen from point `p`."
  [{:keys [probes] :as env} p]
  (set (filter #(sees? env p (probes %)) (range (count probes)))))

(defn visibility-polygon
  "The region seen from point `p`, exactly, as a polygon: a ray toward
  each corner of the boundary and just either side of it, each stopped at
  the first edge it meets, in order of angle."
  [{:keys [polygon edges]} p]
  (let [reach (* 2.0 (reduce max (map #(polygon/distance p %) polygon)))
        hit (fn [a]
              (let [q (polygon/add p [(* reach (math/cos a)) (* reach (math/sin a))])
                    t (reduce min 1.0 (keep (fn [[c d]]
                                              (when-let [[t u] (polygon/segment-intersection p q c d)]
                                                (when (and (> t 1e-9) (<= -1e-9 u (+ 1.0 1e-9))) t)))
                                            edges))]
                (polygon/lerp p q t)))]
    (mapv hit (sort (for [v polygon :when (not= v p) d [-1e-6 0.0 1e-6]]
                      (+ (polygon/angle-of (polygon/sub v p)) d))))))

(defn shadow
  "The shadow from point `p`, the region it cannot see, as its connected
  components: each a sorted vector of probes, ordered by their first."
  [{:keys [probes neighbors] :as env} p]
  (let [seen (visible env p)
        hidden-neighbors (fn [k] (remove seen (neighbors k)))
        fill (fn [k] (-> (search/breadth-first-seq {:start k :goal? (constantly false) :neighbors hidden-neighbors})
                         last :came-from keys sort vec))]
    ;; each probe not yet in a component starts the next
    (second (reduce (fn [[found comps :as acc] k]
                      (if (found k) acc (let [c (fill k)] [(into found c) (conj comps c)])))
                    [#{} []]
                    (remove seen (range (count probes)))))))

(defn contaminated-probes
  "The probes that might hide the evader in information state `s`."
  [{:keys [shadow contaminated]}]
  (into #{} (mapcat shadow) contaminated))

(defn relabel
  "Which of the components `to` are contaminated, coming from components
  `from` with those `contaminated`: those holding or touching a probe of
  one. It is the book's rules at once -- a component that appears holds
  none and is cleared, a split's parts each hold some of a contaminated
  one, a merge holds some of each."
  [{:keys [neighbors]} from contaminated to]
  (let [dirty (into #{} (mapcat from) contaminated)
        touched (into dirty (mapcat neighbors) dirty)]
    (set (keep-indexed (fn [i c] (when (some touched c) i)) to))))

(defn start
  "The information state at point `p` before anything is known: every
  component of its shadow contaminated."
  [env p]
  (let [s (shadow env p)]
    {:at p :shadow s :contaminated (set (range (count s)))}))

(defn clearing-step
  "The information state after the pursuer moves from where `state` has
  it to point `p`, a step short enough that no component jumps a probe."
  [env {from :shadow dirty :contaminated} p]
  (let [to (shadow env p)]
    {:at p :shadow to :contaminated (relabel env from dirty to)}))

(defn clearing-seq
  "The information states as the pursuer follows the points of `path`
  from the first, every component contaminated there."
  [env path]
  (reductions (partial clearing-step env) (start env (first path)) (rest path)))

(defn cleared?
  "Whether information state `s` proves the region empty of evaders."
  [s]
  (empty? (:contaminated s)))

(defn events
  "The critical events between components `from` and `to`, the shadow
  before and after a step: `[:appear j]`, `[:disappear i]`, `[:split i
  js]` and `[:merge is j]`, by which components overlap."
  [from to]
  (let [owner (into {} (for [[j c] (map-indexed vector to) k c] [k j]))
        links (set (for [[i c] (map-indexed vector from) k c :let [j (owner k)] :when j] [i j]))
        succ (group-by first links) pred (group-by second links)]
    (concat (for [j (range (count to)) :when (not (pred j))] [:appear j])
            (for [i (range (count from)) :when (not (succ i))] [:disappear i])
            (for [[i ls] succ :when (> (count ls) 1)] [:split i (sort (map second ls))])
            (for [[j ls] pred :when (> (count ls) 1)] [:merge (sort (map first ls)) j]))))

;; ---------------------------------------------------------------------------
;; The complete algorithm: search the information space

(defn- graph-neighbors
  "The vertices one pursuer step from `[k dirty]`, the pursuer at probe
  `k` with the components `dirty` contaminated: each neighbor probe, with
  the components there that the step leaves contaminated."
  [{:keys [neighbors] :as env} shadow-of [k dirty]]
  (for [n (neighbors k)]
    [n (relabel env (shadow-of k) dirty (shadow-of n))]))

(defn search-seq
  "The search of the information space from point `p`, breadth first, a
  vertex `[probe contaminated]` expanded a step, until one whose shadow is
  all cleared: `allgo.search`'s search states. Each probe's shadow is
  computed once, the first time the search reaches it."
  [{:keys [probes] :as env} p]
  (let [shadow-of (memoize (fn [k] (shadow env (probes k))))
        k (nearest-probe env p)]
    (search/breadth-first-seq {:start [k (set (range (count (shadow-of k))))]
                               :goal? (comp empty? second)
                               :neighbors (partial graph-neighbors env shadow-of)})))

(defn plan
  "The pursuer's path from point `p` -- the probes it visits, as points --
  that leaves no place for an evader to hide, or nil when there is none:
  the search has run out of information states without reaching one
  with the shadow all cleared."
  [{:keys [probes] :as env} p]
  (some->> (search-seq env p) search/goal-state search/path (mapv (comp probes first))))
