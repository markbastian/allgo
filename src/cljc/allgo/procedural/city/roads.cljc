(ns allgo.procedural.city.roads
  "Highways and streets, grown by Parish and Mueller's extended
  L-system.

  ## The L-system, and what is moved out of it

  Their road L-system has only a handful of productions. A road module
  `R` proposes the next segment, an *ideal successor*, with its
  parameters left unassigned; branch modules `B` wait out a delay and
  then become roads of their own. Everything that makes a road go one
  way rather than another is moved out of the rules into two functions
  the rules call:

      global goals       what the road would like: highways head for
                         population, streets fill in where people live,
                         and both follow the street pattern of their
                         district
      local constraints  what the world allows: no building in the
                         water, and a road that runs into, or nearly
                         into, another one joins it

  so a new pattern or constraint is a change to one function and none to
  the rules. Their delays order the growth -- a branch waits `del` steps
  before it starts -- and that ordering is all the productions decide,
  so here they are run as what they amount to: a queue of proposed
  segments ordered by when each is due. Popping a proposal and checking
  it against the local constraints is the paper's insertion query `?I`;
  if it succeeds the segment is laid, and the global goals propose its
  successors, each due some steps later.

  ## Global goals

  Highways connect centers of population. From the end of a highway,
  rays are cast in a fan ahead of it, the population along each summed
  with weights falling off with distance, and the road turns toward the
  best. Streets do not chase population; they fill in the street
  pattern, and stop where nobody lives. The pattern rules narrow the
  angles and lengths a road may take:

      :grid     New York: along the local grid's two axes, a block apart
      :radial   Paris: toward or around a center
      :organic  the paper's basic rule: no imposed pattern
      :terrain  San Francisco: highways take the route of least climb,
                streets go straight up the slope or along it, and short

  and where the land's pattern map gives more than one a weight, the
  directions each proposes are summed with those weights.

  ## Local constraints, and self-sensitivity

  A segment that ends in the water is pruned back to the shore, or
  rotated until it fits, which is what makes roads run along a coast;
  a highway may instead cross a short stretch as a bridge. Then the
  segment is checked against the roads already there, the paper's
  self-sensitive L-system: if it crosses one, it stops at the crossing;
  if it ends near an existing junction, it is pulled to it; and if it
  ends just short of a road, it is extended to meet it. Each of these
  makes a junction, and a road that has joined another stops growing.
  This is what turns a branching tree into a network of loops.

  Not done here: the paper's smoothing of sharp bends by subdivision,
  and its streets lowering the population around them (from Open
  L-systems); the junction rules already keep streets from piling up."
  (:require [allgo.geometry.polygon :as poly]
            [allgo.procedural.city.land :as land]
            [allgo.random :as random]
            [clojure.math :as math]))

(def defaults
  {:seed 1
   :highway-length 170.0
   ;; A block, more or less.
   :street-length 95.0
   ;; Chance a highway segment sends a highway off to one side.
   :highway-branch 0.1
   ;; Chance, per side, that a highway segment sends a street off.
   :highway-streets 0.8
   ;; Chance, per side, that a street segment branches.
   :street-branch 0.75
   ;; Steps a street waits before it starts: highways first.
   :street-delay 6
   ;; Below this, streets stop growing.
   :min-population 0.1
   ;; How far a highway looks ahead for people, and how wide.
   :ray-length 700.0
   :fan-angle 0.35
   ;; Junction rules, as fractions of the segment's own length.
   :snap 0.45
   :extend 0.55
   ;; No two roads leave a junction closer than this, radians.
   :min-angle 0.52
   :max-bridge 320.0
   :max-segments 5000})

;; ---------------------------------------------------------------------------
;; The network so far: segments, and a grid of which segments touch which
;; cells, so a new segment is checked against its neighbors only.

(def ^:private cell 80.0)

(defn- cells-of
  "The grid cells under the bounding box of `a`-`b`, a little padded."
  [[ax ay] [bx by] pad]
  (let [pad (double pad)
        i0 (long (math/floor (/ (- (min ax bx) pad) cell)))
        i1 (long (math/floor (/ (+ (max ax bx) pad) cell)))
        j0 (long (math/floor (/ (- (min ay by) pad) cell)))
        j1 (long (math/floor (/ (+ (max ay by) pad) cell)))]
    (for [i (range i0 (inc i1)) j (range j0 (inc j1))] [i j])))

(defn- index-segment [index id a b]
  (reduce (fn [idx c] (update idx c (fnil conj #{}) id)) index (cells-of a b 0.0)))

(defn- unindex-segment [index id a b]
  (reduce (fn [idx c] (update idx c disj id)) index (cells-of a b 0.0)))

(defn- nearby
  "Ids of segments that might come within `pad` of `a`-`b`."
  [{:keys [index]} a b pad]
  (into #{} (mapcat #(get index %)) (cells-of a b pad)))

(defn- add-segment [state seg]
  (let [id (count (:segments state))
        {:keys [a b]} seg]
    (-> state
        (update :segments conj seg)
        (update :index index-segment id a b)
        (update-in [:nodes a] (fnil conj #{}) id)
        (update-in [:nodes b] (fnil conj #{}) id))))

(defn- split-segment
  "Cuts segment `id` at point `p`: it keeps the part from its start, and
  the rest becomes a new segment of the same kind."
  [state id p]
  (let [{:keys [a b] :as seg} (get-in state [:segments id])]
    (-> state
        (update :index unindex-segment id a b)
        (update-in [:nodes b] disj id)
        (assoc-in [:segments id] (assoc seg :b p))
        (update :index index-segment id a p)
        (update-in [:nodes p] (fnil conj #{}) id)
        (add-segment (assoc seg :a p :b b)))))

(defn- directions-at
  "Unit directions of the segments leaving node `p`."
  [state p]
  (for [id (get-in state [:nodes p])
        :let [{:keys [a b]} (get-in state [:segments id])
              q (if (= a p) b a)]]
    (poly/unit (poly/sub q p))))

(defn- crowded?
  "Whether a road leaving `p` along `d` would run too close beside one
  already there."
  [state p d min-angle]
  (let [c (math/cos min-angle)]
    (some #(> (poly/dot % d) c) (directions-at state p))))

;; ---------------------------------------------------------------------------
;; Global goals

(defn- perp [[x y]] [(- (double y)) (double x)])

(defn- closest-of
  "Of the unit vectors `cands`, the one most nearly along `d`."
  [cands d]
  (apply max-key #(poly/dot % d) cands))

(defn- grid-axes [angle]
  (let [u (poly/heading angle) v (perp u)]
    [u v (poly/scale u -1.0) (poly/scale v -1.0)]))

(defn- radial-axes [{:keys [center]} p]
  (let [r (poly/unit (poly/sub p center))
        r (if (zero? (poly/norm r)) [1.0 0.0] r)
        t (perp r)]
    [r t (poly/scale r -1.0) (poly/scale t -1.0)]))

(defn- fan
  "Directions spread over `half-width` radians either side of `d`."
  [d half-width n]
  (let [a (poly/angle-of d)]
    (for [i (range n)]
      (poly/heading (+ a (* half-width (- (/ (* 2.0 i) (dec n)) 1.0)))))))

(defn- most-populous
  "The paper's ray casting: of a fan of directions ahead of `p`, the one
  along which the population, weighted by the inverse distance, sums
  highest -- or straight on, once there is nobody ahead worth turning
  for, so that a highway leaving town heads out rather than circling
  back to the people behind it."
  [{:keys [population]} p d {:keys [ray-length fan-angle min-population]}]
  (let [samples 8
        ray-length (double ray-length)
        score (fn [u]
                (reduce + (for [k (range 1 (inc samples))
                                :let [r (* ray-length (/ k (double samples)))
                                      [x y] (poly/add p (poly/scale u r))]]
                            (/ (double (population x y)) r))))
        best (apply max-key score (fan d (double fan-angle) 9))
        ;; What the score would be with the threshold population all
        ;; along the ray.
        floor (* (double min-population) (/ samples ray-length)
                 (reduce + (map #(/ 1.0 %) (range 1 (inc samples)))))]
    (if (> (score best) floor) best d)))

(defn- least-climb
  "Of a fan of directions ahead of `p`, the one whose next `len` meters
  change elevation least."
  [{:keys [elevation]} p d len]
  (let [[px py] p e0 (double (elevation px py))]
    (apply min-key
           (fn [u]
             (let [[x y] (poly/add p (poly/scale u len))]
               (+ (abs (- (double (elevation x y)) e0))
                  ;; A tie goes to the straighter.
                  (* 0.5 (- 1.0 (poly/dot u d))))))
           (fan d 0.6 9))))

(defn- slope-axes
  "Up the slope, down it and along it either way, at `p`."
  [{:keys [elevation]} [px py]]
  (let [h 10.0
        gx (- (double (elevation (+ px h) py)) (double (elevation (- px h) py)))
        gy (- (double (elevation px (+ py h))) (double (elevation px (- py h))))
        g (poly/unit [gx gy])]
    (if (zero? (poly/norm g))
      [[1.0 0.0] [0.0 1.0] [-1.0 0.0] [0.0 -1.0]]
      (let [t (perp g)] [g t (poly/scale g -1.0) (poly/scale t -1.0)]))))

(defn- pattern-direction
  "The direction the pattern rule `rule` would have a road leaving `p`
  take, given that it would like to go along `d`."
  [land rule kind p d opts rng]
  (let [[px py] p]
    (case rule
      :grid (closest-of (grid-axes ((:grid-angle land) px py)) d)
      :radial (closest-of (radial-axes land p) d)
      :organic (if (= kind :highway)
                 d
                 (poly/heading (+ (poly/angle-of d) (random/gaussian rng 0.0 0.12))))
      :terrain (if (= kind :highway)
                 (least-climb land p d (:highway-length opts))
                 (closest-of (slope-axes land p) d)))))

(defn- blended-direction
  "The directions of every active pattern rule at `p`, summed with the
  weights of the land's pattern map."
  [land kind p d opts rng]
  (let [[px py] p
        weights ((:pattern land) px py)
        v (reduce (fn [acc [rule w]]
                    (if (pos? (double w))
                      (poly/add acc (poly/scale (pattern-direction land rule kind p d opts rng) w))
                      acc))
                  [0.0 0.0] weights)]
    (if (< (poly/norm v) 1e-6) d (poly/unit v))))

(defn- street-length [land p opts rng]
  (let [[px py] p
        {:keys [organic terrain] :or {organic 0.0 terrain 0.0}} ((:pattern land) px py)
        base (double (:street-length opts))]
    ;; A grid keeps to its block; the basic rule and the hills vary.
    (* base (+ 1.0 (* (+ (double organic) (* 0.5 (double terrain)))
                      (random/uniform rng -0.35 0.3))))))

(defn- global-goals
  "The successors of a segment just laid, from `a` to `b` along `d`:
  where it would like to continue, and what it would like to branch,
  each with how many steps it waits."
  [land {:keys [kind a b]} opts rng]
  (let [d (poly/unit (poly/sub b a))
        [bx by] b
        pop-here (double ((:population land) bx by))
        {:keys [highway-length highway-branch highway-streets street-branch
                street-delay min-population]} opts
        side (fn [s] (poly/scale (perp d) s))]
    (if (= kind :highway)
      (let [ahead (blended-direction land :highway b (most-populous land b d opts) opts rng)]
        (concat
         [{:kind :highway :a b :dir ahead :len highway-length :delay 1}]
         (for [s [1.0 -1.0]
               :when (< (double (rng)) (* (double highway-branch) (+ 0.5 pop-here)))]
           {:kind :highway :a b :dir (blended-direction land :highway b (side s) opts rng)
            :len highway-length :delay 1})
         (when (> pop-here (double min-population))
           (for [s [1.0 -1.0]
                 :when (< (double (rng)) (double highway-streets))]
             {:kind :street :a b :dir (blended-direction land :street b (side s) opts rng)
              :len (street-length land b opts rng) :delay street-delay}))))
      (when (> pop-here (double min-population))
        (concat
         [{:kind :street :a b :dir (blended-direction land :street b d opts rng)
           :len (street-length land b opts rng) :delay 1}]
         (for [s [1.0 -1.0]
               :when (< (double (rng)) (double street-branch))]
           {:kind :street :a b :dir (blended-direction land :street b (side s) opts rng)
            :len (street-length land b opts rng) :delay 1}))))))

;; ---------------------------------------------------------------------------
;; Local constraints

(defn- legal-run
  "Whether every point along `a`-`b`, every few meters, is legal."
  [land a b]
  (let [n (max 2 (long (math/ceil (/ (poly/distance a b) 8.0))))]
    (every? #(land/legal? land (poly/lerp a b (/ % (double n)))) (range 1 (inc n)))))

(defn- legalize
  "The paper's first local constraint: a segment that ends in, or
  crosses, illegal ground is pruned back to where it leaves the land,
  or turned until it fits -- or, for a highway, carried across as a
  bridge if land is near enough on the far side. nil if none works."
  [land {:keys [kind a dir len] :as prop} {:keys [max-bridge]}]
  (let [len (double len)
        b (poly/add a (poly/scale dir len))]
    (if (legal-run land a b)
      (assoc prop :b b)
      (let [step 8.0
            first-bad (first (filter #(not (land/legal? land (poly/add a (poly/scale dir %))))
                                     (range step (+ len step) step)))
            bridge (when (and (= kind :highway) first-bad)
                     (first (filter #(land/legal? land (poly/add a (poly/scale dir %)))
                                    (range first-bad (+ first-bad (double max-bridge)) step))))]
        (cond
          ;; Across: the bridge ends a little way up the far bank.
          (and bridge (legal-run land (poly/add a (poly/scale dir (+ bridge 4.0)))
                                 (poly/add a (poly/scale dir (+ bridge 30.0)))))
          (assoc prop :b (poly/add a (poly/scale dir (+ bridge 30.0))) :bridge? true)

          ;; Pruned, if enough is left to be worth a road.
          (and first-bad (> (- first-bad step) (* 0.55 len)))
          (assoc prop :b (poly/add a (poly/scale dir (- first-bad step))) :dead-end? true)

          ;; Turned, a little at a time, alternately each way.
          :else
          (some (fn [da]
                  (let [d2 (poly/heading (+ (poly/angle-of dir) da))
                        b2 (poly/add a (poly/scale d2 len))]
                    (when (legal-run land a b2) (assoc prop :dir d2 :b b2))))
                (for [k (range 1 7) s [1.0 -1.0]] (* s k 0.09))))))))

(defn- first-crossing
  "The existing segment `a`-`b` first runs into, as `[t id p]`, ignoring
  those that meet it at `a`."
  [state a b]
  (let [hits (for [id (nearby state a b 0.0)
                   :let [{c :a d :b} (get-in state [:segments id])]
                   :when (and (not= c a) (not= d a))
                   :let [[t u] (poly/segment-intersection a b c d)]
                   :when (and t (> t 1e-9) (<= t 1.0) (<= -1e-9 u (+ 1.0 1e-9)))]
               [t id (poly/lerp a b t) u])]
    (when (seq hits) (apply min-key first hits))))

(defn- nearest-node
  "The existing junction or road end nearest `p` within `r`, other than
  `a`."
  [state a p r]
  (let [cands (for [id (nearby state p p r)
                    :let [{c :a d :b} (get-in state [:segments id])]
                    q [c d]
                    :when (and (not= q a) (< (poly/distance q p) r))]
                q)]
    (when (seq cands) (apply min-key #(poly/distance % p) cands))))

(defn- join-at
  "Ends the proposal at an existing segment `id`, at fraction `u` of its
  length, point `p`: at one of its ends if that is close, otherwise by
  splitting it there. Returns `[state end]`, or nil if the crossing is
  too shallow to be a junction."
  [state id u p {:keys [min-angle]} dir]
  (let [{c :a d :b} (get-in state [:segments id])
        along (poly/unit (poly/sub d c))
        snap 6.0]
    (cond
      (< (* (double u) (poly/distance c d)) snap) [state c]
      (< (* (- 1.0 (double u)) (poly/distance c d)) snap) [state d]
      (< (abs (poly/cross along dir)) (math/sin (* 0.7 (double min-angle)))) nil
      :else [(split-segment state id p) p])))

(defn- self-sensitive
  "The paper's second local constraint: the segment against the roads
  already there. Returns `[state segment stopped?]` or nil."
  [land state {:keys [a b dir len bridge?] :as prop} {:keys [snap extend min-angle] :as opts}]
  (let [len (double len)
        finish (fn [[st end]]
                 ;; Pulled or stretched to a junction, the end has moved:
                 ;; the new stretch has to be dry as well.
                 (when (and end (> (poly/distance a end) (* 0.2 len))
                            (or bridge? (legal-run land a end))
                            (not (crowded? st end (poly/unit (poly/sub a end)) min-angle)))
                   [st (assoc prop :b end) true]))]
    (if-let [[_ id p u] (first-crossing state a b)]
      ;; It runs into a road: a crossing, and it stops there.
      (finish (join-at state id u p opts dir))
      (if-let [q (nearest-node state a b (* (double snap) len))]
        ;; It ends near a junction: pulled to it, if the way there is
        ;; clear.
        (when-not (first-crossing state a q)
          (finish [state q]))
        (let [reach (poly/add b (poly/scale dir (* (double extend) len)))]
          (if-let [[_ id p u] (first-crossing state a reach)]
            ;; It ends just short of a road: extended to meet it.
            (finish (join-at state id u p opts dir))
            [state prop false]))))))

(defn- local-constraints
  [land state prop opts]
  (when-let [prop (legalize land prop opts)]
    (when-not (crowded? state (:a prop) (:dir prop) (:min-angle opts))
      (self-sensitive land state prop opts))))

;; ---------------------------------------------------------------------------
;; Growth

(defn- enqueue [state prop t]
  (let [n (:proposals state)]
    (-> state
        (update :queue assoc [(long t) n] prop)
        (assoc :proposals (inc n)))))

(defn start
  "The axiom: two highways leaving the land's main center in opposite
  directions along its grid."
  ([land] (start land {}))
  ([land opts]
   (let [opts (merge defaults opts)
         c (:center land)
         [cx cy] c
         d (poly/heading ((:grid-angle land) cx cy))]
     (-> {:opts opts :segments [] :index {} :nodes {} :queue (sorted-map)
          :proposals 0 :steps 0 :time 0}
         (enqueue {:kind :highway :a c :dir d :len (:highway-length opts)} 0)
         (enqueue {:kind :highway :a c :dir (poly/scale d -1.0) :len (:highway-length opts)} 0)))))

(defn done? [{:keys [queue segments opts]}]
  (or (empty? queue) (>= (count segments) (long (:max-segments opts)))))

(defn step
  "One proposal out of the queue: checked, laid if it can be, and its
  successors proposed."
  [land state]
  (if (done? state)
    state
    (let [[[t n :as k] prop] (first (:queue state))
          {:keys [opts]} state
          rng (random/rng (:seed opts) n 0 0)
          state (-> state (update :queue dissoc k) (update :steps inc) (assoc :time t))]
      (if-let [[state laid stopped?] (local-constraints land state prop opts)]
        (let [seg (select-keys laid [:kind :a :b :bridge?])
              state (add-segment state (assoc seg :t t))]
          (if (or stopped? (:dead-end? laid))
            state
            (reduce (fn [st succ] (enqueue st succ (+ (long t) (long (:delay succ)))))
                    state
                    (global-goals land seg opts rng))))
        state))))

(defn grow
  "The whole road network of `land`."
  ([land] (grow land {}))
  ([land opts]
   (loop [state (start land opts)]
     (if (done? state) state (recur (step land state))))))

(defn segments
  "The laid segments, `{:kind :a :b :bridge? :t}`, each `:t` the step it
  was laid at."
  [state]
  (:segments state))

(defn graph
  "The network as a plane graph: `:points`, and `:edges` as
  `[i j kind]` into them. Junctions are shared exactly, because every
  crossing was made by splitting both roads at one point."
  [state]
  (let [segs (:segments state)
        points (vec (distinct (mapcat (juxt :a :b) segs)))
        index (zipmap points (range))]
    {:points points
     :edges (into [] (comp (map (fn [{:keys [a b kind]}] [(index a) (index b) kind]))
                           (remove (fn [[i j]] (= i j))))
                  segs)}))
