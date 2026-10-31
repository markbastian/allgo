(ns allgo.procedural.city.agents
  "A town grown by agents, after Lechner, Watson, Wilensky and Felsen,
  \"Procedural City Modeling\".

  Where Parish and Mueller lay roads first and fill them with buildings,
  this grows both together, from a seed, and the town is what simple
  agents leave behind. The world is a grid of patches; agents never talk
  to each other, only change the patches and react to what others have
  changed. Four kinds:

      extenders   roam near developed land looking for a patch the roads
                  do not serve. Every road stamps its distance on the
                  patches around it, so from an unserved patch an
                  extender follows those values back down to the network,
                  avoiding steep steps, and offers the path it took as a
                  road -- if it passes the tests: not too much road
                  around it already, not too close to another junction,
                  not too crooked
      connectors  wander the roads. Each now and then picks a road patch
                  nearby and searches for it along the roads, breadth
                  first; if the way round is more than `k` times the
                  straight Manhattan distance, it offers a road straight
                  there, subject to the same tests. Low `k` gives a tight
                  net of loops, high `k` cul-de-sacs
      developers  residential and commercial. Each roams near the roads,
                  values the empty patch it lands on, and builds if that
                  raises the value; landing on a building of its own kind
                  it builds denser, if that raises the value instead.
                  Commercial value rises with road density and nearness
                  to roads; residential value falls with road density;
                  both like their own kind nearby, water, and ground a
                  little above the average, and dislike a neighbor much
                  denser or sparser than themselves
      speculators the paper's future work, in one paragraph: \"pseudo-
                  developers\" that reserve large tracts of land to sell
                  to developers later, for residential districts,
                  universities and parks. It says no more, so the rules
                  here are this code's own. A speculator roams the edge
                  of town, and on land the roads do not serve yet it
                  reserves a compact tract and holds it: nobody builds
                  there. Rough or waterside land it holds for a park,
                  a little flat land for a campus, the rest for housing.
                  When the roads reach a quarter of the tract, and it
                  has been held a while, it sells. A residential
                  district then takes houses only; a park takes nothing,
                  not even a road, and housing values it as it does
                  water; a campus is built out by an institution, a
                  building at a time with grounds between. A tract the
                  roads never reach is given up

  The paper's grid constraint is `:deviation`: a road may only be laid on
  patches within that many of a grid line, every `:spacing` patches. Near
  zero the town is gridded; large, it is organic. The paper ran in
  NetLogo; this is the same idea in plain arrays, mutated in place, so
  that a browser can step it every frame."
  (:require [allgo.procedural.city.land :as land]
            [allgo.random :as random]
            [clojure.math :as math]))

(def defaults
  {:seed 1
   ;; Patches across, and meters a patch.
   :size 120
   :patch 12.0
   :extenders 12
   :connectors 6
   :residential 40
   :commercial 20
   :speculators 3
   ;; Patches a speculator's tract aims for, the ticks it is held before
   ;; it may sell, and how many campuses a town gets.
   :tract-size 100
   :hold 120
   :campuses 1
   ;; A patch is served by a road this many patches away or nearer.
   :service 3.0
   ;; How far a road stamps its distance.
   :reach 9
   ;; The grid roads keep to, and how far they may stray from it.
   :spacing 6
   :deviation 6
   ;; Connectors build when the way round is more than k times the
   ;; Manhattan distance.
   :k 2.5
   ;; Most road a neighborhood of radius 2 may hold.
   :max-road-density 0.42
   ;; Closest two junctions may be, patches.
   :junction-gap 2.5})

;; Uses of a patch.
(def empty-use 0)
(def residential 1)
(def commercial 2)
(def institution 3)

(defn- idx ^long [w i j] (+ (long i) (* (long w) (long j))))

(defn- inside? [w i j] (and (< -1 (long i) (long w)) (< -1 (long j) (long w))))

(def ^:private neighbors8 [[1 0] [-1 0] [0 1] [0 -1] [1 1] [1 -1] [-1 1] [-1 -1]])
(def ^:private neighbors4 [[1 0] [-1 0] [0 1] [0 -1]])

;; ---------------------------------------------------------------------------
;; The world

(defn- mean-elevation ^double [{:keys [elev water]}]
  (let [n (alength ^doubles elev)
        [s c] (reduce (fn [[s c] k] (if (zero? (aget ^ints water k)) [(+ s (aget ^doubles elev k)) (inc c)] [s c]))
                      [0.0 0] (range n))]
    (if (pos? c) (/ s c) 0.0)))

(defn- stamp!
  "Marks `(i, j)` as road, clearing whatever was built there, and lowers
  the distance to the nearest road of every patch within reach."
  [{:keys [w road dist use density opts]} i j]
  (let [r (long (:reach opts))
        k (idx w i j)]
    (aset ^ints road k 1)
    (aset ^ints use k 0)
    (aset ^ints density k 0)
    (doseq [di (range (- r) (inc r)) dj (range (- r) (inc r))
            :let [a (+ (long i) di) b (+ (long j) dj)]
            :when (inside? w a b)]
      (let [k (idx w a b)
            d (math/sqrt (+ (* di di) (* dj dj)))]
        (when (< d (aget ^doubles dist k)) (aset ^doubles dist k d))))))

(defn- road? [{:keys [w road]} i j]
  (and (inside? w i j) (= 1 (aget ^ints road (idx w i j)))))

(defn- lay!
  "Lays a road along `path`, a sequence of `[i j]`, and records its ends
  as junctions."
  [world path]
  (doseq [[i j] path] (stamp! world i j))
  (swap! (:junctions world) conj (first path) (last path))
  (swap! (:roads world) conj (vec path))
  world)

(defn world
  "A new world on `land-opts`' land, with a short road through the
  middle and agents scattered along it."
  ([] (world {}))
  ([opts]
   (let [opts (merge defaults opts)
         {:keys [seed size patch]} opts
         w (long size)
         half (* 0.5 w (double patch))
         l (land/land {:seed seed :size half :river? true :coast? true :pattern :organic})
         n (* w w)
         elev (double-array n)
         ;; Every int array is given its zero: ClojureScript's
         ;; `(int-array n)` is n nils.
         water (int-array n 0)
         rng (random/rng seed 61 0 0)]
     (doseq [i (range w) j (range w)]
       (let [x (- (* (+ i 0.5) patch) half) y (- (* (+ j 0.5) patch) half)
             k (idx w i j)]
         (aset elev k (double ((:elevation l) x y)))
         (aset water k (if ((:water? l) x y) 1 0))))
     (let [world {:w w :opts opts :land l :rng rng
                  :elev elev :water water
                  :mean-elevation (mean-elevation {:elev elev :water water})
                  :road (int-array n 0)
                  :dist (double-array n ##Inf)
                  :use (int-array n 0)
                  :density (int-array n 0)
                  ;; The tract each patch is reserved in, 0 for none,
                  ;; and the tracts by id.
                  :tract (int-array n 0)
                  :tracts (atom (sorted-map))
                  :junctions (atom #{})
                  :roads (atom [])
                  :agents (atom [])
                  :ticks (atom 0)}
           ;; The seed: the dry patch nearest the main center, and a road
           ;; across it.
           [cx cy] (:center l)
           sp (long (:spacing opts))
           ;; On a grid line, so that a grid can grow from it.
           ci (* sp (math/round (/ (/ (+ cx half) patch) sp)))
           cj (* sp (math/round (/ (/ (+ cy half) patch) sp)))
           seed-path (vec (for [i (range (- ci 5) (+ ci 6))
                                :when (and (inside? w i cj) (zero? (aget water (idx w i cj))))]
                            [i cj]))]
       (lay! world seed-path)
       (reset! (:agents world)
               (vec (for [[kind n] [[:extender (:extenders opts)] [:connector (:connectors opts)]
                                    [:residential (:residential opts)] [:commercial (:commercial opts)]
                                    [:speculator (:speculators opts)]]
                          _ (range n)]
                      {:kind kind :at (random/pick rng seed-path)})))
       world))))

;; ---------------------------------------------------------------------------
;; What a patch is worth

(defn- around
  "The patches within `r` of `(i, j)`, inside the world."
  [w i j r]
  (for [di (range (- r) (inc r)) dj (range (- r) (inc r))
        :let [a (+ (long i) di) b (+ (long j) dj)]
        :when (and (inside? w a b) (<= (+ (* di di) (* dj dj)) (* r r)))]
    [a b]))

(defn road-density
  "The fraction of the patches within 2 of `(i, j)` that are road."
  ^double [{:keys [w] :as world} i j]
  (let [ps (around w i j 2)]
    (/ (count (filter (fn [[a b]] (road? world a b)) ps)) (double (count ps)))))

(defn- held
  "The tract `(i, j)` is reserved in, or nil."
  [{:keys [w tract tracts]} i j]
  (let [id (aget ^ints tract (idx w i j))]
    (when (pos? id) (get @tracts id))))

(defn- park? [world i j]
  (= :park (:kind (held world i j))))

(defn- zoned-for?
  "Whether a developer of `kind` may build on `(i, j)`: anywhere no
  speculator holds, and in a residential district once it is sold, if
  it builds houses."
  [world kind i j]
  (let [t (held world i j)]
    (or (nil? t)
        (and (:sold? t) (= :residential (:kind t)) (= :residential kind)))))

(defn value
  "What `(i, j)` is worth to a developer of `kind`, built to `density`."
  [{:keys [w use water elev dist] :as world} kind i j density]
  (let [ps (around w i j 3)
        k (idx w i j)
        own (if (= kind :residential) residential commercial)
        other (if (= kind :residential) commercial residential)
        frac (fn [u] (/ (count (filter (fn [[a b]] (= u (aget ^ints use (idx w a b)))) ps)) (double (count ps))))
        near-water (if (some (fn [[a b]] (= 1 (aget ^ints water (idx w a b)))) ps) 1.0 0.0)
        ;; Housing likes a park as it likes the water.
        near-park (if (and (= kind :residential) (some (fn [[a b]] (park? world a b)) ps)) 1.0 0.0)
        rise (- (aget ^doubles elev k) (double (:mean-elevation world)))
        ;; A little above the average is best; well above is a climb.
        r (/ rise 12.0)
        height (- (* 0.25 (max -1.0 (min 1.0 r))) (* 0.2 (max 0.0 (- r 1.5))))
        rd (road-density world i j)
        d (aget ^doubles dist k)
        ;; Neighbors of the same kind, and how far their density is from
        ;; this one's.
        same (filter (fn [[a b]] (= own (aget ^ints use (idx w a b)))) ps)
        disparity (if (seq same)
                    (abs (- (double density)
                            (/ (reduce + (map (fn [[a b]] (aget ^ints (:density world) (idx w a b))) same))
                               (double (count same)))))
                    0.0)]
    (+ 0.2 (* 0.4 near-water) (* 0.4 near-park) height
       (* 0.9 (frac own)) (* -0.4 (frac other))
       (* -0.08 disparity)
       ;; Building denser pays less with every storey: more where the
       ;; roads are thick, for commerce.
       (let [dn (double density)]
         (if (= kind :commercial)
           (- (* 1.5 rd dn) (* 0.04 dn dn))
           (- (* 0.12 dn) (* 0.02 dn dn))))
       (if (= kind :commercial)
         (+ (* 2.2 rd) (* -0.25 d))
         (+ (* -2.5 (max 0.0 (- rd 0.12))) (* -0.05 (abs (- d 2.0))))))))

;; ---------------------------------------------------------------------------
;; Tests every road must pass

(defn- on-grid?
  "Whether every patch of `path` lies within the deviation allowed of one
  and the same grid line, across or up."
  [{:keys [opts]} path]
  (let [{:keys [spacing deviation]} opts
        spacing (long spacing) deviation (long deviation)
        near (fn [coords]
               ;; Some line, a multiple of the spacing, that all of them
               ;; are near.
               (let [lo (reduce min coords) hi (reduce max coords)]
                 (some (fn [line] (<= (- hi deviation) line (+ lo deviation)))
                       (map #(* spacing %) (range (quot (- lo deviation spacing) spacing)
                                                  (inc (quot (+ hi deviation spacing) spacing)))))))]
    (or (>= (* 2 deviation) spacing)
        (near (map first path))
        (near (map second path)))))

(defn- organic? [{:keys [opts]}]
  (>= (* 2 (long (:deviation opts))) (long (:spacing opts))))

(defn- acceptable?
  "The paper's constraint tests on a proposed road `path`: dry, on the
  grid within the deviation allowed, not into a neighborhood already
  thick with road, and not ending too near another junction -- nor, a
  rule of this code's, through a park."
  [{:keys [w water opts junctions] :as world} path]
  (let [{:keys [max-road-density junction-gap]} opts
        [i j] (nth path (quot (count path) 2))
        end (last path)
        start (first path)
        ;; Where the new road meets the network.
        joins (filter (fn [[a b]] (road? world a b)) [start end])]
    (and (>= (count path) 2)
         (every? (fn [[a b]] (and (inside? w a b) (zero? (aget ^ints water (idx w a b)))
                                  (not (park? world a b))))
                 path)
         (on-grid? world path)
         (< (road-density world i j) (double max-road-density))
         (not-any? (fn [[[a b] [x y]]]
                     (let [d (math/sqrt (+ (math/pow (- a x) 2) (math/pow (- b y) 2)))]
                       (and (pos? d) (< d (double junction-gap)))))
                   (for [jn @junctions p joins] [jn p]))
         (not= start end))))

;; ---------------------------------------------------------------------------
;; Agents

(defn- wander
  "A step to a random neighbor that `ok?` allows, or stay put."
  [rng [i j] ok?]
  (let [cands (filter (fn [[a b]] (ok? a b)) (map (fn [[di dj]] [(+ (long i) di) (+ (long j) dj)]) neighbors8))]
    (if (seq cands) (random/pick rng (vec cands)) [i j])))

(defn- trace-to-road
  "From `(i, j)`, down the stamped distances to the road network: each
  step to the neighbor nearest a road, the gentler climb breaking ties.
  The path, ending on a road patch, or nil."
  [{:keys [w dist elev water] :as world} i j]
  (loop [p [i j] path [[i j]]]
    (let [[a b] p]
      (cond
        (road? world a b) path
        (> (count path) 30) nil
        :else
        (let [here (aget ^doubles dist (idx w a b))
              e0 (aget ^doubles elev (idx w a b))
              ;; On a grid, only straight steps: the distances then
              ;; lead straight across to the nearest road.
              steps (for [[di dj] (if (organic? world) neighbors8 neighbors4)
                          :let [x (+ (long a) di) y (+ (long b) dj)]
                          :when (and (inside? w x y) (zero? (aget ^ints water (idx w x y))))
                          :let [k (idx w x y) d (aget ^doubles dist k)]
                          :when (< d here)]
                      [[x y] (+ d (* 0.08 (abs (- (aget ^doubles elev k) e0))))])]
          (when (seq steps)
            (let [[q] (apply min-key second steps)]
              (recur q (conj path q)))))))))

(defn- developed-near? [{:keys [w use]} i j r]
  (some (fn [[a b]] (pos? (aget ^ints use (idx w a b)))) (around w i j r)))

(defn- extender [{:keys [w dist opts rng water] :as world} {:keys [at] :as agent}]
  (let [{:keys [service reach]} opts
        ok? (fn [a b] (and (inside? w a b) (zero? (aget ^ints water (idx w a b)))
                           (< (aget ^doubles dist (idx w a b)) (double reach))))
        [i j] (wander rng at ok?)
        d (aget ^doubles dist (idx w i j))]
    ;; An unserved patch, near enough to the town to be worth serving.
    (when (and (> d (double service)) (developed-near? world i j 5))
      (when-let [path (trace-to-road world i j)]
        ;; Not too crooked: no longer than half again the straight line.
        (let [[ei ej] (last path)
              straight (math/sqrt (+ (math/pow (- ei i) 2) (math/pow (- ej j) 2)))]
          (when (and (<= (count path) (+ 2 (* 1.5 straight))) (acceptable? world (vec (reverse path))))
            (lay! world (vec (reverse path)))))))
    (assoc agent :at [i j])))

(defn- road-path
  "Whether `to` is reachable from `from` along road patches in at most
  `limit` steps, by breadth-first search."
  [world from to limit]
  (loop [frontier [from] seen #{from} depth 0]
    (cond
      (contains? seen to) true
      (or (empty? frontier) (>= depth (long limit))) false
      :else
      (let [nxt (for [[a b] frontier [di dj] neighbors4
                      :let [q [(+ (long a) di) (+ (long b) dj)]]
                      :when (and (not (seen q)) (road? world (first q) (second q)))]
                  q)
            nxt (vec (distinct nxt))]
        (recur nxt (into seen nxt) (inc depth))))))

(defn- line
  "The patches on the straight line from `a` to `b`, Bresenham's, four-
  connected so a road along it has no diagonal gaps."
  [[i0 j0] [i1 j1]]
  (let [di (abs (- i1 i0)) dj (abs (- j1 j0))
        si (if (< i0 i1) 1 -1) sj (if (< j0 j1) 1 -1)]
    (loop [i i0 j j0 err (- di dj) out [[i0 j0]]]
      (if (and (= i i1) (= j j1))
        out
        (let [e2 (* 2 err)]
          ;; One axis per step.
          (if (> e2 (- dj))
            (recur (+ i si) j (- err dj) (conj out [(+ i si) j]))
            (recur i (+ j sj) (+ err di) (conj out [i (+ j sj)]))))))))

(defn- connector [{:keys [w opts rng] :as world} {:keys [at] :as agent}]
  (let [[i j] (wander rng at (fn [a b] (road? world a b)))]
    (when (< (double (rng)) 0.3)
      (let [r 8
            ti (+ i (- (random/below rng (inc (* 2 r))) r))
            tj (+ j (- (random/below rng (inc (* 2 r))) r))
            manhattan (+ (abs (- ti i)) (abs (- tj j)))]
        (when (and (inside? w ti tj) (road? world ti tj) (> manhattan 3)
                   (not (road-path world [i j] [ti tj] (* (double (:k opts)) manhattan))))
          (let [path (line [i j] [ti tj])
                ;; Only the stretch off the roads is new.
                fresh (filterv (fn [[a b]] (not (road? world a b))) path)]
            (when (and (seq fresh) (acceptable? world path))
              (lay! world path))))))
    (assoc agent :at [i j])))

(defn- developer [{:keys [w dist use density water opts rng] :as world} {:keys [kind at] :as agent}]
  (let [service (double (:service opts))
        ok? (fn [a b] (and (inside? w a b) (zero? (aget ^ints water (idx w a b)))
                           (<= (aget ^doubles dist (idx w a b)) (+ 2.0 service))))
        ;; A few steps at a time: developers cover ground.
        [i j] (nth (iterate #(wander rng % ok?) at) 3)
        k (idx w i j)
        own (if (= kind :residential) residential commercial)
        u (aget ^ints use k)
        dn (aget ^ints density k)]
    (when (and (not (road? world i j)) (<= (aget ^doubles dist k) service) (zoned-for? world kind i j))
      (cond
        (and (= u empty-use) (pos? (value world kind i j 1)))
        (do (aset ^ints use k (int own)) (aset ^ints density k 1))

        (and (= u own) (< dn 5) (> (value world kind i j (inc dn)) (value world kind i j dn)))
        (aset ^ints density k (int (inc dn)))))
    (assoc agent :at [i j])))

;; ---------------------------------------------------------------------------
;; Speculators, and the institution that builds a campus

(defn- free? [{:keys [w water road use tract]} i j]
  (let [k (idx w i j)]
    (and (zero? (aget ^ints water k)) (zero? (aget ^ints road k))
         (zero? (aget ^ints use k)) (zero? (aget ^ints tract k)))))

(defn- survey
  "The tract a speculator at `(i, j)` could reserve: the free patches
  joined to it within the disc of the tract size, or nil if they are
  too few."
  [{:keys [w opts] :as world} i j]
  (let [size (long (:tract-size opts))
        r (math/sqrt (/ size math/PI))
        in? (fn [a b] (<= (+ (* (- a i) (- a i)) (* (- b j) (- b j))) (* r r)))]
    (when (free? world i j)
      (loop [frontier [[i j]] seen #{[i j]}]
        (if (seq frontier)
          (let [nxt (vec (distinct (for [[a b] frontier [di dj] neighbors4
                                         :let [x (+ (long a) di) y (+ (long b) dj)]
                                         :when (and (inside? w x y) (not (seen [x y])) (in? x y) (free? world x y))]
                                     [x y])))]
            (recur nxt (into seen nxt)))
          (when (>= (count seen) (* 0.6 size)) (vec (sort seen))))))))

(defn- zoning
  "What a speculator holds `patches` for: a park on rough ground or the
  waterside, half the time; a campus on flat ground while the town has
  fewer than it allows; otherwise a residential district."
  [{:keys [w water elev opts tracts rng]} patches]
  (let [es (map (fn [[i j]] (aget ^doubles elev (idx w i j))) patches)
        relief (- (double (reduce max es)) (double (reduce min es)))
        wet? (some (fn [[i j]] (some (fn [[a b]] (= 1 (aget ^ints water (idx w a b)))) (around w i j 2))) patches)
        campuses (count (filter #(= :campus (:kind %)) (vals @tracts)))]
    (cond
      (and (or wet? (> relief 9.0)) (< (double (rng)) 0.5)) :park
      (and (< campuses (long (:campuses opts))) (< relief 6.0)) :campus
      :else :residential)))

(defn- release! [{:keys [w tract tracts]} id]
  (doseq [[i j] (:patches (@tracts id))
          :when (= id (aget ^ints tract (idx w i j)))]
    (aset ^ints tract (idx w i j) 0))
  (swap! tracts dissoc id))

(defn- speculator [{:keys [w dist water opts rng tract tracts ticks] :as world} {:keys [at holding] :as agent}]
  (if holding
    ;; Holding: sell once the roads reach a quarter of it and it has been
    ;; held long enough; give it up if they never do.
    (let [{:keys [patches since]} (@tracts holding)
          held-for (- (long @ticks) (long since))
          served (count (filter (fn [[a b]] (<= (aget ^doubles dist (idx w a b)) (double (:service opts)))) patches))]
      (cond
        (and (>= held-for (long (:hold opts))) (>= served (* 0.25 (count patches))))
        (do (swap! tracts assoc-in [holding :sold?] true) (dissoc agent :holding))

        (> held-for (* 5 (long (:hold opts))))
        (do (release! world holding) (dissoc agent :holding))

        :else agent))
    ;; Roaming the edge of town, for land the roads do not serve yet.
    (let [ok? (fn [a b] (and (inside? w a b) (zero? (aget ^ints water (idx w a b)))
                             (< (aget ^doubles dist (idx w a b)) (double (:reach opts)))))
          [i j] (nth (iterate #(wander rng % ok?) at) 2)
          patches (when (and (> (aget ^doubles dist (idx w i j)) (double (:service opts)))
                             (< (double (rng)) 0.15))
                    (survey world i j))]
      (if patches
        (let [id (inc (long (or (first (keys (rseq @tracts))) 0)))]
          (doseq [[a b] patches] (aset ^ints tract (idx w a b) (int id)))
          (swap! tracts assoc id {:id id :kind (zoning world patches) :patches patches
                                  :since @ticks :sold? false})
          (assoc agent :at [i j] :holding id))
        (assoc agent :at [i j])))))

(defn- institution!
  "A sold campus's builder: a building on a served patch of it, if none
  stands beside it, so that the campus keeps its grounds."
  [{:keys [w use density dist opts rng] :as world} {:keys [patches]}]
  (let [[i j] (random/pick rng patches)
        k (idx w i j)]
    (when (and (not (road? world i j)) (zero? (aget ^ints use k))
               (<= (aget ^doubles dist k) (double (:service opts)))
               (not-any? (fn [[di dj]]
                           (let [a (+ (long i) di) b (+ (long j) dj)]
                             (and (inside? w a b) (= institution (aget ^ints use (idx w a b))))))
                         neighbors8))
      (aset ^ints use k (int institution))
      (aset ^ints density k (int (+ 2 (long (random/below rng 3))))))))

(defn tick!
  "Every agent takes one turn, and every sold campus builds. Returns the
  world."
  [world]
  (swap! (:agents world)
         (fn [agents]
           (mapv (fn [{:keys [kind] :as a}]
                   (case kind
                     :extender (extender world a)
                     :connector (connector world a)
                     :speculator (speculator world a)
                     (developer world a)))
                 agents)))
  (doseq [t (vals @(:tracts world))
          :when (and (:sold? t) (= :campus (:kind t)))]
    (institution! world t))
  (swap! (:ticks world) inc)
  world)

(defn run
  "A world after `n` ticks."
  ([n] (run n {}))
  ([n opts]
   (let [world (world opts)]
     (dotimes [_ n] (tick! world))
     world)))
