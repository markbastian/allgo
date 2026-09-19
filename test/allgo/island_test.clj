(ns allgo.island-test
  (:require [allgo.procedural.island :as island]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded [seed]
  (let [r (java.util.Random. seed)]
    (fn [] (.nextDouble r))))

(defn- generate
  ([] (generate {}))
  ([opts] (island/generate (merge {:points 400 :seed 3 :rng (seeded 3)} opts))))

(def ^:private world (delay (generate)))

(deftest ranges-test
  ;; The invariant that earns its keep. Elevation and moisture are both
  ;; produced by a breadth-first pass that counts steps and a
  ;; redistribution that maps the result onto [0, 1] -- and the
  ;; redistribution covers a *subset* of the corners, so any corner left
  ;; out keeps a raw step count instead. Nothing downstream notices
  ;; except the biome table, whose thresholds are all fractions: it reads
  ;; off the end of itself and returns the last row. The map still draws,
  ;; and it is simply wrong. Both halves of this caught a real bug.
  (let [{:keys [centers corners]} @world]
    (testing "every corner elevation and moisture is a fraction"
      (is (every? #(<= 0.0 (:elevation %) 1.0) corners))
      (is (every? #(<= 0.0 (:moisture %) 1.0) corners)))

    (testing "and so is every cell's, since a cell averages its corners"
      (is (every? #(<= 0.0 (:elevation %) 1.0) centers))
      (is (every? #(<= 0.0 (:moisture %) 1.0) centers)))

    (testing "the ranges are actually used, not squashed into a corner"
      (is (> (reduce max (map :elevation centers)) 0.6))
      (is (> (reduce max (map :moisture centers)) 0.6)))))

(deftest water-test
  (let [{:keys [centers]} @world
        by-id (into {} (map (juxt :id identity)) centers)]

    (testing "there is an island: some land, and more sea than land"
      (let [land (remove :water? centers)]
        (is (> (count land) (* 0.1 (count centers))))
        (is (< (count land) (* 0.6 (count centers))))))

    (testing "the edge of the map is all ocean"
      ;; Border corners are forced to water, so border cells are water,
      ;; and the flood starts there -- an island that ran off the side
      ;; would have no coastline to speak of.
      (is (every? :ocean? (filter :border? centers))))

    (testing "a lake is water the flood could not reach"
      ;; The definition, and the thing a heightmap cannot express. Every
      ;; lake must be water, must not be ocean, and must not touch one --
      ;; if it touched the ocean the fill would have taken it.
      (let [lakes (filter :lake? centers)]
        (is (every? :water? lakes))
        (is (every? #(not (:ocean? %)) lakes))
        (is (every? (fn [l] (not-any? #(:ocean? (by-id %)) (:neighbors l))) lakes))))

    (testing "ocean and lake exhaust the water"
      (is (every? #(or (:ocean? %) (:lake? %)) (filter :water? centers))))

    (testing "a coast is dry land with the sea next door"
      (let [coast (filter :coast? centers)]
        (is (seq coast))
        (is (every? #(not (:water? %)) coast))
        (is (every? (fn [c] (some #(:ocean? (by-id %)) (:neighbors c))) coast))))))

(deftest elevation-test
  (let [{:keys [centers corners]} @world]

    (testing "the sea is at sea level"
      (is (every? #(zero? (:elevation %)) (filter :ocean? corners))))

    (testing "and so is the shoreline"
      (is (every? #(zero? (:elevation %)) (filter :coast? corners))))

    (testing "the land rises away from the sea"
      ;; Elevation here *is* distance from the coast, so this is close to
      ;; a tautology -- but it is the tautology the whole pass exists to
      ;; produce, and it would stop holding if the queue were walked from
      ;; the wrong end.
      (let [land (remove :water? centers)
            coastal (filter :coast? land)
            inland (remove :coast? land)
            mean (fn [xs] (/ (reduce + (map :elevation xs)) (max 1 (count xs))))]
        (is (seq inland))
        (is (> (mean inland) (mean coastal)))))

    (testing "low ground is commoner than high, which is what the curve is for"
      ;; The breadth-first pass gives a roughly uniform spread; the
      ;; redistribution bends it so most land is low. Without it the
      ;; median sits near the middle instead of well below it.
      (let [es (sort (map :elevation (remove :water? centers)))
            median (nth es (quot (count es) 2))]
        (is (< median 0.45) (str "median land elevation " median))))))

(deftest river-test
  (let [{:keys [corners edges]} @world]

    (testing "rivers exist and are not everywhere"
      (let [wet (filter #(pos? (long (:river %))) edges)]
        (is (seq wet))
        (is (< (count wet) (/ (count edges) 3)))))

    (testing "water runs downhill"
      ;; The downslope pointer is the only thing routing a river, so if
      ;; it ever pointed uphill the rivers would climb.
      (is (every? (fn [v]
                    (<= (:elevation (corners (:downslope v))) (:elevation v)))
                  corners)))

    (testing "and gets to the sea"
      ;; Follow the pointers from every river corner. It has to end at
      ;; the coast or in the ocean -- a river that stopped inland would
      ;; mean a pit the downhill pass had left behind.
      (let [drains? (fn [start]
                      (loop [v start seen #{} steps 0]
                        (cond
                          (or (:coast? v) (:ocean? v)) true
                          (> steps 1000) false
                          (contains? seen (:id v)) false
                          (= (:downslope v) (:id v)) false
                          :else (recur (corners (:downslope v))
                                       (conj seen (:id v))
                                       (inc steps)))))]
        (is (every? drains? (filter #(pos? (long (:river %))) corners)))))

    (testing "a trunk carries more than a headwater"
      ;; Flows add where rivers meet, so the busiest edge should be well
      ;; clear of a single stream.
      (is (> (reduce max 0 (map :river edges)) 2)))))

(deftest river-body-test
  (let [{:keys [rivers corners edges]} @world]

    (testing "each river is a run of corners ending at its mouth"
      (is (seq rivers))
      (is (every? #(= (last (:path %)) (:mouth %)) rivers))
      (is (every? #(>= (count (:path %)) 2) rivers))
      (is (every? #(apply distinct? (:path %)) rivers)))

    (testing "a river runs downhill the whole way"
      ;; It is traced up the downhill pointers, so this is the check that
      ;; the walk never took a branch that climbs -- which is what a
      ;; tributary of the *next* river along would be.
      (is (every? (fn [r] (apply >= (map #(:elevation (corners %)) (:path r))))
                  rivers)))

    (testing "it ends in the sea"
      (is (every? (fn [r] (let [c (corners (:mouth r))]
                            (or (:coast? c) (:ocean? c))))
                  rivers)))

    (testing "at every fork it took the fuller branch"
      ;; The rule that decides which stream is the same river as the one
      ;; below, and the whole reason a river has a length to be labelled
      ;; along. Checked by looking at what else drained into each step and
      ;; confirming nothing carried more than the branch taken.
      (let [drains-into (reduce (fn [m v]
                                  (let [d (:downslope v)]
                                    (if (and (not= d (:id v))
                                             (pos? (long (:river v 0))))
                                      (update m d (fnil conj []) (:id v))
                                      m)))
                                {}
                                corners)]
        (is (every? (fn [r]
                      (every? (fn [[up down]]
                                (>= (long (:river (corners up) 0))
                                    (reduce max 0 (map #(long (:river (corners %) 0))
                                                       (get drains-into down [])))))
                              (partition 2 1 (:path r))))
                    rivers))))

    (testing "the edges named are the ones along the path"
      (is (every? (fn [r] (= (count (:edges r)) (dec (count (:path r))))) rivers))
      (is (every? (fn [r]
                    (every? (fn [[a b]]
                              (some (fn [eid]
                                      (= #{a b} (set (:corners (edges eid)))))
                                    (:edges r)))
                            (partition 2 1 (:path r))))
                  rivers)))

    (testing "rivers come largest first, and carry flow"
      (is (apply >= (map :flow rivers)))
      (is (every? #(pos? (:flow %)) rivers)))

    (testing "a tagged edge belongs to the river that claims it"
      (is (every? (fn [e]
                    (or (nil? (:river-id e))
                        (some #{(:id e)} (:edges (nth rivers (:river-id e))))))
                  edges)))))

(deftest biome-test
  (let [{:keys [centers]} @world]

    (testing "every cell has one"
      (is (every? :biome centers))
      (is (every? #(contains? island/biome-colors (:biome %)) centers)))

    (testing "the water cases win over the table"
      (is (every? #(= :ocean (:biome %)) (filter :ocean? centers)))
      (is (every? #(= :beach (:biome %)) (filter :coast? centers)))
      (is (every? #(contains? #{:lake :marsh :ice} (:biome %))
                  (filter :lake? centers))))

    (testing "dry high ground and wet high ground are different places"
      (is (= :snow (island/biome {:elevation 0.9 :moisture 0.9})))
      (is (= :scorched (island/biome {:elevation 0.9 :moisture 0.05})))
      (is (= :tropical-rain-forest (island/biome {:elevation 0.1 :moisture 0.9})))
      (is (= :subtropical-desert (island/biome {:elevation 0.1 :moisture 0.05}))))

    (testing "more than a handful of kinds of place"
      (is (> (count (distinct (map :biome centers))) 5)))))

(deftest shapes-test
  (testing "the radial shape makes one island, so it has no lakes to find"
    ;; Monotone in radius: there is no interior dip for water to sit in.
    ;; Worth pinning, because it is the difference between the two shapes
    ;; and the reason to keep both.
    (let [m (generate {:shape :radial})]
      (is (empty? (filter :lake? (:centers m))))
      (is (seq (filter :coast? (:centers m))))))

  (testing "the noise shape is allowed to break into pieces"
    (let [m (generate {:shape :noise})]
      (is (seq (remove :water? (:centers m))))))

  (testing "a shape can be handed in directly"
    ;; Everything downstream asks the shape one question, so a caller can
    ;; answer it however they like.
    (let [m (generate {:shape (fn [[nx ny]] (< (+ (* nx nx) (* ny ny)) 0.25))})
          land (remove :water? (:centers m))]
      (is (seq land))
      (is (every? (fn [c]
                    (let [[x y] (:point c)
                          nx (- (/ x 500.0) 1.0)
                          ny (- (/ y 500.0) 1.0)]
                      (< (+ (* nx nx) (* ny ny)) 0.45)))
                  land)))))

(deftest determinism-test
  (testing "the same seed and stream give the same island"
    (is (= (map :biome (:centers (generate)))
           (map :biome (:centers (generate))))))

  (testing "a different seed gives a different one"
    (is (not= (map :biome (:centers (generate {:seed 3 :rng (seeded 3)})))
              (map :biome (:centers (generate {:seed 9 :rng (seeded 9)})))))))
