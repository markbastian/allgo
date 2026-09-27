(ns allgo.dungeons-test
  (:require [allgo.graph :as graph]
            [allgo.procedural.dungeons :as d]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(defn- room [id [cx cy] w h]
  {:id id :center [cx cy] :width w :height h})

(defn- overlapping-pairs [rooms]
  (for [a rooms b rooms
        :when (and (< (:id a) (:id b)) (d/overlaps? a b))]
    [(:id a) (:id b)]))

(defn- on-grid? [tile-size room]
  (every? #(zero? (mod % tile-size)) (d/bounds room)))

;; ---------------------------------------------------------------------------

(deftest roundm-test
  (testing "rounds up to the next multiple, as in the original"
    (is (= 8 (d/roundm 5 4)))
    (is (= 4 (d/roundm 4 4)))
    (is (= 4 (d/roundm 1 4)))
    (is (= 0 (d/roundm 0 4))))
  (testing "negatives"
    (is (= 0 (d/roundm -3 4)))
    (is (= -4 (d/roundm -4 4)))))

(deftest random-point-test
  (testing "stays inside the disk"
    (let [r 50.0]
      (is (every? (fn [[x y]] (<= (Math/hypot x y) (+ r 1e-9)))
                  (repeatedly 500 #(d/random-point r))))))

  (testing "is not crowded toward the center"
    ;; Sampling the radius uniformly would put half the points inside
    ;; r/2, which covers only a quarter of the area. A correct sampler
    ;; puts about a quarter of them there.
    (let [r      50.0
          inner  (count (filter (fn [[x y]] (< (Math/hypot x y) (/ r 2)))
                                (repeatedly 4000 #(d/random-point r))))]
      (is (< 0.15 (/ inner 4000.0) 0.35))))

  (testing "snapped variant lands on the grid"
    (is (every? (fn [[x y]] (and (zero? (mod x 4)) (zero? (mod y 4))))
                (repeatedly 200 #(d/random-int-point 40 4))))))

(deftest bounds-test
  (is (= [-10.0 -5.0 10.0 5.0] (d/bounds (room 0 [0 0] 20 10))))
  (is (d/overlaps? (room 0 [0 0] 20 20) (room 1 [10 0] 20 20)))
  (testing "rooms merely touching do not overlap"
    (is (not (d/overlaps? (room 0 [0 0] 20 20) (room 1 [20 0] 20 20))))))

;; ---------------------------------------------------------------------------

(deftest separate-test
  (let [config (merge d/default-config {:tile-size 4 :max-iterations 400})]
    (testing "a deliberate pile-up ends up disjoint and on the grid"
      (let [piled  (mapv #(room % [0 0] 16 16) (range 12))
            spread (d/separate piled config)]
        (is (empty? (overlapping-pairs spread)))
        (is (every? (partial on-grid? 4) spread))
        (is (= 12 (count spread)) "no room is lost")
        (is (= (set (map :id piled)) (set (map :id spread))))))

    (testing "rooms already disjoint are left where they are"
      (let [placed (mapv (fn [i] (room i [(* i 40) 0] 16 16)) (range 5))]
        (is (= (map :center placed) (map :center (d/separate placed config))))))

    (testing "a scattered pile separates for a range of sizes"
      (doseq [n [20 60]]
        (let [cfg    (merge config {:room-count n :radius (* 1.2 n)})
              spread (d/separate (d/scatter cfg) cfg)]
          (is (empty? (overlapping-pairs spread)) (str n " rooms"))
          (is (every? (partial on-grid? 4) spread) (str n " rooms")))))))

(deftest scatter-test
  (let [cfg   (merge d/default-config {:room-count 200 :tile-size 4})
        rooms (d/scatter cfg)]
    (is (= 200 (count rooms)))
    (is (= (range 200) (map :id rooms)) "ids are stable and dense")
    (testing "extents are even multiples of the tile size, so edges stay on grid"
      (is (every? #(zero? (mod (:width %) 8)) rooms))
      (is (every? #(zero? (mod (:height %) 8)) rooms)))
    (testing "no room collapses below the minimum"
      (is (every? #(pos? (:width %)) rooms))
      (is (every? #(pos? (:height %)) rooms)))))

(deftest main-rooms-test
  (testing "a room is a hub only if it clears the threshold on both axes"
    (let [rooms  [(room 0 [0 0] 10 10)
                  (room 1 [0 0] 10 10)
                  (room 2 [0 0] 30 30)   ; big both ways
                  (room 3 [0 0] 30 10)]  ; wide but squat
          tagged (d/main-rooms rooms {:main-threshold 1.25})
          kind   (into {} (map (juxt :id :kind)) tagged)]
      ;; means are 20 and 15, so thresholds are 25 and 18.75
      (is (= :minor (kind 0)))
      (is (= :minor (kind 1)))
      (is (= :main (kind 2)))
      (is (= :minor (kind 3)) "must clear the bar on height too")))

  (testing "hubs pulls out exactly the tagged rooms"
    (let [tagged (d/main-rooms [(room 0 [0 0] 10 10) (room 1 [0 0] 40 40)]
                               {:main-threshold 1.25})]
      (is (= [1] (map :id (d/hubs tagged)))))))

;; ---------------------------------------------------------------------------

(deftest room-graph-test
  (testing "four hubs in a square triangulate into a connected graph"
    (let [hubs  [(room 0 [0 0] 20 20) (room 1 [100 0] 20 20)
                 (room 2 [100 100] 20 20) (room 3 [0 100] 20 20)]
          {:keys [nodes edges]} (d/room-graph hubs)]
      (is (= #{0 1 2 3} nodes))
      ;; A square triangulates to two triangles: four sides and one diagonal.
      (is (= 5 (count edges)))
      (is (graph/connected-graph? nodes edges))
      (is (every? #(= 2 (count %)) edges) "edges are unordered pairs")))

  (testing "fewer than three hubs falls back to connecting them directly"
    (is (= #{#{0 1}} (:edges (d/room-graph [(room 0 [0 0] 10 10)
                                            (room 1 [50 0] 10 10)]))))
    (is (= #{} (:edges (d/room-graph [(room 0 [0 0] 10 10)]))))
    (is (= #{} (:edges (d/room-graph [])))))

  (testing "collinear hubs still get connected"
    ;; Three rooms in a row have no circumcircle, so Delaunay degenerates
    ;; and returns nothing; without a fallback the dungeon has no corridors.
    (let [hubs [(room 0 [-112 -28] 20 20)
                (room 1 [-20 -32] 20 20)
                (room 2 [100 -36] 20 20)]
          {:keys [nodes edges]} (d/room-graph hubs)]
      (is (= 3 (count nodes)))
      (is (pos? (count edges)))
      (is (graph/connected-graph? nodes edges)))
    (testing "exactly collinear, the degenerate case"
      (let [hubs [(room 0 [0 0] 20 20) (room 1 [100 0] 20 20) (room 2 [200 0] 20 20)]
            {:keys [nodes edges]} (d/room-graph hubs)]
        (is (graph/connected-graph? nodes edges)))))

  (testing "hubs sharing a center cannot both be triangulation sites"
    (let [{:keys [nodes]} (d/room-graph [(room 0 [0 0] 10 10)
                                         (room 1 [0 0] 10 10)
                                         (room 2 [50 0] 10 10)])]
      (is (= 2 (count nodes))))))

(deftest connect-test
  (let [hubs  [(room 0 [0 0] 20 20) (room 1 [100 0] 20 20)
               (room 2 [100 100] 20 20) (room 3 [0 100] 20 20)]
        graph (d/room-graph hubs)]
    (testing "with no edges added back the result is exactly a spanning tree"
      (let [edges (d/connect hubs graph {:extra-edge-ratio 0.0})]
        (is (= 3 (count edges)))
        (is (graph/connected-graph? (:nodes graph) edges))))

    (testing "adding all of them back returns the whole Delaunay graph"
      (is (= (:edges graph) (d/connect hubs graph {:extra-edge-ratio 1.0}))))

    (testing "every hub stays reachable whatever the ratio"
      (doseq [r [0.0 0.1 0.5 1.0]]
        (is (graph/connected-graph? (:nodes graph) (d/connect hubs graph {:extra-edge-ratio r}))
            (str "ratio " r))))

    (testing "the tree is always a subset of what comes back"
      (let [tree (d/connect hubs graph {:extra-edge-ratio 0.0})]
        (is (every? (d/connect hubs graph {:extra-edge-ratio 0.5}) tree))))))

;; ---------------------------------------------------------------------------

(deftest corridor-segments-test
  (testing "rooms sharing a vertical band get one vertical run"
    (let [segs (d/corridor-segments (room 0 [0 0] 40 40) (room 1 [0 100] 40 40))]
      (is (= 1 (count segs)))
      (is (= [[0.0 0] [0.0 100]] (first segs)))))

  (testing "rooms sharing a horizontal band get one horizontal run"
    (let [segs (d/corridor-segments (room 0 [0 0] 40 40) (room 1 [100 0] 40 40))]
      (is (= 1 (count segs)))
      (is (= [[0 0.0] [100 0.0]] (first segs)))))

  (testing "rooms sharing neither get an L of two runs"
    (let [segs (d/corridor-segments (room 0 [0 0] 10 10) (room 1 [100 100] 10 10))]
      (is (= 2 (count segs)))
      (is (= [[[0 0] [100 0]] [[100 0] [100 100]]] segs))))

  (testing "every segment produced is axis-aligned"
    (doseq [[a b] [[(room 0 [0 0] 40 40) (room 1 [0 100] 40 40)]
                   [(room 0 [0 0] 40 40) (room 1 [100 0] 40 40)]
                   [(room 0 [0 0] 10 10) (room 1 [100 100] 10 10)]
                   [(room 0 [0 0] 10 10) (room 1 [-70 -90] 10 10)]]]
      (doseq [[[x1 y1] [x2 y2]] (d/corridor-segments a b)]
        (is (or (== x1 x2) (== y1 y2)))))))

(deftest absorb-test
  (let [config {:tile-size 4 :corridor-width 3}
        rooms  [(assoc (room 0 [0 0] 20 20) :kind :main)
                (assoc (room 1 [50 0] 10 10) :kind :minor)   ; on the run
                (assoc (room 2 [50 500] 10 10) :kind :minor) ; far away
                (assoc (room 3 [100 0] 20 20) :kind :main)]
        segs   [[[0 0] [100 0]]]
        out    (into {} (map (juxt :id :kind)) (d/absorb rooms segs config))]
    (is (= :hallway (out 1)) "a minor room the corridor crosses is promoted")
    (is (= :minor (out 2)) "one nowhere near it is not")
    (is (= :main (out 0)) "hubs are never demoted")
    (is (= :main (out 3)))))

(deftest rasterize-test
  (let [config {:tile-size 4 :corridor-width 3}
        rooms  [(assoc (room 0 [0 0] 16 16) :kind :main)
                (assoc (room 1 [100 0] 8 8) :kind :hallway)
                (assoc (room 2 [500 500] 8 8) :kind :minor)]
        grid   (d/rasterize rooms [[[0 0] [100 0]]] config)]
    (testing "a 16x16 room on a 4-tile grid covers 4x4 tiles"
      (is (= 16 (count (filter #{:main} (vals grid))))))
    (testing "leftover minor rooms are not on the map"
      (is (not (contains? grid [125 125]))))
    (testing "rooms win over the corridor running through them"
      (is (= :main (grid [0 0]))))
    (testing "the corridor is laid down between the two"
      (is (some #{:corridor} (vals grid))))
    (testing "tiles are integer indices"
      (is (every? (fn [[x y]] (and (integer? x) (integer? y))) (keys grid))))))

;; ---------------------------------------------------------------------------

(deftest generate-test
  (doseq [n [40 80]]
    (let [{:keys [rooms graph edges corridors grid config]}
          (d/generate {:room-count n :radius (* 1.2 n)})]
      (testing (str n " rooms: the layout is sound")
        (is (= n (count rooms)))
        (is (empty? (overlapping-pairs rooms)) "no two rooms overlap")
        (is (every? (partial on-grid? (:tile-size config)) rooms)))

      (testing (str n " rooms: every hub is reachable")
        (is (graph/connected-graph? (:nodes graph) edges)))

      (testing (str n " rooms: the connected edges come from the Delaunay graph")
        (is (every? (:edges graph) edges)))

      (testing (str n " rooms: there is something to walk around in")
        (is (<= 3 (count (d/hubs rooms))) "the min-hubs floor holds")
        (is (pos? (count corridors)))
        (is (pos? (count grid)))
        (is (contains? (set (vals grid)) :main))
        ;; Hallway rooms appear only where a corridor happens to cross a
        ;; leftover room, so they are not guaranteed on any given run.
        (is (every? #{:main :hallway :corridor} (vals grid))))

      (testing (str n " rooms: only tagged kinds survive")
        (is (every? #{:main :hallway :minor} (map :kind rooms)))))))

(deftest generate-config-test
  (testing "defaults alone produce a dungeon"
    (is (pos? (count (:grid (d/generate))))))

  (testing "an ellipse scatter gives a wider layout than a disk"
    ;; Adonaac's fix for squat dungeons: seed a strip, not a circle.
    (let [aspect (fn [cfg]
                   (let [rooms (:rooms (d/generate cfg))
                         xs    (mapcat (fn [r] [(first (d/bounds r)) (nth (d/bounds r) 2)]) rooms)
                         ys    (mapcat (fn [r] [(second (d/bounds r)) (nth (d/bounds r) 3)]) rooms)]
                     (/ (- (apply max xs) (apply min xs))
                        (double (- (apply max ys) (apply min ys))))))]
      (is (> (aspect {:room-count 60 :ellipse [800 60]})
             (aspect {:room-count 60 :radius 80.0}))))))

(deftest seeded-test
  (testing "the same seed builds the same dungeon"
    (let [build #(d/generate {:room-count 40 :radius 60.0 :rng (random/rng 11)})]
      (is (= (:grid (build)) (:grid (build))))
      (is (= (:edges (build)) (:edges (build)))))))
