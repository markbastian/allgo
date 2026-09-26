(ns allgo.erosion-test
  (:require [allgo.procedural.erosion :as e]
            [allgo.procedural.terrain :as t]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded
  "A repeatable stand-in for `rand`."
  [seed]
  (let [r (java.util.Random. seed)]
    (fn [] (.nextDouble r))))

(defn- flat
  "A level field, `dim` on a side, every cell at `h`."
  [dim h]
  {:heights (double-array (* dim dim) (double h)) :dim dim})

(defn- ramp
  "A field tilting `drop` per cell toward the last row."
  [dim drop]
  {:heights (double-array (for [i (range dim) _ (range dim)] (* (- (dec dim) i) (double drop))))
   :dim dim})

(defn- terrain
  ([] (terrain 7 7))
  ([iterations seed]
   (t/generate {:width 1.0 :iterations iterations :corners [0.0 0.0 0.0 0.0]
                :rng (seeded seed)})))

(defn- heights [{:keys [heights]}] (vec ^doubles heights))

(defn- total [{:keys [^doubles heights]}]
  (areduce heights i s 0.0 (+ s (aget heights i))))

(defn- extremes [{:keys [^doubles heights]}]
  [(reduce min (vec heights)) (reduce max (vec heights))])

(deftest brush-test
  (testing "the weights are a partition: a brush moves what it is given"
    (doseq [r [0 1 2 3 5]]
      (is (< (abs (- 1.0 (reduce + (vec ^doubles (:w (e/brush r)))))) 1e-12)
          (str "radius " r))))

  (testing "radius zero is the single cell under the droplet"
    (let [b (e/brush 0)]
      (is (= 1 (:n b)))
      (is (= [0] (vec ^ints (:dr b))))
      (is (= [0] (vec ^ints (:dc b))))))

  (testing "it is a disk, not a square, and rim cells that take nothing are dropped"
    ;; A 7x7 square would be 49; the disk of radius 3 is 29, and the 8
    ;; cells at distance exactly 3 carry zero weight and are not visited.
    (is (= 29 (:n (e/brush 3))))
    (is (every? pos? (vec ^doubles (:w (e/brush 3))))))

  (testing "weight falls with distance from the center"
    (let [{:keys [n dr dc w]} (e/brush 3)
          by-distance (sort-by first
                               (for [k (range n)]
                                 [(+ (* (aget ^ints dr k) (aget ^ints dr k))
                                     (* (aget ^ints dc k) (aget ^ints dc k)))
                                  (aget ^doubles w k)]))]
      (is (apply >= (map second by-distance))))))

(deftest flat-terrain-test
  ;; Every cut is clamped by the drop the droplet just fell, and on level
  ;; ground there is no drop. Nothing to erode, and nothing picked up to
  ;; deposit -- so this is exact, not approximate.
  (testing "level ground is left exactly level"
    (let [f (e/erode (flat 64 0.5) {:rng (seeded 1)})]
      (is (every? #(== 0.5 %) (heights f)))))

  (testing "and the droplets still record where they went"
    (let [f (e/erode (flat 64 0.5) {:rng (seeded 1)})]
      (is (pos? (second (e/drainage-bounds f)))))))

(deftest purity-test
  ;; The field holds a mutable array, so this is worth stating outright.
  (testing "erode leaves the field it was given alone"
    (let [g (terrain)
          before (heights g)]
      (dotimes [_ 3] (e/erode g {:rng (seeded 2)}))
      (is (= before (heights g)))))

  (testing "erode! is the one that writes, and reports the same field back"
    (let [g (terrain)
          before (heights g)
          after (e/erode! g {:rng (seeded 2)})]
      (is (identical? (:heights g) (:heights after)))
      (is (not= before (heights g))))))

(deftest determinism-test
  (testing "the same seed erodes the same way"
    (is (= (heights (e/erode (terrain) {:rng (seeded 42)}))
           (heights (e/erode (terrain) {:rng (seeded 42)})))))

  (testing "and a different one does not"
    (is (not= (heights (e/erode (terrain) {:rng (seeded 42)}))
              (heights (e/erode (terrain) {:rng (seeded 43)}))))))

(deftest finiteness-test
  (testing "nothing runs away, however hard it rains"
    ;; The border taper is what this is really checking. Without it the
    ;; cells beside the edge are mined by every droplet that leaves and
    ;; refilled by none, and the map grows a trench that deepens without
    ;; bound -- it reached fifty times the depth of the terrain here.
    (let [g (terrain 6 11)
          [lo0 hi0] (extremes g)
          relief (- hi0 lo0)
          rng (seeded 5)
          f (assoc g :heights (double-array (heights g)))]
      (dotimes [_ 12]
        (e/erode! f {:rng rng}))
      (let [[lo hi] (extremes f)]
        (is (every? #(and (not (Double/isNaN %)) (not (Double/isInfinite %))) (heights f)))
        (is (< (- lo0 relief) lo) (str "cut to " lo " under a floor of " lo0))
        (is (<= hi (+ hi0 relief)) (str "piled to " hi " over a ceiling of " hi0))))))

(deftest mass-test
  (testing "material is not invented"
    (let [g (terrain)
          f (e/erode g {:rng (seeded 3)})]
      (is (<= (total f) (+ (total g) 1e-9)))))

  (testing "and barely any leaves, because only the border exports it"
    ;; What a droplet picks up it must put down again, unless it walks off
    ;; the map. One droplet per cell over a 129-cell square should lose a
    ;; few percent, not a few fifths -- the latter is what happens when a
    ;; droplet that dries up is allowed to take its load with it.
    (let [g (terrain)
          f (e/erode g {:rng (seeded 3)})
          lost (/ (- (total g) (total f)) (total g))]
      (is (< 0.0 lost 0.10) (str "lost " lost " of the material"))))

  (testing "and what leaves, leaves at the border"
    ;; Sediment is exported by walking off the map and by no other route,
    ;; so how much goes is set by how many droplets get that far -- and a
    ;; droplet that may take only four steps mostly cannot. Loss falls
    ;; away with the lifetime, to a twentieth of a percent and then to
    ;; nothing much at all.
    (let [g (terrain 7 19)
          lost (fn [lifetime]
                 (/ (- (total g) (total (e/erode g {:rng (seeded 3)
                                                    :droplets 2000
                                                    :lifetime lifetime})))
                    (total g)))]
      (is (apply < (map lost [4 10 30])))
      (is (< (lost 4) 1e-3) (str "a four-step droplet still lost " (lost 4))))))

(deftest carving-test
  (testing "erosion cuts valleys into a slope rather than smoothing it"
    ;; A uniform ramp has no valleys. After rain it does: the droplets
    ;; converge, and the cells they converge on are cut below their
    ;; neighbors, so the spread of heights within a row goes from nothing
    ;; to something.
    (let [dim 96
          r (ramp dim 0.02)
          row (fn [{:keys [^doubles heights]} i]
                (mapv #(aget heights (+ (* i dim) %)) (range dim)))
          spread (fn [f i] (let [v (row f i)] (- (reduce max v) (reduce min v))))]
      (is (zero? (spread r 48)))
      (let [f (e/erode r {:rng (seeded 8) :droplets (* 4 dim dim)})]
        (is (pos? (spread f 48)) "the ramp is no longer flat across")))))

(deftest drainage-test
  (let [dim 96
        f (e/erode (ramp dim 0.02) {:rng (seeded 8) :droplets (* 2 dim dim)})
        row-flux (fn [i] (reduce + (map #(e/drainage f i %) (range dim))))]

    (testing "flux accumulates downstream"
      ;; The ramp falls toward the last row, so every droplet that passes
      ;; a high row must also pass the rows below it.
      (is (< (row-flux 10) (row-flux 40) (row-flux 70))))

    (testing "flux accumulates across calls, so erosion can be run in rounds"
      ;; `erode!` returns the field it was handed, so the total has to be
      ;; read out before the second round rather than after it.
      (let [sum-flux (fn [g] (reduce + (for [i (range dim) j (range dim)]
                                         (e/drainage g i j))))
            once (e/erode (ramp dim 0.02) {:rng (seeded 8) :droplets 500})
            after-one (sum-flux once)
            twice (e/erode! once {:rng (seeded 9) :droplets 500})]
        (is (pos? after-one))
        (is (> (sum-flux twice) after-one))))))

(deftest channel-test
  (testing "flux is heavily skewed, which is what a channel network is"
    ;; If water spread evenly there would be no rivers, only sheet wash.
    ;; A ramp is the wrong place to look -- every droplet on it runs
    ;; straight down and nothing converges -- so this asks real terrain.
    (let [g (terrain)
          dim (:dim g)
          f (e/erode g {:rng (seeded 8) :droplets (* 2 dim dim)})
          vs (vec (sort (for [i (range dim) j (range dim)] (e/drainage f i j))))
          n (count vs)
          mean (/ (reduce + vs) n)]
      (is (zero? (first (e/drainage-bounds f))))
      (is (> (peek vs) (* 4.0 mean))
          (str "busiest " (peek vs) " against a mean of " mean))
      (is (> (vs (int (* 0.99 n))) (* 2.0 (vs (quot n 2))))
          "the top percentile carries multiples of the median"))))

(deftest border-test
  (testing "the border is cut far less than the interior"
    ;; The taper zeroes erosion *centered* on an edge cell, not erosion
    ;; that reaches one: a droplet a few cells in still swings its disk
    ;; over the edge. So the outermost cells are worn lightly rather than
    ;; not at all, and it is the ratio that matters -- an untapered run
    ;; cuts the border harder than the interior, not a fifth as hard, and
    ;; then keeps going until it has dug a trench.
    (let [g (terrain 6 23)
          dim (:dim g)
          f (e/erode g {:rng (seeded 4) :droplets (* 4 dim dim)})
          cut (fn [i j] (abs (- (aget ^doubles (:heights g) (+ (* i dim) j))
                                (aget ^doubles (:heights f) (+ (* i dim) j)))))
          band (fn [rows] (/ (reduce + (for [i rows j (range dim)] (cut i j)))
                             (* (count rows) dim)))
          edge (band [0])
          interior (band (range (quot dim 3) (* 2 (quot dim 3))))]
      (is (pos? interior))
      (is (< edge (* 0.5 interior))
          (str "edge wore " edge " against an interior " interior)))))

(deftest shape-test
  (testing "a non-square field works if it says so"
    (let [f (e/erode {:heights (double-array (for [i (range 40) _ (range 20)] (* 0.02 (- 39 i))))
                      :rows 40 :cols 20}
                     {:rng (seeded 6) :droplets 500})]
      (is (= 800 (alength ^doubles (:heights f))))
      (is (every? #(and (not (Double/isNaN %)) (not (Double/isInfinite %))) (heights f)))))

  (testing "a field too small to interpolate across is refused"
    (is (thrown? clojure.lang.ExceptionInfo
                 (e/erode {:heights (double-array 1) :dim 1} {:rng (seeded 1)})))))
