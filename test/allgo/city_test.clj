(ns allgo.city-test
  (:require [allgo.geometry.polygon :as poly]
            [allgo.procedural.city.agents :as agents]
            [allgo.procedural.city.buildings :as buildings]
            [allgo.procedural.city.facade :as facade]
            [allgo.procedural.city.land :as land]
            [allgo.procedural.city.lots :as lots]
            [allgo.procedural.city.roads :as roads]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(def ^:private the-land (land/land {:seed 2}))
(def ^:private the-roads (roads/grow the-land))
(def ^:private the-graph (roads/graph the-roads))
(def ^:private the-blocks (lots/blocks the-graph))
(def ^:private the-lots (lots/lots the-blocks (:zone the-land)))

(deftest land-test
  (let [{:keys [population water? center pattern size]} the-land]
    (testing "nobody lives on the water, and the main center is dry and busy"
      (is (not (water? (first center) (second center))))
      (is (> (population (first center) (second center)) 0.6))
      (doseq [x (range -1500 1500 97) y (range -1500 1500 89)]
        (when (water? x y) (is (zero? (population x y))))))
    (testing "the mixed pattern is radial at the center and a grid or organic further out"
      (is (> (:radial (pattern (first center) (second center))) 0.99))
      (is (< (:radial (pattern (+ (first center) size) (second center))) 0.01)))
    (testing "a pattern asked for is the only one"
      (is (= {:grid 1.0} ((:pattern (land/land {:pattern :grid})) 0.0 0.0))))
    (testing "the same seed is the same place"
      (is (= (:center (land/land {:seed 2})) center)))))

(deftest roads-test
  (let [segs (roads/segments the-roads)]
    (testing "it grows a city, highways and streets, and stops"
      (is (roads/done? the-roads))
      (is (> (count segs) 500))
      (is (= #{:highway :street} (set (map :kind segs)))))
    (testing "no road is built on the water, except as a bridge"
      (doseq [{:keys [a b bridge?]} segs :when (not bridge?)
              t [0.0 0.25 0.5 0.75 1.0]]
        (is (land/legal? the-land (poly/lerp a b t)))))
    (testing "self-sensitivity: no two roads cross except at a shared junction"
      (let [v (vec segs)]
        (doseq [i (range (count v)) j (range (inc i) (count v))
                :let [{a :a b :b} (v i) {c :a d :b} (v j)]
                :when (empty? (filter #{a b} [c d]))
                :let [[t u] (poly/segment-intersection a b c d)]]
          (is (not (and t (< 1e-6 t (- 1 1e-6)) (< 1e-6 u (- 1 1e-6))))
              [i j]))))
    (testing "and no two leave a junction at too close an angle"
      (let [{:keys [points edges]} the-graph
            out (reduce (fn [m [i j]] (-> m (update i conj j) (update j conj i))) {} edges)]
        (doseq [[n nbrs] out
                :let [p (points n)
                      ds (mapv #(poly/unit (poly/sub (points %) p)) (distinct nbrs))]
                i (range (count ds)) j (range (inc i) (count ds))]
          (is (< (poly/dot (ds i) (ds j)) (Math/cos 0.35))))))
    (testing "deterministic"
      (is (= segs (roads/segments (roads/grow the-land)))))))

(deftest lots-test
  (testing "faces of a square of four roads around a square of four more"
    (let [g {:points [[0.0 0.0] [100.0 0.0] [100.0 100.0] [0.0 100.0] [200.0 0.0] [200.0 100.0] [50.0 -40.0]]
             :edges [[0 1 :street] [1 2 :street] [2 3 :street] [3 0 :street]
                     [1 4 :street] [4 5 :street] [5 2 :highway]
                     ;; A dead end, which bounds nothing.
                     [0 6 :street]]}
          fs (lots/faces g)]
      (is (= 2 (count fs)))
      (is (every? #(< (Math/abs (- 10000.0 (poly/area (:points %)))) 1e-9) fs))
      (is (some #(some #{:highway} (:kinds %)) fs))))
  (testing "the blocks are inside their faces, smaller by the roads around them"
    (is (> (count the-blocks) 100))
    (is (every? #(pos? (poly/signed-area (:points %))) the-blocks)))
  (testing "every lot is convex, faces a street, and is no larger than its zone allows"
    (is (> (count the-lots) 1000))
    (doseq [{:keys [points street zone area]} the-lots]
      (is (poly/convex? points))
      (is (seq street))
      (is (<= area (get-in lots/defaults [:lot-area zone]))))
    (is (= #{:downtown :commercial :residential} (set (map :zone the-lots)))))
  (testing "subdivision covers the block exactly, before the lots are culled"
    (let [sq [[0.0 0.0] [60.0 0.0] [60.0 40.0] [0.0 40.0]]
          pieces (lots/subdivide sq [:s :s :s :s] 300.0 (random/rng 3) {})]
      (is (< (Math/abs (- 2400.0 (reduce + (map #(poly/area (:points %)) pieces)))) 1e-6))
      (is (every? #(<= (poly/area (:points %)) 300.0) pieces)))))

(deftest buildings-test
  (let [stages (mapv #(buildings/rounds % (:max-height the-land)) the-lots)]
    (testing "the axiom is the bounding box, and each round stays inside it"
      (doseq [s (take 300 stages)
              :let [[{:keys [footprint h]}] (first s)
                    top (fn [m] (+ (:z m) (or (:h m) 0.0)))]]
        (is (= 1 (count (first s))))
        (doseq [string (rest s) m string]
          (if (= :antenna (:op m))
            (is (poly/contains-point? footprint (:at m)))
            (do (is (<= (top m) (+ h 1e-6)))
                (is (every? #(poly/contains-point? (poly/inset footprint -0.01) %)
                            (or (:footprint m) (:rect m)))))))))
    (testing "and the last round has nothing left to rewrite"
      (is (every? #(not-any? :sym (peek %)) stages)))
    (testing "towers downtown, houses with roofs in the suburbs"
      (let [by-zone (group-by :zone the-lots)
            styles (fn [zone] (set (mapcat (fn [l] (keep :style (buildings/building l (:max-height the-land))))
                                           (by-zone zone))))]
        (is (contains? (styles :downtown) :tower))
        (is (contains? (styles :residential) :house))
        (is (some #(= :roof (:op %)) (mapcat #(buildings/building % (:max-height the-land)) (by-zone :residential))))))
    (testing "the same lot gives the same building"
      (let [l (first the-lots)]
        (is (= (buildings/building l (:max-height the-land)) (buildings/building l (:max-height the-land))))))))

(deftest facade-test
  (testing "interval groups: a grid, scaled, and the pulse that finds a cell"
    (let [g (facade/scale (facade/grid 3 9.0) 0.5)]
      (is (= [[0.75 2.25] [3.75 5.25] [6.75 8.25]] g))
      (is (= [1 0.5] (facade/locate g 4.5)))
      (is (nil? (facade/locate g 3.0)))))
  (testing "a layer is AND or OR of its two pulses, and a stack falls through to the wall"
    (let [s (facade/grid 2 2.0) t (facade/grid 2 2.0)
          window (facade/layer (facade/scale s 0.5) (facade/scale t 0.5) :and (constantly :glass))
          band (facade/layer [[0.0 0.1]] [[0.0 0.1]] :or (constantly :band))
          wall (constantly :wall)]
      (is (= :glass (facade/evaluate [window band] wall 0.5 0.5)))
      (is (= :band (facade/evaluate [window band] wall 0.05 1.0)))
      (is (= :wall (facade/evaluate [window band] wall 1.0 1.0)))))
  (testing "functions between layers: a house's lintels sit only over its openings"
    (let [p (facade/painter :house 10.0 6.0 3.0 4)]
      (is (every? #(= 3 (count %)) (for [s (range 0.1 10.0 0.7) t (range 0.1 6.0 0.4)] (p s t))))))
  (testing "every style paints every point a color"
    (doseq [style (keys facade/styles)]
      (let [p (facade/painter style 12.0 20.0 3.5 1)]
        (doseq [s (range 0.05 12.0 1.3) t (range 0.05 20.0 1.7)]
          (is (every? #(<= 0.0 % 1.0) (p s t)) style))))))

(deftest agents-test
  (let [w (agents/run 800 {:seed 4})
        n (* (:w w) (:w w))
        built (count (filter pos? (vec (:use w))))]
    (testing "a town grows: roads, and both kinds of building"
      (is (> (count (filter #{1} (vec (:road w)))) 100))
      (is (> built 300))
      (is (= #{0 1 2} (set (vec (:use w))))))
    (testing "nothing is built on a road or in the water"
      (doseq [k (range n)]
        (when (pos? (aget ^ints (:use w) k))
          (is (zero? (aget ^ints (:road w) k)))
          (is (zero? (aget ^ints (:water w) k))))))
    (testing "commerce sits where the roads are thicker than where people live"
      (let [wd (:w w)
            mean-rd (fn [u] (let [ks (filter #(= u (aget ^ints (:use w) %)) (range n))]
                              (/ (reduce + (map #(agents/road-density w (mod % wd) (quot % wd)) ks)) (count ks))))]
        (is (> (mean-rd agents/commercial) (mean-rd agents/residential)))))
    (testing "with no deviation allowed, every road lies along a grid line"
      (let [g (agents/run 600 {:seed 4 :deviation 0})
            sp (:spacing agents/defaults)]
        (doseq [path @(:roads g)]
          (is (or (apply = (map first path)) (apply = (map second path))))
          (is (or (every? #(zero? (mod (first %) sp)) path) (every? #(zero? (mod (second %) sp)) path))))))))
