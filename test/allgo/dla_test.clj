(ns allgo.dla-test
  (:require [allgo.procedural.dla :as dla]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded [seed]
  (let [r (java.util.Random. seed)]
    (fn [] (.nextDouble r))))

(defn- grow
  ([] (grow {}))
  ([opts] (dla/aggregate (merge {:dim 97 :density 0.05 :rng (seeded 3)} opts))))

(def ^:private cluster (delay (grow)))

(deftest tree-test
  (let [{:keys [order seed dim] :as c} @cluster
        ^ints parent (:parent c)
        occupied (set order)]

    (testing "it grew"
      (is (> (count order) 100))
      (is (= seed (first order))))

    (testing "every cell arrived once"
      (is (apply distinct? order)))

    (testing "the seed is the root and nothing else is"
      (is (= -1 (aget parent seed)))
      (is (every? #(not= -1 (aget parent %)) (remove #{seed} order))))

    (testing "every cell stuck to something already there"
      ;; What makes it a tree rather than a set of cells, and the whole
      ;; reason it can be measured for height.
      (is (every? #(contains? occupied (aget parent %)) (remove #{seed} order))))

    (testing "and stuck to a cell it was touching"
      ;; Four-connected: a walker that could slip diagonally between two
      ;; occupied cells would leave holes inside what should be trunk.
      (is (every? (fn [cell]
                    (let [[x y] (dla/xy c cell)
                          [px py] (dla/xy c (aget parent cell))]
                      (= 1 (+ (abs (- x px)) (abs (- y py))))))
                  (remove #{seed} order))))

    (testing "a parent always arrives before its children"
      ;; `subtree-sizes` walks `:order` backward and depends on it.
      (let [arrival (into {} (map-indexed (fn [i c] [c i])) order)]
        (is (every? #(< (arrival (aget parent %)) (arrival %))
                    (remove #{seed} order)))))

    (testing "every cell is on the grid"
      (is (every? #(< -1 % (* dim dim)) order)))))

(deftest sizes-test
  (let [{:keys [order seed] :as c} @cluster
        ^ints parent (:parent c)
        ^ints sizes (dla/subtree-sizes c)]

    (testing "the seed carries the whole cluster"
      (is (= (count order) (aget sizes seed))))

    (testing "a cell carries itself plus everything hanging off it"
      (let [children (group-by #(aget parent %) (remove #{seed} order))]
        (is (every? (fn [cell]
                      (= (aget sizes cell)
                         (inc (reduce + 0 (map #(aget sizes %) (get children cell []))))))
                    order))))

    (testing "tips carry one, and there are plenty of them"
      ;; The branching is the point: a cluster that grew as a blob would
      ;; have few tips.
      (let [tips (filter #(= 1 (aget sizes %)) order)]
        (is (> (count tips) (* 0.1 (count order))))))))

(deftest heightmap-test
  (let [c @cluster]

    (testing "a heightmap comes out the shape the rest of the repo takes"
      (let [{:keys [heights dim]} (dla/heightmap c {:refine 0})]
        (is (= (:dim c) dim))
        (is (= (* dim dim) (alength ^doubles heights)))))

    (testing "refining doubles the resolution each time"
      (is (= 193 (:dim (dla/heightmap c {:refine 1}))))
      (is (= 385 (:dim (dla/heightmap c {:refine 2})))))

    (testing "every height is a finite fraction"
      (let [{:keys [heights]} (dla/heightmap c {:refine 1})
            vs (map #(aget ^doubles heights %) (range (alength ^doubles heights)))]
        (is (every? #(<= 0.0 % 1.0) vs))
        (is (every? #(and (not (Double/isNaN %)) (not (Double/isInfinite %))) vs))))

    (testing "the trunk is higher than the tips"
      ;; The claim the whole height model rests on: height is read off the
      ;; tree, not off the picture, so a cell carrying half the cluster
      ;; stands above one carrying nothing. Checked before blurring,
      ;; which is what mixes the two together.
      (let [{:keys [heights dim]} (dla/heightmap c {:refine 0 :blur 0})
            ^ints sizes (dla/subtree-sizes c)
            at (fn [cell] (aget ^doubles heights (long cell)))
            trunk (filter #(> (aget sizes %) 20) (:order c))
            tips (filter #(= 1 (aget sizes %)) (:order c))
            mean (fn [xs] (/ (reduce + (map at xs)) (max 1 (count xs))))]
        (is (seq trunk))
        (is (= dim (:dim c)))
        (is (> (mean trunk) (* 2.0 (mean tips))))))))

(deftest options-test
  (testing "density sets how much of the grid fills up"
    (let [sparse (count (:order (grow {:density 0.02})))
          dense (count (:order (grow {:density 0.08})))]
      (is (< sparse dense))))

  (testing "a bigger grid grows a bigger cluster"
    (is (< (count (:order (grow {:dim 65})))
           (count (:order (grow {:dim 129})))))))

(deftest determinism-test
  (testing "the same stream grows the same cluster"
    (is (= (:order (grow)) (:order (grow)))))

  (testing "a different one does not"
    (is (not= (:order (grow {:rng (seeded 1)}))
              (:order (grow {:rng (seeded 2)}))))))
