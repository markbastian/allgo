(ns allgo.vision.segmentation-test
  (:require [allgo.graph :as graph]
            [allgo.random :as random]
            [allgo.vision.image :as img]
            [allgo.vision.segmentation :as seg]
            [clojure.test :refer [deftest is testing]]))

;; Segmentation, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 25, on a drawn scene of three
;; colored regions with noise.

(def ^:private scene
  (img/render 48 36 [0.2 0.3 0.6] [{:kind :rect :box [4 4 22 30] :value [0.8 0.3 0.2]}
                                   {:kind :disk :center [34 18] :radius 9 :value [0.3 0.8 0.3]}]))
(def ^:private noisy (img/add-noise scene 0.05 (random/rng 1)))

(defn- purity
  "The fraction of pixels whose region's most common true color is theirs."
  [labels]
  (let [groups (group-by first (map vector (:px labels) (:px scene)))]
    (/ (double (reduce + (map (fn [[_ ps]] (apply max (vals (frequencies (map second ps))))) groups)))
       (count (:px scene)))))

(deftest graph-based-test
  (testing "the regions follow the scene's: few, and nearly pure"
    (let [labels (seg/graph-based (img/smooth noisy 0.8) 2.0 20)]
      (is (<= 3 (seg/region-count labels) 6))
      (is (> (purity labels) 0.95))))
  (testing "a step an edge; regions only ever merge"
    (let [counts (map (fn [st] (count (distinct (map #(graph/canonical (:forest st) %) (keys (:size st))))))
                      (take-nth 500 (seg/graph-based-seq (img/smooth noisy 0.8) 2.0)))]
      (is (apply >= counts)))))

(deftest k-means-test
  (testing "k-means on the colors, with restarts, gives the three regions exactly"
    (let [labels (seg/k-means noisy 3 (random/rng 2))]
      (is (= 3 (seg/region-count labels)))
      (is (> (purity labels) 0.999)))))

(deftest normalized-cuts-test
  (doseq [seed (range 1 7)]
    (let [rng (random/rng seed)
          shapes [{:kind :rect :box [(+ 2 (random/below rng 4)) (+ 2 (random/below rng 4)) (+ 14 (random/below rng 4)) (+ 20 (random/below rng 5))] :value 0.75}
                  {:kind :disk :center [(+ 20 (random/below rng 4)) (+ 12 (random/below rng 6))] :radius (+ 5 (random/below rng 2)) :value 0.45}]
          scene (img/add-noise (img/render 30 28 0.15 shapes) 0.05 rng)
          truth (img/render 30 28 0 (map-indexed (fn [i s] (assoc s :value (inc i))) shapes))
          graph (seg/superpixel-graph scene 0.05 8 0.1)
          labels (seg/normalized-cuts scene {})]
      (testing "the over-segmentation leaves a few dozen superpixels where there were 840 pixels"
        (is (<= 20 (count (:nodes graph)) 50)))
      (testing "and recursive normalized cuts on them find the two shapes and the background, each region all one of them"
        (is (= 3 (seg/region-count labels)))
        (is (every? #(= 1 (count (set %)))
                    (vals (group-by first (map vector (:px labels) (:px truth))))))))))
