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
