(ns allgo.vision.edges-test
  (:require [allgo.random :as random]
            [allgo.vision.edges :as edges]
            [allgo.vision.image :as img]
            [clojure.test :refer [deftest is testing]]))

;; Edges and texture, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 25, on drawn images whose edges
;; are known.

(def ^:private square
  (img/add-noise (img/render 48 48 0.2 [{:kind :rect :box [12 14 36 34] :value 0.8}]) 0.03 (random/rng 1)))

(defn- boundary-distance [[x y]]
  (let [cx (+ x 0.5) cy (+ y 0.5)]
    (min (abs (- cx 12)) (abs (- cx 36)) (abs (- cy 14)) (abs (- cy 34)))))

(deftest canny-test
  (let [stages (edges/canny-seq square)
        e (:edges (last stages))]
    (testing "the stages in order: smoothing, gradient, thinning, then hysteresis growing"
      (is (= [:smoothed :gradient :thinned] (map :stage (take 3 stages))))
      (is (apply <= (map #(count (:edges %)) (drop 3 stages)))))
    (testing "on a noisy square: nearly the whole perimeter, every pixel within half a pixel of the true boundary"
      (is (> (count e) 80))
      (is (every? #(<= (boundary-distance %) 0.5) e)))
    (testing "and one pixel thick: no two edge pixels side by side across the top"
      (is (not-any? #(and (e [% 13]) (e [% 14])) (range 48))))
    (testing "noise alone, without the square, makes no edges"
      (is (empty? (edges/canny (img/add-noise (img/image 48 48 (constantly 0.2)) 0.03 (random/rng 2))))))))

(deftest texture-test
  (testing "vertical stripes vote for one orientation, horizontal ones for the one across it"
    (let [stripes (fn [vertical?] (img/render 32 32 0.0 (for [i (range 0 32 4)]
                                                          {:kind :rect :value 1.0
                                                           :box (if vertical? [i 0 (+ i 2) 32] [0 i 32 (+ i 2)])})))
          h (fn [im] (get (edges/orientation-histograms im 16 4) [0 0]))]
      (is (= [1.0 0.0 0.0 0.0] (h (stripes true))))
      (is (= [0.0 0.0 1.0 0.0] (h (stripes false))))
      (is (= [0.0 0.0 0.0 0.0] (h (img/image 32 32 (constantly 0.5)))) "and a flat patch has none"))))
