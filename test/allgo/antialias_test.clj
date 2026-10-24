(ns allgo.antialias-test
  (:require [allgo.procedural.antialias :as aa]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as shaping]
            [clojure.test :refer [deftest is testing]]))

(defn- brute
  "The box average of `f` over the `w` square, the slow way."
  [f x y w]
  (aa/supersample f x y w 200 :grid [0 0]))

(def ^:private footprints
  [[0.1 0.2 0.3] [0.37 0.81 0.9] [1.3 -2.2 2.5] [0.49 0.49 0.05] [-3.7 5.1 7.3]])

(deftest integration-test
  (testing "a filtered stripe is the fraction of the footprint that is on"
    (doseq [[x y w] footprints]
      (is (< (abs (- (aa/filtered-pulse-train 0.5 x w)
                     (brute (fn [x _] (shaping/pulse-train 0.5 1.0 x)) x y w)))
             0.01))))
  (testing "a filtered checker is the board's average over the footprint"
    (doseq [[x y w] footprints]
      (is (< (abs (- (aa/filtered-checker x y w) (brute aa/checker x y w))) 0.01))))
  (testing "a footprint of nothing is the point"
    (doseq [[x y] [[0.1 0.2] [0.7 0.3] [-1.2 0.6]]]
      (is (= (aa/checker x y) (aa/filtered-checker x y 0.0)))))
  (testing "seen from far enough, a pattern is its average"
    (is (< (abs (- 0.5 (aa/filtered-checker 0.3 0.6 1000.3))) 0.001))
    (is (< (abs (- 0.25 (aa/filtered-pulse-train 0.25 0.3 1000.7))) 0.001))))

(deftest clamped-fbm-test
  (let [f (aa/clamped-fbm (noise/gradient-basis {:seed 3}) {:octaves 8})
        spread (fn [w] (let [vs (for [i (range 400)] (f (* 0.37 i) (* 0.11 i) 0.5 w))]
                         (- (apply max vs) (apply min vs))))]
    (testing "a footprint wider than every octave leaves nothing to alias"
      (is (= 0.0 (f 1.3 2.7 0.5 100.0))))
    (testing "the wider the footprint, the less variation is left"
      (is (> (spread 0.0) (spread 0.3) (spread 1.0))))))

(deftest sampling-test
  (testing "supersampling converges on the average, and jitter gets there too"
    (let [[x y w] [0.37 0.81 0.9]
          exact (aa/filtered-checker x y w)]
      (is (< (abs (- exact (aa/supersample aa/checker x y w 16 :grid [0 0]))) 0.02))
      (is (< (abs (- exact (aa/supersample aa/checker x y w 16 :jittered [4 7]))) 0.03))))
  (testing "jitter is repeatable for a pixel, and differs between pixels"
    (let [s #(aa/supersample aa/checker 0.37 0.81 0.9 3 :jittered %)]
      (is (= (s [1 2]) (s [1 2])))
      (is (not= (s [1 2]) (s [2 1])))))
  (testing "colors average channel by channel"
    (is (= [0.5 0.25] (aa/supersample (fn [x _] (if (< x 0.0) [1.0 0.5] [0.0 0.0])) 0.0 0.0 1.0 2 :grid [0 0]))))
  (testing "adaptive sampling spends little on flat color and more on an edge"
    (let [[flat n-flat] (aa/adaptive aa/checker 0.25 0.25 0.1 0.1 3)
          [edge n-edge] (aa/adaptive aa/checker 0.5 0.25 0.1 0.1 3)]
      (is (= (aa/checker 0.25 0.25) flat))
      (is (= 5 n-flat))
      (is (> n-edge 5))
      (is (< (abs (- 0.5 edge)) 0.1)))))
