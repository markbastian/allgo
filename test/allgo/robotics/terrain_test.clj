(ns allgo.robotics.terrain-test
  (:require [allgo.random :as random]
            [allgo.robotics.terrain :as terrain]
            [clojure.test :refer [deftest is testing]]))

;; Self-supervised terrain learning, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 26.4.3: a road
;; through grass with rocks in it.

(def ^:private colors {:road [0.5 0.5 0.52] :grass [0.25 0.6 0.2] :rock [0.55 0.4 0.25]})

(defn- scene [rng road-left road-right]
  (let [kind (fn [x y] (cond (and (<= road-left x) (< x road-right)) :road
                             (zero? (mod (+ (* 7 x) (* 13 y)) 23)) :rock
                             :else :grass))]
    {:image {:w 40 :h 30 :px (vec (for [y (range 30) x (range 40)] (mapv #(random/gaussian rng % 0.03) (colors (kind x y)))))}
     :kind kind :flat? (fn [x y] (not= :rock (kind x y)))}))

(deftest self-supervised-test
  (let [rng (random/rng 1)
        on-road (scene rng 12 28)
        ;; turned off the road: grass and rocks ahead
        off-road (scene rng 0 2)
        [first-model second-model] (map :model (terrain/adapt-seq [on-road off-road] [15 25 25 30] 2 20 2000 rng))
        test (scene rng 5 20)
        share (fn [model k] (let [ps (for [y (range 30) x (range 40) :when (= k ((:kind test) x y))] ((:px (:image test)) (+ (* y 40) x)))]
                              (/ (count (filter #(terrain/drivable? model %) ps)) (double (count ps)))))]
    (testing "taught by the laser on the road just ahead, the camera finds the road drivable everywhere, and not the grass"
      (is (> (share first-model :road) 0.95))
      (is (< (share first-model :grass) 0.05)))
    (testing "driven onto the grass, the laser finds it flat, and now the camera takes grass for drivable too"
      (is (> (share second-model :grass) 0.95))
      (is (> (share second-model :road) 0.95)))
    (testing "rocks the laser never finds flat, and the camera never takes them for drivable"
      (is (zero? (share first-model :rock)))
      (is (zero? (share second-model :rock))))))
