(ns allgo.robotics.fastslam-test
  (:require [allgo.random :as random]
            [allgo.robotics.fastslam :as fs]
            [allgo.robotics.localization :as loc]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; FastSLAM, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, chapter 26, on the lap of `slam-test`: checked
;; against the truth the simulation knows.

(def ^:private lap
  (let [w (/ (/ math/PI 2) 1.6)]
    (concat (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w])
            (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w]))))

(defn- dist [[x y] [a b]] (math/sqrt (+ (* (- x a) (- x a)) (* (- y b) (- y b)))))

(deftest fastslam-test
  (let [world (loc/landmark-field)
        start [3.0 3.5 0.0]
        {:keys [truth observations]} (loc/drive-and-observe world start lap (random/rng 4))
        sets (vec (fs/fastslam-seq world (random/rng 7) (fs/initial-particles start 50) lap observations))
        final (peek sets)
        dead (reduce (fn [p [v w]] (loc/drive p v w 0.1)) start lap)]
    (testing "every landmark is mapped, near where it is"
      (is (= 12 (count (fs/landmarks (fs/best final)))))
      (doseq [[id m] (fs/landmarks (fs/best final))]
        (is (< (dist m ((:landmarks world) id)) 0.3) id)))
    (testing "the pose is tracked far better than dead reckoning manages"
      (is (< (apply max (map #(dist (fs/pose %1) %2) sets truth)) 0.3))
      (is (< (dist (fs/pose final) (peek truth)) 0.1))
      (is (> (dist dead (peek truth)) 0.2)))
    (testing "each particle keeps a small filter per landmark, not one Gaussian over the map"
      (let [l (get-in (fs/best final) [:landmarks 0])]
        (is (= [2 2] [(count (:cov l)) (count (first (:cov l)))]))))
    (testing "a landmark's covariance shrinks as it is seen again"
      (let [first-seen (ffirst (first (remove empty? observations)))
            k (count (take-while empty? observations))
            trace (fn [ps] (let [[[a _] [_ d]] (get-in (fs/best ps) [:landmarks first-seen :cov])] (+ a d)))]
        (is (< (trace (sets (+ k 40))) (trace (sets (inc k)))))))))
