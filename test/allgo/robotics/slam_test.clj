(ns allgo.robotics.slam-test
  (:require [allgo.random :as random]
            [allgo.robotics.localization :as loc]
            [allgo.robotics.slam :as slam]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; EKF-SLAM, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, chapter 26, on a simulated lap of a field of
;; landmarks: checked against the truth the simulation knows.

(def ^:private lap
  (let [w (/ (/ math/PI 2) 1.6)]
    (concat (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w])
            (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w]))))

(defn- dist [[x y] [a b]] (Math/hypot (- x a) (- y b)))

(deftest slam-test
  (let [world (loc/landmark-field)
        start [3.0 3.5 0.0]
        {:keys [truth observations]} (loc/drive-and-observe world start lap (random/rng 4))
        beliefs (vec (slam/slam-seq world (slam/initial-belief start 12) lap observations))
        final (peek beliefs)
        dead (reduce (fn [p [v w]] (loc/drive p v w 0.1)) start lap)]
    (testing "every landmark is mapped, near where it is"
      (is (= 12 (count (slam/landmarks final))))
      (doseq [[id m] (slam/landmarks final)]
        (is (< (dist m ((:landmarks world) id)) 0.3) id)))
    (testing "the pose is tracked far better than dead reckoning manages"
      (is (< (apply max (map #(dist (slam/pose %1) %2) beliefs truth)) 0.3))
      (is (< (dist (slam/pose final) (peek truth)) 0.1))
      (is (> (dist dead (peek truth)) 0.2)))
    (testing "closing the loop shrinks the uncertainty of the landmarks first seen"
      (let [first-seen (ffirst (first (remove empty? observations)))
            var-of (fn [b] (let [[[a _] [_ d]] (slam/landmark-cov b first-seen)] (+ a d)))
            far (apply max-key #(var-of (beliefs %)) (filter #(contains? (:seen (beliefs %)) first-seen) (range (count beliefs))))]
        (is (< (var-of final) (var-of (beliefs far))))))
    (testing "the covariance stays symmetric"
      (let [P (:cov final)]
        (is (every? #(< (abs (- (get-in P %) (get-in P (reverse %)))) 1e-9)
                    (for [i (range (count P)) j (range i)] [i j])))))))
