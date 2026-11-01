(ns allgo.robotics.localization-test
  (:require [allgo.random :as random]
            [allgo.robotics.localization :as loc]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Monte Carlo localization, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 26.

(deftest resample-test
  (testing "low-variance resampling keeps a particle of weight w about w n times"
    (let [ps (vec (range 10))
          ws [0.0 0.0 0.5 0.0 0.25 0.25 0.0 0.0 0.0 0.0]
          counts (frequencies (loc/resample ps ws (random/rng 1)))]
      (is (= 10 (reduce + (vals counts))))
      (is (= #{2 4 5} (set (keys counts))))
      (is (= 5 (counts 2)))
      (is (every? #(<= 2 (counts %) 3) [4 5])))))

(deftest raycast-test
  (testing "rays meet the walls where geometry says"
    (let [{:keys [walls max-range]} (loc/room)]
      (is (< (abs (- 5.5 (loc/raycast walls [1.5 4.5] 0.0 max-range))) 1e-9) "east to the pillar's west face, at x = 7")
      (is (< (abs (- 8.5 (loc/raycast walls [1.5 5.0] 0.0 max-range))) 1e-9) "and above it, on to the wall at x = 10")
      (is (< (abs (- 1.5 (loc/raycast walls [1.5 4.5] Math/PI max-range))) 1e-9) "west to x = 0")
      (is (< (abs (- 0.5 (loc/raycast walls [7.5 3.0] (/ Math/PI 2) max-range))) 1e-9) "north to the pillar"))))

(deftest localization-test
  (testing "from not knowing where it is, the robot's particles find it: within a meter at the end of the route, every time, and usually within half"
    (let [errors (for [seed [5 6 7 8 9 10]]
                   (let [world (assoc (loc/room) :rng (random/rng seed))
                         controls (concat (repeat 40 [1.0 0.0]) (repeat 16 [0.5 -1.0]) (repeat 15 [1.0 0.0]))
                         {:keys [truth particles]} (last (loc/mcl-seq world [1.5 4.5 0.0] controls 1500))
                         [[x y]] (loc/estimate particles)]
                     (Math/hypot (- x (truth 0)) (- y (truth 1)))))]
      (is (every? #(< % 1.0) errors))
      (is (< (nth (sort errors) 3) 0.5)))))

(def ^:private loop-controls
  ;; A lap of the landmark field: east, up, west, down, quarter turns.
  (let [w (/ (/ math/PI 2) 1.6)]
    (concat (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w])
            (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w]))))

(deftest sensor-model-test
  (testing "range and bearing, and their partials against central differences"
    (let [pose [1.0 2.0 0.3] m [4.0 6.0]
          [r b] (loc/range-bearing pose m)
          J (loc/range-bearing-jacobian pose m)
          h 1e-6
          numeric (fn [k] (let [bump (fn [d] (let [v (vec (concat pose m))
                                                   v (update v k + d)]
                                               (loc/range-bearing (subvec v 0 3) (subvec v 3 5))))
                                [r+ b+] (bump h) [r- b-] (bump (- h))]
                            [(/ (- r+ r-) (* 2 h)) (/ (- b+ b-) (* 2 h))]))]
      (is (< (abs (- r 5.0)) 1e-12))
      (is (< (abs (- b (- (math/atan2 4.0 3.0) 0.3))) 1e-12))
      (doseq [k (range 5) :let [[dr db] (numeric k)]]
        (is (< (abs (- dr (get-in J [0 k]))) 1e-6) k)
        (is (< (abs (- db (get-in J [1 k]))) 1e-6) k)))))

(deftest ekf-localization-test
  (testing "the EKF tracks the robot round a lap to within centimeters, its error consistent with its covariance"
    (doseq [seed [1 2 3]]
      (let [world (loc/landmark-field)
            start [3.0 3.5 0.0]
            {:keys [truth observations]} (loc/drive-and-observe world start loop-controls (random/rng seed))
            beliefs (vec (loc/ekf-localization-seq world [start [[0.01 0.0 0.0] [0.0 0.01 0.0] [0.0 0.0 0.001]]]
                                                   loop-controls observations))
            errors (map (fn [[m] t] (Math/hypot (- (m 0) (t 0)) (- (m 1) (t 1)))) beliefs truth)
            ;; the normalized estimation error squared of the position: 2 on
            ;; average, for a filter whose covariance is honest
            nees (map (fn [[m P] t]
                        (let [e0 (- (m 0) (t 0)) e1 (- (m 1) (t 1))
                              a (get-in P [0 0]) b (get-in P [0 1]) d (get-in P [1 1])]
                          (/ (+ (* d e0 e0) (* -2 b e0 e1) (* a e1 e1)) (- (* a d) (* b b)))))
                      (rest beliefs) (rest truth))
            dead (reduce (fn [p [v w]] (loc/drive p v w 0.1)) start loop-controls)]
        (is (< (apply max errors) 0.15) seed)
        (is (< 1.0 (/ (reduce + nees) (count nees)) 3.5) seed)
        (is (> (Math/hypot (- (dead 0) ((peek truth) 0)) (- (dead 1) ((peek truth) 1))) 0.2)
            "where dead reckoning alone drifts")))))
