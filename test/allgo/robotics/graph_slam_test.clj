(ns allgo.robotics.graph-slam-test
  (:require [allgo.random :as random]
            [allgo.robotics.graph-slam :as gs]
            [allgo.robotics.localization :as loc]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Graph SLAM, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 26, on the lap of `slam-test`:
;; checked against the truth the simulation knows.

(def ^:private lap
  (let [w (/ (/ math/PI 2) 1.6)]
    (concat (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w])
            (repeat 120 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 50 [1.0 0.0]) (repeat 16 [0.5 w]))))

(defn- dist [[x y] [a b]] (math/sqrt (+ (* (- x a) (- x a)) (* (- y b) (- y b)))))

(deftest graph-slam-test
  (let [world (loc/landmark-field)
        start [3.0 3.5 0.0]
        {:keys [truth observations]} (loc/drive-and-observe world start lap (random/rng 4))
        states (vec (gs/graph-slam-seq world start lap observations))
        {:keys [graph estimate]} (peek states)
        kept (:poses graph)
        path-error (fn [{:keys [poses]}] (apply max (map (fn [p k] (dist p (truth k))) poses kept)))]
    (testing "Gauss-Newton never raises the error, and converges in a few iterations"
      (is (apply >= (map :error states)))
      (is (< (count states) 15)))
    (testing "the weighted squared error ends near the measurements' degrees of freedom, as it should"
      (let [residuals (+ (* 3 (count (:odometry graph))) (* 2 (count (:sightings graph))))
            unknowns (+ (* 3 (dec (count kept))) (* 2 (count (:landmarks estimate))))
            dof (- residuals unknowns)]
        (is (< (* 0.6 dof) (:error (peek states)) (* 1.4 dof)))))
    (testing "every landmark is mapped, near where it is"
      (is (= 12 (count (:landmarks estimate))))
      (doseq [[id m] (:landmarks estimate)]
        (is (< (dist m ((:landmarks world) id)) 0.3) id)))
    (testing "the whole path is corrected, its start too, from dead reckoning's drift"
      (is (> (path-error (:estimate (first states))) 0.5))
      (is (< (path-error estimate) 0.25)))
    (testing "the first pose stays where it was put"
      (is (= start (first (:poses estimate)))))))
