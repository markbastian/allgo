(ns allgo.learning.unsupervised-test
  (:require [allgo.learning.unsupervised :as un]
            [allgo.probability.temporal :as tm]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Learning with hidden variables, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 20.

(def ^:private points
  (let [rng (random/rng 4)
        blob (fn [cx cy] (vec (repeatedly 100 (fn [] [(+ cx (* 0.6 (- (rng) 0.5))) (+ cy (* 0.6 (- (rng) 0.5)))]))))]
    (vec (concat (blob 0 0) (blob 3 0) (blob 0 3)))))

(deftest clustering-test
  (testing "k-means: the cost never rises, and the centers find the three blobs"
    (let [states (take 12 (rest (un/k-means-seq points 3 (random/rng 4))))]
      (is (apply >= (map :cost states)))
      (is (= #{[0 0] [3 0] [0 3]} (set (map (fn [c] (mapv #(Math/round (double %)) c)) (:centers (last states))))))))
  (testing "EM: the log likelihood never falls, and the weights come out a third each"
    (let [states (take 30 (rest (un/em-gmm-seq points 3 (random/rng 4))))]
      (is (every? true? (map #(>= %2 (- %1 1e-9)) (map :log-likelihood states) (rest (map :log-likelihood states)))))
      (is (every? #(< (abs (- % (/ 1.0 3))) 0.01) (:weights (last states)))))))

(deftest baum-welch-test
  (testing "Baum-Welch on umbrella observations: the likelihood never falls, and the models stay distributions"
    (let [obs [true true false true true false false false true false true true]
          states (take 30 (un/baum-welch-seq tm/umbrella obs))
          lls (map #(un/log-likelihood % obs) states)]
      (is (every? true? (map #(>= %2 (- %1 1e-9)) lls (rest lls))))
      (is (> (last lls) (first lls)))
      (doseq [row (:transition (last states))] (is (< (abs (- 1.0 (reduce + row))) 1e-9))))))
