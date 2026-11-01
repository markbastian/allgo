(ns allgo.probability.temporal-test
  (:require [allgo.probability.temporal :as tm]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The umbrella world of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 14, and the answers aima-python's
;; (MIT) tests expect. State 0 is rain.

(def ^:private u tm/umbrella)

(defn- brute
  "P(X_k | evidence) for each k by summing over every state sequence: an
  account of smoothing independent of the messages."
  [hmm ev]
  (let [t (count ev)
        seqs (reduce (fn [acc _] (for [a acc s [0 1]] (conj a s))) [[]] (range (inc t)))
        weight (fn [xs] (reduce * (get (:prior hmm) (xs 0))
                                (for [k (range 1 (inc t))]
                                  (* (get-in hmm [:transition (xs (dec k)) (xs k)]) (get-in hmm [:sensor (nth ev (dec k)) (xs k)])))))
        total (reduce + (map weight seqs))]
    (vec (for [k (range 1 (inc t))]
           (mapv (fn [s] (/ (reduce + (map weight (filter #(= s (% k)) seqs))) total)) [0 1])))))

(defn- close? [a b tol] (every? true? (map #(< (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(def ^:private evidences [[true true false true true] [true false true false true]])

(deftest smoothing-test
  (testing "forward-backward equals smoothing by brute force"
    (doseq [ev evidences] (is (close? (brute u ev) (tm/forward-backward u ev) 1e-12))))
  (testing "filtering: 0.818 after one umbrella, 0.883 after two"
    (let [[_ f1 f2] (tm/filtering-seq u [true true])]
      (is (close? 0.8182 (f1 0) 1e-4))
      (is (close? 0.8834 (f2 0) 1e-4))))
  (testing "fixed-lag smoothing equals smoothing d steps back, given everything so far"
    (doseq [ev evidences d [1 2 3]
            [t est] (map-indexed vector (tm/fixed-lag-smoothing-seq u ev d))]
      (if (< t d)
        (is (nil? est))
        (is (close? ((brute u (subvec ev 0 (inc t))) (- t d)) est 1e-9) [ev d t]))))
  (testing "aima-python's fixed-lag value: [0.2259 0.7741]"
    (is (close? [0.2259 0.7741] (nth (tm/fixed-lag-smoothing-seq u (second evidences) 2) 3) 1e-4))))

(deftest viterbi-test
  (testing "the paths and messages aima-python's test expects"
    (let [{:keys [path probabilities]} (tm/viterbi u (first evidences))]
      (is (= [0 0 1 0 0] path))
      (is (close? [0.8182 0.5155 0.1237 0.0334 0.0210] probabilities 1e-4)))
    (let [{:keys [path probabilities]} (tm/viterbi u (second evidences))]
      (is (= [0 1 1 1 0] path))
      (is (close? [0.8182 0.1964 0.0275 0.0154 0.0042] probabilities 1e-4)))))

(deftest particles-and-kalman-test
  (testing "5000 particles come near the exact filtered distribution"
    (let [ps (last (tm/particle-filtering-seq u [true true false] 5000 (random/rng 2)))]
      (is (close? (last (tm/filtering-seq u [true true false])) (tm/particle-estimate ps 2) 0.03))))
  (testing "the Kalman filter's closed form, from aima-python's test: mean 1.875, variance 0.75"
    (let [[[mean cov]] (tm/kalman-seq {:F [[1.0]] :Q [[2.0]] :H [[1.0]] :R [[1.0]]} [[0.0] [[1.0]]] [[2.5]])]
      (is (close? 1.875 (first mean) 1e-12))
      (is (close? 0.75 (ffirst cov) 1e-12))))
  (testing "observing only position, it infers a positive velocity"
    (let [[mean] (last (tm/kalman-seq {:F [[1.0 1.0] [0.0 1.0]] :Q [[0.01 0.0] [0.0 0.01]] :H [[1.0 0.0]] :R [[1.0]]}
                                      [[0.0 0.0] [[1.0 0.0] [0.0 1.0]]] (map vector [1.0 2.0 3.0 4.0 5.0])))]
      (is (pos? (second mean))))))
