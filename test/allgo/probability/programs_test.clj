(ns allgo.probability.programs-test
  (:require [allgo.probability.programs :as pp]
            [allgo.probability.temporal :as tm]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Probabilistic programs, after the idea of Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 15;
;; the examples are this library's, each checked against an exact answer.

(def ^:private letters
  {:letters [:a :b :c] :start {:a 0.5 :b 0.3 :c 0.2}
   :next {:a {:a 0.1 :b 0.6 :c 0.3} :b {:a 0.5 :b 0.2 :c 0.3} :c {:a 0.4 :b 0.4 :c 0.2}}
   :accuracy 0.8})

(deftest markov-letters-test
  (testing "likelihood weighting on the program agrees with HMM smoothing on the same model"
    (let [seen [:a :c :b :b :a]
          T [[0.1 0.6 0.3] [0.5 0.2 0.3] [0.4 0.4 0.2]]
          ;; As an HMM: a state before the first letter whose one transition
          ;; is the start distribution.
          hmm {:prior [0.0 0.0 0.0 1.0]
               :transition (conj (mapv #(conj % 0.0) T) [0.5 0.3 0.2 0.0])
               :sensor (into {} (for [[l i] {:a 0 :b 1 :c 2}] [l (conj (mapv #(if (= % i) 0.8 0.1) (range 3)) 0.0)]))}
          exact ((tm/forward-backward hmm seen) 2)
          lw (:estimate (nth (pp/likelihood-weighting-seq (pp/markov-letters letters seen 2) (random/rng 3)) 40000))]
      (doseq [[l i] {:a 0 :b 1 :c 2}]
        (is (< (abs (- (exact i) (get lw l 0.0))) 0.02) l)))))

(deftest aircraft-test
  (let [m {:lambda 2.0 :n-max 12 :detect 0.7 :clutter 1.0}]
    (testing "the posterior on how many aircraft there are matches summing over the detections, however many blips"
      (doseq [k [0 2 5]
              :let [lw (:estimate (nth (pp/likelihood-weighting-seq (pp/aircraft-count m k) (random/rng 4)) 40000))
                    ex (pp/aircraft-count-exact m k)]
              n (range 13)]
        (is (< (abs (- (get lw n 0.0) (ex n))) 0.01) [k n])))
    (testing "more blips, more aircraft believed in"
      (let [mean (fn [k] (reduce + (map (fn [[n p]] (* n p)) (pp/aircraft-count-exact m k))))]
        (is (< (mean 0) (mean 2) (mean 5)))))))
