(ns allgo.language.seq2seq-test
  (:require [allgo.language.seq2seq :as s2s]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Sequence to sequence with attention, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 24,
;; learning to reverse strings of digits -- whose attention should run
;; back along the source.

(deftest reversal-test
  (let [rng (random/rng 1)
        digits #(vec (repeatedly (+ 3 (random/below rng 3)) (fn [] (random/below rng 5))))
        pairs (vec (repeatedly 500 (fn [] (let [d (digits)] [d (vec (reverse d))]))))
        m (:model (nth (s2s/train-seq (s2s/model :rnn (s2s/vocabulary (range 5)) 8 16 (random/rng 2)) pairs :lr 0.01) 1500))
        tests (repeatedly 30 digits)]
    (testing "it reverses strings it has not seen, greedily and by beam search"
      (is (>= (count (filter #(= (vec (reverse %)) (:target (s2s/greedy m %))) tests)) 27))
      (is (>= (count (filter #(= (vec (reverse %)) (s2s/beam m % 3)) (take 10 tests))) 9)))
    (testing "and each output token attends most to the source token it copies: the anti-diagonal"
      (let [{:keys [attention]} (s2s/greedy m [0 1 2 3 4])]
        (is (= [4 3 2 1 0] (mapv #(apply max-key % (range 5)) (take 5 attention))))))))

(deftest basic-test
  (let [trained (fn [cell attention lo steps]
                  (let [rng (random/rng 1)
                        digits #(vec (repeatedly (+ lo (random/below rng 3)) (fn [] (random/below rng 5))))
                        pairs (vec (repeatedly 500 (fn [] (let [d (digits)] [d (vec (reverse d))]))))
                        tests (vec (repeatedly 30 digits))
                        m (:model (nth (s2s/train-seq (s2s/model cell (s2s/vocabulary (range 5)) 8 16 (random/rng 2) :attention attention) pairs :lr 0.01) steps))]
                    (count (filter #(= (vec (reverse %)) (:target (s2s/greedy m %))) tests))))]
    (testing "without attention, the decoder knows the source by the encoder's last state alone: enough for three to five digits"
      (is (>= (trained :lstm false 3 1000) 20)))
    (testing "but not for six to eight -- the fixed context size the chapter describes"
      (is (zero? (trained :lstm false 6 1000))))
    (testing "with attention, even the plain RNN, given half as long again, reverses them all"
      (is (= 30 (trained :rnn true 6 1500))))))
