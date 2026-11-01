(ns allgo.language.recurrent-test
  (:require [allgo.language.parsing :as parsing]
            [allgo.language.recurrent :as rec]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Recurrent language models, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 24, on sentences
;; drawn from a grammar -- whose parser then judges what they write.

(deftest language-model-test
  (let [rng (random/rng 1)
        sents (vec (repeatedly 400 #(mapv first (parsing/generate parsing/ambiguous-grammar rng))))
        [train test] [(subvec sents 0 300) (subvec sents 300)]
        m (:model (nth (rec/train-seq (rec/model :rnn (rec/vocabulary sents) 8 16 (random/rng 2)) train :lr 0.01) 600))]
    (testing "far less surprised by unseen sentences than a model of word frequencies alone"
      (is (< (rec/perplexity m test) (* 0.5 (rec/unigram-perplexity train test)))))
    (testing "and what it writes is mostly grammatical: its sentences parse"
      (is (>= (count (filter #(parsing/parse-tree % parsing/ambiguous-grammar) (repeatedly 50 #(rec/sample m rng)))) 35)))
    (testing "a low temperature writes the likeliest sentence, a high one more of them"
      (let [distinct-at (fn [t] (count (distinct (repeatedly 20 #(rec/sample m rng :temperature t)))))]
        (is (< (distinct-at 0.05) (distinct-at 1.5)))))))

(deftest long-dependency-test
  (testing "the last word must match the first, six words back: the LSTM learns it, the plain RNN does not"
    (let [rng (random/rng 2)
          data (vec (repeatedly 200 (fn [] (let [f (random/pick rng ["a" "b"])] (vec (concat [f] (repeat 6 "x") [f]))))))
          voc (rec/vocabulary data)
          ;; the probability put on the right word of the two, at the end
          score (fn [m] (/ (reduce + (for [f ["a" "b"]]
                                       (let [d (rec/next-distribution m (vec (concat [f] (repeat 6 "x"))))
                                             o ({"a" "b" "b" "a"} f)]
                                         (/ (d f) (+ (d f) (d o))))))
                           2.0))
          trained (fn [cell] (:model (nth (rec/train-seq (rec/model cell voc 4 8 (random/rng 2)) data :lr 0.05) 800)))]
      (is (> (score (trained :lstm)) 0.9))
      (is (< (score (trained :rnn)) 0.6)))))
