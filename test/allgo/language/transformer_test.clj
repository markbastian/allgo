(ns allgo.language.transformer-test
  (:require [allgo.language.parsing :as parsing]
            [allgo.language.seq2seq :as s2s]
            [allgo.language.transformer :as tf]
            [allgo.learning.autodiff :as ad]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The transformer, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 24: tagging parts of speech in
;; sentences drawn from a grammar whose words are often more than one --
;; "the duck will duck" -- so that the tag needs the context.

(def ^:private data
  (delay (let [rng (random/rng 1)
               tagged (vec (repeatedly 500 #(parsing/generate parsing/ambiguous-grammar rng)))]
           {:train (subvec tagged 0 400) :test (subvec tagged 400)
            :vocab (tf/vocabulary (map #(mapv first %) tagged))})))

(def ^:private ambiguous #{"fish" "duck" "watch" "saw" "can"})
(def ^:private tags [:Det :N :V :Adj :Modal :Adv])

(deftest positional-encoding-test
  (testing "sines and cosines: position 0 is (0 1 0 1 ...), and each pair has length 1"
    (let [pe (tf/positional-encoding 5 8)]
      (is (= [0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0] (first pe)))
      (is (every? #(< (abs (- 1.0 (+ (* (% 0) (% 0)) (* (% 1) (% 1))))) 1e-12) pe)))))

(deftest tagging-test
  (let [{:keys [train test vocab]} @data
        m (:model (nth (tf/train-tagger-seq (tf/tagger (tf/encoder vocab 16 2 32 1 (random/rng 2)) tags (random/rng 3)) train) 800))
        base (tf/most-frequent-tag train)
        base-acc (fn [only?] (let [rs (for [s test [w t] s :when (only? w)] (= t (base w)))]
                               (/ (double (count (filter true? rs))) (count rs))))]
    (testing "the most frequent tag is often wrong on the ambiguous words"
      (is (< (base-acc ambiguous) 0.7)))
    (testing "the transformer, reading the context, gets nearly all of them"
      (is (> (tf/tag-accuracy m test ambiguous) 0.95))
      (is (> (tf/tag-accuracy m test) 0.97)))
    (testing "and its attention is a distribution over the sentence, per head"
      (let [{:keys [attention]} (tf/tag m ["the" "duck" "will" "duck" "a" "dog"])]
        (is (= 2 (count (first attention))))
        (is (every? #(< (abs (- 1.0 (reduce + %))) 1e-9) (apply concat (first attention))))))))

(deftest masked-lm-test
  (testing "pretrained on unlabeled sentences, the masked model fills each hole with a word of the right part of speech, mostly"
    (let [{:keys [train test vocab]} @data
          sents (map #(mapv first %) train)
          m (:model (nth (tf/train-mlm-seq (tf/masked-lm (tf/encoder vocab 16 2 32 1 (random/rng 2)) (random/rng 4)) sents (random/rng 5)) 1500))
          rng (random/rng 9)
          tags-of (fn [w] (set (map first (get-in parsing/ambiguous-grammar [:lexicon w]))))
          results (for [s test
                        :let [words (mapv first s) [masked pos] (tf/mask-words words 0.2 rng) got (tf/fill-in m masked)]
                        [i g] (map vector pos got)]
                    (contains? (tags-of g) (second (s i))))]
      (is (> (/ (double (count (filter true? results))) (count results)) 0.7)))))

(deftest learned-positions-test
  (testing "with a learned embedding for each position, as the chapter has it, in place of the sinusoids: every word tagged right"
    (let [{:keys [train test vocab]} @data
          m (:model (nth (tf/train-tagger-seq (tf/tagger (tf/encoder vocab 16 2 32 1 (random/rng 2) :max-len 8) tags (random/rng 3)) train) 800))]
      (is (= 1.0 (tf/tag-accuracy m test))))))

(deftest causal-mask-test
  (testing "each position may look at itself and those before it"
    (is (= [[0.0 -1e9 -1e9] [0.0 0.0 -1e9] [0.0 0.0 0.0]] (tf/causal-mask 3))))
  (testing "so the decoder's outputs at a position are untouched by the tokens after it"
    (let [m (tf/encoder-decoder (s2s/vocabulary (range 5)) 8 2 16 1 (random/rng 1))
          P (ad/consts (:params m))
          [H] (tf/encode m P [2 3 4])
          out #(ad/value (first (#'tf/decode m P H %)))
          a (out [0 2 3 4])
          b (out [0 2 5 6])]
      (is (= (take 2 a) (take 2 b)))
      (is (not= (drop 2 a) (drop 2 b))))))

(deftest encoder-decoder-test
  (testing "the whole transformer, encoder and decoder, reverses strings of six to eight digits it has not seen"
    (let [rng (random/rng 1)
          digits #(vec (repeatedly (+ 6 (random/below rng 3)) (fn [] (random/below rng 5))))
          pairs (vec (repeatedly 500 (fn [] (let [d (digits)] [d (vec (reverse d))]))))
          tests (vec (repeatedly 30 digits))
          m (:model (nth (tf/train-translator-seq (tf/encoder-decoder (s2s/vocabulary (range 5)) 16 2 32 1 (random/rng 2)) pairs) 2000))]
      (is (>= (count (filter #(= (vec (reverse %)) (:target (tf/translate m %))) tests)) 29))
      (is (= 10 (count (filter #(= (vec (reverse %)) (:target (tf/translate m % :beam 3))) (take 10 tests))))))))
