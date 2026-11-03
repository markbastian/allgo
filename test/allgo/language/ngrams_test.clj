(ns allgo.language.ngrams-test
  (:require [allgo.language.ngrams :as ng]
            [allgo.language.parsing :as parsing]
            [allgo.procedural.naming :as naming]
            [allgo.random :as random]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

;; N-gram models, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 23: worked by hand on tiny
;; corpora, then on sentences drawn from a grammar, and on the words of
;; languages drawn by `allgo.procedural.naming`.

(def ^:private grammar
  "The ambiguous grammar's rules over sixty words."
  (let [words (fn [tag ws] (into {} (map (fn [w] [w [[tag (/ 1.0 (count ws))]]]) ws)))]
    {:rules (:rules parsing/ambiguous-grammar)
     :lexicon (merge (words :N ["dog" "cat" "horse" "cow" "bird" "farmer" "child" "teacher" "boat" "house"
                                "tree" "apple" "river" "stone" "letter" "song" "garden" "window" "road" "city"])
                     (words :V ["sees" "likes" "finds" "keeps" "paints" "follows" "carries" "builds" "hears" "visits"
                                "chases" "loves" "brings" "moves" "washes" "sells"])
                     (words :Adj ["young" "old" "big" "small" "red" "green" "happy" "quiet"])
                     (words :Adv ["quickly" "well" "slowly" "often" "today" "again"])
                     {"the" [[:Det 0.6]] "a" [[:Det 0.4]] "will" [[:Modal 0.5]] "can" [[:Modal 0.3] [:N 0.02]] "must" [[:Modal 0.2]]
                      "fish" [[:N 0.03] [:V 0.03]] "duck" [[:N 0.03] [:V 0.03]] "watch" [[:N 0.03] [:V 0.03]] "saw" [[:N 0.03] [:V 0.03]]})}))

(deftest counting
  (let [m (ng/model 2 [["a" "b"] ["a" "a"]])]
    (testing "bigram counts, the sequences padded with a start before and an end after"
      (is (= {"a" 2} (get-in m [:counts 1 [ng/start]])))
      (is (= {"b" 1 "a" 1 ng/end 1} (get-in m [:counts 1 ["a"]])) "a is followed by b, by a, and by the end")
      (is (= #{"a" "b" ng/end} (:vocab m))))
    (testing "the maximum-likelihood estimate, and Laplace's add-one"
      (is (= (/ 1.0 3.0) (ng/probability m :none ["a"] "b")))
      (is (= 0.0 (ng/probability m :none ["b"] "a")) "never seen after b")
      ;; after "a": counts b 1, a 1, end 1, over 3, with 3 tokens: (1 + 1) / (3 + 3)
      (is (= (/ 2.0 6.0) (ng/probability m :laplace ["a"] "b")))
      (is (= (/ 1.0 4.0) (ng/probability m :laplace ["b"] "a"))))))

(deftest perplexity
  (testing "a model of equally likely tokens has the perplexity of their number: here a, b and the end"
    (let [m (ng/model 1 (repeat 10 ["a" "b"]))]
      (is (< (abs (- 3.0 (ng/perplexity m :none [["a" "b"] ["b" "a"]]))) 1e-9))))
  (let [rng (random/rng 7)
        sents (vec (repeatedly 3000 #(mapv first (parsing/generate grammar rng))))
        [train held test] [(subvec sents 0 2000) (subvec sents 2000 2500) (subvec sents 2500)]
        fitted (fn [n] (let [m (ng/model n train)
                             em (vec (take 30 (ng/interpolation-seq m held)))]
                         {:m m :em em :lambdas (:lambdas (peek em))}))
        [uni bi tri] (map fitted [1 2 3])
        ppl (fn [{:keys [m lambdas]}] (ng/perplexity m [:interpolate lambdas] test))]
    (testing "EM for the interpolation weights: the held-out likelihood never falls, and the weights sum to one"
      (doseq [{:keys [em lambdas]} [uni bi tri]]
        (is (apply <= (map :log-likelihood em)))
        (is (< (abs (- 1.0 (reduce + lambdas))) 1e-9))))
    (testing "on unseen sentences, the bigram model far less perplexed than the unigram; the trigram no worse, the grammar's words depending mostly on the one before"
      (is (< (* 2 (ppl bi)) (ppl uni)))
      (is (<= (ppl tri) (+ (ppl bi) 0.05))))
    (testing "Laplace's add-one spreads too much over the trigrams never seen: worse than interpolation, worse than its own bigram"
      (is (> (ng/perplexity (:m tri) :laplace test) (* 1.5 (ppl tri))))
      (is (> (ng/perplexity (:m tri) :laplace test) (ng/perplexity (:m bi) :laplace test))))
    (testing "samples are made of the words the model knows"
      (let [s (ng/sample (:m tri) [:interpolate (:lambdas tri)] rng)]
        (is (seq s))
        (is (every? (:vocab (:m tri)) s))))))

(deftest unknown-words
  (testing "each word's first appearance becomes <UNK>"
    (is (= [["<UNK>" "<UNK>" "a"] ["b" "<UNK>"]] (ng/with-unknown [["a" "b" "a"] ["b" "c"]]))))
  (testing "a word never seen takes <UNK>'s probability, so unseen text has a finite perplexity"
    (let [rng (random/rng 3)
          sents (vec (repeatedly 400 #(mapv first (parsing/generate grammar rng))))
          m (ng/model 2 (ng/with-unknown (subvec sents 0 300)))
          lambdas (:lambdas (nth (ng/interpolation-seq m (map #(ng/known m %) (subvec sents 300 350))) 20))
          test (map #(ng/known m %) (conj (subvec sents 350) ["the" "zebra" "sees" "a" "dog"]))]
      (is (some #(= ng/unknown %) (last test)))
      (is (< (ng/perplexity m [:interpolate lambdas] test) 30.0)))))

(deftest identification
  (testing "the language of a single word, by character trigrams, among four languages drawn by the naming code"
    (let [langs (mapv naming/language [3 11 29 47])
          words (fn [lang from to] (mapv #(mapv str (str/lower-case (naming/word lang %))) (range from to)))
          models (into {} (map-indexed (fn [i l] [i [(ng/model 3 (words l 0 300)) [:interpolate [0.01 0.09 0.3 0.6]]]]) langs))
          results (for [[i l] (map-indexed vector langs) w (words l 1000 1100)] (= i (ng/identify models w)))]
      (is (> (/ (count (filter true? results)) (double (count results))) 0.95)))))
