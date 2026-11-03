(ns allgo.language.tagging-test
  (:require [allgo.language.parsing :as parsing]
            [allgo.language.tagging :as tagging]
            [allgo.language.transformer :as tf]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Part-of-speech tagging with an HMM, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 23,
;; on sentences drawn from grammars whose words are often two parts of
;; speech.

(def ^:private ambiguous #{"fish" "duck" "watch" "saw" "can"})

(deftest hmm-tagging
  (let [rng (random/rng 1)
        tagged (vec (repeatedly 500 #(parsing/generate parsing/ambiguous-grammar rng)))
        [train test] [(subvec tagged 0 400) (subvec tagged 400)]
        m (tagging/train train)
        base (tf/most-frequent-tag train)
        base-acc (let [rs (for [s test [w t] s :when (ambiguous w)] (= t (base w)))]
                   (/ (double (count (filter true? rs))) (count rs)))]
    (testing "each word's most frequent tag gets the ambiguous words wrong four times in ten"
      (is (< base-acc 0.6)))
    (testing "the HMM, weighing the tags around each, gets every word"
      (is (= 1.0 (tagging/accuracy m test)))
      (is (= 1.0 (tagging/accuracy m test (fn [w _] (ambiguous w))))))
    (testing "an unknown word is tagged as <UNK> is, by its neighbors"
      (is (= [:Det :N :V :Det :N] (tagging/tag m ["the" "zebra" "sees" "a" "dog"]))))))

(deftest viterbi-looks-ahead
  (testing "where only the next word decides the tag, Viterbi -- weighing the whole sentence -- still gets it"
    (let [g {:rules [[:S :A :Adv 1.0] [:A :N :Vi 0.5] [:A :V :Obj 0.5]]
             :lexicon {"duck" [[:N 0.25] [:V 0.25]] "fish" [[:N 0.25] [:V 0.25]] "dogs" [[:N 0.5]] "feed" [[:V 0.5]]
                       "swim" [[:Vi 0.5]] "sleep" [[:Vi 0.5]] "stones" [[:Obj 0.5]] "bread" [[:Obj 0.5]]
                       "today" [[:Adv 0.5]] "again" [[:Adv 0.5]]}}
          rng (random/rng 1)
          data (vec (repeatedly 300 #(parsing/generate g rng)))
          m (tagging/train (subvec data 0 200))]
      (is (= 1.0 (tagging/accuracy m (subvec data 200) (fn [w i] (and (zero? i) (#{"duck" "fish"} w))))))
      (is (= [:N :Vi :Adv] (tagging/tag m ["duck" "swim" "today"])))
      (is (= [:V :Obj :Adv] (tagging/tag m ["duck" "stones" "today"]))))))
