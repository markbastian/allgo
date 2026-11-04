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

;; Tagging by logistic regression, section 23.1.6.

(def ^:private look-ahead-grammar
  "Sentences whose first word, duck or fish, is a noun or a verb as the
  next word is an intransitive verb or an object."
  {:rules [[:S :A :Adv 1.0] [:A :N :Vi 0.5] [:A :V :Obj 0.5]]
   :lexicon {"duck" [[:N 0.25] [:V 0.25]] "fish" [[:N 0.25] [:V 0.25]] "dogs" [[:N 0.5]] "feed" [[:V 0.5]]
             "swim" [[:Vi 0.5]] "sleep" [[:Vi 0.5]] "stones" [[:Obj 0.5]] "bread" [[:Obj 0.5]]
             "today" [[:Adv 0.5]] "again" [[:Adv 0.5]]}})

(defn- logistic-accuracy [tag tagged]
  (let [rs (for [s tagged [[_ t] g] (map vector s (tag (mapv first s)))] (= t g))]
    (/ (double (count (filter true? rs))) (count rs))))

(deftest logistic-tagging
  (testing "the chapter's features: the word, its spelling, its neighbors, the tags before -- the ambiguous grammar's words all right"
    (let [rng (random/rng 1)
          tagged (vec (repeatedly 500 #(parsing/generate parsing/ambiguous-grammar rng)))
          m (nth (tagging/logistic-tagger-seq (subvec tagged 0 400)) 3)]
      (is (= 1.0 (logistic-accuracy #(tagging/greedy-tag m %) (subvec tagged 400))))))
  (let [rng (random/rng 1)
        data (vec (repeatedly 300 #(parsing/generate look-ahead-grammar rng)))
        [train test] [(subvec data 0 200) (subvec data 200)]
        ;; only the word and the tag before: no look at the next word
        local (fn [words tags i] [[:bias] [:w (words i)] [:c-1 (get tags (dec i) :start)]])
        m (nth (tagging/logistic-tagger-seq train local 0.1) 10)]
    (testing "where only what comes after decides a tag, greedy search commits too soon; a beam keeps the alternative until the next word settles it"
      (is (< (logistic-accuracy #(tagging/greedy-tag m %) test) 0.9))
      (is (> (logistic-accuracy #(tagging/beam-tag m % 2) test) 0.92))
      (is (= [:V :Obj :Adv] (tagging/greedy-tag m ["duck" "swim" "today"])))
      (is (= [:N :Vi :Adv] (tagging/beam-tag m ["duck" "swim" "today"] 2))))
    (testing "the classifier may look at any word of the sentence, and with the next word among its features the greedy search needs no beam"
      (is (= 1.0 (logistic-accuracy #(tagging/greedy-tag (nth (tagging/logistic-tagger-seq train) 10) %) test))))
    (testing "a beam of b keeps b sequences, a step a word"
      (is (= [1 2 2 2] (map count (tagging/beam-tag-seq m ["duck" "swim" "today"] 2))))))
  (testing "what an HMM cannot be told: any word ending in -ous is likely an adjective -- words never seen get :Adj by their spelling"
    (let [g {:rules [[:S :NP :VP 1.0] [:NP :Det :N 0.5] [:NP :Adj :N 0.5] [:VP :V :NP 1.0]]
             :lexicon {"the" [[:Det 1.0]] "famous" [[:Adj 0.25]] "nervous" [[:Adj 0.25]] "dangerous" [[:Adj 0.25]] "curious" [[:Adj 0.25]]
                       "dog" [[:N 0.25]] "cat" [[:N 0.25]] "bus" [[:N 0.25]] "house" [[:N 0.25]]
                       "sees" [[:V 0.5]] "likes" [[:V 0.5]]}}
          rng (random/rng 3)
          m (nth (tagging/logistic-tagger-seq (vec (repeatedly 200 #(parsing/generate g rng)))) 5)]
      (is (= [:Adj :N :V :Adj :N] (tagging/greedy-tag m ["jealous" "dog" "likes" "glorious" "house"]))))))
