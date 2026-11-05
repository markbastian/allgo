(ns allgo.language.window-test
  (:require [allgo.language.parsing :as parsing]
            [allgo.language.window :as win]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The feed-forward window tagger of Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 24.1.

(deftest windows-test
  (testing "each word's window, centered and padded"
    (is (= [["<pad>" "<pad>" "they" "cut" "it"]
            ["<pad>" "they" "cut" "it" "<pad>"]
            ["they" "cut" "it" "<pad>" "<pad>"]]
           (win/windows ["they" "cut" "it"] 5))))
  (testing "the vocabulary keeps the words seen more than min-count times, sorted, after the markers"
    (is (= ["<pad>" "<unk>" "a" "the"]
           (:words (win/vocabulary [["the" "cat"] ["a" "dog" "the"] ["a"]] 1))))))

(def ^:private right-grammar
  "A sentence's first word is a noun when an intransitive verb follows
  it, \"duck swim today\", and a verb when an object does, \"duck stones
  today\": only the word after decides."
  {:rules [[:S :A :Adv 1.0] [:A :N :Vi 0.5] [:A :V :Obj 0.5]]
   :lexicon {"duck" [[:N 0.2] [:V 0.2]] "fish" [[:N 0.2] [:V 0.2]] "watch" [[:N 0.2] [:V 0.2]] "saw" [[:N 0.2] [:V 0.2]]
             "dogs" [[:N 0.2]] "feed" [[:V 0.2]]
             "swim" [[:Vi 0.34]] "sleep" [[:Vi 0.33]] "run" [[:Vi 0.33]]
             "stones" [[:Obj 0.34]] "boats" [[:Obj 0.33]] "bread" [[:Obj 0.33]]
             "today" [[:Adv 0.34]] "again" [[:Adv 0.33]] "here" [[:Adv 0.33]]}})

(deftest window-tagging-test
  (let [rng (random/rng 1)
        data (vec (repeatedly 300 #(parsing/generate right-grammar rng)))
        [train test] [(subvec data 0 200) (subvec data 200)]
        voc (win/vocabulary (map #(mapv first %) train) 0)
        trained (fn [w] (:model (nth (win/train-seq (win/tagger voc [:N :V :Vi :Obj :Adv] w 8 16 16 (random/rng 2)) train) 600)))
        first-ambiguous (fn [w i] (and (zero? i) (#{"duck" "fish" "watch" "saw"} w)))
        one (trained 1)
        three (trained 3)]
    (testing "a window of one word is the word alone: the first word's tag a guess"
      (is (< 0.4 (win/accuracy one test first-ambiguous) 0.6)))
    (testing "a window of three sees the next word, which decides it"
      (is (= 1.0 (win/accuracy three test first-ambiguous)))
      (is (= 1.0 (win/accuracy three test))))))

(def ^:private big-grammar
  "The ambiguous grammar's rules over 61 words."
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

(deftest unseen-words-test
  (let [rng (random/rng 1)
        data (vec (repeatedly 330 #(parsing/generate big-grammar rng)))
        [train test] [(subvec data 0 30) (subvec data 30)]
        seen (set (mapcat #(map first %) train))
        unseen? (fn [w _] (not (seen w)))
        trained (fn [min-count]
                  (:model (nth (win/train-seq (win/tagger (win/vocabulary (map #(mapv first %) train) min-count)
                                                          [:Det :N :V :Adj :Modal :Adv] 5 8 16 16 (random/rng 2))
                                              train)
                               600)))
        all-words (trained 0)]
    (testing "thirty labeled sentences leave many test words never labeled, and the window tags them from their company"
      (is (< 150 (count (for [s test [w] s :when (not (seen w))] w))))
      (is (> (win/accuracy all-words test unseen?) 0.9)))
    (testing "reading the rare training words as <unk> too does not help here"
      (is (< (win/accuracy (trained 1) test unseen?) (win/accuracy all-words test unseen?))))))
