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

;; Stacked and bidirectional networks, contextual embeddings and
;; pretraining, after the same chapter.

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

(deftest bidirectional-test
  (let [rng (random/rng 1)
        data (vec (repeatedly 300 #(parsing/generate right-grammar rng)))
        [train test] [(subvec data 0 200) (subvec data 200)]
        voc (rec/vocabulary (map #(mapv first %) data))
        trained (fn [layers bi] (:model (nth (rec/train-tagger-seq (rec/tagger :lstm voc [:N :V :Vi :Obj :Adv] 8 12 layers bi (random/rng 2)) train) 600)))
        first-ambiguous (fn [w i] (and (zero? i) (#{"duck" "fish" "watch" "saw"} w)))]
    (testing "reading left to right, the first word's tag is a guess: the word that decides has not come yet"
      (is (< 0.4 (rec/tag-accuracy (trained 1 false) test first-ambiguous) 0.6)))
    (testing "reading both ways, every one right"
      (is (= 1.0 (rec/tag-accuracy (trained 1 true) test first-ambiguous))))
    (testing "stacking a second layer does not give a one-way network the words it has not read; stacked both ways, all right"
      (is (< 0.4 (rec/tag-accuracy (trained 2 false) test first-ambiguous) 0.6))
      (is (= 1.0 (rec/tag-accuracy (trained 2 true) test))))))

(deftest contextual-test
  (testing "a word's contextual embedding follows its use: \"duck\" the noun nearer \"duck\" the noun than \"duck\" the verb"
    (let [rng (random/rng 1)
          tagged (vec (repeatedly 600 #(parsing/generate parsing/ambiguous-grammar rng)))
          sents (mapv #(mapv first %) tagged)
          lm (:model (nth (rec/train-bilm-seq (rec/bilm :lstm (rec/vocabulary sents) 8 12 (random/rng 2)) sents) 1200))
          uses (for [s (take 300 tagged) :let [vs (rec/contextual lm (mapv first s))]
                     [[w t] v] (map vector s vs) :when (= w "duck")]
                 [t v])
          cosine (fn [a b] (/ (reduce + (map * a b)) (Math/sqrt (* (reduce + (map * a a)) (reduce + (map * b b))))))
          mean (fn [xs] (/ (reduce + xs) (count xs)))
          same (mean (for [[t1 a] uses [t2 b] uses :when (and (= t1 t2) (not (identical? a b)))] (cosine a b)))
          other (mean (for [[t1 a] uses [t2 b] uses :when (not= t1 t2)] (cosine a b)))]
      (is (= #{:N :V} (set (map first uses))))
      (is (> same (+ other 0.2))))))

(def ^:private big-grammar
  "The ambiguous grammar's rules over a vocabulary four times its size:
  twenty nouns, sixteen verbs, eight adjectives, six adverbs."
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

(deftest pretraining-test
  (let [rng (random/rng 1)
        data (vec (repeatedly 1710 #(parsing/generate big-grammar rng)))
        unlabeled (mapv #(mapv first %) (subvec data 10 1510))
        test (subvec data 1510)
        voc (rec/vocabulary (map #(mapv first %) data))
        tags [:Det :N :V :Adj :Modal :Adv]
        lm (:model (nth (rec/train-bilm-seq (rec/bilm :lstm voc 16 16 (random/rng 3)) unlabeled) 3000))
        run (fn [n]
              (let [labeled (subvec data 0 n)
                    seen (set (mapcat #(map first %) labeled))
                    unseen? (fn [w _] (not (seen w)))
                    features (fn [s] (rec/contextual lm (mapv first s)))
                    examples (mapv (fn [s] [(features s) (mapv second s)]) labeled)
                    probe (:model (nth (rec/train-probe-seq (rec/probe 48 tags (random/rng 4)) examples) 300))
                    probe-acc (fn [only?] (let [rs (for [s test [[w t] g] (map vector s (rec/probe-tags probe (features s))) :when (only? w nil)] (= t g))]
                                            (/ (double (count (filter true? rs))) (max 1 (count rs)))))
                    scratch (:model (nth (rec/train-tagger-seq (rec/tagger :lstm voc tags 16 16 1 true (random/rng 2)) labeled) 300))
                    tuned (:model (nth (rec/train-tagger-seq (rec/fine-tune lm tags (random/rng 5)) labeled) 300))]
                {:unseen (count (for [s test [w] s :when (unseen? w nil)] w))
                 :scratch [(rec/tag-accuracy scratch test) (rec/tag-accuracy scratch test unseen?)]
                 :probe [(probe-acc (constantly true)) (probe-acc unseen?)]
                 :tuned [(rec/tag-accuracy tuned test) (rec/tag-accuracy tuned test unseen?)]}))
        few (run 10)
        many (run 100)]
    (testing "ten labeled sentences leave a third of the test's words unseen: there, the language model's pretraining pays"
      (is (< 300 (:unseen few)))
      (is (< (second (:scratch few)) 0.7))
      (is (> (second (:probe few)) (+ 0.15 (second (:scratch few)))))
      (is (> (second (:tuned few)) (+ 0.1 (second (:scratch few)))))
      (is (> (first (:tuned few)) (+ 0.05 (first (:scratch few))))))
    (testing "a hundred cover every word, and the network trained from scratch needs no help"
      (is (zero? (:unseen many)))
      (is (= 1.0 (first (:scratch many))))
      (is (>= (first (:scratch many)) (first (:probe many)))))))

(defn- review
  "A review whose verdict is its fourth word, followed by `tail` to
  `tail` + 2 words that say nothing: `[words class]`."
  [rng tail]
  (let [pick #(nth % (random/below rng (count %)))
        good? (zero? (random/below rng 2))
        verdict (if good? (pick ["great" "moving" "superb" "fine" "clever"]) (pick ["awful" "dull" "poor" "weak" "clumsy"]))
        filler ["and" "the" "film" "was" "shot" "in" "a" "town" "by" "river" "with" "its" "cast" "of" "actors" "over" "summer"]]
    [(into ["this" "movie" "was" verdict] (repeatedly (+ tail (random/below rng 3)) #(pick filler)))
     (if good? :positive :negative)]))

(deftest sentence-classification-test
  (let [run (fn [tail pooling seed]
              (let [rng (random/rng 1)
                    data (vec (repeatedly 400 #(review rng tail)))
                    m (rec/classifier :rnn (rec/vocabulary (map first data)) [:positive :negative] 8 12 pooling (random/rng seed))]
                (rec/class-accuracy (:model (nth (rec/train-classifier-seq m (subvec data 0 300)) 300)) (subvec data 300))))]
    (testing "with the verdict near the end, the last state classifies every review"
      (is (= 1.0 (run 2 :last 3))))
    (testing "ten words further back, the last state has lost it -- the bias toward the end the chapter warns of"
      (is (< (run 10 :last 3) 0.6))
      (is (< (run 10 :last 4) 0.6)))
    (testing "average pooling weighs every word's state alike, and still classifies every one"
      (is (= 1.0 (run 10 :mean 3)))
      (is (= 1.0 (run 10 :mean 4))))))
