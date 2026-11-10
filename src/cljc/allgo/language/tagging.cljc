(ns allgo.language.tagging
  "Part-of-speech tagging, section 23.1.6 of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition: by a hidden
  Markov model, generative, and by logistic regression, discriminative.

  The logistic tagger has a model for each tag over binary features of a
  word in its context -- the word, its spelling, its neighbors, the tags
  already given the words before -- and tags a sentence left to right:
  greedily, each word's likeliest tag taken and never revisited, or by a
  beam search keeping the b likeliest tag sequences, greedy being b = 1.
  It can say what an HMM cannot, that a word ending in \"ous\" is likely an
  adjective; the HMM, below, weighs the whole sentence.

  The tags are the hidden states and the words the evidence: P(tag_1),
  P(tag | the tag before), and P(word | tag), each counted from tagged
  sentences with Laplace's add-one, and the unknown word handled as the
  n-gram models handle it -- each word's first appearance counted as
  `<UNK>`, with its tag. A sentence's tags are then the most likely
  sequence of states given its words: chapter 14's Viterbi algorithm,
  `allgo.probability.temporal`'s, which weighs the whole sentence -- the
  words after a tag as well as before. Chapter 14's HMM starts from
  P(X0) and moves a step before the first evidence, so the sentence's
  start is a state of its own, as the n-gram models' start marker is: all
  of the prior, never a word's, and its transitions are P(first tag)."
  (:require [allgo.language.ngrams :as ngrams]
            [allgo.math :as am]
            [allgo.probability.temporal :as temporal]
            [clojure.math :as math]
            [clojure.string :as str]))

(defn- normalized [counts keys]
  (let [total (reduce + (map #(+ 1.0 (get counts % 0)) keys))]
    (mapv #(/ (+ 1.0 (get counts % 0)) total) keys)))

(defn train
  "An HMM tagger counted from `tagged` sentences, each `[[word tag] ...]`:
  `{:states :words :hmm}`, the HMM chapter 14's form, `:states` the start
  then the tags."
  [tagged]
  (let [words (ngrams/with-unknown (map #(mapv first %) tagged))
        tagged (map (fn [ws s] (mapv vector ws (map second s))) words tagged)
        tags (vec (sort (distinct (map second (apply concat tagged)))))
        vocab (vec (sort (distinct (map first (apply concat tagged)))))
        firsts (frequencies (map (comp second first) tagged))
        follows (reduce (fn [m [[_ a] [_ b]]] (update-in m [a b] (fnil inc 0))) {} (mapcat #(partition 2 1 %) tagged))
        emits (reduce (fn [m [w t]] (update-in m [t w] (fnil inc 0))) {} (apply concat tagged))
        emission (into {} (map (fn [t] [t (zipmap vocab (normalized (emits t) vocab))]) tags))]
    {:states (into [:start] tags)
     :words (set vocab)
     :hmm {:prior (into [1.0] (repeat (count tags) 0.0))
           :transition (into [(into [0.0] (normalized firsts tags))]
                             (map #(into [0.0] (normalized (follows %) tags)) tags))
           :sensor (fn [w] (into [0.0] (map #(get-in emission [% w]) tags)))}}))

(defn tag
  "The tags of `words`, the likeliest sequence by Viterbi; words never
  seen in training read as `<UNK>`."
  [{:keys [states words hmm]} sentence]
  (mapv states (:path (temporal/viterbi hmm (mapv #(if (words %) % ngrams/unknown) sentence)))))

(defn accuracy-of
  "The fraction of the words of `tagged` sentences, each `[[word tag]
  ...]`, that `tag-fn` (fn [words]) -> tags gets right, or of just those
  `(only? word position)` picks out: every tagger's accuracy."
  ([tag-fn tagged] (accuracy-of tag-fn tagged (constantly true)))
  ([tag-fn tagged only?]
   (let [rs (for [s tagged
                  [i [w t] g] (map vector (range) s (tag-fn (mapv first s)))
                  :when (only? w i)]
              (= t g))]
     (/ (double (count (filter true? rs))) (max 1 (count rs))))))

(defn accuracy
  "The fraction of the words of `tagged` sentences the HMM tags right, or
  of just those `(only? word position)` picks out."
  ([m tagged] (accuracy-of #(tag m %) tagged))
  ([m tagged only?] (accuracy-of #(tag m %) tagged only?)))

;; ---------------------------------------------------------------------------
;; Logistic regression

(defn book-features
  "Binary features of the word at `i` of `words`, given the tags assigned
  so far, `tags` -- after the chapter's table of POS-tagging features:
  the word itself, its neighbors one and two away, the tag before, the
  word two back with the tag before, the words either side together, and
  the word's spelling -- its last three and two letters, its first two,
  a hyphen, a digit, all capitals. Each feature a vector, its value 1."
  [words tags i]
  (let [w (get words i)
        at (fn [j] (get words j "<s>"))
        tag-at (fn [j] (get tags j :start))
        lower (str/lower-case w)
        n (count lower)]
    (cond-> [[:bias]
             [:w lower]
             [:w-1 (str/lower-case (at (- i 1)))]
             [:w+1 (str/lower-case (at (+ i 1)))]
             [:w-2 (str/lower-case (at (- i 2)))]
             [:c-1 (tag-at (- i 1))]
             [:w-2-c-1 (str/lower-case (at (- i 2))) (tag-at (- i 1))]
             [:w-1-w+1 (str/lower-case (at (- i 1))) (str/lower-case (at (+ i 1)))]
             [:suffix3 (subs lower (max 0 (- n 3)))]
             [:suffix2 (subs lower (max 0 (- n 2)))]
             [:prefix2 (subs lower 0 (min n 2))]]
      (str/includes? w "-") (conj [:hyphen])
      (re-find #"[0-9]" w) (conj [:digit])
      (and (re-matches #"[A-Z]+" w) (> n 1)) (conj [:upper]))))

(defn- score [weights features] (reduce + (map #(get weights % 0.0) features)))

(defn tag-distribution
  "P(c | words, the tags before) for the word at `i`: each tag's logistic
  model's probability that the word is of that tag, normalized over the
  tags."
  [{:keys [tags weights features]} words prior-tags i]
  (let [fs (features words prior-tags i)
        ps (map #(am/sigmoid (score (weights %) fs)) tags)
        z (reduce + ps)]
    (zipmap tags (map #(/ % z) ps))))

(defn logistic-step
  "One epoch of training the one-vs-rest models on `tagged` sentences:
  for each word, the features its correct context gives -- the tags
  before it the true ones -- and for each tag, its weights moved by
  `alpha` up the gradient of the log likelihood, (y - P) times each
  feature, y 1 for the word's tag and 0 for every other."
  [tagged alpha {:keys [tags features] :as model}]
  (reduce (fn [model sentence]
            (let [words (mapv first sentence) gold (mapv second sentence)]
              (reduce (fn [model i]
                        (let [fs (features words gold i)]
                          (reduce (fn [model t]
                                    (let [w (get-in model [:weights t])
                                          err (- (if (= t (gold i)) 1.0 0.0) (am/sigmoid (score w fs)))]
                                      (assoc-in model [:weights t] (reduce #(update %1 %2 (fnil + 0.0) (* alpha err)) w fs))))
                                  model tags)))
                      model (range (count words)))))
          model tagged))

(defn logistic-tagger-seq
  "A POS tagger by logistic regression, section 23.1.6: one model for
  each tag, a weight for each binary feature, trained by stochastic
  gradient descent on `tagged` sentences, an epoch a step. `features` is
  `(fn [words tags i])` -> the features of word `i`, `book-features` by
  default."
  ([tagged] (logistic-tagger-seq tagged book-features 0.1))
  ([tagged features alpha]
   (let [tags (vec (sort (distinct (map second (apply concat tagged)))))]
     (iterate (partial logistic-step tagged alpha) {:tags tags :features features :weights (zipmap tags (repeat {}))}))))

(defn beam-tag-seq
  "Tagging `words` by beam search, a word a step: each state the `b`
  likeliest tag sequences so far, `[[tags log-p] ...]`, each extended by
  every tag and the `b` best kept. With b = 1 it is the greedy search,
  each word's tag argmax P(c | words, the tags before) and never
  revisited."
  [model words b]
  (reductions (fn [beam i]
                (->> (for [[tags lp] beam [c p] (tag-distribution model words tags i)]
                       [(conj tags c) (+ lp (math/log p))])
                     (sort-by (comp - second))
                     (take b)
                     vec))
              [[[] 0.0]]
              (range (count words))))

(defn beam-tag
  "The tags of `words`, the best of a beam of `b`."
  [model words b]
  (ffirst (last (beam-tag-seq model words b))))

(defn greedy-tag
  "The tags of `words` by greedy search: a beam of 1."
  [model words]
  (beam-tag model words 1))
