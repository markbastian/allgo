(ns allgo.language.tagging
  "Part-of-speech tagging with a hidden Markov model: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 23.

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
            [allgo.probability.temporal :as temporal]))

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

(defn accuracy
  "The fraction of the words of `tagged` sentences tagged right, or of just
  those `(only? word position)` picks out."
  ([m tagged] (accuracy m tagged (constantly true)))
  ([m tagged only?]
   (let [rs (for [s tagged
                  [i [w t] g] (map vector (range) s (tag m (mapv first s)))
                  :when (only? w i)]
              (= t g))]
     (/ (double (count (filter true? rs))) (max 1 (count rs))))))
