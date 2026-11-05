(ns allgo.language.window
  "A feed-forward part-of-speech tagger over a window of words, learning
  its word embeddings as it learns to tag: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, section
  24.1, the procedure the section gives in steps.

    1. A window of w words, w odd, centered on the word to tag; each
       sentence padded so every word has one.
    2. A vocabulary of the words seen more than `min-count` times, the
       rest read as `<unk>`.
    3. Sorted.
    4-5. An embedding matrix E, a row of d for each word, started at
       random.
    6. z1 = sigmoid(W1 x), z2 = sigmoid(W2 z1), y = softmax(Wout z2).
    7. x the window's embeddings concatenated, w d long: a word's
       embedding is the same wherever it is in the window, but each
       position meets its own part of W1, so the network knows where.
    8. E and the weights trained together by gradient descent -- here
       Adam, `allgo.learning.autodiff` differentiating, a sentence a step.

  Unlike a one-way recurrent tagger it sees the words after the one it
  tags, as far as the window reaches; unlike a bidirectional one, no
  further."
  (:require [allgo.language.recurrent :as rec]
            [allgo.learning.autodiff :as ad]))

(def pad "<pad>")
(def unknown "<unk>")

(defn vocabulary
  "The words of `sentences` seen more than `min-count` times, sorted,
  after `<pad>` and `<unk>`: `{:words :index}`."
  [sentences min-count]
  (let [counts (frequencies (apply concat sentences))
        ws (vec (concat [pad unknown] (sort (keep (fn [[w n]] (when (> n min-count) w)) counts))))]
    {:words ws :index (zipmap ws (range))}))

(defn windows
  "Each word's window of `w` in `words`, centered on it and padded past
  either end."
  [words w]
  (let [half (quot w 2) padded (vec (concat (repeat half pad) words (repeat half pad)))]
    (mapv #(subvec padded % (+ % w)) (range (count words)))))

(defn tagger
  "The tagger: a window of `w`, embeddings of `d`, hidden layers of `h1`
  and `h2`, and an output for each of `tags`."
  [vocab tags w d h1 h2 rng]
  (let [tags (vec tags) zeros #(vector (vec (repeat % 0.0)))]
    {:vocab vocab :tags tags :tag-index (zipmap tags (range)) :w w
     :params {:E (ad/init (count (:words vocab)) d rng)
              :W1 (ad/init (* w d) h1 rng) :b1 (zeros h1)
              :W2 (ad/init h1 h2 rng) :b2 (zeros h2)
              :Wout (ad/init h2 (count tags) rng) :bout (zeros (count tags))}}))

(defn- logits
  "A row of tag scores for each word of `words`."
  [{:keys [vocab w]} P words]
  (let [{:keys [index]} vocab
        x (ad/stack-rows (for [win (windows words w)]
                           (ad/concat-cols (map #(ad/row (:E P) (get index % (index unknown))) win))))
        z1 (ad/sigmoid (ad/add-row (ad/matmul x (:W1 P)) (:b1 P)))
        z2 (ad/sigmoid (ad/add-row (ad/matmul z1 (:W2 P)) (:b2 P)))]
    (ad/add-row (ad/matmul z2 (:Wout P)) (:bout P))))

(defn loss
  "The mean cross-entropy of the tags of `tagged`, `[[word tag] ...]`."
  [m P tagged]
  (ad/cross-entropy (logits m P (mapv first tagged)) (mapv #((:tag-index m) (second %)) tagged)))

(defn train-seq
  "The tagger trained on `tagged` sentences by Adam, a sentence a step:
  each state `{:model :loss}`."
  [m tagged & {:keys [lr] :or {lr 0.01}}]
  (rec/adam-seq m tagged loss lr nil))

(defn tag
  "The tags `m` gives `words`."
  [m words]
  (let [P (into {} (map (fn [[k v]] [k (ad/const v)]) (:params m)))]
    (mapv (fn [row] ((:tags m) (apply max-key row (range (count row)))))
          (ad/value (logits m P words)))))

(defn accuracy
  "The fraction of the words of `tagged` sentences tagged right, or of
  just those `(only? word i)` picks out."
  ([m tagged] (accuracy m tagged (constantly true)))
  ([m tagged only?]
   (let [rs (for [s tagged
                  [i [w t] g] (map vector (range) s (tag m (mapv first s)))
                  :when (only? w i)]
              (= t g))]
     (/ (double (count (filter true? rs))) (max 1 (count rs))))))
