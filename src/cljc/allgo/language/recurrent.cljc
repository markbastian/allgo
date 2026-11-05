(ns allgo.language.recurrent
  "Recurrent neural networks for language: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 24.

  A language model reads a sentence a word at a time and, at each, gives
  a distribution over the next. A recurrent network carries what it has
  read in a hidden vector: each word's embedding and the previous hidden
  vector make the next, and the hidden vector makes the prediction. Two
  cells:

    :rnn    h = tanh(x Wx + h Wh + b), Elman's -- whose gradients over a
            long sentence shrink or blow up as they pass back through
            each step
    :lstm   the long short-term memory: a cell vector c beside h, kept or
            forgotten by gates, h read out of it -- gradients flow back
            along c without being squeezed at every step

  Trained by back-propagation through time -- `allgo.learning.autodiff`
  differentiating the unrolled network -- with Adam and gradient clipping,
  a sentence a step. Then `sample` writes sentences, a word drawn at a
  time from the model's distribution (sharpened or flattened by a
  temperature), and `perplexity` measures how surprised the model is by
  sentences it has not seen: e to the mean negative log probability of
  each word.

  Sentences are vectors of words; the model adds `<s>` before and `</s>`
  after each.

  Below the language models: taggers of stacked layers, each layer
  bidirectional or not; a bidirectional language model whose states give
  each word an embedding of its context, as ELMo's do; and a probe, a
  softmax layer trained on frozen features, to see what those hold."
  (:require [allgo.learning.autodiff :as ad]
            [clojure.math :as math]))

(def start "<s>")
(def end "</s>")

(defn vocabulary
  "The words of `sentences` and the two markers: `{:words :index}`."
  [sentences]
  (let [ws (vec (concat [start end] (sort (distinct (apply concat sentences)))))]
    {:words ws :index (zipmap ws (range))}))

(defn model
  "A language model of `cell` (`:rnn` or `:lstm`) over `vocab`: embeddings
  of `dim`, a hidden vector of `hidden`. `{:cell :vocab :params}`."
  [cell vocab dim hidden rng]
  (let [v (count (:words vocab))]
    {:cell cell :vocab vocab :hidden hidden
     :params (merge {:E (ad/init v dim rng) :Wy (ad/init hidden v rng) :by [(vec (repeat v 0.0))]}
                    (case cell
                      :rnn {:Wx (ad/init dim hidden rng) :Wh (ad/init hidden hidden rng) :b [(vec (repeat hidden 0.0))]}
                      ;; the forget gate's bias starts at 1: remember, until
                      ;; taught otherwise
                      :lstm {:W (ad/init (+ dim hidden) (* 4 hidden) rng)
                             :b [(vec (concat (repeat hidden 0.0) (repeat hidden 1.0) (repeat (* 2 hidden) 0.0)))]}))}))

(defn step
  "One step of the cell from state `[h c]` on input row `x`: `[h c]`, c nil
  for the plain RNN."
  [cell P x [h c]]
  (case cell
    :rnn [(ad/tanh (ad/add-row (ad/add (ad/matmul x (:Wx P)) (ad/matmul h (:Wh P))) (:b P))) nil]
    :lstm (let [n (count (first (ad/value h)))
                z (ad/add-row (ad/matmul (ad/concat-cols [x h]) (:W P)) (:b P))
                gate (fn [k f] (f (ad/slice-cols z (* k n) (* (inc k) n))))
                i (gate 0 ad/sigmoid) f (gate 1 ad/sigmoid) o (gate 2 ad/sigmoid) g (gate 3 ad/tanh)
                c (ad/add (ad/mul f c) (ad/mul i g))]
            [(ad/mul o (ad/tanh c)) c])))

(defn initial-state [hidden]
  (let [z (ad/const [(vec (repeat hidden 0.0))])] [z z]))

(defn- ids [{:keys [index]} sentence] (mapv index (concat [start] sentence [end])))

(defn run
  "The logits after each input of `tokens` (indices), and the final
  state: `[logits-node state]`."
  [{:keys [cell hidden]} P tokens]
  (let [[outs st] (reduce (fn [[outs st] t]
                            (let [[h c] (step cell P (ad/row (:E P) t) st)]
                              [(conj outs (ad/add-row (ad/matmul h (:Wy P)) (:by P))) [h c]]))
                          [[] (initial-state hidden)] tokens)]
    [(ad/stack-rows outs) st]))

(defn loss
  "The mean cross-entropy of predicting each word of `sentence`, and the
  end, from what came before."
  [m P sentence]
  (let [t (ids (:vocab m) sentence)
        [logits] (run m P (pop t))]
    (ad/cross-entropy logits (subvec t 1))))

(defn train-seq
  "The model trained on `sentences` by Adam, a sentence a step, the
  gradient clipped at `:clip` (default 5): each state `{:model :loss}`,
  the loss that sentence's before its step. The sentences are taken in
  turn, again and again; shuffle them first."
  [m sentences & {:keys [lr clip] :or {lr 0.01 clip 5.0}}]
  (let [sentences (vec sentences)]
    (map first
         (iterate (fn [[_ m opt k]]
                    (let [s (sentences (mod k (count sentences)))
                          {:keys [value grads]} (ad/value-and-grad #(loss m % s) (:params m))
                          [params opt] (ad/adam-step (:params m) (ad/clip grads clip) opt :lr lr)
                          m (assoc m :params params)]
                      [{:model m :loss value} m opt (inc k)]))
                  [{:model m :loss nil} m nil 0]))))

(defn- distribution [m P st t]
  (let [[h c] (step (:cell m) P (ad/row (:E P) t) st)]
    [(ad/softmax (first (ad/value (ad/add-row (ad/matmul h (:Wy P)) (:by P))))) [h c]]))

(defn- consts [params] (into {} (map (fn [[k v]] [k (ad/const v)]) params)))

(defn next-distribution
  "The model's distribution over the word after `prefix`, `{word p}`."
  [{:keys [vocab hidden params] :as m} prefix]
  (let [P (consts params)
        t (pop (ids vocab prefix))
        [p] (reduce (fn [[_ st] tok] (distribution m P st tok)) [nil (initial-state hidden)] t)]
    (zipmap (:words vocab) p)))

(defn sample
  "A sentence drawn from the model, a word at a time, until `</s>` or
  `max-len` words; `temperature` (default 1) divides the logs of the
  probabilities, below 1 sharpening toward the likeliest word."
  [{:keys [vocab hidden params] :as m} rng & {:keys [temperature max-len] :or {temperature 1.0 max-len 20}}]
  (let [P (consts params) words (:words vocab)]
    (loop [st (initial-state hidden) t ((:index vocab) start) out []]
      (let [[p st] (distribution m P st t)
            p (ad/softmax (mapv #(/ (math/log (max 1e-300 %)) temperature) p))
            r (double (rng))
            nxt (loop [i 0 acc 0.0] (if (or (= i (dec (count p))) (< r (+ acc (p i)))) i (recur (inc i) (+ acc (p i)))))]
        (if (or (= (words nxt) end) (>= (count out) max-len))
          out
          (recur st nxt (conj out (words nxt))))))))

(defn log-probability
  "The natural log of the probability the model gives `sentence`, the
  end included."
  [{:keys [vocab hidden params] :as m} sentence]
  (let [P (consts params) t (ids vocab sentence)]
    (first (reduce (fn [[lp st] [a b]] (let [[p st] (distribution m P st a)] [(+ lp (math/log (max 1e-300 (p b)))) st]))
                   [0.0 (initial-state hidden)] (partition 2 1 t)))))

(defn perplexity
  "e to the mean negative log probability per predicted word -- the end
  of each sentence counted as a word -- over `sentences`."
  [m sentences]
  (math/exp (/ (- (reduce + (map #(log-probability m %) sentences)))
               (reduce + (map #(inc (count %)) sentences)))))

(defn unigram-perplexity
  "The perplexity of the model that knows only how often each word is
  used -- the baseline a language model must beat."
  [train test]
  (let [counts (frequencies (mapcat #(conj (vec %) end) train))
        total (reduce + (vals counts))
        p #(/ (get counts % 0.5) total)]
    (math/exp (/ (- (reduce + (map #(math/log (p %)) (mapcat #(conj (vec %) end) test))))
                 (reduce + (map #(inc (count %)) test))))))

;; ---------------------------------------------------------------------------
;; Stacked and bidirectional networks: tagging

(defn- cell-params
  "The parameters of one cell reading `in`-wide inputs into `hidden`,
  named with `prefix`."
  [cell in hidden prefix rng]
  (let [k #(keyword (str prefix (name %)))]
    (case cell
      :rnn {(k :Wx) (ad/init in hidden rng) (k :Wh) (ad/init hidden hidden rng) (k :b) [(vec (repeat hidden 0.0))]}
      :lstm {(k :W) (ad/init (+ in hidden) (* 4 hidden) rng)
             (k :b) [(vec (concat (repeat hidden 0.0) (repeat hidden 1.0) (repeat (* 2 hidden) 0.0)))]})))

(defn- cell-of
  "The parameters named with `prefix`, under the names `step` knows."
  [P cell prefix]
  (into {} (for [n (case cell :rnn [:Wx :Wh :b] :lstm [:W :b])] [n (P (keyword (str prefix (name n))))])))

(defn- direction-prefix [layer dir] (str "l" layer "-" (name dir) "-"))

(defn- run-cells
  "The hidden rows of a cell run over the row nodes `xs` in order."
  [cell P hidden xs]
  (first (reduce (fn [[hs st] x] (let [[h c] (step cell P x st)] [(conj hs h) [h c]]))
                 [[] (initial-state hidden)] xs)))

(defn encode-sequence
  "The rows `xs` through `layers` stacked recurrent layers, each reading
  the outputs of the one below; with `bidirectional`, each layer two
  networks, one reading left to right and one right to left, their
  states side by side -- so each position's output has heard the whole
  sentence on both sides of it. The rows of the top layer."
  [{:keys [cell hidden layers bidirectional]} P xs]
  (reduce (fn [xs l]
            (let [fwd (run-cells cell (cell-of P cell (direction-prefix l :f)) hidden xs)]
              (if bidirectional
                (let [bwd (vec (rseq (run-cells cell (cell-of P cell (direction-prefix l :b)) hidden (vec (rseq xs)))))]
                  (mapv #(ad/concat-cols [%1 %2]) fwd bwd))
                fwd)))
          xs (range layers)))

(defn tagger
  "A sequence tagger: embeddings of `dim`, `layers` recurrent layers of
  `cell` with `hidden` units -- each `bidirectional` or not -- and a layer
  from the top one's output to `tags`. The sentence is read between `<s>`
  and `</s>`, as the language models read it."
  [cell vocab tags dim hidden layers bidirectional rng]
  (let [dirs (if bidirectional [:f :b] [:f])
        width (* hidden (count dirs))
        tags (vec tags)]
    {:cell cell :vocab vocab :hidden hidden :layers layers :bidirectional bidirectional
     :tags tags :tag-index (zipmap tags (range))
     :params (merge {:E (ad/init (count (:words vocab)) dim rng)
                     :Wt (ad/init width (count tags) rng) :bt [(vec (repeat (count tags) 0.0))]}
                    (apply merge (for [l (range layers) d dirs]
                                   (cell-params cell (if (zero? l) dim width) hidden (direction-prefix l d) rng))))}))

(defn- word-ids [{:keys [index]} words] (mapv #(get index % (index end)) words))

(defn- tag-logits [m P words]
  (let [t (word-ids (:vocab m) (concat [start] words [end]))
        hs (encode-sequence m P (mapv #(ad/row (:E P) %) t))]
    (ad/add-row (ad/matmul (ad/stack-rows (subvec hs 1 (inc (count words)))) (:Wt P)) (:bt P))))

(defn tag-loss
  "The mean cross-entropy of the tags of `tagged`, `[[word tag] ...]`."
  [m P tagged]
  (ad/cross-entropy (tag-logits m P (mapv first tagged)) (mapv #((:tag-index m) (second %)) tagged)))

(defn adam-seq
  "Adam on `(loss m P example)` over `examples` in turn, clipped at 5,
  `only` the named parameters moved if given: `{:model :loss}` each step,
  the loss that example's before its step."
  [m examples loss lr only]
  (let [examples (vec examples)]
    (map first
         (iterate (fn [[_ m opt k]]
                    (let [ex (examples (mod k (count examples)))
                          {:keys [value grads]} (ad/value-and-grad #(loss m % ex) (:params m))
                          grads (if only (select-keys grads only) grads)
                          [params opt] (ad/adam-step (:params m) (ad/clip grads 5.0) opt :lr lr)
                          m (assoc m :params params)]
                      [{:model m :loss value} m opt (inc k)]))
                  [{:model m :loss nil} m nil 0]))))

(defn train-tagger-seq
  "The tagger trained on `tagged` sentences, a sentence a step: each state
  `{:model :loss}`. With `:train-only`, a set of parameter names, the rest
  stay as they are."
  [m tagged & {:keys [lr train-only] :or {lr 0.01}}]
  (adam-seq m tagged tag-loss lr train-only))

(defn tag
  "The tags `m` gives `words`."
  [m words]
  (mapv (fn [row] ((:tags m) (apply max-key row (range (count row)))))
        (ad/value (tag-logits m (consts (:params m)) words))))

(defn tag-accuracy
  "The fraction of words of `tagged` sentences tagged right, or of just
  those `only?` picks out -- by word and position, `(only? word i)`."
  ([m tagged] (tag-accuracy m tagged (constantly true)))
  ([m tagged only?]
   (let [rs (for [s tagged
                  [i [w t] g] (map vector (range) s (tag m (mapv first s)))
                  :when (only? w i)]
              (= t g))]
     (/ (double (count (filter true? rs))) (max 1 (count rs))))))

;; ---------------------------------------------------------------------------
;; Classifying whole sentences

(defn classifier
  "A sentence classifier, section 24.2.2: embeddings of `dim`, a
  recurrent layer of `cell` with `hidden` units, its hidden vectors
  pooled into one for the sentence -- `:last`, the state after the last
  word, which has read it all but most recently its end, or `:mean`,
  average pooling, every word's state weighed alike -- then a
  feed-forward layer of `hidden` and a softmax over `classes`."
  [cell vocab classes dim hidden pooling rng]
  (let [classes (vec classes) zeros #(vector (vec (repeat % 0.0)))]
    {:cell cell :vocab vocab :hidden hidden :layers 1 :bidirectional false :pooling pooling
     :classes classes :class-index (zipmap classes (range))
     :params (merge {:E (ad/init (count (:words vocab)) dim rng)
                     :Wf (ad/init hidden hidden rng) :bf (zeros hidden)
                     :Wc (ad/init hidden (count classes) rng) :bc (zeros (count classes))}
                    (cell-params cell dim hidden (direction-prefix 0 :f) rng))}))

(defn- class-logits [m P words]
  (let [hs (encode-sequence m P (mapv #(ad/row (:E P) %) (word-ids (:vocab m) words)))
        z (case (:pooling m)
            :last (peek hs)
            :mean (ad/scale (reduce ad/add hs) (/ 1.0 (count hs))))]
    (ad/add-row (ad/matmul (ad/tanh (ad/add-row (ad/matmul z (:Wf P)) (:bf P))) (:Wc P)) (:bc P))))

(defn class-loss
  "The cross-entropy of the class of one example, `[words class]`."
  [m P [words class]]
  (ad/cross-entropy (class-logits m P words) [((:class-index m) class)]))

(defn train-classifier-seq
  "The classifier trained on `examples`, each `[words class]`, an example
  a step: `{:model :loss}`."
  [m examples & {:keys [lr] :or {lr 0.01}}]
  (adam-seq m examples class-loss lr nil))

(defn classify
  "The class `m` gives `words`."
  [m words]
  (let [row (first (ad/value (class-logits m (consts (:params m)) words)))]
    ((:classes m) (apply max-key row (range (count row))))))

(defn class-accuracy
  "The fraction of `examples`, `[words class]`, classified right."
  [m examples]
  (/ (double (count (filter (fn [[ws c]] (= c (classify m ws))) examples))) (max 1 (count examples))))

;; ---------------------------------------------------------------------------
;; A bidirectional language model: contextual embeddings

(defn bilm
  "A bidirectional language model after ELMo's (Peters et al., \"Deep
  Contextualized Word Representations\", NAACL 2018): two recurrent
  networks of `cell` over the same embeddings, one predicting each word
  from the words before it and one from the words after, sharing the
  embeddings and the output layer as ELMo does. Its parameters are named
  as a one-layer bidirectional tagger's, so a tagger can start from
  them."
  [cell vocab dim hidden rng]
  (let [v (count (:words vocab))]
    {:cell cell :vocab vocab :hidden hidden :layers 1 :bidirectional true
     :params (merge {:E (ad/init v dim rng) :Wy (ad/init hidden v rng) :by [(vec (repeat v 0.0))]}
                    (cell-params cell dim hidden (direction-prefix 0 :f) rng)
                    (cell-params cell dim hidden (direction-prefix 0 :b) rng))}))

(defn- bilm-states [{:keys [cell hidden vocab]} P words]
  (let [t (word-ids vocab (concat [start] words [end]))
        xs (mapv #(ad/row (:E P) %) t)]
    [t xs
     (run-cells cell (cell-of P cell (direction-prefix 0 :f)) hidden xs)
     (vec (rseq (run-cells cell (cell-of P cell (direction-prefix 0 :b)) hidden (vec (rseq xs)))))]))

(defn bilm-loss
  "The two directions' mean cross-entropies, summed: each word predicted
  from the forward state before it and from the backward state after."
  [m P words]
  (let [[t _ fwd bwd] (bilm-states m P words)
        n (count t)
        out #(ad/add-row (ad/matmul (ad/stack-rows %) (:Wy P)) (:by P))]
    (ad/add (ad/cross-entropy (out (subvec fwd 0 (dec n))) (subvec t 1))
            (ad/cross-entropy (out (subvec bwd 1)) (subvec t 0 (dec n))))))

(defn train-bilm-seq
  "The bidirectional language model trained on unlabeled `sentences`, a
  sentence a step: `{:model :loss}`."
  [m sentences & {:keys [lr] :or {lr 0.01}}]
  (adam-seq m sentences bilm-loss lr nil))

(defn contextual
  "The contextual embedding of each of `words`: its own embedding beside
  the forward state after reading it and the backward state after reading
  it from the other end -- so the same word gets a different vector in
  each sentence, as its company differs. ELMo learns a weighting of its
  layers for each task; here they are kept side by side."
  [m words]
  (let [[_ xs fwd bwd] (bilm-states m (consts (:params m)) words)]
    (mapv (fn [i] (vec (mapcat #(first (ad/value (% i))) [xs fwd bwd])))
          (range 1 (inc (count words))))))

(defn fine-tune
  "A one-layer bidirectional tagger to `tags` started from the
  bidirectional language model `m`: its embeddings and both networks, a
  new output layer."
  [m tags rng]
  (let [{:keys [cell vocab hidden]} m
        fresh (tagger cell vocab tags 1 hidden 1 true rng)]
    (assoc fresh :params (merge (dissoc (:params m) :Wy :by) (select-keys (:params fresh) [:Wt :bt])))))

;; ---------------------------------------------------------------------------
;; Probing frozen features

(defn probe
  "A softmax layer from features of `width` to `tags`, to be trained on
  frozen features -- how much of a tag a representation already holds."
  [width tags rng]
  (let [tags (vec tags)]
    {:tags tags :tag-index (zipmap tags (range))
     :params {:Wt (ad/init width (count tags) rng) :bt [(vec (repeat (count tags) 0.0))]}}))

(defn- probe-loss [m P [features tags]]
  (ad/cross-entropy (ad/add-row (ad/matmul (ad/const features) (:Wt P)) (:bt P)) (mapv (:tag-index m) tags)))

(defn train-probe-seq
  "The probe trained on `examples`, each `[features tags]` for a sentence
  -- a row of features and a tag per word -- a sentence a step."
  [m examples & {:keys [lr] :or {lr 0.01}}]
  (adam-seq m examples probe-loss lr nil))

(defn probe-tags
  "The tags the probe gives the rows of `features`."
  [m features]
  (mapv (fn [row] ((:tags m) (apply max-key row (range (count row)))))
        (ad/value (ad/add-row (ad/matmul (ad/const features) (ad/const (:Wt (:params m)))) (ad/const (:bt (:params m)))))))
