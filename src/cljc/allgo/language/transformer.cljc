(ns allgo.language.transformer
  "The transformer: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 24.

  No recurrence: every token looks at every other at once. In
  self-attention each token's vector makes a query, a key and a value by
  three learned maps; its new vector is the values averaged with weights
  from the softmax of its query's dot products with every key, over the
  square root of their size. Several heads do this side by side on
  slices of the vectors, each free to attend to something different,
  and their outputs are joined and mapped back. A layer adds that to its
  input (a residual connection) and normalizes, then does the same with a
  small feed-forward network applied to each position alone. Since
  attention has no notion of order, each token's embedding is added to a
  positional encoding: sines and cosines of its position at wavelengths
  in geometric progression.

  Two uses here, as the chapter's:

    tagging          a label for every token -- a part of speech -- read
                     off each position of the last layer
    masked language  pretraining without labels: some tokens hidden
    modeling         behind a `<mask>`, and the encoder trained to say
                     what they were from the rest; its encoder can then
                     start a tagger -- `fine-tune` -- the chapter's
                     transfer learning. (On the small grammar here, where
                     twenty labeled sentences already tag perfectly,
                     pretraining measurably learns the parts of speech
                     but does not make a better tagger; NOTES has the
                     numbers.)

  Differentiated by `allgo.learning.autodiff`, trained by Adam, a
  sentence a step. `encode` also returns every head's attention weights,
  to see what each looked at."
  (:require [allgo.learning.autodiff :as ad]
            [clojure.math :as math]))

(def mask-token "<mask>")

(defn vocabulary
  "The words of `sentences`, `<mask>` first."
  [sentences]
  (let [ws (vec (cons mask-token (sort (distinct (apply concat sentences)))))]
    {:words ws :index (zipmap ws (range))}))

(defn positional-encoding
  "The sinusoidal positional encoding of `n` positions in `d` dimensions:
  PE(p, 2i) = sin(p / 10000^(2i/d)), PE(p, 2i+1) = cos(the same)."
  [n d]
  (vec (for [p (range n)]
         (vec (for [j (range d)]
                (let [angle (/ p (math/pow 10000.0 (/ (* 2 (quot j 2)) d)))]
                  (if (even? j) (math/sin angle) (math/cos angle))))))))

(defn- zeros-row [n] [(vec (repeat n 0.0))])
(defn- ones-row [n] [(vec (repeat n 1.0))])

(defn encoder
  "The parameters of an encoder over `vocab`: embeddings of `d`, `layers`
  layers of `heads` heads and a feed-forward width of `ff`."
  [vocab d heads ff layers rng]
  {:d d :heads heads :layers layers :vocab vocab
   :params (into {:E (ad/init (count (:words vocab)) d rng)}
                 (for [l (range layers)
                       [k m] {:Wq (ad/init d d rng) :Wk (ad/init d d rng) :Wv (ad/init d d rng) :Wo (ad/init d d rng)
                              :g1 (ones-row d) :b1 (zeros-row d)
                              :W1 (ad/init d ff rng) :c1 (zeros-row ff) :W2 (ad/init ff d rng) :c2 (zeros-row d)
                              :g2 (ones-row d) :b2 (zeros-row d)}]
                   [(keyword (str "l" l "-" (name k))) m]))})

(defn- p [P l k] (P (keyword (str "l" l "-" (name k)))))

(defn- self-attention
  "Multi-head self-attention of `X` (n x d): `[output weights]`, the
  weights one n x n matrix per head."
  [P l X heads]
  (let [d (count (first (ad/value X))) dk (quot d heads)
        [Q K V] (map #(ad/matmul X (p P l %)) [:Wq :Wk :Wv])
        outs (for [h (range heads)]
               (let [sl #(ad/slice-cols % (* h dk) (* (inc h) dk))
                     A (ad/softmax-rows (ad/scale (ad/matmul (sl Q) (ad/transpose (sl K))) (/ 1.0 (math/sqrt dk))))]
                 [(ad/matmul A (sl V)) A]))]
    [(ad/matmul (ad/concat-cols (map first outs)) (p P l :Wo)) (mapv (comp ad/value second) outs)]))

(defn encode
  "The encoder run on token indices `ids`: `[X attention]`, X the last
  layer's vectors (n x d) and the attention weights by layer and head."
  [{:keys [d heads layers]} P ids]
  (let [X0 (ad/add (ad/rows (:E P) ids) (ad/const (positional-encoding (count ids) d)))]
    (reduce (fn [[X att] l]
              (let [[a w] (self-attention P l X heads)
                    X (ad/add-row (ad/mul-row (ad/layer-norm-rows (ad/add X a)) (p P l :g1)) (p P l :b1))
                    f (ad/add-row (ad/matmul (ad/relu (ad/add-row (ad/matmul X (p P l :W1)) (p P l :c1))) (p P l :W2)) (p P l :c2))
                    X (ad/add-row (ad/mul-row (ad/layer-norm-rows (ad/add X f)) (p P l :g2)) (p P l :b2))]
                [X (conj att w)]))
            [X0 []] (range layers))))

(defn- ids [{:keys [vocab]} words] (mapv #(get (:index vocab) % 0) words))

;; ---------------------------------------------------------------------------
;; Tagging

(defn tagger
  "A tagger: `encoder` (fresh, or a pretrained one's) and an output layer
  to `tags`."
  [enc tags rng]
  (let [tags (vec tags)]
    (-> enc
        (assoc :tags tags :tag-index (zipmap tags (range)))
        (update :params merge {:Wt (ad/init (:d enc) (count tags) rng) :bt (zeros-row (count tags))}))))

(defn- tag-logits [m P words]
  (let [[X att] (encode m P (ids m words))]
    [(ad/add-row (ad/matmul X (:Wt P)) (:bt P)) att]))

(defn tag-loss [m P tagged]
  (let [[logits] (tag-logits m P (mapv first tagged))]
    (ad/cross-entropy logits (mapv #((:tag-index m) (second %)) tagged))))

(defn- consts [params] (into {} (map (fn [[k v]] [k (ad/const v)]) params)))

(defn tag
  "The tags `m` gives `words`, and its attention: `{:tags :attention}`."
  [m words]
  (let [[logits att] (tag-logits m (consts (:params m)) words)]
    {:tags (mapv (fn [row] ((:tags m) (apply max-key row (range (count row))))) (ad/value logits))
     :attention att}))

(defn- train
  "Adam on `loss` over `examples` in turn; only the parameters in
  `only` (default all) moved -- the rest frozen."
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
  "The tagger trained on `tagged` sentences, each `[[word tag] ...]`, a
  sentence a step: each state `{:model :loss}`. With `:train-only`, a set
  of parameter names, the rest stay as they are -- `#{:Wt :bt}` trains
  just the output layer on a pretrained encoder's representations."
  [m tagged & {:keys [lr train-only] :or {lr 0.01}}]
  (train m tagged tag-loss lr train-only))

(defn tag-accuracy
  "The fraction of words `m` tags right, over `tagged` sentences, or over
  just the words `only?` picks out."
  ([m tagged] (tag-accuracy m tagged (constantly true)))
  ([m tagged only?]
   (let [results (for [s tagged
                       :let [got (:tags (tag m (mapv first s)))]
                       [[w t] g] (map vector s got)
                       :when (only? w)]
                   (= t g))]
     (/ (double (count (filter true? results))) (max 1 (count results))))))

(defn most-frequent-tag
  "The baseline: each word given the tag it most often has in `tagged`,
  unknown words the commonest tag. A function of a word."
  [tagged]
  (let [by-word (group-by first (apply concat tagged))
        best (into {} (map (fn [[w ps]] [w (key (apply max-key val (frequencies (map second ps))))]) by-word))
        common (key (apply max-key val (frequencies (map second (apply concat tagged)))))]
    (fn [w] (get best w common))))

;; ---------------------------------------------------------------------------
;; Masked language modeling

(defn masked-lm
  "A masked language model: `encoder` and an output layer back to the
  words."
  [enc rng]
  (update enc :params merge {:Wm (ad/init (:d enc) (count (:words (:vocab enc))) rng)
                             :bm (zeros-row (count (:words (:vocab enc))))}))

(defn mask-words
  "`words` with each hidden behind `<mask>` with probability `rate` (at
  least one): `[masked positions]`."
  [words rate rng]
  (let [chosen (vec (keep-indexed (fn [i _] (when (< (double (rng)) rate) i)) words))
        chosen (if (seq chosen) chosen [(long (* (rng) (count words)))])]
    [(reduce #(assoc %1 %2 mask-token) (vec words) chosen) chosen]))

(defn mlm-loss
  "The cross-entropy of saying the hidden words, `example` being `[words
  masked positions]`."
  [m P [words masked positions]]
  (let [[X] (encode m P (ids m masked))
        logits (ad/add-row (ad/matmul (ad/rows X positions) (:Wm P)) (:bm P))]
    (ad/cross-entropy logits (mapv #((:index (:vocab m)) (words %)) positions))))

(defn train-mlm-seq
  "The masked language model trained on `sentences`, each masked afresh
  at `rate` (default 0.2) every time it comes round: `{:model :loss}`."
  [m sentences rng & {:keys [lr rate] :or {lr 0.01 rate 0.2}}]
  (let [sentences (vec sentences)]
    (map first
         (iterate (fn [[_ m opt k]]
                    (let [words (sentences (mod k (count sentences)))
                          [masked positions] (mask-words words rate rng)
                          {:keys [value grads]} (ad/value-and-grad #(mlm-loss m % [words masked positions]) (:params m))
                          [params opt] (ad/adam-step (:params m) (ad/clip grads 5.0) opt :lr lr)
                          m (assoc m :params params)]
                      [{:model m :loss value} m opt (inc k)]))
                  [{:model m :loss nil} m nil 0]))))

(defn fill-in
  "The word `m` puts at each `<mask>` in `masked`."
  [m masked]
  (let [positions (vec (keep-indexed (fn [i w] (when (= w mask-token) i)) masked))
        P (consts (:params m))
        [X] (encode m P (ids m masked))
        logits (ad/value (ad/add-row (ad/matmul (ad/rows X positions) (:Wm P)) (:bm P)))]
    (mapv (fn [row] ((:words (:vocab m)) (apply max-key row (range (count row))))) logits)))

(defn fine-tune
  "A tagger started from a pretrained model's encoder."
  [pretrained tags rng]
  (tagger (update pretrained :params dissoc :Wm :bm) tags rng))
