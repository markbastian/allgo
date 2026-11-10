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
  (:require [allgo.language.recurrent :as rec]
            [allgo.language.seq2seq :as s2s]
            [allgo.language.tagging :as tagging]
            [allgo.learning.autodiff :as ad]
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

(defn- layer-params
  "One layer's parameters, named with `prefix`: attention's four maps,
  the feed-forward network's two layers, and a gain and bias for each
  normalization -- with `cross`, a second attention, to another
  sequence, and its normalization."
  [prefix d ff cross rng]
  (into {}
        (map (fn [[k m]] [(keyword (str prefix (name k))) m]))
        (merge {:Wq (ad/init d d rng) :Wk (ad/init d d rng) :Wv (ad/init d d rng) :Wo (ad/init d d rng)
                :g1 (ones-row d) :b1 (zeros-row d)
                :W1 (ad/init d ff rng) :c1 (zeros-row ff) :W2 (ad/init ff d rng) :c2 (zeros-row d)
                :g2 (ones-row d) :b2 (zeros-row d)}
               (when cross
                 {:xWq (ad/init d d rng) :xWk (ad/init d d rng) :xWv (ad/init d d rng) :xWo (ad/init d d rng)
                  :g3 (ones-row d) :b3 (zeros-row d)}))))

(defn encoder
  "The parameters of an encoder over `vocab`: embeddings of `d`, `layers`
  layers of `heads` heads and a feed-forward width of `ff`. Positions
  are told by the sinusoidal encoding, or, given `:max-len`, by a
  learned embedding for each position up to it -- the chapter's way."
  [vocab d heads ff layers rng & {:keys [max-len]}]
  {:d d :heads heads :layers layers :vocab vocab
   :params (merge {:E (ad/init (count (:words vocab)) d rng)}
                  (when max-len {:Pos (ad/init max-len d rng)})
                  (apply merge (for [l (range layers)] (layer-params (str "l" l "-") d ff false rng))))})

(defn- p [P prefix k] (P (keyword (str prefix (name k)))))

(defn- attention
  "Multi-head attention from the rows of `Xq` to those of `Xkv`, by the
  maps named `[q k v o]` under `prefix`: `[output weights]`, the weights
  an n x m matrix per head. `mask`, if given, is added to the scores
  before the softmax: -1e9 where a row may not look."
  [P prefix [q k v o] Xq Xkv heads mask]
  (let [d (count (first (ad/value Xq))) dk (quot d heads)
        Q (ad/matmul Xq (p P prefix q))
        [K V] (map #(ad/matmul Xkv (p P prefix %)) [k v])
        outs (for [h (range heads)]
               (let [sl #(ad/slice-cols % (* h dk) (* (inc h) dk))
                     scores (ad/scale (ad/matmul (sl Q) (ad/transpose (sl K))) (/ 1.0 (math/sqrt dk)))
                     A (ad/softmax-rows (if mask (ad/add scores mask) scores))]
                 [(ad/matmul A (sl V)) A]))]
    [(ad/matmul (ad/concat-cols (map first outs)) (p P prefix o)) (mapv (comp ad/value second) outs)]))

(defn- add-norm
  "The residual `X + Y`, normalized, with the gain and bias `g` and `b`."
  [P prefix X Y g b]
  (ad/add-row (ad/mul-row (ad/layer-norm-rows (ad/add X Y)) (p P prefix g)) (p P prefix b)))

(defn- feed-forward [P prefix X]
  (ad/add-row (ad/matmul (ad/relu (ad/add-row (ad/matmul X (p P prefix :W1)) (p P prefix :c1))) (p P prefix :W2)) (p P prefix :c2)))

(defn- embed
  "The embeddings of `ids` plus their positions'."
  [d P ids]
  (ad/add (ad/rows (:E P) ids)
          (if (:Pos P)
            (ad/rows (:Pos P) (range (count ids)))
            (ad/const (positional-encoding (count ids) d)))))

(defn encode
  "The encoder run on token indices `ids`: `[X attention]`, X the last
  layer's vectors (n x d) and the attention weights by layer and head."
  [{:keys [d heads layers]} P ids]
  (reduce (fn [[X att] l]
            (let [pre (str "l" l "-")
                  [a w] (attention P pre [:Wq :Wk :Wv :Wo] X X heads nil)
                  X (add-norm P pre X a :g1 :b1)
                  X (add-norm P pre X (feed-forward P pre X) :g2 :b2)]
              [X (conj att w)]))
          [(embed d P ids) []] (range layers)))

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

(defn tag
  "The tags `m` gives `words`, and its attention: `{:tags :attention}`."
  [m words]
  (let [[logits att] (tag-logits m (ad/consts (:params m)) words)]
    {:tags (mapv (fn [row] ((:tags m) (apply max-key row (range (count row))))) (ad/value logits))
     :attention att}))

(defn train-tagger-seq
  "The tagger trained on `tagged` sentences, each `[[word tag] ...]`, a
  sentence a step: each state `{:model :loss}`. With `:train-only`, a set
  of parameter names, the rest stay as they are -- `#{:Wt :bt}` trains
  just the output layer on a pretrained encoder's representations."
  [m tagged & {:keys [lr train-only] :or {lr 0.01}}]
  (ad/adam-seq m tagged tag-loss :lr lr :only train-only))

(defn tag-accuracy
  "The fraction of words `m` tags right, over `tagged` sentences, or over
  just the words `only?` picks out."
  ([m tagged] (tagging/accuracy-of #(:tags (tag m %)) tagged))
  ([m tagged only?] (tagging/accuracy-of #(:tags (tag m %)) tagged (fn [w _] (only? w)))))

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
  (ad/adam-seq m sentences mlm-loss :lr lr :prepare (fn [words] (into [words] (mask-words words rate rng)))))

(defn fill-in
  "The word `m` puts at each `<mask>` in `masked`."
  [m masked]
  (let [positions (vec (keep-indexed (fn [i w] (when (= w mask-token) i)) masked))
        P (ad/consts (:params m))
        [X] (encode m P (ids m masked))
        logits (ad/value (ad/add-row (ad/matmul (ad/rows X positions) (:Wm P)) (:bm P)))]
    (mapv (fn [row] ((:words (:vocab m)) (apply max-key row (range (count row))))) logits)))

(defn fine-tune
  "A tagger started from a pretrained model's encoder."
  [pretrained tags rng]
  (tagger (update pretrained :params dissoc :Wm :bm) tags rng))

;; ---------------------------------------------------------------------------
;; The encoder-decoder

(defn encoder-decoder
  "The whole transformer, encoder and decoder, over `vocab` -- a
  `seq2seq/vocabulary`, its markers included -- for source and target
  alike: embeddings of `d` shared by both, `layers` layers each of
  `heads` heads and a feed-forward width of `ff`, and an output layer to
  the vocabulary. Each decoder layer has the encoder's parts and a second
  attention, to the encoder's output."
  [vocab d heads ff layers rng]
  (let [enc (encoder vocab d heads ff layers rng)
        v (count (:words vocab))]
    (update enc :params merge
            (apply merge (for [l (range layers)] (layer-params (str "d" l "-") d ff true rng)))
            {:Wout (ad/init d v rng) :bout (zeros-row v)})))

(defn causal-mask
  "The n x n mask that lets each position attend only to itself and the
  positions before it: 0 on and below the diagonal, -1e9 above."
  [n]
  (vec (for [i (range n)] (vec (for [j (range n)] (if (> j i) -1e9 0.0))))))

(defn- decode
  "The decoder run on target indices `ids`, the encoder's output `H`
  beside it: `[logits cross]`, a row of logits per position and the
  attention to the source by layer and head."
  [{:keys [d heads layers]} P H ids]
  (let [mask (ad/const (causal-mask (count ids)))
        [X cross] (reduce (fn [[X cross] l]
                            (let [pre (str "d" l "-")
                                  [a] (attention P pre [:Wq :Wk :Wv :Wo] X X heads mask)
                                  X (add-norm P pre X a :g1 :b1)
                                  [c w] (attention P pre [:xWq :xWk :xWv :xWo] X H heads nil)
                                  X (add-norm P pre X c :g3 :b3)
                                  X (add-norm P pre X (feed-forward P pre X) :g2 :b2)]
                              [X (conj cross w)]))
                          [(embed d P ids) []] (range layers))]
    [(ad/add-row (ad/matmul X (:Wout P)) (:bout P)) cross]))

(defn translation-loss
  "The mean cross-entropy of `target` given `source`, teacher forced: the
  decoder reads `<s>` and the target, and at each position must say the
  token after it -- the mask keeping it from reading ahead."
  [{:keys [vocab] :as m} P [source target]]
  (let [index (:index vocab)
        tgt (mapv index (concat [rec/start] target [rec/end]))
        [H] (encode m P (mapv index source))
        [logits] (decode m P H (pop tgt))]
    (ad/cross-entropy logits (subvec tgt 1))))

(defn train-translator-seq
  "The encoder-decoder trained on `pairs` of `[source target]`, a pair a
  step by Adam: `{:model :loss}`. The rate, by default 0.005, is half
  the taggers': at theirs, reversing strings, it never learned."
  [m pairs & {:keys [lr] :or {lr 0.005}}]
  (ad/adam-seq m pairs translation-loss :lr lr))

(defn- next-distribution
  "The distribution over the token after `prefix` (indices, `<s>`
  first), and the attention to the source so far."
  [m P H prefix]
  (let [[logits cross] (decode m P H prefix)]
    [(ad/softmax (peek (ad/value logits))) cross]))

(defn translate
  "The target `m` writes for `source`: greedily, the likeliest token each
  step -- `{:target :attention}`, the attention the last layer's to the
  source, averaged over the heads, a row per token written -- or, given
  `:beam`, by a beam search of that width, `{:target}`. At most
  `max-len` tokens, by default twice the source's plus 2."
  [{:keys [vocab params] :as m} source & {:keys [beam max-len]}]
  (let [{:keys [index words]} vocab
        P (ad/consts params)
        [H] (encode m P (mapv index source))
        max-len (or max-len (+ 2 (* 2 (count source))))]
    (if beam
      {:target (s2s/beam-search vocab #(first (next-distribution m P H %)) beam max-len)}
      (loop [prefix [(index rec/start)]]
        (let [[p cross] (next-distribution m P H prefix)
              nxt (apply max-key p (range (count p)))]
          (if (or (= nxt (index rec/end)) (> (count prefix) max-len))
            {:target (mapv words (rest prefix))
             :attention (let [heads (peek cross)]
                          (apply mapv (fn [& rows] (apply mapv (fn [& xs] (/ (reduce + xs) (count xs))) rows)) heads))}
            (recur (conj prefix nxt))))))))
