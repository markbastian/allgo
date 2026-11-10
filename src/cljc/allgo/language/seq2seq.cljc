(ns allgo.language.seq2seq
  "Sequence-to-sequence models with attention: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 24.

  An encoder reads the source a token at a time, a recurrent network
  (`allgo.language.recurrent`'s cells); a decoder, started from the
  encoder's last state, writes the target a token at a time, each step
  fed the token before. Squeezing a whole source into one vector is the
  bottleneck the chapter describes, so the decoder also attends: at each
  step it scores every encoder state against its own by their dot
  product, turns the scores into weights by a softmax, and reads the
  weighted sum of encoder states -- the context -- beside its own state
  when choosing the next token. The weights say where it looked, and are
  returned.

  Decoding: greedy, the likeliest token each step; or beam search, which
  is `allgo.search`'s beam search on prefixes, each step's cost the
  negative log probability of the token added -- best first, the frontier
  cut to the beam's width, the first complete sentence off it the answer.

  Tokens are any values; the vocabulary adds `<s>` and `</s>`."
  (:require [allgo.language.recurrent :as rec]
            [allgo.learning.autodiff :as ad]
            [allgo.search :as search]
            [clojure.math :as math]
            [clojure.string :as str]))

(defn vocabulary [tokens]
  (let [ws (vec (concat [rec/start rec/end] (sort-by str (distinct tokens))))]
    {:words ws :index (zipmap ws (range))}))

(defn model
  "An encoder-decoder over `vocab` with `cell` (`:rnn` or `:lstm`),
  embeddings of `dim` and states of `hidden`. With `:attention false`,
  the basic model of the chapter's Figure 24.6: the decoder knows the
  source only by the encoder's last state, which it starts from."
  [cell vocab dim hidden rng & {:keys [attention] :or {attention true}}]
  (let [v (count (:words vocab))
        enc (:params (rec/model cell vocab dim hidden rng))
        dec (:params (rec/model cell vocab dim hidden rng))]
    {:cell cell :vocab vocab :hidden hidden :attention attention
     :params (merge (into {} (map (fn [[k m]] [(keyword (str "enc-" (name k))) m]) (dissoc enc :Wy :by)))
                    (into {} (map (fn [[k m]] [(keyword (str "dec-" (name k))) m]) (dissoc dec :Wy :by)))
                    {:Wo (ad/init (if attention (* 2 hidden) hidden) v rng) :bo [(vec (repeat v 0.0))]})}))

(defn- sub-params [P prefix]
  (into {} (keep (fn [[k v]] (let [n (name k)] (when (str/starts-with? n prefix) [(keyword (subs n (count prefix))) v]))) P)))

(defn- encode
  "The encoder's states, stacked, and its last state."
  [{:keys [cell hidden]} P source-ids]
  (let [E (sub-params P "enc-")
        [hs st] (reduce (fn [[hs st] t] (let [[h c] (rec/step cell E (ad/row (:E E) t) st)] [(conj hs h) [h c]]))
                        [[] (rec/initial-state hidden)] source-ids)]
    [(ad/stack-rows hs) st]))

(defn- decode-step
  "From decoder state `st` on token `t`: `[logits weights st']`, the
  weights nil without attention."
  [{:keys [cell attention]} P H st t]
  (let [D (sub-params P "dec-")
        [h c] (rec/step cell D (ad/row (:E D) t) st)]
    (if attention
      ;; the encoder states scored against h, softmaxed
      (let [weights (ad/softmax-rows (ad/transpose (ad/matmul H (ad/transpose h))))
            context (ad/matmul weights H)]
        [(ad/add-row (ad/matmul (ad/concat-cols [h context]) (:Wo P)) (:bo P)) weights [h c]])
      [(ad/add-row (ad/matmul h (:Wo P)) (:bo P)) nil [h c]])))

(defn- ids [{:keys [index]} tokens] (mapv index tokens))

(defn loss
  "The mean cross-entropy of the target given the source, teacher forced:
  each decoder step fed the true token before."
  [{:keys [vocab] :as m} P source target]
  (let [src (ids vocab source)
        tgt (ids vocab (concat [rec/start] target [rec/end]))
        [H st] (encode m P src)
        [outs] (reduce (fn [[outs st] t] (let [[l _ st] (decode-step m P H st t)] [(conj outs l) st]))
                       [[] st] (pop tgt))]
    (ad/cross-entropy (ad/stack-rows outs) (subvec tgt 1))))

(defn train-seq
  "The model trained on `pairs` of `[source target]`, a pair a step by
  Adam: each state `{:model :loss}`."
  [m pairs & {:keys [lr clip] :or {lr 0.01 clip 5.0}}]
  (ad/adam-seq m pairs (fn [m P [s t]] (loss m P s t)) :lr lr :clip clip))

(defn- prefix-distribution
  "After the decoder has read `prefix` (target ids, `<s>` first): the
  distribution over the next token, and the attention weights of each
  step so far."
  [m P H st0 prefix]
  (let [[_ ws p] (reduce (fn [[st ws _] t]
                           (let [[l w st] (decode-step m P H st t)]
                             [st (cond-> ws w (conj (first (ad/value w)))) (ad/softmax (first (ad/value l)))]))
                         [st0 [] nil] prefix)]
    [p ws]))

(defn greedy
  "The target the model writes for `source`, the likeliest token each
  step, at most `max-len` (default twice the source plus 2): `{:target
  :attention}`, the attention a row of weights per output token."
  [{:keys [vocab params] :as m} source & {:keys [max-len]}]
  (let [P (ad/consts params) [H st] (encode m P (ids vocab source))
        max-len (or max-len (+ 2 (* 2 (count source))))
        end ((:index vocab) rec/end)]
    (loop [prefix [((:index vocab) rec/start)]]
      (let [[p ws] (prefix-distribution m P H st prefix)
            nxt (apply max-key p (range (count p)))]
        (if (or (= nxt end) (> (count prefix) max-len))
          {:target (mapv (:words vocab) (rest prefix)) :attention ws}
          (recur (conj prefix nxt)))))))

(defn beam-search
  "The target by beam search of `width`: `allgo.search/beam-seq` over
  prefixes of token indices from `<s>`, `(dist prefix)` the probability
  of each next token, a step's cost minus the log probability of the
  token it adds, a prefix ending `</s>` -- or `max-len` long -- the goal.
  The target's tokens, the markers removed."
  [vocab dist width max-len]
  (let [end ((:index vocab) rec/end)
        dist (memoize dist)
        found (search/path (search/goal-state
                            (search/beam-seq {:start [((:index vocab) rec/start)]
                                              :beam-width width
                                              :neighbors (fn [prefix] (map #(conj prefix %) (range (count (dist prefix)))))
                                              :cost (fn [a b] (- (math/log (max 1e-300 ((dist a) (peek b))))))
                                              :goal? (fn [prefix] (or (= end (peek prefix)) (> (count prefix) max-len)))})))]
    (when found
      (mapv (:words vocab) (remove #{end} (rest (peek found)))))))

(defn beam
  "The target `m` writes for `source` by a beam search of `width`, at
  most `max-len` long (default twice the source plus 2)."
  [{:keys [vocab params] :as m} source width & {:keys [max-len]}]
  (let [P (ad/consts params) [H st] (encode m P (ids vocab source))]
    (beam-search vocab #(first (prefix-distribution m P H st %)) width
                 (or max-len (+ 2 (* 2 (count source)))))))
