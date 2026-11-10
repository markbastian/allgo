(ns allgo.language.ngrams
  "N-gram language models: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 23 -- the probability of a
  sequence of words, or of characters, as a Markov chain of order n - 1:
  each token's probability given the n - 1 before it, counted from a
  corpus.

    smoothing        the counts leave unseen n-grams at probability zero.
                     Laplace's add-one counts every possible one once more;
                     linear interpolation mixes the n-gram, (n-1)-gram, ...
                     and unigram estimates with weights lambda, which are
                     learned by EM on held-out text -- an iteration a step
    perplexity       P(tokens)^(-1/N): the model's average branching
                     factor, lower the better it predicts unseen text
    unknown words    the infrequent words of the training text replaced
                     by <UNK> (the 3rd edition replaced each word's first
                     appearance instead), so that <UNK> is counted like any
                     word, and an unknown word in a test sentence takes its
                     probability
    sampling         a sequence drawn token by token
    identification   a character model for each language, and a text's
                     language the one making it likeliest, by Bayes' rule

  A sequence is a vector of tokens -- words, or characters -- and is padded
  with n - 1 start markers in front and an end marker behind, so that its
  first tokens have a context and its length is part of what is modeled."
  (:require [clojure.math :as math]))

(def start "<s>")
(def end "</s>")
(def unknown "<UNK>")

(defn- padded [n tokens] (vec (concat (repeat (dec (long n)) start) tokens [end])))

(defn- context-of [v i k] (subvec v (- (long i) (long k)) i))

(defn model
  "An n-gram model of `sequences`: for each order k from 0 to n - 1, how
  often each token follows each k tokens of context. `{:n :vocab
  :counts}`, `:counts` a vector by k of `{context {token count}}`."
  [n sequences]
  (let [n (long n)
        padded (map #(padded n %) sequences)]
    {:n n
     :vocab (disj (set (mapcat identity padded)) start)
     :counts (vec (for [k (range n)]
                    (reduce (fn [m [ctx t]] (update-in m [ctx t] (fnil inc 0)))
                            {}
                            (for [v padded i (range (dec n) (count v))] [(context-of v i k) (v i)]))))}))

(defn- mle
  "The maximum-likelihood estimate of P(token | context): its count over
  the context's; 0 for a context never seen."
  [{:keys [counts]} context token]
  (let [follow (get-in counts [(count context) context])]
    (if follow (/ (double (get follow token 0)) (reduce + (vals follow))) 0.0)))

(defn probability
  "P(`token` | the `context` before it -- the last n - 1 tokens, start
  markers in front if fewer) under `smoothing`:

    :none            the bare counts
    :laplace         add one: (count + 1) / (context count + V), V the
                     size of the vocabulary
    [:interpolate l] the weights `l` mixing the estimates of each order:
                     l0 / V + l1 P(w) + l2 P(w | w-1) + ... -- the
                     first a uniform estimate, below the unigram, so that
                     no token has probability zero"
  [{:keys [n vocab] :as m} smoothing context token]
  (let [n (long n)
        ctx (vec (take-last (dec n) (concat (repeat (dec n) start) context)))]
    (cond
      (= smoothing :none) (mle m ctx token)
      (= smoothing :laplace)
      (let [follow (get-in m [:counts (dec n) ctx])]
        (/ (+ 1.0 (get follow token 0)) (+ (reduce + 0 (vals follow)) (count vocab))))
      (= :interpolate (first smoothing))
      (let [[l0 & ls] (second smoothing)]
        (+ (/ (double l0) (count vocab))
           (reduce + (map-indexed (fn [k l] (* (double l) (mle m (subvec ctx (- (dec n) k)) token))) ls)))))))

(defn- predictions
  "Each token of each sequence, the end included, with its context:
  `[context token]`."
  [n sequences]
  (for [s sequences :let [v (padded n s)] i (range (dec (long n)) (count v))]
    [(context-of v i (dec (long n))) (v i)]))

(defn log-probability
  "The natural log of the probability of `tokens`, the end included."
  [m smoothing tokens]
  (reduce + (for [[ctx t] (predictions (:n m) [tokens])] (math/log (probability m smoothing ctx t)))))

(defn perplexity
  "P(tokens)^(-1/N) over `sequences`, N the tokens predicted -- each
  sequence's end counted as one: e to the mean negative log probability."
  [m smoothing sequences]
  (let [ps (predictions (:n m) sequences)]
    (math/exp (- (/ (reduce + (for [[ctx t] ps] (math/log (probability m smoothing ctx t)))) (count ps))))))

(defn- interpolation-state
  "The weights and the held-out log likelihood they give."
  [table lambdas]
  {:lambdas lambdas
   :log-likelihood (reduce + (map (fn [row] (math/log (reduce + (map * lambdas row)))) table))})

(defn- interpolation-em-step
  "Each weight made its order's share of the held-out tokens."
  [table {:keys [lambdas]}]
  (let [shares (reduce (fn [acc row]
                         (let [parts (map * lambdas row) total (reduce + parts)]
                           (mapv + acc (map #(/ % total) parts))))
                       (vec (repeat (count lambdas) 0.0)) table)]
    (interpolation-state table (mapv #(/ % (count table)) shares))))

(defn interpolation-seq
  "The interpolation weights fitted to `held-out` sequences by EM, an
  iteration a step, from equal weights: the E step asks, of each held-out
  token, how much of its probability each order's estimate supplied; the
  M step makes each weight its order's share. The weights are the
  uniform estimate's, then the unigram's and up. Each state `{:lambdas
  :log-likelihood}`, the held-out log likelihood of the weights it
  holds, which never falls."
  [m held-out]
  (let [n (long (:n m))
        ps (vec (predictions n held-out))
        ;; each held-out token's estimate at each order, the uniform
        ;; first, computed once
        table (mapv (fn [[ctx t]] (into [(/ 1.0 (count (:vocab m)))] (map (fn [k] (mle m (subvec ctx (- (dec n) k)) t)) (range n)))) ps)]
    (iterate (partial interpolation-em-step table) (interpolation-state table (vec (repeat (inc n) (/ 1.0 (inc n))))))))

(defn with-unknown
  "`sequences` with words replaced by `<UNK>`, so that the model learns
  how often unknown words come and where. Given `min-count`, the 4th
  edition's way: every word seen fewer than `min-count` times in all of
  them, the infrequent ones. Without, the 3rd edition's: the first
  appearance of each word."
  ([sequences min-count]
   (let [counts (frequencies (apply concat sequences))]
     (mapv (fn [s] (mapv #(if (< (counts %) min-count) unknown %) s)) sequences)))
  ([sequences]
   (first (reduce (fn [[out seen] s]
                    (let [[s' seen] (reduce (fn [[s' seen] w] (if (seen w) [(conj s' w) seen] [(conj s' unknown) (conj seen w)]))
                                            [[] seen] s)]
                      [(conj out s') seen]))
                  [[] #{}] sequences))))

(defn known
  "`tokens` with each the model has not seen replaced by `<UNK>`."
  [{:keys [vocab]} tokens]
  (mapv #(if (contains? vocab %) % unknown) tokens))

(defn sample
  "A sequence drawn from `m` under `smoothing`, a token at a time from
  its distribution given those before, until the end or `max-len`."
  [{:keys [vocab] :as m} smoothing rng & {:keys [max-len] :or {max-len 30}}]
  (let [choices (sort (map str vocab))]
    (loop [out []]
      (if (>= (count out) (long max-len))
        out
        (let [r (double (rng))
              t (loop [[c & more] choices acc 0.0]
                  (let [acc (+ acc (probability m smoothing out c))]
                    (if (or (nil? more) (< r acc)) c (recur more acc))))]
          (if (= t end) out (recur (conj out t))))))))

(defn identify
  "The language of `text` (a sequence of tokens), of the `models` -- `{lang
  [model smoothing]}` -- the one with the greatest log prior plus log
  likelihood: argmax P(l) P(text | l). `priors` default to equal."
  ([models text] (identify models nil text))
  ([models priors text]
   (key (apply max-key (fn [[l [m smoothing]]]
                         (+ (math/log (get priors l (/ 1.0 (count models)))) (log-probability m smoothing text)))
               (sort-by key models)))))
