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
  after each."
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
