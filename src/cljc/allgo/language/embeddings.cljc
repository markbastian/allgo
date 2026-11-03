(ns allgo.language.embeddings
  "Word embeddings: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 24.

  A word is known by the company it keeps: a vector for each word, learned
  so that words used in the same contexts get vectors pointing the same
  way. Skip-gram with negative sampling, word2vec's method: each word is
  trained to predict the words within `window` of it -- the dot product of
  its vector with a true neighbor's context vector pushed up through a
  sigmoid, and with a few words drawn from the unigram distribution to the
  3/4 power, which are probably not neighbors, pushed down. An epoch a
  step.

  What comes out: `similarity`, the cosine of two words' vectors;
  `nearest`, the words most like one; `analogy`, the word d with a : b ::
  c : d, the nearest to b - a + c -- which works when the contexts that
  make b differ from a are the ones that make d differ from c; and
  `project-2d`, the vectors on their two principal components, by
  `allgo.numerics.eigen`'s Jacobi method.

  The gradient of this objective is two lines, so it is followed directly
  rather than through `allgo.learning.autodiff`: an epoch is many thousand
  small updates, and that is where the time goes.

  GloVe, below, gets at the same vectors from the other side: the table
  of how often each pair of words occurs together, fitted whole."
  (:require [allgo.learning.autodiff :as ad]
            [allgo.numerics.eigen :as eigen]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn vocabulary
  "The words of `sentences` (each a vector of words), most frequent first:
  `{:words [...] :index {word i} :counts [...]}`."
  [sentences]
  (let [freqs (sort-by (fn [[w c]] [(- c) (str w)]) (frequencies (apply concat sentences)))]
    {:words (mapv first freqs) :index (zipmap (map first freqs) (range)) :counts (mapv second freqs)}))

(defn pairs
  "Every (word, context) index pair within `window` of each other."
  [sentences {:keys [index]} window]
  (vec (for [s sentences
             :let [ids (mapv index s)]
             i (range (count ids)) j (range (max 0 (- i window)) (min (count ids) (+ i window 1)))
             :when (not= i j)]
         [(ids i) (ids j)])))

(defn- sigmoid [x] (/ 1.0 (+ 1.0 (math/exp (- x)))))
(defn- dot [a b] (reduce + (map * a b)))

(defn- noise-table
  "Word indices in proportion to count^(3/4), to draw negatives from."
  [counts]
  (let [ws (mapv #(math/pow % 0.75) counts) total (reduce + ws)]
    (vec (mapcat (fn [i w] (repeat (max 1 (math/round (* 1000.0 (/ w total)))) i)) (range) ws))))

(defn- update-pair
  "One word, one context, label 1 for a true pair and 0 for noise: both
  vectors moved along the gradient of the log-sigmoid."
  [{:keys [in out] :as model} w c label lr]
  (let [vw (in w) vc (out c)
        g (* lr (- label (sigmoid (dot vw vc))))]
    (-> model
        (assoc-in [:in w] (mapv #(+ %1 (* g %2)) vw vc))
        (assoc-in [:out c] (mapv #(+ %1 (* g %2)) vc vw)))))

(defn sgns-seq
  "Skip-gram with negative sampling on `sentences`: the model after each
  epoch -- or, with `:chunk`, after every that many pairs, for a caller
  that must not stall -- `{:vocab :in :out}`, the word vectors in `:in`.
  Options `:dim` (default 16), `:window` (2), `:negatives` (4), `:lr`
  (0.05)."
  [sentences rng & {:keys [dim window negatives lr chunk] :or {dim 16 window 2 negatives 4 lr 0.05}}]
  (let [vocab (vocabulary sentences)
        n (count (:words vocab))
        ps (pairs sentences vocab window)
        noise (noise-table (:counts vocab))
        learn (fn [model batch]
                (reduce (fn [m [w c]]
                          (reduce (fn [m _] (update-pair m w (random/pick rng noise) 0.0 lr))
                                  (update-pair m w c 1.0 lr)
                                  (range negatives)))
                        model batch))
        epochs (fn epochs [] (lazy-seq (concat (random/shuffle rng ps) (epochs))))]
    (reductions learn
                {:vocab vocab :in (ad/init n dim rng) :out (vec (repeat n (vec (repeat dim 0.0))))}
                (partition (or chunk (count ps)) (epochs)))))

(defn word-vector [{:keys [vocab in]} word] (in ((:index vocab) word)))

(defn cosine [a b]
  (let [na (math/sqrt (dot a a)) nb (math/sqrt (dot b b))]
    (if (or (zero? na) (zero? nb)) 0.0 (/ (dot a b) (* na nb)))))

(defn similarity
  "The cosine of the angle between two words' vectors."
  [model a b]
  (cosine (word-vector model a) (word-vector model b)))

(defn nearest-to
  "The `k` words whose vectors are nearest `v` by cosine, `[word cosine]`,
  leaving out `exclude`."
  [{:keys [vocab in]} v k exclude]
  (->> (map vector (:words vocab) in)
       (remove (fn [[w]] (contains? exclude w)))
       (map (fn [[w u]] [w (cosine u v)]))
       (sort-by (comp - second))
       (take k)))

(defn nearest
  "The `k` words most like `word`."
  [model word k]
  (nearest-to model (word-vector model word) k #{word}))

(defn analogy
  "The word d for a : b :: c : d -- the nearest to b - a + c, a b and c
  left out."
  [model a b c]
  (let [[va vb vc] (map #(word-vector model %) [a b c])]
    (ffirst (nearest-to model (mapv #(+ (- %2 %1) %3) va vb vc) 1 #{a b c}))))

(defn project-2d
  "`vectors` on their two principal components: centered, then onto the
  covariance's two leading eigenvectors."
  [vectors]
  (let [n (count vectors) d (count (first vectors))
        mean (mapv #(/ % n) (apply mapv + vectors))
        X (mapv #(mapv - % mean) vectors)
        C (vec (for [i (range d)] (vec (for [j (range d)] (/ (reduce + (map #(* (% i) (% j)) X)) n)))))
        {:keys [vectors]} (eigen/jacobi C)
        pc (fn [k] (mapv #(nth % k) vectors))
        [p1 p2] [(pc (dec d)) (pc (- d 2))]]
    (mapv (fn [x] [(dot x p1) (dot x p2)]) X)))

;; ---------------------------------------------------------------------------
;; GloVe

(defn cooccurrences
  "How much each pair of words occurs together: for every two words
  within `window` of each other, d apart, 1 / d added to their count --
  near neighbors counting more, as GloVe weighs them. `{[i j] x}`, by
  index, both ways round."
  [sentences {:keys [index]} window]
  (reduce (fn [m [i j d]] (update m [i j] (fnil + 0.0) (/ 1.0 (double d))))
          {}
          (for [s sentences
                :let [ids (mapv index s)]
                a (range (count ids))
                b (range (max 0 (- a window)) (min (count ids) (+ a window 1)))
                :when (not= a b)]
            [(ids a) (ids b) (abs (- a b))])))

(defn- glove-weight
  "f(x) = (x / x-max)^alpha, 1 from x-max up: rare pairs, whose counts
  are noise, count for less, and the commonest no more than the rest."
  [x x-max alpha]
  (if (< x x-max) (math/pow (/ x x-max) alpha) 1.0))

(defn- adagrad
  "`v` moved against gradient `g`, each component by `lr` over the root
  of its own running sum of squared gradients, `acc`: `[v' acc']`."
  [v acc g lr]
  (let [acc (mapv #(+ %1 (* %2 %2)) acc g)]
    [(mapv (fn [x gi a] (- x (/ (* lr gi) (math/sqrt a)))) v g acc) acc]))

(defn glove-seq
  "GloVe (Pennington, Socher and Manning, \"GloVe: Global Vectors for Word
  Representation\", EMNLP 2014): rather than predict each word's
  neighbors a pair at a time, fit the whole table of how often words
  occur together. Each pair's word vector w_i and context vector w~_j,
  with biases, are fitted so that w_i . w~_j + b_i + b~_j is the log of
  their count x_ij, by weighted least squares -- the weight `f(x)` --
  over every pair that occurs at all, with AdaGrad as the paper does.
  The model after each epoch, the pairs taken in a new random order each
  time: `{:vocab :in :w :wt :cost}`, `:in` the vectors the similarities
  and analogies use and `:cost` the epoch's weighted squared error. The
  paper takes w + w~ as a word's vector, finding the sum a little better
  than either; `:vectors :sum` does that. The default, `:word`, is w
  alone: where words keep to their roles -- a person never a context word
  about a person -- the sum mixes in what each is as a context, and on
  such corpora the analogies suffered. With `:chunk`, a state after every
  that many pairs rather than every epoch, `:cost` theirs, for a caller
  that must not stall. Options `:dim` (default 16),
  `:window` (2), `:x-max` (100), `:alpha` (0.75) and `:lr` (0.05)."
  [sentences rng & {:keys [dim window x-max alpha lr vectors chunk] :or {dim 16 window 2 x-max 100.0 alpha 0.75 lr 0.05 vectors :word}}]
  (let [vocab (vocabulary sentences)
        n (count (:words vocab))
        entries (vec (cooccurrences sentences vocab window))
        ones (fn [] (vec (repeat n (vec (repeat dim 1.0)))))
        with-in (fn [{:keys [w wt] :as m}] (assoc m :in (if (= vectors :sum) (mapv #(mapv + %1 %2) w wt) w)))
        start {:vocab vocab :w (ad/init n dim rng) :wt (ad/init n dim rng)
               :b (vec (repeat n 0.0)) :bt (vec (repeat n 0.0))
               ;; AdaGrad's sums start at 1, as the paper's code does
               :gw (ones) :gwt (ones) :gb (vec (repeat n 1.0)) :gbt (vec (repeat n 1.0))}
        fit (fn [m batch]
              (with-in
                (reduce (fn [{:keys [w wt b bt gw gwt gb gbt] :as m} [[i j] x]]
                          (let [wi (w i) wj (wt j)
                                diff (- (+ (dot wi wj) (b i) (bt j)) (math/log x))
                                fd (* (glove-weight x x-max alpha) diff)
                                [wi' gwi] (adagrad wi (gw i) (mapv #(* fd %) wj) lr)
                                [wj' gwj] (adagrad wj (gwt j) (mapv #(* fd %) wi) lr)
                                [[bi'] [gbi]] (adagrad [(b i)] [(gb i)] [fd] lr)
                                [[bj'] [gbj]] (adagrad [(bt j)] [(gbt j)] [fd] lr)]
                            (-> m
                                (assoc-in [:w i] wi') (assoc-in [:gw i] gwi)
                                (assoc-in [:wt j] wj') (assoc-in [:gwt j] gwj)
                                (assoc-in [:b i] bi') (assoc-in [:gb i] gbi)
                                (assoc-in [:bt j] bj') (assoc-in [:gbt j] gbj)
                                (update :cost + (* 0.5 fd diff)))))
                        (assoc m :cost 0.0)
                        batch)))
        epochs (fn epochs [] (lazy-seq (concat (random/shuffle rng entries) (epochs))))]
    (reductions fit (assoc (with-in start) :cost nil) (partition (or chunk (count entries)) (epochs)))))
