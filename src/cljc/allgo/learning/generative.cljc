(ns allgo.learning.generative
  "Unsupervised deep learning: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 21.7.1 -- models
  that learn a representation z of data x, or a distribution new x can be
  sampled from, or both.

    probabilistic PCA     z ~ N(0, I), x ~ N(W z, sigma^2 I), so P(x) =
                          N(0, W W^T + sigma^2 I), equation 21.16; W and
                          sigma^2 learned by EM (here about the data's
                          mean, where the chapter's model has mean zero),
                          z^ = E[z | x] the features, a low P(x) an
                          anomaly
    linear autoencoder    z^ = W x, x = W^T z^, W trained to minimize the
                          squared reconstruction error -- it comes to span
                          the data's m principal components
    deep autoencoder      any network with a narrow middle layer, trained
                          to give back its input: the layers to the
                          middle the encoder f, those after the decoder g
    variational           an encoder giving Q(z), a Gaussian with diagonal
    autoencoder           covariance, a decoder giving the mean of P(x |
                          z), and both trained to maximize the variational
                          lower bound L = H(Q) + E_{z~Q} log P(z, x): H(Q)
                          analytic, the expectation estimated by sampling
                          z from Q
    autoregressive        each x_t predicted from x_{t-k} ... x_{t-1}: the
                          classical linear-Gaussian AR model of order k,
                          fitted by the Yule-Walker equations, or a deep
                          one, a network in place of the linear model; an
                          n-gram model is one of order n - 1
    GAN                   a generator h from z ~ N(0, I) to x and a
                          discriminator classifying x as real or the
                          generator's, trained at once, the generator to
                          fool the discriminator and the discriminator not
                          to be fooled

  Data are vectors of numbers, a data set a vector of them."
  (:require [allgo.learning.autodiff :as ad]
            [allgo.learning.neural :as nn]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn- mean-of [X] (let [n (count X)] (mapv #(/ % n) (apply mapv + X))))

(defn covariance
  "The data's covariance matrix about their mean, over n."
  [X]
  (let [mu (mean-of X) C (mapv #(lin/sub % mu) X)]
    (lin/mat-scale (lin/mat-mul (lin/transpose C) C) (/ 1.0 (count X)))))

(defn- trace [M] (reduce + (map-indexed (fn [i row] (row i)) M)))

(defn- gaussians [rng n] (vec (repeatedly n #(random/gaussian rng 0.0 1.0))))

(defn- gaussian-log-density
  "log N(x; mu, C), C symmetric positive definite."
  [mu C x]
  (let [L (lin/cholesky C) r (lin/sub x mu)]
    (* -0.5 (+ (* (count x) (math/log (* 2.0 math/PI)))
               (* 2.0 (reduce + (map-indexed (fn [i row] (math/log (row i))) L)))
               (lin/dot r (lin/cholesky-solve C r))))))

;; ---------------------------------------------------------------------------
;; Probabilistic PCA

(defn- ppca-covariance [{:keys [W var]}]
  (lin/mat-add (lin/mat-mul W (lin/transpose W)) (lin/mat-scale (lin/eye (count W)) var)))

(defn ppca-em-step
  "One EM step for probabilistic PCA, given the data's covariance `S`,
  in the form Tipping and Bishop (\"Probabilistic Principal Component
  Analysis\", J. R. Statist. Soc. B 61, 1999) give it, the E step folded
  into the M step's sums: with M = W^T W + sigma^2 I,

    W'        = S W (sigma^2 I + M^-1 W^T S W)^-1
    sigma^2'  = tr(S - S W M^-1 W'^T) / d"
  [S {:keys [W var] :as model}]
  (let [m (count (first W)) d (count W)
        Wt (lin/transpose W)
        Minv (ls/inverse (lin/mat-add (lin/mat-mul Wt W) (lin/mat-scale (lin/eye m) var)))
        SW (lin/mat-mul S W)
        W' (lin/mat-mul SW (ls/inverse (lin/mat-add (lin/mat-scale (lin/eye m) var) (lin/mat-mul Minv (lin/mat-mul Wt SW)))))]
    (assoc model
           :W W'
           :var (/ (trace (lin/mat-sub S (lin/mat-mul SW (lin/mat-mul Minv (lin/transpose W'))))) d))))

(defn ppca-em-seq
  "EM for a probabilistic PCA model of `X` with `m` latent dimensions, a
  step at a time: `{:mean :W :var}`, the mean the data's, W `d` by `m`
  and random to start, sigma^2 1."
  [X m rng]
  (let [S (covariance X) d (count S)]
    (iterate (partial ppca-em-step S)
             {:mean (mean-of X) :W (vec (repeatedly d (fn [] (gaussians rng m)))) :var 1.0})))

(defn ppca-log-density
  "log P(x) under the model, equation 21.16: N(x; mean, W W^T + sigma^2
  I). Data the model finds very improbable are potential anomalies."
  [model x]
  (gaussian-log-density (:mean model) (ppca-covariance model) x))

(defn ppca-log-likelihood [model X] (reduce + (map #(ppca-log-density model %) X)))

(defn ppca-encode
  "The features of `x`: z^ = E[z | x] = M^-1 W^T (x - mean)."
  [{:keys [W var mean]} x]
  (let [Wt (lin/transpose W)
        M (lin/mat-add (lin/mat-mul Wt W) (lin/mat-scale (lin/eye (count Wt)) var))]
    (lin/mat-vec (ls/inverse M) (lin/mat-vec Wt (lin/sub x mean)))))

(defn ppca-sample
  "A new x: z from its prior N(0, I), then x from N(W z, sigma^2 I),
  about the mean."
  [{:keys [W var mean]} rng]
  (let [z (gaussians rng (count (first W)))]
    (mapv (fn [mu x e] (+ mu x (* (math/sqrt var) e))) mean (lin/mat-vec W z) (gaussians rng (count W)))))

;; ---------------------------------------------------------------------------
;; Autoencoders

(defn- linear-reconstruction-loss
  "The mean over the rows of X of || x - W^T W x ||^2, data a row each:
  Z = X W^T, X^ = Z W -- W used twice, its gradient the sum of both
  uses'."
  [{:keys [W]} X]
  (let [R (ad/sub X (ad/matmul (ad/matmul X (ad/transpose W)) W))]
    (ad/scale (ad/sum (ad/mul R R)) (/ 1.0 (count (ad/value X))))))

(defn linear-autoencoder-seq
  "A linear autoencoder for `X` with an `m`-dimensional z, trained by
  Adam on the whole data set, a step at a time: `{:mean :W :loss}`, W
  `m` by `d`, the data centered on their mean first, as PCA does."
  [X m rng & {:keys [lr] :or {lr 0.01}}]
  (let [mu (mean-of X)
        Xc (ad/const (mapv #(lin/sub % mu) X))]
    (map #(assoc (:params %) :mean mu :loss (:loss %))
         (iterate (fn [{:keys [params opt]}]
                    (let [{:keys [value grads]} (ad/value-and-grad #(linear-reconstruction-loss % Xc) params)
                          [params opt] (ad/adam-step params grads opt :lr lr)]
                      {:params params :opt opt :loss value}))
                  {:params {:W (ad/init m (count (first X)) rng)} :opt nil}))))

(defn linear-encode [{:keys [W mean]} x] (lin/mat-vec W (lin/sub x mean)))
(defn linear-decode [{:keys [W mean]} z] (lin/add mean (lin/mat-vec (lin/transpose W) z)))

(defn autoencoder
  "A deep autoencoder from `d` dimensions to `m` and back, through a tanh
  layer of `hidden` each way: the encoder its first two layers, the
  decoder its last two."
  [d hidden m rng]
  [(nn/dense d hidden :tanh rng) (nn/dense hidden m :identity rng)
   (nn/dense m hidden :tanh rng) (nn/dense hidden d :identity rng)])

(defn autoencoder-seq
  "The autoencoder trained to give back `X`, an epoch a step: `train-seq`
  with each example its own target."
  [network X opts]
  (nn/train-seq network X X opts))

(defn encode [network x] (nn/predict (subvec network 0 (quot (count network) 2)) x))
(defn decode [network z] (nn/predict (subvec network (quot (count network) 2)) z))

(defn reconstruction-error
  "The mean squared distance between each x and `reconstruct` of it."
  [reconstruct X]
  (/ (reduce + (map #(lin/distance-squared % (reconstruct %)) X)) (count X)))

;; ---------------------------------------------------------------------------
;; Variational autoencoders

(defn vae
  "A variational autoencoder from `d` dimensions to `k` latent ones: an
  encoder through a tanh layer of `hidden` to the mean and log variance
  of Q(z), and a decoder through a tanh layer of `decoder-hidden` -- or,
  with 0, a linear one -- to the mean of P(x | z), whose variance is
  `noise-var`. Row vectors, as `allgo.learning.autodiff` has them."
  [d hidden k decoder-hidden noise-var rng]
  (let [row (fn [n] [(vec (repeat n 0.0))])]
    {:noise-var noise-var
     :params (merge {:We (ad/init d hidden rng) :be (row hidden)
                     :Wm (ad/init hidden k rng) :bm (row k)
                     :Wv (ad/init hidden k rng) :bv (row k)}
                    (if (pos? decoder-hidden)
                      {:Wd (ad/init k decoder-hidden rng) :bd (row decoder-hidden)
                       :Wo (ad/init decoder-hidden d rng) :bo (row d)}
                      {:Wo (ad/init k d rng) :bo (row d)}))}))

(defn- affine [P W b X] (ad/add-row (ad/matmul X (P W)) (P b)))

(defn- q-nodes
  "Q(z) for each row of X: the means and log variances."
  [P X]
  (let [H (ad/tanh (affine P :We :be X))] [(affine P :Wm :bm H) (affine P :Wv :bv H)]))

(defn- decoder-node [P Z]
  (if (P :Wd) (affine P :Wo :bo (ad/tanh (affine P :Wd :bd Z))) (affine P :Wo :bo Z)))

(defn- elbo-node
  "The variational lower bound summed over the rows of X, L = H(Q) +
  E_{z~Q} log P(z, x), with one z drawn from Q per row as mu + sigma
  eps, `E` the eps -- so that the sample is a function of the encoder's
  outputs and the gradient reaches them (Kingma and Welling's
  reparameterization). `[node constant]`, the bound the node's value
  plus the constant, which the gradient does not need:

    H(Q)      = sum over z's units of (1 + log 2 pi + log sigma^2) / 2
    log P(z)  = -(|z|^2 + k log 2 pi) / 2
    log P(x | z) = -(|x - g(z)|^2 / noise + d log (2 pi noise)) / 2"
  [P noise-var X E]
  (let [[Mu Lv] (q-nodes P X)
        Z (ad/add Mu (ad/mul (ad/exp (ad/scale Lv 0.5)) E))
        R (ad/sub X (decoder-node P Z))
        n (count (ad/value X)) k (count (first (ad/value Mu))) d (count (first (ad/value X)))]
    [(ad/add (ad/add (ad/scale (ad/sum Lv) 0.5) (ad/scale (ad/sum (ad/mul Z Z)) -0.5))
             (ad/scale (ad/sum (ad/mul R R)) (/ -0.5 noise-var)))
     (* n (- (* 0.5 k) (* 0.5 d (math/log (* 2.0 math/PI noise-var)))))]))

(defn- noise-rows [rng n k] (vec (repeatedly n #(gaussians rng k))))

(defn vae-step
  "One step of Adam on a minibatch `X`, both networks moved up the
  gradient of the variational lower bound."
  [{:keys [params noise-var opt] :as model} X rng lr]
  (let [k (count (first (:bm params)))
        E (ad/const (noise-rows rng (count X) k))
        {:keys [grads]} (ad/value-and-grad #(ad/scale (first (elbo-node % noise-var (ad/const X) E)) (/ -1.0 (count X))) params)
        [params opt] (ad/adam-step params grads opt :lr lr)]
    (assoc model :params params :opt opt)))

(defn vae-seq
  "The variational autoencoder trained on `X`, an epoch a step, the data
  shuffled and taken `batch` at a time."
  [model X rng & {:keys [batch lr] :or {batch 32 lr 0.005}}]
  (iterate (fn [m] (reduce #(vae-step %1 (vec %2) rng lr) m (partition-all batch (random/shuffle rng X))))
           model))

(defn elbo
  "The variational lower bound on log P(x), estimated with `samples`
  draws of z from Q."
  [{:keys [params noise-var]} x rng samples]
  (let [k (count (first (:bm params)))
        X (ad/const (vec (repeat samples x)))
        [node c] (elbo-node (into {} (map (fn [[n m]] [n (ad/const m)]) params)) noise-var X (ad/const (noise-rows rng samples k)))]
    (/ (+ (ffirst (ad/value node)) c) samples)))

(defn vae-encode
  "Q(z) for `x`: `{:mean :var}`."
  [{:keys [params]} x]
  (let [P (into {} (map (fn [[n m]] [n (ad/const m)]) params))
        [Mu Lv] (q-nodes P (ad/const [x]))]
    {:mean (first (ad/value Mu)) :var (mapv math/exp (first (ad/value Lv)))}))

(defn vae-decode
  "The mean of P(x | z)."
  [{:keys [params]} z]
  (first (ad/value (decoder-node (into {} (map (fn [[n m]] [n (ad/const m)]) params)) (ad/const [z])))))

(defn vae-sample
  "A new x: z from its prior N(0, I), decoded -- the mean of P(x | z)."
  [{:keys [params] :as model} rng]
  (vae-decode model (gaussians rng (count (first (:bm params))))))

(defn vae-log-likelihood
  "With a linear decoder, x = z Wo + bo plus noise, the exact log P(x) the
  bound is below: N(x; bo, Wo^T Wo + noise I)."
  [{:keys [params noise-var]} x]
  (let [{:keys [Wo bo]} params]
    (gaussian-log-density (first bo) (lin/mat-add (lin/mat-mul (lin/transpose Wo) Wo) (lin/mat-scale (lin/eye (count x)) noise-var)) x)))

;; ---------------------------------------------------------------------------
;; Autoregressive models

(defn autocovariance
  "The series' autocovariance at `lag`, about its mean, over n."
  [xs lag]
  (let [n (count xs) mu (/ (reduce + xs) n) v (vec xs)]
    (/ (reduce + (map (fn [t] (* (- (v t) mu) (- (v (+ t lag)) mu))) (range (- n lag)))) n)))

(defn yule-walker
  "The classical AR model of order `k`, x_t = mean + sum over j of phi_j
  (x_{t-j} - mean) + noise of variance sigma^2, fitted by the Yule-Walker
  equations: R phi = r, R the k by k matrix of autocovariances r_|i-j|,
  r those at lags 1 to k; sigma^2 = r_0 - phi . r. `{:mean :phi :var}`."
  [xs k]
  (let [r (mapv #(autocovariance xs %) (range (inc k)))
        R (mapv (fn [i] (mapv (fn [j] (r (abs (- i j)))) (range k))) (range k))
        phi (ls/gauss R (subvec r 1))]
    {:mean (/ (reduce + xs) (count xs)) :phi phi :var (- (r 0) (lin/dot phi (subvec r 1)))}))

(defn ar-predict
  "The AR model's mean for the next value after `window`, the last `k`
  values, oldest first."
  [{:keys [mean phi]} window]
  (+ mean (reduce + (map (fn [p x] (* p (- x mean))) phi (rseq (vec window))))))

(defn windows
  "The series as examples for a deep AR model of order `k`: `[X Y]`, each
  x the `k` values before, its y the next value."
  [xs k]
  (let [v (vec xs)]
    [(mapv #(subvec v % (+ % k)) (range (- (count v) k)))
     (mapv #(vector (v (+ % k))) (range (- (count v) k)))]))

(defn deep-ar-predict [network window] (first (nn/predict network window)))

(defn ar-generate-seq
  "A new series from an AR model, a value a step: from `init`, the first
  `k` values, each next one `predict` of the `k` before it -- `(partial
  ar-predict model)` or `(partial deep-ar-predict network)` -- plus, with
  `sd` and `rng`, Gaussian noise of that standard deviation."
  ([predict init] (ar-generate-seq predict init 0.0 nil))
  ([predict init sd rng]
   (concat init
           (map peek (rest (iterate (fn [w] (conj (subvec w 1) (+ (predict w) (if rng (random/gaussian rng 0.0 sd) 0.0))))
                                    (vec init)))))))

;; ---------------------------------------------------------------------------
;; Generative adversarial networks

(defn gan
  "A generator from `kz` dimensions of noise to `d`, and a discriminator
  from `d` to the logit that its input is real, each through a layer of
  `hidden`: tanh in the generator, relu in the discriminator -- with tanh
  there it saturated, and missed more of what the generator got wrong."
  [kz d hidden rng]
  (let [row (fn [n] [(vec (repeat n 0.0))])]
    {:kz kz
     :generator {:G1 (ad/init kz hidden rng) :g1 (row hidden) :G2 (ad/init hidden d rng) :g2 (row d)}
     :discriminator {:D1 (ad/init d hidden rng) :d1 (row hidden) :D2 (ad/init hidden 1 rng) :d2 (row 1)}}))

(defn- consts [P] (into {} (map (fn [[n m]] [n (ad/const m)]) P)))

(defn- generator-node [G Z] (affine G :G2 :g2 (ad/tanh (affine G :G1 :g1 Z))))

(defn- real-or-fake
  "The discriminator's logits for the rows of X as the two classes fake
  and real, `[0 logit]`, so softmax gives P(real) = sigmoid(logit)."
  [D X]
  (let [L (affine D :D2 :d2 (ad/relu (affine D :D1 :d1 X)))]
    (ad/concat-cols [(ad/const (vec (repeat (count (ad/value X)) [0.0]))) L])))

(defn gan-step
  "One round of the game on a minibatch of real data `X`: the
  discriminator a step of Adam toward classifying `X` real and as many
  of the generator's samples fake, then the generator a step toward the
  discriminator's calling its samples real (minimizing -log D(h(z)),
  which pushes hardest where the discriminator is surest)."
  [{:keys [kz generator discriminator d-opt g-opt] :as model} X rng lr]
  (let [n (count X)
        fake (ad/value (generator-node (consts generator) (ad/const (noise-rows rng n kz))))
        {dg :grads} (ad/value-and-grad #(ad/cross-entropy (real-or-fake % (ad/const (into (vec X) fake)))
                                                          (into (vec (repeat n 1)) (repeat n 0)))
                                       discriminator)
        [discriminator d-opt] (ad/adam-step discriminator dg d-opt :lr lr :b1 0.5)
        Z (ad/const (noise-rows rng n kz))
        D (consts discriminator)
        {gg :grads} (ad/value-and-grad #(ad/cross-entropy (real-or-fake D (generator-node % Z)) (vec (repeat n 1))) generator)
        [generator g-opt] (ad/adam-step generator gg g-opt :lr lr :b1 0.5)]
    (assoc model :generator generator :discriminator discriminator :d-opt d-opt :g-opt g-opt)))

(defn gan-seq
  "The game played on `X`, a minibatch of `batch` drawn at random a step."
  [model X rng & {:keys [batch lr] :or {batch 64 lr 0.002}}]
  (let [v (vec X)]
    (iterate #(gan-step % (vec (repeatedly batch (fn [] (random/pick rng v)))) rng lr) model)))

(defn gan-sample
  "`n` new x from the generator."
  [{:keys [kz generator]} rng n]
  (ad/value (generator-node (consts generator) (ad/const (noise-rows rng n kz)))))

(defn discriminate
  "The discriminator's probability that each of `X` is real."
  [{:keys [discriminator]} X]
  (mapv #(second (ad/softmax %)) (ad/value (real-or-fake (consts discriminator) (ad/const (vec X))))))
