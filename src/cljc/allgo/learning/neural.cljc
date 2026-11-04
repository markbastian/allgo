(ns allgo.learning.neural
  "Neural networks: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 21 -- BACK-PROP-LEARNING (figure
  21's back-propagation; the 4th edition gives it in the text), with
  dense layers, convolutional and pooling ones, residual layers and
  batch normalization.

    dense         in_j = sum over i of w_ij a_i + b_j, a_j = g(in_j)
    convolution   2-D kernels slid over a grid, stride 1, no padding: the
                  layers of a convolutional network, flattened into the
                  dense ones after
    pooling       section 21.3.1: each window of l by l of every feature
                  map, every s units, summarized by its maximum -- a
                  feature present somewhere in the window -- or its
                  average, which is a convolution with the uniform kernel
                  [1/l ... 1/l]; with s = l, the maps downsampled by s
    residual      section 21.3.3, equation 21.10: z' = g_r(z + f(z)), the
                  residual f(z) = V g(W z) perturbing the representation
                  rather than replacing it; with V = 0 and relus, the
                  layer passes relu outputs through unchanged
    batch norm    section 21.4.2: each unit's values over a minibatch
                  standardized, then scaled by a learned gamma and shifted
                  by a learned beta, z^ = gamma (z - mu) / sqrt(eps +
                  sigma^2) + beta
    training      stochastic gradient descent, an epoch a step, the
                  examples shuffled each epoch and taken a minibatch of
                  `:batch` at a time (default 1), the weights moved
                  against the minibatch's mean gradient; the deltas
                  propagated back from the output exactly as the figure
                  does
    weight decay  section 21.5.3: the penalty lambda sum W_ij^2 added to
                  the loss (`:lambda`), the weights' (and kernels'), not
                  the biases' -- MAP learning with a zero-mean Gaussian
                  prior on the weights
    dropout       section 21.5.4: each step a thinned network, every
                  hidden unit (`:keep p`) and input unit (`:keep-input`)
                  kept with probability p and its output multiplied by
                  1/p, or else zero; at test time no dropout

  A network is a vector of layers:

    {:type :dense :W [[...]] :b [...] :g :sigmoid|:tanh|:relu|:identity}
    {:type :conv :kernels [[[...]]] :b [...] :g ... :shape [h w]}
    {:type :pool :op :max|:average :size l :stride s :shape [h w] :channels c}
    {:type :residual :W [[...]] :bw [...] :V [[...]] :b [...] :inner g :g g_r}
    {:type :batch-norm :gamma [...] :beta [...] :mean [...] :var [...] :g ...}

  a conv layer taking a flat input it reads as an `h` by `w` grid and
  giving its feature maps flattened, one after another, and a pooling
  layer reading `c` maps of `h` by `w` so. The loss is the squared error,
  or with a sigmoid output, `:cross-entropy`."
  (:require [allgo.random :as random]
            [clojure.math :as math]))

(defn- g [kind z]
  (case kind
    :sigmoid (/ 1.0 (+ 1.0 (math/exp (- z))))
    :tanh (math/tanh z)
    :relu (max 0.0 z)
    :identity z))

(defn- g' [kind z a]
  (case kind
    :sigmoid (* a (- 1.0 a))
    :tanh (- 1.0 (* a a))
    :relu (if (pos? z) 1.0 0.0)
    :identity 1.0))

(defn- dot [u v] (reduce + (map * u v)))

;; ---------------------------------------------------------------------------
;; Layers

(defn dense
  "A dense layer from `n-in` to `n-out`, weights small and random."
  [n-in n-out g rng]
  (let [scale (/ 1.0 (math/sqrt n-in))]
    {:type :dense :g g
     :W (vec (repeatedly n-out (fn [] (vec (repeatedly n-in #(* scale (- (* 2.0 (rng)) 1.0)))))))
     :b (vec (repeat n-out 0.0))}))

(defn conv
  "A convolution layer of `n` kernels of `k` by `k` over an input grid of
  `shape`."
  [shape k n g rng]
  {:type :conv :g g :shape shape
   :kernels (vec (repeatedly n (fn [] (vec (repeatedly k (fn [] (vec (repeatedly k #(* (/ 1.0 k) (- (* 2.0 (rng)) 1.0))))))))))
   :b (vec (repeat n 0.0))})

(defn pool
  "A pooling layer over `channels` feature maps of `shape`, as a conv
  layer gives them: each `size` by `size` window, every `stride` units,
  summarized by its maximum (`:max`) or its mean (`:average`) -- `size`
  and `stride` numbers, or `[rows cols]`. Nothing in it is learned."
  [shape channels size stride op]
  {:type :pool :op op :shape shape :channels channels :size size :stride stride :g :identity})

(defn residual
  "A residual layer of `n` units whose residual goes through `k` hidden
  ones: z' = relu(z + V relu(W z)), relus for both g and g_r as the
  chapter's argument has them."
  [n k rng]
  (let [{W :W bw :b} (dense n k :relu rng)
        {V :W b :b} (dense k n :relu rng)]
    {:type :residual :W W :bw bw :V V :b b :inner :relu :g :relu}))

(defn batch-norm
  "Batch normalization of `n` units, then `g`: gamma 1 and beta 0 to
  start. In training each unit is standardized by its minibatch's mean
  and variance; at test time there is no minibatch, and it is
  standardized by running averages of the minibatches' (each step
  `rho` of the old one and 1 - `rho` of the minibatch's) -- the usual
  way, which the chapter does not go into."
  [n g]
  {:type :batch-norm :g g :eps 1e-5 :rho 0.9
   :gamma (vec (repeat n 1.0)) :beta (vec (repeat n 0.0))
   :mean (vec (repeat n 0.0)) :var (vec (repeat n 1.0))})

(defn- conv-out-shape [{:keys [shape kernels]}]
  (let [[h w] shape k (count (first kernels))] [(inc (- h k)) (inc (- w k))]))

(defn- rows-cols [x] (if (vector? x) x [x x]))

(defn- pool-out-shape [{[h w] :shape :keys [size stride]}]
  (let [[lh lw] (rows-cols size) [sh sw] (rows-cols stride)]
    [(inc (quot (- h lh) sh)) (inc (quot (- w lw) sw))]))

(defn layer-out-size [{:keys [type W kernels channels b gamma] :as layer}]
  (case type
    :dense (count W)
    :conv (let [[oh ow] (conv-out-shape layer)] (* (count kernels) oh ow))
    :pool (let [[oh ow] (pool-out-shape layer)] (* channels oh ow))
    :residual (count b)
    :batch-norm (count gamma)))

(defn- layer-in-size [{:keys [type W shape channels] :as layer}]
  (case type
    :dense (count (first W))
    :conv (apply * shape)
    :pool (* channels (apply * shape))
    (:residual :batch-norm) (layer-out-size layer)))

(defn- pool-windows
  "For each of a pooling layer's outputs, the indices of the inputs in
  its window."
  [{[h w] :shape :keys [channels size stride] :as layer}]
  (let [[oh ow] (pool-out-shape layer) [lh lw] (rows-cols size) [sh sw] (rows-cols stride)]
    (for [c (range channels) r (range oh) q (range ow)]
      (vec (for [i (range lh) j (range lw)]
             (+ (* c h w) (* (+ (* r sh) i) w) (* q sw) j))))))

(defn- dot-window
  "The kernel row `row` dotted with the `k` inputs from `base` on."
  [row a base k]
  (loop [j 0 acc 0.0] (if (= j k) acc (recur (inc j) (+ acc (* (row j) (a (+ base j))))))))

(defn- scatter-window!
  "`df` times kernel row `row` added into the transient `acc` at the `k`
  inputs from `base` on."
  [acc row df base k]
  (loop [j 0 acc acc]
    (if (= j k) acc (recur (inc j) (assoc! acc (+ base j) (+ (nth acc (+ base j)) (* df (row j))))))))

(defn- standardize
  "Batch normalization of one example `a` by means `mu` and variances
  `var`."
  [{:keys [gamma beta eps]} a mu var]
  (mapv (fn [z gm bt m v] (+ bt (/ (* gm (- z m)) (math/sqrt (+ eps v))))) a gamma beta mu var))

(defn- pre-activation
  "in, before the activation function -- for a batch-norm layer, as at
  test time."
  [{:keys [type W b kernels shape] :as layer} a]
  (case type
    :dense (mapv (fn [row bj] (+ bj (reduce + (map * row a)))) W b)
    :conv (let [[_ w] shape [oh ow] (conv-out-shape layer) k (count (first kernels))]
            ;; loops rather than nested `for` and `get-in`: this is the
            ;; inner loop of every convolutional network's training
            (persistent!
             (reduce (fn [out [f r c]]
                       (let [K (kernels f)]
                         (conj! out (loop [i 0 acc (b f)]
                                      (if (= i k)
                                        acc
                                        (recur (inc i) (+ acc (dot-window (K i) a (+ (* (+ r i) w) c) k))))))))
                     (transient [])
                     (for [f (range (count kernels)) r (range oh) c (range ow)] [f r c]))))
    :pool (let [op (:op layer)]
            (mapv (fn [idx] (let [xs (map a idx)] (if (= op :max) (reduce max xs) (/ (reduce + xs) (count idx)))))
                  (pool-windows layer)))
    :residual (let [{:keys [bw V inner]} layer
                    h (mapv (fn [row bj] (g inner (+ bj (dot row a)))) W bw)]
                (mapv (fn [ai row bi] (+ ai bi (dot row h))) a V b))
    :batch-norm (standardize layer a (:mean layer) (:var layer))))

(defn forward
  "Every layer's `[in a]`, input first: `[[nil x] [in1 a1] ...]`. With
  `masks`, a vector or nil a layer, each layer's activations are
  multiplied by its mask: dropout's thinned network."
  ([network x] (forward network x nil))
  ([network x masks]
   (reductions (fn [[_ a] [layer m]]
                 (let [z (pre-activation layer a) a (mapv #(g (:g layer) %) z)]
                   [z (if m (mapv * a m) a)]))
               [nil (vec x)] (map vector network (or masks (repeat nil))))))

(defn predict [network x] (second (last (forward network x))))

(defn- batch-statistics
  "Each unit's mean and variance over the minibatch `A`: `[mu var]`."
  [A]
  (let [m (count A)
        mu (mapv #(/ % m) (apply mapv + A))]
    [mu (mapv #(/ % m) (apply mapv + (map (fn [a] (mapv #(let [d (- %1 %2)] (* d d)) a mu)) A)))]))

(defn forward-batch
  "`forward` over a minibatch `X` as training sees it, every layer's
  `[ins as]`: a batch-norm layer standardizing by the minibatch's own
  statistics."
  [network X masks]
  (reductions (fn [[_ A] [layer m]]
                (let [Z (if (= :batch-norm (:type layer))
                          (let [[mu var] (batch-statistics A)] (mapv #(standardize layer % mu var) A))
                          (mapv #(pre-activation layer %) A))]
                  [Z (mapv (fn [z] (let [a (mapv #(g (:g layer) %) z)] (if m (mapv * a m) a))) Z)]))
              [nil (vec X)] (map vector network (or masks (repeat nil)))))

;; ---------------------------------------------------------------------------
;; Back-propagation

(defn- layer-back
  "The gradients of `layer` given the deltas at its output and the
  activations `a` it took in, and the deltas at its input (before the
  previous layer's activation derivative)."
  [{:keys [type W kernels shape] :as layer} delta a]
  (case type
    :dense {:dW (mapv (fn [dj] (mapv #(* dj %) a)) delta)
            :db delta
            :back (mapv (fn [i] (reduce + (map (fn [row dj] (* (row i) dj)) W delta))) (range (count a)))}
    :conv (let [[_ w] shape [oh ow] (conv-out-shape layer) k (count (first kernels))
                d (fn [f r c] (delta (+ (* f oh ow) (* r ow) c)))
                cells (for [r (range oh) c (range ow)] [r c])
                sum-over (fn [g] (loop [[[r c] & more] cells acc 0.0] (if r (recur more (+ acc (g r c))) acc)))]
            {:dK (mapv (fn [f]
                         (mapv (fn [i] (mapv (fn [j] (sum-over (fn [r c] (* (d f r c) (a (+ (* (+ r i) w) c j)))))) (range k)))
                               (range k)))
                       (range (count kernels)))
             :db (mapv (fn [f] (sum-over (fn [r c] (d f r c)))) (range (count kernels)))
             ;; each input's share of the deltas, accumulated in place
             :back (persistent!
                    (reduce (fn [acc [f r c]]
                              (let [df (d f r c) K (kernels f)]
                                (if (zero? df)
                                  acc
                                  (loop [i 0 acc acc]
                                    (if (= i k)
                                      acc
                                      (recur (inc i) (scatter-window! acc (K i) df (+ (* (+ r i) w) c) k)))))))
                            (transient (vec (repeat (count a) 0.0)))
                            (for [f (range (count kernels)) r (range oh) c (range ow)] [f r c])))})
    ;; a window's maximum passes its delta to the input that was largest
    ;; (the first, on a tie); its average, an equal share to each
    :pool (let [op (:op layer)]
            {:back (persistent!
                    (reduce (fn [acc [idx d]]
                              (if (= op :max)
                                (let [i (reduce #(if (> (a %2) (a %1)) %2 %1) idx)] (assoc! acc i (+ (nth acc i) d)))
                                (let [share (/ d (count idx))] (reduce #(assoc! %1 %2 (+ (nth %1 %2) share)) acc idx))))
                            (transient (vec (repeat (count a) 0.0)))
                            (map vector (pool-windows layer) delta)))})
    :residual (let [{:keys [bw V inner]} layer
                    u (mapv (fn [row bj] (+ bj (dot row a))) W bw)
                    h (mapv #(g inner %) u)
                    du (mapv (fn [i ui hi] (* (g' inner ui hi) (reduce + (map (fn [row dj] (* (row i) dj)) V delta))))
                             (range) u h)]
                {:dV (mapv (fn [dj] (mapv #(* dj %) h)) delta)
                 :db delta
                 :dW (mapv (fn [di] (mapv #(* di %) a)) du)
                 :dbw du
                 ;; the identity's delta, and the residual's
                 :back (mapv (fn [i dj] (+ dj (reduce + (map (fn [row di] (* (row i) di)) W du)))) (range) delta)})))

(defn- batch-norm-back
  "A batch-norm layer's gradients over a whole minibatch, its deltas `D`
  and inputs `A` an example a row -- every example's output depending on
  every other's through the mean and variance -- with its minibatch's
  statistics, for the running averages."
  [{:keys [gamma eps]} D A]
  (let [m (count A)
        [mu var] (batch-statistics A)
        sd (mapv #(math/sqrt (+ eps %)) var)
        Xh (mapv (fn [a] (mapv #(/ (- %1 %2) %3) a mu sd)) A)
        md (mapv #(/ % m) (apply mapv + D))
        mdx (mapv #(/ % m) (apply mapv + (map #(mapv * %1 %2) D Xh)))]
    {:dgamma (apply mapv + (map #(mapv * %1 %2) D Xh))
     :dbeta (apply mapv + D)
     :batch-mean mu :batch-var var
     :back (mapv (fn [d xh] (mapv (fn [di xi gm s mdi mdxi] (* (/ gm s) (- di mdi (* xi mdxi)))) d xh gamma sd md mdx)) D Xh)}))

(defn- tree-map [f a b] (if (vector? a) (mapv #(tree-map f %1 %2) a b) (f a b)))

(defn- add-gradients [a b] (merge-with #(tree-map + %1 %2) a b))

(defn batch-gradients
  "The gradient of the mean loss over the minibatch `X` `Y` with respect
  to every weight, a layer at a time: the output deltas, then each
  layer's deltas from the next one's, as BACK-PROP-LEARNING propagates
  them, summed over the examples -- through a batch-norm layer, over the
  whole minibatch at once. With dropout's `masks` (`dropout-masks`), the
  gradient of the thinned network's: a dropped unit passes back nothing,
  a kept one its delta times its mask."
  ([network X Y loss] (batch-gradients network X Y loss nil))
  ([network X Y loss masks]
   (let [acts (vec (forward-batch network X masks))
         n (count network)
         m (count X)
         [Z-out A-out] (acts n)
         out-layer (network (dec n))
         ;; With a sigmoid output and cross-entropy the g' cancels.
         ce? (and (= loss :cross-entropy) (= :sigmoid (:g out-layer)))
         D (mapv (fn [z a y]
                   (mapv #(/ % m) (if ce?
                                    (mapv - a y)
                                    (mapv (fn [zi ai ti] (* 2.0 (- ai ti) (g' (:g out-layer) zi ai))) z a y))))
                 Z-out A-out Y)]
     (loop [l (dec n) D D grads ()]
       (let [layer (network l)
             A-in (second (acts l))
             {:keys [back] :as gr} (if (= :batch-norm (:type layer))
                                     (batch-norm-back layer D A-in)
                                     (let [gs (mapv #(layer-back layer %1 %2) D A-in)]
                                       (assoc (reduce add-gradients (map #(dissoc % :back) gs)) :back (mapv :back gs))))
             grads (cons (dissoc gr :back) grads)]
         (if (zero? l)
           (vec grads)
           (let [[Z A] (acts l) gl (:g (network (dec l))) mk (get masks (dec l))]
             (recur (dec l)
                    (mapv (fn [b z a]
                            (if mk
                              ;; a = m g(in): the unmasked g(in) for g'
                              (mapv (fn [bi zi mi] (if (zero? mi) 0.0 (* bi mi (g' gl zi (g gl zi))))) b z mk)
                              (mapv (fn [bi zi ai] (* bi (g' gl zi ai))) b z a)))
                          back Z A)
                    grads))))))))

(defn gradients
  "The gradient of the loss on one example: `batch-gradients` of a
  minibatch of one."
  ([network x y loss] (batch-gradients network [x] [y] loss nil))
  ([network x y loss masks] (batch-gradients network [x] [y] loss masks)))

(def ^:private parameter-of
  {:dW :W :dK :kernels :dV :V :db :b :dbw :bw :dgamma :gamma :dbeta :beta})

(def ^:private weights #{:W :kernels :V})

(defn- apply-gradients
  "Every parameter moved against its gradient by `alpha` -- a weight's
  with weight decay's 2 lambda w added -- and a batch-norm layer's
  running mean and variance moved toward its minibatch's."
  [network grads alpha lambda]
  (mapv (fn [layer gr]
          (reduce-kv (fn [layer gk G]
                       (if-let [pk (parameter-of gk)]
                         (update layer pk #(tree-map (if (and lambda (weights pk))
                                                       (fn [p gi] (- p (* alpha (+ gi (* 2.0 lambda p)))))
                                                       (fn [p gi] (- p (* alpha gi))))
                                                     % G))
                         layer))
                     (if-let [mu (:batch-mean gr)]
                       (let [rho (:rho layer) run (fn [old new] (mapv #(+ (* rho %1) (* (- 1.0 rho) %2)) old new))]
                         (-> layer (update :mean run mu) (update :var run (:batch-var gr))))
                       layer)
                     gr))
        network grads))

(defn- example-loss [a y kind]
  (if (= kind :cross-entropy)
    (- (reduce + (map (fn [ai yi] (+ (* yi (math/log (max 1e-12 ai))) (* (- 1 yi) (math/log (max 1e-12 (- 1 ai)))))) a y)))
    (reduce + (map #(let [d (- %1 %2)] (* d d)) a y))))

(defn loss
  "The mean loss of `network` on examples `X` `Y`."
  [network X Y kind]
  (/ (reduce + (map (fn [x y] (example-loss (predict network x) y kind)) X Y))
     (count X)))

(defn batch-loss
  "The mean loss on the minibatch `X` `Y` as training sees it, batch
  normalization by the minibatch's own statistics: what
  `batch-gradients` is the gradient of."
  [network X Y kind]
  (/ (reduce + (map #(example-loss %1 %2 kind) (second (last (forward-batch network X nil))) Y))
     (count X)))

(defn weight-penalty
  "Weight decay's penalty, lambda times the sum of the squares of every
  weight and kernel entry."
  [network lambda]
  (* lambda (reduce + (for [layer network k weights :let [P (k layer)] :when P
                            x (if (= k :kernels) (flatten P) (apply concat P))]
                        (* x x)))))

;; ---------------------------------------------------------------------------
;; Dropout

(defn- unit-mask
  "`n` units' mask: each kept with probability `p`, its output then
  multiplied by 1/p, or else fixed at zero."
  [n p rng]
  (let [up (/ 1.0 p)] (vec (repeatedly n #(if (< (rng) p) up 0.0)))))

(defn dropout-masks
  "The masks for one step of dropout, `{:input mask :layers [mask ...]}`:
  the input's from the first layer's `:keep-input p`, and each hidden
  layer's outputs' from its `:keep p` -- every unit kept with probability
  p and multiplied by 1/p, or fixed at zero; nil where a layer has no
  `:keep`, and always for the output layer. Multiplying the kept units by
  1/p leaves each unit's expected output what it is with every unit
  active, so at test time the network runs with no dropout at all."
  [network rng]
  (let [last-l (dec (count network))
        p-in (:keep-input (first network))]
    {:input (when (and p-in (< p-in 1.0)) (unit-mask (layer-in-size (first network)) p-in rng))
     :layers (vec (map-indexed (fn [l {p :keep :as layer}]
                                 (when (and p (< p 1.0) (< l last-l)) (unit-mask (layer-out-size layer) p rng)))
                               network))}))

(defn with-dropout
  "`network` trained with dropout: every hidden unit kept with probability
  `hidden`, and with `input`, every input unit with that -- the book's
  0.5 and 0.8."
  ([network hidden] (with-dropout network hidden nil))
  ([network hidden input]
   (cond-> (into (mapv #(assoc % :keep hidden) (pop network)) [(peek network)])
     input (assoc-in [0 :keep-input] input))))

(defn- dropout? [network]
  (or (some :keep network) (:keep-input (first network))))

;; ---------------------------------------------------------------------------
;; Training

(defn minibatch-step
  "One step of stochastic gradient descent on the minibatch `X` `Y`: the
  network moved against the minibatch's mean gradient by `alpha`, with
  weight decay `lambda` if given. Given `rng`, a network with dropout
  takes the step as a thinned network, its units deactivated by fresh
  `dropout-masks` -- the same ones for the whole minibatch."
  [network X Y {:keys [alpha loss rng lambda] :or {alpha 0.1 loss :squared}}]
  (if (and rng (dropout? network))
    (let [{:keys [input layers]} (dropout-masks network rng)]
      (apply-gradients network (batch-gradients network (if input (mapv #(mapv * % input) X) X) Y loss layers) alpha lambda))
    (apply-gradients network (batch-gradients network X Y loss) alpha lambda)))

(defn sgd-step
  "One stochastic gradient step on one example: `minibatch-step` on a
  minibatch of one. `train-seq` is this over every example an epoch; a
  caller that must not stall -- a demo training between frames -- can
  take it a few examples at a time."
  ([network x y alpha loss] (sgd-step network x y alpha loss nil))
  ([network x y alpha loss rng] (minibatch-step network [x] [y] {:alpha alpha :loss loss :rng rng})))

(defn train-seq
  "Stochastic gradient descent, an epoch a step: the examples shuffled,
  and the weights moved against each minibatch's gradient by `alpha` --
  with dropout, each step through a newly thinned network. Options
  `:alpha` (0.1), `:loss` (`:squared`), `:batch` (1), `:lambda` (no
  weight decay), `:rng`."
  [network X Y {:keys [rng batch] :or {batch 1} :as opts}]
  (let [rng (or rng (random/rng 1))
        opts (assoc opts :rng rng)
        pairs (vec (map vector X Y))]
    (iterate (fn [net]
               (reduce (fn [net chunk] (minibatch-step net (mapv first chunk) (mapv second chunk) opts))
                       net (partition-all batch (random/shuffle rng pairs))))
             network)))
