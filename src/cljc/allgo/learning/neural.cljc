(ns allgo.learning.neural
  "Neural networks: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 21 -- BACK-PROP-LEARNING (figure
  21's back-propagation; the 4th edition gives it in the text), with
  dense layers and a convolutional one.

    dense         in_j = sum over i of w_ij a_i + b_j, a_j = g(in_j)
    convolution   2-D kernels slid over a grid, stride 1, no padding: the
                  layers of a convolutional network, flattened into the
                  dense ones after
    training      stochastic gradient descent, an epoch a step, the
                  examples shuffled each epoch; the deltas propagated back
                  from the output exactly as the figure does
    dropout       section 21.5.4: each step a thinned network, every
                  hidden unit (`:keep p`) and input unit (`:keep-input`)
                  kept with probability p and its output multiplied by
                  1/p, or else zero; at test time no dropout

  A network is a vector of layers:

    {:type :dense :W [[...]] :b [...] :g :sigmoid|:tanh|:relu|:identity}
    {:type :conv :kernels [[[...]]] :b [...] :g ... :shape [h w]}

  a conv layer taking a flat input it reads as an `h` by `w` grid and
  giving its feature maps flattened, one after another. The loss is the
  squared error, or with a sigmoid output, `:cross-entropy`."
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

(defn- conv-out-shape [{:keys [shape kernels]}]
  (let [[h w] shape k (count (first kernels))] [(inc (- h k)) (inc (- w k))]))

(defn layer-out-size [{:keys [type W kernels] :as layer}]
  (if (= type :dense) (count W) (let [[oh ow] (conv-out-shape layer)] (* (count kernels) oh ow))))

(defn- layer-in-size [{:keys [type W shape]}]
  (if (= type :dense) (count (first W)) (apply * shape)))

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

(defn- pre-activation
  "in, before the activation function."
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
                     (for [f (range (count kernels)) r (range oh) c (range ow)] [f r c]))))))

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

;; ---------------------------------------------------------------------------
;; Back-propagation

(defn- back
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
                            (for [f (range (count kernels)) r (range oh) c (range ow)] [f r c])))})))

(defn gradients
  "The gradient of the loss on one example with respect to every weight,
  a layer at a time: the output deltas, then each layer's deltas from the
  next one's, as BACK-PROP-LEARNING propagates them. With dropout's
  `masks` (`dropout-masks`), the gradient of the thinned network's: a
  dropped unit passes back nothing, a kept one its delta times its mask."
  ([network x y loss] (gradients network x y loss nil))
  ([network x y loss masks]
   (let [acts (vec (forward network x masks))
         n (count network)
         [z-out a-out] (acts n)
         out-layer (network (dec n))
        ;; With a sigmoid output and cross-entropy the g' cancels.
         delta (if (and (= loss :cross-entropy) (= :sigmoid (:g out-layer)))
                 (mapv - a-out y)
                 (mapv (fn [z a t] (* 2.0 (- a t) (g' (:g out-layer) z a))) z-out a-out y))]
     (loop [l (dec n) delta delta grads ()]
       (if (neg? l)
         (vec grads)
         (let [{:keys [back] :as gr} (back (network l) delta (second (acts l)))
               grads (cons (dissoc gr :back) grads)]
           (if (zero? l)
             (vec grads)
             (let [[z a] (acts l) gl (:g (network (dec l))) m (get masks (dec l))]
               (recur (dec l)
                      (if m
                       ;; a = m g(in): the unmasked g(in) for g'
                        (mapv (fn [bi zi mi] (if (zero? mi) 0.0 (* bi mi (g' gl zi (g gl zi))))) back z m)
                        (mapv (fn [bi zi ai] (* bi (g' gl zi ai))) back z a))
                      grads)))))))))

(defn- apply-gradients [network grads alpha]
  (mapv (fn [layer {:keys [dW dK db]}]
          (cond-> (update layer :b (fn [b] (mapv #(- %1 (* alpha %2)) b db)))
            dW (update :W (fn [W] (mapv (fn [row drow] (mapv #(- %1 (* alpha %2)) row drow)) W dW)))
            dK (update :kernels (fn [K] (mapv (fn [kk dk] (mapv (fn [row drow] (mapv #(- %1 (* alpha %2)) row drow)) kk dk)) K dK)))))
        network grads))

(defn loss
  "The mean loss of `network` on examples `X` `Y`."
  [network X Y kind]
  (/ (reduce + (map (fn [x y]
                      (let [a (predict network x)]
                        (if (= kind :cross-entropy)
                          (- (reduce + (map (fn [ai yi] (+ (* yi (math/log (max 1e-12 ai))) (* (- 1 yi) (math/log (max 1e-12 (- 1 ai)))))) a y)))
                          (reduce + (map #(let [d (- %1 %2)] (* d d)) a y)))))
                    X Y))
     (count X)))

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

(defn sgd-step
  "One stochastic gradient step on one example: the network moved
  against its gradient by `alpha`. `train-seq` is this over every example
  an epoch; a caller that must not stall -- a demo training between
  frames -- can take it a few examples at a time. Given `rng`, a network
  with dropout takes the step as a thinned network, its units deactivated
  by fresh `dropout-masks`."
  ([network x y alpha loss] (sgd-step network x y alpha loss nil))
  ([network x y alpha loss rng]
   (if (and rng (dropout? network))
     (let [{:keys [input layers]} (dropout-masks network rng)]
       (apply-gradients network (gradients network (if input (mapv * x input) x) y loss layers) alpha))
     (apply-gradients network (gradients network x y loss) alpha))))

(defn train-seq
  "Stochastic gradient descent, an epoch a step: the examples shuffled,
  and the weights moved against each one's gradient by `alpha` -- with
  dropout, each step through a newly thinned network."
  [network X Y {:keys [alpha loss rng] :or {alpha 0.1 loss :squared}}]
  (let [rng (or rng (random/rng 1))
        pairs (vec (map vector X Y))]
    (iterate (fn [net]
               (reduce (fn [net [x y]] (sgd-step net x y alpha loss rng))
                       net (random/shuffle rng pairs)))
             network)))
