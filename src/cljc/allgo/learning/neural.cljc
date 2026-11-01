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

(defn- pre-activation
  "in, before the activation function."
  [{:keys [type W b kernels shape] :as layer} a]
  (case type
    :dense (mapv (fn [row bj] (+ bj (reduce + (map * row a)))) W b)
    :conv (let [[_ w] shape [oh ow] (conv-out-shape layer) k (count (first kernels))]
            (vec (for [f (range (count kernels)) r (range oh) c (range ow)]
                   (+ (b f) (reduce + (for [i (range k) j (range k)]
                                        (* (get-in kernels [f i j]) (a (+ (* (+ r i) w) c j)))))))))))

(defn forward
  "Every layer's `[in a]`, input first: `[[nil x] [in1 a1] ...]`."
  [network x]
  (reductions (fn [[_ a] layer]
                (let [z (pre-activation layer a)] [z (mapv #(g (:g layer) %) z)]))
              [nil (vec x)] network))

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
                d (fn [f r c] (delta (+ (* f oh ow) (* r ow) c)))]
            {:dK (vec (for [f (range (count kernels))]
                        (vec (for [i (range k)]
                               (vec (for [j (range k)]
                                      (reduce + (for [r (range oh) c (range ow)] (* (d f r c) (a (+ (* (+ r i) w) c j)))))))))))
             :db (vec (for [f (range (count kernels))] (reduce + (for [r (range oh) c (range ow)] (d f r c)))))
             :back (let [acc (volatile! (vec (repeat (count a) 0.0)))]
                     (doseq [f (range (count kernels)) r (range oh) c (range ow) i (range k) j (range k)]
                       (vswap! acc update (+ (* (+ r i) w) c j) + (* (d f r c) (get-in kernels [f i j]))))
                     @acc)})))

(defn gradients
  "The gradient of the loss on one example with respect to every weight,
  a layer at a time: the output deltas, then each layer's deltas from the
  next one's, as BACK-PROP-LEARNING propagates them."
  [network x y loss]
  (let [acts (vec (forward network x))
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
            (let [[z a] (acts l) gl (:g (network (dec l)))]
              (recur (dec l) (mapv (fn [bi zi ai] (* bi (g' gl zi ai))) back z a) grads))))))))

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

(defn train-seq
  "Stochastic gradient descent, an epoch a step: the examples shuffled,
  and the weights moved against each one's gradient by `alpha`."
  [network X Y {:keys [alpha loss rng] :or {alpha 0.1 loss :squared}}]
  (let [rng (or rng (random/rng 1))
        pairs (vec (map vector X Y))]
    (iterate (fn [net]
               (reduce (fn [net [x y]] (apply-gradients net (gradients net x y loss) alpha))
                       net (random/shuffle rng pairs)))
             network)))
