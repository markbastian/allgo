(ns allgo.learning.autodiff
  "Reverse-mode automatic differentiation over matrices: what training a
  recurrent network, an attention model or a transformer needs, after
  Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
  edition, chapters 21 and 24 -- back-propagation, done once for every
  architecture rather than derived by hand for each.

  A computation is built from nodes, each a matrix value -- a vector of
  row vectors, a scalar being 1 x 1 -- and how to send a gradient back
  to the nodes it came from. `value-and-grad` runs a function of named
  parameters forward, then walks the graph from its scalar result
  backwards, each node's gradient the sum of what its consumers send it:

    linear      matmul, add, add-row (a row added to every row), sub,
                mul and mul-row (elementwise), scale, transpose
    nonlinear   tanh, sigmoid, relu, exp
    structure   rows (a gather: embedding lookup), row, stack-rows,
                concat-cols, slice-cols
    reductions  sum, mean, softmax-rows, layer-norm-rows
    losses      cross-entropy: rows of logits against target indices,
                softmax and log fused, as is usual, for stability

  and training: `adam-step`, and `clip` by the gradients' global norm.

  Row vectors throughout: an input is a row, a layer x W + b."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def ^:private ids (atom 0))

(defn- node [value parents]
  {:id (swap! ids inc) :value value :parents parents})

(defn const
  "A node holding `m`, a matrix, which nothing differentiates."
  [m]
  (node m []))

(defn value [n] (:value n))

(defn- zeros [r c] (vec (repeat r (vec (repeat c 0.0)))))
(defn- shape [m] [(count m) (count (first m))])
(defn- emap [f m] (mapv (fn [row] (mapv f row)) m))
(defn- emap2 [f a b] (mapv (fn [r s] (mapv f r s)) a b))
(defn- col-sums [m] (apply mapv + m))

;; ---------------------------------------------------------------------------
;; Operations: each a value and, per parent, what a gradient G of the
;; result sends back to it

(defn- row-dot
  "The dot product of two rows of length `n`."
  [a b n]
  (loop [i 0 acc 0.0]
    (if (= i n) acc (recur (inc i) (+ acc (* (double (nth a i)) (double (nth b i))))))))

(defn- mat-mul
  "A B, in loops: the product every layer makes, and most of the time a
  model takes to train."
  [A B]
  (let [Bt (lin/transpose B) n (count (first A))]
    (mapv (fn [a] (mapv #(row-dot a % n) Bt)) A)))

(defn matmul [a b]
  (let [A (value a) B (value b)]
    (node (mat-mul A B)
          [[a #(mat-mul % (lin/transpose B))]
           [b #(mat-mul (lin/transpose A) %)]])))

(defn add [a b] (node (emap2 + (value a) (value b)) [[a identity] [b identity]]))
(defn sub [a b] (node (emap2 - (value a) (value b)) [[a identity] [b #(emap - %)]]))

(defn add-row
  "`r`, one row, added to every row of `a`."
  [a r]
  (let [R (first (value r))]
    (node (mapv #(mapv + % R) (value a)) [[a identity] [r #(vector (col-sums %))]])))

(defn mul
  "Elementwise product."
  [a b]
  (let [A (value a) B (value b)]
    (node (emap2 * A B) [[a #(emap2 * % B)] [b #(emap2 * % A)]])))

(defn mul-row
  "Every row of `a` multiplied elementwise by the row `r`."
  [a r]
  (let [A (value a) R (first (value r))]
    (node (mapv #(mapv * % R) A)
          [[a (fn [G] (mapv #(mapv * % R) G))]
           [r (fn [G] [(col-sums (emap2 * G A))])]])))

(defn scale [a k] (node (emap #(* k %) (value a)) [[a (fn [G] (emap #(* k %) G))]]))

(defn transpose [a] (node (lin/transpose (value a)) [[a lin/transpose]]))

(defn- pointwise [a f df]
  ;; df of the input x and the output y
  (let [A (value a) Y (emap f A)]
    (node Y [[a (fn [G] (mapv (fn [g x y] (mapv (fn [gi xi yi] (* gi (df xi yi))) g x y)) G A Y))]])))

(defn tanh [a] (pointwise a math/tanh (fn [_ y] (- 1.0 (* y y)))))
(defn sigmoid [a] (pointwise a #(/ 1.0 (+ 1.0 (math/exp (- %)))) (fn [_ y] (* y (- 1.0 y)))))
(defn relu [a] (pointwise a #(max 0.0 %) (fn [x _] (if (pos? x) 1.0 0.0))))
(defn exp [a] (pointwise a math/exp (fn [_ y] y)))

(defn rows
  "The rows of `a` at `indices`, in order -- an embedding lookup. The
  gradient is added back into those rows."
  [a indices]
  (let [A (value a) [r c] (shape A)]
    (node (mapv A indices)
          [[a (fn [G] (reduce (fn [acc [i g]] (update acc i #(mapv + % g))) (zeros r c) (map vector indices G)))]])))

(defn row [a i] (rows a [i]))

(defn stack-rows
  "The nodes' rows one under another."
  [ns]
  (let [counts (mapv #(count (value %)) ns)
        starts (vec (reductions + 0 counts))]
    (node (vec (mapcat value ns))
          (mapv (fn [n s c] [n (fn [G] (subvec G s (+ s c)))]) ns starts counts))))

(defn concat-cols
  "The nodes side by side."
  [ns]
  (let [widths (mapv #(count (first (value %))) ns)
        starts (vec (reductions + 0 widths))]
    (node (apply mapv (fn [& rs] (vec (apply concat rs))) (map value ns))
          (mapv (fn [n s w] [n (fn [G] (mapv #(subvec % s (+ s w)) G))]) ns starts widths))))

(defn slice-cols
  "Columns `from` (inclusive) to `to` (exclusive) of `a`."
  [a from to]
  (let [A (value a) c (count (first A))]
    (node (mapv #(subvec % from to) A)
          [[a (fn [G] (mapv (fn [g] (vec (concat (repeat from 0.0) g (repeat (- c to) 0.0)))) G))]])))

(defn sum [a] (let [A (value a) [r c] (shape A)] (node [[(reduce + (map #(reduce + %) A))]] [[a (fn [G] (vec (repeat r (vec (repeat c (ffirst G))))))]])))
(defn mean [a] (let [[r c] (shape (value a))] (scale (sum a) (/ 1.0 (* r c)))))

(defn- softmax-row [xs]
  (let [m (reduce max xs) es (mapv #(math/exp (- % m)) xs) s (reduce + es)]
    (mapv #(/ % s) es)))

(defn softmax-rows
  "Each row turned into a probability distribution."
  [a]
  (let [Y (mapv softmax-row (value a))]
    (node Y [[a (fn [G] (mapv (fn [g y] (let [d (reduce + (map * g y))] (mapv (fn [gi yi] (* yi (- gi d))) g y))) G Y))]])))

(defn layer-norm-rows
  "Each row shifted to mean 0 and scaled to variance 1 (plus `eps`,
  default 1e-5) -- a gain and bias, if wanted, applied after by `mul-row`
  and `add-row`."
  ([a] (layer-norm-rows a 1e-5))
  ([a eps]
   (let [A (value a) n (count (first A))
         stats (mapv (fn [x] (let [mu (/ (reduce + x) n)
                                   var (/ (reduce + (map #(let [d (- % mu)] (* d d)) x)) n)]
                               [mu (math/sqrt (+ var eps))]))
                     A)
         Y (mapv (fn [x [mu sd]] (mapv #(/ (- % mu) sd) x)) A stats)]
     (node Y [[a (fn [G] (mapv (fn [g y [_ sd]]
                                 (let [mg (/ (reduce + g) n) mgy (/ (reduce + (map * g y)) n)]
                                   (mapv (fn [gi yi] (/ (- gi mg (* yi mgy)) sd)) g y)))
                               G Y stats))]]))))

(defn cross-entropy
  "The mean over rows of -log softmax(logits)[target], `targets` an index
  per row. The gradient is softmax minus the one-hot target, over the
  number of rows."
  [logits targets]
  (let [L (value logits) n (count L)
        P (mapv softmax-row L)
        loss (/ (reduce + (map (fn [p t] (- (math/log (max 1e-300 (p t))))) P targets)) n)]
    (node [[loss]]
          [[logits (fn [G] (let [g (ffirst G)]
                             (mapv (fn [p t] (mapv (fn [j pj] (* g (/ (- pj (if (= j t) 1.0 0.0)) n))) (range) p)) P targets)))]])))

;; ---------------------------------------------------------------------------
;; Differentiation

(defn- topological
  "The nodes `out` depends on, each after everything it depends on."
  [out]
  (loop [stack [[out false]] seen #{} order []]
    (if-let [[n done?] (peek stack)]
      (let [stack (pop stack)]
        (cond
          done? (recur stack seen (conj order n))
          (seen (:id n)) (recur stack seen order)
          :else (recur (into (conj stack [n true]) (map (fn [[p]] [p false]) (:parents n)))
                       (conj seen (:id n)) order)))
      order)))

(defn backward
  "The gradient of the scalar node `out` with respect to every node it
  depends on: `{id gradient}`."
  [out]
  (reduce (fn [grads n]
            (if-let [G (grads (:id n))]
              (reduce (fn [grads [p back]]
                        (update grads (:id p) (fn [acc] (let [g (back G)] (if acc (emap2 + acc g) g)))))
                      grads (:parents n))
              grads))
          {(:id out) [[1.0]]}
          (reverse (topological out))))

(defn value-and-grad
  "`f`, a function of a map of parameter nodes returning a scalar node,
  at the parameters `params` -- `{name matrix}` -- and its gradient with
  respect to each: `{:value :grads}`."
  [f params]
  (let [leaves (into {} (map (fn [[k m]] [k (node m [])]) params))
        out (f leaves)
        g (backward out)]
    {:value (ffirst (value out))
     :grads (into {} (map (fn [[k n]] [k (or (g (:id n)) (let [[r c] (shape (value n))] (zeros r c)))]) leaves))}))

;; ---------------------------------------------------------------------------
;; Training

(defn clip
  "The gradients scaled down, all together, so that their global norm is
  at most `max-norm` -- the usual guard against an RNN's exploding
  gradients."
  [grads max-norm]
  (let [norm (math/sqrt (reduce + (for [[_ G] grads row G x row] (* x x))))]
    (if (<= norm max-norm)
      grads
      (into {} (map (fn [[k G]] [k (emap #(* % (/ max-norm norm)) G)]) grads)))))

(defn adam-step
  "One step of Adam: each parameter moved against its gradient by `lr`,
  scaled by running averages of the gradient and its square. `state` is
  `{:m :v :t}`, nil to start. `[params' state']`."
  [params grads state & {:keys [lr b1 b2 eps] :or {lr 0.01 b1 0.9 b2 0.999 eps 1e-8}}]
  (let [{:keys [m v t] :or {t 0}} state
        t (inc t)
        m (into {} (map (fn [[k G]] [k (emap2 #(+ (* b1 %1) (* (- 1 b1) %2)) (get m k (emap (constantly 0.0) G)) G)]) grads))
        v (into {} (map (fn [[k G]] [k (emap2 #(+ (* b2 %1) (* (- 1 b2) %2 %2)) (get v k (emap (constantly 0.0) G)) G)]) grads))
        c1 (- 1 (math/pow b1 t)) c2 (- 1 (math/pow b2 t))
        params (into {} (map (fn [[k P]]
                               [k (if-let [mk (m k)]
                                    (mapv (fn [p mr vr] (mapv (fn [pi mi vi] (- pi (/ (* lr (/ mi c1)) (+ (math/sqrt (/ vi c2)) eps)))) p mr vr)) P mk (v k))
                                    P)])
                             params))]
    [params {:m m :v v :t t}]))

(defn init
  "A `rows` by `cols` matrix of small random weights, scaled by 1 /
  sqrt(rows) as for a layer of that many inputs."
  [rows cols rng]
  (let [s (/ 1.0 (math/sqrt rows))]
    (vec (repeatedly rows (fn [] (vec (repeatedly cols #(* s (- (* 2.0 (rng)) 1.0)))))))))

(defn softmax
  "A row of numbers as a probability distribution."
  [xs]
  (softmax-row xs))
