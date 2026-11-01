(ns allgo.learning.models
  "Parametric and nonparametric models: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapters 19 and 20,
  whose text teaches these without pseudocode.

    linear regression     the exact least-squares weights, by
                          `allgo.numerics.regression/multiple`: one
                          step, no learning rate. The chapter also
                          teaches gradient descent on it; that is
                          `gradient-descent-seq`, kept for the models
                          that have no closed form
    perceptron            the threshold classifier's learning rule
    logistic regression   gradient descent on the log loss
    k-nearest neighbors   the majority of the k nearest, by any distance
    support vector        the dual problem by sequential minimal
    machines              optimization (Platt's SMO, simplified), with a
                          kernel -- linear, polynomial, or Gaussian
    naive Bayes           counts with Laplace smoothing, and its Gaussian
                          form for real-valued attributes

  Inputs are vectors of numbers unless a model says otherwise; labels for
  the classifiers here are 0 and 1, or -1 and 1 for the SVM."
  (:require [allgo.numerics.regression :as regression]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn- dot [a b] (reduce + (map * a b)))

;; ---------------------------------------------------------------------------
;; Linear models

(defn linear-regression
  "The least-squares weights of y = w0 + w1 x1 + ...: `[w0 w1 ...]`."
  [X ys]
  (vec (:a (regression/multiple X ys))))

(defn- with-bias [x] (into [1.0] x))

(defn gradient-descent-seq
  "The weights after each epoch of batch gradient descent: each step w <-
  w - alpha times the mean gradient of `loss-gradient` (fn [w x y]) over
  the examples, x with a 1 in front for the bias."
  [X ys alpha loss-gradient]
  (let [Xb (mapv with-bias X)
        n (count ys)]
    (iterate (fn [w]
               (let [g (reduce (fn [acc [x y]] (mapv + acc (loss-gradient w x y))) (vec (repeat (count w) 0.0)) (map vector Xb ys))]
                 (mapv (fn [wi gi] (- wi (* alpha (/ gi n)))) w g)))
             (vec (repeat (inc (count (first X))) 0.0)))))

(defn predict-linear [w x] (dot w (with-bias x)))

(defn- sigmoid [z] (/ 1.0 (+ 1.0 (math/exp (- z)))))

(defn logistic-regression-seq
  "Gradient descent on the log loss of the logistic function."
  [X ys alpha]
  (gradient-descent-seq X ys alpha (fn [w x y] (let [err (- (sigmoid (dot w x)) y)] (mapv #(* err %) x)))))

(defn predict-logistic
  "P(y = 1 | x) under the logistic model `w`."
  [w x]
  (sigmoid (dot w (with-bias x))))

(defn perceptron-seq
  "The perceptron learning rule, an epoch a step: for each example in
  turn, w <- w + alpha (y - h(x)) x, h the hard threshold. It converges
  when the data are linearly separable, and not otherwise."
  [X ys alpha]
  (let [Xb (mapv with-bias X)]
    (iterate (fn [w]
               (reduce (fn [w [x y]]
                         (let [h (if (>= (dot w x) 0.0) 1 0)]
                           (mapv (fn [wi xi] (+ wi (* alpha (- y h) xi))) w x)))
                       w (map vector Xb ys)))
             (vec (repeat (inc (count (first X))) 0.0)))))

(defn predict-threshold [w x] (if (>= (dot w (with-bias x)) 0.0) 1 0))

;; ---------------------------------------------------------------------------
;; Nearest neighbors

(defn euclidean [a b] (math/sqrt (reduce + (map #(let [d (- %1 %2)] (* d d)) a b))))

(defn knn
  "The k-nearest-neighbors classifier on `X` `ys`: (fn [x]) -> the label
  most of the `k` nearest have, ties to the nearer."
  ([X ys k] (knn X ys k euclidean))
  ([X ys k distance]
   (fn [x]
     (let [near (take k (sort-by first (map (fn [xi yi] [(distance x xi) yi]) X ys)))
           counts (frequencies (map second near))
           top (apply max (vals counts))]
       (second (first (filter #(= top (counts (second %))) near)))))))

;; ---------------------------------------------------------------------------
;; Support vector machines

(defn linear-kernel [a b] (dot a b))
(defn polynomial-kernel [degree] (fn [a b] (math/pow (+ 1.0 (dot a b)) degree)))
(defn gaussian-kernel [gamma] (fn [a b] (math/exp (- (* gamma (reduce + (map #(let [d (- %1 %2)] (* d d)) a b)))))))

(defn- svm-output [{:keys [X ys alphas b kernel]} x]
  (+ b (reduce + (map (fn [xi yi ai] (if (zero? ai) 0.0 (* ai yi (kernel xi x)))) X ys alphas))))

(defn svm-seq
  "Simplified SMO (Platt 1998, as it is usually taught): each step a pass
  over the examples, each that breaks the KKT conditions by more than
  `tol` optimized jointly with another chosen at random, the pair's
  multipliers moved analytically within [0, C]. Labels are -1 and 1.
  Each state is the model `{:alphas :b ...}`; `svm-predict` classifies
  with it."
  [X ys {:keys [C kernel tol rng] :or {C 1.0 tol 1e-3}}]
  (let [kernel (or kernel linear-kernel)
        rng (or rng (random/rng 1))
        n (count ys)
        X (vec X) ys (vec ys)
        K (vec (for [i (range n)] (vec (for [j (range n)] (kernel (X i) (X j))))))
        f (fn [alphas b i] (+ b (reduce + (map (fn [j] (* (alphas j) (ys j) ((K j) i))) (range n)))))]
    (iterate
     (fn [{:keys [alphas b] :as model}]
       (let [[alphas b changed]
             (reduce
              (fn [[alphas b changed] i]
                (let [Ei (- (f alphas b i) (ys i))
                      yi (ys i) ai (alphas i)]
                  (if-not (or (and (< (* yi Ei) (- tol)) (< ai C)) (and (> (* yi Ei) tol) (> ai 0.0)))
                    [alphas b changed]
                    (let [j (let [r (random/below rng (dec n))] (if (>= r i) (inc r) r))
                          Ej (- (f alphas b j) (ys j))
                          yj (ys j) aj (alphas j)
                          [L H] (if (not= yi yj)
                                  [(max 0.0 (- aj ai)) (min C (+ C (- aj ai)))]
                                  [(max 0.0 (- (+ ai aj) C)) (min C (+ ai aj))])
                          eta (- (* 2.0 ((K i) j)) ((K i) i) ((K j) j))]
                      (if (or (>= L H) (>= eta 0.0))
                        [alphas b changed]
                        (let [aj' (min H (max L (- aj (/ (* yj (- Ei Ej)) eta))))]
                          (if (< (abs (- aj' aj)) 1e-7)
                            [alphas b changed]
                            (let [ai' (+ ai (* yi yj (- aj aj')))
                                  b1 (- b Ei (* yi (- ai' ai) ((K i) i)) (* yj (- aj' aj) ((K i) j)))
                                  b2 (- b Ej (* yi (- ai' ai) ((K i) j)) (* yj (- aj' aj) ((K j) j)))
                                  b' (cond (< 0.0 ai' C) b1 (< 0.0 aj' C) b2 :else (/ (+ b1 b2) 2.0))]
                              [(assoc alphas i ai' j aj') b' (inc changed)]))))))))
              [alphas b 0] (range n))]
         (assoc model :alphas alphas :b b :changed changed :passes (inc (long (:passes model))))))
     {:X X :ys ys :kernel kernel :alphas (vec (repeat n 0.0)) :b 0.0 :changed nil :passes 0})))

(defn svm-train
  "An SVM trained until `max-quiet` (default 5) passes in a row change
  nothing, or `max-passes` (default 200)."
  [X ys opts]
  (let [{:keys [max-quiet max-passes] :or {max-quiet 5 max-passes 200}} opts]
    (loop [[m & more] (rest (svm-seq X ys opts)) quiet 0]
      (let [quiet (if (zero? (long (:changed m))) (inc quiet) 0)]
        (if (or (>= quiet max-quiet) (>= (long (:passes m)) max-passes)) m (recur more quiet))))))

(defn svm-predict [model x] (if (>= (svm-output model x) 0.0) 1 -1))

(defn support-vectors [{:keys [X alphas]}] (keep-indexed (fn [i a] (when (> a 1e-8) (X i))) alphas))

;; ---------------------------------------------------------------------------
;; Naive Bayes

(defn naive-bayes
  "The naive Bayes classifier on examples, maps of discrete attributes:
  P(class) and each P(attribute = value | class) counted, with Laplace
  smoothing of `k` (default 1). (fn [example]) -> `{:class :posterior}`."
  ([examples attributes target] (naive-bayes examples attributes target 1.0))
  ([examples attributes target k]
   (let [classes (distinct (map target examples))
         values (into {} (for [a attributes] [a (set (map a examples))]))
         by-class (group-by target examples)
         n (count examples)]
     (fn [example]
       (let [log-score (into {} (for [c classes
                                      :let [exs (by-class c) nc (count exs)]]
                                  [c (+ (math/log (/ (double nc) n))
                                        (reduce + (for [a attributes :when (contains? example a)]
                                                    (math/log (/ (+ k (count (filter #(= (a example) (a %)) exs)))
                                                                 (+ nc (* k (count (values a)))))))))]))
             top (apply max (vals log-score))
             unnorm (into {} (map (fn [[c s]] [c (math/exp (- s top))]) log-score))
             total (reduce + (vals unnorm))]
         {:class (key (apply max-key val log-score))
          :posterior (into {} (map (fn [[c u]] [c (/ u total)]) unnorm))})))))

(defn gaussian-naive-bayes
  "Naive Bayes for real-valued inputs, each attribute Gaussian within a
  class: (fn [x]) -> the likeliest label."
  [X ys]
  (let [by (group-by second (map vector X ys))
        n (count ys)
        stats (into {} (for [[c pairs] by
                             :let [xs (map first pairs) m (count xs) d (count (first xs))
                                   mean (mapv (fn [j] (/ (reduce + (map #(nth % j) xs)) m)) (range d))
                                   var (mapv (fn [j] (+ 1e-9 (/ (reduce + (map #(let [v (- (nth % j) (mean j))] (* v v)) xs)) m))) (range d))]]
                         [c {:prior (/ (double m) n) :mean mean :var var}]))]
    (fn [x]
      (key (apply max-key
                  (fn [[_ {:keys [prior mean var]}]]
                    (+ (math/log prior)
                       (reduce + (map (fn [xi mu v] (- (* -0.5 (math/log (* 2 math/PI v))) (/ (* (- xi mu) (- xi mu)) (* 2 v)))) x mean var))))
                  stats)))))
