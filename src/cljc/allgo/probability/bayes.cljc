(ns allgo.probability.bayes
  "Bayesian networks: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 13.

    exact          ENUMERATION-ASK (figure 13.11), and ELIMINATION-ASK
                   (figure 13.13) by factors summed out in reverse
                   topological order
    sampling       PRIOR-SAMPLE (13.16), REJECTION-SAMPLING (13.17),
                   LIKELIHOOD-WEIGHTING (13.18) and GIBBS-ASK (13.20),
                   each a seq of estimates, one sample a step

  A network is `{:variables :parents :domains :cpt}`: the variables in
  topological order, each one's parents, its values, and its conditional
  table, a map from the vector of its parents' values to a distribution
  `{value p}`. `boolean-net` builds one the way aima-python (MIT) writes
  them, with `burglary` and `sprinkler` the chapter's two networks as it
  states them. A distribution is a map from values to probabilities."
  (:require [allgo.random :as random]))

;; ---------------------------------------------------------------------------
;; Networks

(defn boolean-net
  "A network of true/false variables from `[name parents table]` rows,
  parents first: `table` is P(name = true) given the parents' values --
  a number for no parents, else a map from the tuple of parent values
  (one value alone for one parent) to it."
  [rows]
  {:variables (mapv first rows)
   :parents (into {} (map (fn [[n ps]] [n (vec ps)]) rows))
   :domains (into {} (map (fn [[n]] [n [true false]]) rows))
   :cpt (into {} (map (fn [[n ps t]]
                        [n (into {} (for [vals (if (empty? ps)
                                                 [[]]
                                                 (reduce (fn [acc _] (for [a acc v [true false]] (conj a v))) [[]] ps))
                                          :let [p (cond (number? t) t
                                                        (= 1 (count ps)) (t (first vals))
                                                        :else (t vals))]]
                                      [vals {true p false (- 1.0 p)}]))])
                      rows))})

(def burglary
  "The chapter's burglary network, as aima-python's `burglary`."
  (boolean-net [[:Burglary [] 0.001]
                [:Earthquake [] 0.002]
                [:Alarm [:Burglary :Earthquake] {[true true] 0.95 [true false] 0.94 [false true] 0.29 [false false] 0.001}]
                [:JohnCalls [:Alarm] {true 0.90 false 0.05}]
                [:MaryCalls [:Alarm] {true 0.70 false 0.01}]]))

(def sprinkler
  "The chapter's sprinkler network, as aima-python's `sprinkler`."
  (boolean-net [[:Cloudy [] 0.5]
                [:Sprinkler [:Cloudy] {true 0.10 false 0.50}]
                [:Rain [:Cloudy] {true 0.80 false 0.20}]
                [:WetGrass [:Sprinkler :Rain] {[true true] 0.99 [true false] 0.90 [false true] 0.90 [false false] 0.0}]]))

(defn p
  "P(var = value | the values `event` gives its parents)."
  [{:keys [parents cpt]} var value event]
  (get-in cpt [var (mapv event (parents var)) value] 0.0))

(defn normalize [dist]
  (let [total (reduce + (vals dist))]
    (if (pos? total) (into {} (map (fn [[k v]] [k (/ v total)]) dist)) dist)))

;; ---------------------------------------------------------------------------
;; Exact inference

(defn- enumerate-all [bn vars e]
  (if-let [y (first vars)]
    (if (contains? e y)
      (* (p bn y (e y) e) (enumerate-all bn (rest vars) e))
      (reduce + (for [v ((:domains bn) y)]
                  (* (p bn y v e) (enumerate-all bn (rest vars) (assoc e y v))))))
    1.0))

(defn enumeration-ask
  "ENUMERATION-ASK: P(`x` | evidence `e`), summing over every hidden
  variable's every value, depth first."
  [x e bn]
  (normalize (into {} (for [v ((:domains bn) x)] [v (enumerate-all bn (:variables bn) (assoc e x v))]))))

;; A factor is `{:vars [...] :table {[values...] number}}`.

(defn- make-factor [bn var e]
  (let [vars (vec (remove #(contains? e %) (cons var ((:parents bn) var))))
        assignments (reduce (fn [acc v] (for [a acc x ((:domains bn) v)] (assoc a v x))) [e] vars)]
    {:vars vars
     :table (into {} (for [a assignments] [(mapv a vars) (p bn var (a var) a)]))}))

(defn- pointwise-product [bn f g]
  (let [vars (vec (distinct (concat (:vars f) (:vars g))))
        assignments (reduce (fn [acc v] (for [a acc x ((:domains bn) v)] (assoc a v x))) [{}] vars)]
    {:vars vars
     :table (into {} (for [a assignments]
                       [(mapv a vars) (* ((:table f) (mapv a (:vars f))) ((:table g) (mapv a (:vars g))))]))}))

(defn- sum-out [bn var factors]
  (let [{with true without false} (group-by #(boolean (some #{var} (:vars %))) factors)
        product (reduce #(pointwise-product bn %1 %2) with)
        i (first (keep-indexed (fn [k v] (when (= v var) k)) (:vars product)))
        vars (vec (concat (subvec (:vars product) 0 i) (subvec (:vars product) (inc i))))]
    (conj (vec without)
          {:vars vars
           :table (reduce (fn [t [k v]] (update t (vec (concat (subvec k 0 i) (subvec k (inc i)))) (fnil + 0.0) v))
                          {} (:table product))})))

(defn elimination-ask
  "ELIMINATION-ASK: P(`x` | `e`) by variable elimination: a factor per
  variable, in reverse topological order, each hidden variable summed out
  as soon as its factor is made."
  [x e bn]
  (let [factors (reduce (fn [fs var]
                          (let [fs (conj fs (make-factor bn var e))]
                            (if (or (= var x) (contains? e var)) fs (sum-out bn var fs))))
                        [] (reverse (:variables bn)))
        product (reduce #(pointwise-product bn %1 %2) factors)]
    (normalize (into {} (map (fn [[k v]] [(first k) v]) (:table product))))))

;; ---------------------------------------------------------------------------
;; Sampling

(defn- draw [rng dist]
  (let [r (double (rng))]
    (loop [[[v pr] & more] (seq dist) acc 0.0]
      (if (or (empty? more) (< r (+ acc pr))) v (recur more (+ acc pr))))))

(defn prior-sample
  "PRIOR-SAMPLE: an event drawn from the network's joint distribution."
  [bn rng]
  (reduce (fn [x var] (assoc x var (draw rng (get-in (:cpt bn) [var (mapv x ((:parents bn) var))]))))
          {} (:variables bn)))

(defn- estimates
  "A seq of normalized estimates, `{:estimate :n}`, from a seq of
  `[value weight]` samples (weight 0 rejects)."
  [samples]
  (rest (reductions (fn [{:keys [counts n]} [v w]]
                      (let [counts (if (pos? w) (update counts v (fnil + 0.0) w) counts)]
                        {:counts counts :n (inc (long n)) :estimate (normalize counts)}))
                    {:counts {} :n 0}
                    samples)))

(defn rejection-sampling-seq
  "REJECTION-SAMPLING's estimate of P(`x` | `e`) after each prior sample:
  samples that contradict the evidence are thrown away."
  [x e bn rng]
  (estimates (repeatedly (fn [] (let [s (prior-sample bn rng)]
                                  [(s x) (if (every? (fn [[k v]] (= v (s k))) e) 1.0 0.0)])))))

(defn weighted-sample
  "WEIGHTED-SAMPLE: an event with the evidence fixed and the rest drawn,
  and its weight, the likelihood of the evidence given its parents."
  [bn e rng]
  (reduce (fn [[x w] var]
            (if (contains? e var)
              [(assoc x var (e var)) (* w (p bn var (e var) x))]
              [(assoc x var (draw rng (get-in (:cpt bn) [var (mapv x ((:parents bn) var))]))) w]))
          [{} 1.0] (:variables bn)))

(defn likelihood-weighting-seq
  "LIKELIHOOD-WEIGHTING's estimate of P(`x` | `e`) after each weighted
  sample."
  [x e bn rng]
  (estimates (repeatedly (fn [] (let [[s w] (weighted-sample bn e rng)] [(s x) w])))))

(defn- markov-blanket-sample
  "A value for `var` drawn from P(var | its Markov blanket), with the rest
  of `x` fixed: its parents through its own table, its children through
  theirs."
  [bn var x rng]
  (let [children (filter #(some #{var} ((:parents bn) %)) (:variables bn))
        dist (normalize (into {} (for [v ((:domains bn) var)
                                       :let [x' (assoc x var v)]]
                                   [v (reduce * (p bn var v x') (map #(p bn % (x' %) x') children))])))]
    (draw rng dist)))

(defn gibbs-seq
  "GIBBS-ASK's estimate of P(`x` | `e`) after each step: one nonevidence
  variable resampled from its Markov blanket, chosen in turn (the 4th
  edition lets the order be any distribution; round robin is one)."
  [x e bn rng]
  (let [hidden (vec (remove #(contains? e %) (:variables bn)))
        start (merge (into {} (map (fn [v] [v (random/pick rng ((:domains bn) v))]) hidden)) e)]
    (estimates (map (fn [state] [(state x) 1.0])
                    (rest (reductions (fn [state k]
                                        (let [var (hidden (mod k (count hidden)))]
                                          (assoc state var (markov-blanket-sample bn var state rng))))
                                      start (range)))))))
