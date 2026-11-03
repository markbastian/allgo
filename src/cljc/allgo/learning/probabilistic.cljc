(ns allgo.learning.probabilistic
  "Learning probabilistic models: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 20.

    Bayesian learning   the posterior over hypotheses, P(h | d) = alpha P(d
                        | h) P(h), a datum a step; predictions averaged
                        over every hypothesis by its posterior; the MAP and
                        maximum-likelihood hypotheses beside it
    parameters by ML    with complete data, a Bayes net's tables are the
                        observed frequencies of each value given each
                        setting of its parents; a Gaussian's mean and
                        standard deviation are the sample's; a linear
                        Gaussian's line is least squares
    parameters by       a beta prior over a Boolean parameter -- the
    Bayes               conjugate prior: each observation adds one to a
                        virtual count -- and the Dirichlet for many values
    structure           hill climbing over networks, adding, deleting or
                        reversing a link a step, never making a cycle,
                        scored by the log likelihood less a penalty for
                        the parameters (BIC's half log N each): without
                        the penalty, more links always fit better
    nonparametric       density estimation by k nearest neighbors and by
                        Gaussian kernels
    hidden variables    EM for a mixture of naive Bayes models, the class
                        hidden: the expected counts of each class from
                        the posterior over it, an iteration a step

  Networks are `allgo.probability.bayes`'s: `{:variables :parents :domains
  :cpt}`. Data are events, maps from variables to values."
  (:require [allgo.numerics.special :as special]
            [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Bayesian learning

(defn- normalize [ps] (let [t (reduce + ps)] (mapv #(/ % t) ps)))

(defn posterior-seq
  "Bayesian learning over `hypotheses`, each a function giving P(datum |
  h), from `prior`, a vector of their probabilities: the posterior
  after each of `data` in turn, the prior first."
  [hypotheses prior data]
  (reductions (fn [post d] (normalize (mapv #(* %1 (%2 d)) post hypotheses))) (vec prior) data))

(defn predict
  "P(x | d), each hypothesis's prediction weighted by its posterior."
  [hypotheses posterior x]
  (reduce + (map #(* %1 (%2 x)) posterior hypotheses)))

(defn map-hypothesis
  "The index of the most probable hypothesis under `posterior`."
  [posterior]
  (apply max-key posterior (range (count posterior))))

(defn ml-hypothesis
  "The index of the hypothesis making `data` likeliest: MAP with a
  uniform prior."
  [hypotheses data]
  (map-hypothesis (last (posterior-seq hypotheses (vec (repeat (count hypotheses) 1.0)) data))))

;; ---------------------------------------------------------------------------
;; Parameters by maximum likelihood

(defn- settings
  "Every tuple of values of `vars`."
  [domains vars]
  (reduce (fn [acc v] (for [a acc x (domains v)] (conj a x))) [[]] vars))

(defn learn-parameters
  "A network's tables from complete `data`: each variable's distribution
  given each setting of its parents the observed frequencies, `pseudo`
  (default 0) added to every count first -- 1 is Laplace's, which keeps
  a value never seen from probability zero. A setting of the parents
  never seen gets the uniform distribution."
  [{:keys [variables parents domains] :as structure} data & {:keys [pseudo] :or {pseudo 0.0}}]
  (let [counts (frequencies (for [e data v variables] [v (mapv e (parents v)) (e v)]))]
    (assoc structure :cpt
           (into {} (for [v variables]
                      [v (into {} (for [u (settings domains (parents v))
                                        :let [cs (mapv #(+ (double pseudo) (get counts [v u %] 0)) (domains v))
                                              total (reduce + cs)]]
                                    [u (zipmap (domains v) (if (pos? total) (mapv #(/ % total) cs)
                                                               (repeat (/ 1.0 (count (domains v))))))]))])))))

(defn log-likelihood
  "The natural log of the probability of `data` under the network."
  [{:keys [variables parents cpt]} data]
  (reduce + (for [e data v variables] (math/log (get-in cpt [v (mapv e (parents v)) (e v)])))))

(defn gaussian-ml
  "The maximum-likelihood Gaussian for `xs`: the sample's mean, and the
  root of its mean squared deviation -- divided by N, not N - 1.
  `{:mean :sd}`."
  [xs]
  (let [n (double (count xs)) mu (/ (reduce + xs) n)]
    {:mean mu :sd (math/sqrt (/ (reduce + (map #(let [d (- % mu)] (* d d)) xs)) n))}))

(defn linear-gaussian-ml
  "The maximum-likelihood linear Gaussian y = theta1 x + theta2 + N(0,
  sigma^2): the line least squares gives, and the root mean square of
  its residuals. `{:theta1 :theta2 :sigma}`."
  [xs ys]
  (let [n (double (count xs))
        mx (/ (reduce + xs) n) my (/ (reduce + ys) n)
        sxy (reduce + (map #(* (- %1 mx) (- %2 my)) xs ys))
        sxx (reduce + (map #(let [d (- % mx)] (* d d)) xs))
        t1 (/ sxy sxx) t2 (- my (* t1 mx))]
    {:theta1 t1 :theta2 t2
     :sigma (math/sqrt (/ (reduce + (map #(let [r (- %2 (+ (* t1 %1) t2))] (* r r)) xs ys)) n))}))

;; ---------------------------------------------------------------------------
;; Parameters by Bayesian learning

(defn beta-density
  "beta[a, b](theta) = alpha theta^(a-1) (1 - theta)^(b-1), alpha the
  1 / B(a, b) that makes it integrate to 1."
  [a b theta]
  (math/exp (+ (- (special/log-gamma (+ a b)) (special/log-gamma a) (special/log-gamma b))
               (* (dec (double a)) (math/log theta)) (* (dec (double b)) (math/log (- 1.0 theta))))))

(defn beta-mean [[a b]] (/ (double a) (+ a b)))

(defn beta-seq
  "The beta posterior over a Boolean parameter after each of
  `observations` (true or false), from the prior `[a b]`: a true adds one
  to a, a false one to b -- the beta family is closed under update, so
  its hyperparameters are virtual counts."
  [[a b] observations]
  (reductions (fn [[a b] o] (if o [(inc a) b] [a (inc b)])) [a b] observations))

(defn dirichlet-update
  "The Dirichlet posterior, `{value alpha}`, after `observations`: each
  adds one to its value's alpha."
  [alphas observations]
  (reduce (fn [m o] (update m o (fnil inc 0))) alphas observations))

(defn dirichlet-predictive
  "P(next = each value) under the Dirichlet `alphas`: each alpha over
  their sum."
  [alphas]
  (let [t (double (reduce + (vals alphas)))] (into {} (map (fn [[k a]] [k (/ a t)]) alphas))))

;; ---------------------------------------------------------------------------
;; Structure

(defn parameter-count
  "The free parameters of a network with `parents`: for each variable,
  (values - 1) for each setting of its parents."
  [domains parents]
  (reduce + (for [[v ps] parents] (* (dec (count (domains v))) (reduce * (map #(count (domains %)) ps))))))

(defn structure-score
  "The score of the structure `parents` on `data`: the maximum log
  likelihood it gives -- each table fitted to the data -- less `penalty`
  per free parameter, by default BIC's half log N."
  ([domains parents data] (structure-score domains parents data (* 0.5 (math/log (count data)))))
  ([domains parents data penalty]
   (let [net (learn-parameters {:variables (vec (sort (keys domains))) :parents parents :domains domains} data)
        ;; tables fitted by bare counts give an unseen value 0, but no
        ;; unseen value is ever evaluated: the data is what fitted them
         ll (log-likelihood net data)]
     (- ll (* (double penalty) (parameter-count domains parents))))))

(defn- cyclic? [parents]
  (letfn [(reaches? [from to seen]
            (some (fn [p] (or (= p to) (and (not (seen p)) (reaches? p to (conj seen p))))) (parents from)))]
    (some (fn [v] (reaches? v v #{})) (keys parents))))

(defn- moves
  "Every network one link away from `parents`: a link added, deleted or
  reversed, no cycles, at most `max-parents` parents a variable."
  [variables parents max-parents]
  (for [a variables b variables :when (not= a b)
        :let [has? (some #{a} (parents b))]
        [move ps] (if has?
                    [[[:delete a b] (update parents b #(vec (remove #{a} %)))]
                     [[:reverse a b] (-> parents (update b #(vec (remove #{a} %))) (update a conj b))]]
                    [[[:add a b] (update parents b conj a)]])
        :when (and (<= (count (ps b)) (long max-parents)) (<= (count (ps a)) (long max-parents)) (not (cyclic? ps)))]
    [move ps]))

(defn- random-structure
  "A random network over `variables`: a random order, and each variable
  given each earlier one as a parent with probability 0.3, at most
  `max-parents` -- acyclic by construction."
  [variables max-parents rng]
  (let [order (random/shuffle rng variables)]
    (into {} (map-indexed (fn [i v] [v (vec (take max-parents (filter (fn [_] (< (double (rng)) 0.3)) (subvec order 0 i))))])
                          order))))

(defn structure-search-seq
  "Hill climbing over network structures for `data` over `domains`: each
  step makes the one move -- a link added, deleted or reversed -- that
  raises the score most, until none does. Climbing stops at a local
  maximum, so with `:restarts` n (and an `:rng`) it climbs again from n
  random networks after the first climb from no links, as the chapter
  suggests. The score is the fitted log likelihood less `:penalty` per
  free parameter: by default BIC's, half the log of the number of
  examples. Each state `{:parents :score :move :restart :best}`, `:best`
  the best `{:parents :score}` met so far."
  [domains data & {:keys [penalty max-parents restarts rng] :or {max-parents 3 restarts 0}}]
  (let [variables (vec (sort (keys domains)))
        penalty (if (some? penalty) penalty (* 0.5 (math/log (count data))))
        memo (atom {})
        score-of (fn [ps] (or (@memo ps) (let [s (structure-score domains ps data penalty)] (swap! memo assoc ps s) s)))
        starts (cons (zipmap variables (repeat [])) (repeatedly restarts #(random-structure variables max-parents rng)))]
    (letfn [(climb [{:keys [parents score] :as st}]
              (lazy-seq
               (cons st
                     (let [[move ps s] (reduce (fn [best [m ps]] (let [s (score-of ps)] (if (> s (nth best 2)) [m ps s] best)))
                                               [nil nil score] (moves variables parents max-parents))]
                       (when move (climb {:parents ps :score s :move move}))))))]
      (rest (reductions (fn [{:keys [best]} st]
                          (assoc st :best (if (or (nil? best) (> (double (:score st)) (double (:score best))))
                                            (select-keys st [:parents :score])
                                            best)))
                        {:best nil}
                        (mapcat (fn [k start] (map #(assoc % :restart k) (climb {:parents start :score (score-of start) :move nil})))
                                (range) starts))))))

(defn learn-structure
  "The best network structure `structure-search-seq` finds: `{:parents
  :score}`."
  [domains data & opts]
  (:best (last (apply structure-search-seq domains data opts))))

;; ---------------------------------------------------------------------------
;; Nonparametric density estimation

(defn- distance [a b] (math/sqrt (reduce + (map #(let [d (- %1 %2)] (* d d)) a b))))

(defn- ball-volume
  "The volume of a ball of radius `r` in `d` dimensions."
  [d r]
  (* (/ (math/pow math/PI (/ d 2.0)) (math/exp (special/log-gamma (inc (/ d 2.0))))) (math/pow r d)))

(defn knn-density
  "The density at `x` by its `k` nearest of `points`: k over N times the
  volume of the smallest ball round x that holds them."
  [points k x]
  (let [r (nth (sort (map #(distance x %) points)) (dec (long k)))]
    (/ (double k) (* (count points) (ball-volume (count x) r)))))

(defn kernel-density
  "The density at `x` as the average of a spherical Gaussian of standard
  deviation `w` round each of `points`."
  [points w x]
  (let [d (count x) w (double w)
        norm (math/pow (* w w 2.0 math/PI) (/ d 2.0))]
    (/ (reduce + (map #(math/exp (- (/ (let [r (distance x %)] (* r r)) (* 2.0 w w)))) points))
       (* (count points) norm))))

;; ---------------------------------------------------------------------------
;; Hidden variables: EM for a mixture of naive Bayes models

(defn- class-posteriors
  "P(class | features) for one example, and P(features)."
  [{:keys [prior features]} x]
  (let [joint (mapv (fn [p fs] (reduce * p (map (fn [theta v] (if v theta (- 1.0 theta))) fs x))) prior features)
        total (reduce + joint)]
    [(mapv #(/ % total) joint) total]))

(defn mixture-log-likelihood
  "The log likelihood of `data` -- `[[features count] ...]`, the
  features Booleans -- under the mixture `params`."
  [params data]
  (reduce + (map (fn [[x n]] (* (double n) (math/log (second (class-posteriors params x))))) data)))

(defn naive-bayes-em-seq
  "EM for a mixture of naive Bayes models with Boolean features, the class
  hidden: `params` `{:prior [P(class)] :features [[P(feature true |
  class) ...] ...]}`, `data` `[[features count] ...]`. The E step gives
  each example's posterior over the classes, by Bayes' rule and the
  features' independence given the class; the M step makes each
  parameter its expected count over its class's. An iteration a step,
  each state `{:params :log-likelihood}`, which never falls."
  [params data]
  (let [n (double (reduce + (map second data)))
        state (fn [ps] {:params ps :log-likelihood (mixture-log-likelihood ps data)})]
    (iterate (fn [{:keys [params]}]
               (let [posts (map (fn [[x c]] [x c (first (class-posteriors params x))]) data)
                     k (count (:prior params))
                     totals (mapv (fn [j] (reduce + (map (fn [[_ c post]] (* c (post j))) posts))) (range k))]
                 (state {:prior (mapv #(/ % n) totals)
                         :features (mapv (fn [j]
                                           (mapv (fn [f] (/ (reduce + (map (fn [[x c post]] (if (x f) (* c (post j)) 0.0)) posts)) (totals j)))
                                                 (range (count (ffirst data)))))
                                         (range k))})))
             (state params))))
