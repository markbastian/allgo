(ns allgo.learning.unsupervised
  "Learning with hidden variables: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 20, whose text
  gives these without pseudocode. Each is an iteration, a step a pass.

    k-means         assign each point to its nearest center, move each
                    center to its points' mean
    EM for          the E-step weighs each point's membership of each
    Gaussian        Gaussian; the M-step refits each Gaussian's weight,
    mixtures        mean and (diagonal) variances to them
    Baum-Welch      EM for a hidden Markov model of
                    `allgo.probability.temporal`: smoothing gives the
                    expected transitions and emissions, which become the
                    new model"
  (:require [allgo.probability.temporal :as temporal]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn- sq-dist [a b] (reduce + (map #(let [d (- %1 %2)] (* d d)) a b)))
(defn- mean [xs] (let [n (count xs)] (apply mapv (fn [& vs] (/ (reduce + vs) n)) xs)))

(defn k-means-seq
  "k-means on `points` from `k` of them chosen at random: each state
  `{:centers :assignment :cost}`, cost the sum of squared distances to
  the centers -- which never rises."
  [points k rng]
  (let [points (vec points)]
    (iterate (fn [{:keys [centers]}]
               (let [assignment (mapv (fn [p] (apply min-key #(sq-dist p (centers %)) (range k))) points)
                     groups (group-by assignment (range (count points)))
                     centers (mapv (fn [c] (if-let [is (groups c)] (mean (map points is)) (centers c))) (range k))]
                 {:centers centers :assignment assignment
                  :cost (reduce + (map (fn [p a] (sq-dist p (centers a))) points assignment))}))
             {:centers (vec (take k (random/shuffle rng points))) :assignment nil :cost nil})))

(defn- gaussian-density [x mu var]
  (reduce * (map (fn [xi m v] (/ (math/exp (- (/ (* (- xi m) (- xi m)) (* 2 v)))) (math/sqrt (* 2 math/PI v)))) x mu var)))

(defn em-gmm-seq
  "EM for a mixture of `k` Gaussians with diagonal covariances, from
  centers chosen at random: each state `{:weights :means :vars
  :log-likelihood}` -- the log likelihood never falls."
  [points k rng]
  (let [points (vec points) n (count points) d (count (first points))
        all-var (mapv (fn [j] (/ (reduce + (map #(let [v (- (nth % j) (nth (mean points) j))] (* v v)) points)) n)) (range d))]
    (iterate (fn [{:keys [weights means vars]}]
               (let [;; E-step: each point's responsibilities.
                     resp (mapv (fn [x] (let [ps (mapv (fn [c] (* (weights c) (gaussian-density x (means c) (vars c)))) (range k))
                                              total (reduce + ps)]
                                          (if (pos? total) (mapv #(/ % total) ps) (vec (repeat k (/ 1.0 k))))))
                                points)
                     ;; M-step.
                     nk (mapv (fn [c] (reduce + (map #(% c) resp))) (range k))
                     means (mapv (fn [c] (mapv (fn [j] (/ (reduce + (map (fn [x r] (* (r c) (nth x j))) points resp)) (max 1e-12 (nk c)))) (range d))) (range k))
                     vars (mapv (fn [c] (mapv (fn [j] (+ 1e-6 (/ (reduce + (map (fn [x r] (let [v (- (nth x j) ((means c) j))] (* (r c) v v))) points resp))
                                                                 (max 1e-12 (nk c)))))
                                              (range d)))
                                (range k))
                     weights (mapv #(/ % n) nk)]
                 {:weights weights :means means :vars vars
                  :log-likelihood (reduce + (map (fn [x] (math/log (reduce + (map (fn [c] (* (weights c) (gaussian-density x (means c) (vars c)))) (range k))))) points))}))
             {:weights (vec (repeat k (/ 1.0 k))) :means (vec (take k (random/shuffle rng points)))
              :vars (vec (repeat k all-var)) :log-likelihood nil})))

(defn log-likelihood
  "log P(evidence) under the hidden Markov model."
  [{:keys [prior sensor] :as hmm} evidence]
  (first (reduce (fn [[ll f] e]
                   (let [pred (temporal/forward hmm f e)
                         ;; The normalizer of the forward step.
                         raw (let [{:keys [transition]} hmm n (count f)]
                               (mapv (fn [j] (* ((sensor e) j) (reduce + (map (fn [i] (* (f i) (get-in transition [i j]))) (range n))))) (range n)))]
                     [(+ ll (math/log (reduce + raw))) pred]))
                 [0.0 prior] evidence)))

(defn baum-welch-seq
  "Baum-Welch for a hidden Markov model given one sequence of `evidence`:
  each step re-estimates the transition and sensor models (and the prior)
  from the expected counts that smoothing gives. The log likelihood of
  the evidence never falls; `log-likelihood` measures it."
  [hmm evidence]
  (let [values (vec (distinct evidence))
        t (count evidence)]
    (iterate
     (fn [{:keys [prior transition sensor] :as m}]
       (let [n (count prior)
             fv (vec (temporal/filtering-seq m evidence))
             ;; Backward messages b_k for k = 1..t, b_t = 1.
             bv (vec (reverse (reductions (fn [b e] (temporal/backward m b e)) (vec (repeat n 1.0)) (reverse (rest evidence)))))
             gamma (mapv (fn [k] (let [g (mapv * (fv (inc k)) (bv k)) s (reduce + g)] (mapv #(/ % s) g))) (range t))
             ;; Expected transitions, from X0 -> X1 (k = -1) to X_t-1 -> X_t:
             ;; the model's first transition comes before any evidence.
             xi (vec (for [k (range -1 (dec t))]
                       (let [e (nth evidence (inc k))
                             raw (vec (for [i (range n)] (vec (for [j (range n)] (* ((fv (inc k)) i) (get-in transition [i j]) ((sensor e) j) ((bv (inc k)) j))))))
                             s (reduce + (flatten raw))]
                         (mapv (fn [row] (mapv #(/ % s) row)) raw))))
             transition' (vec (for [i (range n)]
                                (let [row (mapv (fn [j] (reduce + (map #(get-in % [i j]) xi))) (range n))
                                      s (reduce + row)]
                                  (if (pos? s) (mapv #(/ % s) row) (transition i)))))
             sensor' (into {} (for [v values]
                                [v (mapv (fn [i] (/ (reduce + (map (fn [g e] (if (= e v) (g i) 0.0)) gamma evidence))
                                                    (reduce + (map #(% i) gamma))))
                                         (range n))]))]
         ;; The prior is P(X0), before any evidence: smoothed, it is
         ;; P(X0) times the backward message from step 1.
         (assoc m :prior (let [g (mapv * prior (temporal/backward m (bv 0) (first evidence))) s (reduce + g)] (mapv #(/ % s) g))
                :transition transition' :sensor sensor')))
     hmm)))
