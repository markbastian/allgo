(ns allgo.ethics.privacy
  "Privacy: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, section 27.3.2 -- what can be learned from data
  about people without learning about any one of them.

    k-anonymity        records stripped of names can still be picked out by
                       what is left -- birth date, sex, zip code -- so those
                       fields are generalized, a year for a date, a range
                       for an age, until every record looks like at least
                       k - 1 others; a step a generalization of whichever
                       field still splits the records most
    aggregate queries  the records kept private and only summaries
                       answered, none for a group smaller than n -- which
                       two overlapping queries can still undo: the
                       difference of two averages is one person's value
    differential       noise added to each answer, Laplace noise of scale
    privacy            sensitivity / epsilon, so that one person's being in
                       the data or not changes the log probability of any
                       answer by at most epsilon
    federated          no central data at all: each user trains on their
    learning           own and sends only parameters, which the server
                       averages -- here logistic regression, a round of
                       local gradient steps then the average
    secure             the server learning only the average: each pair of
    aggregation        users agrees a random mask, one adding it and the
                       other subtracting, so every upload looks random and
                       the masks cancel in the sum

  A record is a map of field to value."
  (:require [allgo.learning.models :as models]
            [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; k-anonymity

(defn classes
  "The records grouped by their values of the `quasi` fields -- those an
  attacker could link to another database -- `{values [records]}`."
  [records quasi]
  (group-by (fn [r] (mapv r quasi)) records))

(defn anonymity
  "The k for which `records` are k-anonymous over the `quasi` fields: the
  size of the smallest group of records that agree on them all."
  [records quasi]
  (reduce min (map count (vals (classes records quasi)))))

(defn- outliers
  "The records in groups smaller than `k` over the `quasi` fields."
  [records quasi k]
  (mapcat identity (filter #(< (count %) k) (vals (classes records quasi)))))

(defn- release
  "The `records` at `levels`, done if at most `suppress` stand out -- and
  those withheld."
  [records levels quasi k suppress]
  (let [out (outliers records quasi k)]
    (if (<= (count out) suppress)
      {:records (vec (remove (set out) records)) :levels levels :done true}
      {:records records :levels levels})))

(defn- anonymize-step
  "A rung climbed on the field with the most distinct values left, the
  `original` records generalized to the new levels."
  [original hierarchies k suppress {:keys [records levels done] :as st}]
  (if done
    st
    (let [quasi (vec (keys hierarchies))
          open (filter #(< (levels %) (count (hierarchies %))) quasi)
          f (apply max-key #(count (distinct (map % records))) open)
          levels (update levels f inc)
          generalized (mapv (fn [r] (reduce (fn [r f]
                                              (let [l (levels f)]
                                                (if (pos? l) (update r f (nth (hierarchies f) (dec l))) r)))
                                            r quasi))
                            original)]
      (release generalized levels quasi k suppress))))

(defn anonymize-seq
  "Generalization until `records` are `k`-anonymous over the fields of
  `hierarchies`, a step a generalization: each field has a ladder of
  functions from a value to a coarser one, the last mapping everything to
  `*` (the field deleted); a step climbs one rung of the field with the
  most distinct values left. As in Sweeney's Datafly (\"Achieving
  k-Anonymity Privacy Protection Using Generalization and Suppression\",
  2002), up to `suppress` records (default 0) that still stand out may be
  withheld rather than every record generalized further for their sake.
  Each state `{:records :levels :done}`, `:records` the ones released."
  ([records hierarchies k] (anonymize-seq records hierarchies k 0))
  ([records hierarchies k suppress]
   (let [quasi (vec (keys hierarchies))]
     (iterate (partial anonymize-step records hierarchies k suppress)
              (release (vec records) (zipmap quasi (repeat 0)) quasi k suppress)))))

(defn anonymize
  "The released records `anonymize-seq` ends with, and the levels it
  reached: `{:records :levels}`."
  ([records hierarchies k] (anonymize records hierarchies k 0))
  ([records hierarchies k suppress]
   (dissoc (first (filter :done (anonymize-seq records hierarchies k suppress))) :done)))

;; ---------------------------------------------------------------------------
;; Aggregate queries

(defn aggregate
  "The mean of `field` over the records `where` selects, and how many
  they are -- `[mean count]` -- or nil if fewer than `n`: a summary is
  given only for a group big enough to hide in."
  [records where field n]
  (let [rs (filter where records)]
    (when (>= (count rs) n)
      [(/ (reduce + (map field rs)) (double (count rs))) (count rs)])))

(defn difference-attack
  "One person's value from two aggregate answers whose groups differ by
  that person: `[mean count]` for the larger and the smaller. The two
  totals' difference."
  [[mean-big n-big] [mean-small n-small]]
  (- (* mean-big n-big) (* mean-small n-small)))

;; ---------------------------------------------------------------------------
;; Differential privacy

(defn laplace
  "A draw from the Laplace distribution of location 0 and scale `b`."
  [rng b]
  (let [u (- (double (rng)) 0.5)]
    (* -1.0 b (math/signum u) (math/log (- 1.0 (* 2.0 (abs u)))))))

(defn laplace-mechanism
  "The `value` of a query with Laplace noise of scale `sensitivity` /
  `epsilon` added: epsilon-differentially private when one person can
  change the query's answer by at most `sensitivity`."
  [value sensitivity epsilon rng]
  (+ value (laplace rng (/ sensitivity epsilon))))

(defn private-count
  "How many `records` satisfy `pred`, with epsilon-differential privacy:
  one person changes a count by at most 1."
  [records pred epsilon rng]
  (laplace-mechanism (count (filter pred records)) 1.0 epsilon rng))

(defn laplace-log-density
  "The log density of the Laplace mechanism's answer `y` when the true
  value is `value`."
  [value sensitivity epsilon y]
  (let [b (/ sensitivity epsilon)]
    (- (- (/ (abs (- y value)) b)) (math/log (* 2 b)))))

;; ---------------------------------------------------------------------------
;; Federated learning and secure aggregation

(defn local-update
  "A user's model after `steps` epochs of gradient descent on their own
  data `[X ys]` from the shared weights `w`: logistic regression,
  `allgo.learning.models`'."
  [w [X ys] alpha steps]
  (nth (models/logistic-regression-seq X ys alpha w) steps))

(defn- federated-round
  "Every user's local update of `w`, averaged by how much data each
  holds."
  [datasets alpha local-steps w]
  (let [total (reduce + (map (comp count second) datasets))
        updates (map #(local-update w % alpha local-steps) datasets)]
    (apply mapv + (map (fn [u [_ ys]] (mapv #(* % (/ (count ys) (double total))) u)) updates datasets))))

(defn federated-seq
  "Federated averaging (after McMahan et al., \"Communication-Efficient
  Learning of Deep Networks from Decentralized Data\", AISTATS 2017): a
  round a step, every user in `datasets` trains the shared weights on
  their own data and sends back only the weights, and the server averages
  them, weighted by how much data each holds. Each state the shared
  weights."
  [datasets dims alpha local-steps]
  (iterate (partial federated-round datasets alpha local-steps) (vec (repeat (inc dims) 0.0))))

(defn masked-uploads
  "Secure aggregation's uploads (after Bonawitz et al., \"Practical Secure
  Aggregation for Privacy-Preserving Machine Learning\", CCS 2017, without
  its dropout handling): for every pair of users i < j a random mask
  vector, drawn from a seed the pair shares, that i adds and j subtracts.
  Each user's upload alone is noise; the uploads' sum is the vectors'."
  [vectors seed scale]
  (let [n (count vectors) d (count (first vectors))
        mask (fn [i j] (let [rng (random/rng (+ (* 1000003 seed) (* 1009 i) j))]
                         (vec (repeatedly d #(* scale (- (* 2.0 (rng)) 1.0))))))]
    (vec (for [i (range n)]
           (reduce (fn [v j]
                     (cond (< i j) (mapv + v (mask i j))
                           (> i j) (mapv - v (mask j i))
                           :else v))
                   (vectors i) (range n))))))
