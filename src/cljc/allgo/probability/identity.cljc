(ns allgo.probability.identity
  "Identity uncertainty: which observations are of the same unknown
  object, when how many objects there are is unknown too. Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  section 15.2.3's citation matching -- are two citation strings, written
  differently and with errors, citations of the same paper?

  An open-universe model, simplified from the book's figure 15.5:

    papers              how many is unknown: the citations fall into
                        papers by the Chinese restaurant process, each
                        citation of a paper already cited in proportion
                        to its citations so far, of a new one in
                        proportion to `alpha`
    a paper's words     each word of the vocabulary in its title or
                        author list with probability `q`, independently
    a citation's words  each of its paper's words kept with probability
                        `keep`; each word not its paper's put in by
                        error with probability `noise`

  Given the grouping, the papers' words sum out word by word, so a
  group's likelihood is exact and cheap; the grouping itself is
  uncertain. `exact-posterior` enumerates every grouping (the Bell number
  of them, so for a handful of citations); `gibbs-seq` samples groupings,
  each step moving one citation to the paper, old or new, it most
  probably belongs to -- the MCMC the book uses for such models. The
  answer to a query like \"are citations i and j of the same paper?\" is
  the posterior probability of the groupings in which they share a
  group.

  A citation is the words of its author and title fields, the parse
  taken as given: the book's model parses as part of its inference, and
  its figure's paper has only a title and an author, so the venue plays
  no part. Whole strings, venue and all, mislead -- as the section says
  text similarity does.

  As the section says, the inference is collective: a third citation of
  a paper makes the first two more surely the same."
  (:require [allgo.random :as random]
            [clojure.math :as math]
            [clojure.string :as str]))

(defn words
  "The words of a citation string: lowercase, letters and digits only."
  [s]
  (set (re-seq #"[a-z0-9]+" (str/lower-case s))))

(defn- log-sum-exp [xs]
  (let [m (reduce max ##-Inf xs)]
    (if (= m ##-Inf) m (+ m (math/log (reduce + (map #(math/exp (- % m)) xs)))))))

(defn group-log-likelihood
  "log p(the citations in `group` | they cite one paper), the paper's
  words summed out: for each word of the vocabulary, the paper has it
  (probability q) and each citation kept it or not, or lacks it and each
  citation put it in by error or not."
  [{:keys [q keep noise vocabulary]} citations group]
  (reduce + (for [w vocabulary
                  :let [has (count (filter #(contains? (citations %) w) group))
                        lacks (- (count group) has)]]
              (log-sum-exp [(+ (math/log q) (* has (math/log keep)) (* lacks (math/log (- 1.0 keep))))
                            (+ (math/log (- 1.0 q)) (* has (math/log noise)) (* lacks (math/log (- 1.0 noise))))]))))

(defn- crp-log-prior
  "log p(grouping) under the Chinese restaurant process: alpha^k
  prod (n_i - 1)! / (alpha (alpha + 1) ... (alpha + n - 1)) for k groups
  of sizes n_i."
  [alpha groups]
  (let [n (reduce + (map count groups))
        log-fact (fn [m] (reduce + 0.0 (map #(math/log %) (range 1 (inc m)))))]
    (- (+ (* (count groups) (math/log alpha)) (reduce + (map #(log-fact (dec (count %))) groups)))
       (reduce + (map #(math/log (+ alpha %)) (range n))))))

(defn model
  "A model over the word sets `citations`: `:alpha` the Chinese
  restaurant process's concentration, `:q` `:keep` `:noise` as above, the
  vocabulary every word cited."
  [citations {:keys [alpha q keep noise] :or {alpha 1.0 q 0.1 keep 0.8 noise 0.01}}]
  {:alpha alpha :q q :keep keep :noise noise :vocabulary (reduce into #{} citations)})

(defn log-joint
  "log p(grouping, citations): the prior times each group's likelihood."
  [m citations groups]
  (+ (crp-log-prior (:alpha m) groups) (reduce + (map #(group-log-likelihood m citations %) groups))))

(defn- partitions
  "Every way to split the items `xs` into nonempty groups."
  [xs]
  (if (empty? xs)
    [[]]
    (let [[x & more] xs]
      (for [p (partitions more)
            option (cons (conj p [x]) (map-indexed (fn [i g] (assoc p i (conj g x))) p))]
        option))))

(defn exact-posterior
  "The posterior over every grouping of the citations, as a map from a
  set of groups (each a set of citation indices) to its probability."
  [m citations]
  (let [groups (map (fn [p] (set (map set p))) (partitions (range (count citations))))
        lj (map #(log-joint m citations %) groups)
        z (log-sum-exp lj)]
    (zipmap groups (map #(math/exp (- % z)) lj))))

(defn same-paper
  "The probability under `posterior` (groupings to probabilities) that
  citations `i` and `j` are of the same paper."
  [posterior i j]
  (reduce + (for [[groups p] posterior :when (some #(and (% i) (% j)) groups)] p)))

(defn- gibbs-step
  "One sweep: each citation in turn taken out of its group and put back
  into an existing group or a new one of its own, in proportion to the
  joint probability of the grouping that results."
  [m citations rng groups]
  (reduce (fn [groups i]
            (let [rest-groups (vec (remove empty? (map #(disj % i) groups)))
                  options (conj (mapv (fn [k] (update rest-groups k conj i)) (range (count rest-groups)))
                                (conj rest-groups #{i}))
                  lj (mapv #(log-joint m citations %) options)
                  z (log-sum-exp lj)]
              (options (random/weighted rng (zipmap (range (count options)) (map #(math/exp (- % z)) lj))))))
          groups
          (range (count citations))))

(defn gibbs-seq
  "Groupings sampled by Gibbs sweeps from every citation alone, each a
  vector of sets of citation indices."
  [m citations rng]
  (iterate (partial gibbs-step m citations rng) (mapv hash-set (range (count citations)))))

(defn sampled-posterior
  "The posterior over groupings estimated from `samples`."
  [samples]
  (update-vals (frequencies (map set samples)) #(/ (double %) (count samples))))
