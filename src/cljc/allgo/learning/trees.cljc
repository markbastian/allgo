(ns allgo.learning.trees
  "Learning from examples, by trees and their ensembles: Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 19.

    decision trees      DECISION-TREE-LEARNING (figure 19.5), choosing
                        each test by information gain, with an optional
                        depth limit
    model selection     MODEL-SELECTION and CROSS-VALIDATION (figure
                        19.8): k-fold cross-validation over sizes, keeping
                        the size with the least validation error
    decision lists      DECISION-LIST-LEARNING (figure 19.11), tests of
                        up to `k` attribute values
    AdaBoost            ADABOOST (figure 19.25), a round a step
    random forests      the chapter's bagged trees, each split chosen from
                        a random subset of the attributes

  An example is a map from attributes to values, the classification under
  `target`. A tree is a leaf value, or `{:attribute :branches {value
  subtree} :default}`."
  (:require [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Decision trees

(defn plurality-value
  "The most common classification, ties broken by first appearance."
  [examples target]
  (let [counts (frequencies (map target examples))]
    (apply max-key counts (distinct (map target examples)))))

(defn entropy
  "The entropy, in bits, of the classifications of `examples`."
  [examples target]
  (let [n (count examples)]
    (- (reduce + (for [[_ c] (frequencies (map target examples)) :let [p (/ (double c) n)]]
                   (* p (/ (math/log p) (math/log 2))))))))

(defn information-gain
  "The entropy of `examples` less the expected entropy after testing `a`."
  [a examples target]
  (let [n (double (count examples))]
    (- (entropy examples target)
       (reduce + (for [[_ exs] (group-by a examples)] (* (/ (count exs) n) (entropy exs target)))))))

(defn decision-tree-learning
  "DECISION-TREE-LEARNING: a tree that classifies `examples` by `target`,
  testing first the attribute of greatest information gain. `values` maps
  each attribute to all its values, so a branch exists even for a value
  no example has, labeled with its parent's plurality. `:max-depth`
  limits the depth."
  ([examples attributes target values] (decision-tree-learning examples attributes target values {}))
  ([examples attributes target values {:keys [max-depth choose]}]
   (letfn [(learn [exs attrs parent depth]
             (cond
               (empty? exs) (plurality-value parent target)
               (apply = (map target exs)) (target (first exs))
               (or (empty? attrs) (and max-depth (>= depth max-depth))) (plurality-value exs target)
               :else
               (let [candidates (if choose (choose attrs) attrs)
                     a (apply max-key #(information-gain % exs target) (reverse candidates))]
                 {:attribute a
                  :default (plurality-value exs target)
                  :branches (into {} (for [v (values a)]
                                       [v (learn (filter #(= v (a %)) exs) (remove #{a} attrs) exs (inc depth))]))})))]
     (learn examples (vec attributes) examples 0))))

(defn classify
  "The classification `tree` gives `example`."
  [tree example]
  (if (map? tree)
    (classify (get (:branches tree) ((:attribute tree) example) (:default tree)) example)
    tree))

(defn tree-size
  "The number of tests in `tree`."
  [tree]
  (if (map? tree) (inc (reduce + (map tree-size (vals (:branches tree))))) 0))

(defn error-rate
  "The fraction of `examples` hypothesis `h` gets wrong."
  [h examples target]
  (if (empty? examples) 0.0 (/ (double (count (remove #(= (target %) (h %)) examples))) (count examples))))

;; ---------------------------------------------------------------------------
;; Model selection

(defn partition-folds
  "The `k` folds of `examples`, as `[training validation]` pairs."
  [examples k]
  (let [v (vec examples) n (count v)]
    (for [fold (range k)
          :let [lo (quot (* fold n) k) hi (quot (* (inc fold) n) k)]]
      [(into (subvec v 0 lo) (subvec v hi)) (subvec v lo hi)])))

(defn cross-validation
  "CROSS-VALIDATION: the mean training and validation error rates of
  `learner`, a function of a size and training examples returning a
  hypothesis, over `k` folds: `[train-error validation-error]`."
  [learner size examples k target]
  (let [errs (for [[train valid] (partition-folds examples k)
                   :let [h (learner size train)]]
               [(error-rate h train target) (error-rate h valid target)])]
    [(/ (reduce + (map first errs)) k) (/ (reduce + (map second errs)) k)]))

(defn model-selection
  "MODEL-SELECTION: the size, from 1 to `max-size`, of least validation
  error by `k`-fold cross-validation, and the learner retrained at that
  size on every example: `{:size :hypothesis :errors}`, `:errors` each
  size's `[train validation]`."
  [learner examples k target max-size]
  (let [errors (into (sorted-map) (for [size (range 1 (inc max-size))] [size (cross-validation learner size examples k target)]))
        best (apply min-key #(second (errors %)) (keys errors))]
    {:size best :hypothesis (learner best examples) :errors errors}))

;; ---------------------------------------------------------------------------
;; Decision lists

(defn- tests
  "Every conjunction of up to `k` attribute = value literals."
  [values attributes k]
  (letfn [(combos [attrs k]
            (if (zero? k)
              [[]]
              (concat [[]]
                      (for [[i a] (map-indexed vector attrs) v (values a) rest (combos (drop (inc i) attrs) (dec k))]
                        (cons [a v] rest)))))]
    (distinct (remove empty? (combos attributes k)))))

(defn decision-list-learning
  "DECISION-LIST-LEARNING: tests of up to `k` attribute values, each with
  the outcome its examples all share, the smallest test first that
  covers a nonempty set of examples of one class -- or nil if there is
  none. A decision list is `[[test outcome] ...]`, a test a seq of
  `[attribute value]`; the list ends with the `default`."
  [examples attributes target values k]
  (let [ts (sort-by count (tests values attributes k))
        matches? (fn [t e] (every? (fn [[a v]] (= v (a e))) t))]
    (loop [exs examples dl []]
      (if (empty? exs)
        dl
        (if-let [[t covered] (some (fn [t] (let [c (filter #(matches? t %) exs)]
                                             (when (and (seq c) (apply = (map target c))) [t c])))
                                   ts)]
          (recur (remove (set covered) exs) (conj dl [t (target (first covered))]))
          nil)))))

(defn decision-list-classify [dl default example]
  (or (some (fn [[t o]] (when (every? (fn [[a v]] (= v (a example))) t) o)) dl) default))

;; ---------------------------------------------------------------------------
;; AdaBoost

(defn- weighted-stump
  "The one-test tree of least weighted error: the weighted learner AdaBoost
  boosts here."
  [examples weights attributes target values]
  (let [stump (fn [a]
                (let [by (group-by #(a (first %)) (map vector examples weights))
                      branch (into {} (for [v (values a)]
                                        [v (let [ws (by v)]
                                             (if (seq ws)
                                               (key (apply max-key val (reduce (fn [m [e w]] (update m (target e) (fnil + 0.0) w)) {} ws)))
                                               (plurality-value examples target)))]))]
                  {:attribute a :branches branch :default (plurality-value examples target)}))
        err (fn [tree] (reduce + (map (fn [e w] (if (= (target e) (classify tree e)) 0.0 w)) examples weights)))]
    (apply min-key err (map stump (reverse attributes)))))

(defn adaboost-seq
  "ADABOOST, a round a step: each round a hypothesis from the weighted
  learner (default a decision stump), its weight log((1 - error) /
  error), and the examples it got right made lighter. Each state has the
  `:hypotheses` and their `:z` so far."
  ([examples attributes target values] (adaboost-seq examples attributes target values {}))
  ([examples attributes target values {:keys [learner]}]
   (let [n (count examples)
         learner (or learner #(weighted-stump %1 %2 attributes target values))]
     (iterate (fn [{:keys [w hypotheses z]}]
                (let [h (learner examples w)
                      error (min (- 1.0 1e-10) (max 1e-10 (reduce + (map (fn [e wj] (if (= (target e) (classify h e)) 0.0 wj)) examples w))))
                      w (mapv (fn [e wj] (if (= (target e) (classify h e)) (* wj (/ error (- 1.0 error))) wj)) examples w)
                      total (reduce + w)]
                  {:w (mapv #(/ % total) w)
                   :hypotheses (conj hypotheses h)
                   :z (conj z (math/log (/ (- 1.0 error) error)))}))
              {:w (vec (repeat n (/ 1.0 n))) :hypotheses [] :z []}))))

(defn weighted-majority
  "WEIGHTED-MAJORITY: the class the hypotheses vote for, each vote weighted
  by its z."
  [{:keys [hypotheses z]} example]
  (key (apply max-key val (reduce (fn [m [h zk]] (update m (classify h example) (fnil + 0.0) zk)) {} (map vector hypotheses z)))))

;; ---------------------------------------------------------------------------
;; Random forests

(defn random-forest
  "`n` trees, each learned from a bootstrap sample of `examples` with
  every split chosen among `m` attributes drawn at random (default the
  square root of their number): `(fn [example]) -> class`, by vote."
  [examples attributes target values n rng & {:keys [m]}]
  (let [m (or m (max 1 (int (math/round (math/sqrt (count attributes))))))
        v (vec examples)
        trees (vec (for [_ (range n)]
                     (let [sample (vec (repeatedly (count v) #(random/pick rng v)))]
                       (decision-tree-learning sample attributes target values
                                               {:choose (fn [attrs] (take m (random/shuffle rng (vec attrs))))}))))]
    (fn [example]
      (key (apply max-key val (frequencies (map #(classify % example) trees)))))))
