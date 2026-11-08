(ns allgo.ethics.fairness
  "Fairness and bias: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 27.3.3.

  A decision about people -- a loan, a release before trial -- judged by
  the section's criteria, each comparing groups:

    demographic parity  each group decided for at the same rate
    equal opportunity   of those who truly deserve it, each group decided
                        for at the same rate (equal true positive rates),
                        and, as COMPAS was faulted for, wrongly decided
                        against at the same rate (equal false positives)
    calibration         of those given the same score, the same share in
                        each group turn out positive
    equal impact        the decision made where expected utility says,
                        the same costs and benefits for everyone

  Kleinberg, Mullainathan and Raghavan (2016) proved that when the groups'
  base rates differ, no score both calibrated and equal in opportunity
  exists; `kleinberg-scores` makes calibrated scores for two groups and
  shows the false positive rates part. And two of the section's sources
  of bias: fairness through unawareness fails, since other features
  predict the deleted one (`proxy-accuracy`); and a group with few
  examples is fitted worse, which oversampling it -- SMOTE, new examples
  between its neighbors -- can redress.

  An example is a map; `pred` and `label` are booleans, `group` any
  value, `score` a number."
  (:require [allgo.random :as random]
            [clojure.math :as math]))

(defn- rate [xs pred] (if (seq xs) (/ (count (filter pred xs)) (double (count xs))) ##NaN))

(defn group-rates
  "For each group of `examples` (by `group`): the share decided for, the
  true positive rate and the false positive rate of the decision `pred`
  against the outcome `label`. `{group {:selection :tpr :fpr :n}}`."
  [examples group pred label]
  (into {} (for [[g xs] (group-by group examples)]
             [g {:n (count xs)
                 :selection (rate xs pred)
                 :tpr (rate (filter label xs) pred)
                 :fpr (rate (remove label xs) pred)}])))

(defn gap
  "The largest difference between groups in the rate `k` of `rates`."
  [rates k]
  (let [vs (map k (vals rates))] (- (reduce max vs) (reduce min vs))))

(defn calibration
  "For each group, for each score bin of width `width`: the share of
  examples whose `label` is true. `{group {bin share}}` -- equal across
  groups in each bin if the score is calibrated."
  [examples group score label width]
  (into {} (for [[g xs] (group-by group examples)]
             [g (into (sorted-map) (for [[b ys] (group-by #(* width (math/floor (/ (score %) width))) xs)]
                                     [b (rate ys label)]))])))

(defn kleinberg-scores
  "A calibrated score for two groups whose `base-rates` differ: each
  person's score is drawn from a beta-like spread round their group's
  rate, and the outcome then drawn true with probability the score -- so,
  by construction, a score s means an s chance of the outcome in either
  group. `n` people a group. Each example `{:group :score :label}`."
  [base-rates n rng]
  (vec (for [[g rate] base-rates
             _ (range n)
             :let [s (max 0.01 (min 0.99 (+ rate (* 0.25 (- (* 2.0 (rng)) 1.0)))))]]
         {:group g :score s :label (< (double (rng)) s)})))

(defn threshold-decision
  "A decision by thresholds on a score, one per group (`thresholds`, a
  map of group to threshold) or one for all (a number)."
  [score group thresholds]
  (fn [x] (>= (score x) (if (map? thresholds) (thresholds (group x)) thresholds))))

(defn equal-opportunity-thresholds
  "For each group, the threshold on `score` whose true positive rate is
  closest to `target`: the decisions then give each group's truly
  positive members the same chance."
  [examples group score label target]
  (into {} (for [[g xs] (group-by group examples)]
             (let [pos (filter label xs)
                   ts (sort (distinct (map score xs)))]
               [g (apply min-key #(abs (- target (rate pos (fn [x] (>= (score x) %))))) ts)]))))

(defn utility-threshold
  "Equal impact: decide for a person when their probability of the good
  outcome `p` makes deciding for them worth more than against -- p
  `benefit` - (1 - p) `cost` > 0 -- the threshold cost / (benefit +
  cost), the same for everyone with the same p, whatever their group."
  [benefit cost]
  (/ cost (+ benefit cost)))

;; ---------------------------------------------------------------------------
;; Sources of bias

(defn proxy-accuracy
  "Fairness through unawareness, tested: how well the protected `group`
  can be predicted from the other `features` by 1-nearest-neighbor, left
  out one at a time -- near 1 when, deleted or not, the attribute is
  still there in what correlates with it."
  [examples features group]
  (let [xs (vec examples)
        d (fn [a b] (reduce + (map #(let [u (- (double (a %)) (double (b %)))] (* u u)) features)))]
    (rate (range (count xs))
          (fn [i] (let [x (xs i)
                        nearest (apply min-key #(if (= % i) ##Inf (d x (xs %))) (range (count xs)))]
                    (= (group x) (group (xs nearest))))))))

(defn smote
  "SMOTE (Chawla et al., \"SMOTE: Synthetic Minority Over-sampling
  Technique\", JAIR 16, 2002): `n` new examples for a scarce group, each
  on the segment between one of its `examples` and one of that example's
  `k` nearest neighbors in the group, at a random point -- new examples
  that fill in the group's region rather than repeating its members.
  `features` are the numeric fields interpolated; the rest are copied."
  [examples features k n rng]
  (let [xs (vec examples)
        d (fn [a b] (reduce + (map #(let [u (- (double (a %)) (double (b %)))] (* u u)) features)))]
    (vec (repeatedly n (fn []
                         (let [x (random/pick rng xs)
                               near (take k (sort-by #(d x %) (remove #(identical? % x) xs)))
                               y (random/pick rng (vec near))
                               t (double (rng))]
                           (reduce (fn [z f] (assoc z f (+ (x f) (* t (- (y f) (x f)))))) x features)))))))
