(ns allgo.ethics.fairness-test
  (:require [allgo.ethics.fairness :as fair]
            [allgo.learning.models :as models]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Fairness and bias, after Russell and Norvig, *Artificial Intelligence:
;; A Modern Approach*, 4th edition, section 27.3.3.

(deftest criteria-test
  (let [xs [{:g :a :pred true :label true} {:g :a :pred false :label true} {:g :a :pred true :label false} {:g :a :pred false :label false}
            {:g :b :pred true :label true} {:g :b :pred true :label true} {:g :b :pred false :label false} {:g :b :pred false :label false}]
        r (fair/group-rates xs :g :pred :label)]
    (testing "both groups decided for at the same rate, a half: demographic parity"
      (is (zero? (fair/gap r :selection))))
    (testing "but not equal opportunity: group a's deserving get it half the time, b's always"
      (is (= 0.5 (get-in r [:a :tpr])))
      (is (= 1.0 (get-in r [:b :tpr])))
      (is (= 0.5 (get-in r [:a :fpr])))
      (is (= 0.0 (get-in r [:b :fpr]))))))

(deftest kleinberg-test
  (let [xs (fair/kleinberg-scores {:a 0.3 :b 0.6} 4000 (random/rng 1))]
    (testing "a calibrated score: within a score bin both groups turn out positive at about the same rate"
      (let [c (fair/calibration xs :group :score :label 0.1)]
        (doseq [b [0.4 0.5]] (is (< (abs (- (get-in c [:a b]) (get-in c [:b b]))) 0.06)))))
    (testing "but the base rates differ, so one threshold wrongly flags the group with the higher rate far more often"
      (let [r (fair/group-rates xs :group (fair/threshold-decision :score :group 0.5) :label)]
        (is (> (fair/gap r :fpr) 0.4))))
    (testing "thresholds set per group to equalize the true positive rate still leave the false positive rates apart: the impossibility"
      (let [th (fair/equal-opportunity-thresholds xs :group :score :label 0.8)
            r (fair/group-rates xs :group (fair/threshold-decision :score :group th) :label)]
        (is (< (fair/gap r :tpr) 0.01))
        (is (> (fair/gap r :fpr) 0.05)))))
  (testing "equal impact: decide for whoever's chance of the good outcome makes it worth it, cost / (benefit + cost)"
    (is (= 0.25 (fair/utility-threshold 3.0 1.0)))))

(deftest unawareness-test
  (testing "the group deleted, zip code still gives it away: fairness through unawareness fails"
    (let [rng (random/rng 2)
          xs (vec (for [_ (range 400)] (let [g (random/pick rng [:a :b])]
                                         {:group g :zip (random/gaussian rng (if (= g :a) 0.0 3.0) 1.0) :income (random/gaussian rng 5.0 1.0)})))]
      (is (> (fair/proxy-accuracy xs [:zip :income] :group) 0.85))
      (is (< (fair/proxy-accuracy xs [:income] :group) 0.6)))))

(deftest sample-size-disparity-test
  (let [rng (random/rng 3)
        make (fn [g n] (vec (for [_ (range n)] (let [x1 (random/gaussian rng 0.0 1.0) x2 (random/gaussian rng 0.0 1.0)]
                                                 {:x1 x1 :x2 x2 :y (if (= g :major) (> (+ x1 (* 0.3 x2)) 0) (> (- x2 (* 0.3 x1)) 0))}))))
        major (make :major 900) minor (make :minor 100)
        test-major (make :major 1000) test-minor (make :minor 1000)
        fit (fn [train] (nth (models/logistic-regression-seq (mapv (juxt :x1 :x2) train) (mapv #(if (:y %) 1 0) train) 1.0) 300))
        accuracy (fn [w xs] (/ (count (filter #(= (:y %) (> (models/predict-logistic w [(:x1 %) (:x2 %)]) 0.5)) xs)) (double (count xs))))
        pooled (fit (into major minor))
        balanced (fit (into (into major minor) (fair/smote minor [:x1 :x2] 5 800 rng)))]
    (testing "a group with a different rule and a tenth of the data: one linear model fits the majority and fails the minority"
      (is (> (accuracy pooled test-major) 0.95))
      (is (< (accuracy pooled test-minor) 0.6)))
    (testing "SMOTE's synthetic examples give the minority its share of the fit -- at the majority's expense: no linear model fits both"
      (is (> (accuracy balanced test-minor) 0.7))
      (is (< (accuracy balanced test-major) (accuracy pooled test-major))))))
