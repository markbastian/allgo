(ns allgo.learning.trees-test
  (:require [allgo.learning.models :as mo]
            [allgo.learning.trees :as tr]
            [allgo.random]
            [clojure.test :refer [deftest is testing]]))

;; The restaurant examples of Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 19, as the
;; aima-data repository's restaurant.csv has them, read by aima-python
;; (MIT): will the customer wait for a table?

(def ^:private attrs [:Alternate :Bar :FriSat :Hungry :Patrons :Price :Raining :Reservation :Type :WaitEstimate])

(def ^:private rows
  (mapv #(zipmap (conj attrs :WillWait) %)
        [["Yes" "No" "No" "Yes" "Some" "$$$" "No" "Yes" "French" "0-10" "Yes"]
         ["Yes" "No" "No" "Yes" "Full" "$" "No" "No" "Thai" "30-60" "No"]
         ["No" "Yes" "No" "No" "Some" "$" "No" "No" "Burger" "0-10" "Yes"]
         ["Yes" "No" "Yes" "Yes" "Full" "$" "No" "No" "Thai" "10-30" "Yes"]
         ["Yes" "No" "Yes" "No" "Full" "$$$" "No" "Yes" "French" ">60" "No"]
         ["No" "Yes" "No" "Yes" "Some" "$$" "Yes" "Yes" "Italian" "0-10" "Yes"]
         ["No" "Yes" "No" "No" "None" "$" "Yes" "No" "Burger" "0-10" "No"]
         ["No" "No" "No" "Yes" "Some" "$$" "Yes" "Yes" "Thai" "0-10" "Yes"]
         ["No" "Yes" "Yes" "No" "Full" "$" "Yes" "No" "Burger" ">60" "No"]
         ["Yes" "Yes" "Yes" "Yes" "Full" "$$$" "No" "Yes" "Italian" "10-30" "No"]
         ["No" "No" "No" "No" "None" "$" "No" "No" "Thai" "0-10" "No"]
         ["Yes" "Yes" "Yes" "Yes" "Full" "$" "No" "No" "Burger" "30-60" "Yes"]]))

(def ^:private values (into {} (for [a attrs] [a (vec (distinct (map a rows)))])))

(deftest decision-tree-test
  (testing "Patrons has the greatest information gain, 0.541 bits"
    (is (= :Patrons (apply max-key #(tr/information-gain % rows :WillWait) attrs)))
    (is (< (abs (- 0.541 (tr/information-gain :Patrons rows :WillWait))) 5e-4)))
  (let [tree (tr/decision-tree-learning rows attrs :WillWait values)]
    (testing "the tree tests Patrons, then Hungry, Type and Fri/Sat, and fits every example"
      (is (= :Patrons (:attribute tree)))
      (is (= :Hungry (get-in tree [:branches "Full" :attribute])))
      (is (= :Type (get-in tree [:branches "Full" :branches "Yes" :attribute])))
      (is (= :FriSat (get-in tree [:branches "Full" :branches "Yes" :branches "Thai" :attribute])))
      (is (zero? (tr/error-rate #(tr/classify tree %) rows :WillWait))))))

(deftest model-selection-test
  (testing "cross-validation over tree depth picks a depth, and the depth-limited trees err more on training data the shallower they are"
    (let [learner (fn [size exs] (let [t (tr/decision-tree-learning exs attrs :WillWait values {:max-depth size})] #(tr/classify t %)))
          {:keys [size errors]} (tr/model-selection learner rows 4 :WillWait 4)]
      (is (<= 1 size 4))
      (is (apply >= (map first (vals errors)))))))

(deftest decision-list-test
  (testing "a decision list of single tests fits the examples"
    (let [dl (tr/decision-list-learning rows attrs :WillWait values 1)]
      (is (some? dl))
      (is (every? #(= (:WillWait %) (tr/decision-list-classify dl "No" %)) rows)))))

(deftest ensembles-test
  (testing "AdaBoost on decision stumps: the training error of the weighted majority reaches zero"
    (let [states (tr/adaboost-seq rows attrs :WillWait values)
          err (fn [st] (tr/error-rate #(tr/weighted-majority st %) rows :WillWait))]
      (is (pos? (err (nth states 1))) "one stump is not enough")
      (is (zero? (err (nth states 10))))))
  (testing "a random forest fits the examples it was grown on"
    (let [f (tr/random-forest rows attrs :WillWait values 25 (allgo.random/rng 3))]
      (is (< (tr/error-rate f rows :WillWait) 0.2))))
  (testing "naive Bayes: posteriors sum to one, and most training examples come out right"
    (let [nb (mo/naive-bayes rows attrs :WillWait)]
      (is (every? #(< (abs (- 1.0 (reduce + (vals (:posterior (nb %)))))) 1e-12) rows))
      (is (>= (count (filter #(= (:WillWait %) (:class (nb %))) rows)) 9)))))
