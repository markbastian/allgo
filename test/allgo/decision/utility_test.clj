(ns allgo.decision.utility-test
  (:require [allgo.decision.utility :as ut]
            [clojure.test :refer [deftest is testing]]))

;; The value of information, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 16: n blocks,
;; one holding oil worth C, each sold for C/n, and a survey of block 0.

(defn- oil [n C]
  (let [states (range n)]
    {:prior (zipmap states (repeat (/ 1.0 n)))
     :actions (cons :none states)
     :utility (fn [a s] (cond (= a :none) 0.0 (= a s) (- C (/ C n)) :else (- (/ C n))))}))

(def ^:private survey
  {:name :survey :cost 0.1 :values [true false] :likelihood (fn [e s] (if (= e (= s 0)) 1.0 0.0))})

(deftest vpi-test
  (testing "the survey is worth C/n: worked out, (1/n)(C - C/n) + ((n-1)/n)(C/(n-1) - C/n)"
    (doseq [n [2 3 4 10] C [1.0 7.0]]
      (is (< (abs (- (/ C n) (ut/vpi (oil n C) (:prior (oil n C)) survey))) 1e-12))))
  (testing "evidence that says nothing is worth nothing"
    (is (< (abs (ut/vpi (oil 4 1.0) (:prior (oil 4 1.0)) {:values [true false] :likelihood (constantly 0.5)})) 1e-12))))

(deftest information-gathering-test
  (testing "the agent buys a survey worth more than it costs, then acts on what it says"
    (let [[agent m] (ut/information-gathering-agent (oil 4 1.0) [survey])
          [a1 m1] (agent m nil)
          [a2] (agent m1 [:survey true])]
      (is (= [:request :survey] a1))
      (is (= [:act 0] a2))))
  (testing "but not one that costs more than it is worth"
    (let [[agent m] (ut/information-gathering-agent (oil 4 1.0) [(assoc survey :cost 0.5)])]
      (is (= :act (first (first (agent m nil))))))))
