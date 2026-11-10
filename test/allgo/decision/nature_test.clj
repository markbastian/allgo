(ns allgo.decision.nature-test
  (:require [allgo.decision.nature :as nature]
            [clojure.test :refer [deftest is testing]]))

;; Games against nature, after LaValle, *Planning Algorithms*, Cambridge
;; University Press, 2006, section 9.2, with Examples 9.8 to 9.10's cost
;; matrix.

(def ^:private costs
  {:u1 {:t1 1 :t2 -1 :t3 0}
   :u2 {:t1 -1 :t2 2 :t3 -2}
   :u3 {:t1 2 :t2 -1 :t3 1}})

(deftest nature-test
  (testing "Example 9.9, nature unknown: u1, whose worst case is 1"
    (is (= [:u1 1] (nature/worst-case costs))))
  (testing "nature's distribution 1/5, 1/5, 3/5: expected costs 0, -1 and 4/5, so u2"
    (let [[a c] (nature/expected costs {:t1 0.2 :t2 0.2 :t3 0.6})]
      (is (= :u2 a))
      (is (< (abs (+ 1.0 c)) 1e-12))))
  (testing "1/10, 4/5, 1/10: u1"
    (is (= :u1 (first (nature/expected costs {:t1 0.1 :t2 0.8 :t3 0.1})))))
  (testing "Example 9.10's regret matrix, and the least worst-case regret: u1, 2"
    (is (= {:u1 {:t1 2 :t2 0 :t3 2} :u2 {:t1 0 :t2 3 :t3 0} :u3 {:t1 3 :t2 0 :t3 3}} (nature/regrets costs)))
    (is (= [:u1 2] (nature/minimax-regret costs))))
  (testing "the least expected regret with nature uniform: the example says u1, at 1, but u1's is 4/3 -- it is u2's that is 1"
    (let [[a c] (nature/expected (nature/regrets costs) {:t1 (/ 1.0 3) :t2 (/ 1.0 3) :t3 (/ 1.0 3)})]
      (is (= :u2 a))
      (is (< (abs (- 1.0 c)) 1e-12))))
  (testing "Hurwicz: fully pessimistic is minimax; fully optimistic takes the best best case, u2's -2"
    (is (= :u1 (first (nature/hurwicz costs 1.0))))
    (is (= [:u2 -2.0] (nature/hurwicz costs 0.0))))
  (testing "an observation that rules out t2 leaves u2 best in the worst case"
    (is (= :u2 (first (nature/worst-case (nature/observe costs [:t1 :t3])))))))
