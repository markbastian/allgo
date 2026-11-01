(ns allgo.logic.description-test
  (:require [allgo.logic.description :as dl]
            [clojure.test :refer [deftest is testing]]))

;; The worked examples of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 3rd edition, section 12.5.2 (4th edition, chapter 10)

(def ^:private men-with-sons
  "The chapter's men with at least three sons who are all unemployed and
  married to doctors, and at most two daughters who are all professors in
  physics or math departments."
  [:and :Man [:at-least 3 :Son] [:at-most 2 :Daughter]
   [:all :Son [:and :Unemployed :Married [:all :Spouse :Doctor]]]
   [:all :Daughter [:and :Professor [:fills :Department :Physics :Math]]]])

(def ^:private tbox
  {:Bachelor [:and :Unmarried :Adult :Male]
   :Man [:and :Adult :Male]
   :Father [:and :Man [:at-least 1 :Son]]
   :Busy men-with-sons})

(deftest subsumption
  (testing "Bachelor = And(Unmarried, Adult, Male)"
    (is (dl/subsumes? tbox :Man :Bachelor))
    (is (not (dl/subsumes? tbox :Bachelor :Man)))
    (is (dl/equivalent? tbox :Bachelor [:and :Male [:and :Adult :Unmarried]])))
  (testing "what follows from the men with sons"
    (is (dl/subsumes? tbox :Father :Busy))
    (is (dl/subsumes? tbox [:at-least 2 :Son] :Busy))
    (is (dl/subsumes? tbox [:all :Son :Married] :Busy))
    (is (dl/subsumes? tbox [:all :Son [:all :Spouse :Doctor]] :Busy))
    (is (dl/subsumes? tbox [:all :Daughter [:fills :Department :Math]] :Busy))
    (is (dl/subsumes? tbox [:at-most 5 :Daughter] :Busy)))
  (testing "and what does not"
    (is (not (dl/subsumes? tbox [:at-least 1 :Daughter] :Busy)))
    (is (not (dl/subsumes? tbox [:at-most 1 :Daughter] :Busy)))
    (is (not (dl/subsumes? tbox [:all :Son :Doctor] :Busy)))
    (is (not (dl/subsumes? tbox :Busy :Father))))
  (testing "fills and one-of"
    (is (dl/subsumes? tbox [:at-least 2 :Department] [:fills :Department :Physics :Math]))
    (is (dl/subsumes? tbox [:one-of :Physics :Math :Biology] [:one-of :Math :Physics]))
    (is (not (dl/subsumes? tbox [:one-of :Math] [:one-of :Math :Physics]))))
  (testing "an incoherent concept is subsumed by everything"
    (let [none [:and [:at-least 3 :Son] [:at-most 2 :Son]]]
      (is (dl/incoherent? (dl/normalize tbox none)))
      (is (dl/subsumes? tbox :Bachelor none))
      (is (dl/incoherent? (dl/normalize tbox [:and [:fills :Department :Physics] [:all :Department [:one-of :Math]]])))
      (is (not (dl/incoherent? (dl/normalize tbox men-with-sons)))))))

(deftest classification
  (testing "the names placed in a taxonomy, each under its most specific subsumers"
    (is (= {:Adult #{} :Male #{} :Man #{:Adult :Male} :Bachelor #{:Man} :Father #{:Man} :Busy #{:Father}}
           (dl/taxonomy tbox :Adult :Male))))
  (testing "an individual is a Bachelor from what it is told; fillers taken as all it has"
    (let [abox {:Bob {:prims [:Adult :Male :Unmarried]}
                :Al {:prims [:Adult :Male] :fills {:Son [:S1 :S2 :S3]}}
                :S1 {:prims [:Unemployed :Married]} :S2 {:prims [:Unemployed :Married]} :S3 {:prims [:Married]}}]
      (is (dl/instance-of? tbox abox :Bachelor :Bob))
      (is (dl/instance-of? tbox abox :Father :Al))
      (is (not (dl/instance-of? tbox abox [:all :Son :Unemployed] :Al)))
      (is (dl/instance-of? tbox abox [:all :Son :Married] :Al)))))
