(ns allgo.logic.modal-test
  (:require [allgo.logic.modal :as md]
            [clojure.test :refer [deftest is testing]]))

;; Figure 12.4 of Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 3rd edition (4th edition, chapter 10): I is "Superman's
;; secret identity is Clark Kent", R "the weather report for tomorrow is
;; rain". The arrows are as the text describes them -- the figure's own
;; are not machine-readable -- and every world is accessible from itself.

(defn- both-ways [pairs]
  (reduce (fn [m [a b]] (-> m (update a (fnil conj #{}) b) (update b (fnil conj #{}) a))) {} pairs))

(defn- all-of [ws] (for [x ws y ws] [x y]))

(def ^:private valuations
  {:w0 [:I :R] :w1 [:R] :w2 [:I] :w3 [] :w4 [:I :R] :w5 [:R] :w6 [:I] :w7 []})

(def ^:private top-left
  "Superman knows his identity; neither has seen the report."
  (md/model (select-keys valuations [:w0 :w1 :w2 :w3])
            {:Superman (both-ways [[:w0 :w2] [:w1 :w3]])
             :Lois (both-ways (all-of [:w0 :w1 :w2 :w3]))}))

(def ^:private top-right
  "Lois has seen the report, and everyone knows it."
  (md/model (select-keys valuations [:w4 :w5 :w6 :w7])
            {:Superman (both-ways [[:w4 :w6] [:w5 :w7]])
             :Lois (both-ways [[:w4 :w5] [:w6 :w7]])}))

(def ^:private bottom
  "The two joined: Lois may or may not have seen it, and Superman can't
  tell which."
  (md/model valuations
            {:Superman (both-ways (concat (all-of [:w0 :w2 :w4 :w6]) (all-of [:w1 :w3 :w5 :w7])))
             :Lois (both-ways (concat (all-of [:w0 :w1 :w2 :w3]) [[:w4 :w5] [:w6 :w7]]))}))

(defn- knows-whether [a p] [:or [:knows a p] [:knows a [:not p]]])

(deftest figure-12-4
  (testing "top left: Lois doesn't know whether Clark is Superman, but knows Superman knows"
    (is (md/valid? top-left [:knows :Lois (knows-whether :Superman :I)]))
    (is (empty? (md/worlds-where top-left (knows-whether :Lois :I))))
    (is (= #{:w0 :w2} (md/worlds-where top-left [:knows :Superman :I]))))
  (testing "top right: in w4 Lois knows R, in w6 she knows ~R; Superman knows she knows"
    (is (md/holds? top-right :w4 [:knows :Lois :R]))
    (is (md/holds? top-right :w6 [:knows :Lois [:not :R]]))
    (is (md/valid? top-right [:knows :Superman (knows-whether :Lois :R)]))
    (is (not (md/holds? top-right :w4 [:knows :Superman :R]))))
  (testing "bottom: in w0 Superman knows I but not R, nor whether Lois knows R"
    (is (md/holds? bottom :w0 [:knows :Superman :I]))
    (is (not (md/holds? bottom :w0 (knows-whether :Superman :R))))
    (is (not (md/holds? bottom :w0 (knows-whether :Superman (knows-whether :Lois :R)))))
    (is (md/holds? bottom :w4 [:knows :Lois :R]))))

(deftest bond-and-the-spy
  (testing "K_Bond (exists x Spy(x)) holds and (exists x K_Bond Spy(x)) does not"
    (let [m (md/model {:u [[:Spy :A]] :v [[:Spy :B]]} {:Bond (both-ways [[:u :v]])} :domain [:A :B])]
      (is (md/valid? m '[:knows :Bond [:exists [x] [:Spy x]]]))
      (is (empty? (md/worlds-where m '[:exists [x] [:knows :Bond [:Spy x]]]))))))

(deftest the-axioms
  (let [p-and-q [:and [:knows :a :P] [:knows :a [:implies :P :Q]]]
        refl (md/model {:u [:P :Q] :v [:P] :x [:Q]} {:a {:u #{:v} :v #{:x}}})
        not-refl (md/model {:u [] :v [:P]} {:a {:u #{:v}}} :reflexive? false)]
    (testing "an agent knows what follows from what it knows, in any model"
      (is (md/valid? refl [:implies p-and-q [:knows :a :Q]]))
      (is (md/valid? not-refl [:implies p-and-q [:knows :a :Q]])))
    (testing "K_a P => P needs reflexive accessibility"
      (is (md/reflexive? refl :a))
      (is (md/valid? refl [:implies [:knows :a :P] :P]))
      (is (not (md/reflexive? not-refl :a)))
      (is (not (md/valid? not-refl [:implies [:knows :a :P] :P]))))
    (testing "K_a P => K_a K_a P needs transitive accessibility"
      (let [m (md/model {:u [:P] :v [:P] :x []} {:a {:u #{:v} :v #{:x}}})]
        (is (not (md/transitive? m :a)))
        (is (not (md/holds? m :u [:implies [:knows :a :P] [:knows :a [:knows :a :P]]]))))
      (is (md/valid? top-left [:implies [:knows :Lois :I] [:knows :Lois [:knows :Lois :I]]])))
    (testing "every agent knows P | ~P; not every agent knows P or knows ~P"
      (is (md/valid? top-left [:knows :Lois [:or :R [:not :R]]]))
      (is (not (md/valid? top-left (knows-whether :Lois :R)))))))

(deftest labeling
  (testing "a step a subsentence, innermost first"
    (let [steps (md/check-seq top-left [:knows :Lois (knows-whether :Superman :I)])]
      (is (= :I (:sentence (first steps))))
      (is (= #{:w0 :w1 :w2 :w3} (:worlds (last steps)))))))
