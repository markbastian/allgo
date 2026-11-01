(ns allgo.logic.categories-test
  (:require [allgo.logic.categories :as cat]
            [allgo.logic.first-order :as fo]
            [clojure.test :refer [deftest is testing]]))

;; The worked examples of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 3rd edition, sections 12.2 and 12.5 (4th edition,
;; chapter 10)

(def ^:private figure-12-5
  "Four objects -- John, Mary, 1 and 2 -- and four categories."
  {:member {:Mary #{:FemalePersons} :John #{:MalePersons}}
   :subset {:FemalePersons #{:Persons} :MalePersons #{:Persons} :Persons #{:Mammals}}
   :links {:Mary {:SisterOf :John} :John {:Legs 1}}
   :boxed {:Persons {:Legs 2}}
   :double-boxed {:Persons {:HasMother :FemalePersons}}})

(deftest inheritance
  (testing "Mary has two legs, by being a person"
    (is (= 2 (cat/inherit figure-12-5 :Mary :Legs)))
    (is (= {:Persons 2} (:values (last (cat/inherit-seq figure-12-5 :Mary :Legs))))))
  (testing "John's own one leg overrides the default"
    (is (= 1 (cat/inherit figure-12-5 :John :Legs))))
  (testing "or a OneLeggedPersons category below Persons does"
    (let [net (-> figure-12-5
                  (assoc-in [:member :John] #{:MalePersons :OneLeggedPersons})
                  (assoc-in [:subset :OneLeggedPersons] #{:Persons})
                  (update :links dissoc :John)
                  (assoc-in [:boxed :OneLeggedPersons :Legs] 1))]
      (is (= 1 (cat/inherit net :John :Legs)))))
  (testing "a person's mother is a female person: a link of the members, not of the category"
    (is (= :FemalePersons (cat/value-category figure-12-5 :Mary :HasMother)))
    (is (nil? (cat/inherit figure-12-5 :Persons :HasMother))))
  (testing "inheritance walks up a node at a time"
    (is (= [:Mary :FemalePersons :Persons :Mammals] (map :node (cat/inherit-seq figure-12-5 :Mary :Legs)))))
  (testing "membership through subsets"
    (is (cat/member? figure-12-5 :Mary :Mammals))
    (is (not (cat/member? figure-12-5 :Mary :MalePersons)))))

(deftest the-nixon-diamond
  (testing "multiple inheritance with values that disagree is ambiguous"
    (let [net {:member {:Nixon #{:Quakers :Republicans}}
               :boxed {:Quakers {:Pacifist true} :Republicans {:Pacifist false}}}]
      (is (= {:ambiguous {:Quakers true :Republicans false}} (cat/inherit net :Nixon :Pacifist)))))
  (testing "a more specific category's value wins even when it is further away"
    (let [net {:member {:Opus #{:Penguins :Pets}}
               :subset {:Penguins #{:Seabirds} :Seabirds #{:Birds} :Pets #{:Animals} :Birds #{:Animals}}
               :boxed {:Birds {:Flies true} :Penguins {:Flies false}}}]
      (is (false? (cat/inherit net :Opus :Flies))))))

(deftest reification
  (testing "Fly(Shankar, NewYork, NewDelhi, Yesterday) as figure 12.6's Fly17"
    (let [net (cat/reify-event :Fly17 :FlyEvents [:Agent :Origin :Destination :During]
                               [:Shankar :NewYork :NewDelhi :Yesterday])]
      (is (= :NewDelhi (cat/inherit net :Fly17 :Destination)))
      (is (cat/member? net :Fly17 :FlyEvents)))))

(deftest decompositions
  (let [males #{:Al :Bo} females #{:Cy :Di} animals #{:Al :Bo :Cy :Di}
        americans #{:Al :Bo} canadians #{:Bo :Cy} mexicans #{:Di}]
    (testing "Partition({Males, Females}, Animals)"
      (is (cat/partition? [males females] animals)))
    (testing "an exhaustive decomposition with a dual citizen is no partition"
      (is (cat/exhaustive-decomposition? [americans canadians mexicans] animals))
      (is (not (cat/partition? [americans canadians mexicans] animals))))
    (testing "Disjoint({Animals, Vegetables})"
      (is (cat/disjoint? [animals #{:Kale}])))))

(deftest part-of
  (testing "PartOf is transitive: PartOf(Bucharest, Earth), by allgo.logic.first-order"
    (let [kb (map fo/definite ["PartOf(Bucharest, Romania)" "PartOf(Romania, EasternEurope)"
                               "PartOf(EasternEurope, Europe)" "PartOf(Europe, Earth)"
                               "PartOf(x, y) & PartOf(y, z) => PartOf(x, z)"])]
      (is (some #{{'x :Earth}} (fo/fol-fc-ask kb (fo/parse "PartOf(Bucharest, x)")))))))
