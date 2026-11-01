(ns allgo.agents-test
  (:require [allgo.agents :as agents]
            [clojure.test :refer [deftest is testing]]))

;; The vacuum world of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 2, scored as aima-python (MIT)
;; scores it.

(def ^:private worlds
  (for [loc [:A :B] a [:clean :dirty] b [:clean :dirty]]
    (agents/vacuum-world loc {:A a :B b})))

(defn- final [agent world steps] (:world (last (agents/run-seq agent world steps))))

(deftest vacuum-test
  (testing "the reflex agent cleans both squares from every start"
    (doseq [w worlds]
      (is (= {:A :clean :B :clean} (:status (final agents/reflex-vacuum-agent w 4))))))
  (testing "the table-driven agent does what the reflex agent does, as far as its table reaches"
    (doseq [w worlds]
      (is (= (:status (final agents/reflex-vacuum-agent w 3))
             (:status (final (agents/table-driven-agent agents/vacuum-table) w 3))))))
  (testing "a simple reflex agent of rules is the reflex agent again"
    (let [rules [[(fn [[_ s]] (= s :dirty)) :suck] [(fn [[l _]] (= l :A)) :right] [(constantly true) :left]]
          agent (agents/simple-reflex-agent identity rules)]
      (doseq [w worlds]
        (is (= (map :action (agents/run-seq agents/reflex-vacuum-agent w 6))
               (map :action (agents/run-seq agent w 6)))))))
  (testing "the model-based agent stops once it knows both are clean, and so scores better than shuttling"
    (doseq [w worlds]
      (let [model (final agents/model-based-vacuum-agent w 10)
            reflex (final agents/reflex-vacuum-agent w 10)]
        (is (= {:A :clean :B :clean} (:status model)))
        (is (> (:score model) (:score reflex))))))
  (testing "scores: +10 for dirt sucked up, -1 a move -- suck, right, suck: 10 - 1 + 10"
    (is (= 19 (:score (final agents/reflex-vacuum-agent (agents/vacuum-world :A {:A :dirty :B :dirty}) 3))))))
