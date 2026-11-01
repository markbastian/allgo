(ns allgo.logic.tms-test
  (:require [allgo.logic.tms :as tms]
            [clojure.test :refer [deftest is testing]]))

;; The cases of Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 3rd edition, section 12.6.2 (4th edition, chapter 10)

(def ^:private p=>q {:antecedents #{:P :P=>Q} :consequent :Q})
(def ^:private p|r=>q-by-p {:antecedents #{:P :PvR=>Q} :consequent :Q})
(def ^:private p|r=>q-by-r {:antecedents #{:R :PvR=>Q} :consequent :Q})

(defn- told [& justifications] (tms/enable (tms/jtms justifications) :P=>Q :PvR=>Q :R :P))

(deftest jtms
  (testing "TELL(P) with P => Q in the KB brings Q in"
    (is (tms/in? (told p=>q) :Q)))
  (testing "Q's only justification {P, P => Q}: retracting P takes Q out"
    (is (not (tms/in? (tms/retract (told p=>q) :P) :Q))))
  (testing "with {P, P | R => Q} as well it still goes"
    (is (not (tms/in? (tms/retract (told p=>q p|r=>q-by-p) :P) :Q))))
  (testing "with {R, P | R => Q} as well it is spared"
    (is (tms/in? (tms/retract (told p=>q p|r=>q-by-p p|r=>q-by-r) :P) :Q)))
  (testing "a sentence out is kept, and comes back in when its justification does"
    (let [t (tms/retract (told p=>q) :P)]
      (is (= 1 (count (:justifications t))))
      (is (tms/in? (tms/enable t :P) :Q))))
  (testing "two sentences that justify each other do not hold each other up"
    (let [t (tms/enable (tms/jtms [{:antecedents #{:P} :consequent :A}
                                   {:antecedents #{:A} :consequent :B}
                                   {:antecedents #{:B} :consequent :A}]) :P)]
      (is (= #{:P :A :B} (:in t)))
      (is (= #{} (:in (tms/retract t :P))))))
  (testing "rules from definite clauses, each its own node"
    (let [t (tms/enable (tms/jtms (tms/rules "A & B => C" "C => D")) "A & B => C" "C => D" :A :B)]
      (is (tms/in? t :D))
      (is (not (tms/in? (tms/retract t "C => D") :D)))))
  (testing "retraction a step at a time: marked out, then supported back"
    (let [phases (map :phase (tms/retract-seq (told p=>q p|r=>q-by-p p|r=>q-by-r) :P))]
      (is (= :out (first phases)))
      (is (some #{:support} phases)))))

(def ^:private car
  "Why won't the car start? The chapter's example, its knowledge of how
  cars work our own: a dead battery or a broken starter means no
  cranking, no cranking or no fuel means no start, and a dead battery
  can't light the lights."
  (tms/atms [:BatteryDead :StarterBroken :NoFuel] [:LightsOn]
            (concat (tms/horn "BatteryDead => NoCrank" "StarterBroken => NoCrank"
                              "NoCrank => WontStart" "NoFuel => WontStart")
                    [{:antecedents #{:LightsOn :BatteryDead} :consequent :contradiction}])))

(deftest atms
  (testing "the label of WontStart is its minimal explanations"
    (let [end (last (tms/atms-seq car))]
      (is (= #{#{:StarterBroken} #{:NoFuel}} ((:labels end) :WontStart)))
      (is (= #{#{:BatteryDead}} (:nogoods end)) "the lights rule the battery out")
      (is (= #{#{:StarterBroken}} ((:labels end) :NoCrank)))))
  (testing "every hypothetical world at once"
    (let [end (last (tms/atms-seq car))]
      (is (tms/holds-under? end :WontStart #{:NoFuel}))
      (is (tms/holds-under? end :WontStart #{:NoFuel :StarterBroken}))
      (is (not (tms/holds-under? end :WontStart #{})))
      (is (not (tms/holds-under? end :WontStart #{:BatteryDead :NoFuel})) "a nogood")))
  (testing "minimal environments: a sentence explained by A alone is not also by A and B"
    (let [a (tms/atms [:A :B] [] (tms/horn "A => X" "A & B => X"))]
      (is (= #{#{:A}} ((tms/labels a) :X))))))
