(ns allgo.logic.defaults-test
  (:require [allgo.logic.defaults :as df]
            [clojure.test :refer [deftest is testing]]))

;; The worked examples of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 3rd edition, section 12.6.1 (4th edition, chapter 10)

(def ^:private nixon
  ["Republican(Nixon) & Quaker(Nixon)"
   "Republican(x) & ~Abnormal2(x) => ~Pacifist(x)"
   "Quaker(x) & ~Abnormal3(x) => Pacifist(x)"])

(defn- k [s] (keyword s))

(defn- true-atoms [m] (set (keep (fn [[k v]] (when v k)) m)))

(deftest circumscription
  (testing "birds fly: Flies(Tweety) from Bird(Tweety), but not once Abnormal1(Tweety) is told"
    (let [kb ["Bird(Tweety)" "Bird(x) & ~Abnormal1(x) => Flies(x)"]]
      (is (df/circumscription-entails? kb #{:Abnormal1} "Flies(Tweety)"))
      (is (not (df/circumscription-entails? (conj kb "Abnormal1(Tweety)") #{:Abnormal1} "Flies(Tweety)")))))
  (testing "the Nixon diamond has two preferred models, and the reasoner stays agnostic"
    (is (= #{#{(k "Republican(Nixon)") (k "Quaker(Nixon)") (k "Abnormal2(Nixon)") (k "Pacifist(Nixon)")}
             #{(k "Republican(Nixon)") (k "Quaker(Nixon)") (k "Abnormal3(Nixon)")}}
           (set (map true-atoms (df/preferred-models nixon #{:Abnormal2 :Abnormal3})))))
    (is (not (df/circumscription-entails? nixon #{:Abnormal2 :Abnormal3} "Pacifist(Nixon)")))
    (is (not (df/circumscription-entails? nixon #{:Abnormal2 :Abnormal3} "~Pacifist(Nixon)"))))
  (testing "prioritized: religion over politics -- Abnormal3 minimized first -- makes him a pacifist"
    (is (df/circumscription-entails? nixon [#{:Abnormal3} #{:Abnormal2}] "Pacifist(Nixon)"))
    (is (df/circumscription-entails? nixon [#{:Abnormal2} #{:Abnormal3}] "~Pacifist(Nixon)")))
  (testing "every model, a step at a time"
    (is (= 4 (count (df/models-seq (df/propositionalize nixon)))))))

(def ^:private nixon-defaults
  ["Republican(x) : ~Pacifist(x) / ~Pacifist(x)"
   "Quaker(x) : Pacifist(x) / Pacifist(x)"])

(deftest default-logic
  (testing "a default rule P : J1, ..., Jn / C"
    (is (= {:prerequisite [:Bird 'x] :justifications [[:Flies 'x]] :conclusion [:Flies 'x]}
           (df/default-rule "Bird(x) : Flies(x) / Flies(x)"))))
  (testing "the Nixon diamond: one fact, two defaults, two extensions"
    (let [exts (df/extensions ["Republican(Nixon) & Quaker(Nixon)"] nixon-defaults)]
      (is (= #{#{(k "Pacifist(Nixon)")} #{[:not (k "Pacifist(Nixon)")]}} (set (map :conclusions exts))))
      (is (= #{true false} (set (map #(df/believes? % "Pacifist(Nixon)") exts))))))
  (testing "Tweety flies by default; a penguin is an exception that stops it"
    (let [exts (df/extensions ["Bird(Tweety)" "Penguin(Opus)" "Penguin(x) => Bird(x) & ~Flies(x)"]
                              ["Bird(x) : Flies(x) / Flies(x)"])]
      (is (= 1 (count exts)))
      (is (df/believes? (first exts) "Flies(Tweety)"))
      (is (df/believes? (first exts) "~Flies(Opus)"))))
  (testing "nonmonotonic: told more, it believes less"
    (is (df/believes? (first (df/extensions ["Bird(Tweety)"] ["Bird(x) : Flies(x) / Flies(x)"])) "Flies(Tweety)"))
    (is (not (df/believes? (first (df/extensions ["Bird(Tweety)" "~Flies(Tweety)"] ["Bird(x) : Flies(x) / Flies(x)"]))
                           "Flies(Tweety)"))))
  (testing "a process that applies a rule whose justification a later conclusion denies fails"
    (let [states (df/extensions-seq ["P"] ["P : A / B" "P : ~B / C"])]
      (is (some :failed? states) "C first, then B denies ~B")
      (is (= [#{:B}] (map :conclusions (:extensions (last states)))) "B first, and ~B is no longer consistent")))
  (testing "the process tree: the root, then each order of the two rules"
    (is (= 3 (count (df/extensions-seq ["Republican(Nixon) & Quaker(Nixon)"] nixon-defaults))))))
