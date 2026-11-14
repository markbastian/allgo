(ns allgo.logic.temporal-test
  "Linear temporal logic, after Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 10.4.1: the four
  operators on behaviors, the identities that tie them together, and a
  mutual-exclusion protocol checked for safety and liveness."
  (:require [allgo.logic.temporal :as ltl]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(deftest the-four-operators
  (let [;; p, then p, then q, then p and q forever after
        trace {:prefix [#{:p} #{:p} #{:q}] :loop [#{:p :q}]}]
    (testing "X P: P in the next step"
      (is (ltl/holds? trace [:next :p]))
      (is (not (ltl/holds? trace [:next [:next :p]]))))
    (testing "F P: P eventually"
      (is (ltl/holds? trace [:eventually :q]))
      (is (not (ltl/holds? trace [:eventually :r]))))
    (testing "G P: P always -- not p, which fails at the third step, but eventually always p"
      (is (not (ltl/holds? trace [:always :p])))
      (is (ltl/holds? trace [:eventually [:always :p]])))
    (testing "P U Q: P remains true until Q occurs"
      (is (ltl/holds? trace [:until :p :q]))
      (is (not (ltl/holds? {:prefix [#{:p} #{}] :loop [#{:q}]} [:until :p :q]))))
    (testing "and a loop never reaching Q: P U Q fails though P holds forever, P R Q holds"
      (is (not (ltl/holds? {:prefix [] :loop [#{:p}]} [:until :p :q])))
      (is (ltl/holds? {:prefix [] :loop [#{:q}]} [:release :p :q])))))

(defn- random-lasso [rng]
  (let [state #(set (filter (fn [_] (< (rng) 0.5)) [:p :q]))]
    {:prefix (vec (repeatedly (random/below rng 4) state))
     :loop (vec (repeatedly (inc (random/below rng 3)) state))}))

(deftest identities
  (testing "the identities the operators obey, on 400 random behaviors"
    (let [rng (random/rng 10)
          same (fn [l a b] (= (ltl/holds? l a) (ltl/holds? l b)))]
      (doseq [l (repeatedly 400 #(random-lasso rng))]
        (is (same l [:eventually :p] [:until [:or :p [:not :p]] :p]) "F P is true U P")
        (is (same l [:always :p] [:not [:eventually [:not :p]]]) "G P is not F not P")
        (is (same l [:not [:until :p :q]] [:release [:not :p] [:not :q]]) "not (P U Q) is (not P) R (not Q)")
        (is (same l [:until :p :q] [:or :q [:and :p [:next [:until :p :q]]]]) "P U Q unfolds one step")
        (is (same l [:always :p] [:and :p [:next [:always :p]]]) "G P unfolds one step")
        (is (same l [:eventually [:always :p]] [:not [:always [:eventually [:not :p]]]]) "F G P is not G F not P")))))

;; Two processes, each idle, trying or critical, and either may step, or
;; neither. The naive protocol checks nothing; the turn-taking one lets a
;; process into its critical section only on its turn, passing the turn
;; when it leaves.
(defn- mutex [checked?]
  {:initial [{:a :idle :b :idle :turn :a}]
   :successors (fn [{:keys [turn] :as s}]
                 ;; either process steps, or nothing happens: a process
                 ;; may stay idle as long as it likes
                 (cons s (for [p [:a :b]
                               :let [other (if (= p :a) :b :a)
                                     next-state (case (s p) :idle :trying :trying :critical :critical :idle)]
                               :when (or (not checked?) (not= next-state :critical) (= turn p))]
                           (cond-> (assoc s p next-state)
                             (= (s p) :critical) (assoc :turn other)))))
   :atoms (fn [s] (set (for [p [:a :b] :let [st (s p)]] (keyword (str (name p) "-" (name st))))))})

(deftest mutual-exclusion
  (let [safety [:always [:not [:and :a-critical :b-critical]]]
        liveness [:always [:implies :a-trying [:eventually :a-critical]]]]
    (testing "the naive protocol lets both in at once: a counterexample, which indeed falsifies safety"
      (let [ce (ltl/counterexample (mutex false) safety 12)]
        (is (some? ce))
        (is (not (ltl/holds? ce safety)))))
    (testing "taking turns keeps them apart: no counterexample among every behavior of up to twelve states"
      (is (nil? (ltl/counterexample (mutex true) safety 12))))
    (testing "but taking turns alone does not promise each a turn: b may idle forever with the turn, a trying in vain"
      (let [ce (ltl/counterexample (mutex true) liveness 12)]
        (is (some? ce))
        (is (not (ltl/holds? ce liveness)))))))
