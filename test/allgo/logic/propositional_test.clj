(ns allgo.logic.propositional-test
  (:require [allgo.logic.propositional :as pl]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; ---------------------------------------------------------------------------
;; The worked problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 7

(def ^:private wumpus-kb
  "The chapter's knowledge base after the agent has felt nothing in [1,1]
  and a breeze in [2,1]: R1 to R5."
  (into [:and] (map pl/parse ["~P1,1"
                              "B1,1 <=> (P1,2 | P2,1)"
                              "B2,1 <=> (P1,1 | P2,2 | P3,1)"
                              "~B1,1"
                              "B2,1"])))

(def ^:private horn-kb
  "The chapter's definite clauses for forward and backward chaining."
  (map pl/definite ["P => Q" "L & M => P" "B & L => M" "A & P => L" "A & B => L" "A" "B"]))

(deftest wumpus-kb-test
  (testing "no pit in [1,2], by model checking, resolution and DPLL alike"
    (let [q (pl/parse "~P1,2")]
      (is (pl/tt-entails? wumpus-kb q))
      (is (pl/pl-resolution wumpus-kb q))
      (is (pl/entails? wumpus-kb q))))
  (testing "but about [2,2] nothing: neither a pit nor none is entailed"
    (doseq [q ["P2,2" "~P2,2"]]
      (is (not (pl/tt-entails? wumpus-kb (pl/parse q))))
      (is (not (pl/pl-resolution wumpus-kb (pl/parse q))))
      (is (not (pl/entails? wumpus-kb (pl/parse q))))))
  (testing "model checking goes through all 2^7 = 128 models of the seven symbols to prove it"
    (is (= 128 (:checked (last (pl/tt-entails-seq wumpus-kb (pl/parse "~P1,2"))))))))

(deftest horn-test
  (testing "forward chaining proves Q, inferring A B L M P Q in that order"
    (is (pl/pl-fc-entails? horn-kb :Q))
    (is (= [:A :B :L :M :P :Q] (vec (keep :popped (pl/pl-fc-entails-seq horn-kb :Q))))))
  (testing "backward chaining proves it too, and both fail on a symbol nothing concludes"
    (is (pl/pl-bc-entails? horn-kb :Q))
    (is (not (pl/pl-fc-entails? horn-kb :Z)))
    (is (not (pl/pl-bc-entails? horn-kb :Z))))
  (testing "backward chaining does not loop on P => Q, Q => P"
    (is (not (pl/pl-bc-entails? (map pl/definite ["P => Q" "Q => P"]) :Q)))))

;; ---------------------------------------------------------------------------
;; Sentences and clauses

(deftest sentences-test
  (testing "parse and show round-trip, with the precedence of not, and, or, implies, iff"
    (is (= [:iff :A [:or [:and :B [:not :C]] :D]] (pl/parse "A <=> B & ~C | D")))
    (is (= [:implies :A [:implies :B :C]] (pl/parse "A => B => C")))
    (doseq [s ["~(A | B) => C" "A <=> (B & C)" "P1,2 | ~P2,1"]]
      (is (= (pl/parse s) (pl/parse (pl/show (pl/parse s)))))))
  (testing "CNF keeps a sentence's models, on every model of its symbols"
    (doseq [text ["A <=> (B | C)" "~(A & (B => C))" "(A => B) & (B => C) => (A => C)" "A & ~A"]
            :let [s (pl/parse text) cs (pl/cnf s) syms (vec (pl/symbols s))]
            k (range (bit-shift-left 1 (count syms)))
            :let [m (zipmap syms (map #(odd? (bit-shift-right k %)) (range (count syms))))]]
      (is (= (boolean (pl/pl-true? s m)) (every? #(true? (pl/clause-true? % m)) cs)) text)))
  (testing "a tautology has no clauses left; a contradiction has an unsatisfiable one"
    (is (empty? (pl/cnf (pl/parse "(A => B) & (B => C) => (A => C)"))))
    (is (not (pl/satisfiable? (pl/parse "A & ~A"))))))

;; ---------------------------------------------------------------------------
;; DPLL and WalkSAT

(defn- brute-satisfiable? [clauses]
  (let [syms (vec (distinct (map pl/symbol-of (mapcat identity clauses))))
        n (count syms)]
    (boolean (some (fn [k] (let [m (zipmap syms (map #(odd? (bit-shift-right k %)) (range n)))]
                             (every? #(true? (pl/clause-true? % m)) clauses)))
                   (range (bit-shift-left 1 n))))))

(deftest sat-test
  (let [rng (random/rng 9)]
    (testing "DPLL agrees with trying every model, and its models satisfy"
      (doseq [_ (range 200)
              :let [cs (pl/random-3-sat rng 9 (+ 10 (random/below rng 40)))
                    m (pl/dpll cs)]]
        (is (= (brute-satisfiable? cs) (some? m)))
        (when m (is (every? #(true? (pl/clause-true? % m)) cs)))))
    (testing "WalkSAT's models satisfy, and it finds one for easy problems; DPLL proves the hard ones unsatisfiable"
      (doseq [_ (range 20)
              :let [cs (pl/random-3-sat rng 30 60)
                    m (pl/walksat cs {:rng rng})]]
        (is (some? m))
        (is (every? #(true? (pl/clause-true? % m)) cs)))
      (is (nil? (pl/dpll (pl/random-3-sat rng 20 200)))))))
