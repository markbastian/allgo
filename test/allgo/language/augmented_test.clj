(ns allgo.language.augmented-test
  (:require [allgo.language.augmented :as aug]
            [allgo.logic.first-order :as fo]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

;; Augmented grammars, after Russell and Norvig, *Artificial Intelligence:
;; A Modern Approach*, 4th edition, section 23.4.

(def ^:private r aug/rule)

(defn- words [s] (str/split s #" "))

(defn- renamed
  "`term` with its variables renamed in order of appearance: equal for two
  terms exactly when one is a renaming of the other."
  [term]
  (fo/subst (into {} (map-indexed (fn [i v] [v (symbol (str "v" i))]) (fo/variables term))) term))

(def ^:private arithmetic aug/arithmetic)

(deftest arithmetic-test
  (testing "figure 23.11: 3 + (4 ÷ 2) is Exp(5), and every node of its tree has the figure's value"
    (let [{:keys [term tree]} (aug/parse (words "3 + ( 4 ÷ 2 )") arithmetic '[:Exp x])]
      (is (= [:Exp 5] term))
      (is (= [[:Exp 5]
              [[:Exp 3] [[:Number 3] [[:Digit 3] "3"]]]
              [[:Operator :+] "+"]
              [[:Exp 2] "(" [[:Exp 2] [[:Exp 4] [[:Number 4] [[:Digit 4] "4"]]]
                             [[:Operator :÷] "÷"]
                             [[:Exp 2] [[:Number 2] [[:Digit 2] "2"]]]] ")"]]
             tree))))
  (testing "Number(10 x1 + x2) -> Number(x1) Digit(x2): numbers of more than one digit"
    (is (= [:Exp 123] (:term (aug/parse (words "1 2 3") arithmetic '[:Exp x])))))
  (testing "the grammar knows no precedence: 1 + 2 × 3 has two parses, 7 and 9"
    (is (= #{[:Exp 7] [:Exp 9]} (set (map :term (aug/parses (words "1 + 2 × 3") arithmetic '[:Exp x]))))))
  (testing "a step a span, shortest first: four words, ten spans"
    (is (= 11 (count (aug/chart-seq (words "3 + 4 5") arithmetic))))))

;; Figure 23.12's grammar, its semantics sentences of first-order logic.

(def ^:private loves aug/ali-loves-bo)

(deftest semantics-test
  (testing "figure 23.12: Ali loves Bo is S(Loves(Ali, Bo)), its VP lambda x Loves(x, Bo), its verb lambda y lambda x Loves(x, y)"
    (let [{:keys [term tree]} (aug/parse (words "Ali loves Bo") loves)
          [_ np [vp-term verb _]] tree]
      (is (= [:S [:Loves :Ali :Bo]] term))
      (is (= [[:NP :Ali] [[:Name :Ali] "Ali"]] np))
      (is (= (renamed '[:VP [:lambda x [:Loves x :Bo]]]) (renamed vp-term)))
      (is (= (renamed '[:Verb [:lambda y [:lambda x [:Loves x y]]]]) (renamed (first verb))))))
  (testing "and the figure's three other sentences"
    (is (= [[:S [:Loves :Bo :Ali]] [:S [:Loves :Ali :Ali]] [:S [:Loves :Bo :Bo]]]
           (map #(:term (aug/parse (words %) loves)) ["Bo loves Ali" "Ali loves Ali" "Bo loves Bo"]))))
  (testing "beta reduction leaves a lambda that binds its own variable alone"
    (is (= '[:lambda x x] (aug/beta-reduce '[[:lambda y [:lambda x x]] :A])))
    (is (= '[:P :A [:lambda x x]] (aug/beta-reduce '[[:lambda x [:P x [:lambda x x]]] :A])))))

;; Agreement, after figure 23.9: case, and person and number. The figure's
;; VP rule gives the object the verb's person and number too; here it has
;; its own, pn2, so that "I see them" parses.

(def ^:private agreement
  {:rules [(r '[:S v] '[:NP :Sbj pn n] '[:VP pn v])
           (r '[:NP c pn n] '[:Pronoun c pn n])
           (r '[:VP pn v] '[:Verb pn v] '[:NP :Obj pn2 n])
           (r [:Pronoun :Sbj :1S "I"] "I") (r [:Pronoun :Sbj :1P "we"] "we") (r [:Pronoun :Sbj :3S "she"] "she")
           (r [:Pronoun :Obj :1S "me"] "me") (r [:Pronoun :Obj :3P "them"] "them") (r [:Pronoun :Obj :3S "her"] "her")
           (r [:Verb :1S "see"] "see") (r [:Verb :1P "see"] "see") (r [:Verb :3S "see"] "sees")
           (r '[:Verb pn "saw"] "saw")]})

(deftest agreement-test
  (testing "the subject in the subjective case, the object in the objective, and subject and verb agreeing"
    (doseq [s ["I see them" "we see her" "she sees me" "I saw her" "she saw them"]]
      (is (aug/parse (words s) agreement) s)))
  (testing "what E0 overgenerates -- \"I saw she\", \"I sees her\" -- and \"Me see them\" do not parse"
    (doseq [s ["I saw she" "I sees her" "me see them" "she see them"]]
      (is (nil? (aug/parse (words s) agreement)) s)))
  (testing "the sentence's head is the verb's: S(v)"
    (is (= [:S "see"] (:term (aug/parse (words "I see them") agreement))))))

;; A lexicalized PCFG: rules' probabilities functions of the phrases' head
;; words. The chapter gives P1(v, n) and P4(n, c, m) no values; these are
;; ours, with its two lexical rules' probabilities.

(def ^:private p1 {["ate" "banana"] 0.1 ["ate" "bandanna"] 0.0001})
(def ^:private p4 {["spaghetti" "and" "meatballs"] 0.1 ["spaghetti" "or" "lasagna"] 0.1 ["meatballs" "or" "lasagna"] 0.01
                   ["pie" "or" "cake"] 0.1 ["spaghetti" "or" "cake"] 0.01})

(def ^:private lexicalized
  {:rules [{:lhs '[:VP v] :rhs '[[:Verb v] [:NP n]] :p (fn [{:syms [v n]}] (get p1 [v n] 0.001))}
           (r '[:NP n] '[:Article a] '[:Noun n])
           (r '[:NP n] '[:Noun n])
           {:lhs '[:NP n] :rhs '[[:NP n] [:Conjunction c] [:NP m]] :p (fn [{:syms [n c m]}] (get p4 [n c m] 0.05))}
           {:lhs [:Verb "ate"] :rhs ["ate"] :p 0.002}
           {:lhs [:Noun "banana"] :rhs ["banana"] :p 0.0007}
           {:lhs [:Noun "bandanna"] :rhs ["bandanna"] :p 0.0007}
           (r [:Article "a"] "a")
           (r [:Conjunction "and"] "and") (r [:Conjunction "or"] "or")
           (r [:Noun "spaghetti"] "spaghetti") (r [:Noun "meatballs"] "meatballs") (r [:Noun "lasagna"] "lasagna")
           (r [:Noun "pie"] "pie") (r [:Noun "cake"] "cake")]})

(defn- bracketing [tree]
  (if (string? tree) tree (let [kids (map bracketing (rest tree))] (if (= 1 (count kids)) (first kids) (vec kids)))))

(deftest lexicalized-test
  (testing "P1(ate, banana) > P1(ate, bandanna), the nouns themselves equally likely: ate a banana is the likelier VP"
    (let [p #(:p (aug/parse (words %) lexicalized '[:VP v]))]
      (is (< (* 100 (p "ate a bandanna")) (p "ate a banana")))))
  (testing "the head words choose the bracketing: [[spaghetti and meatballs] or lasagna], but [spaghetti and [pie or cake]]"
    (is (= [["spaghetti" "and" "meatballs"] "or" "lasagna"]
           (bracketing (:tree (aug/parse (words "spaghetti and meatballs or lasagna") lexicalized '[:NP n])))))
    (is (= ["spaghetti" "and" ["pie" "or" "cake"]]
           (bracketing (:tree (aug/parse (words "spaghetti and pie or cake") lexicalized '[:NP n])))))
    (is (= 2 (count (aug/parses (words "spaghetti and pie or cake") lexicalized '[:NP n]))))))
