(ns allgo.logic.first-order-test
  (:require [allgo.logic.first-order :as fo]
            [clojure.test :refer [deftest is testing]]))

;; ---------------------------------------------------------------------------
;; The worked problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 9, as aima-python (MIT) states
;; them

(def ^:private crime-text
  "The chapter's crime: is Colonel West a criminal? aima-python's
  `crime_kb`."
  ["American(x) & Weapon(y) & Sells(x, y, z) & Hostile(z) => Criminal(x)"
   "Owns(Nono, M1)"
   "Missile(M1)"
   "Missile(x) & Owns(Nono, x) => Sells(West, x, Nono)"
   "Missile(x) => Weapon(x)"
   "Enemy(x, America) => Hostile(x)"
   "American(West)"
   "Enemy(Nono, America)"])

(def ^:private crime (map fo/definite crime-text))

(def ^:private rabbits
  "aima-python's `test_kb`, which it marks as its own, not the book's."
  (map fo/definite ["Farmer(Mac)" "Rabbit(Pete)" "Mother(MrsMac, Mac)" "Mother(MrsRabbit, Pete)"
                    "Rabbit(r) & Farmer(f) => Hates(f, r)" "Mother(m, c) => Loves(m, c)"
                    "Mother(m, r) & Rabbit(r) => Rabbit(m)" "Farmer(f) => Human(f)"
                    "Mother(m, h) & Human(h) => Human(m)"]))

(def ^:private append
  "The chapter's Append, with lists built from Cons and Nil."
  (map fo/definite ["Append(Nil, y, y)" "Append(x, y, z) => Append(Cons(a, x), y, Cons(a, z))"]))

(defn- answers [ask kb q] (set (ask kb (fo/parse q))))

(deftest crime-test
  (testing "forward and backward chaining both find West the criminal, and Nono the enemy"
    (doseq [ask [fo/fol-fc-ask fo/fol-bc-ask]]
      (is (= #{{'x :West}} (answers ask crime "Criminal(x)")))
      (is (= #{{'x :Nono}} (answers ask crime "Enemy(x, America)")))))
  (testing "forward chaining infers it in two rounds: Sells, Weapon and Hostile, then Criminal"
    (let [states (fo/fol-fc-ask-seq crime (fo/parse "Criminal(x)"))
          by-round (group-by :round (filter #(seq (:inferred %)) states))]
      (is (= #{"Sells(West, M1, Nono)" "Weapon(M1)" "Hostile(Nono)"}
             (set (map fo/show (mapcat :inferred (by-round 1))))))
      (is (= ["Criminal(West)"] (map fo/show (mapcat :inferred (by-round 2)))))))
  (testing "resolution refutes ~Criminal(West) from the sentences, quantifiers and all"
    (let [kb (map fo/parse ["forall x, y, z: American(x) & Weapon(y) & Sells(x, y, z) & Hostile(z) => Criminal(x)"
                            "exists x: Owns(Nono, x) & Missile(x)"
                            "forall x: Missile(x) & Owns(Nono, x) => Sells(West, x, Nono)"
                            "forall x: Missile(x) => Weapon(x)"
                            "forall x: Enemy(x, America) => Hostile(x)"
                            "American(West)" "Enemy(Nono, America)"])
          end (last (take 500 (fo/resolution-seq kb (fo/parse "Criminal(West)"))))
          lines (fo/proof-lines (:proof end))]
      (is (true? (:proved? end)))
      (is (empty? (:clause (peek lines))) "the proof ends in the empty clause")
      (is (:goal? (first lines)) "and starts from the negated query")
      (is (not (true? (fo/resolution kb (fo/parse "Criminal(Nono)") 300)))))))

(deftest rabbits-test
  (testing "aima-python's answers, by both chainings"
    (doseq [ask [fo/fol-fc-ask fo/fol-bc-ask]]
      (is (= #{{'x :Mac}} (answers ask rabbits "Farmer(x)")))
      (is (= #{{'x :Mac} {'x :MrsMac}} (answers ask rabbits "Human(x)")))
      (is (= #{{'x :Pete} {'x :MrsRabbit}} (answers ask rabbits "Rabbit(x)")))))
  (testing "the conjuncts in the other order send backward chaining down forever; the depth limit stops it, answers found"
    (let [loopy (concat (butlast rabbits) [(fo/definite "Human(h) & Mother(m, h) => Human(m)")])]
      (is (contains? (set (take 2 (fo/fol-bc-ask loopy (fo/parse "Human(x)") {:max-depth 12}))) {'x :Mac})))))

(deftest append-test
  (testing "Append(x, y, [1, 2]) has three answers, in the order depth-first search finds them"
    (is (= [["Nil" "Cons(1, Cons(2, Nil))"] ["Cons(1, Nil)" "Cons(2, Nil)"] ["Cons(1, Cons(2, Nil))" "Nil"]]
           (map (fn [a] [(fo/show (a 'x)) (fo/show (a 'y))])
                (fo/fol-bc-ask append (fo/parse "Append(x, y, Cons(1, Cons(2, Nil)))")))))))

;; ---------------------------------------------------------------------------
;; Unification and clausal form

(defn- same? [theta a b] (= (fo/subst theta a) (fo/subst theta b)))

(deftest unify-test
  (testing "aima-python's unification tests, compared by what the substitution does"
    (doseq [[a b] [["A(x)" "A(B)"]
                   ["American(x) & Weapon(B)" "American(A) & Weapon(y)"]
                   ["P(F(x,z), G(u, z))" "P(F(y,a), y)"]
                   ["P(A, x, F(G(y)))" "P(z, F(z), F(u))"]
                   ["P(x, A, F(G(y)))" "P(F(z), z, F(u))"]]
            :let [x (fo/parse a) y (fo/parse b) th (fo/unify x y)]]
      (is (some? th) a)
      (is (same? th x y) a))
    (is (= {'x :B} (fo/unify (fo/parse "A(x)") (fo/parse "A(B)")))))
  (testing "the occurs check: x and F(x) do not unify"
    (is (nil? (fo/unify 'x (fo/parse "F(x)")))))
  (testing "two sentences that share a variable fail until standardized apart"
    (let [a (fo/parse "Knows(John, x)") b (fo/parse "Knows(x, Elizabeth)")]
      (is (nil? (fo/unify a b)))
      (is (some? (fo/unify a (fo/standardize-apart b)))))))

(deftest clausal-form-test
  (testing "an existential inside a universal becomes a function of it"
    (let [[c & more] (fo/clausal-form (fo/parse "forall x: exists y: Loves(y, x)"))
          [lit] (vec c)]
      (is (empty? more))
      (is (vector? (second lit)) "the lover is a Skolem function of x")
      (is (= (nth lit 2) (second (second lit))))))
  (testing "parse and show round-trip"
    (doseq [t (concat crime-text ["forall x: exists y: Loves(y, x)" "~(A(x) | B(x)) <=> C"])]
      (is (= (fo/parse t) (fo/parse (fo/show (fo/parse t))))))))
