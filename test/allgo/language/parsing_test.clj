(ns allgo.language.parsing-test
  (:require [allgo.language.parsing :as lang]
            [allgo.random :as random]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

;; Parsing, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, chapter 23, with aima-python's (MIT) grammar
;; E_Prob_Chomsky.

(defn- words [s] (str/split s #" "))

(deftest cyk-test
  (testing "the robot is good: S -> NP VP, each part by its likeliest rule -- 1 x 0.6 x 0.5 x 0.4 x 0.5 x 0.5 x 0.5"
    (let [[tree p] (lang/parse-tree (words "the robot is good") lang/e-prob-chomsky)]
      (is (= [:S [:NP [:Article "the"] [:Noun "robot"]] [:VP [:Verb "is"] [:Adjective "good"]]] tree))
      (is (< (abs (- p 0.015)) 1e-12))))
  (testing "a verb phrase with an object"
    (let [[tree p] (lang/parse-tree (words "the robot say the fence") lang/e-prob-chomsky)]
      (is (= ["the" "robot" "say" "the" "fence"] (lang/leaves tree)))
      (is (< (abs (- p (* 0.6 0.5 0.4 0.5 0.3 0.6 0.5 0.2))) 1e-12))))
  (testing "no parse for what the grammar cannot derive"
    (is (nil? (lang/parse-tree (words "a sad sheep say") lang/e-prob-chomsky))))
  (testing "a step a span: the lexicon, then the 10 spans of 2 to 5 words"
    (is (= 11 (count (lang/cyk-seq (words "the robot say the fence") lang/e-prob-chomsky))))))

(deftest treebank-test
  (testing "the chapter's Penn Treebank sentence, read and written back"
    (let [text "[[S [NP-SBJ-2 Her eyes] [VP were [VP glazed [NP *-2] [SBAR-ADV as if [S [NP-SBJ she] [VP did n't [VP [VP hear [NP *-1]] or [VP [ADVP even] see [NP *-1]] [NP-1 him]]]]]]]] .]"
          tree (lang/read-tree text)]
      (is (= text (lang/write-tree tree)))
      (is (= ["Her" "eyes" "were" "glazed" "*-2" "as" "if" "she" "did" "n't" "hear" "*-1" "or" "even" "see" "*-1" "him" "."]
             (lang/leaves tree))))))

(deftest generation-test
  (testing "sampled sentences are the grammar's: every one parses, and every word's tag is one the lexicon gives it"
    (let [rng (allgo.random/rng 1)
          sentences (repeatedly 50 #(lang/generate lang/ambiguous-grammar rng))]
      (is (every? #(lang/parse-tree (mapv first %) lang/ambiguous-grammar) sentences))
      (is (every? (fn [[w t]] (some #(= t (first %)) (get-in lang/ambiguous-grammar [:lexicon w])))
                  (apply concat sentences))))))

;; Learning a grammar, after the same chapter: from a treebank by
;; counting, and from unparsed sentences by the inside-outside algorithm.

(def ^:private g lang/ambiguous-grammar)

(defn- flat
  "A grammar's probabilities by rule: `{[X Y Z] p}` and `{[X word] p}`."
  [{:keys [rules lexicon]}]
  (merge (into {} (for [[x y z p] rules] [[x y z] p]))
         (into {} (for [[w xs] lexicon [x p] xs] [[x w] p]))))

(defn- worst-difference [a b]
  (let [fa (flat a) fb (flat b)]
    (apply max (for [k (into (set (keys fa)) (keys fb))] (abs (- (get fa k 0.0) (get fb k 0.0)))))))

(defn- with-probabilities
  "`binary` and `lexical` rules, `{[X Y Z] weight}` and `{[X word]
  weight}`, made a grammar: each left-hand side's weights normalized."
  [binary lexical]
  (let [totals (reduce (fn [m [[x] c]] (update m x (fnil + 0.0) c)) {} (concat binary lexical))]
    {:rules (vec (for [[[x y z] c] binary] [x y z (/ c (totals x))]))
     :lexicon (reduce (fn [m [[x w] c]] (update m w (fnil conj []) [x (/ c (totals x))])) {} lexical)}))

(deftest generated-trees
  (testing "a sampled tree's tagged leaves are the sentence the same draws give"
    (is (= (lang/generate g (random/rng 5)) (lang/tagged-leaves (lang/generate-tree g (random/rng 5)))))))

(deftest treebank-learning
  (let [rng (random/rng 1)
        trees (vec (repeatedly 3000 #(lang/generate-tree g rng)))
        learned (lang/learn-grammar trees)]
    (testing "counting the rules of three thousand trees gives every probability to within 0.02"
      (is (< (worst-difference g learned) 0.02)))
    (testing "and the learned grammar parses new sentences as the true one does"
      (let [test (repeatedly 50 #(lang/leaves (lang/generate-tree g rng)))]
        (is (<= 48 (count (filter #(= (first (lang/parse-tree % g)) (first (lang/parse-tree % learned))) test))))))))

(deftest inside-outside
  (let [rng (random/rng 2)]
    (testing "the rules known, their probabilities not: from uniform ones, EM climbs to the sample's own maximum"
      (let [sentences (vec (repeatedly 400 #(lang/leaves (lang/generate-tree g rng))))
            uniform (with-probabilities (into {} (for [[x y z] (:rules g)] [[x y z] 1.0]))
                      (into {} (for [[w xs] (:lexicon g) [x] xs] [[x w] 1.0])))
            states (vec (take 21 (lang/inside-outside-seq uniform sentences)))
            truth (:log-likelihood (first (lang/inside-outside-seq g sentences)))]
        (is (apply <= (map :log-likelihood states)) "the likelihood never falls")
        (is (> (worst-difference g uniform) (* 5 (worst-difference g (:grammar (peek states))))) "far nearer the truth")
        (is (< (worst-difference g (:grammar (peek states))) 0.05))
        (is (>= (:log-likelihood (peek states)) truth) "as likely as the true grammar, or more, on this sample")))
    (testing "every rule allowed among the categories, from random starts: it climbs, but to local maxima, as the chapter warns"
      (let [sentences (vec (repeatedly 100 #(lang/leaves (lang/generate-tree g rng))))
            heads [:S :NP :VP :VB] tags [:Det :N :V :Adj :Modal :Adv] cats (into heads tags)
            truth (:log-likelihood (first (lang/inside-outside-seq g sentences)))
            run (fn [seed]
                  (let [r (random/rng seed)
                        start (with-probabilities (into {} (for [x heads y cats z cats] [[x y z] (+ 0.1 (r))]))
                                (into {} (for [t tags w (keys (:lexicon g))] [[t w] (+ 0.1 (r))])))]
                    (mapv :log-likelihood (take 31 (lang/inside-outside-seq start sentences)))))
            [a b] (map run [1 3])]
        (is (apply <= (map #(- % 1e-6) a)))
        (is (apply <= (map #(- % 1e-6) b)))
        (is (< (peek a) (- truth 20)) "one start stalls well short of the true grammar")
        (is (> (peek b) (peek a)) "another, from elsewhere, does better")))))

;; Dependency parses, section 23.3.1 and figure 23.7.

(def ^:private detect
  "Figure 23.7's phrase structure parse of \"I detect the smelly wumpus
  near me\"."
  [:S [:NP [:Pronoun "I"]]
   [:VP [:Verb "detect"]
    [:NP [:NP [:Article "the"] [:Adjs [:Adjective "smelly"]] [:Noun "wumpus"]]
     [:PP [:Prep "near"] [:NP [:Pronoun "me"]]]]]])

(def ^:private figure-heads
  "The heads figure 23.7's dependency parse implies: a sentence's its verb
  phrase's, a verb phrase's its verb, a noun phrase's its noun or pronoun
  or inner noun phrase, a prepositional phrase's its noun phrase."
  {:S [:VP] :VP [:Verb] :NP [:Noun :Pronoun :NP] :PP [:NP] :Adjs [:Adjective]})

(deftest dependency-test
  (testing "figure 23.7: the phrase structure tree, annotated with heads, gives the dependency tree -- I and wumpus depend on detect; the, smelly and me on wumpus; near on me"
    (is (= [1 nil 4 4 1 6 4] (lang/dependencies detect figure-heads))))
  (testing "and back: a phrase structure tree, its categories made up, from which the same dependencies come"
    (let [tree (lang/dependency-tree (lang/leaves detect) [1 nil 4 4 1 6 4])]
      (is (= [:XP [:X "I"] [:H "detect"] [:XP [:X "the"] [:X "smelly"] [:H "wumpus"] [:XP [:X "near"] [:H "me"]]]] tree))
      (is (= (lang/leaves detect) (lang/leaves tree)))
      (is (= [1 nil 4 4 1 6 4] (lang/dependencies tree {:XP [:H]}))))))

;; A* and shift-reduce parsing, section 23.3.

(deftest a-star-parsing-test
  (testing "the first parse A* finds is the most probable: CYK's, on sampled sentences of the ambiguous grammar"
    (let [rng (random/rng 2)]
      (doseq [ws (repeatedly 30 #(mapv first (lang/generate lang/ambiguous-grammar rng)))]
        (let [[t1 p1] (lang/parse-tree ws lang/ambiguous-grammar)
              [t2 p2] (lang/a-star-parse ws lang/ambiguous-grammar)]
          (is (= t1 t2) ws)
          (is (< (abs (- p1 p2)) 1e-12) ws)))))
  (testing "the robot is good: 0.015, as CYK has it"
    (let [[tree p] (lang/a-star-parse (words "the robot is good") lang/e-prob-chomsky)]
      (is (= [:S [:NP [:Article "the"] [:Noun "robot"]] [:VP [:Verb "is"] [:Adjective "good"]]] tree))
      (is (< (abs (- p 0.015)) 1e-12))))
  (testing "no parse where the grammar has none"
    (is (nil? (lang/a-star-parse (words "a sad sheep say") lang/e-prob-chomsky)))))

(deftest shift-reduce-test
  (let [rng (random/rng 2)
        sentences (vec (repeatedly 30 #(mapv first (lang/generate lang/ambiguous-grammar rng))))
        outcome (fn [b] (frequencies (for [ws sentences
                                           :let [best (first (lang/parse-tree ws lang/ambiguous-grammar))
                                                 found (lang/shift-reduce-parse ws lang/ambiguous-grammar b)]]
                                       (cond (nil? found) :none (= best (first found)) :best :else :other))))]
    (testing "a deterministic parser, b = 1, commits to its likeliest step and mostly finds no parse at all; wider beams find more, and with 16 every parse is the most probable"
      (is (< (get (outcome 1) :best 0) 10))
      (is (< (get (outcome 1) :best 0) (get (outcome 4) :best 0)))
      (is (= {:best 30} (outcome 16))))
    (testing "every parse takes 2n - 1 actions: the shifts of the n words, as their parts of speech, and n - 1 reductions"
      (is (= 10 (count (lang/shift-reduce-seq (words "the dog sees a cat") lang/ambiguous-grammar 4)))))))
