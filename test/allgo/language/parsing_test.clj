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
