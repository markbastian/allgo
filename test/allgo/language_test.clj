(ns allgo.language-test
  (:require [allgo.language :as lang]
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
