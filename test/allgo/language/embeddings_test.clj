(ns allgo.language.embeddings-test
  (:require [allgo.language.embeddings :as emb]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Word embeddings, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 24, on a corpus built so that the
;; analogies have answers: each person word is a gender, a rank and an age,
;; each said by context words of its own.

(defn- corpus [rng n]
  (let [people {"king" [:m :royal :old] "queen" [:f :royal :old] "prince" [:m :royal :young] "princess" [:f :royal :young]
                "man" [:m :common :old] "woman" [:f :common :old] "boy" [:m :common :young] "girl" [:f :common :young]}
        verbs {:royal ["rules" "crowns" "reigns"] :common ["works" "plows" "sells"]}
        things {:royal ["throne" "crown" "kingdom"] :common ["field" "cart" "goods"]}
        pron {:m ["his" "he"] :f ["her" "she"]}
        age {:old ["grown" "elder"] :young ["little" "small"]}]
    (vec (repeatedly n (fn [] (let [[w [g r a]] (random/pick rng (vec people))]
                                [(random/pick rng (age a)) w (random/pick rng (verbs r)) (random/pick rng (pron g)) (random/pick rng (things r))]))))))

(def ^:private model
  (delay (let [rng (random/rng 1)] (nth (emb/sgns-seq (corpus rng 2000) rng :dim 12 :window 2) 8))))

(deftest similarity-test
  (testing "words used alike point alike"
    (is (> (emb/similarity @model "rules" "reigns") 0.8))
    (is (> (emb/similarity @model "rules" "reigns") (emb/similarity @model "rules" "plows")))
    (is (contains? #{"queen" "prince" "man"} (ffirst (emb/nearest @model "king" 1)))
        "the king's nearest differs from him in one way only")))

(deftest analogy-test
  (testing "a : b :: c : d by vector arithmetic, d the nearest to b - a + c"
    (is (= "queen" (emb/analogy @model "man" "woman" "king")))
    (is (= "girl" (emb/analogy @model "king" "queen" "boy")))
    (is (= "queen" (emb/analogy @model "man" "king" "woman")))
    (is (= "boy" (emb/analogy @model "woman" "girl" "man")))
    (is (= "queen" (emb/analogy @model "prince" "king" "princess")))
    (is (= "princess" (emb/analogy @model "boy" "prince" "girl")))))

(deftest projection-test
  (testing "the projection onto two principal components keeps the points' spread along them"
    (let [pts (emb/project-2d [[1.0 0.0 0.0] [-1.0 0.0 0.0] [0.0 0.5 0.0] [0.0 -0.5 0.0]])]
      (is (= 4 (count pts)))
      (is (< (abs (- 1.0 (apply max (map #(abs (first %)) pts)))) 1e-9))
      (is (< (abs (- 0.5 (apply max (map #(abs (second %)) pts)))) 1e-9)))))
