(ns allgo.probability.identity-test
  "Citation matching, after Russell and Norvig, *Artificial Intelligence:
  A Modern Approach*, 4th edition, section 15.2.3, on the section's own
  two citations of Lashkari, Metral and Maes -- which, it says, are of
  the same paper."
  (:require [allgo.probability.identity :as identity]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(defn- fields
  "A citation's author and title fields, parsed, as words."
  [authors title]
  (identity/words (str authors " " title)))

;; the section's two citations, their author and title fields
(def ^:private lashkari-1 (fields "Yezdi Lashkari, Max Metral, and Pattie Maes" "Collaborative Interface Agents"))
(def ^:private lashkari-2 (fields "Metral M. Lashkari, Y. and P. Maes" "Collaborative interface agents"))
;; and two others
(def ^:private maes (fields "Maes, P." "Agents that reduce work and information overload"))
(def ^:private aima (fields "Russell, S. and Norvig, P." "Artificial Intelligence: A Modern Approach"))

(defn- posterior [citations] (identity/exact-posterior (identity/model citations {}) citations))

(deftest groupings
  (testing "the exact posterior covers every grouping, the Bell numbers of them, and sums to one"
    (doseq [[n bell] [[1 1] [2 2] [3 5] [4 15] [5 52]]]
      (let [cs (vec (take n [lashkari-1 lashkari-2 maes aima (fields "Maes" "Agents")]))
            p (posterior cs)]
        (is (= bell (count p)))
        (is (< (abs (- 1.0 (reduce + (vals p)))) 1e-9))))))

(deftest the-sections-citations
  (let [p (posterior [lashkari-1 lashkari-2 maes aima])]
    (testing "the two Lashkari citations, written differently, are of the same paper"
      (is (> (identity/same-paper p 0 1) 0.99)))
    (testing "a paper sharing an author and a word of title is not that paper"
      (is (< (identity/same-paper p 0 2) 0.01))
      (is (< (identity/same-paper p 2 3) 0.01))))
  (testing "the whole strings, venue and publisher and all, mislead: the section's point that text similarity fails, and that the model must parse"
    (let [whole [(identity/words "[Lashkari et al 94] Collaborative Interface Agents, Yezdi Lashkari, Max Metral, and Pattie Maes, Proceedings of the Twelfth National Conference on Articial Intelligence, MIT Press, Cambridge, MA, 1994.")
                 (identity/words "Metral M. Lashkari, Y. and P. Maes. Collaborative interface agents. In Conference of the American Association for Artificial Intelligence, Seattle, WA, August 1994.")]]
      (is (< (identity/same-paper (posterior whole) 0 1) 0.01)))))

(deftest collective-disambiguation
  (testing "\"the more citations for a given paper, the more accurately each of them is parsed\": a third citation sharing words with each of two makes them more surely one paper"
    (let [a (fields "Lashkari" "Collaborative interface agents")
          b (fields "Lashkari Metral Maes" "Interface agents learning")
          bridge (fields "Lashkari Metral Maes" "Collaborative interface agents learning")
          alone (identity/same-paper (posterior [a b maes]) 0 1)
          together (identity/same-paper (posterior [a b bridge maes]) 0 1)]
      (is (< alone 0.9))
      (is (> together 0.99)))))

(deftest gibbs
  (testing "Gibbs sweeps sample the exact posterior"
    (let [cs [lashkari-1 lashkari-2 maes aima (fields "Lashkari Metral Maes" "Collaborative agents for interfaces")]
          m (identity/model cs {})
          exact (identity/exact-posterior m cs)
          sampled (identity/sampled-posterior (take 4000 (drop 100 (identity/gibbs-seq m cs (random/rng 2)))))]
      (doseq [[i j] [[0 1] [0 4] [1 4] [0 2] [2 3]]]
        (is (< (abs (- (identity/same-paper exact i j) (identity/same-paper sampled i j))) 0.02) (str i " " j))))))
