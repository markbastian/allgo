(ns allgo.probability.bayes-test
  (:require [allgo.probability.bayes :as bn]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The burglary and sprinkler networks of Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 13, and the
;; answers aima-python's (MIT) tests expect of them.

(defn- close? [a b tol] (< (abs (- a b)) tol))

(deftest exact-test
  (testing "P(Burglary | JohnCalls, MaryCalls) = 0.284, by enumeration and by elimination"
    (doseq [ask [bn/enumeration-ask bn/elimination-ask]]
      (is (close? 0.284 ((ask :Burglary {:JohnCalls true :MaryCalls true} bn/burglary) true) 5e-4))))
  (testing "the two agree on every query of the sprinkler network"
    (doseq [x [:Cloudy :Sprinkler :Rain :WetGrass]
            e [{} {:WetGrass true} {:Sprinkler true :WetGrass true} {:Cloudy false}]
            :when (not (contains? e x))
            :let [a (bn/enumeration-ask x e bn/sprinkler) b (bn/elimination-ask x e bn/sprinkler)]]
      (is (close? (a true) (b true) 1e-12) [x e]))))

(deftest sampling-test
  (let [exact ((bn/enumeration-ask :Rain {:Sprinkler true} bn/sprinkler) true)]
    (testing "every sampler's estimate of P(Rain | Sprinkler) nears the exact 0.3"
      (doseq [ask [bn/rejection-sampling-seq bn/likelihood-weighting-seq bn/gibbs-seq]]
        (is (close? exact ((:estimate (nth (ask :Rain {:Sprinkler true} bn/sprinkler (random/rng 3)) 20000)) true) 0.02))))
    (testing "Gibbs sampling finds P(Cloudy | Rain) = 0.8, as aima-python's test does"
      (is (close? 0.8 ((:estimate (nth (bn/gibbs-seq :Cloudy {:Rain true} bn/sprinkler (random/rng 4)) 30000)) true) 0.02)))
    (testing "prior samples have the network's marginals"
      (let [rng (random/rng 5) xs (repeatedly 20000 #(bn/prior-sample bn/sprinkler rng))]
        (is (close? 0.5 (/ (count (filter :Rain xs)) 20000.0) 0.02))))))
