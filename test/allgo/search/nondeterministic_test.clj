(ns allgo.search.nondeterministic-test
  (:require [allgo.random :as random]
            [allgo.search.nondeterministic :as nd]
            [clojure.test :refer [deftest is testing]]))

;; The slippery and local-sensing vacuum worlds of Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 4.
;; A state is [where left-dirty? right-dirty?].

(defn- result [[w l r] a]
  (case a :left [:left l r] :right [:right l r] :suck (if (= w :left) [w false r] [w l false])))

(defn- clean? [[_ l r]] (not (or l r)))

(def ^:private slippery
  "Moving sometimes fails: the agent stays where it was."
  {:initial [:left true true] :actions (constantly [:suck :left :right])
   :results (fn [s a] (if (#{:left :right} a) (hash-set (result s a) s) #{(result s a)}))
   :goal? clean?})

(deftest cyclic-test
  (testing "no acyclic plan: moving may fail every time"
    (is (nil? (nd/and-or slippery))))
  (testing "a cyclic one: suck, keep trying right until it works, suck"
    (let [plan (nd/cyclic-and-or slippery)]
      (is (= [:suck {[:left false true] [:right {[:left false true] [:loop [:left false true]]
                                                 [:right false true] [:suck {[:right false false] []}]}]}]
             plan))
      (testing "and it gets there however the world slips"
        (let [rng (random/rng 2)]
          (doseq [_ (range 200)]
            (is (clean? (peek (nd/execute slippery plan (:initial slippery) #(random/pick rng %))))))))))
  (testing "a plan that could only loop is not a plan"
    (is (nil? (nd/cyclic-and-or (assoc slippery :actions (constantly [:left])))))))

(deftest partially-observable-test
  (testing "local sensing, the right square unknown: suck, go right, and suck there only if it is dirty"
    (let [problem (nd/partially-observable
                   {:initial [[:left true true] [:left true false]]
                    :actions (constantly [:suck :left :right])
                    :results (fn [s a] #{(result s a)})
                    :percept (fn [[w l r]] [w (if (= w :left) l r)])
                    :goal? clean?})]
      (is (= [:suck {#{[:left false false] [:left false true]}
                     [:right {#{[:right false false]} []
                              #{[:right false true]} [:suck {#{[:right false false]} []}]}]}]
             (nd/and-or problem))))))
