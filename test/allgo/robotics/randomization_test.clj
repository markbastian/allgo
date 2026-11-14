(ns allgo.robotics.randomization-test
  "Domain randomization, after Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 26.8: gains
  searched in one nominal simulation against gains searched over many
  randomized ones, both tried on robots that are not the nominal one."
  (:require [allgo.random :as random]
            [allgo.robotics.randomization :as dr]
            [clojure.test :refer [deftest is testing]]))

(def ^:private nominal {:mass 1.0 :friction 0.0 :delay 3})

(deftest sim-to-real
  (let [tuned (dr/search-policy [nominal])
        robust (dr/search-policy (dr/randomized (random/rng 3) 30 [0.2 2.0] 0.5 3))]
    (testing "on the robot it was tuned for, the nominal policy is the better"
      (is (< (dr/cost nominal tuned) (dr/cost nominal robust))))
    (testing "on real robots lighter than modeled, it shakes itself apart where the randomized policy settles"
      (doseq [mass [0.3 0.5]]
        (let [real {:mass mass :friction 0.2 :delay 3}]
          (is (> (dr/cost real tuned) (* 5 (dr/cost real robust))) (str mass)))))
    (testing "and its worst case over the range of masses is far the worse"
      (let [worst (fn [gains] (apply max (for [m [0.3 0.5 1.0 1.5 2.0]] (dr/cost {:mass m :friction 0.2 :delay 3} gains))))]
        (is (> (worst tuned) (* 5 (worst robust))))))))
