(ns allgo.robotics.reactive-test
  (:require [allgo.robotics.reactive :as reactive]
            [clojure.test :refer [deftest is testing]]))

;; Reactive control, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, section 26.9: a six-legged walker.

(defn- grounded [state] (count (filter reactive/grounded? (vals (:legs state)))))

(deftest gait-test
  (testing "on flat ground the tripods take turns, a stride every four ticks, and never fewer than three legs stand"
    (let [states (vec (take 201 (reactive/walk-seq (reactive/walker 0.1) (constantly 0.0) 0.1 0.05)))]
      (is (= 50 (:strides (peek states))))
      (is (every? #(>= (grounded %) 3) states))))
  (testing "a rock 0.27 high in the right front leg's way: its swing is blocked, it retracts and lifts 0.05 higher each time -- four times, until 0.3 clears it"
    (let [rock (fn [leg stride] (if (and (= leg :rf) (= stride 4)) 0.27 0.0))
          states (vec (take 201 (reactive/walk-seq (reactive/walker 0.1) rock 0.1 0.05)))]
      (is (= 4 (get-in (peek states) [:legs :rf :retries])))
      (is (every? #(>= (grounded %) 3) states))
      (testing "-- its tripod waiting for it, standing; each retry three ticks, blocked, lift and swing, so twelve ticks, three strides, lost"
        (is (= 47 (:strides (peek states))))))))

(deftest subsumption-test
  (let [layers [reactive/walk-layer reactive/avoid-layer]]
    (testing "with nothing touched, the walking layer's command goes through"
      (is (= :walk (reactive/subsumption layers {:bumper nil}))))
    (testing "the bumper's layer, above, suppresses it"
      (is (= :back-up-turn-right (reactive/subsumption layers {:bumper :left})))
      (is (= :back-up-turn-left (reactive/subsumption layers {:bumper :right}))))))
