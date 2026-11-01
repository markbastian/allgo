(ns allgo.decision.mdp-test
  (:require [allgo.decision.mdp :as mdp]
            [clojure.test :refer [deftest is testing]]))

;; The 4 x 3 world of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 17, and the expected answers, as
;; aima-python (MIT) and its tests give them.

(deftest four-by-three-test
  (testing "value iteration's utilities are aima-python's"
    (let [U (mdp/value-iteration mdp/four-by-three 0.01)
          expect {[3 2] 1.0 [3 1] -1.0 [3 0] 0.12958868267972745 [0 1] 0.39810203830605462
                  [0 2] 0.50928545646220924 [1 0] 0.25348746162470537 [0 0] 0.29543540628363629
                  [1 2] 0.64958064617168676 [2 0] 0.34461306281476806 [2 1] 0.48643676237737926
                  [2 2] 0.79536093684710951}]
      (doseq [[s v] expect] (is (< (abs (- v (U s))) 1e-12) s))))
  (testing "policy iteration's policy is aima-python's"
    (is (= {[0 0] [0 1] [0 1] [0 1] [0 2] [1 0] [1 0] [1 0] [1 2] [1 0] [2 0] [0 1]
            [2 1] [0 1] [2 2] [1 0] [3 0] [-1 0] [3 1] nil [3 2] nil}
           (mdp/policy-iteration mdp/four-by-three))))
  (testing "and the two agree"
    (is (= (mdp/policy-iteration mdp/four-by-three)
           (mdp/best-policy mdp/four-by-three (mdp/value-iteration mdp/four-by-three 0.0001))))))

(def ^:private pomdp
  "aima-python's POMDP test example: two states, three actions."
  {:transition [[[0.65 0.35] [0.65 0.35]] [[0.65 0.35] [0.65 0.35]] [[1.0 0.0] [0.0 1.0]]]
   :sensor [[[0.5 0.5] [0.5 0.5]] [[0.5 0.5] [0.5 0.5]] [[0.8 0.2] [0.3 0.7]]]
   :reward [[5 -10] [-20 5] [-1 -1]]
   :gamma 0.95})

(defn- lookahead
  "Expectimax over beliefs to depth `d`: an independent account of what the
  plans of depth `d` are worth."
  [{:keys [transition sensor reward gamma] :as p} b d]
  (if (zero? d)
    0.0
    (apply max (for [a (range 3)]
                 (+ (reduce + (map * b (reward a)))
                    (* gamma (reduce + (for [e (range 2)
                                             :let [pe (reduce + (for [s (range 2) s2 (range 2)]
                                                                  (* (b s) (get-in transition [a s s2]) (get-in sensor [a s2 e]))))]
                                             :when (pos? pe)]
                                         (* pe (lookahead p (mdp/update-belief p b a e) (dec d)))))))))))

(deftest pomdp-test
  (testing "belief update, as aima-python's test has it"
    (let [[b0 b1] (mdp/update-belief pomdp [0.5 0.5] 2 0)]
      (is (< (abs (- b0 (/ 0.4 0.55))) 1e-12))
      (is (< (abs (- b1 (/ 0.15 0.55))) 1e-12))))
  (testing "the plans of each depth are worth what expectimax over beliefs says, everywhere"
    (doseq [{:keys [plans depth]} (take 6 (mdp/pomdp-value-iteration-seq pomdp 0.1))
            k (range 21)
            :let [b [(/ k 20.0) (- 1.0 (/ k 20.0))]]]
      (is (< (abs (- (mdp/pomdp-value plans b) (lookahead pomdp b depth))) 1e-9) [depth b]))))
