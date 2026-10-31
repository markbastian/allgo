(ns allgo.logic.wumpus-test
  (:require [allgo.logic.propositional :as pl]
            [allgo.logic.wumpus :as wumpus]
            [clojure.test :refer [deftest is testing]]))

(defn- truth
  "The model the true cave makes of the knowledge base's symbols."
  [{:keys [pits] :as w}]
  (into {} (for [x (range 1 5) y (range 1 5)
                 :let [sq [x y] nb (wumpus/neighbors sq)]
                 [k v] [[(wumpus/pit sq) (contains? pits sq)]
                        [(wumpus/wumpus sq) (= sq (:wumpus w))]
                        [(wumpus/breeze sq) (boolean (some pits nb))]
                        [(wumpus/stench sq) (boolean (some #{(:wumpus w)} nb))]]]
             [k v])))

(deftest world-test
  (testing "the physics is true of every cave, and the percepts agree with it"
    (doseq [seed (range 1 30)
            :let [w (wumpus/world seed) m (truth w)]]
      (is (every? #(true? (pl/clause-true? % m)) wumpus/physics))
      (doseq [x (range 1 5) y (range 1 5)
              :let [p (wumpus/percept (assoc (wumpus/start w) :pos [x y]))]]
        (is (= (:breeze p) (m (wumpus/breeze [x y]))))
        (is (= (:stench p) (m (wumpus/stench [x y])))))))
  (testing "the scores: -1 an action, -1000 for a pit, +1000 for leaving with the gold"
    (let [w {:pits #{[2 1]} :wumpus [4 4] :gold [1 2]}
          s (wumpus/start w)]
      (is (= -1001 (:score (wumpus/execute s :forward))))
      (is (= -5 (:score (reduce wumpus/execute s [:left :forward :grab :right :right]))))
      (is (= (+ -7 1000) (:score (reduce wumpus/execute s [:left :forward :grab :right :right :forward :climb]))))
      (is (= -1 (:score (wumpus/execute s :climb))) "climbing out empty-handed scores nothing more"))))

(deftest agent-test
  (let [games (for [seed (range 1 41)]
                (let [ss (vec (wumpus/run-seq (wumpus/world seed)))]
                  {:world (wumpus/world seed) :states ss :end (peek ss)}))]
    (testing "every game ends, with the agent out of the cave or dead"
      (is (every? #(wumpus/over? (:state (:end %))) games)))
    (testing "what the agent knows is true of the cave: every square it proved safe is safe, every pit it proved is a pit"
      (doseq [{:keys [world end]} games
              :let [{:keys [safe pits]} (:memory end)]]
        (is (not-any? (:pits world) safe))
        (is (every? (:pits world) pits))))
    (testing "it never dies in a square it had proved safe: deaths are risks it chose"
      (doseq [{:keys [states end]} games
              :when (not (:alive? (:state end)))
              :let [before (states (- (count states) 2))]]
        (is (not (contains? (:safe (:memory before)) (:pos (:state end)))))))
    (testing "and it brings the gold out of some of them"
      (is (some #(and (:gold? (:state (:end %))) (:out? (:state (:end %)))) games)))))
