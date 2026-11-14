(ns allgo.planning.satplan-test
  (:require [allgo.planning.satplan :as satplan]
            [clojure.test :refer [deftest is testing]]))

(def ^:private cake
  "The 'have cake and eat cake too' problem of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, as their pseudocode
  repository (aimacode/aima-pseudocode) states it."
  {:fluents [:Have :Eaten]
   :init #{:Have}
   :goal #{:Have :Eaten}
   :actions [{:name :Eat :pre #{:Have} :effect #{[:not :Have] :Eaten}}
             {:name :Bake :pre #{[:not :Have]} :effect #{:Have}}]})

(def ^:private rooms
  "Three rooms in a row and a locked door: the key is in the first."
  {:fluents [:InA :InB :InC :Key :Open]
   :init #{:InA}
   :goal #{:InC}
   :actions [{:name :GoAB :pre #{:InA} :effect #{[:not :InA] :InB}}
             {:name :GoBA :pre #{:InB} :effect #{[:not :InB] :InA}}
             {:name :TakeKey :pre #{:InA} :effect #{:Key}}
             {:name :Unlock :pre #{:InB :Key} :effect #{:Open}}
             {:name :GoBC :pre #{:InB :Open} :effect #{[:not :InB] :InC}}]})

(deftest satplan-test
  (testing "eat the cake, then bake another"
    (is (= [:Eat :Bake] (satplan/satplan cake 5))))
  (testing "the shortest plan: no plan of three steps exists, one of four does"
    (is (= [:TakeKey :GoAB :Unlock :GoBC] (satplan/satplan rooms 8)))
    (is (nil? (satplan/satplan rooms 3))))
  (testing "a goal no plan reaches"
    (is (nil? (satplan/satplan (assoc rooms :actions (pop (:actions rooms))) 6)))))
