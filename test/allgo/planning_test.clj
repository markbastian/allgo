(ns allgo.planning-test
  (:require [allgo.logic.first-order :as fo]
            [allgo.planning :as planning]
            [clojure.test :refer [deftest is testing]]))

;; The planning problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 11, as aima-pseudocode (CC BY
;; 4.0) states them.

(deftest classical-test
  (testing "the blocks world: the three moves the figure gives"
    (is (= ["MoveToTable(C, A)" "Move(B, Table, C)" "Move(A, Table, B)"]
           (mapv fo/show (planning/plan planning/blocks-world)))))
  (testing "every problem's plan is executable and reaches its goal, in the fewest actions"
    (doseq [[p n] [[planning/blocks-world 3] [planning/spare-tire 3] [planning/air-cargo 6]]
            :let [plan (planning/plan p)]]
      (is (= n (count plan)))
      (is (planning/holds? (planning/execute p plan) (:goal p))))))

(deftest scheduling-test
  (let [{:keys [actions capacity]} planning/car-assembly]
    (testing "without resources the critical path is engine 2, its wheels and inspection: 85"
      (let [{:keys [start length]} (planning/critical-path actions)]
        (is (= 85 length))
        (is (= #{:AddEngine2 :AddWheels2 :Inspect2} (set (keep (fn [[a [es ls]]] (when (= es ls) a)) start))))
        (is (= [0 15] (start :AddEngine1)) "engine 1 has 15 of slack")))
    (testing "with one engine hoist the best is 115 -- engine 1 first -- where minimum slack, starting with engine 2, takes 130"
      (is (= 115 (:length (planning/optimal-schedule actions capacity))))
      (is (= 130 (:length (planning/schedule actions capacity)))))))

;; Hierarchies: getting to the airport, and navigating a grid.

(def ^:private airport
  {:init {:at :Home :cash 0 :car? true}
   :goal? #(= :SFO (:at %))
   :step (fn [s [act from to]]
           (when (= from (:at s))
             (case act
               :Drive (when (:car? s) (assoc s :at to))
               :Shuttle (assoc s :at to)
               :Taxi (when (>= (:cash s) 50) (-> s (assoc :at to) (update :cash - 50))))))})

(def ^:private go-hierarchy
  {:Go (fn [[from to] _] [[[:Drive from :LongTermParking] [:Shuttle :LongTermParking to]]
                          [[:Taxi from to]]])})

(def ^:private grid
  {:init [0 0]
   :goal? #(= % [2 1])
   :step (fn [[x y] [a]]
           (let [[nx ny] (case a :Left [(dec x) y] :Right [(inc x) y] :Up [x (inc y)] :Down [x (dec y)] [x y])]
             (when (and (<= 0 nx 2) (<= 0 ny 2)) [nx ny])))})

(def ^:private navigate
  ;; Navigate([x y] [gx gy]): done, or a step toward the goal and the rest.
  {:Navigate (fn [[[gx gy]] [x y]]
               (cond-> []
                 (= [x y] [gx gy]) (conj [])
                 (< x gx) (conj [[:Right] [:Navigate [gx gy]]])
                 (> x gx) (conj [[:Left] [:Navigate [gx gy]]])
                 (< y gy) (conj [[:Up] [:Navigate [gx gy]]])
                 (> y gy) (conj [[:Down] [:Navigate [gx gy]]])))})

(deftest hierarchical-test
  (testing "with a car and no cash, the airport by drive and shuttle"
    (is (= [[:Drive :Home :LongTermParking] [:Shuttle :LongTermParking :SFO]]
           (planning/hierarchical-search airport go-hierarchy [:Go :Home :SFO]))))
  (testing "with cash and no car, by taxi"
    (is (= [[:Taxi :Home :SFO]]
           (planning/hierarchical-search (assoc airport :init {:at :Home :cash 80 :car? false}) go-hierarchy [:Go :Home :SFO]))))
  (testing "a recursive refinement: navigating the grid, three moves"
    (let [plan (planning/hierarchical-search grid navigate [:Navigate [2 1]])]
      (is (= 3 (count plan)))
      (is (= [2 1] (reduce (:step grid) (:init grid) plan))))))

(deftest angelic-test
  (testing "angelic search commits to the plan that surely works, and decomposes it to the same primitives"
    (is (= [[:Drive :Home :LongTermParking] [:Shuttle :LongTermParking :SFO]]
           (planning/angelic-search airport go-hierarchy [[:Go :Home :SFO]] 3)))
    (is (nil? (planning/angelic-search (assoc airport :init {:at :Home :cash 0 :car? false}) go-hierarchy [[:Go :Home :SFO]] 3))
        "with neither car nor cash, no plan")
    (let [plan (planning/angelic-search grid navigate [[:Navigate [2 1]]] 6)]
      (is (= [2 1] (reduce (:step grid) (:init grid) plan))))))
