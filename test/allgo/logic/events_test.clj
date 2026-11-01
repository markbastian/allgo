(ns allgo.logic.events-test
  (:require [allgo.logic.events :as ev]
            [clojure.test :refer [deftest is testing]]))

;; The worked examples of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 3rd edition, section 12.3 (4th edition, chapter 10)

(deftest the-time-scale
  (testing "the chapter's numbers hold on Pacific Standard Time, eight hours behind Greenwich"
    (let [ad2001 (ev/year-interval 2001 :utc-offset -8)]
      (is (= [3187324800 3218860800] ad2001))
      (is (= 31536000 (ev/duration ad2001)))
      (is (= (* 8 3600) (first (ev/year-interval 1900 :utc-offset -8))) "the scale's zero stays midnight at Greenwich")))
  (testing "Date(0, 20, 21, 24, 1, 1995) = Seconds(3000000000): seconds first, then minutes, then hours"
    (is (= 3000000000 (ev/date 0 20 21 24 1 1995 :utc-offset -8)))
    (is (not= 3000000000 (ev/date 0 20 21 24 1 1995)) "not on Greenwich time, which the text names"))
  (testing "on Greenwich time AD2001 starts eight hours earlier"
    (is (= (- 3187324800 (* 8 3600)) (first (ev/year-interval 2001)))))
  (testing "leap years"
    (is (= (* 366 86400) (ev/duration (ev/year-interval 2000))))
    (is (= (* 365 86400) (ev/duration (ev/year-interval 1900))))))

(deftest interval-relations
  (testing "the chapter's definitions"
    (is (= #{:meet} (ev/relations-of [0 5] [5 9])))
    (is (= #{:before} (ev/relations-of [0 4] [5 9])))
    (is (= #{:after} (ev/relations-of [5 9] [0 4])))
    (is (= #{:during} (ev/relations-of [2 3] [0 9])))
    (is (= #{:overlap} (ev/relations-of [0 4] [2 9])))
    (is (= #{} (ev/relations-of [2 9] [0 4])) "Overlap only one way round: i begins first")
    (is (= #{:begins} (ev/relations-of [0 4] [0 9])))
    (is (= #{:finishes} (ev/relations-of [5 9] [0 9])))
    (is (= #{:begins :finishes :equals} (ev/relations-of [0 9] [0 9]))))
  (testing "the reign of Elizabeth II immediately followed that of George VI"
    (let [george [1936 1952] elizabeth [1952 2022]]
      (is (ev/meet? george elizabeth)))))

(def ^:private presidents
  "Figure 12.3's President(USA): Washington from 1789 to 1797, Adams from
  1797 to 1801, Jefferson from 1801 -- years as the time points."
  {:events [{:id :inauguration1 :interval [1789 1789]}
            {:id :inauguration2 :interval [1797 1797]}
            {:id :inauguration3 :interval [1801 1801]}]
   :initiates (fn [e _] (case (:id e)
                          :inauguration1 [[:Equals :President-USA :GeorgeWashington]]
                          :inauguration2 [[:Equals :President-USA :JohnAdams]]
                          :inauguration3 [[:Equals :President-USA :ThomasJefferson]]))
   :terminates (fn [e _] (case (:id e)
                           :inauguration1 []
                           :inauguration2 [[:Equals :President-USA :GeorgeWashington]]
                           :inauguration3 [[:Equals :President-USA :JohnAdams]]))})

(deftest event-calculus
  (testing "T(Equals(President(USA), GeorgeWashington), AD1790)"
    (is (ev/holds-over presidents [:Equals :President-USA :GeorgeWashington] [1790 1791]))
    (is (false? (ev/holds presidents [:Equals :President-USA :GeorgeWashington] 1798)))
    (is (true? (ev/holds presidents [:Equals :President-USA :JohnAdams] 1798)))
    (is (nil? (ev/holds presidents [:Equals :President-USA :JohnAdams] 1790)) "nothing yet says either way"))
  (testing "the wumpus agent's arrow: Initiates only at Start, Terminates by a Shooting"
    (let [n {:events [{:id :Start :interval [0 0]} {:id :Shoot1 :category :Shootings :interval [5 6]}]
             :initiates (fn [e _] (if (= :Start (:id e)) [:HaveArrow] []))
             :terminates (fn [e _] (if (= :Shootings (:category e)) [:HaveArrow] []))}]
      (is (= [nil true true true true true false false] (mapv #(ev/holds n :HaveArrow %) (range 8))))
      (is (ev/clipped? n :HaveArrow [0 6]))
      (is (not (ev/clipped? n :HaveArrow [0 5])))
      (testing "the timeline gives the same answers a step at a time"
        (is (= (mapv #(ev/holds n :HaveArrow %) (range 1 12))
               (mapv (fn [{:keys [holds not-holds]}] (cond (holds :HaveArrow) true (not-holds :HaveArrow) false))
                     (take 11 (rest (ev/timeline-seq n 0)))))))))
  (testing "an initiation and a termination at once leave the fluent undecided"
    (let [n {:events [{:id :a :interval [0 1]} {:id :b :interval [0 1]}]
             :initiates (fn [e _] (if (= :a (:id e)) [:F] []))
             :terminates (fn [e _] (if (= :b (:id e)) [:F] []))}]
      (is (nil? (ev/holds n :F 3)))
      (is (contains? (:undecided (nth (ev/timeline-seq n 0) 3)) :F)))))

(deftest processes
  (testing "a process happens over every subinterval; a discrete event only over its own"
    (let [flight {:id :E1 :category :Flyings :interval [0 10]}]
      (is (ev/happens? flight [3 7] #{:Flyings}))
      (is (not (ev/happens? flight [3 7])))
      (is (ev/happens? flight [0 10]))
      (is (not (ev/happens? flight [3 12] #{:Flyings}))))))
