(ns allgo.almanac-test
  (:require [allgo.astro.almanac :as al]
            [allgo.astro.time :as t]
            [clojure.test :refer [deftest is testing]]))

(defn- on-day? [[y m d] mjd]
  (let [[y' m' d'] (t/mjd->calendar mjd)]
    (= [y m d] [y' m' (long d')])))

(deftest next-and-previous
  (testing "the eclipses of 2024 and 2025"
    (is (on-day? [2024 4 8] (al/next-event :solar-eclipse (t/calendar->mjd 2024 1 1))))
    (is (on-day? [2025 3 14] (al/next-event :lunar-eclipse (t/calendar->mjd 2025 1 1))))
    (is (on-day? [2024 4 8] (al/previous-event :solar-eclipse (t/calendar->mjd 2024 5 1)))))
  (testing "a new moon either side of a date"
    (is (on-day? [2024 1 11] (al/previous-event :new-moon (t/calendar->mjd 2024 1 15))))
    (is (on-day? [2024 2 9] (al/next-event :new-moon (t/calendar->mjd 2024 1 15)))))
  (testing "the seasons, across the turn of a year"
    (is (on-day? [2025 3 20] (al/next-event :march-equinox (t/calendar->mjd 2024 6 1))))
    (is (on-day? [2024 3 20] (al/previous-event :march-equinox (t/calendar->mjd 2024 6 1)))))
  (testing "Mars's opposition of 2025"
    (is (on-day? [2025 1 16] (al/next-event :mars-opposition (t/calendar->mjd 2024 6 1)))))
  (testing "every event can be found in both directions, in order"
    (doseq [k (keys al/events)]
      (let [now (t/calendar->mjd 2030 7 1)
            n (al/next-event k now)
            p (al/previous-event k now)]
        (is (and n p (< p now n)) (str k))
        (is (< (- n p) (* 1.5 (max 800 (:period (al/events k))))) (str k))))))
