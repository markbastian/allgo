(ns allgo.phenomena-search-test
  "`phenomena/search`, which finds the planets' phenomena from VSOP87,
  against the rows of Meeus's table 36.B that were checked against a
  second transcription -- which cover every kind it finds -- and then the
  phenomena those rows do not have, by how they must fall around the
  ones they do."
  (:require [allgo.astro.phenomena :as ph]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private checked
  "Meeus's verified rows, what `search` calls each, and how far apart the
  two may be, in minutes: Meeus's series are good to minutes for Mercury
  and Venus and hours for the outer planets, worst at a station, where
  the longitude is flat."
  {:mercury-inferior-conjunction [:mercury :inferior-conjunction 20]
   :mercury-superior-conjunction [:mercury :superior-conjunction 40]
   :mercury-greatest-east-elongation [:mercury :greatest-east-elongation 30]
   :mercury-greatest-west-elongation [:mercury :greatest-west-elongation 30]
   :venus-inferior-conjunction [:venus :inferior-conjunction 30]
   :mars-opposition [:mars :opposition 240]
   :mars-second-station [:mars :second-station 330]
   :jupiter-opposition [:jupiter :opposition 270]
   :saturn-opposition [:saturn :opposition 120]
   :saturn-conjunction [:saturn :conjunction 150]
   :uranus-opposition [:uranus :opposition 270]
   :neptune-opposition [:neptune :opposition 150]})

(deftest against-meeus-test
  (testing "the same phenomenon as each verified row of table 36.B, 1000-3000, within its accuracy"
    (doseq [[event [planet kind minutes]] checked
            y (range 1000 3001 41)]
      (let [r (ph/phenomenon event (double y))
            t (if (vector? r) (first r) r)
            s (ph/search planet kind (- t 20.0))]
        (is (some? s) (str event " " y))
        (when s
          (is (< (* 1440.0 (abs (- s t))) minutes) (str event " " y))))))
  (testing "Mercury's greatest elongations, the angle as Meeus has it within 0.005 degrees"
    (doseq [[event kind] [[:mercury-greatest-east-elongation :greatest-east-elongation]
                          [:mercury-greatest-west-elongation :greatest-west-elongation]]
            y (range 1000 3001 41)]
      (let [[t e] (ph/phenomenon event (double y))
            s (ph/search :mercury kind (- t 20.0))]
        (is (< (math/to-degrees (abs (- e (ph/elongation :mercury s)))) 0.005) (str event " " y))))))

(defn- days [a b] (- b a))

(deftest the-rest-test
  (doseq [y (range 1900 2101 11)
          :let [t0 (ph/phenomenon :venus-inferior-conjunction (double y))]]
    (testing "Venus: stations three weeks either side of inferior conjunction, greatest elongations
              ten weeks either side at 45-47 degrees, superior conjunction half a synodic period away"
      (let [first-station (ph/search :venus :first-station (- t0 60.0))
            second-station (ph/search :venus :second-station t0)
            east (ph/search :venus :greatest-east-elongation (- t0 120.0))
            west (ph/search :venus :greatest-west-elongation t0)
            superior (ph/search :venus :superior-conjunction t0)]
        (is (< 18 (days first-station t0) 24) (str y))
        (is (< 18 (days t0 second-station) 24) (str y))
        (is (< 66 (days east t0) 76) (str y))
        (is (< 66 (days t0 west) 76) (str y))
        (doseq [e [east west]]
          (is (< 45.0 (math/to-degrees (ph/elongation :venus e)) 47.5) (str y)))
        (is (< 280 (days t0 superior) 305) (str y)))))
  (doseq [y (range 1900 2101 7)]
    (testing "Mercury's stations either side of inferior conjunction"
      (let [t0 (ph/phenomenon :mercury-inferior-conjunction (double y))]
        (is (< 5 (days (ph/search :mercury :first-station (- t0 30.0)) t0) 17) (str y))
        (is (< 5 (days t0 (ph/search :mercury :second-station t0)) 17) (str y))))
    (testing "Mars's first station before opposition, its conjunction about a year from it"
      (let [t0 (ph/phenomenon :mars-opposition (double y))]
        (is (< 25 (days (ph/search :mars :first-station (- t0 60.0)) t0) 45) (str y))
        (is (< 330 (days t0 (ph/search :mars :conjunction t0)) 450) (str y))))
    (testing "the outer planets' conjunctions about half a synodic period after their oppositions"
      (doseq [[planet synodic] {:jupiter 398.88 :uranus 369.66 :neptune 367.49 :saturn 378.09}]
        (let [t0 (ph/phenomenon (keyword (str (name planet) "-opposition")) (double y))
              half (* 0.5 synodic)]
          (is (< (abs (- (days t0 (ph/search planet :conjunction t0)) half)) 20.0) (str planet " " y))))))
  (testing "backward finds what forward does"
    (let [t (ph/search :venus :greatest-west-elongation 60000.0)]
      (is (< (abs (- t (ph/search :venus :greatest-west-elongation (+ t 10.0) -1))) 1e-6)))))
