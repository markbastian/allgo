(ns allgo.meeus-misc-test
  "Meeus chapters 15 and 56 to 58 -- rising and setting, magnitudes,
  binary stars and sundials -- against the book's worked examples."
  (:require [allgo.astro.rise :as rise]
            [allgo.astro.stars :as stars]
            [allgo.astro.sundial :as dial]
            [allgo.astro.time :as t]
            [allgo.meeus-support :refer [->arcsec ->deg close? deg dms hms]]
            [clojure.test :refer [deftest is testing]]))

(def ^:private boston [(dms 42 20 0) (dms -1 71 5 0)])

(deftest rising-and-setting
  (testing "15.a: Venus from Boston, 1988 March 20, the first approximation"
    (let [{:keys [rise transit set]} (rise/approximate boston (:star rise/standard-altitude)
                                                       (hms 11 50 58.1)
                                                       [(hms 2 46 55.51) (dms 18 26 27.3)])]
      (is (close? 0.51816 rise 1e-5))
      (is (close? 0.81965 transit 1e-5))
      (is (close? 0.12113 set 1e-5))))
  (testing "15.a: corrected for Venus's motion"
    (let [{:keys [rise transit set]}
          (rise/times boston (deg -0.5667) (hms 11 50 58.1) 56.0
                      [[(hms 2 42 43.25) (dms 18 2 51.4)]
                       [(hms 2 46 55.51) (dms 18 26 27.3)]
                       [(hms 2 51 7.69) (dms 18 49 38.7)]])]
      (is (close? 0.51766 rise 1e-5))
      (is (close? 0.81980 transit 1e-5))
      (is (close? 0.12130 set 1e-5))))
  (testing "and from the ephemeris alone"
    (let [{:keys [rise transit set]} (rise/planet :venus boston (t/calendar->mjd 1988 3 20))]
      (is (close? 0.51766 rise 2e-5))
      (is (close? 0.81980 transit 2e-5))
      (is (close? 0.12130 set 2e-5))))
  (testing "the Sun at the equinox is up for about half the day"
    (let [{:keys [rise set]} (rise/sun [0.0 0.0] (t/calendar->mjd 2024 3 20))]
      (is (close? 0.5 (- set rise) 0.01))))
  (testing "a circumpolar star never sets"
    (is (nil? (rise/approximate boston (:star rise/standard-altitude) 0.0 [0.0 (deg 80)])))))

(deftest magnitudes
  (testing "56.a: Castor's two components"
    (is (close? 1.58 (stars/combined-magnitude 1.96 2.89) 0.005)))
  (testing "56.b: a triple"
    (is (close? 3.93 (stars/combined-magnitude 4.73 5.22 5.60) 0.005)))
  (testing "56.c: a cluster"
    (is (close? 2.02 (apply stars/combined-magnitude
                            (concat (repeat 4 5) (repeat 14 6) (repeat 23 7) (repeat 38 8)))
                0.005)))
  (testing "56.d, 56.e: ratios"
    (is (close? 6.19 (stars/brightness-ratio 0.14 2.12) 0.005))
    (is (close? 6.75 (stars/magnitude-difference 500) 0.005)))
  (testing "absolute magnitude is the apparent one at 10 parsecs"
    (is (close? 4.0 (stars/absolute-magnitude 4.0 10.0) 1e-12))
    (is (close? (stars/absolute-magnitude 1.0 20.0)
                (stars/absolute-magnitude-from-parallax 1.0 (dms 0 0 0.05)) 1e-9))))

(deftest binary-stars
  (let [eta-cor {:P 41.623 :T 1934.008 :e 0.2763 :a (dms 0 0 0.907)
                 :i (deg 59.025) :om (deg 23.717) :w (deg 219.907)}]
    (testing "57.a: eta Coronae Borealis in 1980"
      (let [[theta rho] (stars/binary-position eta-cor 1980.0)]
        (is (close? 318.4 (->deg theta) 0.05))
        (is (close? 0.411 (->arcsec rho) 5e-4))))
    (testing "57.b: its apparent eccentricity"
      (is (close? 0.860 (stars/apparent-eccentricity 0.2763 (deg 59.025) (deg 219.907)) 5e-4)))))

(deftest sundials
  (testing "58.a: a declining, inclined dial at latitude 40"
    (let [{:keys [lines center psi]} (dial/general (deg 40) (deg 70) 1.0 (deg 50))
          by-hour (into {} (map (juxt :hour :points)) lines)]
      (is (= [9 10 11 12 13 14 15 16 17 18 19] (mapv :hour lines)))
      (is (every? true? (map #(close? %1 %2 5e-5) [-2.0007 -1.1069] (get-in by-hour [11 2]))))
      (is (every? true? (map #(close? %1 %2 5e-5) [-0.0390 -0.3615] (get-in by-hour [14 6]))))
      (is (every? true? (map #(close? %1 %2 5e-5) [3.3880 -3.1102] center)))
      (is (close? 12.2672 (->deg psi) 5e-5))))
  (testing "58.b: a vertical dial in the southern hemisphere"
    (let [{:keys [lines center psi]} (dial/general (deg -35) (deg 160) 1.0 (deg 90))
          by-hour (into {} (map (juxt :hour :points)) lines)]
      (is (every? true? (map #(close? %1 %2 5e-5) [0.3640 -0.7410] (get-in by-hour [12 5]))))
      (is (every? true? (map #(close? %1 %2 5e-5) [-0.8439 -0.9298] (get-in by-hour [15 3]))))
      (is (every? true? (map #(close? %1 %2 5e-5) [0.3640 0.7451] center)))
      (is (close? 50.3315 (->deg psi) 5e-5))
      (testing "and it agrees with the vertical dial's own formulae"
        (is (= (map :hour lines) (map :hour (:lines (dial/vertical (deg -35) (deg 160) 1.0))))))))
  (testing "58.c: which hours a steep north-facing dial sees"
    (is (= [5 6 13 14 15 16 17 18 19] (mapv :hour (:lines (dial/general (deg 40) (deg 160) 1.0 (deg 75)))))))
  (testing "a horizontal dial is the general one with z = 0"
    (let [h (dial/horizontal (deg 40) 1.0)
          g (dial/general (deg 40) 0.0 1.0 0.0)]
      (is (= (map :hour (:lines h)) (map :hour (:lines g))))
      (is (every? true? (map #(close? %1 %2 1e-12)
                             (flatten (map :points (:lines h))) (flatten (map :points (:lines g)))))))))
