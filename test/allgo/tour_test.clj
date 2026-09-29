(ns allgo.tour-test
  "Gravity-assist tours: the powered flyby flown through its hyperbolas,
  and three real tours flown on their real dates -- Galileo, Cassini and
  Voyager 2 -- against what NASA published of them.

  Dates of launch and closest approach, UTC, and the distances: Galileo
  from the PDS mission catalog (go_1002/catalog/mission.cat) and JPL's
  press releases (Venus 16,000 km, Earth 960 km and 303 km); Cassini from
  JPL's VVEJGA trajectory sheet (JPL 400-856D; the December 1998 deep
  space maneuver, 'close to 450 meters per second' in NASA's accounts)
  and NASA's Cassini timeline (Venus 284 km and about 600 km, Earth about
  1,100 km, Jupiter about 10 million km); Voyager 2 from NSSDCA
  (1977-076A) and NASA Science and History (Jupiter 350,000 miles and
  Uranus 81,500 km above the cloud tops).

  The patched conics, on mean planetary elements, fly each of them nearly
  ballistically -- small burns at the flybys stand in for the missions'
  trajectory corrections and the model's own approximations -- and pass
  close to the published distances where the geometry fixes them."
  (:require [allgo.astro.interplanetary :as ip]
            [allgo.astro.time :as time]
            [allgo.astro.tour :as tour]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- utc [y m d h mi] (time/utc->tt (+ (time/calendar->mjd y m d) (/ (+ h (/ mi 60.0)) 24.0))))

(defn- close? [a b tol] (every? true? (map #(< (abs (- %1 %2)) tol) a b)))

(deftest powered-flyby
  (testing "flown: arrive on one hyperbola, burn at periapsis, leave on the other"
    (let [mu (tour/gm :earth)
          v-in [3.0 -7.5 1.2]
          ;; the leaving excess velocity: faster, and turned 40 degrees in
          ;; some plane holding the arriving one
          h (v3/normalize (v3/cross v-in [0.2 0.3 1.0]))
          turn (math/to-radians 40.0)
          dir (v3/add (v3/scale (v3/normalize v-in) (math/cos turn))
                      (v3/scale (v3/cross h (v3/normalize v-in)) (math/sin turn)))
          v-out (v3/scale dir 8.6)
          {:keys [rp dv]} (tour/powered-flyby mu v-in v-out)
          [r v] (ip/approach-hyperbola mu rp v-in h)
          v+ (v3/scale v (/ (+ (v3/length v) dv) (v3/length v)))
          [_ v-far] (u/propagate mu r v+ (* 3.0e7))]
      (is (< 0.0 dv 1.0))
      ;; far out, as nearly as the 1/r pull allows
      (is (close? v-far v-out 2e-3))))
  (testing "unpowered: the same speed turned by what the periapsis gives"
    (let [mu (tour/gm :venus)
          rp 6400.0
          v-in [0.0 6.0 0.0]
          {:keys [v-out]} (ip/flyby mu [0.0 0.0 0.0] v-in rp [0.0 0.0 1.0])
          {r :rp dv :dv} (tour/powered-flyby mu v-in v-out)]
      (is (< (abs (- r rp)) 1e-3))
      (is (< dv 1e-9)))))

(def ^:private galileo
  [[:earth :venus :earth :earth :jupiter]
   [(utc 1989 10 18 16 54) (utc 1990 2 10 5 59) (utc 1990 12 8 20 35) (utc 1992 12 8 15 9) (utc 1995 12 7 22 0)]
   {:long? [false true false true]}])

(def ^:private cassini
  [[:earth :venus :venus :earth :jupiter :saturn]
   [(utc 1997 10 15 8 43) (utc 1998 4 26 13 45) (utc 1999 6 24 20 30) (utc 1999 8 18 3 28)
    (utc 2000 12 30 10 5) (utc 2004 7 1 2 48)]
   {:long? [true true false false false] :dsm [false true false false false]}])

(def ^:private voyager-2
  [[:earth :jupiter :saturn :uranus :neptune]
   [(utc 1977 8 20 14 29) (utc 1979 7 9 22 29) (utc 1981 8 26 3 24) (utc 1986 1 24 17 59) (utc 1989 8 25 3 56)]
   {}])

(defn- fly [[bodies dates opts]] (tour/fly bodies dates opts))

(deftest galileo-veega
  (let [{:keys [flybys legs feasible?]} (fly galileo)
        [venus earth-1 earth-2] flybys]
    (is feasible?)
    (testing "nearly ballistic"
      (is (every? #(< (:dv %) 0.15) flybys)))
    (testing "Venus at 16,000 km"
      (is (< (abs (- (:altitude venus) 16000.0)) 800.0)))
    (testing "both Earth passes close, as they were (960 and 303 km)"
      (is (every? #(< (:altitude %) 2000.0) [earth-1 earth-2])))
    (testing "the two years between Earth passes flown as a resonant orbit, whose
              return point the Earth has left by less than its sphere of influence"
      (let [leg (legs 2)]
        (is (:resonant? leg))
        (is (< (:miss leg) (ip/sphere-of-influence :earth)))))))

(deftest cassini-vvejga
  (let [{:keys [flybys legs feasible? dsm-dv]} (fly cassini)
        [venus-1 venus-2 earth jupiter] flybys
        {:keys [t]} (:dsm (legs 1))]
    (is feasible?)
    (testing "the deep space maneuver: close to 450 m/s, in December 1998"
      (is (< 0.40 dsm-dv 0.50))
      (is (< (abs (- t (time/calendar->mjd 1998 12 3))) 30.0)))
    (testing "and with it the flybys nearly ballistic"
      (is (every? #(< (:dv %) 0.2) flybys)))
    (testing "Venus twice under a thousand kilometers (284 and about 600)"
      (is (every? #(< (:altitude %) 1000.0) [venus-1 venus-2])))
    (testing "Earth at about 1,100 km"
      (is (< 300.0 (:altitude earth) 2000.0)))
    (testing "Jupiter at about 10 million km"
      (is (< 7e6 (:altitude jupiter) 13e6)))))

(deftest voyager-2-grand-tour
  (let [{:keys [flybys feasible?]} (fly voyager-2)
        [jupiter _ uranus] flybys]
    (is feasible?)
    (testing "ballistic from Jupiter to Neptune"
      (is (every? #(< (:dv %) 0.1) flybys)))
    (testing "Jupiter 350,000 miles above the cloud tops"
      (is (< (abs (- (:altitude jupiter) 563000.0)) (* 0.15 563000.0))))
    (testing "Uranus 81,500 km above them"
      (is (< (abs (- (:altitude uranus) 81500.0)) 800.0)))))

(deftest optimizing
  (testing "from Voyager's dates put off by days, a cheaper tour that can be flown"
    (let [[bodies dates] voyager-2
          start (mapv + dates [5.0 -10.0 12.0 -8.0 6.0])
          opts {:arrival :none}
          before (tour/cost (tour/fly bodies start {}) opts)
          {:keys [cost feasible?]} (tour/optimize bodies start opts)]
      (is feasible?)
      (is (< cost before)))))
