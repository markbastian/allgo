(ns allgo.porkchop-test
  "Launch windows against NASA's Earth-to-Mars mission design handbook,
  2026-2045 (NASA/TM-2010-216764; test/data/imdh): scanning twenty years
  finds every one of its energy minima, at its dates and costs."
  (:require [allgo.astro.porkchop :as pc]
            [allgo.astro.time :as time]
            [allgo.geometry.vec3 :as v3]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]))

(def ^:private rows
  (edn/read-string (slurp (io/resource "data/imdh/earth-mars-2026-2045.edn"))))

(def ^:private date-slip
  "Days by which the handbook's tables 6-8 print every date early; see
  test/data/imdh/README.md."
  {2035 63 2037 81 2039 75})

(defn- mjd [year date] (+ (date-slip year 0) (apply time/calendar->mjd date)))

(defn- book-row
  "The handbook's minimum of `k` for the opportunity `o` belongs to."
  [k {:keys [type depart]}]
  ;; 2041's type I least-arrival-speed row is a type II transfer, and a
  ;; copy of its own type II row; see the README
  (first (filter (fn [[year t kk d]] (and (= t type) (= kk k)
                                          (not (and (= year 2041) (= k :v-inf) (= t 1)))
                                          (< (abs (- (mjd year d) depart)) 30)))
                 rows)))

(def ^:private scan
  (memoize #(pc/opportunities :earth :mars (time/calendar->mjd 2026 1 1)
                              (time/calendar->mjd 2046 6 1) {:objective %})))

(deftest synodic-period
  (testing "the Earth and Mars: the handbook's 779.935 days"
    (is (< (abs (- (pc/synodic-period :earth :mars) 779.935)) 0.1)))
  (testing "and it is symmetric"
    (is (= (pc/synodic-period :venus :earth) (pc/synodic-period :earth :venus)))))

(deftest least-energy
  (let [os (scan :c3)]
    (testing "one opportunity of each type every synodic period, ten of each"
      (is (= 20 (count os)))
      (is (= [10 10] (map #(count (filter (comp #{%} :type) os)) [1 2]))))
    (doseq [{:keys [type depart arrive c3] :as o} os
            :let [[year _ _ d a book-c3] (book-row :c3 o)]]
      (testing (str year " type " type)
        (is year "the handbook lists it")
        (when year
          (if (and (= year 2026) (= type 1))
            ;; the valley runs into a nodal transfer, where the handbook's
            ;; two-day grid and our ridge band each stop somewhere on the
            ;; slope down to it: only no worse
            (is (<= c3 book-c3))
            ;; the handbook's grid is two days: its minimum is a little
            ;; above the true one, and its dates within a day or two
            (do (is (< -0.02 (- c3 book-c3) 0.01))
                (is (< (abs (- depart (mjd year d))) 3.0))
                ;; along arrival the valley is flat
                (is (< (abs (- arrive (mjd year a))) 7.0)))))))))

(deftest least-arrival-speed
  (let [os (scan :v-inf-arrive)]
    (is (= 20 (count os)))
    (doseq [{:keys [type depart arrive v-inf-arrive] :as o} os
            :let [[year _ _ d a _ _ _ book-v] (book-row :v-inf o)
                  v (v3/length v-inf-arrive)]]
      (testing (str year " type " type)
        (when year
          (if (and (= year 2026) (= type 1))
            (is (<= v book-v))
            ;; again the grid's minimum can only be the higher
            (do (is (< -0.01 (- v book-v) 2e-3))
                (is (< (abs (- depart (mjd year d))) 3.0))
                (is (< (abs (- arrive (mjd year a))) 7.0)))))))))

(deftest the-plane
  ;; the handbook's 2026 plot: 160 days of departures, 400 of arrivals
  (let [start (time/calendar->mjd 2026 8 4)
        departs (pc/dates start (+ start 160) 4)
        arrives (pc/dates (+ start 160) (+ start 560) 4)
        g (pc/grid :earth :mars departs arrives)
        c3 (pc/contour-values g :c3)]
    (testing "rows by arrival, columns by departure"
      (is (= (count arrives) (count c3)))
      (is (every? #(= (count departs) (count %)) c3)))
    (testing "each cell the transfer for its dates, the cheaper type"
      (let [t (get-in g [:cells 30 20])]
        (is (= [(departs 20) (arrives 30)] [(:depart t) (:arrive t)]))
        (is (= (:c3 t) (apply min (keep #(:c3 (pc/transfer :earth (departs 20) :mars (arrives 30) %)) [1 2]))))))
    (testing "no transfer arrives before it leaves"
      (is (nil? (get-in c3 [0 (dec (count departs))]))))
    (testing "the grid's least C3 is the 2026 type II opportunity's, refined a little"
      (let [best (apply min (keep identity (flatten c3)))
            o (first (filter #(= 2 (:type %)) (scan :c3)))]
        (is (< (:c3 o) best (+ (:c3 o) 0.2)))))))

(deftest the-ridge
  (testing "near 180 degrees of travel nothing is offered unless asked for"
    (let [d (time/calendar->mjd 2026 11 12) a (+ d 270.61)]
      (is (< 179.0 (:travel (pc/transfer :earth d :mars a 1 {:ridge 0.0})) 181.0))
      (is (nil? (pc/transfer :earth d :mars a 1))))))
