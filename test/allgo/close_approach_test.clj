(ns allgo.close-approach-test
  "ANCAS against the truth found the slow way: two satellites set up to
  pass 150 m apart at a known moment, and every close approach of a pair
  of crossing orbits over a day checked against a dense search of the
  true distance."
  (:require [allgo.astro.close-approach :as ca]
            [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kep]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.roots :as roots]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)

(defn- ephemeris [[r v]] (fn [t] (u/propagate mu r v t)))

(def ^:private a-state (kep/elements->state mu {:a 7000.0 :e 0.001 :i 1.7 :raan 0.2 :argp 0.0 :M 0.0}))

(deftest a-set-up-encounter
  (testing "B passes 150 m from A an hour and a quarter in, crossing at 14 km/s"
    (let [tca 4500.0
          [ra va] (u/propagate mu (first a-state) (second a-state) tca)
          ;; B at tca: 150 m off A, normal to both velocities, on a
          ;; retrograde-ish orbit through there
          vb (v3/scale (v3/normalize (v3/add (v3/scale va -0.2) (v3/cross ra va))) 7.5)
          off (v3/scale (v3/normalize (v3/cross va vb)) 0.150)
          [rb0 vb0] (u/propagate mu (v3/add ra off) vb (- tca))
          eph-a (ephemeris a-state) eph-b (ephemeris [rb0 vb0])
          found (ca/close-approaches eph-a eph-b 0.0 7200.0 120.0 {:within 10.0})]
      (is (= 1 (count found)))
      (let [{:keys [t miss t-fit miss-fit]} (first found)]
        (is (< (abs (- t tca)) 1e-4))
        (is (< (abs (- miss 0.150)) 1e-4))
        (testing "and the fitted cubics alone already close"
          (is (< (abs (- t-fit tca)) 0.5))
          (is (< (abs (- miss-fit 0.150)) 0.5)))))))

(deftest a-day-of-passes
  (testing "every local minimum of the true distance over a day, found and matched"
    (let [b-state (kep/elements->state mu {:a 7010.0 :e 0.002 :i 1.2 :raan 0.9 :argp 1.0 :M 2.0})
          eph-a (ephemeris a-state) eph-b (ephemeris b-state)
          dist (fn [t] (v3/distance (first (eph-a t)) (first (eph-b t))))
          ;; the slow way: every second, then golden section
          ts (range 0.0 86400.0 1.0)
          truth (for [[x y z] (partition 3 1 ts)
                      :when (and (< (dist y) (dist x)) (<= (dist y) (dist z)))]
                  (roots/minimize dist x z {:tol 1e-6}))
          found (ca/close-approaches eph-a eph-b 0.0 86400.0 300.0)]
      (is (= (count truth) (count found)))
      (doseq [[t {ft :t miss :miss}] (map vector truth found)]
        (is (< (abs (- t ft)) 1e-3))
        (is (< (abs (- (dist t) miss)) 1e-6))))))
