(ns allgo.range-iod-test
  "Initial orbits from ranges and range rates (Vallado, section 7.4),
  against a known orbit seen at once from several turning ground
  stations, placed about its ground track so that they see it."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.iod :as iod]
            [allgo.astro.kepler :as kep]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)

(def ^:private truth
  (kep/elements->state mu {:a 12000.0 :e 0.1 :i 0.9 :raan 0.4 :argp 0.3 :M 0.2}))

(defn- at [t] (kep/propagate mu (first truth) (second truth) t))

(defn- station
  "A station at latitude and longitude (degrees) on a spherical Earth
  turning from 0 at t = 0: its inertial position and velocity."
  [lat lon t]
  (let [la (math/to-radians lat) th (+ (math/to-radians lon) (* c/omega-earth t))
        r (v3/scale [(* (math/cos la) (math/cos th)) (* (math/cos la) (math/sin th)) (math/sin la)] c/R-earth)]
    [r (v3/cross [0.0 0.0 c/omega-earth] r)]))

(def ^:private places
  "Stations around the satellite's ground point at t = 1000 s, which see
  it overhead-ish then and through the twenty minutes either side."
  (let [[[x y z]] (at 1000.0)
        lat (math/to-degrees (math/asin (/ z (v3/length [x y z]))))
        lon (- (math/to-degrees (math/atan2 y x)) (math/to-degrees (* c/omega-earth 1000.0)))]
    (mapv (fn [[dl dn]] [(+ lat dl) (+ lon dn)]) [[10.0 -10.0] [-8.0 12.0] [3.0 15.0] [-12.0 -9.0]])))

(defn- seen
  "At time `t`, the sites' states and the exact ranges and range rates."
  [t n]
  (let [[r v] (at t)
        sts (map (fn [[la lo]] (station la lo t)) (take n places))]
    {:sites (mapv first sts) :site-v (mapv second sts)
     :ranges (mapv (fn [[s]] (v3/distance r s)) sts)
     :rates (mapv (fn [[s sv]] (let [d (v3/sub r s)] (/ (v3/dot d (v3/sub v sv)) (v3/length d)))) sts)}))

(defn- close? [a b tol] (every? #(< (abs %) tol) (map - a b)))

(deftest trilateration
  (let [[r] (at 1000.0)]
    (testing "three stations: the point above their plane is the satellite"
      (let [{:keys [sites ranges]} (seen 1000.0 3)
            [above below] (iod/trilaterate sites ranges)]
        (is (close? above r 1e-6))
        (is (not (close? below r 1.0)))))
    (testing "four: the one point that fits them all"
      (let [{:keys [sites ranges]} (seen 1000.0 4)]
        (is (close? (iod/trilaterate sites ranges) r 1e-8))))))

(deftest range-rates
  (testing "with the position, three range rates give the velocity"
    (let [[r v] (at 1000.0)
          {:keys [sites site-v rates]} (seen 1000.0 3)]
      (is (close? (iod/range-rate-velocity r sites site-v rates) v 1e-9)))))

(deftest ranges-alone
  (doseq [[label dt tol] [["a minute apart, by Herrick-Gibbs" 60.0 1e-6]
                          ["twenty minutes apart, by Gibbs" 1200.0 1e-6]]]
    (testing label
      (let [ts [(- 1000.0 dt) 1000.0 (+ 1000.0 dt)]
            obs (mapv (fn [t] (let [{:keys [sites ranges]} (seen t 3)] [t sites ranges])) ts)
            [r2 v2] (iod/range-only obs)
            [r v] (at 1000.0)]
        (is (close? r2 r 1e-6))
        (is (close? v2 v tol))))))
