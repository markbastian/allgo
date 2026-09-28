(ns allgo.constellation-test
  "Walker patterns and coverage against their own geometry: the pattern's
  spacing as defined, the coverage angle's edge at the minimum
  elevation, one satellite covering its spherical cap, the mean in view
  the sum of the caps, and a GPS-like 24/6/1 covering the globe."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.constellation :as con]
            [allgo.astro.kepler :as kep]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest walker-pattern
  (testing "24/6/1: nodes 60 degrees apart, satellites 90 apart in each plane, 15 degrees phase between planes"
    (let [els (con/walker 26560.0 (math/to-radians 55.0) 24 6 1)
          by-plane (group-by :raan els)]
      (is (= 24 (count els)))
      (is (= 6 (count by-plane)))
      (is (every? #(= 4 (count %)) (vals by-plane)))
      (let [nodes (sort (keys by-plane))]
        (is (every? #(< (abs (- % (/ math/PI 3))) 1e-12) (map - (rest nodes) nodes))))
      (let [first-of (fn [raan] (apply min (map :M (by-plane raan))))
            nodes (sort (keys by-plane))]
        (is (< (abs (- (am/wrap-angle (- (first-of (second nodes)) (first-of (first nodes)))) (math/to-radians 15.0))) 1e-12))))))

(deftest coverage-angle
  (testing "a ground point at the coverage angle sees the satellite at exactly the minimum elevation"
    (doseq [[r eps] [[7000.0 0.0] [7000.0 0.2] [26560.0 (math/to-radians 10.0)]]]
      (let [lam (con/coverage-angle r eps)
            sat [r 0.0 0.0]
            ground (v3/scale [(math/cos lam) (math/sin lam) 0.0] c/R-earth)
            to-sat (v3/sub sat ground)
            elevation (- (/ math/PI 2) (math/acos (/ (v3/dot to-sat ground) (* (v3/length to-sat) c/R-earth))))]
        (is (< (abs (- elevation eps)) 1e-12) (str r " " eps))))))

(deftest coverage
  (let [eps (math/to-radians 10.0)]
    (testing "one satellite covers its spherical cap, (1 - cos lambda)/2 of the globe"
      (let [r 12000.0
            lam (con/coverage-angle r eps)
            {:keys [fraction]} (con/coverage [[0.0 0.0 r]] eps {:bands 300})]
        (is (< (abs (- fraction (* 0.5 (- 1.0 (math/cos lam))))) 2e-3))))
    (testing "the mean number in view is the sum of the caps"
      (let [els (con/walker 20000.0 (math/to-radians 50.0) 12 4 2)
            sats (con/positions els 1234.0)
            lam (con/coverage-angle 20000.0 eps)
            {:keys [mean]} (con/coverage sats eps {:bands 200})]
        (is (< (abs (- mean (* 12 0.5 (- 1.0 (math/cos lam))))) 5e-3))))
    (testing "a GPS-like 24/6/1 at 55 degrees sees everywhere, from four satellites or more, at any moment"
      (let [els (con/walker 26560.0 (math/to-radians 55.0) 24 6 1)]
        (doseq [dt [0.0 3000.0 20000.0]]
          (let [{:keys [fraction min]} (con/coverage (con/positions els dt) (math/to-radians 5.0))]
            (is (< (abs (- fraction 1.0)) 1e-12) (str dt))
            (is (>= min 4) (str dt))))))
    (testing "and a sparse pattern leaves gaps, which the gap count finds"
      (let [els (con/walker 7500.0 (math/to-radians 60.0) 6 3 1)
            {:keys [covered longest]} (con/gaps els [1.0 0.0 0.0] eps (* 2 (kep/period c/GM-earth 7500.0)) 30.0)]
        (is (< 0.0 covered 1.0))
        (is (pos? longest))))))
