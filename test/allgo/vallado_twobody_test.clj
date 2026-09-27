(ns allgo.vallado-twobody-test
  "Vallado chapter 2 -- universal variables, the anomalies of every conic,
  time of flight and the J2 drift -- checked against what each must agree
  with: the series the Stumpff functions sum, the elliptic propagator, the
  conics' own Kepler equations and geometry, and orbits whose timing is
  known."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))
(defn- all-close? [as bs tol] (every? true? (map #(close? %1 %2 tol) as bs)))
(defn- same-angle? [a b tol] (< (abs (am/wrap-angle (- a b))) tol))

(defn- series
  "c2 and c3 by their power series: sum (-psi)^k / (2k+2)! and / (2k+3)!."
  [psi]
  (let [terms (fn [start] (loop [k 0 term (/ 1.0 (reduce * (range 1 (inc start)))) sum 0.0]
                            (if (or (> k 80) (and (> k 5) (< (abs term) 1e-18)))
                              sum
                              (recur (inc k)
                                     (/ (* term (- psi)) (* (+ start (* 2 k) 1) (+ start (* 2 k) 2)))
                                     (+ sum term)))))]
    [(terms 2) (terms 3)]))

(deftest stumpff-functions
  (doseq [psi [-39.5 -4.0 -1e-8 0.0 1e-8 0.57 4.0 39.5]]
    (let [[a b] (u/stumpff psi)
          [c2 c3] (series psi)]
      (is (close? a c2 1e-10) (str psi))
      (is (close? b c3 1e-10) (str psi)))))

(deftest universal-kepler
  (testing "agrees with the elliptic propagator, and runs backward"
    (doseq [[r0 v0 dt] [[[7000.0 1000.0 -500.0] [0.5 7.2 1.1] 4000.0]
                        [[-6518.0 -2404.0 -22.0] [2.6 -7.1 -0.26] 120.0]
                        [[7000.0 0.0 0.0] [0.0 9.5 3.0] 40000.0]]]
      (let [[r v] (u/propagate mu r0 v0 dt)
            [r' v'] (kep/propagate mu r0 v0 dt)
            [rb vb] (u/propagate mu r v (- dt))]
        (is (all-close? r r' 1e-9))
        (is (all-close? v v' 1e-9))
        (is (all-close? rb r0 1e-9))
        (is (all-close? vb v0 1e-9)))))
  (testing "a hyperbola: where the hyperbolic Kepler equation puts it"
    (let [r0 [7000.0 0.0 0.0] v0 [0.0 12.0 1.0] dt 20000.0
          [r v] (u/propagate mu r0 v0 dt)
          h (v3/cross r0 v0)
          a (/ 1.0 (- (/ 2.0 (v3/length r0)) (/ (v3/length-squared v0) mu)))
          e (math/sqrt (- 1.0 (/ (v3/length-squared h) (* mu a))))
          n (math/sqrt (/ mu (- (* a a a))))
          ;; starting at periapsis: M = 0
          H (kep/kepler-hyperbolic (* n dt) e)]
      (is (close? (v3/length r) (* a (- 1.0 (* e (math/cosh H)))) 1e-10))
      (is (close? (kep/specific-energy mu r0 v0) (kep/specific-energy mu r v) 1e-10))
      (is (all-close? h (v3/cross r v) 1e-10))))
  (testing "near a parabola, a billion kilometers out: conserved, and reversible"
    (let [r0 [15912.0 1.352e9 0.0] v0 [-0.000784 31889.0 0.0]
          [r v] (u/propagate mu r0 v0 100000.0)
          [rb vb] (u/propagate mu r v -100000.0)]
      (is (all-close? (v3/cross r0 v0) (v3/cross r v) 1e-10))
      (is (all-close? rb r0 1e-8))
      (is (all-close? vb v0 1e-8))))
  (testing "the Taylor series matches for a short step"
    (let [r0 [-6518.0 -2404.0 -22.0] v0 [2.6 -7.1 -0.26]
          dt 60.0
          uv (u/chi mu r0 v0 dt)]
      (is (all-close? (u/f-and-g mu (v3/length r0) dt uv)
                      (u/f-and-g-series mu r0 v0 dt) 1e-9)))))

(deftest anomalies-of-every-conic
  (testing "each anomaly, its mean anomaly and its true anomaly round-trip"
    (doseq [[e x] [[0.0 0.74] [0.4 5.84] [0.9 2.0] [1.05 2.15] [1.5 -0.8] [3.0 1.2] [1.0 0.58] [1.0 -1.3]]]
      (let [[M nu] (kep/anomaly->mean-and-true e x)
            [x' nu'] (kep/mean->anomaly-and-true e M)
            [x'' M''] (kep/true->anomaly-and-mean e nu)]
        (is (close? x' (if (< e 1.0) (am/wrap-2pi x) x) 1e-10) (str e))
        (is (same-angle? nu' nu 1e-10) (str e))
        (is (same-angle? x'' x 1e-10) (str e))
        (is (same-angle? M'' M 1e-10) (str e)))))
  (testing "and put the body at the same distance: r = p / (1 + e cos nu)"
    (let [p 10000.0]
      (doseq [[e x] [[0.4 5.84] [0.9 2.0] [1.05 2.15] [1.5 -0.8] [1.0 0.58]]]
        (let [[_ nu] (kep/anomaly->mean-and-true e x)
              r (/ p (+ 1.0 (* e (math/cos nu))))
              r' (cond (< e 1.0) (* (/ p (- 1.0 (* e e))) (- 1.0 (* e (math/cos x))))
                       (> e 1.0) (* (/ p (- 1.0 (* e e))) (- 1.0 (* e (math/cosh x))))
                       :else (* 0.5 p (+ 1.0 (* x x))))]
          (is (close? r r' 1e-10) (str e))))))
  (testing "past periapsis, the hyperbolic anomaly, mean and true anomalies share a sign"
    (let [[M nu] (kep/anomaly->mean-and-true 1.05 2.15)]
      (is (pos? M))
      (is (pos? nu))))
  (testing "circular: all three are the same angle"
    (let [[M nu] (kep/anomaly->mean-and-true 0.0 0.74)]
      (is (close? M 0.74 1e-15))
      (is (close? nu 0.74 1e-15)))))

(deftest time-of-flight
  (testing "between two points of orbits whose timing is known"
    (let [at (fn [p e nu] (let [r (/ p (+ 1.0 (* e (math/cos nu))))]
                            [(* r (math/cos nu)) (* r (math/sin nu)) 0.0]))
          elapsed (fn [p e nu1 nu2]
                    (cond
                      (< e 1.0) (let [a (/ p (- 1.0 (* e e)))
                                      n (math/sqrt (/ mu (* a a a)))]
                                  (/ (- (kep/true->mean nu2 e) (kep/true->mean nu1 e)) n))
                      (> e 1.0) (let [a (/ p (- 1.0 (* e e)))
                                      n (math/sqrt (/ mu (- (* a a a))))
                                      m (fn [nu] (second (kep/true->anomaly-and-mean e nu)))]
                                  (/ (- (m nu2) (m nu1)) n))
                      :else (let [m (fn [nu] (second (kep/true->anomaly-and-mean 1.0 nu)))]
                              (* 0.5 (math/sqrt (/ (* p p p) mu)) (- (m nu2) (m nu1))))))]
      (doseq [[p e nu1 nu2] [[9000.0 0.3 0.2 1.9] [8000.0 0.0 0.0 1.0] [12000.0 1.4 -0.5 0.9]
                             [14000.0 1.0 0.1 1.2]]]
        (is (close? (u/time-of-flight mu (at p e nu1) (at p e nu2) p) (elapsed p e nu1 nu2) 1e-8)
            (str [p e])))))
  (is (NaN? (u/time-of-flight mu [7000.0 0.0 0.0] [14000.0 0.0 0.0] 14000.0))
      "a = 0: no such orbit"))

(deftest hitting-the-earth
  (is (= :at-an-end (u/hits-earth? [6378.0 0.0 0.0] [0.0 7.8 0.0] [7000.0 0.0 0.0] [0.0 7.0 0.0] 0 100.0)))
  (is (= :during-revolutions (u/hits-earth? [7000.0 0.0 0.0] [0.0 7.0 0.0] [8000.0 0.0 0.0] [0.0 6.5 0.0] 1 100.0)))
  (is (nil? (u/hits-earth? [8000.0 0.0 0.0] [0.0 7.0 0.0] [9000.0 0.0 0.0] [0.0 6.5 0.0] 0 100.0))))

(deftest j2-drift
  (testing "the node and perigee move at J2's averaged rates"
    (let [r0 [-6518.1083 -2403.8479 -22.1722] v0 [2.604057 -7.105717 -0.263218]
          dt 86400.0
          el0 (kep/state->elements c/GM-earth r0 v0)
          [r v] (u/propagate-j2 r0 v0 dt)
          el1 (kep/state->elements c/GM-earth r v)
          {:keys [a e i]} el0
          n (kep/mean-motion c/GM-earth a)
          p (* a (- 1 (* e e)))
          k (/ (* 1.5 n geo/J2 c/R-earth c/R-earth) (* p p))]
      (is (same-angle? (- (:raan el1) (:raan el0)) (* (- k) (math/cos i) dt) 1e-9))
      (is (same-angle? (- (:argp el1) (:argp el0)) (* k (- 2 (* 2.5 (math/sin i) (math/sin i))) dt) 1e-7))
      (is (close? (:a el1) a 1e-9) "no drag, no decay"))))
