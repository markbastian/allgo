(ns allgo.vallado-iod-test
  "Vallado chapter 7 -- Gibbs, Herrick-Gibbs, Gauss and Lambert -- each
  given positions, sightings or times taken from an orbit flown with the
  two-body propagator, and asked to recover it; every Lambert solution is
  flown as well, to check it arrives."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.iod :as iod]
            [allgo.astro.kepler :as kep]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))
(defn- rel-err [a b] (/ (v3/distance a b) (v3/length b)))

(def ^:private truth
  "An orbit to recover: 12000 km, e 0.2, inclined 30 degrees."
  (kep/elements->state mu {:a 12000.0 :e 0.2 :i (math/to-radians 30) :raan (math/to-radians 40)
                           :argp (math/to-radians 60) :M 0.3}))

(defn- at [t] (let [[r v] truth] (kep/propagate mu r v t)))
(defn- angle [a b] (math/acos (/ (v3/dot a b) (* (v3/length a) (v3/length b)))))

(deftest three-positions
  (testing "Gibbs, from three positions tens of degrees apart"
    (let [[r1] (at -1800.0) [r2 v2] (at 0.0) [r3] (at 2400.0)
          {:keys [theta12 theta23 copa] v :v2} (iod/gibbs r1 r2 r3)]
      (is (< (rel-err v v2) 1e-10))
      (is (close? theta12 (angle r1 r2) 1e-12))
      (is (close? theta23 (angle r2 r3) 1e-12))
      (is (< (abs copa) 1e-12))))
  (testing "Gibbs sees a position off the plane"
    (let [[r1] (at -1800.0) [r2] (at 0.0) [r3] (at 2400.0)]
      (is (> (abs (:copa (iod/gibbs (v3/add r1 [0.0 0.0 50.0]) r2 r3))) 1e-3))))
  (testing "Herrick-Gibbs, from three positions a minute apart"
    (let [[r1] (at -60.0) [r2 v2] (at 0.0) [r3] (at 75.0)
          {v :v2} (iod/herrick-gibbs r1 r2 r3 -60.0 0.0 75.0)]
      ;; a Taylor expansion, good to its truncation
      (is (< (rel-err v v2) 1e-6)))))

(def ^:private omega 7.292115e-5)

(defn- site
  "A site at latitude `lat` on a spherical Earth, turned by the Earth's
  rotation `t` seconds from when it was on the x-z plane."
  [lat t]
  (let [th (* omega t)]
    (v3/scale [(* (math/cos lat) (math/cos th)) (* (math/cos lat) (math/sin th)) (math/sin lat)]
              c/R-earth)))

(deftest angles-only
  (testing "Gauss, from three sightings a few minutes apart"
    (let [ts [-480.0 0.0 480.0]
          sites (map #(site (math/to-radians 40) %) ts)
          obs (map (fn [t s]
                     (let [[x y z] (v3/normalize (v3/sub (first (at t)) s))]
                       [(math/atan2 y x) (math/asin z)]))
                   ts sites)
          {:keys [r2 v2]} (iod/gauss obs ts sites)
          [r-true v-true] (at 0.0)]
      ;; Gauss's method truncates the f and g series, so it is close, not exact
      (is (< (rel-err r2 r-true) 0.02))
      (is (< (rel-err v2 v-true) 0.05)))))

(defn- sightings
  "Sightings of the true orbit at `ts` from a site at 40 degrees north."
  [ts]
  (let [sites (map #(site (math/to-radians 40) %) ts)]
    [(map (fn [t s]
            (let [[x y z] (v3/normalize (v3/sub (first (at t)) s))]
              [(math/atan2 y x) (math/asin z)]))
          ts sites)
     sites]))

(deftest double-r
  (testing "refines Gauss to the orbit itself"
    (doseq [ts [[-480.0 0.0 480.0] [-900.0 0.0 1200.0]]]
      (let [[obs sites] (sightings ts)
            {:keys [r2 v2]} (iod/double-r obs ts sites)
            gauss (iod/gauss obs ts sites)
            [r-true v-true] (at 0.0)]
        (is (< (rel-err r2 r-true) 1e-9) (str ts))
        (is (< (rel-err v2 v-true) 1e-9) (str ts))
        (is (> (rel-err (:r2 gauss) r-true) 1e-4) "where Gauss alone was off"))))
  (testing "from a guess well away from the answer"
    (let [ts [-480.0 0.0 480.0]
          [obs sites] (sightings ts)
          [r-true] (at 0.0)
          {:keys [r2]} (iod/double-r c/GM-earth obs ts sites [9000.0 9000.0])]
      (is (< (rel-err r2 r-true) 1e-9)))))

(deftest laplace
  (testing "Laplace's method, from three sightings a few minutes apart"
    (let [ts [-480.0 0.0 480.0]
          [obs sites] (sightings ts)
          {:keys [r2 v2]} (iod/laplace obs ts sites)
          [r-true v-true] (at 0.0)]
      ;; the derivatives are a quadratic's, so it is close, not exact:
      ;; 0.4% in position here, where Gauss's single pass is 1.1% off
      (is (< (rel-err r2 r-true) 0.01))
      (is (< (rel-err v2 v-true) 0.03))))
  (testing "and closer as the sightings close up, as its derivatives do"
    (let [err (fn [dt] (let [ts [(- dt) 0.0 dt] [obs sites] (sightings ts)]
                         (rel-err (:r2 (iod/laplace obs ts sites)) (first (at 0.0)))))]
      ;; a quarter the error for half the arc
      (is (< (* 3.0 (err 120.0)) (err 240.0)))
      (is (< (* 3.0 (err 240.0)) (err 480.0))))))

(defn- arrives? [r1 v1 r2 dt]
  (let [[r _] (u/propagate mu r1 v1 dt)]
    (< (v3/distance r r2) (* 1e-6 (v3/length r2)))))

(deftest lambert
  (let [[r1 v1] (at 0.0)
        period (kep/period mu 12000.0)]
    (testing "the short way: the velocities of the orbit the two points came from"
      (let [dt 2400.0 [r2 v2] (at dt)
            [a b] (iod/lambert r1 r2 dt {})]
        (is (< (rel-err a v1) 1e-8))
        (is (< (rel-err b v2) 1e-8))))
    (testing "the long way"
      (let [dt (* 0.7 period) [r2 v2] (at dt)
            [a b] (iod/lambert r1 r2 dt {:long? true})]
        (is (< (rel-err a v1) 1e-8))
        (is (< (rel-err b v2) 1e-8))))
    (testing "with a revolution on the way, one of the two orbits is the true one, and both arrive"
      (let [dt (+ period 2400.0) [r2 v2] (at dt)
            sols (for [high? [false true]] (iod/lambert r1 r2 dt {:revs 1 :high? high?}))]
        (is (some (fn [[a b]] (and (< (rel-err a v1) 1e-7) (< (rel-err b v2) 1e-7))) sols))
        (doseq [[a] sols] (is (arrives? r1 a r2 dt)))))
    (testing "all four one-revolution transfers arrive"
      (let [r2 (first (at 3000.0)) dt (* 1.6 period)]
        (doseq [long? [false true] high? [false true]
                :let [[a] (iod/lambert r1 r2 dt {:long? long? :revs 1 :high? high?})]]
          (is (arrives? r1 a r2 dt) (str {:long? long? :high? high?})))))
    (testing "too little time for a revolution"
      (is (nil? (iod/lambert r1 (first (at 3000.0)) 4560.0 {:revs 1}))))))

(deftest lambert-limits
  (let [r1 (first (at 0.0)) r2 (first (at 3000.0))
        chord (v3/distance r1 r2)
        s (* 0.5 (+ (v3/length r1) (v3/length r2) chord))]
    (testing "the minimum-energy orbit's semi-major axis is half the semiperimeter"
      (is (close? (:a (iod/lambert-min-energy r1 r2 false 0)) (* 0.5 s) 1e-12)))
    (testing "flying it for its time arrives, the short way and the long"
      (doseq [long? [false true]]
        (let [{:keys [v1 tof]} (iod/lambert-min-energy r1 r2 long? 0)]
          (is (arrives? r1 v1 r2 tof) (str long?)))))
    (testing "and Lambert, asked for that time, finds that orbit"
      (let [{:keys [v1 tof]} (iod/lambert-min-energy r1 r2 false 0)
            [a] (iod/lambert r1 r2 tof {})]
        (is (< (rel-err a v1) 1e-7))))
    (testing "the long way with a revolution takes longer than the short"
      (is (< (:tof (iod/lambert-min-energy r1 r2 false 1)) (:tof (iod/lambert-min-energy r1 r2 true 1)))))
    (testing "no ellipse beats the parabola, and just above its time the orbit is nearly one"
      (let [{:keys [tof tof-parabolic]} (iod/lambert-min-energy r1 r2 false 0)
            [a] (iod/lambert r1 r2 (* tof-parabolic 1.0001) {})
            energy (kep/specific-energy mu r1 a)]
        (is (< tof-parabolic tof))
        (is (< (abs energy) (* 1e-3 (/ mu (v3/length r1)))))))))
