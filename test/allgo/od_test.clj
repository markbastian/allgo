(ns allgo.od-test
  "Orbit determination against a known orbit observed from a turning
  ground site: the transition matrix against differences of the flown
  trajectory; the sequential batch against the batch taken whole; the
  batch iterated to the orbit itself from exact data and within its
  covariance from noisy; and the extended and unscented Kalman filters
  converging from a start kilometers out, their errors within their
  covariances."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kep]
            [allgo.astro.od :as od]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(defn- two-body [_ r _] (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))))

(def ^:private truth
  (kep/elements->state mu {:a 7200.0 :e 0.01 :i 1.0 :raan 0.5 :argp 0.3 :M 0.0}))

(defn- site
  "A station at 35 degrees north on a spherical Earth, turning with it."
  [t]
  (let [th (+ 0.2 (* 7.292115e-5 t)) lat (math/to-radians 35.0)]
    (v3/scale [(* (math/cos lat) (math/cos th)) (* (math/cos lat) (math/sin th)) (math/sin lat)] c/R-earth)))

(defn- measure
  "Range (km), and topocentric right ascension and declination (rad)."
  [t r _]
  (let [rho (v3/sub r (site t))
        [x y z] rho
        d (v3/length rho)]
    [d (math/atan2 y x) (math/asin (/ z d))]))

(def ^:private times
  "Every 20 seconds of a ten-minute pass over the station an orbit after
  the epoch, the satellite above its horizon throughout."
  (let [[r v] truth
        ts (vec (range 6100.0 6700.0 20.0))]
    (assert (every? (fn [t] (let [[rt] (kep/propagate mu r v t) s (site t)] (pos? (v3/dot (v3/sub rt s) s)))) ts))
    ts))

(def ^:private sigma [0.01 2e-5 2e-5])

(defn- observations
  "The observations at `times`, noise from a seeded generator scaled by
  `noise` standard deviations."
  [noise]
  (let [rng (java.util.Random. 20260928)
        [r v] truth]
    (mapv (fn [t]
            (let [[rt vt] (kep/propagate mu r v t)]
              {:t t :model measure :sigma sigma
               :z (mapv #(+ %1 (* noise %2 (.nextGaussian rng))) (measure t rt vt) sigma)}))
          times)))

(deftest transition-matrix
  (testing "the integrated transition matrix is the flown trajectory's derivative"
    (let [[r v] truth
          x0 (vec (concat r v))
          fly (fn [x] (let [{r1 :r v1 :v} (od/propagate-with-stm two-body 0.0 [(subvec x 0 3) (subvec x 3 6)] 3000.0)]
                        (vec (concat r1 v1))))
          {:keys [phi]} (od/propagate-with-stm two-body 0.0 truth 3000.0)
          numeric (diff/jacobian fly x0 {:steps (vec (concat (repeat 3 1e-3) (repeat 3 1e-6)))})]
      (is (every? true? (for [i (range 6) j (range 6)]
                          (< (abs (- (get-in phi [i j]) (get-in numeric [i j])))
                             (* 1e-6 (max 1.0 (abs (get-in numeric [i j])))))))))))

(deftest sequential-batch
  (testing "the batches taken in turn give what the batch taken whole gives"
    (let [rs (od/batch-rows two-body 0.0 truth (observations 1.0))
          whole (peek (od/sequential-batch nil [rs]))
          split (od/sequential-batch nil (partition-all 17 rs))
          prior {:x [0.1 -0.2 0.05 1e-4 0.0 -2e-4] :P (lin/mat-scale (lin/eye 6) 1.0)}
          whole-p (peek (od/sequential-batch prior [rs]))
          split-p (od/sequential-batch prior (partition-all 17 rs))]
      (is (< (lin/distance (:x whole) (:x (peek split))) 1e-9))
      (is (< (lin/distance (:x whole-p) (:x (peek split-p))) 1e-9))
      (testing "and the uncertainty shrinks with each batch"
        (let [traces (map #(reduce + (map (fn [i] (get-in (:P %) [i i])) (range 3))) split-p)]
          (is (apply > traces)))))))

(deftest batch-least-squares
  (let [start [(v3/add (first truth) [3.0 -2.0 1.0]) (v3/add (second truth) [0.002 0.001 -0.003])]]
    (testing "from exact observations, the orbit itself"
      (let [{:keys [r v]} (od/batch two-body 0.0 start (observations 0.0) {})]
        (is (< (v3/distance r (first truth)) 1e-6))
        (is (< (v3/distance v (second truth)) 1e-9))))
    (testing "from noisy ones, within three standard deviations the covariance gives"
      (let [{:keys [r v P]} (od/batch two-body 0.0 start (observations 1.0) {})
            err (vec (concat (v3/sub r (first truth)) (v3/sub v (second truth))))]
        (doseq [i (range 6)]
          (is (< (abs (err i)) (* 3.0 (math/sqrt (get-in P [i i])))) (str i)))))
    (testing "and taken in batches, the same orbit"
      (let [whole (od/batch two-body 0.0 start (observations 1.0) {})
            split (od/batch two-body 0.0 start (observations 1.0) {:batches 20})]
        (is (< (v3/distance (:r whole) (:r split)) 1e-8))))))

(deftest extended-kalman-filter
  (testing "from a start 5 km and 5 m/s out, noisy observations bring it within its covariance"
    (let [start [(v3/add (first truth) [5.0 -3.0 2.0]) (v3/add (second truth) [0.005 -0.003 0.002])]
          P0 (vec (for [i (range 6)] (vec (for [j (range 6)] (if (= i j) (if (< i 3) 100.0 1e-4) 0.0)))))
          run (od/ekf two-body start P0 0.0 (observations 1.0) {:q 1e-15})
          {:keys [t r v P]} (peek run)
          [rt vt] (kep/propagate mu (first truth) (second truth) t)
          err (vec (concat (v3/sub r rt) (v3/sub v vt)))]
      (is (< (v3/length (subvec err 0 3)) 0.05) "tens of meters")
      (doseq [i (range 6)]
        (is (< (abs (err i)) (* 3.0 (math/sqrt (get-in P [i i])))) (str i)))
      (testing "and its residuals settle to the noise"
        (let [late (map :residuals (drop (quot (count run) 2) run))
              rms (math/sqrt (/ (reduce + (map #(* (first %) (first %)) late)) (count late)))]
          (is (< rms (* 2.0 (first sigma)))))))))

(deftest unscented-kalman-filter
  (testing "the same start and observations as the extended filter: within its covariance, and as close as it"
    (let [start [(v3/add (first truth) [5.0 -3.0 2.0]) (v3/add (second truth) [0.005 -0.003 0.002])]
          P0 (vec (for [i (range 6)] (vec (for [j (range 6)] (if (= i j) (if (< i 3) 100.0 1e-4) 0.0)))))
          obs (observations 1.0)
          u (peek (od/ukf two-body start P0 0.0 obs {:q 1e-15}))
          e (peek (od/ekf two-body start P0 0.0 obs {:q 1e-15}))
          [rt vt] (kep/propagate mu (first truth) (second truth) (:t u))
          err (vec (concat (v3/sub (:r u) rt) (v3/sub (:v u) vt)))]
      (is (< (v3/length (subvec err 0 3)) 0.05) "tens of meters")
      (doseq [i (range 6)]
        (is (< (abs (err i)) (* 3.0 (math/sqrt (get-in (:P u) [i i])))) (str i)))
      (is (< (v3/distance (:r u) (:r e)) 0.02) "the two filters agree to meters"))))
