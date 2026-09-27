(ns allgo.three-body-test
  (:require [allgo.astro.three-body :as tb]
            [allgo.geometry.vec3 :as v3]
            [clojure.test :refer [deftest is testing]]))

(defn- run [id method t]
  (first (tb/advance (tb/integrator (tb/by-id id) method) t 1000000)))

(defn- drift
  "Relative energy error of an integrator over the run."
  [id s]
  (let [{:keys [masses pos vel]} (tb/by-id id)
        e0 (tb/energy masses pos vel)]
    (abs (/ (- (tb/energy masses (tb/positions s) (tb/velocities s)) e0) e0))))

(defn- return-error
  "How far the bodies are from where they started."
  [id s]
  (reduce max (map v3/distance (tb/bodies (tb/positions s)) (tb/bodies (:pos (tb/by-id id))))))

(deftest acceleration-test
  (testing "equal and opposite: the pulls on the system sum to zero"
    (let [{:keys [masses pos]} (tb/by-id :pythagorean)
          acc (tb/bodies ((tb/acceleration masses) 0.0 pos))]
      (is (< (v3/length (reduce v3/add (map v3/scale acc masses))) 1e-12))))

  (testing "two bodies alone feel m / r^2"
    (let [[ax] ((tb/acceleration [1.0 2.0]) 0.0 [0.0 0.0 0.0 2.0 0.0 0.0])]
      (is (< (abs (- ax 0.5)) 1e-12)))))

(deftest scenarios-test
  (testing "every start is at rest as a whole, at the origin"
    (doseq [{:keys [id masses pos vel]} tb/scenarios]
      (is (< (v3/length (tb/momentum masses vel)) 1e-9) (name id))
      (is (< (v3/length (tb/center-of-mass masses pos)) 1e-7) (name id))))

  (testing "free fall is repeatable by seed and different across seeds"
    (is (= (tb/free-fall 3) (tb/free-fall 3)))
    (is (not= (:pos (tb/free-fall 3)) (:pos (tb/free-fall 4))))))

(deftest periodic-orbits-test
  (testing "the figure eight closes after one period, by every adaptive method"
    (doseq [m ["DOPRI5(4)" "RKF4(5)" "RKN4" "GBS8-2"]]
      (let [s (run :figure-eight m (:period (tb/by-id :figure-eight)))]
        (is (< (return-error :figure-eight s) 1e-5) m)
        (is (< (drift :figure-eight s) 1e-7) m))))

  (testing "Butterfly I closes, to the six digits its start is known to"
    (let [s (run :butterfly "DOPRI5(4)" (:period (tb/by-id :butterfly)))]
      (is (< (return-error :butterfly s) 5e-3))))

  (testing "the Lagrange triangle turns rigidly: a third of a period is a
            third of a turn, which swaps the corners round"
    (let [{:keys [pos period]} (tb/by-id :lagrange)
          [a b c] (tb/bodies (tb/positions (run :lagrange "DOPRI5(4)" (/ period 3.0))))
          [a0 b0 c0] (tb/bodies pos)]
      (is (< (v3/distance a b0) 1e-8))
      (is (< (v3/distance b c0) 1e-8))
      (is (< (v3/distance c a0) 1e-8)))))

(deftest conservation-test
  (testing "momentum and angular momentum hold through close encounters"
    (let [{:keys [masses pos vel]} (tb/by-id :pythagorean)
          s (run :pythagorean "DOPRI5(4)" 30.0)]
      (is (< (v3/length (tb/momentum masses (tb/velocities s))) 1e-8))
      (is (< (v3/distance (tb/angular-momentum masses pos vel)
                          (tb/angular-momentum masses (tb/positions s) (tb/velocities s)))
             1e-6))
      (is (< (drift :pythagorean s) 1e-6)))))

(deftest pythagorean-test
  (testing "ends as Szebehely and Peters found: the lightest body thrown
            out, the other two bound to each other"
    (let [s (run :pythagorean "DOPRI5(4)" 100.0)
          [a b c] (tb/bodies (tb/positions s))
          [_ vb vc] (tb/bodies (tb/velocities s))
          pair-com (v3/scale (v3/add (v3/scale b 4.0) (v3/scale c 5.0)) (/ 1.0 9.0))
          pair-energy (- (* 0.5 (/ 20.0 9.0) (v3/length-squared (v3/sub vb vc)))
                         (/ 20.0 (v3/distance b c)))]
      (is (> (v3/distance a pair-com) 20.0))
      (is (neg? pair-energy)))))

(deftest advance-test
  (testing "lands exactly on the time asked for, and records each step"
    (let [[s path] (tb/advance (tb/integrator (tb/by-id :figure-eight) "DOPRI5(4)") 0.37 10000)]
      (is (== 0.37 (:t s)))
      (is (seq path))
      (is (= (tb/positions s) (peek path)))))

  (testing "stops when its step budget runs out"
    (let [[s path] (tb/advance (tb/integrator (tb/by-id :figure-eight) "Verlet") 1.0 10)]
      (is (= 10 (count path)))
      (is (< (:t s) 1.0))))

  (testing "a nudge moves only the first body"
    (let [sc (tb/by-id :figure-eight)]
      (is (= (subvec (:pos sc) 1) (subvec (:pos (tb/nudge sc 1e-6)) 1))))))
