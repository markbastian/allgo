(ns allgo.gauss-jackson-test
  "Gauss-Jackson against Kepler's orbit, known exactly, and against the
  coefficients its series are known to give."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kep]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.gauss-jackson :as gj]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest coefficients
  (testing "the Cowell and Adams-Moulton-type series, exactly"
    (let [{:keys [g v]} (gj/coefficients 7)]
      (is (= [1 -1 1/12 0 -1/240 -1/240 -221/60480] g))
      (is (= [1 -1/2 -1/12 -1/24 -19/720 -3/160 -863/60480] v)))))

(def ^:private mu c/GM-earth)
(def ^:private s0 (kep/elements->state mu {:a 8000.0 :e 0.1 :i 0.5 :raan 0.3 :argp 0.7 :M 0.0}))
(defn- kepler-accel [_ r] (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))))

(defn- error-after [h t]
  (let [steps (long (/ t h))
        [[tn y]] [(peek (gj/integrate kepler-accel 0.0 (first s0) (second s0) h steps))]
        [r] (kep/propagate mu (first s0) (second s0) tn)]
    (v3/distance y r)))

(deftest against-kepler
  (let [period (kep/period mu 8000.0)]
    (testing "ten orbits at a minute a step to centimeters, at half a minute to a tenth of a millimeter"
      (is (< (error-after 60.0 (* 10 period)) 2e-5))
      (is (< (error-after 30.0 (* 10 period)) 1e-7)))
    (testing "tenth order: halving the step cuts the error by hundreds"
      (let [e1 (error-after 240.0 (* 2 period)) e2 (error-after 120.0 (* 2 period))]
        (is (> (/ e1 e2) 200.0) (str e1 " " e2))))))
