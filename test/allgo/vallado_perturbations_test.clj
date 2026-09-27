(ns allgo.vallado-perturbations-test
  "Vallado chapters 8 and 9's general perturbations, each checked against
  what it summarizes: Gauss's variational equations against the elements
  of an orbit integrated under the force, and each closed-form secular
  rate against Gauss's rates averaged numerically round the orbit."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.perturbations :as pt]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1e-300 (abs b)))))

(defn- rk4
  "Integrate two-body motion plus `(accel r v)` from `[r v]` for `t`
  seconds in `n` steps."
  [accel [r v] t n]
  (let [h (/ t n)
        f (fn [[r v]] [v (v3/add (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))) (accel r v))])
        add (fn [[r v] [dr dv] k] [(v3/add-scaled r dr k) (v3/add-scaled v dv k)])]
    (loop [s [r v] k 0]
      (if (= k n)
        s
        (let [k1 (f s) k2 (f (add s k1 (* 0.5 h))) k3 (f (add s k2 (* 0.5 h))) k4 (f (add s k3 h))]
          (recur (-> s (add k1 (/ h 6)) (add k2 (/ h 3)) (add k3 (/ h 3)) (add k4 (/ h 6))) (inc k)))))))

(def ^:private el0 {:a 8000.0 :e 0.1 :i (math/to-radians 50) :raan 0.7 :argp 1.1 :M 0.4})

(deftest gauss-variational-equations
  (testing "the element rates are the rates of the elements of the orbit the force bends"
    (let [f [2e-7 -3e-7 4e-7]
          accel (fn [_ _] f)
          s0 (kep/elements->state mu el0)
          dt 2.0
          elements (fn [[r v]] (kep/state->elements mu r v))
          e+ (elements (rk4 accel s0 dt 20))
          e- (elements (rk4 accel s0 (- dt) 20))
          rates (pt/gauss-rates el0 (pt/rsw s0 f))
          n (kep/mean-motion mu (:a el0))]
      (doseq [k [:a :e :i :raan :argp]]
        (let [d (- (e+ k) (e- k)) d (if (#{:raan :argp} k) (am/wrap-angle d) d)]
          (is (close? (rates k) (/ d (* 2 dt)) 1e-5) (str k))))
      (is (close? (+ n (:M rates)) (/ (am/wrap-angle (- (:M e+) (:M e-))) (* 2 dt)) 1e-7) "M"))))

(defn- j2-accel [r _]
  (geo/acceleration {:GM mu :R c/R-earth :normalized? false :C {[0 0] 1.0 [2 0] (- geo/J2)} :S {}} r 2))

(defn- minus-central [accel] (fn [r v] (v3/add (accel r v) (v3/scale r (/ mu (math/pow (v3/length r) 3))))))

(deftest j2-secular
  (testing "the node and perigee rates are J2's own acceleration averaged round the orbit"
    (doseq [el [el0 (assoc el0 :i (math/to-radians 98) :e 0.001) (assoc el0 :e 0.3 :i (math/to-radians 20))]]
      (let [{:keys [a e i]} el
            closed (pt/j2-secular a e i)
            averaged (pt/averaged-rates el (minus-central j2-accel))]
        ;; first order in J2: they agree but for J2^2, a part in a thousand
        (is (close? (:raan averaged) (:raan closed) 2e-3) (str el))
        (is (close? (:argp averaged) (:argp closed) 2e-3) (str el))
        (doseq [k [:a :e :i]]
          (is (< (abs (k averaged)) (* 1e-3 (abs (:raan closed)) (if (= k :a) a 1.0))) (str k " has no secular part"))))))
  (testing "the critical inclination freezes the perigee"
    (is (< (abs (:argp (pt/j2-secular 8000.0 0.1 (math/asin (math/sqrt 0.8))))) 1e-20)))
  (testing "a sun-synchronous orbit turns its node once a year"
    (let [i (pt/sun-synchronous-inclination 7078.0 0.0)]
      (is (> i (/ math/PI 2)) "retrograde")
      (is (close? (:raan (pt/j2-secular 7078.0 0.0 i)) (/ (* 2 math/PI) (* 365.2421897 86400.0)) 1e-12)))
    (is (nil? (pt/sun-synchronous-inclination 30000.0 0.0)) "too high for J2 to turn it that fast")))
