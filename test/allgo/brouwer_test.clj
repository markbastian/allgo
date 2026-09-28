(ns allgo.brouwer-test
  "Brouwer's short-period terms against J2's motion integrated
  numerically: what they leave is second order in J2 -- a tenth of J2,
  a hundredth of the error -- at every eccentricity and inclination,
  where Spacetrack Report No. 3's shortened terms leave first-order
  errors growing with e."
  (:require [allgo.astro.brouwer :as br]
            [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.perturbations :as pt]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(def ^:private R c/R-earth)

(defn- j2-perturbation
  "J2's acceleration less the central term."
  [J2]
  (fn [r]
    (v3/add (geo/acceleration {:GM mu :R R :normalized? false :C {[0 0] 1.0 [2 0] (- J2)} :S {}} r 2)
            (v3/scale r (/ mu (math/pow (v3/length r) 3))))))

(defn- rk4 [accel [r v] t n]
  (let [h (/ t n)
        f (fn [[r v]] [v (v3/add (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))) (accel r))])
        add (fn [[r v] [dr dv] k] [(v3/add-scaled r dr k) (v3/add-scaled v dv k)])]
    (loop [s [r v] k 0]
      (if (= k n)
        s
        (let [k1 (f s) k2 (f (add s k1 (* 0.5 h))) k3 (f (add s k2 (* 0.5 h))) k4 (f (add s k3 h))]
          (recur (-> s (add k1 (/ h 6)) (add k2 (/ h 3)) (add k3 (/ h 3)) (add k4 (/ h 6))) (inc k)))))))

(defn- error-after-three-orbits
  "How far J2's motion, integrated from the osculating state of mean
  elements `el`, is from the theory's after three orbits, km."
  [J2 el]
  (let [t (* 3 (kep/period mu (:a el)))
        [r] (rk4 (j2-perturbation J2) (br/osculating-state mu R J2 el) t 3000)]
    (v3/distance r (first (br/propagate mu R J2 el t)))))

(def ^:private orbits
  [{:a 7000.0 :e 0.001 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7000.0 :e 0.01 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 8000.0 :e 0.1 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 14000.0 :e 0.5 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7000.0 :e 0.01 :i 0.001 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7200.0 :e 0.05 :i 2.8 :raan 0.3 :argp 1.0 :M 0.2}])

(deftest semi-major-axis
  (testing "the change in a is the closed form n dW1/dl gives"
    (doseq [{:keys [a e i argp M] :as el} orbits]
      (let [[dL] (br/short-period el)
            da (/ (* 2.0 (math/sqrt (* mu a)) dL) mu)
            gamma2 (/ (* 0.5 geo/J2 R R) (* a a))
            E (kep/kepler-equation M e)
            f (kep/eccentric->true E e)
            ar3 (math/pow (/ 1.0 (- 1.0 (* e (math/cos E)))) 3)
            eta3 (math/pow (- 1.0 (* e e)) 1.5)
            th2 (math/pow (math/cos i) 2)
            closed (* a gamma2 (+ (* (- (* 3.0 th2) 1.0) (- ar3 (/ 1.0 eta3)))
                                  (* 3.0 (- 1.0 th2) ar3 (math/cos (* 2.0 (+ argp f))))))]
        (is (< (abs (- da closed)) (* 1e-9 (+ (abs closed) (* a gamma2)))) (str el))))))

(deftest against-numerical-j2
  (testing "what is left is J2^2: a tenth of J2 leaves a hundredth of the error"
    (doseq [el orbits]
      (let [full (error-after-three-orbits geo/J2 el)
            tenth (error-after-three-orbits (* 0.1 geo/J2) el)]
        ;; J2^2's along-track drift, some 70 m an orbit low and inclined,
        ;; more near the equator where the node and perigee turn fastest
        (is (< full 2.0) (str el " " full))
        (is (< 90.0 (/ full tenth) 110.0) (str el " " full " " tenth)))))
  (testing "where the report's shortened terms, missing J2 e, are nearly a hundred times further out at e = 0.1"
    (let [el (orbits 2)
          t (* 3 (kep/period mu (:a el)))
          [r] (rk4 (j2-perturbation geo/J2) (pt/j2-osculating el) t 3000)
          report (v3/distance r (first (pt/j2-propagate el t)))]
      (is (> report (* 50.0 (error-after-three-orbits geo/J2 el)))))))

(deftest osculating-to-mean
  (testing "mean to osculating and back"
    (doseq [{:keys [a e i raan argp M] :as el} orbits]
      (let [m (br/mean-elements (br/osculating-state el))]
        (is (< (abs (- (:a m) a)) 1e-8) (str el))
        (is (< (abs (- (:e m) e)) 1e-10) (str el))
        (is (< (abs (- (:i m) i)) 1e-11) (str el))
        (is (< (abs (am/wrap-angle (- (:raan m) raan))) 1e-11) (str el))
        (is (< (abs (am/wrap-angle (- (+ (:argp m) (:M m)) argp M))) 1e-11) (str el))))))
