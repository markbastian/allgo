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

(defn- exponential-drag
  "Drag in still air falling off exponentially from rho-p at radius rp."
  [B rho-p rp H]
  (fn [r v]
    (let [rho (* rho-p (math/exp (- (/ (- (v3/length r) rp) H))))]
      (v3/scale v (* -0.5 B rho 1e3 (v3/length v))))))

(deftest drag
  (testing "King-Hele's decay per revolution is drag averaged round the orbit"
    (let [B 0.01 rho-p 3e-11 H 60.0 rp (+ c/R-earth 300.0)]
      (doseq [e [0.002 0.01 0.03]]
        (let [a (/ rp (- 1.0 e))
              el {:a a :e e :i 0.9 :raan 0.1 :argp 0.3 :M 0.0}
              {:keys [a-rate e-rate]} (pt/drag-secular a e B rho-p H)
              averaged (pt/averaged-rates mu el (exponential-drag B rho-p rp H) 3600)]
          ;; first order in e: they part by e^2
          (is (close? (:a averaged) a-rate (* 3 e e)) (str "a, e = " e))
          (is (close? (:e averaged) e-rate (+ 1e-3 (* 3 e))) (str "e, e = " e))))))
  (testing "a circular orbit's decay rate is Gauss's equation with along-track drag"
    (let [a 6778.0 B 0.01 rho 3e-12
          el {:a a :e 1e-9 :i 0.9 :raan 0.1 :argp 0.3 :M 0.0}
          averaged (pt/averaged-rates mu el (exponential-drag B rho a 1e9) 72)]
      (is (close? (:a averaged) (pt/circular-decay-rate a B rho) 1e-9))))
  (testing "the lifetime is the decay rate integrated in time"
    (let [B 0.01 density (fn [h] (* 3e-11 (math/exp (- (/ (- h 300.0) 55.0)))))
          a0 (+ c/R-earth 350.0) a1 (+ c/R-earth 250.0)
          t (pt/circular-lifetime a0 a1 B density)
          ;; step the orbit down with RK4 in time until it reaches a1
          rate (fn [a] (pt/circular-decay-rate a B (density (- a c/R-earth))))
          dt 600.0
          t' (loop [a a0 t 0.0]
               (let [k1 (rate a) k2 (rate (+ a (* 0.5 dt k1))) k3 (rate (+ a (* 0.5 dt k2))) k4 (rate (+ a (* dt k3)))
                     a' (+ a (* (/ dt 6) (+ k1 (* 2 k2) (* 2 k3) k4)))]
                 (if (<= a' a1)
                   (+ t (* dt (/ (- a a1) (- a a'))))
                   (recur a' (+ t dt)))))]
      (is (close? t t' 1e-4)))))

(deftest steady-push
  (let [el {:a 26000.0 :e 0.3 :i 0.6 :raan 0.4 :argp 1.3 :M 0.0}
        F [3e-10 -1e-10 2e-10]
        s (kep/elements->state mu el)
        {he :e hh :h} (pt/steady-push-secular s F)
        ;; the exact rates, averaged numerically round the orbit
        n 720
        states (for [k (range n)] (kep/elements->state mu (assoc el :M (* 2 math/PI (/ k n)))))
        mean (fn [f] (v3/scale (reduce v3/add (map f states)) (/ 1.0 n)))
        dh (mean (fn [[r _]] (v3/cross r F)))
        de (mean (fn [[r v]] (v3/scale (v3/add (v3/cross F (v3/cross r v)) (v3/cross v (v3/cross r F))) (/ 1.0 mu))))]
    (testing "the closed forms are the exact rates averaged"
      (is (< (v3/distance hh dh) (* 1e-10 (v3/length dh))))
      (is (< (v3/distance he de) (* 1e-10 (v3/length de)))))
    (testing "and agree with Gauss's equations averaged: no change in a, and e's rate"
      (let [averaged (pt/averaged-rates mu el (constantly F) 720)
            e-hat (v3/normalize (pt/eccentricity-vector s))]
        (is (< (abs (:a averaged)) 1e-15))
        (is (close? (:e averaged) (v3/dot e-hat he) 1e-9)))))
  (testing "radiation pressure pushes from the Sun"
    (let [s (kep/elements->state mu {:a 26000.0 :e 0.3 :i 0.6 :raan 0.4 :argp 1.3 :M 0.0})
          {:keys [e]} (pt/srp-secular s [c/AU 0.0 0.0] 0.02 1.3)
          F (pt/steady-push-secular s [-1.0 0.0 0.0])]
      ;; the push is along -x, so the rates are those of a push along -x, scaled
      (is (< (v3/length (v3/cross e (:e F))) (* 1e-12 (v3/length e) (v3/length (:e F)))))
      (is (pos? (v3/dot e (:e F)))))))
