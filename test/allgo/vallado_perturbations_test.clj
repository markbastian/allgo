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

(deftest j2-short-period
  (let [per (kep/period mu 7000.0)
        t (* 3 per)
        error (fn [e with?]
                (let [mean {:a 7000.0 :e e :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
                      s0 (pt/j2-osculating mean)
                      num (rk4 (minus-central j2-accel) s0 t 2400)
                      rates (pt/j2-secular 7000.0 e 0.9)
                      mean-t (assoc mean :raan (+ 0.3 (* (:raan rates) t)) :argp (+ 1.0 (* (:argp rates) t))
                                    :M (+ 0.2 (* (:M rates) t)))
                      model (if with? (pt/j2-propagate mean t) (kep/elements->state mu mean-t))]
                  (v3/distance (first num) (first model))))]
    (testing "three orbits of J2's motion, integrated, against the mean elements with the short-period terms"
      ;; J2^2's along-track drift, some 45 m an orbit, is what first order leaves
      (is (< (error 0.0 true) 0.2))
      (is (< (error 0.001 true) 0.2)))
    (testing "and without them, the orbit is kilometers out"
      (is (> (error 0.0 false) (* 10 (error 0.0 true))))
      (is (> (error 0.001 false) 2.0)))
    (testing "the report drops the terms of order J2 e: at 0.01 they show"
      (is (< (error 0.01 true) 2.0))))
  (testing "osculating to mean and back"
    (let [mean {:a 7000.0 :e 0.01 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
          m (pt/osculating->mean (pt/j2-osculating mean))]
      (doseq [k [:a :e :i :raan]]
        (is (< (abs (- (m k) (mean k))) (* 1e-9 (max 1.0 (mean k)))) (str k)))
      ;; the perigee and anomaly, each ill-set by a small eccentricity, together
      (is (< (abs (- (+ (:argp m) (:M m)) (+ (:argp mean) (:M mean)))) 1e-9)))))

(deftest lagranges-equations
  (testing "the averaged J2 disturbing function gives J2's secular rates exactly"
    (doseq [[a e i] [[8000.0 0.1 0.9] [7078.0 0.001 1.71] [26000.0 0.7 1.1]]]
      (let [R-bar (fn [{:keys [a e i]}]
                    (* (/ (* mu geo/J2 c/R-earth c/R-earth) (* 4.0 a a a (math/pow (- 1.0 (* e e)) 1.5)))
                       (- 2.0 (* 3.0 (math/pow (math/sin i) 2)))))
            rates (pt/lagrange-rates-of {:a a :e e :i i :raan 0.4 :argp 1.2 :M 0.5} R-bar)
            closed (pt/j2-secular a e i)]
        (doseq [k [:raan :argp :M]]
          (is (close? (rates k) (closed k) 1e-8) (str k " " [a e i])))
        (doseq [k [:a :e :i]]
          (is (< (abs (rates k)) 1e-15) (str k))))))
  (testing "J2's instantaneous disturbing function gives what Gauss's equations give from its acceleration"
    (let [el {:a 8000.0 :e 0.1 :i 0.9 :raan 0.4 :argp 1.2 :M 0.5}
          R-j2 (fn [el'] (let [[[_ _ z :as r]] (kep/elements->state mu el')
                               rm (v3/length r)]
                           (* (/ (* -0.5 mu geo/J2 c/R-earth c/R-earth) (* rm rm rm))
                              (- (* 3.0 (/ (* z z) (* rm rm))) 1.0))))
          lag (pt/lagrange-rates-of el R-j2)
          s (kep/elements->state mu el)
          gauss (pt/gauss-rates el (pt/rsw s ((minus-central j2-accel) (first s) (second s))))
          n (kep/mean-motion mu 8000.0)]
      (doseq [k [:a :e :i :raan :argp]]
        (is (close? (lag k) (gauss k) 1e-6) (str k)))
      (is (close? (- (:M lag) n) (:M gauss) 1e-6) "M"))))

(deftest third-bodies
  (let [el {:a 8000.0 :e 0.2 :i 0.9 :raan 0.4 :argp 1.2 :M 0.0}]
    (testing "the orbit's second moment, averaged numerically"
      (let [n 720
            M (reduce (fn [acc k] (let [[r] (kep/elements->state mu (assoc el :M (* 2 math/PI (/ k n))))]
                                    (mapv (fn [row x] (mapv #(+ %1 (/ (* x %2) n)) row r)) acc r)))
                      [[0.0 0.0 0.0] [0.0 0.0 0.0] [0.0 0.0 0.0]] (range n))]
        (is (every? #(< (abs %) 1e-6) (flatten (mapv (fn [a b] (mapv - a b)) M (pt/orbit-moment el)))))))
    (testing "the doubly averaged rates are the tidal acceleration itself averaged over both orbits"
      (let [mu3 c/GM-moon r3 384400.0
            h3 (v3/normalize [0.0 (- (math/sin 0.4)) (math/cos 0.4)])
            u3 [1.0 0.0 0.0] w3 (v3/cross h3 u3)
            secular (pt/third-body-secular el mu3 r3 h3)
            ;; the third body at 24 places round its orbit, each averaged over the satellite's
            k 24
            tidal (fn [n3] (fn [r _] (v3/scale (v3/sub (v3/scale n3 (* 3.0 (v3/dot r n3))) r) (/ mu3 (* r3 r3 r3)))))
            brute (let [rs (for [j (range k) :let [th (* 2 math/PI (/ j k))
                                                   n3 (v3/add (v3/scale u3 (math/cos th)) (v3/scale w3 (math/sin th)))]]
                             (pt/averaged-rates mu el (tidal n3) 360))]
                    (into {} (for [key [:e :i :raan :argp]] [key (/ (reduce + (map key rs)) k)])))]
        (doseq [key [:e :i :raan :argp]]
          (is (< (abs (- (secular key) (brute key))) (* 1e-4 (apply max (map #(abs (brute %)) [:raan :argp])))) (str key))))))
  (testing "for a circular orbit, the classical closed forms and the familiar coefficients"
    (let [a 7000.0 i 0.9 n (kep/mean-motion mu a)
          eps (math/to-radians 23.44)
          ;; the Moon's pull, as the Earth-Moon system feels it
          {:keys [raan argp]} (pt/third-body-node-perigee n i c/GM-moon 384400.0 eps)
          rev-per-day (/ (* n 86400.0) (* 2 math/PI))
          deg-per-day #(math/to-degrees (* % 86400.0))]
      (is (< (abs (- (deg-per-day raan) (/ (* -0.00338 (math/cos i)) rev-per-day))) 1e-5))
      (is (< (abs (- (deg-per-day argp) (/ (* 0.00169 (- 4.0 (* 5.0 (math/pow (math/sin i) 2)))) rev-per-day))) 1e-5)))))
