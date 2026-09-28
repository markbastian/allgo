(ns allgo.semianalytic-test
  "The semi-analytic theory against what it must reproduce: the
  equinoctial elements' round trip; the averaged J2 rates the closed-form
  secular rates; the mean-to-osculating map and its inverse; and three
  days of J2 and J3, and two of drag, against the motion integrated
  numerically -- meters with the second-order averages, kilometers
  without."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.perturbations :as pt]
            [allgo.astro.semianalytic :as sa]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)

(defn- zonal [n J2 J3]
  (let [field {:GM mu :R c/R-earth :normalized? false :C {[0 0] 1.0 [2 0] (- J2) [3 0] J3} :S {}}]
    (fn [_ r _] (v3/add (geo/acceleration field r n) (v3/scale r (/ mu (math/pow (v3/length r) 3)))))))

(defn- integrate
  "The motion under two-body gravity plus `pert`, by Dormand-Prince at
  tight tolerance: the state at `t1`."
  [pert [r v] t1]
  (let [f (fn [t y] (let [r (subvec y 0 3) v (subvec y 3 6)]
                      (into v (v3/add (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))) (pert t r v)))))
        y (:y (core/step-until (rk/integrator rk/dopri54 f 0.0 (into r v) 30.0 {:tol-abs 1e-11 :tol-rel 1e-11}) t1))]
    [(subvec y 0 3) (subvec y 3 6)]))

(def ^:private s0 (kep/elements->state mu {:a 7200.0 :e 0.01 :i 0.9 :raan 0.4 :argp 1.1 :M 2.0}))

(deftest elements
  (testing "equinoctial elements and back"
    (doseq [el [{:a 7200.0 :e 0.05 :i 0.9 :raan 0.4 :argp 1.1 :M 2.0}
                {:a 42164.0 :e 0.0001 :i 0.001 :raan 3.0 :argp 0.5 :M 5.0}
                {:a 26000.0 :e 0.7 :i 1.1 :raan 5.0 :argp 4.0 :M 0.3}]]
      (let [[r v :as s] (kep/elements->state mu el)
            [r' v'] (sa/elements->state (sa/state->elements s))]
        (is (< (v3/distance r r') 1e-8))
        (is (< (v3/distance v v') 1e-11))))))

(deftest averaged-rates
  (testing "J2's rates averaged round the orbit are its closed-form secular rates"
    (let [x (sa/state->elements s0)
          {:keys [a e i]} (kep/state->elements mu (first s0) (second s0))
          [_ af ag chi psi] x
          [_ daf dag dchi dpsi dlam] (sa/mean-rates (zonal 2 geo/J2 0.0) mu 0.0 x 64 false)
          closed (pt/j2-secular a e i)
          node (/ (- (* psi dchi) (* chi dpsi)) (+ (* chi chi) (* psi psi)))
          perigee-longitude (/ (- (* af dag) (* ag daf)) (+ (* af af) (* ag ag)))]
      (is (< (abs (- node (:raan closed))) (* 1e-8 (abs (:raan closed)))))
      (is (< (abs (- perigee-longitude (+ (:raan closed) (:argp closed)))) (* 1e-8 (abs (:argp closed)))))
      (is (< (abs (- dlam (+ (:raan closed) (:argp closed) (:M closed)))) (* 1e-10 (:M closed)))))))

(deftest mean-and-osculating
  (testing "the mean elements of an osculating state give it back"
    (let [pert (zonal 3 geo/J2 2.53e-6)
          x (sa/mean-elements pert 0.0 s0)
          [r v] (sa/osculating pert 0.0 x)]
      (is (< (v3/distance r (first s0)) 1e-6))
      (is (< (v3/distance v (second s0)) 1e-9)))))

(defn- errors
  "Along the numerically integrated trajectory's radial, along-track and
  cross-track axes, the semi-analytic state's error after `t` seconds."
  [pert t second?]
  (let [[y v] (integrate pert s0 t)
        x (sa/propagate-mean pert mu (sa/mean-elements pert 0.0 s0) 0.0 t 21600.0 second?)
        [r] (sa/osculating pert t x)
        d (v3/sub r y)
        R (v3/normalize y) W (v3/normalize (v3/cross y v)) S (v3/cross W R)]
    (mapv #(abs (v3/dot d %)) [R S W])))

(deftest zonal-harmonics
  (testing "three days of J2 and J3 at six-hour steps: meters with the second-order averages"
    (let [pert (zonal 3 geo/J2 2.53e-6)
          [radial along cross] (errors pert (* 3 86400.0) true)
          [_ along1] (errors pert (* 3 86400.0) false)]
      (is (< radial 0.01))
      (is (< along 0.03))
      (is (< cross 0.005))
      (testing "and kilometers along track without them"
        (is (> along1 1.0))))))

(deftest drag
  (testing "two days of drag in an exponential atmosphere, the decay the averaged drag gives"
    (let [pert (fn [_ r v]
                 (let [h (- (v3/length r) c/R-earth)
                       rho (* 3e-12 (math/exp (/ (- 400.0 h) 60.0)))            ; kg/m^3
                       ;; B = 0.02 m^2/kg; the air turning with the Earth
                       vrel (v3/sub v (v3/cross [0.0 0.0 7.292115e-5] r))]
                   (v3/scale vrel (* -0.5 0.02 rho 1e3 (v3/length vrel)))))
          s0 (kep/elements->state mu {:a (+ c/R-earth 400.0) :e 0.001 :i 0.9 :raan 0.4 :argp 1.1 :M 2.0})
          t (* 2 86400.0)
          [y] (integrate pert s0 t)
          x (sa/propagate-mean pert mu (sa/mean-elements pert 0.0 s0) 0.0 t 21600.0 true)
          [r] (sa/osculating pert t x)
          a-decay (- (first (sa/state->elements s0)) (first (sa/state->elements (integrate pert s0 t))))]
      (is (> a-decay 0.05) "the orbit does decay")
      (is (< (v3/distance r y) 0.05)))))
