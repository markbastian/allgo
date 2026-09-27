(ns allgo.vallado-encke-test
  "Encke's method against Cowell's: the same J2-perturbed orbit integrated
  as a deviation from a conic and as a whole, by the same integrator, must
  land in the same place -- and Encke's, integrating something small and
  slow, should get there in fewer steps."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.encke :as encke]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)

(defn- j2 [_ r _]
  (v3/add (geo/acceleration {:GM mu :R c/R-earth :normalized? false :C {[0 0] 1.0 [2 0] (- geo/J2)} :S {}} r 2)
          (v3/scale r (/ mu (math/pow (v3/length r) 3)))))

(def ^:private s0 (kep/elements->state mu {:a 7200.0 :e 0.05 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}))

(defn- cowell
  "The whole motion, two-body plus J2, by the same integrator."
  [[r0 v0] t]
  (let [f (fn [t y] (let [r (subvec y 0 3) v (subvec y 3 6)]
                      (into v (v3/add (v3/scale r (- (/ mu (math/pow (v3/length r) 3)))) (j2 t r v)))))
        integ (rk/integrator rk/dopri54 f 0.0 (into r0 v0) 60.0 {:tol-abs 1e-9 :tol-rel 1e-12})
        steps (count (take-while #(< (:t %) t) (core/trajectory integ)))
        {:keys [y]} (core/step-until integ t)]
    {:state [(subvec y 0 3) (subvec y 3 6)] :steps steps}))

(deftest battins-f
  (testing "f(q) is (1+q)^3/2 - 1"
    (doseq [q [-0.3 -1e-3 1e-3 0.5 2.0]]
      (is (< (abs (- (encke/f-of-q q) (- (math/pow (+ 1 q) 1.5) 1.0))) (* 1e-14 (max 1.0 (abs q)))))))
  (testing "and keeps its digits where the subtraction would lose them"
    (let [q 1e-13]
      ;; to first order f = 3q/2, with the next term 3q^2/8
      (is (< (abs (- (encke/f-of-q q) (+ (* 1.5 q) (* 0.375 q q)))) (* 1e-15 q))))))

(deftest no-perturbation
  (testing "with nothing to integrate, the reference conic is the answer"
    (let [[r0 v0] s0
          {:keys [state rectifications]} (encke/propagate mu r0 v0 (fn [_ _ _] [0.0 0.0 0.0]) 20000.0 {})
          [r] (u/propagate mu r0 v0 20000.0)]
      (is (< (v3/distance (first state) r) 1e-9))
      (is (zero? rectifications)))))

(deftest against-cowell
  (let [[r0 v0] s0 t 86400.0
        {e-state :state e-steps :steps rects :rectifications} (encke/propagate mu r0 v0 j2 t {})
        {c-state :state c-steps :steps} (cowell s0 t)]
    (testing "a day of J2: the two methods land together"
      (is (< (v3/distance (first e-state) (first c-state)) 1e-4))
      (is (< (v3/distance (second e-state) (second c-state)) 1e-7)))
    (testing "rectifying as the deviation grows"
      (is (pos? rects)))
    (testing "in fewer steps, integrating something small and slow"
      (is (< e-steps c-steps) (str e-steps " against " c-steps)))))
