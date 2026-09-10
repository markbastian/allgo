(ns procedurals.numerics-test
  (:require [procedurals.numerics.core :as core]
            [procedurals.numerics.rk :as rk]
            [clojure.test :refer [deftest is testing]]))

;; y' = y, y(0) = 1  =>  y(1) = e. A closed form to measure against.
(defn- exponential [_ y] y)
(def ^:private target (Math/exp 1.0))

(defn- global-error [method h]
  (let [s (core/step-until (rk/integrator method exponential 0.0 [1.0] h {:adaptive? false}) 1.0)]
    (abs (- (first (:y s)) target))))

(defn- observed-order
  "Halving the step should divide the global error by 2^p."
  [method]
  (let [es (mapv #(global-error method %) [0.1 0.05 0.025 0.0125])]
    (mapv (fn [[a b]] (/ (Math/log (/ a b)) (Math/log 2.0))) (partition 2 1 es))))

(deftest tableaux-are-consistent
  ;; M&G eq. 4.13: a stage's weights sum to its node, and the solution
  ;; weights to one. Cheap, and it catches a mistyped coefficient -- which
  ;; is exactly how the 125/192 in Dormand-Prince announced itself.
  (doseq [{:keys [name c a b b-hat]} rk/catalog]
    (testing name
      (doseq [i (range (count c))]
        (is (< (abs (- (reduce + 0.0 (nth a i)) (nth c i))) 1e-12)
            (str name " row " i " must sum to c")))
      (is (< (abs (- (reduce + 0.0 b) 1.0)) 1e-12) (str name " b must sum to 1"))
      (when b-hat
        (is (< (abs (- (reduce + 0.0 b-hat) 1.0)) 1e-12) (str name " b-hat must sum to 1"))))))

(deftest methods-converge-at-their-stated-order
  (doseq [method rk/catalog]
    (testing (:name method)
      (doseq [p (observed-order method)]
        ;; Generous below, tight above: a method may not beat its order, but
        ;; roundoff and pre-asymptotic behaviour cost a little beneath it.
        (is (< (- (:order method) 0.35) p (+ (:order method) 0.35))
            (str (:name method) " claims order " (:order method) ", observed " p))))))

(deftest higher-order-is-more-accurate-at-equal-step
  (let [errs (mapv #(global-error % 0.01) rk/catalog)]
    (is (apply > errs) "error falls monotonically as order rises")))

(deftest step-until-lands-exactly-on-the-end-time
  ;; Overshooting by up to a step makes any error comparison meaningless.
  (doseq [h [0.1 0.03 0.007 (/ 1.0 3.0)]]
    (let [s (core/step-until (rk/integrator rk/rk4 exponential 0.0 [1.0] h {:adaptive? false}) 1.0)]
      (is (< (abs (- (:t s) 1.0)) 1e-12) (str "h=" h " ended at " (:t s))))))

(deftest adaptive-stepping-honours-its-tolerance
  (doseq [method [rk/rkf45 rk/dopri54]]
    (testing (:name method)
      (doseq [tol [1e-4 1e-6 1e-8 1e-10]]
        (let [s (core/step-until
                 (rk/integrator method exponential 0.0 [1.0] 0.1 {:tol-abs tol :tol-rel tol}) 1.0)
              err (abs (- (first (:y s)) target))]
          (is (< err (* 200 tol))
              (str (:name method) " tol " tol " gave error " err)))))))

(deftest adaptive-stepping-adjusts-the-step
  (testing "a stiffening problem forces the step down and it recovers after"
    ;; y' = -100 y for a while, then benign. The step must shrink and regrow.
    (let [f     (fn [t y] (if (< t 0.2) (core/v* y -100.0) (core/v* y -0.1)))
          integ (rk/integrator rk/dopri54 f 0.0 [1.0] 0.05 {:tol-abs 1e-8 :tol-rel 1e-8})
          hs    (->> (core/trajectory integ)
                     (take-while #(< (:t %) 1.0))
                     (mapv :h))]
      (is (seq hs))
      (is (< (apply min hs) 0.02) "step shrinks through the fast phase")
      (is (> (apply max hs) (apply min hs)) "and grows again once it can"))))

(deftest rejected-steps-do-not-advance-time
  (let [integ (assoc (rk/integrator rk/dopri54 exponential 0.0 [1.0] 5.0
                                    {:tol-abs 1e-12 :tol-rel 1e-12}) :h 5.0)
        after (core/step integ)]
    (is (not (:accepted after)) "an absurd step must be rejected")
    (is (= 0.0 (:t after)) "time stands still on a rejection")
    (is (= [1.0] (:y after)) "and so does the state")
    (is (< (:h after) 5.0) "only the step shrinks")))

(deftest a-system-of-equations
  (testing "the harmonic oscillator as a first-order pair keeps its energy"
    (let [f (fn [_ [x v]] [v (- x)])
          s (core/step-until (rk/integrator rk/dopri54 f 0.0 [1.0 0.0] 0.01
                                            {:tol-abs 1e-12 :tol-rel 1e-12})
                             (* 4 Math/PI))
          [x v] (:y s)]
      (is (< (abs (- x 1.0)) 1e-8) "back to where it started after two periods")
      (is (< (abs v) 1e-8))
      (is (< (abs (- (+ (* x x) (* v v)) 1.0)) 1e-8) "energy conserved"))))
