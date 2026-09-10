(ns procedurals.numerics-test
  (:require [procedurals.numerics.core :as core]
            [procedurals.numerics.rk :as rk]
            [procedurals.numerics.rkn :as rkn]
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

;; ------------------------------------------------------- Runge-Kutta-Nystrom

;; y'' = -y, y(0)=1, y'(0)=0  =>  y = cos t. Measured away from a period
;; boundary: at a whole period the oscillator's symmetry cancels error and
;; every method reads a full order or two better than it really is.
(defn- oscillator [_ y] (core/v* y -1.0))

(defn- rkn-order [method which t-end]
  (let [exact (if (= which :y) #(Math/cos %) #(- (Math/sin %)))
        errs  (mapv (fn [h]
                      (let [s (core/step-until (rkn/integrator method oscillator 0.0 [1.0] [0.0] h) t-end)]
                        (max 1e-18 (abs (- (first (which s)) (exact t-end))))))
                    [0.1 0.05 0.025 0.0125])]
    (mapv (fn [[a b]] (/ (Math/log (/ a b)) (Math/log 2.0))) (partition 2 1 errs))))

(deftest nystrom-tableaux-are-consistent
  ;; The Nystrom conditions differ from Runge-Kutta's: a stage's weights sum
  ;; to half the square of its node, position weights to 1/2, velocity
  ;; weights to 1.
  (doseq [{:keys [name c a b b-dot]} rkn/catalog]
    (testing name
      (doseq [i (range (count c))]
        (is (< (abs (- (reduce + 0.0 (nth a i)) (/ (* (nth c i) (nth c i)) 2.0))) 1e-14)
            (str name " row " i)))
      (is (< (abs (- (reduce + 0.0 b) 0.5)) 1e-14) (str name " b sums to 1/2"))
      (is (< (abs (- (reduce + 0.0 b-dot) 1.0)) 1e-14) (str name " b-dot sums to 1")))))

(deftest nystrom-methods-converge-at-their-stated-order
  (doseq [method rkn/catalog]
    (testing (:name method)
      (doseq [which [:y :dy]
              p     (rkn-order method which 1.0)]
        (is (< (- (:order method) 0.35) p (+ (:order method) 0.35))
            (str (:name method) " " which " claims " (:order method) ", observed " p))))))

(deftest nystrom-on-a-nonlinear-system
  (testing "y'' = -y^3, which has no symmetry to flatter the method"
    (let [f    (fn [_ y] (mapv (fn [yi] (- (* yi yi yi))) y))
          ref  (first (:y (core/step-until (rkn/integrator rkn/rkn4 f 0.0 [1.0] [0.5] 1e-5) 1.0)))
          errs (mapv (fn [h]
                       (abs (- (first (:y (core/step-until (rkn/integrator rkn/rkn4 f 0.0 [1.0] [0.5] h) 1.0)))
                               ref)))
                     [0.1 0.05 0.025])
          orders (mapv (fn [[a b]] (/ (Math/log (/ a b)) (Math/log 2.0))) (partition 2 1 errs))]
      (doseq [p orders]
        (is (< 3.6 p 4.4) (str "observed order " p))))))

(deftest nystrom-beats-runge-kutta-per-force-evaluation
  ;; The reason Nystrom methods exist. RK4 on the doubled system spends four
  ;; evaluations a step; RKN4 spends three, because the two midpoint stages
  ;; of RK4 ask for the same force at the same point.
  (testing "same accuracy for fewer evaluations"
    (let [t-end   (* 2.0 Math/PI)
          exact   (Math/cos t-end)
          rk-err  (fn [h] (let [f (fn [_ [x v]] [v (- x)])
                                s (core/step-until (rk/integrator rk/rk4 f 0.0 [1.0 0.0] h {:adaptive? false}) t-end)]
                            (abs (- (first (:y s)) exact))))
          rkn-err (fn [h] (let [s (core/step-until (rkn/integrator rkn/rkn4 oscillator 0.0 [1.0] [0.0] h) t-end)]
                            (abs (- (first (:y s)) exact))))
          ;; Equal work per unit time: RK4 spends 4/h evaluations, RKN4
          ;; spends 3/H, so H = 3h/4 buys Nystrom the same budget.
          h       0.05]
      (is (= 4 (:stages rk/rk4)))
      (is (= 3 (:stages rkn/rkn4)))
      (is (< (rkn-err (* h 0.75)) (rk-err h))
          "for the same force-evaluation budget, RKN4 is the more accurate")
      (is (< (rkn-err h) (* 3.0 (rk-err h)))
          "and at equal step it is within a small factor while doing 3/4 the work"))))

(deftest nystrom-step-doubling-honours-its-tolerance
  (doseq [tol [1e-6 1e-8 1e-10]]
    (let [s (core/step-until (rkn/integrator rkn/rkn4 oscillator 0.0 [1.0] [0.0] 0.1
                                             {:adaptive? true :tol-abs tol :tol-rel tol})
                             10.0)]
      (is (< (abs (- (first (:y s)) (Math/cos 10.0))) (* 100 tol))
          (str "tol " tol)))))
