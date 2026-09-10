(ns procedurals.numerics-test
  (:require [procedurals.numerics.core :as core]
            [procedurals.numerics.multistep :as ms]
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

;; ------------------------------------------------------------- multistep

(defn- ratio-order [errs]
  (mapv (fn [[a b]] (/ (Math/log (/ a b)) (Math/log 2.0))) (partition 2 1 errs)))

(deftest adams-coefficients-match-the-published-tables
  ;; Derived from the defining integral rather than transcribed, so this
  ;; checks the derivation against the numbers everyone prints.
  (letfn [(scaled [cs d] (mapv #(Math/round (* d %)) cs))]
    (testing "Adams-Bashforth"
      (is (= [1] (scaled (ms/adams-bashforth-coefficients 1) 1)))
      (is (= [3 -1] (scaled (ms/adams-bashforth-coefficients 2) 2)))
      (is (= [23 -16 5] (scaled (ms/adams-bashforth-coefficients 3) 12)))
      (is (= [55 -59 37 -9] (scaled (ms/adams-bashforth-coefficients 4) 24)))
      (is (= [1901 -2774 2616 -1274 251] (scaled (ms/adams-bashforth-coefficients 5) 720))))
    (testing "Adams-Moulton"
      (is (= [1] (scaled (ms/adams-moulton-coefficients 1) 1)))
      (is (= [1 1] (scaled (ms/adams-moulton-coefficients 2) 2)))
      (is (= [5 8 -1] (scaled (ms/adams-moulton-coefficients 3) 12)))
      (is (= [9 19 -5 1] (scaled (ms/adams-moulton-coefficients 4) 24)))
      (is (= [251 646 -264 106 -19] (scaled (ms/adams-moulton-coefficients 5) 720))))
    (testing "Stoermer and Cowell"
      (is (= [1] (scaled (ms/stoermer-coefficients 1) 1)))
      (is (= [1 0] (scaled (ms/stoermer-coefficients 2) 1)))
      (is (= [13 -2 1] (scaled (ms/stoermer-coefficients 3) 12)))
      (is (= [1] (scaled (ms/cowell-coefficients 1) 1)))
      (is (= [1 10 1] (scaled (ms/cowell-coefficients 3) 12))))))

(deftest every-adams-weight-set-sums-to-one
  (doseq [k (range 1 9)]
    (is (< (abs (- (reduce + (ms/adams-bashforth-coefficients k)) 1.0)) 1e-12))
    (is (< (abs (- (reduce + (ms/adams-moulton-coefficients k)) 1.0)) 1e-12))
    ;; Stoermer-Cowell integrates twice, so its weights sum to one as well:
    ;; a constant acceleration g must give exactly h^2 g of second difference.
    (is (< (abs (- (reduce + (ms/stoermer-coefficients k)) 1.0)) 1e-12))
    (is (< (abs (- (reduce + (ms/cowell-coefficients k)) 1.0)) 1e-12))))

(deftest adams-methods-converge-at-their-stated-order
  (doseq [method ms/catalog]
    (testing (:name method)
      (let [errs (mapv (fn [h]
                         (max 1e-17 (abs (- (first (:y (core/step-until
                                                        (ms/integrator method exponential 0.0 [1.0] h) 1.0)))
                                            target))))
                       [0.05 0.025 0.0125])]
        (doseq [p (ratio-order errs)]
          (is (< (- (:order method) 0.5) p (+ (:order method) 0.5))
              (str (:name method) " observed " p)))))))

(deftest stoermer-cowell-converges-at-its-stated-order
  (doseq [method ms/catalog-2]
    (testing (:name method)
      (let [errs (mapv (fn [h]
                         (max 1e-16 (abs (- (first (:dy (core/step-until
                                                         (ms/integrator-2 method oscillator 0.0 [1.0] [0.0] h) 1.0)))
                                            (- (Math/sin 1.0))))))
                       [0.05 0.025 0.0125])]
        ;; Velocity rides an Adams sum over the same history, so it carries
        ;; the method's order rather than being capped by a difference.
        (doseq [p (ratio-order errs)]
          (is (< (- (:order method) 0.5) p (+ (:order method) 0.6))
              (str (:name method) " velocity observed " p)))))))

(deftest correcting-beats-predicting-alone
  (testing "PECE has a much smaller error constant than the predictor at equal order"
    (doseq [k [3 4 5]]
      (let [err (fn [m] (abs (- (first (:y (core/step-until
                                            (ms/integrator m exponential 0.0 [1.0] 0.02) 1.0)))
                                target)))]
        (is (< (err (ms/adams-pece k)) (err (ms/adams-bashforth k)))
            (str "ABM" k " should beat AB" k))))))

(deftest multistep-costs-one-or-two-evaluations-a-step
  ;; The reason to tolerate a history and a starting procedure at all.
  (is (every? #(= 1 (:stages %)) (mapv ms/adams-bashforth [2 3 4 5])))
  (is (every? #(= 2 (:stages %)) (mapv ms/adams-pece [2 3 4 5])))
  (is (= 7 (:stages rk/dopri54)) "against seven for a comparable single-step method"))

(deftest the-starting-procedure-does-not-cap-the-order
  ;; A startup one order too coarse silently limits everything after it, and
  ;; would show up as a fifth-order method converging at four.
  (let [errs (mapv (fn [h]
                     (abs (- (first (:y (core/step-until
                                         (ms/integrator (ms/adams-pece 5) exponential 0.0 [1.0] h) 1.0)))
                             target)))
                   [0.05 0.025 0.0125])]
    (doseq [p (ratio-order errs)]
      (is (> p 4.4) (str "order-5 Adams observed " p)))))
