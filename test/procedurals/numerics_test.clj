(ns procedurals.numerics-test
  (:require [procedurals.numerics :as num]
            [procedurals.numerics.core :as core]
            [procedurals.numerics.extrapolation :as ex]
            [procedurals.numerics.multistep :as ms]
            [procedurals.numerics.rk :as rk]
            [procedurals.numerics.rkn :as rkn]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (< (abs (double (- a b))) tol))

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

;; ---------------------------------------------------------- extrapolation

(deftest modified-midpoint-is-second-order-on-its-own
  ;; Before extrapolation it is a crude rule; the point is the shape of its
  ;; error, not its accuracy.
  (let [errs (mapv (fn [n] (abs (- (first (ex/modified-midpoint exponential 0.0 [1.0] 1.0 n)) target)))
                   [8 16 32 64])]
    (doseq [p (ratio-order errs)]
      (is (< 1.7 p 2.3) (str "observed " p)))))

(deftest extrapolation-reaches-twice-the-levels-in-order
  ;; The payoff of an error expansion in even powers of h: each level buys
  ;; two orders, not one.
  (doseq [[levels expected] [[2 4] [3 6] [4 8]]]
    (testing (str levels " levels")
      (let [method (ex/gbs levels)
            errs   (mapv (fn [h]
                           (max 1e-16 (abs (- (first (:y (core/step-until
                                                          (ex/integrator method exponential 0.0 [1.0] h) 1.0)))
                                              target))))
                         [0.5 0.25 0.125])]
        (is (= expected (:order method)))
        (doseq [p (ratio-order errs)]
          (is (< (- expected 0.6) p (+ expected 0.6)) (str "observed " p)))))))

(deftest extrapolation-beats-the-rule-it-is-built-from
  ;; One macro-step across the whole interval, both spending 20 sub-steps.
  (let [raw (abs (- (first (ex/modified-midpoint exponential 0.0 [1.0] 1.0 20)) target))
        ext (abs (- (first (:y (core/step-until (ex/integrator (ex/gbs 4) exponential 0.0 [1.0] 1.0) 1.0))) target))]
    (is (< ext (/ raw 100.0))
        (str "extrapolated " ext " should be orders better than the raw midpoint " raw))))

(deftest extrapolation-handles-a-system
  (let [f (fn [_ [x v]] [v (- x)])
        s (core/step-until (ex/integrator (ex/gbs 4) f 0.0 [1.0 0.0] 0.25) (* 2 Math/PI))
        [x v] (:y s)]
    (is (< (abs (- x 1.0)) 1e-9))
    (is (< (abs v) 1e-9))))

(deftest extrapolation-adapts-its-step
  (doseq [tol [1e-8 1e-10]]
    (let [s (core/step-until (ex/integrator (ex/gbs 4) exponential 0.0 [1.0] 0.25
                                            {:adaptive? true :tol-abs tol :tol-rel tol}) 1.0)]
      (is (< (abs (- (first (:y s)) target)) (* 500 tol)) (str "tol " tol)))))

;; ------------------------------------------------------------------ Kepler

(def ^:private mu 1.0)

(defn- energy-error [a e steps-per-orbit orbits method]
  (let [[r0 v0] (num/periapsis-state mu a e)
        period  (num/orbital-period mu a)
        e0      (num/specific-energy mu r0 v0)
        s       (core/step-until (num/integrator-2 method (num/kepler mu) 0.0 r0 v0
                                                   (/ period steps-per-orbit))
                                 (* orbits period))]
    (abs (/ (- (num/specific-energy mu (:y s) (:dy s)) e0) e0))))

(deftest orbital-elements-agree-with-the-closed-forms
  (doseq [a [0.5 1.0 3.0], e [0.0 0.3 0.7]]
    (let [[r v] (num/periapsis-state mu a e)]
      (is (< (abs (- (num/specific-energy mu r v) (/ (- mu) (* 2.0 a)))) 1e-12)
          "specific energy is -mu/2a")
      (is (< (abs (- (Math/sqrt (reduce + (map * (num/angular-momentum r v)
                                               (num/angular-momentum r v))))
                     (Math/sqrt (* mu a (- 1.0 (* e e))))))
             1e-12)
          "angular momentum is sqrt(mu a (1-e^2))"))))

(deftest an-orbit-closes-on-itself
  (doseq [e [0.0 0.3 0.7]]
    (let [[r0 v0] (num/periapsis-state mu 1.0 e)
          period  (num/orbital-period mu 1.0)
          s       (core/step-until (num/integrator-2 rkn/rkn4 (num/kepler mu) 0.0 r0 v0
                                                     (/ period 2000))
                                   period)
          drift   (Math/sqrt (reduce + (map (fn [x y] (let [d (- x y)] (* d d))) (:y s) r0)))]
      (is (< drift 1e-5) (str "e=" e " returned " drift " from where it started")))))

(deftest conserved-quantities-are-conserved
  (doseq [method num/second-order]
    (testing (:name method)
      (is (< (energy-error 1.0 0.3 400 1 method) 1e-4)
          (str (:name method) " loses energy over a single orbit")))))

(deftest symplectic-methods-hold-their-energy-over-long-runs
  ;; The reason a second-order method can beat a fourth over a long
  ;; propagation, and the whole argument for Verlet in orbit work: its
  ;; energy error oscillates within a bound instead of accumulating, while
  ;; a more accurate but non-symplectic method walks steadily away.
  (let [drift (fn [method orbits] (energy-error 1.0 0.5 400 orbits method))]
    (testing "Verlet's energy error stops growing"
      (let [at-100 (drift rkn/verlet 100)
            at-400 (drift rkn/verlet 400)]
        (is (< at-400 (* 1.5 at-100))
            (str "bounded: " at-100 " at 100 orbits, " at-400 " at 400"))))
    (testing "RKN4 is far more accurate early but drifts linearly"
      (let [at-10  (drift rkn/rkn4 10)
            at-100 (drift rkn/rkn4 100)
            at-400 (drift rkn/rkn4 400)]
        (is (< at-10 (drift rkn/verlet 10)) "more accurate over ten orbits")
        (is (> (/ at-100 at-10) 5.0) "and grows roughly in proportion to time")
        (is (> (/ at-400 at-100) 2.0))))
    (testing "so over four hundred orbits the gap has closed to under two decades"
      (is (< (/ (drift rkn/verlet 400) (drift rkn/rkn4 400)) 200.0)))))

;; ------------------------------------------- second-order extrapolation

(deftest stoermer-extrapolation-reaches-twice-the-levels-in-order
  ;; Checked on Kepler rather than an oscillator, and against the exact
  ;; answer rather than a fine-stepped reference: after exactly one period
  ;; the orbit closes on its starting state, so the truth is known to the
  ;; last bit. A reference produced by the same family plateaus at its own
  ;; error and makes every high-order method look like it stopped
  ;; converging around 1e-8.
  (let [[r0 v0] (num/periapsis-state mu 1.0 0.6)
        period  (num/orbital-period mu 1.0)
        f       (num/kepler mu)
        err     (fn [method steps]
                  (max 1e-16
                       (Math/sqrt (reduce + (map (fn [x y] (let [d (- x y)] (* d d)))
                                                 (:y (core/step-until
                                                      (num/integrator-2 method f 0.0 r0 v0 (/ period steps))
                                                      period))
                                                 r0)))))]
    (doseq [[levels expected] [[2 4] [3 6] [4 8]]]
      (testing (str levels " levels")
        (let [method (ex/gbs-2 levels)
              orders (ratio-order (mapv #(err method %) [12 24 48 96]))]
          (is (= expected (:order method)))
          ;; The coarsest ratio is pre-asymptotic and the finest can touch
          ;; the roundoff floor, so it is the middle that must land.
          (is (< (- expected 1.0) (second orders) (+ expected 1.0))
              (str (:name method) " observed " orders)))))))

(deftest second-order-extrapolation-is-offered-for-kepler
  (is (some #(= :extrapolation-2 (:kind %)) num/second-order)
      "the Kepler demo draws its menu from this list")
  (doseq [method num/second-order]
    (testing (:name method)
      (is (< (energy-error 1.0 0.3 400 1 method) 1e-4)
          "every offered method holds its energy over an orbit"))))

(deftest extrapolation-integrates-second-order-systems-natively
  ;; Not by doubling the system: the Stoermer rule never forms a velocity
  ;; along the way, taking one at the end from the last difference.
  (let [[y dy] (ex/stoermer-midpoint (num/harmonic 1.0) 0.0 [1.0] [0.0] 1.0 64)]
    (is (< (abs (- (first y) (Math/cos 1.0))) 1e-4))
    (is (< (abs (- (first dy) (- (Math/sin 1.0)))) 1e-4))))

;; -------------------------------------------------- variable-step multistep

(deftest variable-coefficients-reduce-to-the-fixed-tables
  ;; The generalisation must contain the special case exactly, not
  ;; approximately: evenly spaced offsets are just particular node values.
  (doseq [k [2 3 4 5 6]]
    (testing (str "order " k)
      (doseq [[a b] (map vector (ms/adams-bashforth-coefficients k)
                         (ms/variable-coefficients (mapv #(- (double %)) (range k))))]
        (is (close? a b 1e-14) "explicit"))
      (doseq [[a b] (map vector (ms/adams-moulton-coefficients k)
                         (ms/variable-coefficients
                          (vec (cons 1.0 (mapv #(- (double %)) (range (dec k)))))))]
        (is (close? a b 1e-14) "implicit")))))

(deftest variable-coefficients-integrate-a-constant-exactly
  ;; Whatever the spacing, a constant derivative must be integrated exactly,
  ;; which means the weights sum to one. It is the one condition that holds
  ;; for every node arrangement.
  (doseq [offs [[0.0 -1.0 -2.5] [0.0 -0.3 -0.9 -2.2] [0.0 -1.7 -1.9 -4.0 -7.1]
                [1.0 0.0 -0.4] [0.0 -5.0] [0.0 -0.01 -0.02]]]
    (is (close? 1.0 (reduce + (ms/variable-coefficients offs)) 1e-11) (str offs))))

(deftest variable-step-adams-honours-its-tolerance
  (doseq [method ms/catalog-variable
          tol [1e-6 1e-9 1e-12]]
    (let [s (core/step-until (ms/variable-integrator method exponential 0.0 [1.0] 0.05
                                                     {:tol-abs tol :tol-rel tol})
                             1.0)]
      (is (< (abs (- (first (:y s)) target)) (* 500 tol))
          (str (:name method) " at tol " tol)))))

(deftest the-step-follows-the-problem
  ;; The reason for the machinery. On an eccentric orbit the satellite races
  ;; through periapsis and loafs at apoapsis, and a step sized for the former
  ;; is wasted on the latter.
  (let [mu 398600.4415
        a  26000.0
        ecc 0.72
        rp (* a (- 1.0 ecc))
        vp (Math/sqrt (/ (* mu (+ 1.0 ecc)) rp))
        f  (fn [_ y] (let [r (subvec (vec y) 0 3)
                           d (Math/sqrt (reduce + (map * r r)))]
                       (into (subvec (vec y) 3 6)
                             (mapv #(* (- (/ mu (* d d d))) %) r))))
        period (* 2.0 Math/PI (Math/sqrt (/ (* a a a) mu)))
        integ  (ms/variable-integrator (ms/adams-variable 5) f 0.0
                                       [rp 0.0 0.0 0.0 vp 0.0] 5.0
                                       {:tol-abs 1e-10 :tol-rel 1e-10})
        traj   (vec (take-while #(< (:t %) period) (core/trajectory integ)))
        radius (fn [s] (Math/sqrt (reduce + (map * (subvec (vec (:y s)) 0 3)
                                                 (subvec (vec (:y s)) 0 3)))))]
    (is (seq traj))
    (testing "the step varies by a large factor over one revolution"
      (is (> (/ (apply max (map :h traj)) (apply min (map :h traj))) 10.0)))
    (testing "and it is small where the motion is fast"
      (is (< (radius (apply min-key :h traj)) (* 0.5 (radius (apply max-key :h traj))))
          "tightest step near periapsis, loosest near apoapsis"))))

(deftest variable-step-beats-fixed-at-equal-work
  ;; Same number of steps, spent where they matter.
  (let [mu 398600.4415, a 26000.0, ecc 0.72
        rp (* a (- 1.0 ecc))
        vp (Math/sqrt (/ (* mu (+ 1.0 ecc)) rp))
        y0 [rp 0.0 0.0 0.0 vp 0.0]
        f  (fn [_ y] (let [r (subvec (vec y) 0 3)
                           d (Math/sqrt (reduce + (map * r r)))]
                       (into (subvec (vec y) 3 6)
                             (mapv #(* (- (/ mu (* d d d))) %) r))))
        period (* 2.0 Math/PI (Math/sqrt (/ (* a a a) mu)))
        err (fn [s] (Math/sqrt (reduce + (map (fn [p q] (let [d (- p q)] (* d d)))
                                              (subvec (vec (:y s)) 0 3) (subvec y0 0 3)))))
        variable (core/step-until (ms/variable-integrator (ms/adams-variable 5) f 0.0 y0 5.0
                                                          {:tol-abs 1e-10 :tol-rel 1e-10})
                                  period)
        fixed    (core/step-until (ms/integrator (ms/adams-pece 5) f 0.0 y0 (/ period 800)) period)]
    (is (< (err variable) (err fixed))
        (str "variable " (err variable) " km vs fixed " (err fixed) " km at comparable step counts"))
    (is (< (err variable) (* 0.1 (err fixed))) "and by a wide margin")))
