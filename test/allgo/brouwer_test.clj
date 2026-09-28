(ns allgo.brouwer-test
  "Brouwer's theory against J2's motion integrated numerically, and his
  second-order secular rates against the mean Hamiltonian averaged from
  first principles. With the rates and the mean L both second order, what
  thirty orbits leave near-circular orbits is third order in J2 -- a
  tenth of J2, a thousandth of the error -- and a few meters; eccentric
  orbits keep a drift of order J2^2 e from the long-period terms left
  out. Spacetrack Report No. 3's shortened terms are kilometers out."
  (:require [allgo.astro.brouwer :as br]
            [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.perturbations :as pt]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.gauss-jackson :as gj]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(def ^:private R c/R-earth)

(defn- j2-acceleration [J2]
  (fn [_ r] (geo/acceleration {:GM mu :R R :normalized? false :C {[0 0] 1.0 [2 0] (- J2)} :S {}} r 2)))

(defn- integrated
  "J2's motion from state `s0`, `orbits` periods of a, by Gauss-Jackson
  at a step fine enough for the perigee: `[t r]` at the end."
  [J2 a e s0 orbits]
  (let [k (long (/ 200 (math/pow (- 1.0 e) 2)))
        [t r] (peek (gj/integrate (j2-acceleration J2) 0.0 (first s0) (second s0)
                                  (/ (kep/period mu a) k) (* k orbits)))]
    [t r]))

(defn- error-after
  "How far J2's motion, integrated from the osculating state of mean
  elements `el`, is from the theory's after `orbits` orbits, km."
  [J2 {:keys [a e] :as el} orbits]
  (let [[t r] (integrated J2 a e (br/osculating-state mu R J2 el) orbits)]
    (v3/distance r (first (br/propagate mu R J2 el t)))))

(def ^:private orbits
  [{:a 7000.0 :e 0.001 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7000.0 :e 0.01 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 8000.0 :e 0.1 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 14000.0 :e 0.5 :i 0.9 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7000.0 :e 0.01 :i 0.001 :raan 0.3 :argp 1.0 :M 0.2}
   {:a 7200.0 :e 0.05 :i 2.8 :raan 0.3 :argp 1.0 :M 0.2}])

(deftest semi-major-axis
  (testing "the first-order change in a is the closed form n dW1/dl gives"
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

;; ------------------------------------------------ second-order secular rates

(defn- disturbing
  "J2's Hamiltonian at Delaunay [L G H l g]."
  [k2 [L G H l g]]
  (let [a (/ (* L L) mu) e (/ (math/sqrt (* (- L G) (+ L G))) L) th2 (math/pow (/ H G) 2)
        E (kep/kepler-equation l e) f (kep/eccentric->true E e) r (* a (- 1.0 (* e (math/cos E))))]
    (* (/ (* mu k2) (* r r r)) (- (* 0.5 (- 1.0 (* 3.0 th2))) (* 1.5 (- 1.0 th2) (math/cos (* 2.0 (+ f g))))))))

(defn- second-order-hamiltonian
  "The mean Hamiltonian's second-order secular part at mean momenta [L G
  H]: 1/2 F0_LL W1_l^2 + F1_L W1_l + F1_G W1_g averaged over l and g, by
  the trapezoid rule, every partial by differences."
  [k2 [L G H]]
  (let [h (min (* 1e-6 L) (* 0.25 (- L G)))
        nl 48 ng 8
        term (fn [l g]
               (let [[[_ _ _ Wl Wg]] (diff/jacobian (fn [x] [(br/generating-function mu k2 x)]) [L G H l g]
                                                    {:steps [h h (* 1e-6 G) 1e-6 1e-6]})
                     [[FL FG]] (diff/jacobian (fn [[L G]] [(disturbing k2 [L G H l g])]) [L G] {:steps [h h]})]
                 (+ (* -1.5 (/ (* mu mu) (math/pow L 4)) Wl Wl) (* FL Wl) (* FG Wg))))]
    (/ (reduce + (for [i (range nl) j (range ng)] (term (* 2 math/PI (/ i nl)) (* 2 math/PI (/ j ng)))))
       (* nl ng))))

(deftest secular-rates
  (testing "the closed forms are the derivatives of the second-order mean Hamiltonian, averaged from first principles"
    (doseq [{:keys [a e i] :as el} [(orbits 2) (orbits 5)]]
      (let [k2 (* 0.5 geo/J2 R R)
            [L G H] (br/delaunay el)
            dLG (min (* 1e-4 L) (* 0.2 (- L G)))
            [[l2 g2 h2]] (diff/jacobian (fn [x] [(second-order-hamiltonian k2 x)]) [L G H]
                                        {:steps [dLG dLG (* 1e-4 G)] :richardson? false})
            closed (br/secular-rates a e i)
            first-order (pt/j2-secular a e i)
            second (fn [k] (- (closed k) (first-order k)))]
        (is (< (abs (- (second :M) l2)) (* 1e-4 (abs l2))) (str el))
        (is (< (abs (- (second :argp) g2)) (* 1e-4 (abs g2))) (str el))
        (is (< (abs (- (second :raan) h2)) (* 1e-4 (abs h2))) (str el)))))
  (testing "for a circular orbit they are SGP4's J2^2 terms"
    (let [a 7000.0 i 1.1
          n (kep/mean-motion mu a)
          th (math/cos i) th2 (* th th)
          k (* n geo/J2 geo/J2 (math/pow (/ R a) 4))
          closed (br/secular-rates a 0.0 i)
          first-order (pt/j2-secular a 0.0 i)
          second (fn [key] (- (closed key) (first-order key)))]
      (is (< (abs (- (second :M) (* (/ 3.0 64.0) k (+ 13.0 (* -78.0 th2) (* 137.0 th2 th2))))) (* 1e-15 n)))
      (is (< (abs (- (second :argp) (* (/ 3.0 64.0) k (+ 7.0 (* -114.0 th2) (* 395.0 th2 th2))))) (* 1e-15 n)))
      (is (< (abs (- (second :raan) (* 0.375 k th (- 4.0 (* 19.0 th2))))) (* 1e-15 n))))))

;; ------------------------------------------------------------ propagation

(deftest against-numerical-j2
  (testing "thirty orbits: meters near-circular, tens of meters eccentric"
    (doseq [[el bound] (map vector orbits [0.005 0.005 0.04 0.05 0.05 0.04])]
      (is (< (error-after geo/J2 el 30) bound) (str el))))
  (testing "near-circular, what is left is third order: a tenth of J2 leaves a thousandth of the error, or near it"
    (doseq [el [(orbits 0) (orbits 4)]]
      (let [full (error-after geo/J2 el 30)
            tenth (error-after (* 0.1 geo/J2) el 30)]
        (is (> (/ full tenth) 400.0) (str el " " full " " tenth)))))
  (testing "the second-order rates are what does it: J2's first-order rates alone leave kilometers"
    (let [{:keys [a e i raan argp M] :as el} (orbits 1)
          [t r] (integrated geo/J2 a e (br/osculating-state el) 30)
          rates (pt/j2-secular a e i)
          first-order (br/osculating-state (assoc el :raan (+ raan (* (:raan rates) t))
                                                  :argp (+ argp (* (:argp rates) t))
                                                  :M (+ M (* (:M rates) t))))]
      (is (> (v3/distance r (first first-order)) 0.5))))
  (testing "and the report's shortened terms, missing J2 e, are kilometers out at e = 0.1"
    (let [{:keys [a e] :as el} (orbits 2)
          [t r] (integrated geo/J2 a e (pt/j2-osculating el) 30)]
      (is (> (v3/distance r (first (pt/j2-propagate el t))) 50.0)))))

(deftest osculating-to-mean
  (testing "mean to osculating and back"
    (doseq [{:keys [a e i raan argp M] :as el} orbits]
      (let [m (br/mean-elements (br/osculating-state el))]
        (is (< (abs (- (:a m) a)) 1e-8) (str el))
        (is (< (abs (- (:e m) e)) 1e-10) (str el))
        (is (< (abs (- (:i m) i)) 1e-11) (str el))
        (is (< (abs (am/wrap-angle (- (:raan m) raan))) 1e-11) (str el))
        (is (< (abs (am/wrap-angle (- (+ (:argp m) (:M m)) argp M))) 1e-11) (str el))))))
