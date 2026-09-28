(ns allgo.vallado-covariance-test
  "Vallado's covariance transformations, checked three ways: the
  numerical Jacobian against one known in closed form, each
  transformation against its inverse, and the transformed covariance
  against the sample covariance of a cloud of states converted one by
  one -- the definition it approximates."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.covariance :as cov]
            [allgo.astro.kepler :as kep]
            [allgo.astro.states :as st]
            [allgo.math :as am]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private s (kep/elements->state c/GM-earth {:a 8000.0 :e 0.12 :i 0.9 :raan 0.6 :argp 1.4 :M 2.2}))

(def ^:private P
  "A plausible Cartesian covariance: tens of meters and centimeters per
  second, correlated."
  (let [L [[0.030 0 0 0 0 0] [0.010 0.020 0 0 0 0] [-0.005 0.004 0.025 0 0 0]
           [1e-5 -2e-5 3e-6 3e-5 0 0] [-4e-6 1e-5 2e-6 5e-6 2e-5 0] [2e-6 3e-6 -1e-5 -2e-6 4e-6 2.5e-5]]]
    (lin/mat-mul L (lin/transpose L))))

(defn- max-rel
  "The largest difference between two matrices, relative to each entry's
  natural scale sqrt(A_ii A_jj)."
  [A B]
  (apply max (for [i (range 6) j (range 6)]
               (/ (abs (- (get-in A [i j]) (get-in B [i j])))
                  (math/sqrt (* (get-in B [i i]) (get-in B [j j])))))))

(deftest round-trips
  (testing "classical and back"
    (let [{:keys [a e i raan argp M]} (kep/state->elements c/GM-earth (first s) (second s))
          Pc (cov/cartesian->classical P s)]
      (is (< (max-rel (cov/classical->cartesian Pc [a e i raan argp M]) P) 1e-7))))
  (testing "equinoctial and back"
    (let [eq (st/state->equinoctial s)]
      (is (< (max-rel (cov/equinoctial->cartesian (cov/cartesian->equinoctial P s) eq) P) 1e-7))))
  (testing "RSW and NTW are rotations: back exactly, and the trace kept"
    (doseq [[to from] [[cov/cartesian->rsw cov/rsw->cartesian] [cov/cartesian->ntw cov/ntw->cartesian]]]
      (let [P' (to P s)]
        (is (< (max-rel (from P' s) P) 1e-12))
        (is (< (abs (- (reduce + (map #(get-in P' [% %]) (range 3)))
                       (reduce + (map #(get-in P [% %]) (range 3)))))
               1e-15)))))
  (testing "the radial variance is the covariance seen along r"
    (let [R (lin/normalize (first s))
          Pr (mapv #(subvec % 0 3) (subvec P 0 3))]
      (is (< (abs (- (get-in (cov/cartesian->rsw P s) [0 0]) (lin/dot R (lin/mat-vec Pr R)))) 1e-15)))))

(defn- sample-covariance [xs]
  (let [n (count xs)
        m (lin/scale (reduce lin/add xs) (/ 1.0 n))
        ds (map #(lin/sub % m) xs)]
    (lin/mat-scale (reduce lin/mat-add (map (fn [d] (mapv (fn [a] (mapv #(* a %) d)) d)) ds)) (/ 1.0 (dec n)))))

(deftest monte-carlo
  (testing "a cloud of states, each converted, has the covariance the transformation gives"
    (let [rng (java.util.Random. 20260927)
          L (lin/cholesky P)
          x0 (vec (concat (first s) (second s)))
          n 20000
          draws (repeatedly n #(lin/add x0 (lin/mat-vec L (vec (repeatedly 6 (fn [] (.nextGaussian rng)))))))
          center (cov/classical-vector x0)
          ;; differences from the center, the angles wrapped
          elements (map (fn [x] (vec (map-indexed (fn [k [a b]] (if (>= k 2) (+ b (am/wrap-angle (- a b))) a))
                                                  (map vector (cov/classical-vector x) center))))
                        draws)
          sampled (sample-covariance elements)
          transformed (cov/cartesian->classical P [(subvec x0 0 3) (subvec x0 3 6)])]
      ;; a sample of 20000 pins a covariance to a percent or two
      (is (< (max-rel sampled transformed) 0.03)))))

(deftest flight-elements
  (testing "the radius's variance is the position covariance seen along r, whatever the Earth's turn"
    (let [Pf (cov/cartesian->flight P s 60580.0 60580.0 {})
          R (lin/normalize (first s))
          Pr (mapv #(subvec % 0 3) (subvec P 0 3))]
      (is (< (abs (- (get-in Pf [0 0]) (lin/dot R (lin/mat-vec Pr R)))) (* 1e-8 (get-in Pf [0 0])))))))

(deftest analytic-jacobians
  (testing "the analytic partials of the state in the classical elements are the numerical ones"
    (doseq [el [[8000.0 0.12 0.9 0.6 1.4 2.2] [7000.0 0.001 1.7 5.9 0.2 4.0] [26000.0 0.7 0.3 2.0 4.5 0.1]]]
      (let [A (cov/classical-partials el)
            ;; steps no smaller than 1e-7, below which the differences are noise
            N (diff/jacobian #(cov/classical->flat c/GM-earth %) el
                             {:steps (mapv #(* 1e-6 (max 0.1 (abs %))) el)})]
        (is (every? true? (for [i (range 6) j (range 6)]
                            (< (abs (- (get-in A [i j]) (get-in N [i j])))
                               (* 1e-7 (max 1e-3 (abs (get-in N [i j])))))))
            (str el)))))
  (testing "and their inverse is the numerical Jacobian of the elements in the state"
    (let [x (vec (concat (first s) (second s)))
          J (lin/inverse-general (cov/classical-partials (cov/classical-vector x)))
          N (diff/jacobian cov/classical-vector x {:angles #{2 3 4 5}})]
      (is (every? true? (for [i (range 6) j (range 6)]
                          (< (abs (- (get-in J [i j]) (get-in N [i j])))
                             (* 1e-7 (max 1e-6 (abs (get-in N [i j]))))))))))
  (testing "so the covariance carried either way round comes back"
    (let [el (cov/classical-vector (vec (concat (first s) (second s))))]
      (is (< (max-rel (cov/classical->cartesian (cov/cartesian->classical P s) el) P) 1e-9))))
  (testing "the analytic partials of the state in the equinoctial elements are the numerical ones, column by column"
    (doseq [el [{:a 8000.0 :e 0.12 :i 0.9 :raan 0.6 :argp 1.4 :M 2.2}
                {:a 7000.0 :e 1e-4 :i 1e-4 :raan 0.3 :argp 1.0 :M 2.2}
                {:a 26000.0 :e 0.7 :i 1.1 :raan 4.0 :argp 5.0 :M 3.0}
                {:a 7200.0 :e 0.05 :i 2.8 :raan 0.3 :argp 1.0 :M 5.2}]]
      (let [eq (st/state->equinoctial (kep/elements->state c/GM-earth el))
            x (mapv eq [:a :af :ag :chi :psi :meanlon])
            f (fn [[a af ag chi psi meanlon]]
                (let [[r v] (st/equinoctial->state {:a a :af af :ag ag :chi chi :psi psi :meanlon meanlon :fr (:fr eq)})]
                  (vec (concat r v))))
            A (cov/equinoctial-partials eq)
            ;; the numerical side goes by way of the classical elements,
            ;; whose round-off near e = i = 0 wants steps no smaller than 1e-8
            N (diff/jacobian f x {:steps (mapv #(* 1e-6 (max 0.01 (abs %))) x)})
            column (fn [M j] (mapv #(nth % j) M))]
        (doseq [j (range 6)]
          (let [a (column A j) n (column N j)]
            (is (< (lin/distance a n) (* 1e-7 (lin/length n))) (str el " column " j)))))))
  (testing "and where the classical elements fail -- circular and equatorial -- the covariance still comes back"
    (let [s0 (kep/elements->state c/GM-earth {:a 7000.0 :e 0.0 :i 0.0 :raan 0.0 :argp 0.0 :M 1.0})]
      (is (< (max-rel (cov/equinoctial->cartesian (cov/cartesian->equinoctial P s0) (st/state->equinoctial s0)) P)
             1e-9)))))
