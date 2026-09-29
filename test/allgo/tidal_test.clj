(ns allgo.tidal-test
  "The tide the Moon raises on the Earth, acting on a satellite (Vallado,
  section 9.5.3): the acceleration against the potential's gradient, the
  averaged potential against the potential averaged round the orbit
  numerically, and its secular rates through Lagrange's equations against
  the acceleration put through Gauss's equations and averaged."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kep]
            [allgo.astro.perturbations :as pt]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private k2 0.30)
(def ^:private R c/R-earth)
(def ^:private r3 384400.0)
(def ^:private n3 (v3/normalize [0.6 -0.7 0.3]))
(def ^:private moon (v3/scale n3 r3))
(def ^:private el {:a 7000.0 :e 0.05 :i 0.9 :raan 0.4 :argp 1.1 :M 0.0})

(defn- potential [r]
  (let [rm (v3/length r) cs (/ (v3/dot r n3) rm)]
    (* k2 (/ c/GM-moon r3) (math/pow (/ R r3) 2.0) (math/pow (/ R rm) 3.0) 0.5 (- (* 3.0 cs cs) 1.0))))

(deftest acceleration-is-the-gradient
  (let [r [5000.0 -3000.0 4200.0] h 1e-3
        grad (mapv (fn [k] (/ (- (potential (update r k + h)) (potential (update r k - h))) (* 2.0 h))) (range 3))
        acc (pt/tidal-acceleration k2 R c/GM-moon moon r)]
    (is (every? #(< (abs %) (* 1e-6 (v3/length acc))) (map - acc grad)))))

(deftest averaged-potential
  (testing "the closed form against the potential averaged over the mean anomaly"
    (let [n 2000
          avg (/ (reduce + (for [k (range n)]
                             (potential (first (kep/elements->state c/GM-earth (assoc el :M (* 2.0 math/PI (/ k n))))))))
                 n)]
      (is (< (abs (- avg (pt/tidal-averaged-potential el k2 R c/GM-moon r3 n3))) (* 1e-9 (abs avg)))))))

(deftest secular-rates
  (testing "Lagrange on the averaged potential against Gauss on the acceleration, averaged"
    (let [lagrange (pt/tidal-secular el k2 R c/GM-moon r3 n3)
          gauss (pt/averaged-rates c/GM-earth el (fn [r _] (pt/tidal-acceleration k2 R c/GM-moon moon r)) 720)]
      (doseq [k [:i :raan :argp]]
        (is (< (abs (- (k lagrange) (k gauss))) (* 1e-3 (abs (k gauss)))) (name k)))
      (testing "the eccentricity steady: averaged, the tide does not depend on the perigee"
        (is (< (abs (:e lagrange)) (* 1e-12 (abs (:raan gauss)))))
        (is (< (abs (:e gauss)) (* 1e-12 (abs (:raan gauss))))))
      (testing "and no secular change in the semi-major axis"
        (is (< (abs (:a gauss)) (* 1e-9 (:a el) (abs (:raan gauss)))))))))
