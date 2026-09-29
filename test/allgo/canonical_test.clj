(ns allgo.canonical-test
  (:require [allgo.astro.canonical :as cu]
            [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest earth
  (let [u (cu/units :earth)]
    (testing "the TU and speed unit as quoted, 806.8 s and 7.905 km/s"
      (is (< (abs (- (:tu u) 806.81)) 0.01))
      (is (< (abs (- (:velocity u) 7.9054)) 1e-4)))
    (testing "an orbit of one DU takes 2 pi TU"
      (is (< (abs (- (kepler/period c/GM-earth c/R-earth) (* 2.0 math/PI (:tu u)))) 1e-9)))))

(deftest sun
  (testing "a TU is 1/k days, k the Gaussian constant, to the eight
            figures GM-sun and the astronomical unit now agree with it"
    (let [u (cu/units :sun)]
      (is (< (abs (- (/ (:tu u) 86400.0) (/ 1.0 c/gaussian-k))) (* 1e-8 (/ 1.0 c/gaussian-k)))))))

(deftest round-trip
  (let [u (cu/units :earth) r [7000.0 -1200.5 300.0] v [1.0 7.2 -0.5]]
    (is (every? #(< (abs %) 1e-12) (map - (cu/<-canonical u :length (cu/->canonical u :length r)) r)))
    (is (every? #(< (abs %) 1e-12) (map - (cu/<-canonical u :velocity (cu/->canonical u :velocity v)) v)))
    (is (< (abs (- (cu/->canonical u :time (:tu u)) 1.0)) 1e-15))))
