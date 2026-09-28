(ns allgo.jacchia-test
  "Jacchia 1971's static models against the report's own tables (SAO
  Special Report 332): the inflection gradients of section 5 and rows of
  Table 6 at exospheric temperatures of 600, 800 and 1800 K."
  (:require [allgo.astro.jacchia :as j]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest inflection-gradient
  (testing "Gx as the report tabulates it, K/km"
    (doseq [[tinf gx] [[500.0 6.23] [1000.0 11.22] [1600.0 14.24] [2000.0 15.66]]]
      (is (< (abs (- (j/inflection-gradient tinf) gx)) 0.006) (str tinf)))))

(def ^:private table-6
  "[T-inf z T {species log10 n} rho] read from Table 6."
  [[800.0 100.0 192.6 {:N2 12.9533 :O2 12.2695 :O 12.0407 :Ar 11.0311 :He 7.8493} 5.528e-10]
   [800.0 400.0 795.6 {:N2 5.7907 :O2 4.1755 :O 7.6858 :Ar 1.0794 :He 6.5320} 1.341e-15]
   [800.0 1000.0 799.9 {:O 2.6813 :He 5.2791 :H 4.8490} 1.395e-18]
   [600.0 200.0 555.8 {:N2 8.8689 :O2 7.6698 :O 9.5096 :Ar 5.4037 :He 7.0459} 1.229e-13]
   [600.0 1000.0 600.0 {:O 0.2617 :He 4.7200 :H 5.7804} 1.358e-18]
   [1800.0 300.0 1702.0 {:N2 8.7539 :O2 7.6065 :O 9.2391 :Ar 5.4432 :He 6.7984} 7.467e-14]
   [1800.0 1000.0 1799.5 {:N2 4.0918 :O2 2.2845 :O 6.5660 :He 6.1207 :H 3.0459} 1.072e-16]])

(deftest static-models
  (testing "temperatures, number densities and mass densities of Table 6"
    (doseq [[tinf z t logs rho] table-6]
      (let [m (j/static tinf z)]
        (is (< (abs (- (:t m) t)) 0.06) (str tinf " " z " T"))
        (doseq [[k v] logs]
          (is (< (abs (- (math/log10 (get-in m [:n k])) v)) 3e-4) (str tinf " " z " " k)))
        (is (< (abs (- (:rho m) rho)) (* 1e-3 rho)) (str tinf " " z " rho"))))))
