(ns allgo.jacchia-test
  "Jacchia 1971 against the report's own tables (SAO Special Report 332):
  the inflection gradients of section 5 and rows of Table 6 at
  exospheric temperatures of 600, 800 and 1800 K for the static models;
  Table 1 for the diurnal variation, Tables 3 and 4 for the semiannual
  and seasonal-latitudinal ones, Table 5 for helium."
  (:require [allgo.astro.jacchia :as j]
            [allgo.astro.solar :as solar]
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

(deftest diurnal
  (testing "Table 1: Tl/Tc x 1000 by local solar time, the Sun at +23.44 degrees"
    (let [dec (math/to-radians 23.44)
          rows {90.0 (repeat 24 1202)
                45.0 [1093 1087 1085 1085 1085 1089 1099 1116 1141 1172 1206 1238
                      1264 1282 1288 1284 1270 1248 1222 1194 1167 1142 1121 1104]
                0.0 [1021 1013 1009 1009 1010 1015 1028 1052 1086 1128 1174 1218
                     1254 1277 1286 1280 1261 1232 1196 1158 1121 1087 1058 1036]
                -45.0 [1016 1010 1008 1008 1008 1012 1021 1037 1060 1089 1120 1150
                       1175 1191 1197 1193 1180 1160 1136 1110 1084 1061 1041 1026]
                -90.0 (repeat 24 1080)}]
      (doseq [[lat row] rows
              [lst printed] (map vector (range) row)]
        (is (<= (abs (- (* 1000.0 (j/local-temperature 1.0 (math/to-radians lat) dec lst)) printed)) 0.51)
            (str lat " " lst))))))

(def ^:private mjd-1971 40952.0)                            ; 1971 January 1, 0h

(deftest semiannual
  (testing "Table 3a: f(z)"
    (doseq [[z f] [[100.0 0.068] [200.0 0.112] [400.0 0.237] [800.0 0.353] [1200.0 0.285]]]
      (is (<= (abs (- (j/semiannual-height z) f)) 0.0005) (str z))))
  (testing "Table 3b: g(t), the phase counted from January 1"
    (doseq [[doy g] [[1 -0.145] [91 0.354] [201 -0.515] [301 0.478]]]
      (is (<= (abs (- (j/semiannual-time (+ 36204.0 (dec doy))) g)) 0.0015) (str doy)))))

(deftest seasonal-latitudinal
  (testing "Table 4: the half-range S at the pole, and the phase P"
    (let [t0 36204.0
          pole (/ math/PI 2)
          ;; the date where P = 1, so the variation at the pole is S itself
          tp (+ t0 (* 365.2422 (/ (- (/ math/PI 2) 1.72) (* 2 math/PI))))]
      (doseq [[z s] [[95.0 0.068] [110.0 0.166] [130.0 0.070] [150.0 0.008]]]
        (is (<= (abs (- (j/seasonal-latitudinal z pole tp) s)) 0.0006) (str z)))
      ;; a calendar date's phase shifts by some 0.003 with the year, and
      ;; the table does not say which it took
      (doseq [[doy p] [[1 0.989] [191 -0.961] [101 -0.297]]]
        (is (<= (abs (- (/ (j/seasonal-latitudinal 110.0 pole (+ t0 (dec doy)))
                           (j/seasonal-latitudinal 110.0 pole tp))
                        p))
                0.003)
            (str doy))))))

(deftest helium
  (testing "Table 5: Delta log10 n(He), the Sun's declination on the day in 1971"
    (doseq [[day row] [[0 {-90.0 -0.226 -60.0 -0.215 -30.0 -0.146 0.0 0.0 30.0 0.189 60.0 0.350 90.0 0.413}]
                       [170 {-90.0 0.420 -60.0 0.356 -30.0 0.192 30.0 -0.148 60.0 -0.218 90.0 -0.230}]
                       [90 {-90.0 0.076 90.0 -0.042}]]
            [lat printed] row]
      (let [dec (second (solar/equatorial-low (+ mjd-1971 day)))]
        (is (<= (abs (- (j/helium (math/to-radians lat) dec) printed)) 0.003) (str day " " lat))))))

(deftest geomagnetic
  (testing "equations (18) and (20)"
    (is (< (abs (- (j/geomagnetic-temperature 4.0 400.0) (+ 112.0 (* 0.03 (math/exp 4.0))))) 1e-12))
    (is (< (abs (- (j/geomagnetic-temperature 4.0 150.0) (+ 56.0 (* 0.02 (math/exp 4.0))))) 1e-12))
    (is (zero? (j/geomagnetic-density 4.0 400.0)))))

(deftest atmosphere
  (testing "with every variation neutral, the static model at the temperature the inputs give"
    (let [in {:mjd 41000.0 :alt 400.0 :lat 0.0 :lst 3.0
              :f107 150.0 :f107a 150.0 :kp 0.0 :dec 0.0}
          r (j/atmosphere in)
          tinf (+ (j/local-temperature (j/night-minimum 150.0 150.0) 0.0 0.0 3.0) 0.03)
          s (j/static tinf 400.0)
          f (math/pow 10.0 (j/semiannual 400.0 (:mjd in)))]
      (is (< (abs (- (:t-exo r) tinf)) 1e-9))
      (is (< (abs (- (:rho r) (* 1e3 (:rho s) f))) (* 1e-9 (:rho r))))
      (is (< (abs (- (:N2 r) (* 1e6 (get-in s [:n :N2]) f))) (* 1e-9 (:N2 r))))))
  (testing "the density rises with the flux and with geomagnetic activity"
    (let [base {:mjd 41000.0 :alt 400.0 :lat 30.0 :lst 14.0 :f107 150.0 :f107a 150.0 :kp 2.0}
          rho #(:rho (j/atmosphere (merge base %)))]
      (is (> (rho {:f107a 200.0 :f107 200.0}) (rho {})))
      (is (> (rho {:kp 6.0}) (rho {}))))))
