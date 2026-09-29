(ns allgo.jb2008-test
  "JB2008, implemented from its papers, against what the papers fix --
  its pieces joining where the papers say they join, Jacchia 1970's
  helium extremes, the storm equations' behavior -- and against an
  independent implementation's densities (Orekit 12.2, test/data/jb2008)
  from 220 to 400 km."
  (:require [allgo.astro.jb2008 :as jb]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(deftest night-minimum
  (testing "with every index at its average, Tc = 392.4 + 3.227 F (JB2008 eq. 2), whatever S10"
    (doseq [f [70.0 150.0 240.0]]
      (is (< (abs (- (jb/night-minimum {:f10 f :f10b f :s10 f :s10b f :m10 f :m10b f :y10 f :y10b f})
                     (+ 392.4 (* 3.227 f))))
             1e-9))))
  (testing "the daily departures each add their own coefficient"
    (let [base {:f10 150.0 :f10b 150.0 :s10 150.0 :s10b 150.0 :m10 150.0 :m10b 150.0 :y10 150.0 :y10b 150.0}
          tc (jb/night-minimum base)]
      (doseq [[k c] [[:f10 0.298] [:s10 2.259] [:m10 0.312] [:y10 0.178]]]
        (is (< (abs (- (jb/night-minimum (update base k + 10.0)) tc (* 10.0 c))) 1e-9) (str k))))))

(deftest local-time-corrections
  (testing "JB2006 eqs. 7 and 8 meet at 250 km, at every local time, latitude and flux -- to the rounding of
            the printed coefficients; a misread one would leave degrees"
    (doseq [lst (range 0.0 24.0 1.5) lat [-80.0 -30.0 0.0 45.0 85.0] f [65.0 120.0 180.0 250.0]]
      (is (< (abs (- (jb/tc-correction 250.0 lst lat f) (jb/tc-correction (- 250.0 1e-9) lst lat f))) 1e-4)
          (str lst " " lat " " f))))
  (testing "and are tens of degrees, as the papers' figures show"
    (doseq [z [220.0 300.0 450.0 650.0] lst (range 0.0 24.0 3.0)]
      (is (< (abs (jb/tc-correction z lst 30.0 150.0)) 60.0))))
  (testing "the inflection correction, below 200 km only"
    (is (zero? (jb/tx-correction 200.0 10.0 150.0)))
    (is (< (abs (jb/tx-correction 170.0 10.0 150.0)) 30.0))))

(deftest high-altitude-factor
  (testing "1 up to 1000 km, and value and slope continuous at 1000 and 1500"
    (doseq [f [70.0 150.0 220.0]]
      (let [F #(jb/high-altitude-factor % f)
            slope (fn [z] (/ (- (F (+ z 1e-3)) (F (- z 1e-3))) 2e-3))]
        (is (== 1.0 (F 900.0)))
        (is (< (abs (- (F 1000.0) (F 1000.001))) 1e-6))
        (is (< (abs (slope 1000.0005)) 1e-5))
        (is (< (abs (- (F 1499.999) (F 1500.001))) 1e-5))
        (is (< (abs (- (slope 1499.99) (slope 1500.01))) 1e-5)))))
  (testing "the 1500-2000 km factors of JB2006's figure 15: about 1.6 and 2.1 at 70 sfu, 0.8 and 1.1 at 220"
    (is (< (abs (- (jb/high-altitude-factor 1500.0 70.0) 1.60)) 0.05))
    (is (< (abs (- (jb/high-altitude-factor 2000.0 70.0) 2.10)) 0.05))
    (is (< (abs (- (jb/high-altitude-factor 1500.0 220.0) 0.80)) 0.05))
    (is (< (abs (- (jb/high-altitude-factor 2000.0 220.0) 1.13)) 0.05))))

(deftest helium
  (testing "Jacchia 1970's extremes at the solstice: 2.3 times the helium over the winter pole, 0.5 over the summer"
    (is (< (abs (- (jb/helium-factor -90.0 23.44) 2.3)) 1e-12))
    (is (< (abs (- (jb/helium-factor 90.0 23.44) 0.5)) 1e-12))))

(deftest storms
  (let [dst [-10.0 -40.0 -90.0 -150.0 -200.0 -180.0 -150.0 -120.0 -100.0 -90.0 -85.0 -80.0]]
    (testing "the temperature rises through the main phase and falls through the recovery, never below zero"
      (let [dtc (vec (jb/storm-dtc dst 4 8 20.0))]
        (is (= (count dst) (count dtc)))
        (is (< (dtc 0) (dtc 4)))
        (is (> (dtc 5) (dtc 8)))
        (is (every? #(>= % 0.0) dtc))))
    (testing "eq. 10's slope, and -1.40 for the greatest storms"
      (is (< (abs (- (jb/main-phase-slope -100.0) (+ (* -1.5050e-5 1e4) 1.0604 -3.20))) 1e-12))
      (is (== -1.40 (jb/main-phase-slope -500.0))))
    (testing "outside storms, Jacchia 1970's ap equation, ap capped at 50"
      (is (< (abs (- (jb/ap-dtc 15.0) (+ 15.0 (* 100.0 (- 1.0 (Math/exp -1.2)))))) 1e-12))
      (is (== (jb/ap-dtc 50.0) (jb/ap-dtc 300.0))))))

(def ^:private orekit
  (for [l (str/split-lines (slurp (io/resource "data/jb2008/orekit-12.2-grid.txt")))]
    (mapv parse-double (str/split (str/trim l) #"\s+"))))

(deftest against-orekit
  (testing "220 to 400 km: 648 points, dates in 2006, three levels of solar activity, storm heating on and off"
    (let [ratios (for [[mjd dec lat lst alt f10 f10b s10 s10b m10 m10b y10 y10b dtc rho] orekit]
                   (/ (:rho (jb/atmosphere {:mjd mjd :dec dec :lat lat :lst lst :alt alt :f10 f10 :f10b f10b :s10 s10
                                            :s10b s10b :m10 m10 :m10b m10b :y10 y10 :y10b y10b :dtc dtc}))
                      rho))]
      (is (= 648 (count ratios)))
      (is (< (abs (- (/ (reduce + ratios) (count ratios)) 1.0)) 0.005) "within 0.5% on average")
      (is (every? #(< (abs (- % 1.0)) 0.07) ratios) "and 7% at worst"))))
