(ns allgo.jacchia-roberts-test
  "Roberts's analytic model against what his paper claims of it, and its
  closed forms against the integrals they replace: identical to
  Jacchia's 1970 models, integrated numerically, from 90 to 125 km; within
  the percentages of his Table II above; the numerical 1970 models
  themselves against rows of Jacchia's SAO Special Report 313."
  (:require [allgo.astro.jacchia :as j]
            [allgo.astro.jacchia-roberts :as r]
            [allgo.numerics.quadrature :as quadrature]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest jacchia-1970
  (testing "the numerical 1970 models against SR 313's 900 K table"
    ;; [z T log n(N2) log n(O2) rho]
    (doseq [[z t n2 o2 rho] [[100.0 195.2 12.9492 12.3138 5.476e-10]
                             [125.0 408.4 11.3263 10.5151 1.332e-11]
                             [200.0 814.6 9.4760 8.4442 2.619e-13]
                             [400.0 895.6 6.4414 4.9837 2.072e-15]]]
      (let [{m-t :t :keys [n] m-rho :rho} (j/static j/j70 900.0 z)]
        (is (< (abs (- m-t t)) 0.06) (str z))
        (is (< (abs (- (math/log10 (:N2 n)) n2)) 3e-4) (str z))
        (is (< (abs (- (math/log10 (:O2 n)) o2)) 3e-4) (str z))
        (is (< (abs (- m-rho rho)) (* 1e-3 rho)) (str z)))))
  (testing "and its 600 K table, where the hydrogen counted below 500 km lowers the mean molecular mass"
    (doseq [[z rho m] [[300.0 3.684e-15 15.84] [350.0 8.505e-16 14.18] [400.0 2.177e-16 11.27]]]
      (let [s (j/static j/j70 600.0 z)]
        (is (< (abs (- (:rho s) rho)) (* 1e-3 rho)) (str z))
        (is (< (abs (- (:m s) m)) 0.015) (str z))))))

(deftest identical-below-125
  (testing "from 90 to 125 km Roberts's closed forms are Jacchia's integrals -- but for the trace of
            hydrogen Jacchia counts there and Roberts leaves out"
    (doseq [tinf [600.0 1000.0 1600.0 2000.0]
            z [92.0 100.0 105.0 110.0 118.0 125.0]]
      (let [a (r/static tinf z) b (j/static j/j70 tinf z)
            rho-b (j/mass-density (dissoc (:n b) :H))]
        (is (< (abs (- (:t a) (:t b))) 1e-9) (str tinf " " z))
        (is (< (abs (- (:rho a) rho-b)) (* 1e-7 rho-b)) (str tinf " " z))
        (doseq [k [:N2 :O2 :O :Ar :He]]
          (is (< (abs (- (get-in a [:n k]) (get-in b [:n k]))) (* 1e-7 (get-in b [:n k]))) (str tinf " " z " " k)))))))

(defn- worst-difference
  "The largest difference, percent, between Roberts's densities and
  Jacchia's at `tinf`, over `heights`."
  [tinf heights]
  (apply max (for [z heights] (* 100.0 (abs (- (/ (:rho (r/static tinf z)) (:rho (j/static j/j70 tinf z))) 1.0))))))

(deftest close-above-125
  (testing "above 125 km, the largest difference from 125 to 1000 km that Roberts's Table II gives"
    (doseq [[tinf pct] [[1300.0 1.08] [2000.0 4.86]]]
      (let [worst (worst-difference tinf (range 130.0 1001.0 10.0))]
        (is (< (abs (- worst pct)) 0.05) (str tinf " K: " worst "%")))))
  (testing "and at 600 K the same, but for just below 500 km, where Roberts leaves out the hydrogen
            Jacchia's tables count -- a few percent of the density at so low a temperature"
    (let [worst (worst-difference 600.0 (remove #(< 430.0 % 501.0) (range 130.0 1001.0 10.0)))]
      (is (< (abs (- worst 3.89)) 0.05) (str worst "%")))))

(deftest closed-forms
  (testing "the temperature profile's roots: the quartic and its factored form agree, and the roots lie where Roberts says"
    (doseq [tinf [600.0 1300.0 2000.0]]
      (let [{:keys [r1 r2 x y]} (#'r/profile tinf)]
        (is (< 164.0 r1 168.0)) (is (< 57.0 r2 66.0)) (is (< 97.0 x 101.0)) (is (< 22.0 y 33.0))
        (doseq [z [90.0 100.0 112.0 125.0]]
          (is (< (abs (- (r/temperature tinf z) (j/temperature j/j70 tinf z))) 1e-9))))))
  (testing "above 125 km each density is the diffusion equation integrated through Roberts's profile"
    (let [tinf 1100.0
          temp #(r/temperature tinf %)
          n125 (r/static tinf 125.0)]
      (doseq [z [200.0 450.0 900.0]
              k [:N2 :He]]
        (let [{:keys [m alpha]} (j/species k)
              ;; ln n(z)/n(125) = -int M g/(R T) dz - (1 + alpha) ln T(z)/T(125), per km
              path (quadrature/simpson #(/ (* m 9.80665 (math/pow (/ 6356.766 (+ 6356.766 %)) 2)) (* 8.31432 (temp %)))
                                       125.0 z 4000)
              expected (* (get-in n125 [:n k]) (math/exp (- (+ path (* (+ 1.0 alpha) (math/log (/ (temp z) (temp 125.0))))))))]
          (is (< (abs (- (get-in (r/static tinf z) [:n k]) expected)) (* 1e-8 expected)) (str z " " k))))))
  (testing "the length l interpolates Table II, and gives the continuous gradient near 1.9 (Ra + Zx)"
    (is (= 11825.0 (r/length-l 600.0)))
    (is (= 14515.0 (r/length-l 2000.0)))
    (is (< 12000.0 (r/length-l 900.0) 13515.0))))

(deftest atmosphere
  (testing "the J71 variations on Roberts's densities: within a few percent of them on J70's numerical ones"
    (let [in {:mjd 41000.0 :alt 400.0 :lat 30.0 :lst 14.0 :f107 150.0 :f107a 150.0 :kp 2.0}
          a (r/atmosphere in)
          b (j/atmosphere in {:static-model #(j/static j/j70 %1 %2)})]
      (is (== (:t-exo a) (:t-exo b)))
      (is (< (abs (- (/ (:rho a) (:rho b)) 1.0)) 0.05)))))
