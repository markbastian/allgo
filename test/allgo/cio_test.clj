(ns allgo.cio-test
  "The IAU 2006/2000A CIO-based reduction against ERFA (the BSD-licensed
  release of the IAU's SOFA library), run at 2020 January 1 7h12m TT:
  eraFal03 and its kin for the fundamental arguments, eraXy06, eraS06,
  eraEra00 and eraSp00, eraC2ixys for the CIP matrix, and eraPom00 with
  eraC2tcio for the whole GCRS-to-ITRS matrix. Then the chain inverted,
  its velocity against finite differences, and against the FK5 chain it
  succeeds."
  (:require [allgo.astro.cio :as cio]
            [allgo.astro.reduction :as rd]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private as (/ math/PI 648000.0))
(def ^:private tt (+ 58849.5 0.3))
(def ^:private ut (- tt (/ 69.184 86400.0)))
(def ^:private polar {:xp (* 0.0714 as) :yp (* 0.2997 as)})

(defn- near? [a b tol] (every? #(< (abs %) tol) (map - (flatten a) (flatten b))))

(deftest against-erfa
  (testing "the fundamental arguments"
    (is (near? (cio/fundamental-arguments (/ (- tt 51544.5) 36525.0))
               [3.1181488125460959 6.2419055013253484 4.4967403835258892 1.307768904704063 -4.5692458099077777
                4.6777186134198061 0.10513121985531626 1.7564418074767758 3.9034348975331667 4.9106155070861206
                5.1401740505367242 0.69380195841796422 6.0745783205284916 0.0048767657548789197]
               1e-14)))
  (testing "the pole's X and Y, to a picoradian"
    (is (near? (cio/xy tt) [0.0019113128678595976 -1.2381535734355221e-05] 1e-15)))
  (testing "s, the Earth rotation angle and s'"
    (let [[x y] (cio/xy tt)]
      (is (< (abs (- (cio/s tt x y) 6.7558895973451562e-10)) 1e-18)))
    (is (< (abs (- (cio/earth-rotation-angle ut) 0.4950628283283578)) 1e-15))
    (is (< (abs (- (cio/s-prime tt) -4.5574357584297403e-11)) 1e-20)))
  (testing "the CIP matrix and the whole GCRS to ITRS matrix"
    (is (near? (cio/c2i-matrix tt)
               [[0.99999817343989239 1.1156916133146788e-08 -0.0019113128678679623]
                [1.2508092820176508e-08 0.99999999992334876 1.2381534443093346e-05]
                [0.0019113128678595976 -1.238153573435522e-05 0.99999817336324104]]
               1e-15))
    (is (near? (cio/c2t-matrix tt ut polar)
               [[0.87993726201928502 0.47508694703296195 -0.001675610009557512]
                [-0.47508606125125197 0.8799388573258623 0.00091748178591994961]
                [0.0019103179777648727 -1.1267451021113248e-05 0.99999817527746926]]
               1e-15))))

(def ^:private gcrs [[5102.5096 6123.01152 6378.1363] [-4.74322 0.790536 5.533756]])

(deftest the-chain
  (let [eop (merge polar {:dx (* 1e-4 as) :dy (* -2e-4 as) :lod 0.0015})]
    (testing "and back"
      (let [[r v] (cio/itrs->gcrs (cio/gcrs->itrs gcrs tt ut eop) tt ut eop)]
        (is (< (v3/distance r (first gcrs)) 1e-9))
        (is (< (v3/distance v (second gcrs)) 1e-12))))
    (testing "the Earth-fixed velocity is the rate of the Earth-fixed position"
      (let [h 0.5
            at (fn [dt] (let [[r v] gcrs d (/ dt 86400.0)]
                          (first (cio/gcrs->itrs [(v3/add-scaled r v dt) v] (+ tt d) (+ ut d) eop))))
            fd (v3/scale (v3/sub (at h) (at (- h))) (/ 1.0 (* 2 h)))]
        (is (< (v3/distance (second (cio/gcrs->itrs gcrs tt ut eop)) fd) 1e-6))))
    (testing "the pole offsets move the pole by their size"
      (let [p0 (first (cio/gcrs->itrs gcrs tt ut polar))
            p1 (first (cio/gcrs->itrs gcrs tt ut (assoc polar :dx (* 1e-3 as))))]
        (is (< 0.5e-3 (/ (v3/distance p0 p1) (* as (v3/length p0))) 1.5e-3))))))

(deftest against-fk5
  ;; the IAU 1976/1980 and 2006/2000A models part by tens of milliarcseconds
  ;; without the IERS corrections: 1.9 m here
  (testing "the FK5 chain and this one put a point within a few meters of each other"
    (is (< (v3/distance (first (cio/gcrs->itrs gcrs tt ut polar))
                        (first (rd/eci->ecef gcrs tt ut polar)))
           5e-3))))
