(ns allgo.vallado-frames-test
  "The FK5 reduction of Vallado chapter 3, stage by stage, against ERFA
  (the BSD-licensed release of the IAU's SOFA library): its IAU 1976
  precession, IAU 1980 nutation with the IERS corrections, 1982 mean
  sidereal time, 1994 equation of the equinoxes and polar motion, run on
  the state and dates below. The velocities are checked against finite
  differences of the positions, and every stage must invert."
  (:require [allgo.astro.reduction :as rd]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private as (/ math/PI 180.0 3600.0))
(def ^:private tt 53101.000754)
(def ^:private ut 53100.999931)
(def ^:private nutation {:ddpsi (* -0.052195 as) :ddeps (* -0.003875 as)})
(def ^:private polar {:xp (* -0.140682 as) :yp (* 0.333309 as)})
(def ^:private eci [[5102.5096 6123.01152 6378.1363] [-4.74322 0.790536 5.533756]])

(defn- near? [a b tol] (< (v3/distance a b) tol))

(deftest fk5-stages
  ;; ERFA: eraPmat76, eraNut80 and eraObl80 through eraNumat with the
  ;; corrections added, eraGmst82 + eraEqeq94, eraPom00
  (testing "mean and true of date, to the micrometer"
    (is (near? (first (rd/eci->mod eci tt)) [5094.0308007762314 6127.8699152122144 6380.2474447706718] 1e-9))
    (is (near? (first (rd/eci->tod eci tt nutation)) [5094.5172789594199 6127.364665291826 6380.3442623764176] 1e-9)))
  (testing "TEME: true of date turned by the equation of the equinoxes without the 1997 terms"
    (is (near? (first (rd/eci->teme eci tt nutation)) [5094.182113484967 6127.6433187427565 6380.3442623764176] 1e-9)))
  ;; to a tenth of a millimeter: the mean sidereal time here is Meeus's
  ;; form of the 1982 expression, which rounds differently by 1e-11 rad
  (testing "pseudo-Earth-fixed and Earth-fixed"
    (is (near? (first (rd/eci->pef eci tt ut {})) [-6473.9396675043063 -4646.1598353521968 6380.3448903172057] 1e-6))
    (is (near? (first (rd/eci->ecef eci tt ut polar)) [-6473.9440191888498 -4646.1701455146076 6380.3329669273753] 1e-6))))

(deftest rotating-velocities
  (testing "the Earth-fixed velocity is the rate of the Earth-fixed position"
    (let [h 0.5
          at (fn [dt] (let [[r v] eci d (/ dt 86400.0)]
                        (first (rd/eci->ecef [(v3/add-scaled r v dt) v] (+ tt d) (+ ut d) polar))))
          fd (v3/scale (v3/sub (at h) (at (- h))) (/ 1.0 (* 2 h)))]
      (is (near? (second (rd/eci->ecef eci tt ut polar)) fd 1e-6)))))

(deftest inverses
  (let [eop (merge nutation polar {:lod 0.0015563})]
    (doseq [[to from] [[#(rd/eci->ecef % tt ut eop) #(rd/ecef->eci % tt ut eop)]
                       [#(rd/eci->pef % tt ut eop) #(rd/pef->eci % tt ut eop)]
                       [#(rd/eci->teme % tt eop) #(rd/teme->eci % tt eop)]
                       [#(rd/eci->tod % tt eop) #(rd/tod->eci % tt eop)]
                       [#(rd/eci->mod % tt) #(rd/mod->eci % tt)]]]
      (let [[r v] (from (to eci))]
        (is (near? r (first eci) 1e-9))
        (is (near? v (second eci) 1e-12))))))
