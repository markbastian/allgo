(ns allgo.halo-test
  "Periodic orbits about L1 and L2 and their manifolds, against JPL's
  catalog of three-body periodic orbits (test/data/jpl-periodic): every
  sampled orbit corrected from its own start to the same orbit, and the
  halo families, continued from their branch points on the Lyapunov
  families, arriving at the catalog's orbits."
  (:require [allgo.astro.cr3bp :as cr3bp]
            [allgo.astro.halo :as h]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private catalog
  (edn/read-string (slurp (io/resource "data/jpl-periodic/earth-moon.edn"))))

(def ^:private mu (:mass-ratio catalog))

(defn- start [{:keys [family] [x z vy] :start}]
  [x 0.0 (if (= family :lyapunov) 0.0 z) 0.0 vy 0.0])

(defn- jacobi [s] (cr3bp/jacobi mu [(subvec s 0 3) (subvec s 3 6)]))

(deftest collinear-points
  (testing "the catalog's L1 and L2"
    (let [ps (cr3bp/lagrange-points mu)]
      (is (< (abs (- (first (:L1 ps)) (:L1 catalog))) 1e-13))
      (is (< (abs (- (first (:L2 ps)) (:L2 catalog))) 1e-13)))))

(deftest the-catalog
  (doseq [{:keys [family point period stability] :as o} (:orbits catalog)]
    (testing (str (name point) " " (name family) " from " (:start o))
      (let [s (start o)
            fixed (h/correct mu s)]
        (is fixed)
        (when fixed
          (testing "corrected from its own start, it is the same orbit"
            (is (< (abs (- ((:state fixed) 4) (s 4))) 1e-9))
            (is (< (abs (- ((:state fixed) 0) (s 0))) 1e-9))
            (is (< (abs (- (:period fixed) period)) 1e-9))
            (is (< (abs (- (:jacobi fixed) (:jacobi o))) 1e-9)))
          (testing "of the same stability"
            (is (< (abs (- (h/stability mu fixed) stability)) (* 1e-5 stability))))
          (testing "and flown a period it closes"
            (let [end (:state (h/propagate mu (:state fixed) (:period fixed)))]
              (is (every? #(< (abs %) 1e-8) (map - end (:state fixed)))))))))))

(deftest halo-families
  (testing "continued from the branch point on the Lyapunov family"
    (doseq [{:keys [point] [x z vy] :start :as o} (:orbits catalog)
            :when (and (= :halo (:family o)) (< 0.04 z 0.06))]
      (let [{[x' _ _ _ vy'] :state :as found} (h/halo mu point z {:step 5e-3})]
        (testing (name point)
          (is (< (abs (- x' x)) 1e-9))
          (is (< (abs (- vy' vy)) 1e-9))
          (is (< (abs (- (:period found) (:period o))) 1e-9)))))))

(deftest manifolds
  (let [o (h/correct mu (start (first (filter #(= :halo (:family %)) (:orbits catalog)))))
        st (h/stability mu o)
        lam (+ st (math/sqrt (- (* st st) 1.0)))
        period (:period o)]
    (doseq [stable? [false true]
            {:keys [tau start end]} (h/manifold mu o period {:stable? stable? :n 3 :epsilon 1e-9})]
      (testing (str (if stable? "stable" "unstable") " manifold from tau " tau)
        (let [on (:state (h/propagate mu (:state o) tau))
              dist (fn [s] (math/sqrt (reduce + (map #(* % %) (map - (subvec s 0 3) (subvec on 0 3))))))]
          (testing "a period on (back, for the stable), the departure has grown by lambda"
            (is (< (abs (- (/ (dist end) (dist start)) lam)) (* 0.005 lam))))
          (testing "and the energy is the orbit's"
            (is (< (abs (- (jacobi end) (:jacobi o))) 1e-10))))))))
