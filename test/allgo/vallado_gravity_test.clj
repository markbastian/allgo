(ns allgo.vallado-gravity-test
  "Vallado chapter 8's gravity-field formulations -- the spherical partials
  and Pines's -- against the library's Montenbruck & Gill recursion and
  the numerical gradient of its potential, on a field of random
  coefficients to degree 20: three routes to one acceleration."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.gravity :as grav]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private degree 20)

(defn- random-field
  "Normalized coefficients of the size a real field's have -- 1e-5/n^2 --
  but drawn at random (from a fixed seed), so every term is exercised."
  [degree]
  (let [rng (java.util.Random. 20260927)
        draw (fn [n] (* (/ 1e-5 (* n n)) (.nextGaussian rng)))]
    {:GM c/GM-earth :R c/R-earth :normalized? true
     :C (into {[0 0] 1.0} (for [n (range 2 (inc degree)) m (range (inc n))] [[n m] (draw n)]))
     :S (into {} (for [n (range 2 (inc degree)) m (range 1 (inc n))] [[n m] (draw n)]))}))

(def ^:private field (random-field degree))

(def ^:private points
  [[6778.0 0.0 0.0] [-4000.0 5200.0 3100.0] [1200.0 -900.0 7100.0] [7000.0 7000.0 -2000.0]])

(defn- perturbing
  "The part of the acceleration beyond the point mass, which the
  comparisons are made on, since the central term would swamp them."
  [a r]
  (v3/add a (v3/scale r (/ c/GM-earth (math/pow (v3/length r) 3)))))

(deftest three-routes
  (doseq [r points
          :let [mg (geo/acceleration field r degree)
                small (* 1e-10 (v3/length (perturbing mg r)))]]
    (testing (str "at " r)
      (is (< (v3/distance (grav/pines-acceleration field r degree) mg) small) "Pines")
      (is (< (v3/distance (grav/spherical-acceleration field r degree) mg) small) "spherical partials")
      (is (< (v3/distance (grav/lear-acceleration field r degree) mg) small) "Lear, normalized")
      (is (< (v3/distance (grav/gottlieb-acceleration field r degree) mg) small) "Gottlieb, normalized"))))

(deftest gradient-of-the-potential
  (testing "Pines's acceleration is the gradient of the potential, numerically"
    (doseq [r points]
      ;; the field without its central term: the whole potential's
      ;; rounding would swamp the differences
      (let [h 1e-2
            f (assoc-in field [:C [0 0]] 0.0)
            U #(geo/potential f % degree)
            grad (mapv (fn [k] (/ (- (U (update r k + h)) (U (update r k - h))) (* 2 h))) (range 3))
            a (grav/pines-acceleration f r degree)]
        (is (< (v3/distance grad a) (* 1e-6 (v3/length a))))))))

(deftest the-poles
  (testing "Pines's, Lear's and Gottlieb's are finite on the polar axis, and agree with the recursion there"
    (doseq [r [[0.0 0.0 7000.0] [0.0 0.0 -7000.0]]
            f [grav/pines-acceleration grav/lear-acceleration grav/gottlieb-acceleration]]
      (let [a (f field r degree)
            mg (geo/acceleration field r degree)]
        (is (every? #(not (NaN? %)) a))
        (is (< (v3/distance a mg) (* 1e-10 (v3/length (perturbing mg r)))))))))

(deftest high-degree
  (testing "normalized, Lear's and Gottlieb's agree with each other to degree 60"
    (let [f60 (random-field 60)]
      (doseq [r points]
        (let [l (grav/lear-acceleration f60 r 60)
              g (grav/gottlieb-acceleration f60 r 60)]
          (is (< (v3/distance l g) (* 1e-10 (v3/length (perturbing g r))))))))))
