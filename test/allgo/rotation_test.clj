(ns allgo.rotation-test
  (:require [allgo.astro.ephemeris :as eph]
            [allgo.astro.rotation :as rot]
            [allgo.astro.time :as time]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private j2000 51544.5)

(defn- close? [a b tol] (< (abs (- a b)) tol))

(deftest frame-test
  (testing "every model gives a rotation: orthonormal, right-handed"
    (doseq [body (keys rot/models) mjd [j2000 (+ j2000 7300.3)]]
      (let [m (rot/body->icrf body mjd)]
        (is (close? 1.0 (v3/dot (lin/mat-vec m [0 0 1]) (lin/mat-vec m [0 0 1])) 1e-12) (name body))
        (is (< (reduce max (map abs (flatten (lin/mat-sub (lin/mat-mul m (lin/transpose m)) (lin/eye 3))))) 1e-12)
            (name body))
        (is (close? 1.0 (v3/dot (v3/cross (lin/mat-vec m [1 0 0]) (lin/mat-vec m [0 1 0]))
                                (lin/mat-vec m [0 0 1]))
                    1e-12)
            (name body)))))

  (testing "the body's z axis is its pole"
    (doseq [body (keys rot/models)]
      (is (< (v3/distance (rot/pole body j2000) (lin/mat-vec (rot/body->icrf body j2000) [0 0 1]))
             1e-12)
          (name body)))))

(deftest earth-test
  (testing "the pole is the ICRF's at J2000"
    (is (< (v3/distance [0.0 0.0 1.0] (rot/pole :earth j2000)) 1e-12)))

  (testing "Greenwich faces where sidereal time says it does, to the
            fraction of a degree the IAU's simple model manages: it starts
            0.31 degrees off at J2000 and drifts, through its slightly slow
            spin and its linear pole, to about 0.65 by 2022"
    (doseq [mjd [j2000 (+ j2000 0.3) (+ j2000 3652.7) (+ j2000 8000.1)]]
      (let [[x y] (lin/mat-vec (rot/body->icrf :earth mjd) [1 0 0])]
        (is (< (abs (am/wrap-angle (- (math/atan2 y x) (time/gmst mjd))))
               (math/to-radians 0.75)))))))

(deftest moon-test
  (testing "the Moon keeps one face to the Earth: over a month the
            sub-Earth point wanders only as far as the librations take it"
    (let [points (for [k (range 60)]
                   (let [mjd (+ j2000 (* 0.5 k))
                         to-earth (v3/normalize (v3/negate (eph/moon mjd)))
                         [x y z] (lin/mat-vec (lin/transpose (rot/body->icrf :moon mjd)) to-earth)]
                     [(math/to-degrees (math/atan2 y x)) (math/to-degrees (math/asin z))]))]
      (is (every? (fn [[lon lat]] (and (< (abs lon) 9.0) (< (abs lat) 8.0))) points))
      (is (< (abs (/ (reduce + (map first points)) (count points))) 2.0)))))

(deftest spin-test
  (testing "day lengths, and which way round"
    (is (close? 23.934 (rot/sidereal-day :earth) 0.001))
    (is (close? 24.623 (rot/sidereal-day :mars) 0.001))
    (is (neg? (rot/sidereal-day :venus)))
    (is (neg? (rot/sidereal-day :uranus)))
    (is (pos? (rot/sidereal-day :neptune)))))
