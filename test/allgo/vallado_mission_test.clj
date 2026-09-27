(ns allgo.vallado-mission-test
  "Vallado chapter 11's mission geometry, each checked against what it
  must agree with: great circles against their own inverse and the
  ellipsoid's distance, the ground point against the place it came from,
  the frozen eccentricity against J2 and J3's acceleration averaged round
  the orbit, and the field of view against a ray traced to the sphere."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.mission :as ms]
            [allgo.astro.perturbations :as pt]
            [allgo.astro.reduction :as rd]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(def ^:private R c/R-earth)
(defn- deg [x] (math/to-radians x))
(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))

(deftest between-places
  (testing "along the equator, a quarter of the way round, due east"
    (let [{:keys [range azimuth]} (ms/range-azimuth [0.0 0.0] [0.0 (deg 90)])]
      (is (close? range (* R 0.5 math/PI) 1e-12))
      (is (close? azimuth (deg 90) 1e-12))))
  (testing "toward the pole is due north"
    (is (< (:azimuth (ms/range-azimuth [(deg 10) (deg 20)] [(deg 50) (deg 20)])) 1e-12)))
  (testing "setting out along the azimuth for the range arrives"
    (doseq [[p1 p2] [[[(deg 39.7) (deg -105.0)] [(deg 51.5) (deg -0.1)]]
                     [[(deg -33.9) (deg 151.2)] [(deg 35.7) (deg 139.7)]]
                     [[(deg 1.0) (deg 2.0)] [(deg -60.0) (deg -170.0)]]]]
      (let [{:keys [range azimuth]} (ms/range-azimuth p1 p2)
            [lat lon] (ms/destination p1 range azimuth)]
        (is (close? lat (first p2) 1e-12))
        (is (< (abs (am/wrap-angle (- lon (second p2)))) 1e-12)))))
  (testing "and the sphere's distance is the ellipsoid's to its flattening"
    (let [p1 [(deg 39.7) (deg -105.0)] p2 [(deg 51.5) (deg -0.1)]]
      (is (close? (:range (ms/range-azimuth p1 p2)) (gd/distance p1 p2) 4e-3)))))

(deftest ground-track
  (testing "the point beneath a satellite is the place it was put over"
    (let [tt 60580.3 ut (- tt (/ 69.0 86400.0))
          [lat lon h] [(deg 40.0) (deg -100.0) 500.0]
          ecef (gd/geodetic->cartesian lat lon h)
          [r] (rd/ecef->eci [ecef [0.0 0.0 0.0]] tt ut {})
          [lat' lon' h'] (ms/ground-point [r [0.0 0.0 0.0]] tt ut)]
      (is (close? lat' lat 1e-10))
      (is (< (abs (am/wrap-angle (- lon' lon))) 1e-10))
      (is (close? h' h 1e-8)))))

(deftest repeat-ground-tracks
  (testing "fifteen revolutions a day, sun-synchronous: the node comes back over the same place"
    (let [a (ms/repeat-ground-track 15 1.0 0.0 (deg 98))
          {:keys [days]} (ms/repeat-period a 0.0 (deg 98) 15)]
      (is (close? days 1.0 1e-12))
      (is (< 450.0 (- a R) 700.0) "a low orbit")))
  (testing "slower than the node's two-body period by what J2 adds"
    (let [a (ms/repeat-ground-track 14 1.0 0.0 (deg 51.6))]
      (is (not (close? (ms/nodal-period a 0.0 (deg 51.6)) (* 2 math/PI (math/sqrt (/ (* a a a) mu))) 1e-5))))))

(defn- zonal-accel
  "J2 and J3's acceleration beyond the point mass."
  [r _]
  (v3/add (geo/acceleration {:GM mu :R R :normalized? false
                             :C {[0 0] 1.0 [2 0] (- geo/J2) [3 0] (- ms/J3)} :S {}} r 3)
          (v3/scale r (/ mu (math/pow (v3/length r) 3)))))

(deftest frozen-orbits
  (let [a 7078.0 i (deg 98)
        ef (ms/frozen-eccentricity a i)
        argp-rate (fn [e] (:argp (pt/averaged-rates mu {:a a :e e :i i :raan 0.3 :argp (/ math/PI 2) :M 0.0}
                                                    zonal-accel 720)))
        j2-alone (:argp (pt/j2-secular a ef i))]
    (testing "about a thousandth, for a low orbit"
      (is (< 5e-4 ef 2e-3)))
    (testing "J2 and J3 together leave the perigee still there"
      (is (< (abs (argp-rate ef)) (* 0.02 (abs j2-alone)))))
    (testing "and turn it opposite ways either side"
      (is (neg? (* (argp-rate (* 0.5 ef)) (argp-rate (* 2.0 ef))))))
    (testing "with the eccentricity itself steady"
      (is (< (abs (:e (pt/averaged-rates mu {:a a :e ef :i i :raan 0.3 :argp (/ math/PI 2) :M 0.0} zonal-accel 720)))
             1e-14)))))

(deftest field-of-view
  (testing "a ray traced from the satellite to the sphere sees what the formulas say"
    (let [h 700.0 S [(+ R h) 0.0 0.0]]
      (doseq [eta [(deg 5) (deg 30) (deg 55)]]
        (let [{:keys [elevation central-angle slant-range swath]} (ms/field-of-view h eta)
              d [(- (math/cos eta)) (math/sin eta) 0.0]
              ;; the near intersection of S + t d with the sphere
              b (v3/dot S d) cc (- (v3/dot S S) (* R R))
              t (- (- b) (math/sqrt (- (* b b) cc)))
              P (v3/add-scaled S d t)
              up (v3/normalize P)
              to-sat (v3/sub S P)]
          (is (close? slant-range t 1e-10))
          (is (close? central-angle (math/atan2 (second P) (first P)) 1e-10))
          (is (close? elevation (math/asin (/ (v3/dot to-sat up) (v3/length to-sat))) 1e-10))
          (is (close? swath (* 2 R central-angle) 1e-12))))))
  (testing "at the horizon the elevation is zero, and past it nothing"
    (let [h 700.0 {:keys [nadir central-angle]} (ms/horizon h)
          fov (ms/field-of-view h (- nadir 1e-12))]
      (is (< (:elevation fov) 1e-5))
      (is (close? (:central-angle fov) central-angle 1e-5))
      (is (nil? (ms/field-of-view h (+ nadir 1e-6)))))))
