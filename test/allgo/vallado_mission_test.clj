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
            [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
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

;; ------------------------------------------- repeat ground-track upkeep

(defn- j2-accel
  "The point mass and J2, with an optional steady drag-like deceleration
  `f` km/s^2 along the velocity."
  [f]
  (fn [_ [x y z vx vy vz]]
    (let [r [x y z] v [vx vy vz]
          a (v3/add (geo/acceleration {:GM mu :R R :normalized? false :C {[0 0] 1.0 [2 0] (- geo/J2)} :S {}} r 2)
                    (v3/scale (v3/normalize v) (- f)))]
      (into v a))))

(defn- nodes
  "The ascending nodes over `duration` seconds of the orbit from `s0`:
  `[[t longitude] ...]`, the longitude east of Greenwich (the Earth turned
  from 0 at t = 0), each crossing refined by Newton's method on z."
  [f s0 duration]
  (let [rhs (j2-accel f)
        opts {:tol-abs 1e-12 :tol-rel 1e-12}
        fly (fn [y dt] (:y (core/step-until (rk/integrator rk/dopri54 rhs 0.0 y (* 0.1 dt) opts) dt)))
        steps (take-while #(< (:t %) duration) (core/trajectory (rk/integrator rk/dopri54 rhs 0.0 s0 30.0 opts)))]
    (vec (for [[a b] (partition 2 1 steps)
               :when (and (neg? (get-in a [:y 2])) (not (neg? (get-in b [:y 2]))))]
           (let [{:keys [t y]} (loop [t (:t a) y (:y a) i 0]
                                 (let [dt (- (/ (y 2) (y 5)))]
                                   (if (or (< (abs (y 2)) 1e-9) (> i 10))
                                     {:t t :y y}
                                     (recur (+ t dt) (if (pos? dt) (fly y dt) y) (inc i)))))]
             [t (am/wrap-angle (- (math/atan2 (y 1) (y 0)) (* c/omega-earth t)))])))))

(defn- circular [a i] [a 0.0 0.0 0.0 (* (math/sqrt (/ mu a)) (math/cos i)) (* (math/sqrt (/ mu a)) (math/sin i))])

(deftest ground-track-drift
  (testing "two orbits a km apart, flown under J2, their nodes compared rev
            for rev, slide apart at K per km"
    (let [a 7078.0 i (deg 98.2) days 2.0
          lo (nodes 0.0 (circular a i) (* days 86400.0))
          hi (nodes 0.0 (circular (+ a 1.0) i) (* days 86400.0))
          n (dec (min (count lo) (count hi)))
          apart (fn [k] (am/wrap-angle (- (second (hi k)) (second (lo k)))))
          measured (/ (- (apart n) (apart 0)) (- (first (lo n)) (first (lo 0))))]
      (is (< (abs (- measured (ms/ground-track-drift a 0.0 i))) (* 0.02 (abs measured)))))))

(deftest ground-track-maintenance
  (testing "started da0 high at the east edge, drag brings the track to the west
            edge and back in one cycle, as flown under J2 against the orbit
            without drag"
    (let [a 7078.0 i (deg 98.2) tol 1.0
          a-rate (/ -2.0 86400.0)                     ; a steep 2 km a day, to keep the flight short
          {:keys [da cycle]} (ms/ground-track-maintenance a 0.0 i a-rate tol)
          ;; the tangential deceleration that lowers a circular orbit at a-rate
          f (/ (* (- a-rate) mu) (* 2.0 a a (math/sqrt (/ mu a))))
          ref (nodes 0.0 (circular a i) (* 1.05 cycle))
          kept (nodes f (circular (+ a da) i) (* 1.05 cycle))
          ;; both start on a node, so rev for rev from there
          err (map (fn [[_ l] [_ m]] (* R (am/wrap-angle (- m l)))) ref kept)]
      (is (< 0.5 da 1.5) "about a kilometer")
      (is (< (abs (+ (apply min err) (* 2.0 tol))) (* 0.1 tol)) "the far edge just reached")
      (is (< (abs (last err)) (* 0.25 tol)) "and back by the cycle's end"))))
