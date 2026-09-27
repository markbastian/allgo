(ns allgo.vallado-states-test
  "Vallado chapters 2 to 4's state representations and local frames, each
  checked by its definition: every conversion inverts, angles point where
  the vectors do, rates match finite differences of the angles, and each
  frame is the orthonormal set its name says. Where the library already
  had an independent version -- classical elements, look angles -- the
  two must agree."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.kepler :as kep]
            [allgo.astro.reduction :as rd]
            [allgo.astro.states :as st]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))
(defn- all-close? [as bs tol] (every? true? (map #(close? %1 %2 tol) (flatten as) (flatten bs))))
(defn- same-angle? [a b tol] (< (abs (am/wrap-angle (- a b))) tol))
(defn- deg [x] (math/to-radians x))

(def ^:private rv [[1526.0 -5867.2 3499.6] [1.483 -7.093 0.9566]])

(defn- rates-by-difference
  "The rates of `f`'s first three outputs, by central differences along
  straight-line motion at the state's velocity."
  [f [r v]]
  (let [h 1e-3
        a (f [(v3/add-scaled r v h) v]) b (f [(v3/add-scaled r v (- h)) v])]
    (map (fn [x y] (/ (am/wrap-angle (- x y)) (* 2 h))) (take 3 a) (take 3 b))))

(defn- orthonormal? [[a b cc]]
  (and (every? #(close? (v3/length %) 1.0 1e-12) [a b cc])
       (every? #(< (abs %) 1e-12) [(v3/dot a b) (v3/dot b cc) (v3/dot a cc)])
       (< (v3/distance (v3/cross a b) cc) 1e-12)))

(deftest spherical-forms
  (testing "range, right ascension and declination: where the vector points, and inverts"
    (let [[rr ra dec :as out] (st/state->radec rv)
          [x y z] (first rv)]
      (is (close? rr (v3/length (first rv)) 1e-12))
      (is (same-angle? ra (math/atan2 y x) 1e-12))
      (is (close? dec (math/asin (/ z rr)) 1e-12))
      (is (all-close? (st/radec->state out) rv 1e-10))
      (is (all-close? (drop 3 out) (rates-by-difference st/state->radec rv) 1e-6))))
  (testing "ecliptic: the equatorial form turned about x by the obliquity"
    (let [[_ lon lat] (st/state->ecliptic rv)
          e (deg 23.439291)
          [x y z] (first rv)
          ;; the vector's ecliptic components
          ye (+ (* y (math/cos e)) (* z (math/sin e)))
          ze (- (* z (math/cos e)) (* y (math/sin e)))]
      (is (same-angle? lon (math/atan2 ye x) 1e-12))
      (is (close? lat (math/asin (/ ze (v3/length (first rv)))) 1e-12))
      (is (all-close? (st/ecliptic->state (st/state->ecliptic rv)) rv 1e-10))))
  (testing "topocentric: the same, from the site"
    (let [site [[-2968.7 3980.6 3992.9] [-0.29 -0.216 0.0]]
          tr (st/state->topocentric-radec rv site)]
      (is (close? (first tr) (v3/distance (first rv) (first site)) 1e-12))
      (is (all-close? (st/topocentric-radec->state tr site) rv 1e-10))))
  (testing "adbar: angle from the vertical and azimuth, and inverts"
    (let [[rm vm _ _ fpav :as out] (st/state->adbar rv)
          [r v] rv]
      (is (close? rm (v3/length r) 1e-12))
      (is (close? vm (v3/length v) 1e-12))
      (is (close? (math/cos fpav) (/ (v3/dot r v) (* rm vm)) 1e-12))
      (is (all-close? (st/adbar->state out) rv 1e-10)))))

(deftest classical-elements
  (testing "agree with the library's own elements, and invert"
    (doseq [s [[[6524.8 6862.9 6448.3] [4.9013 5.5338 -1.9763]] rv
               [[7000.0 0.0 0.0] [0.0 5.0 6.0]]]]
      (let [{:keys [p a e i raan argp nu M arglat] :as el} (st/state->classical s)
            k (kep/state->elements c/GM-earth (first s) (second s))]
        (is (close? a (:a k) 1e-10))
        (is (close? e (:e k) 1e-10))
        (is (same-angle? i (:i k) 1e-10))
        (is (same-angle? raan (:raan k) 1e-10))
        (is (same-angle? argp (:argp k) 1e-10))
        (is (same-angle? nu (:nu k) 1e-10))
        (is (same-angle? M (:M k) 1e-10))
        (is (close? p (* a (- 1.0 (* e e))) 1e-10))
        (is (same-angle? arglat (+ argp nu) 1e-10))
        (is (all-close? (st/classical->state el) s 1e-9)))))
  (testing "the orbit types, and the angles that stand in"
    (let [circular-inclined [[7000.0 0.0 0.0] [0.0 (* (math/sqrt (/ c/GM-earth 7000.0)) (math/cos 0.5))
                                               (* (math/sqrt (/ c/GM-earth 7000.0)) (math/sin 0.5))]]
          circular-equatorial [[0.0 7000.0 0.0] [(- (math/sqrt (/ c/GM-earth 7000.0))) 0.0 0.0]]
          elliptical-equatorial [[7000.0 0.0 0.0] [0.0 8.5 0.0]]]
      (is (= :elliptical-inclined (:type (st/state->classical rv))))
      (let [{:keys [type arglat]} (st/state->classical circular-inclined)]
        (is (= :circular-inclined type))
        (is (same-angle? arglat 0.0 1e-9) "at the ascending node"))
      (let [{:keys [type truelon]} (st/state->classical circular-equatorial)]
        (is (= :circular-equatorial type))
        (is (same-angle? truelon (/ math/PI 2) 1e-9)))
      (let [{:keys [type lonper]} (st/state->classical elliptical-equatorial)]
        (is (= :elliptical-equatorial type))
        (is (same-angle? lonper 0.0 1e-9) "periapsis on x"))))
  (testing "the perifocal axes: orthonormal, P at periapsis, W along the momentum"
    (let [[i raan argp] [(deg 87.87) (deg 227.89) (deg 53.38)]
          [P _ W :as axes] (st/perifocal-axes i raan argp)
          [r v] (st/classical->state {:p 11000.0 :e 0.8 :i i :raan raan :argp argp :nu 0.0})]
      (is (orthonormal? axes))
      (is (close? (v3/dot P (v3/normalize r)) 1.0 1e-12))
      (is (close? (v3/dot W (v3/normalize (v3/cross r v))) 1.0 1e-12)))))

(deftest equinoctial-elements
  (testing "invert, and agree with the classical elements"
    (doseq [s [rv [[6524.8 6862.9 6448.3] [4.9013 5.5338 -1.9763]] [[7000.0 0.0 0.0] [0.0 5.0 6.0]]]]
      (let [eq (st/state->equinoctial s)
            {:keys [a e raan argp M nu]} (st/state->classical s)]
        (is (close? (:a eq) a 1e-10))
        (is (close? (math/hypot (:af eq) (:ag eq)) e 1e-10))
        (is (same-angle? (:meanlon eq) (+ M argp (* (:fr eq) raan)) 1e-10) "fr -1 for a retrograde orbit")
        (is (same-angle? (:truelon eq) (+ nu argp (* (:fr eq) raan)) 1e-10))
        (is (close? (:n eq) (kep/mean-motion c/GM-earth a) 1e-12))
        (is (all-close? (st/equinoctial->state eq) s 1e-9)))))
  (testing "from elements to a state and back"
    (let [eq {:a 7000.0 :af 0.001 :ag 0.001 :chi 0.001 :psi 0.001 :meanlon (/ math/PI 4) :fr 1.0}
          out (st/state->equinoctial (st/equinoctial->state eq))]
      (doseq [k [:a :af :ag :chi :psi :meanlon]]
        (is (close? (k out) (k eq) 1e-9) (str k))))))

(deftest local-frames
  (testing "RSW: radial, along-track, cross-track"
    (let [axes (st/rsw rv)
          [[rx ry rz] [vx vy vz]] (st/into-frame axes rv)]
      (is (orthonormal? axes))
      (is (close? rx (v3/length (first rv)) 1e-12))
      (is (every? #(< (abs %) 1e-9) [ry rz vz]))
      (is (pos? vy) "moving forward")
      (is (close? (math/hypot vx vy) (v3/length (second rv)) 1e-12))))
  (testing "NTW: the velocity all along T"
    (let [axes (st/ntw rv)
          [[_ _ rz] [vn vt vw]] (st/into-frame axes rv)]
      (is (orthonormal? axes))
      (is (close? vt (v3/length (second rv)) 1e-12))
      (is (every? #(< (abs %) 1e-9) [vn vw rz]))))
  (testing "PQW: in the orbit's plane, periapsis along P"
    (let [[[_ _ z] [_ _ vz] :as pqw] (st/state->pqw rv)
          {:keys [nu]} (st/state->classical rv)]
      (is (every? #(< (abs %) 1e-9) [z vz]))
      (is (same-angle? (apply math/atan2 (reverse (take 2 (first pqw)))) nu 1e-10))))
  (testing "SEZ: south, east and up, the local horizon frame reordered"
    (let [lat (deg 39.007) lon (deg -104.883)
          [S E Z :as axes] (st/sez-axes lat lon)
          [e n u] (gd/east-north-up lat lon)]
      (is (orthonormal? axes))
      (is (< (v3/distance S (v3/negate n)) 1e-12))
      (is (< (v3/distance E e) 1e-12))
      (is (< (v3/distance Z u) 1e-12)))))

(deftest radar
  (let [lat (deg 39.007) lon (deg -104.883) alt 2.1
        target [[-5505.5 56.5 3821.2] [2.3 7.5 -1.0]]
        [rho az el :as razel] (st/ecef->razel target lat lon alt)
        site (gd/geodetic->cartesian lat lon alt)
        {:keys [azimuth elevation range]} (gd/look-angles site (first target))]
    (testing "the same look angles as the library's own"
      (is (close? rho range 1e-9))
      (is (same-angle? az azimuth 1e-9))
      (is (close? el elevation 1e-9)))
    (testing "rates by finite differences"
      (is (all-close? (drop 3 razel)
                      (rates-by-difference #(st/ecef->razel % lat lon alt) target) 1e-6)))
    (testing "and every form inverts"
      (is (all-close? (st/razel->ecef razel lat lon alt) target 1e-9))
      (is (all-close? (st/sez->razel (st/razel->sez razel)) razel 1e-9)))
    (testing "straight up is elevation 90, due north azimuth 0"
      (let [up (v3/add site (v3/scale (nth (gd/east-north-up lat lon) 2) 500.0))
            north (v3/add site (v3/scale (second (gd/east-north-up lat lon)) 500.0))]
        (is (close? (nth (st/ecef->razel [up [0.0 0.0 0.0]] lat lon alt) 2) (/ math/PI 2) 1e-9))
        (is (same-angle? (second (st/ecef->razel [north [0.0 0.0 0.0]] lat lon alt)) 0.0 1e-9))))))

(deftest geodetic
  (testing "latitude and height invert"
    (let [r [6524.834 6862.875 6448.296]
          [latgd lon h] (gd/cartesian->geodetic r)]
      (is (all-close? (gd/geodetic->cartesian latgd lon h) r 1e-9)))))

(deftest flight-elements
  (let [tt 53101.0 ut 53101.0 eop {}
        fl [7000.0 7.546 (/ math/PI 6) (/ math/PI 2) (/ math/PI -6) (/ math/PI 4)]
        s (st/flight->state fl tt ut eop)
        [r v] (rd/eci->ecef s tt ut eop)]
    (testing "magnitudes, and the Earth-fixed position at that latitude and longitude"
      (is (close? (v3/length r) 7000.0 1e-12))
      (is (close? (v3/length v) 7.546 1e-12))
      (is (close? (math/asin (/ (nth r 2) 7000.0)) (/ math/PI 6) 1e-12))
      (is (same-angle? (math/atan2 (second r) (first r)) (/ math/PI 2) 1e-12)))
    (testing "the flight-path angle, below the horizontal"
      (is (close? (math/asin (/ (v3/dot r v) (* 7000.0 7.546))) (/ math/PI -6) 1e-12)))
    (testing "and back"
      (is (all-close? (take 6 (st/state->flight s tt ut eop)) fl 1e-9)))))
