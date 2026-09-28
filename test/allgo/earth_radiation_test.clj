(ns allgo.earth-radiation-test
  "Earth radiation pressure against what a uniform Earth must give: a
  uniform Lambertian emitter's infrared irradiance, (e Phi/4)(R/r)^2 and
  radial; the albedo irradiance over the subsolar point by the integral
  symmetry reduces it to, tending far away to a Lambert sphere's (2/3) a
  Phi (R/r)^2; none from over the midnight point. Then Knocke's models themselves."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.earth-radiation :as er]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.quadrature :as quadrature]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private uniform
  {:albedo (constantly 0.3) :emissivity (constantly 0.7)})

(def ^:private sun [c/AU 0.0 0.0])

(def ^:private flux (* c/solar-pressure 299792458.0))

(deftest infrared
  (testing "a uniform emitter: (e Phi/4)(R/r)^2, radially out, at any height and direction"
    (doseq [r [[7000.0 0.0 0.0] [0.0 -12000.0 5000.0] [-20000.0 3000.0 -1000.0]]]
      (let [{:keys [infrared]} (er/irradiance r sun 50000.0 {:models uniform :rings 200 :sectors 72})
            rm (v3/length r)
            expected (* 0.7 0.25 flux (math/pow (/ c/R-earth rm) 2))]
        (is (< (abs (- (v3/length infrared) expected)) (* 1e-3 expected)) (str r))
        (is (< (v3/length (v3/cross (v3/normalize infrared) (v3/normalize r))) 1e-9))))))

(defn- subsolar-albedo
  "The albedo irradiance over the subsolar point at distance `r`, by the
  one-dimensional integral symmetry leaves:
  2 a Phi R^2 integral of cos t (r cos t - R)(r - R cos t)/d^4 sin t dt,
  d^2 = R^2 + r^2 - 2 R r cos t, out to the horizon."
  [a r]
  (let [R c/R-earth
        th (math/acos (/ R r))
        f (fn [t] (let [ct (math/cos t) d2 (+ (* R R) (* r r) (* -2.0 R r ct))]
                    (/ (* ct (- (* r ct) R) (- r (* R ct)) (math/sin t)) (* d2 d2))))]
    (* 2.0 a flux R R (quadrature/simpson f 0.0 th 20000))))

(deftest albedo
  (testing "over the subsolar point, the integral symmetry leaves"
    (doseq [h [1000.0 10000.0 36000.0]]
      (let [r [(+ c/R-earth h) 0.0 0.0]
            {:keys [albedo]} (er/irradiance r sun 50000.0 {:models uniform :rings 200 :sectors 72})
            expected (subsolar-albedo 0.3 (+ c/R-earth h))]
        (is (< (abs (- (v3/length albedo) expected)) (* 2e-3 expected)) (str h)))))
  (testing "and far away, a Lambert sphere at full phase: (2/3) a Phi (R/r)^2"
    (let [rm 2e6
          {:keys [albedo]} (er/irradiance [rm 0.0 0.0] sun 50000.0 {:models uniform :rings 200 :sectors 72})
          expected (* (/ 2.0 3.0) 0.3 flux (math/pow (/ c/R-earth rm) 2))]
      ;; the limit's own correction is of order R/r, some 0.3% here
      (is (< (abs (- (v3/length albedo) expected)) (* 5e-3 expected)))))
  (testing "nothing reflected to a satellite over the midnight point that sees no daylight"
    (let [{:keys [albedo]} (er/irradiance [-8000.0 0.0 0.0] sun 50000.0 {:models uniform})]
      (is (zero? (v3/length albedo))))))

(deftest knocke-models
  (testing "Knocke's zonal models: at the equator 0.34 - 0.29/2 and 0.68 + 0.18/2, the poles swinging with the season"
    (let [{:keys [albedo emissivity]} er/knocke
          dec22 44960.0 jun22 (+ 44960.0 (/ 365.25 2))]
      (is (< (abs (- (albedo 0.0 dec22) (- 0.34 0.145))) 1e-12))
      (is (< (abs (- (emissivity 0.0 dec22) (+ 0.68 0.09))) 1e-12))
      ;; the southern summer pole, in December, the darker and warmer
      (is (> (albedo (/ math/PI 2) dec22) (albedo (/ math/PI -2) dec22)))
      (is (< (abs (- (albedo (/ math/PI 2) dec22) (albedo (/ math/PI -2) jun22))) 1e-12)))))

(deftest acceleration
  (testing "the acceleration is the irradiance over c, and of the order of 1e-9 km/s^2 low down"
    (let [r [7000.0 0.0 0.0]
          {:keys [albedo infrared]} (er/irradiance r sun 50000.0)
          a (er/acceleration r sun 50000.0 0.02 1.3)]
      (is (< (v3/distance a (v3/scale (v3/add albedo infrared) (/ (* 1.3 0.02 1e-3) 299792458.0))) 1e-22))
      (is (< 1e-11 (v3/length a) 1e-9)))))
