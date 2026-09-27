(ns allgo.vallado-cr3bp-test
  "The circular restricted three-body problem, checked against the
  problem itself: the rotating frame's motion against the same motion
  integrated in the inertial frame with the primaries on their circles;
  the Lagrange points against the potential's gradient and the published
  Earth-Moon places; the Jacobi constant held along a trajectory; and the
  stability of L4 changing at Routh's mass ratio."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.cr3bp :as cr]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private earth-moon (cr/mass-ratio c/GM-earth c/GM-moon))

(defn- integrate [f y0 t]
  (:y (core/step-until (rk/integrator rk/dopri54 f 0.0 y0 1e-3 {:tol-abs 1e-13 :tol-rel 1e-13}) t)))

(def ^:private y0 [0.5 0.3 0.05 0.1 -0.2 0.02])

(deftest the-rotating-frame
  (testing "the motion is the inertial motion about two primaries on circles, turned with them"
    (let [mu earth-moon t 2.0
          rot (integrate (cr/derivative mu) y0 t)
          ;; the same start in the inertial frame, which coincides with the
          ;; rotating one at t = 0: the velocity gains omega x r
          [x y z vx vy vz] y0
          inertial0 [x y z (- vx y) (+ vy x) vz]
          f (fn [t [x y z vx vy vz]]
              (let [c (math/cos t) s (math/sin t)
                    p1 [(* (- mu) c) (* (- mu) s) 0.0]
                    p2 [(* (- 1.0 mu) c) (* (- 1.0 mu) s) 0.0]
                    r [x y z]
                    pull (fn [p m] (let [d (v3/sub p r)] (v3/scale d (/ m (math/pow (v3/length d) 3)))))
                    [ax ay az] (v3/add (pull p1 (- 1.0 mu)) (pull p2 mu))]
                [vx vy vz ax ay az]))
          [ix iy iz] (integrate f inertial0 t)
          c (math/cos t) s (math/sin t)
          ;; the inertial position seen from the frame turned by t
          back [(+ (* c ix) (* s iy)) (- (* c iy) (* s ix)) iz]]
      (is (< (v3/distance back (subvec rot 0 3)) 1e-9)))))

(deftest jacobi-constant
  (testing "held along a trajectory"
    (let [mu earth-moon
          C0 (cr/jacobi mu [(subvec y0 0 3) (subvec y0 3 6)])
          y (integrate (cr/derivative mu) y0 5.0)]
      (is (< (abs (- (cr/jacobi mu [(subvec y 0 3) (subvec y 3 6)]) C0)) 1e-10)))))

(deftest lagrange-points
  (let [mu earth-moon
        {:keys [L1 L2 L3 L4] :as pts} (cr/lagrange-points mu)]
    (testing "each is a stationary point of the potential"
      (doseq [[k p] pts]
        (is (< (v3/length (cr/gradient mu p)) 1e-12) (str k))))
    (testing "the Earth-Moon points where they are usually quoted"
      (is (< (abs (- (first L1) 0.8369)) 1e-3))
      (is (< (abs (- (first L2) 1.1557)) 1e-3))
      (is (< (abs (- (first L3) -1.0051)) 1e-3))
      (is (< (abs (- (cr/jacobi mu [L1 [0.0 0.0 0.0]]) 3.1883)) 1e-3)))
    (testing "the collinear points are unstable, L4 stable for the Earth and Moon"
      (doseq [p [L1 L2 L3]] (is (not (cr/stable? mu p))))
      (is (cr/stable? mu L4)))
    (testing "and L4 turns unstable at Routh's mass ratio"
      (let [below (- cr/routh-critical 1e-4) above (+ cr/routh-critical 1e-4)]
        (is (cr/stable? below (:L4 (cr/lagrange-points below))))
        (is (not (cr/stable? above (:L4 (cr/lagrange-points above)))))))
    (testing "the zero-velocity fence: with L1's constant, the gap at L1 just closes"
      (let [C (cr/jacobi mu [L1 [0.0 0.0 0.0]])]
        (is (cr/forbidden? mu (+ C 1e-6) L1))
        (is (not (cr/forbidden? mu (- C 1e-6) L1)))))))
