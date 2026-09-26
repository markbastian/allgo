(ns allgo.quaternion-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [clojure.test :refer [deftest is testing]]))

(defn- close?
  ([a b] (close? a b 1e-9))
  ([a b eps] (< (v/distance a b) eps)))

(def ^:private half-pi (/ Math/PI 2))

(deftest rotation-test
  (testing "a quarter turn about z takes x to y"
    (is (close? [0.0 1.0 0.0] (q/rotate (q/from-axis-angle [0 0 1] half-pi) [1 0 0]))))

  (testing "and about x takes y to z"
    (is (close? [0.0 0.0 1.0] (q/rotate (q/from-axis-angle [1 0 0] half-pi) [0 1 0]))))

  (testing "no rotation leaves a vector alone"
    (is (close? [1.0 -2.0 3.0] (q/rotate q/identity-q [1 -2 3]))))

  (testing "a rotation about an axis leaves that axis alone"
    (let [axis (v/normalize [1 2 -0.5])]
      (is (close? axis (q/rotate (q/from-axis-angle axis 1.234) axis)))))

  (testing "rotation preserves lengths and angles"
    ;; What makes it a rotation rather than some other transform.
    (let [r (q/from-euler 0.3 -1.1 2.2)
          a [1.0 2.0 3.0] b [-2.0 0.5 1.0]]
      (is (< (abs (- (v/length a) (v/length (q/rotate r a)))) 1e-12))
      (is (< (abs (- (v/dot a b) (v/dot (q/rotate r a) (q/rotate r b)))) 1e-12)))))

(deftest algebra-test
  (testing "multiplication reads right to left, as with matrices"
    (let [a (q/from-axis-angle [1 0 0] 0.7)
          b (q/from-axis-angle [0 1 0] 1.1)
          p [0.3 -0.5 0.9]]
      (is (close? (q/rotate (q/mul a b) p) (q/rotate a (q/rotate b p))))))

  (testing "the inverse undoes"
    (let [a (q/normalize [0.3 -0.2 0.5 0.78])
          p [1.0 2.0 -3.0]]
      (is (close? p (q/rotate (q/inverse a) (q/rotate a p))))
      (is (close? [0.0 0.0 0.0] (v/sub (q/rotate (q/mul a (q/inverse a)) p) p)))))

  (testing "for a unit quaternion the inverse is the conjugate"
    (let [a (q/from-axis-angle [1 2 3] 0.8)]
      (is (< (v/distance (vec (take 3 (q/inverse a))) (vec (take 3 (q/conjugate a)))) 1e-12))))

  (testing "normalize puts it back on the unit sphere"
    (is (< (abs (- 1.0 (q/length (q/normalize [4 3 2 1])))) 1e-12)))

  (testing "degenerate input gives the identity rather than NaN"
    (is (= q/identity-q (q/normalize [0 0 0 0])))
    (is (= q/identity-q (q/inverse [0 0 0 0]))))

  (testing "between gives the rotation from one orientation to another"
    (let [a (q/from-axis-angle [0 1 0] 0.4)
          b (q/from-axis-angle [0 1 0] 1.3)
          p [1.0 0.0 0.0]]
      (is (close? (q/rotate b p) (q/rotate (q/mul (q/between a b) a) p))))))

(deftest euler-test
  (testing "XYZ order composes as Rx*Ry*Rz"
    ;; Both readings of \"XYZ\" are in use and they are opposites. Getting
    ;; it backward gives rotations that are wrong only when two of the
    ;; three angles are nonzero, which is exactly the kind of error a demo
    ;; does not show.
    (let [x 0.3 y -0.7 z 1.2
          composed (q/mul (q/from-axis-angle [1 0 0] x)
                          (q/mul (q/from-axis-angle [0 1 0] y)
                                 (q/from-axis-angle [0 0 1] z)))
          p [0.4 0.9 -0.2]]
      (is (close? (q/rotate (q/from-euler x y z) p) (q/rotate composed p)))
      (is (not (close? (q/rotate (q/from-euler x y z) p)
                       (q/rotate (q/mul (q/from-axis-angle [0 0 1] z)
                                        (q/mul (q/from-axis-angle [0 1 0] y)
                                               (q/from-axis-angle [1 0 0] x)))
                                 p)))
          "and not as Rz*Ry*Rx")))

  (testing "a single angle is a single axis rotation"
    (is (close? (q/rotate (q/from-euler 0 0 half-pi) [1 0 0]) [0.0 1.0 0.0])))

  (testing "all zero is the identity"
    (is (close? [1.0 2.0 3.0] (q/rotate (q/from-euler 0 0 0) [1 2 3])))))

(deftest axis-angle-test
  (testing "round trip"
    (let [axis (v/normalize [1 -2 0.5]) angle 1.234
          [axis' angle'] (q/to-axis-angle (q/from-axis-angle axis angle))]
      (is (close? axis axis'))
      (is (< (abs (- angle angle')) 1e-9))))

  (testing "the identity has no angle and needs no axis"
    (let [[_ angle] (q/to-axis-angle q/identity-q)]
      (is (zero? angle))))

  (testing "angle measures the turn, the short way round"
    (is (< (abs (- 1.0 (q/angle (q/from-axis-angle [0 0 1] 1.0)))) 1e-9))
    (is (< (abs (- (q/angle (q/from-axis-angle [0 0 1] (- 1.0))) 1.0)) 1e-9)
        "and does not care about the sign of the turn")
    (is (zero? (q/angle q/identity-q)))))

(deftest slerp-test
  (let [a (q/from-axis-angle [0 0 1] 0.0)
        b (q/from-axis-angle [0 0 1] 2.0)]
    (testing "the ends are the ends"
      (is (close? (q/rotate (q/slerp a b 0.0) [1 0 0]) (q/rotate a [1 0 0])))
      (is (close? (q/rotate (q/slerp a b 1.0) [1 0 0]) (q/rotate b [1 0 0]))))

    (testing "and it moves at a constant rate between them"
      (is (< (abs (- 1.0 (q/angle (q/slerp a b 0.5)))) 1e-9))
      (is (< (abs (- 0.5 (q/angle (q/slerp a b 0.25)))) 1e-9))))

  (testing "nearly identical orientations do not divide by zero"
    (let [a (q/from-axis-angle [0 0 1] 1.0)
          b (q/from-axis-angle [0 0 1] 1.0000000001)]
      (is (v/finite? (vec (take 3 (q/slerp a b 0.5)))))))

  (testing "it takes the short way round"
    ;; A quaternion and its negation are the same orientation, so without
    ;; picking the nearer one the path can go the long way.
    (let [a (q/from-axis-angle [0 0 1] 0.0)
          b (mapv - (q/from-axis-angle [0 0 1] 0.4))]
      (is (< (q/angle (q/slerp a b 0.5)) 0.3)))))
