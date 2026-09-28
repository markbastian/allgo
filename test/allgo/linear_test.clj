(ns allgo.linear-test
  "The small dense helpers: the 3x3 determinant and solve, and carrying a
  matrix through a linear map."
  (:require [allgo.numerics.linear :as lin]
            [clojure.test :refer [deftest is testing]]))

(def ^:private m [[2.0 -1.0 0.5] [0.3 4.0 -2.0] [1.0 1.0 3.0]])

(deftest three-by-three
  (testing "the determinant is the triple product of the rows"
    (let [[a b c] m
          cross (fn [[x y z] [u v w]] [(- (* y w) (* z v)) (- (* z u) (* x w)) (- (* x v) (* y u))])]
      (is (< (abs (- (lin/det-3 m) (lin/dot a (cross b c)))) 1e-12))))
  (testing "and the solve satisfies the system"
    (let [b [1.0 -2.0 0.25]
          x (lin/solve-3 m b)]
      (is (every? #(< (abs %) 1e-12) (lin/sub (lin/mat-vec m x) b))))))

(deftest congruence
  (testing "A B A^T: a rotation keeps the trace, and a scaling scales it"
    (let [B [[4.0 1.0 0.0] [1.0 2.0 0.5] [0.0 0.5 1.0]]
          c 0.6 s 0.8
          R [[c (- s) 0.0] [s c 0.0] [0.0 0.0 1.0]]
          tr (fn [M] (+ (get-in M [0 0]) (get-in M [1 1]) (get-in M [2 2])))]
      (is (< (abs (- (tr (lin/congruence R B)) (tr B))) 1e-12))
      (is (< (abs (- (tr (lin/congruence (lin/mat-scale (lin/eye 3) 2.0) B)) (* 4.0 (tr B)))) 1e-12)))))

(deftest general-inverse
  (testing "an unsymmetric matrix, and a zero on the diagonal that pivoting must step round"
    (let [A [[0.0 2.0 1.0] [1.0 -1.0 4.0] [3.0 0.5 -2.0]]]
      (is (every? #(< (abs %) 1e-12) (flatten (lin/mat-sub (lin/mat-mul A (lin/inverse-general A)) (lin/eye 3)))))))
  (testing "and nil for a singular one"
    (is (nil? (lin/inverse-general [[1.0 2.0] [2.0 4.0]])))))
