(ns allgo.special-test
  "Special functions against their definitions and closed forms."
  (:require [allgo.numerics.special :as special]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest bessel
  (testing "I_n against its integral, (1/pi) int_0^pi exp(x cos t) cos(n t) dt"
    (doseq [n [0 1 2 3] x [0.0 0.5 3.0 12.0]]
      (let [k 2000 h (/ math/PI k)
            integral (* (/ h math/PI)
                        (reduce + (for [j (range (inc k))
                                        :let [t (* j h) w (if (or (zero? j) (= j k)) 0.5 1.0)]]
                                    (* w (math/exp (* x (math/cos t))) (math/cos (* n t))))))]
        (is (< (abs (- (special/bessel-i n x) integral)) (* 1e-10 (max 1.0 integral))) (str [n x]))))))

(deftest legendre
  (let [x 0.37 s (math/sqrt (- 1.0 (* x x)))
        P (special/associated-legendre x s 3)
        A (special/derived-legendre x 4)]
    (testing "the associated functions against their closed forms"
      (doseq [[k v] {[2 0] (* 0.5 (- (* 3 x x) 1.0)) [2 1] (* 3 x s) [2 2] (* 3 s s)
                     [3 1] (* 1.5 s (- (* 5 x x) 1.0)) [3 3] (* 15 s s s)}]
        (is (< (abs (- (P k) v)) 1e-14) (str k))))
    (testing "the derived functions are the polynomials' derivatives: P_nm = s^m A_nm"
      (doseq [n (range 4) m (range (inc n))]
        (is (< (abs (- (P [n m]) (* (math/pow s m) (A [n m])))) 1e-13) (str [n m]))))))
