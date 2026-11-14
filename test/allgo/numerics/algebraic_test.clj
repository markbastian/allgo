(ns allgo.numerics.algebraic-test
  "Integer polynomials against their factors multiplied out by hand, and
  resultants against eliminating y by hand."
  (:require [allgo.numerics.algebraic :as alg]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private x-1-x-2-x-3 [1.0 -6.0 11.0 -6.0])

(deftest arithmetic-test
  (testing "(x - 1)(x - 2)(x - 3) multiplied out, and divided back"
    (is (= x-1-x-2-x-3 (alg/mul (alg/mul [1.0 -1.0] [1.0 -2.0]) [1.0 -3.0])))
    (is (= [1.0 -3.0 2.0] (alg/quotient x-1-x-2-x-3 [1.0 -3.0]))))
  (testing "a division that does not come out even throws"
    (is (thrown? clojure.lang.ExceptionInfo (alg/quotient x-1-x-2-x-3 [1.0 -4.0]))))
  (testing "the pseudo-remainder of x^3 + 1 by 2x + 1 is 2^3 (-1/8 + 1) = 7"
    (is (= [7.0] (alg/prem [1.0 0.0 0.0 1.0] [2.0 1.0]))))
  (testing "a product past 2^53 throws rather than round"
    (is (thrown? clojure.lang.ExceptionInfo (alg/mul [1e10 1.0] [1e10 1.0])))))

(deftest gcd-test
  (testing "the common factor of (x - 1)(x - 2) and (x - 2)(x + 5), and none of x^2 + 1 and x - 1"
    (is (= [1.0 -2.0] (alg/gcd [1.0 -3.0 2.0] [1.0 3.0 -10.0])))
    (is (= [1.0] (alg/gcd [1.0 0.0 1.0] [1.0 -1.0]))))
  (testing "content and primitive part"
    (is (= 6.0 (alg/content [12.0 -18.0 6.0])))
    (is (= [2.0 -3.0 1.0] (alg/primitive [-12.0 18.0 -6.0]))))
  (testing "the square-free part of (x - 1)^2 (x + 2)^3 is (x - 1)(x + 2)"
    (let [p (reduce alg/mul [[1.0 -1.0] [1.0 -1.0] [1.0 2.0] [1.0 2.0] [1.0 2.0]])]
      (is (= [1.0 1.0 -2.0] (alg/square-free p))))))

(deftest roots-test
  (testing "Sturm's theorem counts x^3 - 6x^2 + 11x - 6's roots in an interval"
    (let [chain (alg/sturm x-1-x-2-x-3)]
      (is (= 3 (alg/count-roots chain -10 10)))
      (is (= 2 (alg/count-roots chain 1.5 3.5)))
      (is (= 0 (alg/count-roots chain 3.5 10)))))
  (testing "integer roots exactly, a double root once, an irrational one to the last bits, none of x^2 + 1"
    (is (= [1.0 2.0 3.0] (alg/real-roots x-1-x-2-x-3)))
    (is (= [-2.0 1.0] (alg/real-roots (reduce alg/mul [[1.0 -1.0] [1.0 -1.0] [1.0 2.0]]))))
    (is (< (abs (- (math/sqrt 2) (second (alg/real-roots [1.0 0.0 -2.0])))) 1e-15))
    (is (= [] (alg/real-roots [1.0 0.0 1.0]))))
  (testing "two roots a thousandth apart, each isolated: (1000 x - 1414)(1000 x - 1415)"
    (let [p (alg/mul [1000.0 -1414.0] [1000.0 -1415.0])]
      (is (= 2 (count (alg/real-roots p))))
      (is (< (abs (- 1.414 (first (alg/real-roots p)))) 1e-12))
      (is (< (abs (- 1.415 (second (alg/real-roots p)))) 1e-12)))))

(deftest resultant-test
  (testing "eliminating y from x^2 + y^2 - 1 and y - x leaves 2x^2 - 1"
    (is (= [2.0 0.0 -1.0] (alg/resultant [[1.0] [] [1.0 0.0 -1.0]] [[1.0] [-1.0 0.0]]))))
  (testing "the circle's discriminant in y, res(y^2 + x^2 - 1, 2y) = 4(x^2 - 1): zero where it folds, at x = +-1"
    (is (= [-1.0 1.0] (alg/real-roots (alg/resultant [[1.0] [] [1.0 0.0 -1.0]] [[2.0] []])))))
  (testing "a 3 x 3 integer determinant, with a row swap"
    (is (= [-3.0] (alg/determinant [[[] [1.0] [2.0]] [[3.0] [4.0] [5.0]] [[6.0] [7.0] [9.0]]])))))
