(ns allgo.vec3-test
  (:require [allgo.geometry.vec3 :as v]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b] (< (v/distance a b) 1e-12))

(deftest arithmetic-test
  (testing "the basics"
    (is (= [4.0 6.0 8.0] (v/add [1 2 3] [3.0 4.0 5.0])))
    (is (= [-2 -2 -2] (v/sub [1 2 3] [3 4 5])))
    (is (= [2.0 4.0 6.0] (v/scale [1 2 3] 2)))
    (is (= [-1 -2 -3] (v/negate [1 2 3])))
    (is (= [3.0 8.0 15.0] (v/mul [1.0 2.0 3.0] [3.0 4.0 5.0]))))

  (testing "add takes more than two"
    (is (= [3.0 3.0 3.0] (v/add [1 1 1] [1.0 1.0 1.0] [1.0 1.0 1.0]))))

  (testing "add-scaled is the integrator step"
    (is (= [1.0 2.0 3.0] (v/add-scaled [1 2 3] [4 5 6] 0.0)))
    (is (= [5.0 7.0 9.0] (v/add-scaled [1 2 3] [4 5 6] 1.0)))))

(deftest products-test
  (testing "dot"
    (is (== 32 (v/dot [1 2 3] [4 5 6])))
    (is (zero? (v/dot [1 0 0] [0 1 0])) "perpendicular vectors give nothing"))

  (testing "cross follows the right hand"
    (is (= [0 0 1] (v/cross [1 0 0] [0 1 0])))
    (is (= [0 0 -1] (v/cross [0 1 0] [1 0 0]))))

  (testing "a vector crossed with itself is zero"
    ;; Relied on in `allgo.physics.rigid`: it is why moving a body before
    ;; measuring its lever arm makes no difference.
    (is (= [0 0 0] (v/cross [1 2 3] [1 2 3])))
    (is (close? [0.0 0.0 0.0] (v/cross [1 2 3] (v/scale [1 2 3] 7.5)))))

  (testing "cross is perpendicular to both"
    (let [a [1.0 -2.0 0.5] b [0.3 4.0 -1.0] c (v/cross a b)]
      (is (< (abs (v/dot c a)) 1e-12))
      (is (< (abs (v/dot c b)) 1e-12)))))

(deftest length-test
  (testing "length and its square"
    (is (== 25 (v/length-squared [3 4 0])))
    (is (== 5.0 (v/length [3 4 0])))
    (is (== 5.0 (v/distance [1 1 0] [4 5 0]))))

  (testing "normalize gives a unit vector"
    (is (< (abs (- 1.0 (v/length (v/normalize [3 -4 12])))) 1e-12)))

  (testing "normalizing nothing gives nothing rather than NaN"
    (is (= v/zero (v/normalize [0 0 0])))
    (is (v/finite? (v/normalize [0 0 0]))))

  (testing "zero-vector? and finite?"
    (is (v/zero-vector? [0 0 0]))
    (is (not (v/zero-vector? [0 0 1e-9])))
    (is (v/finite? [1 2 3]))
    (is (not (v/finite? [1 ##Inf 3])))
    (is (not (v/finite? [##NaN 0 0])))))

(deftest lerp-test
  (testing "lerp hits both ends and the middle"
    (is (= [0.0 0.0 0.0] (v/lerp [0 0 0] [2 4 6] 0.0)))
    (is (= [2.0 4.0 6.0] (v/lerp [0 0 0] [2 4 6] 1.0)))
    (is (= [1.0 2.0 3.0] (v/lerp [0 0 0] [2 4 6] 0.5)))))
