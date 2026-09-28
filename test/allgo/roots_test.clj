(ns allgo.roots-test
  "Bracketing root-finders, event location and golden-section search."
  (:require [allgo.numerics.roots :as roots]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest bisection
  (is (< (abs (- (roots/bisect #(- (* % %) 2.0) 0.0 2.0) (math/sqrt 2.0))) 1e-15))
  (testing "either way round"
    (is (< (abs (- (roots/bisect math/cos 3.0 0.0) (/ math/PI 2))) 1e-15))))

(deftest events
  (testing "each change of a condition, found to the tolerance"
    (let [ts (roots/transitions #(pos? (math/sin %)) 0.5 10.0 0.3 1e-9)]
      ;; each reported with the value it had before: true leaving at pi,
      ;; false entering at 2 pi, true leaving at 3 pi
      (is (= [true false true] (map second ts)))
      (is (every? true? (map #(< (abs (- (first %1) %2)) 1e-9) ts [math/PI (* 2 math/PI) (* 3 math/PI)]))))))

(deftest golden-section
  (testing "a minimum and a maximum, with the bracket either way round"
    (is (< (abs (- (roots/minimize #(math/pow (- % 1.3) 2) 0.0 4.0) 1.3)) 1e-8))
    (is (< (abs (- (roots/minimize #(math/pow (- % -1.3) 2) 0.0 -4.0 {:tol 0.0 :max-iter 100}) -1.3)) 1e-8))
    ;; a peak is flat, so its place is found only to sqrt(epsilon)
    (is (< (abs (- (roots/maximize math/sin 0.0 3.0) (/ math/PI 2))) 1e-7))))
