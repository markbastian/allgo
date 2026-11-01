(ns allgo.learning.neural-test
  (:require [allgo.learning.neural :as nn]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Back-propagation, after Russell and Norvig, *Artificial Intelligence:
;; A Modern Approach*, 4th edition, chapter 21.

(defn- numeric-gradient [net x y loss path]
  (let [h 1e-6 f (fn [n] (nn/loss n [x] [y] loss))]
    (/ (- (f (update-in net path + h)) (f (update-in net path - h))) (* 2 h))))

(deftest gradient-test
  (testing "back-propagated gradients agree with numerical ones, through a convolution and two dense layers"
    (let [rng (random/rng 2)
          net [(nn/conv [4 4] 2 2 :tanh rng) (nn/dense 18 3 :relu rng) (nn/dense 3 2 :sigmoid rng)]
          x (vec (repeatedly 16 rng))
          y [1.0 0.0]]
      (doseq [loss [:squared :cross-entropy]
              :let [g (nn/gradients net x y loss)]
              [path gpath] [[[0 :kernels 1 0 1] [0 :dK 1 0 1]] [[0 :kernels 0 1 1] [0 :dK 0 1 1]] [[0 :b 0] [0 :db 0]]
                            [[1 :W 2 5] [1 :dW 2 5]] [[2 :W 1 2] [2 :dW 1 2]] [[2 :b 1] [2 :db 1]]]]
        (is (< (abs (- (get-in g gpath) (numeric-gradient net x y loss path))) 1e-7) [loss path])))))

(deftest xor-test
  (testing "a hidden layer learns XOR, which no single layer can"
    (let [rng (random/rng 7)
          X [[0.0 0.0] [0.0 1.0] [1.0 0.0] [1.0 1.0]]
          Y [[0.0] [1.0] [1.0] [0.0]]
          net (nth (nn/train-seq [(nn/dense 2 4 :tanh rng) (nn/dense 4 1 :sigmoid rng)] X Y
                                 {:alpha 0.05 :loss :cross-entropy :rng (random/rng 3)})
                   3000)]
      (is (= [0 1 1 0] (map #(Math/round (double (first (nn/predict net %)))) X))))))
