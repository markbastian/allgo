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

;; Dropout, section 21.5.4.

(deftest dropout-test
  (let [rng (random/rng 4)
        net (nn/with-dropout [(nn/dense 6 40 :tanh rng) (nn/dense 40 30 :relu rng) (nn/dense 30 2 :sigmoid rng)] 0.5 0.8)]
    (testing "every unit kept with probability p and its output multiplied by 1/p, else zero: hidden p = 0.5, input 0.8"
      (let [draws (repeatedly 400 #(nn/dropout-masks net rng))
            kept (fn [ms] (/ (count (filter pos? (apply concat ms))) (double (count (apply concat ms)))))]
        (is (= #{0.0 2.0} (set (mapcat #(get-in % [:layers 0]) draws))))
        (is (= #{0.0 1.25} (set (mapcat :input draws))))
        (is (< (abs (- 0.5 (kept (map #(get-in % [:layers 1]) draws)))) 0.02))
        (is (< (abs (- 0.8 (kept (map :input draws)))) 0.02))
        (testing "and never the output layer"
          (is (every? nil? (map #(get-in % [:layers 2]) draws))))))
    (testing "back-propagation through the thinned network: its gradients the numerical ones"
      (let [{:keys [layers]} (nn/dropout-masks net (random/rng 5))
            x (vec (repeatedly 6 rng)) y [1.0 0.0]
            f (fn [n] (let [a (second (last (nn/forward n x layers)))] (reduce + (map #(let [d (- %1 %2)] (* d d)) a y))))
            g (nn/gradients net x y :squared layers)
            kept (first (keep-indexed #(when (pos? %2) %1) (layers 0)))
            dropped (first (keep-indexed #(when (zero? %2) %1) (layers 0)))]
        (doseq [[path gpath] [[[0 :W kept 3] [0 :dW kept 3]] [[0 :b kept] [0 :db kept]] [[1 :W 4 kept] [1 :dW 4 kept]]
                              [[1 :b 7] [1 :db 7]] [[2 :W 1 9] [2 :dW 1 9]]]]
          (is (< (abs (- (get-in g gpath) (let [h 1e-6] (/ (- (f (update-in net path + h)) (f (update-in net path - h))) (* 2 h)))))
                 1e-7)
              path))
        (testing "and a dropped unit learns nothing that step"
          (is (every? zero? (get-in g [0 :dW dropped])))
          (is (zero? (get-in g [0 :db dropped]))))))
    (testing "multiplying by 1/p: in a linear network, the thinned networks' outputs average to the whole network's"
      (let [linear (mapv #(assoc % :g :identity) net)
            x (vec (repeatedly 6 rng))
            mean (/ (reduce + (repeatedly 20000 #(let [{:keys [input layers]} (nn/dropout-masks linear rng)]
                                                   (first (second (last (nn/forward linear (mapv * x input) layers)))))))
                    20000)]
        (is (< (abs (- mean (first (nn/predict linear x)))) 0.01))))
    (testing "p = 1 keeps every unit: training exactly as without dropout"
      (let [X (vec (repeatedly 10 #(vec (repeatedly 6 rng)))) Y (vec (repeatedly 10 #(vector (rng) (rng))))
            plain (mapv #(dissoc % :keep :keep-input) net)
            opts {:alpha 0.1 :rng (random/rng 6)}]
        (is (= (nth (nn/train-seq plain X Y opts) 3)
               (mapv #(dissoc % :keep :keep-input) (nth (nn/train-seq (nn/with-dropout plain 1.0 1.0) X Y (assoc opts :rng (random/rng 6))) 3))))))))

(deftest dropout-fit-test
  (testing "dropout makes the training set harder to fit: on 50 noisy points, a higher training loss"
    (let [rng (random/rng 1)
          label (fn [[x y]] (if (< (+ (* x x) (* y y)) 0.5) 1.0 0.0))
          X (vec (repeatedly 50 (fn [] [(- (* 2 (rng)) 1) (- (* 2 (rng)) 1)])))
          Y (mapv (fn [x] [(if (< (rng) 0.2) (- 1.0 (label x)) (label x))]) X)
          net [(nn/dense 2 20 :tanh rng) (nn/dense 20 20 :tanh rng) (nn/dense 20 1 :sigmoid rng)]
          fit (fn [n] (nn/loss (nth (nn/train-seq n X Y {:alpha 0.1 :loss :cross-entropy :rng (random/rng 11)}) 100) X Y :cross-entropy))]
      (is (< (fit net) (fit (nn/with-dropout net 0.5)))))))
