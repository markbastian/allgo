(ns allgo.learning.autodiff-test
  (:require [allgo.learning.autodiff :as ad]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Reverse-mode differentiation, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapters 21 and 24:
;; every operation's gradient against central differences.

(defn- numeric-grad [f params k i j]
  (let [h 1e-6
        at (fn [d] (:value (ad/value-and-grad f (update-in params [k i j] + d))))]
    (/ (- (at h) (at (- h))) (* 2 h))))

(defn- check [f params]
  (let [{:keys [grads]} (ad/value-and-grad f params)]
    (doseq [[k M] params i (range (count M)) j (range (count (first M)))]
      (is (< (abs (- (get-in grads [k i j]) (numeric-grad f params k i j))) 1e-6) [k i j]))))

(def ^:private rng (random/rng 7))
(defn- m [r c] (ad/init r c rng))

(deftest gradient-test
  (testing "linear operations, nonlinearities and reductions"
    (check (fn [{:keys [A B r]}]
             (ad/sum (ad/tanh (ad/add-row (ad/matmul A B) r))))
           {:A (m 3 4) :B (m 4 2) :r (m 1 2)})
    (check (fn [{:keys [A B]}]
             (ad/mean (ad/mul (ad/sigmoid A) (ad/exp (ad/sub B (ad/scale A 0.5))))))
           {:A (m 2 3) :B (m 2 3)})
    (check (fn [{:keys [A]}] (ad/sum (ad/mul (ad/relu A) (ad/transpose (ad/transpose A)))))
           {:A [[0.5 -0.3 0.2] [-0.7 0.9 0.1]]}))
  (testing "structure: lookup, stacking, concatenation and slicing"
    (check (fn [{:keys [E W]}]
             (let [x (ad/rows E [2 0 2])
                   h (ad/concat-cols [x (ad/slice-cols x 1 3)])]
               (ad/sum (ad/tanh (ad/matmul (ad/stack-rows [h (ad/row h 1)]) W)))))
           {:E (m 4 3) :W (m 5 2)}))
  (testing "softmax, layer norm with gain and bias, and cross-entropy"
    (check (fn [{:keys [A g b]}]
             (ad/sum (ad/mul (ad/softmax-rows A) (ad/add-row (ad/mul-row (ad/layer-norm-rows A) g) b))))
           {:A (m 3 4) :g (m 1 4) :b (m 1 4)})
    (check (fn [{:keys [A W]}] (ad/cross-entropy (ad/matmul A W) [1 0 2]))
           {:A (m 3 4) :W (m 4 3)})))

(deftest training-test
  (testing "Adam on cross-entropy fits a small classification, the loss falling"
    (let [X [[1.0 0.0] [0.0 1.0] [1.0 1.0] [0.0 0.0]]
          f (fn [{:keys [W b]}] (ad/cross-entropy (ad/add-row (ad/matmul (ad/const X) W) b) [0 1 1 0]))
          losses (map first (iterate (fn [[_ params state]]
                                       (let [{:keys [value grads]} (ad/value-and-grad f params)
                                             [p s] (ad/adam-step params grads state :lr 0.1)]
                                         [value p s]))
                                     [nil {:W (m 2 2) :b (m 1 2)} nil]))
          ls (vec (take 200 (rest losses)))]
      (is (< (peek ls) (* 0.2 (first ls))))))
  (testing "clipping scales every gradient by the same factor, to the norm asked"
    (let [g (ad/clip {:a [[3.0 0.0]] :b [[0.0 4.0]]} 1.0)]
      (is (< (abs (- 0.6 (get-in g [:a 0 0]))) 1e-12))
      (is (< (abs (- 0.8 (get-in g [:b 0 1]))) 1e-12)))))
