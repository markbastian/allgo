(ns allgo.learning.models-test
  (:require [allgo.learning.models :as mo]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Models of Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, chapters 19 and 20, on data made to have a
;; known answer.

(def ^:private rng (random/rng 3))
(def ^:private X (vec (repeatedly 200 (fn [] [(* 4 (- (rng) 0.5)) (* 4 (- (rng) 0.5))]))))
(def ^:private linear-labels (mapv (fn [[a b]] (if (> (+ a (* 2 b) -0.5) 0) 1 0)) X))
(def ^:private circle-labels (mapv (fn [[a b]] (if (< (+ (* a a) (* b b)) 1.5) 1 -1)) X))

(defn- accuracy [h labels] (/ (count (filter true? (map #(= (h %1) %2) X labels))) (double (count X))))

(deftest regression-test
  (testing "least squares recovers y = 1 + 2 x1 - 3 x2 through noise of 0.05"
    (let [ys (mapv (fn [[a b]] (+ 1.0 (* 2 a) (* -3 b) (* 0.1 (- (rng) 0.5)))) X)
          w (mo/linear-regression X ys)]
      (is (every? true? (map #(< (abs (- %1 %2)) 0.01) w [1.0 2.0 -3.0])))
      (testing "and gradient descent on the squared loss reaches the same weights"
        (let [gd (nth (mo/gradient-descent-seq X ys 0.1 (fn [w x y] (let [e (- (reduce + (map * w x)) y)] (mapv #(* 2.0 e %) x)))) 500)]
          (is (every? true? (map #(< (abs (- %1 %2)) 1e-6) w gd))))))))

(deftest classifiers-test
  (testing "on linearly separable data, logistic regression and the perceptron are perfect"
    (let [lw (nth (mo/logistic-regression-seq X linear-labels 1.0) 2000)
          pw (nth (mo/perceptron-seq X linear-labels 0.1) 200)]
      (is (= 1.0 (accuracy #(if (> (mo/predict-logistic lw %) 0.5) 1 0) linear-labels)))
      (is (= 1.0 (accuracy #(mo/predict-threshold pw %) linear-labels)))))
  (testing "k nearest neighbors and Gaussian naive Bayes do well"
    (is (> (accuracy (mo/knn X linear-labels 5) linear-labels) 0.95))
    (is (> (accuracy (mo/gaussian-naive-bayes X linear-labels) linear-labels) 0.85)))
  (testing "an SVM with a Gaussian kernel learns a circle a linear one cannot"
    (let [rbf (mo/svm-train X circle-labels {:C 10.0 :kernel (mo/gaussian-kernel 1.0) :rng (random/rng 1)})
          lin (mo/svm-train X circle-labels {:C 10.0 :rng (random/rng 1)})]
      (is (> (accuracy #(mo/svm-predict rbf %) circle-labels) 0.95))
      (is (< (accuracy #(mo/svm-predict lin %) circle-labels) 0.85))
      (is (< (count (mo/support-vectors rbf)) 60) "only some examples are support vectors"))))
