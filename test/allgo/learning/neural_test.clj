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

;; Pooling, section 21.3.1.

(deftest pooling-test
  (testing "the chapter's example: max-pooling [5, 9, 4] gives 9 -- the feature is somewhere in the window"
    (is (= [9] (nn/predict [(nn/pool [1 3] 1 [1 3] 1 :max)] [5 9 4])))
    (is (= [6] (nn/predict [(nn/pool [1 3] 1 [1 3] 1 :average)] [5 9 4]))))
  (testing "average-pooling is a convolution with the uniform kernel [1/l ... 1/l]"
    (let [rng (random/rng 1)
          x (vec (repeatedly 36 rng))
          uniform (assoc (nn/conv [6 6] 3 1 :identity rng) :kernels [(vec (repeat 3 (vec (repeat 3 (/ 1.0 9)))))])]
      (is (every? #(< (abs %) 1e-12) (map - (nn/predict [uniform] x) (nn/predict [(nn/pool [6 6] 1 3 1 :average)] x))))))
  (testing "with s = l, downsampling by s: 10s pixels by 10s become 10 by 10, in each feature map"
    (is (= 200 (count (nn/predict [(nn/pool [30 30] 2 3 3 :max)] (vec (range 1800)))))))
  (testing "back-propagated gradients through both kinds agree with numerical ones"
    (let [rng (random/rng 2)
          net [(nn/conv [6 6] 3 2 :tanh rng) (nn/pool [4 4] 2 2 2 :max) (nn/dense 8 4 :tanh rng)
               (nn/pool [1 4] 1 [1 2] 1 :average) (nn/dense 3 2 :sigmoid rng)]
          x (vec (repeatedly 36 rng)) y [1.0 0.0]
          g (nn/gradients net x y :cross-entropy)]
      (doseq [[path gpath] [[[0 :kernels 1 0 1] [0 :dK 1 0 1]] [[0 :b 0] [0 :db 0]] [[2 :W 3 5] [2 :dW 3 5]] [[4 :W 1 2] [4 :dW 1 2]]]]
        (is (< (abs (- (get-in g gpath) (numeric-gradient net x y :cross-entropy path))) 1e-7) path)))))

;; Residual networks, section 21.3.3.

(defn- norm [M] (Math/sqrt (reduce + (for [r M x r] (* x x)))))

(deftest residual-test
  (let [rng (random/rng 5)]
    (testing "V = 0, with relus: the layer passes its relu inputs through unchanged"
      (let [z (vec (repeatedly 6 #(max 0.0 (- (rng) 0.5))))
            layer (assoc (nn/residual 6 4 rng) :V (vec (repeat 6 (vec (repeat 4 0.0)))))]
        (is (= z (nn/predict [layer] z)))))
    (testing "back-propagated gradients agree with numerical ones"
      (let [net [(nn/dense 3 5 :relu rng) (nn/residual 5 4 rng) (nn/residual 5 2 rng) (nn/dense 5 1 :identity rng)]
            x [0.3 0.9 0.5] y [0.2]
            g (nn/gradients net x y :squared)]
        (doseq [[path gpath] [[[0 :W 2 1] [0 :dW 2 1]] [[1 :W 3 4] [1 :dW 3 4]] [[1 :bw 1] [1 :dbw 1]]
                              [[1 :V 2 3] [1 :dV 2 3]] [[2 :b 0] [2 :db 0]] [[3 :W 0 4] [3 :dW 0 4]]]]
          (is (< (abs (- (get-in g gpath) (numeric-gradient net x y :squared path))) 1e-7) path))))
    (testing "thirty layers deep, the first layer's gradient vanishes through plain sigmoid or relu layers, and a residual network keeps it"
      (let [x (vec (repeatedly 8 rng))
            deep (fn [layer] (vec (concat [(nn/dense 8 8 :relu rng)] (repeatedly 30 layer) [(nn/dense 8 1 :identity rng)])))
            first-grad (fn [net] (norm (:dW (first (nn/gradients net x [1.0] :squared)))))]
        (is (< (first-grad (deep #(nn/dense 8 8 :sigmoid rng))) 1e-9))
        (is (< (first-grad (deep #(nn/dense 8 8 :relu rng))) 1e-9))
        (is (> (first-grad (deep #(nn/residual 8 8 rng))) 0.1))))))

;; Batch normalization, section 21.4.2.

(defn- sd [xs] (let [m (/ (reduce + xs) (count xs))] (Math/sqrt (/ (reduce + (map #(* (- % m) (- % m)) xs)) (count xs)))))
(defn- shrink [layer k] (update layer :W (fn [W] (mapv (fn [r] (mapv #(* k %) r)) W))))

(deftest batch-norm-test
  (let [rng (random/rng 8)
        X (vec (repeatedly 64 (fn [] (vec (repeatedly 10 #(- (* 2 (rng)) 1))))))]
    (testing "in training, each unit's values over the minibatch have mean beta and standard deviation gamma (less a hair for epsilon)"
      (let [layer (assoc (nn/batch-norm 10 :identity) :gamma (vec (range 1.0 11.0)) :beta (vec (range -5.0 5.0)))
            out (second (last (nn/forward-batch [(nn/dense 10 10 :identity rng) layer] X nil)))]
        (doseq [j (range 10) :let [col (map #(% j) out)]]
          (is (< (abs (- (/ (reduce + col) 64) (- j 5.0))) 1e-9))
          (is (< (abs (- (/ (sd col) (inc j)) 1.0)) 1e-3)))))
    (testing "back-propagated gradients -- through every example's dependence on the others -- agree with numerical ones"
      (let [net [(nn/dense 10 6 :identity rng) (assoc (nn/batch-norm 6 :tanh) :gamma [1.2 0.7 1.0 0.9 1.5 1.1] :beta [0.1 0 0 -0.2 0 0.3])
                 (nn/dense 6 2 :sigmoid rng)]
            X (subvec X 0 8) Y (vec (repeatedly 8 (fn [] [(rng) (rng)])))
            g (nn/batch-gradients net X Y :cross-entropy)
            f #(nn/batch-loss % X Y :cross-entropy)]
        (doseq [[path gpath] [[[0 :W 3 4] [0 :dW 3 4]] [[0 :b 2] [0 :db 2]] [[1 :gamma 2] [1 :dgamma 2]]
                              [[1 :beta 5] [1 :dbeta 5]] [[2 :W 1 3] [2 :dW 1 3]]]]
          (is (< (abs (- (get-in g gpath) (let [h 1e-6] (/ (- (f (update-in net path + h)) (f (update-in net path - h))) (* 2 h)))))
                 1e-7)
              path))))
    (testing "weights too small: without it the spread of a deep layer's values decays to near zero; with it, it does not"
      (let [layers (vec (repeatedly 8 #(shrink (nn/dense 10 10 :identity rng) 0.3)))
            spread (fn [acts] (sd (map first (second acts))))]
        (is (< (spread (last (nn/forward-batch (mapv #(assoc % :g :tanh) layers) X nil))) 1e-5))
        (is (< 0.3 (spread (last (nn/forward-batch (vec (mapcat (fn [l] [l (nn/batch-norm 10 :tanh)]) layers)) X nil)))))))
    (testing "the rate of convergence: six small-weighted tanh layers learn with it what they barely start to without; and at test time, by the running averages, it classifies"
      (let [label (fn [[a b]] (if (pos? (* a b)) 1.0 0.0))
            X (vec (repeatedly 256 (fn [] [(- (* 2 (rng)) 1) (- (* 2 (rng)) 1)])))
            Y (mapv #(vector (label %)) X)
            hidden (vec (cons (shrink (nn/dense 2 16 :identity rng) 0.3) (repeatedly 5 #(shrink (nn/dense 16 16 :identity rng) 0.3))))
            out (nn/dense 16 1 :sigmoid rng)
            train #(nth (nn/train-seq % X Y {:alpha 0.1 :batch 16 :loss :cross-entropy :rng (random/rng 2)}) 40)
            plain (train (conj (mapv #(assoc % :g :tanh) hidden) out))
            normed (train (conj (vec (mapcat (fn [l] [l (nn/batch-norm 16 :tanh)]) hidden)) out))]
        (is (> (nn/loss plain X Y :cross-entropy) 0.65))
        (is (< (nn/loss normed X Y :cross-entropy) 0.45))
        (is (> (count (filter true? (map #(= (Math/round (double (first (nn/predict normed %1)))) (long (first %2))) X Y))) 220))))))

;; Weight decay, section 21.5.3.

(deftest weight-decay-test
  (let [rng (random/rng 7)
        xs (vec (repeatedly 40 #(* 4 (rng))))
        X (mapv vector xs)
        Y (mapv (fn [x] [(+ 1.0 (* 2.0 x) (random/gaussian rng 0 1))]) xs)
        mean (fn [vs] (/ (reduce + vs) (count vs)))
        mx (mean xs) my (mean (map first Y))
        vx (mean (map #(* (- % mx) (- % mx)) xs))
        cxy (mean (map #(* (- %1 mx) (- (first %2) my)) xs Y))]
    (testing "one linear unit, the penalty lambda w^2 on its weight but not its bias: gradient descent finds ridge regression's w = cov(x, y) / (var x + lambda)"
      (doseq [lambda [0.0 0.5 2.0]]
        (let [net (nth (nn/train-seq [(nn/dense 1 1 :identity rng)] X Y {:alpha 0.05 :batch 40 :lambda lambda}) 3000)]
          (is (< (abs (- (get-in net [0 :W 0 0]) (/ cxy (+ vx lambda)))) 1e-9) lambda))))
    (testing "the penalty is lambda times the sum of the squared weights, and lambda = 0 is no weight decay"
      (let [net [(nn/dense 1 3 :tanh rng) (nn/dense 3 1 :identity rng)]]
        (is (< (abs (- (nn/weight-penalty net 0.1) (* 0.1 (reduce + (for [l net r (:W l) w r] (* w w)))))) 1e-15))
        (is (= (nth (nn/train-seq net X Y {:alpha 0.01}) 2)
               (nth (nn/train-seq net X Y {:alpha 0.01 :lambda 0.0}) 2)))))))
