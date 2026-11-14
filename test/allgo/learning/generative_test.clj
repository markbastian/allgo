(ns allgo.learning.generative-test
  (:require [allgo.learning.generative :as gen]
            [allgo.learning.neural :as nn]
            [allgo.numerics.eigen :as eigen]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Unsupervised deep learning, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 21.7.1.

(defn- mean [xs] (/ (reduce + xs) (count xs)))
(defn- sd [xs] (let [m (mean xs)] (Math/sqrt (mean (map #(* (- % m) (- % m)) xs)))))
(defn- max-abs [M] (reduce max (map abs (flatten M))))

(def ^:private planar
  "500 points in 4-D, nearly all their variance in a plane."
  (let [rng (random/rng 1)]
    (vec (repeatedly 500 (fn [] (let [a (random/gaussian rng 0 3) b (random/gaussian rng 0 1.5)]
                                  (mapv + [1 2 3 4] (lin/scale [1 1 0 0.5] a) (lin/scale [0 1 -1 0] b)
                                        (repeatedly 4 #(random/gaussian rng 0 0.3)))))))))

(def ^:private ring
  "400 points near the circle of radius 2."
  (let [rng (random/rng 2)]
    (vec (repeatedly 400 (fn [] (let [t (* 2 Math/PI (rng))]
                                  [(+ (* 2 (Math/cos t)) (random/gaussian rng 0 0.05)) (+ (* 2 (Math/sin t)) (random/gaussian rng 0 0.05))]))))))

(deftest probabilistic-pca-test
  (let [S (gen/covariance planar)
        {:keys [values vectors]} (eigen/jacobi S)
        models (gen/ppca-em-seq planar 2 (random/rng 2))]
    (testing "EM never lowers the likelihood"
      (is (apply <= (map #(gen/ppca-log-likelihood % planar) (take 40 models)))))
    (testing "and converges to the maximum-likelihood solution: sigma^2 the mean of the discarded eigenvalues of the covariance, W W^T + sigma^2 I its two largest eigenvalues and their eigenvectors"
      (let [{:keys [var] :as m} (nth models 2000)
            C (lin/mat-add (lin/mat-mul (:W m) (lin/transpose (:W m))) (lin/mat-scale (lin/eye 4) var))
            U (mapv #(subvec % 2) vectors)
            noise (/ (+ (values 0) (values 1)) 2)]
        (is (< (abs (- var noise)) 1e-9))
        (is (< (max-abs (lin/mat-sub C (lin/mat-add (lin/mat-mul (lin/mat-mul U [[(values 2) 0] [0 (values 3)]]) (lin/transpose U))
                                                    (lin/mat-scale (lin/mat-sub (lin/eye 4) (lin/mat-mul U (lin/transpose U))) noise))))
               1e-6))))
    (testing "samples from the model have the data's mean and covariance"
      (let [m (nth models 2000) rng (random/rng 3)
            xs (vec (repeatedly 20000 #(gen/ppca-sample m rng)))]
        (is (< (max-abs (lin/sub (lin/mean xs) (:mean m))) 0.05))
        (is (< (max-abs (lin/mat-sub (gen/covariance xs) S)) 0.3))))
    (testing "a point off the plane is far less probable than one on it: an anomaly"
      (let [m (nth models 2000)]
        (is (< (+ (gen/ppca-log-density m [2 1 2 4]) 10) (gen/ppca-log-density m [1 2 3 4])))))
    (testing "the features z^ = E[z | x]: zero at the mean"
      (is (< (max-abs (gen/ppca-encode (nth models 100) (:mean (nth models 100)))) 1e-12)))))

(deftest linear-autoencoder-test
  (testing "trained on the squared error, W spans the m principal components, as in PCA: W^T W the projection onto the covariance's top two eigenvectors, W W^T = I, and the error the sum of the other eigenvalues"
    (let [{:keys [values vectors]} (eigen/jacobi (gen/covariance planar))
          m (nth (gen/linear-autoencoder-seq planar 2 (random/rng 3) :lr 0.02) 1000)
          U (mapv #(subvec % 2) vectors)]
      (is (< (max-abs (lin/mat-sub (lin/mat-mul (lin/transpose (:W m)) (:W m)) (lin/mat-mul U (lin/transpose U)))) 1e-9))
      (is (< (max-abs (lin/mat-sub (lin/mat-mul (:W m) (lin/transpose (:W m))) (lin/eye 2))) 1e-9))
      (is (< (abs (- (:loss m) (+ (values 0) (values 1)))) 1e-9))
      (is (< (abs (- (gen/reconstruction-error #(gen/linear-decode m (gen/linear-encode m %)) planar) (+ (values 0) (values 1)))) 1e-9)))))

(deftest deep-autoencoder-test
  (testing "data on a curve: one number per point is enough for a deep autoencoder, not for a linear one"
    (let [rng (random/rng 4)
          arc (vec (repeatedly 300 (fn [] (let [t (* Math/PI (rng))]
                                            [(+ (Math/cos t) (random/gaussian rng 0 0.03)) (+ (Math/sin t) (random/gaussian rng 0 0.03))]))))
          linear (nth (gen/linear-autoencoder-seq arc 1 (random/rng 1)) 500)
          deep (nth (gen/autoencoder-seq (gen/autoencoder 2 16 1 (random/rng 2)) arc {:alpha 0.02 :rng (random/rng 3)}) 150)
          linear-error (gen/reconstruction-error #(gen/linear-decode linear (gen/linear-encode linear %)) arc)
          deep-error (gen/reconstruction-error #(gen/decode deep (gen/encode deep %)) arc)]
      (is (= 1 (count (gen/encode deep (first arc)))))
      (is (< (* 20 deep-error) linear-error)))))

(deftest variational-autoencoder-test
  (testing "with a linear decoder the model is a linear Gaussian, its log P(x) exact: the variational lower bound is below it, and rises toward it in training"
    (let [rng (random/rng 1)
          X (vec (repeatedly 300 (fn [] (let [a (random/gaussian rng 0 2)]
                                          [(+ 1 a (random/gaussian rng 0 0.3)) (+ 2 (* 0.5 a) (random/gaussian rng 0 0.3)) (random/gaussian rng 0 0.3)]))))
          models (gen/vae-seq (gen/vae 3 8 1 0 0.09 (random/rng 2)) X (random/rng 3))
          r (random/rng 4)
          bound (fn [m] (mean (map #(gen/elbo m % r 400) (take 30 X))))
          exact (fn [m] (mean (map #(gen/vae-log-likelihood m %) (take 30 X))))]
      (doseq [n [0 5 60] :let [m (nth models n)]]
        (is (< (bound m) (exact m)) n))
      (is (< (bound (nth models 5)) (bound (nth models 60))))
      (is (< (- (exact (nth models 60)) (bound (nth models 60))) 1.0))))
  (testing "H(Q), the bound's entropy term, is the Gaussian's: what -log Q(z) averages to over samples of z"
    (let [m (gen/vae 2 4 2 4 0.1 (random/rng 5))
          {mu :mean var :var} (gen/vae-encode m [0.5 -1.0])
          rng (random/rng 6)
          neg-log-q (fn [] (let [z (mapv #(+ %1 (* (Math/sqrt %2) (random/gaussian rng 0 1))) mu var)]
                             (- (reduce + (map (fn [zi mi vi] (- (- (/ (* (- zi mi) (- zi mi)) (* 2 vi))) (* 0.5 (Math/log (* 2 Math/PI vi))))) z mu var)))))]
      (is (< (abs (- (mean (repeatedly 40000 neg-log-q))
                     (reduce + (map #(* 0.5 (+ 1 (Math/log (* 2 Math/PI %)))) var))))
             0.02))))
  (testing "a deep decoder learns a ring: samples near the circle, where a Gaussian -- probabilistic PCA -- scatters them across the disk"
    (let [m (nth (gen/vae-seq (gen/vae 2 32 2 32 0.01 (random/rng 5)) ring (random/rng 6)) 150)
          ppca (nth (gen/ppca-em-seq ring 2 (random/rng 7)) 100)
          rng (random/rng 8)
          radii (fn [sample] (map lin/length (repeatedly 400 sample)))]
      (is (< (sd (radii #(gen/vae-sample m rng))) 0.6))
      (is (< (abs (- 2.0 (mean (radii #(gen/vae-sample m rng))))) 0.4))
      (is (> (sd (radii #(gen/ppca-sample ppca rng))) 0.8)))))

(deftest autoregressive-test
  (testing "the Yule-Walker equations recover an AR(2) model from 20,000 of its values"
    (let [xs (vec (take 20000 (gen/ar-generate-seq (partial gen/ar-predict {:mean 10.0 :phi [0.6 -0.3]}) [10.0 10.0] 1.0 (random/rng 5))))
          {:keys [phi var] m :mean} (gen/yule-walker xs 2)]
      (is (< (abs (- m 10.0)) 0.05))
      (is (< (abs (- (phi 0) 0.6)) 0.02))
      (is (< (abs (- (phi 1) -0.3)) 0.02))
      (is (< (abs (- var 1.0)) 0.05))))
  (testing "the logistic map x' = 3.9 x (1 - x): no linear AR model predicts it, a deep one nearly exactly"
    (let [xs (vec (take 600 (iterate #(* 3.9 % (- 1 %)) 0.3)))
          [X Y] (gen/windows xs 2)
          net (nth (nn/train-seq [(nn/dense 2 16 :tanh (random/rng 6)) (nn/dense 16 1 :identity (random/rng 7))] X Y {:alpha 0.05 :rng (random/rng 8)}) 200)
          rms (fn [predict] (Math/sqrt (mean (map #(let [d (- (predict %1) (first %2))] (* d d)) X Y))))]
      (is (> (rms (partial gen/ar-predict (gen/yule-walker xs 2))) 0.2))
      (is (< (rms (partial gen/deep-ar-predict net)) 0.02))
      (testing "and from the first two values generates the next few -- only a few, the map being chaotic"
        (is (< (reduce max (map #(abs (- %1 %2)) (take 4 (gen/ar-generate-seq (partial gen/deep-ar-predict net) (subvec xs 0 2))) xs))
               0.1))))))

(deftest gan-test
  (testing "a generator learns, against a discriminator, to turn N(0, 1) noise into the data's N(4, 0.5^2) -- its spread short of the data's, as GANs' tend to be -- until the discriminator can do no better than guess"
    (let [rng (random/rng 1)
          data (vec (repeatedly 1000 #(vector (random/gaussian rng 4 0.5))))
          m (nth (gen/gan-seq (gen/gan 1 1 32 (random/rng 7)) data (random/rng 8) :lr 0.002) 1000)
          xs (map first (gen/gan-sample m (random/rng 9) 2000))]
      (is (< (abs (- (mean xs) 4.0)) 0.15))
      (is (< 0.15 (sd xs) 0.75))
      (is (< (abs (- 0.5 (mean (gen/discriminate m (take 500 data))))) 0.1))
      (is (< (abs (- 0.5 (mean (gen/discriminate m (gen/gan-sample m (random/rng 10) 500))))) 0.1)))))
