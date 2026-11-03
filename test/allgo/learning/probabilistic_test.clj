(ns allgo.learning.probabilistic-test
  (:require [allgo.learning.probabilistic :as pl]
            [allgo.probability.bayes :as bayes]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Learning probabilistic models: Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 20, its worked
;; examples -- the candy bags, the beta family, Bayesian linear regression,
;; and EM on two bags mixed together.

(defn- close? [a b tol] (< (abs (- (double a) (double b))) tol))

;; The five kinds of bag, h1 all cherry to h5 all lime, and the
;; manufacturer's prior over them.
(def ^:private bags
  (mapv (fn [cherry] (fn [d] (if (= d :cherry) cherry (- 1.0 cherry)))) [1.0 0.75 0.5 0.25 0.0]))
(def ^:private bag-prior [0.1 0.2 0.4 0.2 0.1])

(deftest bayesian-learning
  (let [posts (vec (pl/posterior-seq bags bag-prior (repeat 10 :lime)))]
    (testing "Figure 20.1(a): h3 most likely at first and after one lime, h4 after two, h5 after three or more"
      (is (= [2 2 3 4 4 4 4 4 4 4 4] (mapv pl/map-hypothesis posts))))
    (testing "and after ten limes in a row, fairly certain of h5: 0.1 against 0.2 (3/4)^10 for h4, and less"
      (is (close? 0.8956 ((posts 10) 4) 1e-4)))
    (testing "Figure 20.1(b): the Bayesian prediction that the next is lime rises toward 1"
      (let [p (mapv #(pl/predict bags % :lime) posts)]
        (is (close? 0.5 (p 0) 1e-12))
        (is (apply < p))))
    (testing "after three limes the MAP learner says lime for certain; the Bayesian prediction is about 0.8"
      (is (= 1.0 ((nth bags (pl/map-hypothesis (posts 3))) :lime)))
      (is (close? 0.8 (pl/predict bags (posts 3) :lime) 0.01)))
    (testing "maximum likelihood ignores the prior: one lime makes h5 the likeliest"
      (is (= 4 (pl/ml-hypothesis bags [:lime]))))))

(deftest maximum-likelihood
  (testing "the candy network of Figure 20.2(b): each parameter its observed fraction"
    (let [structure {:variables [:Flavor :Wrapper] :parents {:Flavor [] :Wrapper [:Flavor]}
                     :domains {:Flavor [:cherry :lime] :Wrapper [:red :green]}}
          data (concat (repeat 6 {:Flavor :cherry :Wrapper :red}) (repeat 2 {:Flavor :cherry :Wrapper :green})
                       (repeat 1 {:Flavor :lime :Wrapper :red}) (repeat 3 {:Flavor :lime :Wrapper :green}))
          net (pl/learn-parameters structure data)]
      (is (close? (/ 8.0 12) (get-in net [:cpt :Flavor [] :cherry]) 1e-12))
      (is (close? (/ 6.0 8) (get-in net [:cpt :Wrapper [:cherry] :red]) 1e-12))
      (is (close? (/ 1.0 4) (get-in net [:cpt :Wrapper [:lime] :red]) 1e-12))
      (testing "and with too little data, a zero -- which Laplace's added count avoids"
        (let [one (pl/learn-parameters structure [{:Flavor :cherry :Wrapper :red}])
              laplace (pl/learn-parameters structure [{:Flavor :cherry :Wrapper :red}] :pseudo 1)]
          (is (= 0.0 (get-in one [:cpt :Flavor [] :lime])))
          (is (close? (/ 1.0 3) (get-in laplace [:cpt :Flavor [] :lime]) 1e-12))))))
  (testing "equation 20.4: the Gaussian's mean the sample's, its deviation the root mean square about it"
    (is (= {:mean 2.0 :sd (math/sqrt (/ 2.0 3))} (pl/gaussian-ml [1.0 2.0 3.0]))))
  (testing "a linear Gaussian: least squares recovers the line, the residuals the noise"
    (let [rng (random/rng 1)
          xs (vec (repeatedly 2000 rng))
          ys (mapv #(+ (* 0.5 %) 0.2 (random/gaussian rng 0.0 0.05)) xs)
          {:keys [theta1 theta2 sigma]} (pl/linear-gaussian-ml xs ys)]
      (is (close? 0.5 theta1 0.01))
      (is (close? 0.2 theta2 0.01))
      (is (close? 0.05 sigma 0.003)))))

(deftest bayesian-parameters
  (testing "beta[a, b]: its mean a / (a + b); beta[1, 1] the uniform density"
    (is (close? 0.75 (pl/beta-mean [3 1]) 1e-12))
    (is (every? #(close? 1.0 (pl/beta-density 1 1 %) 1e-12) [0.1 0.5 0.9])))
  (testing "it integrates to 1"
    (doseq [[a b] [[2 2] [5 5] [3 1] [30 10]]]
      (is (close? 1.0 (* 0.001 (reduce + (map #(pl/beta-density a b (* 0.001 (+ % 0.5))) (range 1000)))) 1e-3))))
  (testing "closed under update: a cherry adds one to a, a lime one to b -- virtual counts"
    (is (= [[1 1] [2 1] [2 2] [3 2]] (vec (pl/beta-seq [1 1] [true false true])))))
  (testing "Figure 20.5(b): beta[3,1], beta[6,2], beta[30,10] narrowing round 0.75"
    (let [spread (fn [[a b]] (math/sqrt (/ (* a b) (* (+ a b) (+ a b) (+ a b 1.0)))))]
      (is (apply > (map spread [[3 1] [6 2] [30 10]])))
      (is (every? #(close? 0.75 (pl/beta-mean %) 1e-12) [[3 1] [6 2] [30 10]]))))
  (testing "the Dirichlet, the same for many values"
    (is (= {:a 0.6 :b 0.2 :c 0.2} (pl/dirichlet-predictive (pl/dirichlet-update {:a 1 :b 1 :c 1} [:a :a]))))))

(deftest bayesian-linear-regression
  (let [noise 0.2
        xs [-0.5 0.1 0.5 -1.8 1.8] ys [-0.4 0.1 0.2 -0.9 0.9]
        prior {:mean 0.0 :var 4.0}
        {:keys [mean var]} (pl/regression-posterior prior noise xs ys)]
    (testing "the closed form is the posterior P(theta | d) proportional to P(d | theta) P(theta), on a grid"
      (let [grid (map #(* 0.0005 %) (range -4000 4001))
            w (map (fn [t] (math/exp (- (* -0.5 (/ (* t t) (:var prior)))
                                        (* 0.5 (/ (reduce + (map #(let [d (- %2 (* t %1))] (* d d)) xs ys)) noise)))))
                   grid)
            z (reduce + w)
            m (/ (reduce + (map * grid w)) z)
            v (/ (reduce + (map #(* (- %1 m) (- %1 m) %2) grid w)) z)]
        (is (close? m mean 1e-6))
        (is (close? v var 1e-6))))
    (testing "a point at a time, each posterior the next prior, ends at the same posterior"
      (let [post (last (pl/regression-posterior-seq prior noise xs ys))]
        (is (close? mean (:mean post) 1e-12))
        (is (close? var (:var post) 1e-12))))
    (testing "an uninformative prior: the least-squares slope, its variance sigma^2 / sum x^2"
      (let [p (pl/regression-posterior {:mean 0.0 :var 1e12} noise xs ys)
            sxx (reduce + (map * xs xs))]
        (is (close? (/ (reduce + (map * xs ys)) sxx) (:mean p) 1e-9))
        (is (close? (/ noise sxx) (:var p) 1e-9))))
    (testing "data near the origin barely constrain the slope: the posterior variance about the prior's"
      (is (close? 4.0 (:var (pl/regression-posterior prior noise [0.001 -0.002] [0.0 0.0])) 1e-3)))
    (testing "Figure 20.7: sigma^2 = 0.2, and the slope's variance 0.3861 falls to 0.0286 with the two points added at about x = +-1.8"
      (is (close? 0.0286 (:var (pl/regression-posterior {:mean 0.0 :var 0.3861} 0.2 [-1.8 1.8] [-0.9 0.9])) 1e-4)))
    (testing "the prediction's variance the noise's plus sigma_N^2 x^2, its mean theta_N x"
      (is (= {:mean (* mean 2.0) :var (+ noise (* var 4.0))} (pl/regression-predictive {:mean mean :var var} noise 2.0)))
      (is (apply < (map #(:var (pl/regression-predictive {:mean mean :var var} noise %)) [0.0 1.0 2.0 4.0]))))))

(deftest structure-learning
  (let [domains (zipmap (:variables bayes/sprinkler) (repeat [true false]))
        skeleton (fn [parents] (set (for [[v ps] parents p ps] #{v p})))
        right? (fn [parents] (and (= (skeleton (:parents bayes/sprinkler)) (skeleton parents))
                                  (= #{:Sprinkler :Rain} (set (parents :WetGrass)))))
        samples (fn [seed] (let [rng (random/rng seed)] (vec (repeatedly 3000 #(bayes/prior-sample bayes/sprinkler rng)))))
        data (samples 1)
        climb (vec (pl/structure-search-seq domains data))]
    (testing "a move a step, the score rising every one"
      (is (apply < (map :score climb))))
    (testing "one climb from no links can stop at a local maximum, below the sprinkler network's own score"
      (is (not (right? (:parents (peek climb)))))
      (is (< (:score (peek climb)) (pl/structure-score domains (:parents bayes/sprinkler) data))))
    (testing "climbing again from ten random networks finds it: its links, and the two into WetGrass directed -- on each of five samples"
      (doseq [seed (range 1 6)]
        (is (right? (:parents (pl/learn-structure domains (samples seed) :restarts 10 :rng (random/rng (+ 100 seed))))) seed)))
    (testing "with no penalty, more links always fit better: the search adds links the network does not have"
      (is (> (count (skeleton (:parents (peek (vec (pl/structure-search-seq domains data :penalty 0.0))))))
             (count (skeleton (:parents bayes/sprinkler))))))))

(def ^:private blobs
  "Three Gaussian blobs in the unit square, equally likely: `[[center sd]
  ...]`."
  [[[0.25 0.3] 0.1] [[0.7 0.7] 0.12] [[0.75 0.25] 0.09]])

(defn- mixture-sample [rng]
  (let [[[mx my] s] (random/pick rng blobs)]
    [(random/gaussian rng mx s) (random/gaussian rng my s)]))

(defn- mixture-density [[x y]]
  (/ (reduce + (for [[[mx my] s] blobs]
                 (/ (math/exp (- (/ (+ (* (- x mx) (- x mx)) (* (- y my) (- y my))) (* 2 s s)))) (* 2 math/PI s s))))
     (count blobs)))

(deftest density-estimation
  (doseq [seed [2 3 4]]
    (let [rng (random/rng seed)
          train (vec (repeatedly 128 #(mixture-sample rng)))
          probe (vec (repeatedly 300 #(mixture-sample rng)))
          ;; how far off an estimate is: the mean absolute log ratio to the
          ;; true density, at points drawn from it
          error (fn [density] (/ (reduce + (map #(abs (- (math/log (density %)) (math/log (mixture-density %)))) probe)) (count probe)))]
      (testing "Figure 20.10's kernels on 128 points: too narrow is spiky, too wide too smooth -- 0.07 nearest the truth"
        (let [[narrow mid wide] (map (fn [w] (error #(pl/kernel-density train w %))) [0.02 0.07 0.2])]
          (is (< mid narrow) seed)
          (is (< mid wide) seed)))
      (testing "Figure 20.9's nearest neighbors: k = 10 nearer the truth than 3 or 40"
        (let [[few mid many] (map (fn [k] (error #(pl/knn-density train k %))) [3 10 40])]
          (is (< mid few) seed)
          (is (< mid many) seed)))
      (testing "the kernel estimate is a density: it integrates to about 1"
        (is (close? 1.0 (* 0.0001 (reduce + (for [i (range -50 150) j (range -50 150)]
                                              (pl/kernel-density train 0.07 [(* 0.01 (+ i 0.5)) (* 0.01 (+ j 0.5))]))))
                    0.02)
            seed)))))

(deftest em-candy
  ;; Section 20.3.2: two bags mixed, each candy's Flavor (true cherry),
  ;; Wrapper (true red) and Hole; 1000 drawn from bags with theta 0.5,
  ;; 0.8 for bag 1's three features and 0.3 for bag 2's (equation 20.9),
  ;; EM started from equation 20.10's.
  (let [data [[[true true true] 273] [[true true false] 93] [[true false true] 104] [[true false false] 90]
              [[false true true] 79] [[false true false] 100] [[false false true] 94] [[false false false] 167]]
        start {:prior [0.6 0.4] :features [[0.6 0.6 0.6] [0.4 0.4 0.4]]}
        states (vec (take 11 (pl/naive-bayes-em-seq start data)))
        {:keys [prior features]} (:params (states 1))
        truth (pl/mixture-log-likelihood {:prior [0.5 0.5] :features [[0.8 0.8 0.8] [0.3 0.3 0.3]]} data)]
    (testing "equation 20.11: the parameters after one iteration"
      (is (close? 0.6124 (prior 0) 5e-5))
      (is (every? true? (map #(close? %1 %2 5e-5) [0.6684 0.6483 0.6558] (features 0))))
      (is (every? true? (map #(close? %1 %2 5e-5) [0.3887 0.3817 0.3827] (features 1)))))
    (testing "the log likelihood from about -2044 to about -2021 in the first iteration"
      (is (close? -2044 (:log-likelihood (states 0)) 1.0))
      (is (close? -2021 (:log-likelihood (states 1)) 1.0)))
    (testing "never falling; by the tenth iteration a better fit than the model that made the data, L = -1982.214"
      (is (apply <= (map :log-likelihood states)))
      (is (close? -1982.214 truth 1e-3))
      (is (> (:log-likelihood (states 10)) truth)))))
