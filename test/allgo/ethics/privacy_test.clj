(ns allgo.ethics.privacy-test
  (:require [allgo.ethics.privacy :as privacy]
            [allgo.learning.models :as models]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Privacy, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, section 27.3.2.

(deftest difference-attack-test
  (testing "the section's example: the average salary of 13 employees aged 30-41 and of the 12 aged 30-40 give away the one 41-year-old's"
    (is (= 80779 (privacy/difference-attack [81199 13] [81234 12]))))
  (testing "a query over too few records gets no answer"
    (is (nil? (privacy/aggregate [{:a 1} {:a 2}] (constantly true) :a 5)))
    (is (= [1.5 2] (privacy/aggregate [{:a 1} {:a 2}] (constantly true) :a 2)))))

(def ^:private hierarchies
  {:age [#(* 10 (quot % 10)) #(* 20 (quot % 20)) (constantly "*")]
   :zip [#(* 10 (quot % 10)) (constantly "*")]
   :sex [(constantly "*")]})

(deftest k-anonymity-test
  (let [rng (random/rng 1)
        records (vec (repeatedly 300 (fn [] {:age (+ 18 (random/below rng 70)) :sex (random/pick rng [:f :m])
                                             :zip (+ 94700 (random/below rng 30)) :diagnosis (random/pick rng [:flu :none])})))
        quasi [:age :sex :zip]]
    (testing "age, sex and zip code with names stripped: most records are the only ones of their kind"
      (is (= 1 (privacy/anonymity records quasi)))
      (is (> (count (filter #(= 1 (count %)) (vals (privacy/classes records quasi)))) 250)))
    (testing "generalized until 2-anonymous: a few records still stand out after ages to decades and zips to tens, so age goes entirely"
      (let [{:keys [records levels]} (privacy/anonymize records hierarchies 2)]
        (is (>= (privacy/anonymity records quasi) 2))
        (is (= 3 (levels :age)))))
    (testing "withholding up to three records instead, the decades are kept"
      (let [{:keys [records levels]} (privacy/anonymize records hierarchies 2 3)]
        (is (>= (privacy/anonymity records quasi) 2))
        (is (= {:age 1 :zip 1 :sex 0} levels))
        (is (<= 297 (count records) 299))))))

(deftest differential-privacy-test
  (testing "Laplace noise of scale 1 / epsilon on a count: whether one more person is in it changes the log density of any answer by at most epsilon"
    (is (every? #(<= (abs (- (privacy/laplace-log-density 100 1.0 0.5 %) (privacy/laplace-log-density 101 1.0 0.5 %))) (+ 0.5 1e-12))
                (range 80 120 0.25))))
  (testing "the noisy answers are right on average"
    (let [rng (random/rng 2)
          answers (repeatedly 20000 #(privacy/private-count (range 100) even? 0.5 rng))]
      (is (< (abs (- 50.0 (/ (reduce + answers) 20000))) 0.1)))))

(deftest federated-test
  (let [rng (random/rng 3)
        ;; five users, each seeing only part of the space: the truth is x + y > 0
        user (fn [cx cy] (let [X (vec (repeatedly 40 (fn [] [(random/gaussian rng cx 1.0) (random/gaussian rng cy 1.0)])))]
                           [X (mapv (fn [[x y]] (if (pos? (+ x y)) 1 0)) X)]))
        users [(user -2.0 1.5) (user 2.0 -1.5) (user 0.0 0.0) (user 1.5 1.5) (user -1.5 -1.5)]
        test (user 0.0 0.0)
        accuracy (fn [w [X ys]] (/ (count (filter true? (map (fn [x y] (= y (if (> (models/predict-logistic w x) 0.5) 1 0))) X ys))) (double (count ys))))
        shared (nth (privacy/federated-seq users 2 0.5 5) 30)
        alone (privacy/local-update [0.0 0.0 0.0] (first users) 0.5 150)]
    (testing "averaging the users' locally trained weights, no data leaving any of them, learns the whole rule"
      (is (> (accuracy shared test) 0.95)))
    (testing "where one user alone, seeing only their corner, does worse"
      (is (< (accuracy alone test) (accuracy shared test)))))
  (testing "secure aggregation: each masked upload is far from what it hides, and the masks cancel in the sum"
    (let [vs [[1.0 2.0] [3.0 -1.0] [0.5 0.5]]
          up (privacy/masked-uploads vs 7 100.0)]
      (is (every? (fn [[u v]] (> (reduce max (map #(abs (- %1 %2)) u v)) 10.0)) (map vector up vs)))
      (is (every? #(< (abs %) 1e-9) (map - (apply mapv + up) (apply mapv + vs)))))))
