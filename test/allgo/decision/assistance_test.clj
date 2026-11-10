(ns allgo.decision.assistance-test
  (:require [allgo.decision.assistance :as assist]
            [clojure.test :refer [deftest is testing]]))

;; Assistance games, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition: the off-switch game of section 16.7.2 and
;; the paperclip game of section 18.2.5, with the sections' own numbers.

(defn- close? [a b] (< (abs (- a b)) 1e-3))

(deftest off-switch-test
  (testing "the hotel worth anywhere from -40 to +60 to Harriet: booking now is worth +10, deferring to her +18 -- 0.4 x 0 + 0.6 x 30"
    (let [{:keys [act defer best]} (assist/off-switch (assist/uniform-prior -40 60 100001))]
      (is (close? 10.0 act))
      (is (close? 18.0 defer))
      (is (= :defer best))))
  (testing "from -60 to +40 he would switch himself off rather than act -- but deferring is better still, 0.4 x 20"
    (let [{:keys [act defer best]} (assist/off-switch (assist/uniform-prior -60 40 100001))]
      (is (close? -10.0 act))
      (is (close? 8.0 defer))
      (is (= :defer best))))
  (testing "deferring is never worse than acting, and equal only when he is sure she wants it"
    (let [{:keys [act defer]} (assist/off-switch (assist/uniform-prior 5 60 1001))]
      (is (close? act defer))))
  (testing "a cost for her time takes it off deferring's value"
    (is (close? 13.0 (:defer (assist/off-switch (assist/uniform-prior -40 60 100001) 5.0 0.0)))))
  (testing "the more often she decides wrongly, the less deferring is worth: at 0.4 he should act"
    (let [at (fn [e] (assist/off-switch (assist/uniform-prior -40 60 100001) 0.0 e))]
      (is (= :defer (:best (at 0.3))))
      (is (= :act (:best (at 0.4)))))))

(deftest paperclip-test
  (let [game assist/paperclip-game
        prior (assist/uniform-prior 0.0 1.0 100001)
        [greedy second third] (assist/best-response-seq game prior)
        eq (assist/equilibrium game prior)]
    (testing "on her own, a Harriet who values a paperclip at 0.45 makes two staples, worth $1.10"
      (is (= [0 2] ((assist/greedy-human game prior) 0.45)))
      (is (close? 1.10 (assist/payoff 0.45 [0 2]))))
    (testing "seeing greedy Harriet make two paperclips, Robbie believes theta uniform on (0.5, 1] and makes 90 paperclips, worth $67.50; one of each, 50 of each; two staples, 90 staples"
      (is (= {[2 0] [90 0] [1 1] [50 50] [0 2] [0 90]} (:robot greedy)))
      (is (close? 67.5 (assist/payoff 0.75 [90 0]))))
    (testing "knowing that, Harriet makes one of each not only when indifferent but anywhere from about 0.446 to 0.554"
      (let [[lo hi] (assist/signal-range (:human second) [1 1])]
        (is (close? lo (assist/threshold [[1 1] [50 50]] [[0 2] [0 90]])))
        (is (close? hi (assist/threshold [[1 1] [50 50]] [[2 0] [90 0]])))
        (is (close? 0.446 lo))
        (is (close? 0.554 hi))))
    (testing "and Robbie's reply to that is unchanged: the equilibrium, after one round"
      (is (= second third eq)))
    (testing "the Harriet at 0.45 now makes one of each, and gets 50 of each: worth $51, more than the $50.60 of two staples and 90 more"
      (is (= [1 1] ((:human eq) 0.45)))
      (is (close? 51.0 (+ (assist/payoff 0.45 [1 1]) (assist/payoff 0.45 [50 50]))))
      (is (close? 50.6 (+ (assist/payoff 0.45 [0 2]) (assist/payoff 0.45 [0 90])))))))
