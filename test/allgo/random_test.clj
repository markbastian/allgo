(ns allgo.random-test
  (:require [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest generator-test
  (testing "the stream is in range and never falls into the fixed point"
    (loop [i 0 s (random/cell-seed 0 1 2 3)]
      (when (< i 10000)
        (is (pos? s))
        (is (< s 2147483647.0))
        (is (<= 0.0 (random/unit s)))
        (is (< (random/unit s) 1.0))
        (recur (inc i) (random/advance s)))))

  (testing "a cell's stream is a function of the cell, and of nothing else"
    (is (== (random/cell-seed 3 10 -4 7) (random/cell-seed 3 10 -4 7)))
    (is (not= (random/cell-seed 3 10 -4 7) (random/cell-seed 3 10 -4 8)))
    (is (not= (random/cell-seed 3 10 -4 7) (random/cell-seed 4 10 -4 7)))))

(deftest rng-test
  (testing "a seed gives the same numbers every time"
    (is (= (repeatedly 100 (random/rng 42)) (repeatedly 100 (random/rng 42)))))

  (testing "and different seeds, and different keys, give different ones"
    (is (not= (repeatedly 10 (random/rng 42)) (repeatedly 10 (random/rng 43))))
    (is (not= (repeatedly 10 (random/rng 42 1 0 0)) (repeatedly 10 (random/rng 42 2 0 0)))))

  (testing "uniform on [0, 1)"
    (let [xs (repeatedly 20000 (random/rng 7))]
      (is (every? #(and (<= 0.0 %) (< % 1.0)) xs))
      (is (< 0.49 (/ (reduce + xs) 20000) 0.51)))))

(deftest draws-test
  (let [rng (random/rng 5)]
    (testing "below stays in range and reaches both ends"
      (let [ks (repeatedly 5000 #(random/below rng 6))]
        (is (= (set (range 6)) (set ks)))))

    (testing "uniform stays in its interval"
      (is (every? #(and (<= -2.0 %) (< % 3.0)) (repeatedly 1000 #(random/uniform rng -2.0 3.0)))))

    (testing "pick only returns members"
      (is (every? #{:a :b :c} (repeatedly 200 #(random/pick rng [:a :b :c])))))

    (testing "shuffle is a permutation, and a repeatable one"
      (let [xs (range 50)
            s1 (random/shuffle (random/rng 9) xs)]
        (is (vector? s1))
        (is (= (sort s1) xs))
        (is (not= s1 (vec xs)))
        (is (= s1 (random/shuffle (random/rng 9) xs)))
        (is (= [] (random/shuffle rng [])))))

    (testing "sample keeps about the fraction asked for"
      (let [n (count (random/sample rng 0.3 (range 10000)))]
        (is (< 2800 n 3200))))

    (testing "gaussian has the mean and spread asked for"
      (let [xs (repeatedly 20000 #(random/gaussian rng 24.0 8.0))
            mean (/ (reduce + xs) (count xs))
            sd   (math/sqrt (/ (reduce + (map #(let [d (- % mean)] (* d d)) xs))
                               (dec (count xs))))]
        (is (< 23.5 mean 24.5))
        (is (< 7.6 sd 8.4))))))
