(ns allgo.spatial.kdtree-test
  (:require [allgo.random :as random]
            [allgo.spatial.kdtree :as kd]
            [clojure.test :refer [deftest is testing]]))

;; k-d trees, after LaValle, *Planning Algorithms*, Cambridge University
;; Press, 2006, section 5.5.2, checked against measuring every point.

(defn- sq-dist [a b] (reduce + (map #(let [d (- %1 %2)] (* d d)) a b)))

(deftest kdtree-test
  (let [rng (random/rng 1)
        points (vec (repeatedly 2000 (fn [] [(rng) (rng) (rng)])))
        tree (kd/build points)
        queries (repeatedly 100 (fn [] [(rng) (rng) (rng)]))]
    (testing "the nearest point is the one measuring them all finds"
      (is (every? (fn [q] (= (apply min-key #(sq-dist % q) points) (first (kd/nearest tree q)))) queries)))
    (testing "visiting only a small part of the tree: under a tenth of the points on average"
      (is (< (/ (reduce + (map #(:visited (meta (kd/nearest tree %))) queries)) 100.0) 200)))
    (testing "the k nearest, in order"
      (is (every? (fn [q] (= (take 5 (sort-by #(sq-dist % q) points)) (kd/k-nearest tree q 5))) (take 20 queries))))
    (testing "every point within a radius"
      (is (every? (fn [q] (= (set (filter #(<= (sq-dist % q) 0.01) points)) (set (kd/within tree q 0.1)))) (take 20 queries))))))
