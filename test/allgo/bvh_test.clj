(ns allgo.bvh-test
  (:require [allgo.spatial.bvh :as bvh]
            [allgo.spatial.sweep :as sweep]
            [clojure.test :refer [deftest is testing]]))

(defn- scene [seed n {:keys [mixed?] :or {mixed? false}}]
  (let [rng (java.util.Random. seed)
        pos (double-array (repeatedly (* 3 n) #(* 10.0 (.nextDouble rng))))
        half (double-array (repeatedly (* 3 n)
                                       #(if mixed?
                                          (+ 0.05 (* 1.2 (Math/pow (.nextDouble rng) 6)))
                                          0.35)))]
    (into [pos] (sweep/boxes-of pos half n))))

(defn- brute [^doubles mins ^doubles maxs n]
  (vec (for [i (range n) j (range (inc i) n)
             :when (every? (fn [k]
                             (and (<= (aget mins (+ (* 3 i) k)) (aget maxs (+ (* 3 j) k)))
                                  (>= (aget maxs (+ (* 3 i) k)) (aget mins (+ (* 3 j) k)))))
                           (range 3))]
         [i j])))

(deftest structure-test
  (testing "an empty tree is a legitimate thing to build and to query"
    ;; The reference halves an empty range and recurses until the stack
    ;; goes.
    (let [t (bvh/build (double-array 0) (double-array 0) 0)]
      (is (zero? (bvh/node-count t)))
      (is (zero? (bvh/depth t)))
      (is (= [] (bvh/query t [0 0 0] [1 1 1])))
      (is (= [] (bvh/overlapping-pairs t (double-array 0) (double-array 0) 0)))))

  (testing "one box is one leaf"
    (let [t (bvh/build (double-array [0 0 0]) (double-array [1 1 1]) 1)]
      (is (= 1 (bvh/node-count t)))
      (is (= 1 (bvh/depth t)))
      (is (bvh/leaf? t (:root t)))
      (is (= [0] (bvh/query-point t [0.5 0.5 0.5])))
      (is (= [] (bvh/query-point t [9 9 9])))))

  (testing "a halving split gives 2n-1 nodes and a tree as shallow as it can be"
    (doseq [n [2 3 7 8 64 500]]
      (let [[_ mins maxs] (scene 3 n {})
            t (bvh/build mins maxs n)]
        (is (= (dec (* 2 n)) (bvh/node-count t)))
        (is (= (inc (long (Math/ceil (/ (Math/log n) (Math/log 2)))))
               (bvh/depth t))
            (str n " boxes"))))))

(deftest query-test
  (testing "a region query finds exactly the boxes overlapping it"
    (let [n 400
          [_ mins maxs] (scene 5 n {})
          t (bvh/build mins maxs n)
          lo [2.0 2.0 2.0] hi [5.0 5.0 5.0]
          want (set (for [i (range n)
                          :when (every? (fn [k]
                                          (and (<= (aget ^doubles mins (+ (* 3 i) k)) (nth hi k))
                                               (>= (aget ^doubles maxs (+ (* 3 i) k)) (nth lo k))))
                                        (range 3))]
                      i))]
      (is (seq want))
      (is (= want (set (bvh/query t lo hi))))))

  (testing "a query outside everything finds nothing"
    (let [[_ mins maxs] (scene 5 100 {})
          t (bvh/build mins maxs 100)]
      (is (= [] (bvh/query t [100 100 100] [101 101 101]))))))

(deftest pairs-test
  (testing "the tree finds exactly the pairs brute force does"
    (doseq [n [2 5 50 400]]
      (let [[_ mins maxs] (scene 7 n {})
            t (bvh/build mins maxs n)]
        (is (= (sort (brute mins maxs n))
               (sort (bvh/overlapping-pairs t mins maxs n)))
            (str n " boxes")))))

  (testing "including when the boxes are wildly different sizes"
    (let [n 400
          [_ mins maxs] (scene 11 n {:mixed? true})
          t (bvh/build mins maxs n)]
      (is (= (sort (brute mins maxs n))
             (sort (bvh/overlapping-pairs t mins maxs n))))))

  (testing "each pair comes back once, not once from each end"
    ;; The reference queries every box against the whole tree and keeps
    ;; both directions, so its pair count is twice the brute force count
    ;; it is displayed beside.
    (let [n 300
          [_ mins maxs] (scene 13 n {})
          t (bvh/build mins maxs n)
          pairs (bvh/overlapping-pairs t mins maxs n)]
      (is (seq pairs))
      (is (every? (fn [[i j]] (< i j)) pairs))
      (is (= (count pairs) (count (distinct pairs))))
      (is (= (count (brute mins maxs n)) (count pairs)))))

  (testing "and it agrees with the other two broad phases exactly"
    ;; Three different structures, one answer -- which is what makes them
    ;; interchangeable.
    (let [n 400
          [_ mins maxs] (scene 17 n {:mixed? true})
          t (bvh/build mins maxs n)
          s (sweep/sweep n)]
      (is (= (sort (bvh/overlapping-pairs t mins maxs n))
             (sort (sweep/overlapping-pairs s mins maxs n)))))))

(deftest adjacency-test
  (testing "adjacency holds the same pairs as the pair list"
    (let [n 300
          [_ mins maxs] (scene 19 n {})
          t (bvh/build mins maxs n)
          pairs (bvh/overlapping-pairs t mins maxs n)
          {:keys [^ints starts ^ints ids] :as adj} (bvh/adjacency t mins maxs n)]
      (is (= (count pairs) (:pairs adj) (alength ids)))
      (is (= (set pairs)
             (set (for [j (range n)
                        k (range (aget starts j) (aget starts (inc j)))]
                    [(aget ids k) j])))))))
