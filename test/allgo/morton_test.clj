(ns allgo.morton-test
  (:require [allgo.spatial.morton :as m]
            [clojure.test :refer [deftest is testing]]))

(deftest interleaving-test
  (testing "a single axis lands on every third bit"
    (is (= 1 (m/encode-3 1 0 0)))
    (is (= 2 (m/encode-3 0 1 0)))
    (is (= 4 (m/encode-3 0 0 1)))
    (is (= 7 (m/encode-3 1 1 1)))
    (is (= 8 (m/encode-3 2 0 0)))
    (is (= 64 (m/encode-3 4 0 0))))

  (testing "a single axis lands on every other bit in 2D"
    (is (= 1 (m/encode-2 1 0)))
    (is (= 2 (m/encode-2 0 1)))
    (is (= 3 (m/encode-2 1 1))))

  (testing "the whole range fits in 30 bits and stays positive"
    (is (= (dec (bit-shift-left 1 30)) (m/encode-3 m/max-3d m/max-3d m/max-3d)))
    (is (= (dec (bit-shift-left 1 30)) (m/encode-2 m/max-2d m/max-2d)))
    (is (pos? (m/encode-3 m/max-3d m/max-3d m/max-3d))))

  (testing "codes come apart again"
    (doseq [xyz [[0 0 0] [1 2 3] [1023 0 512] [7 7 7] [1023 1023 1023]]]
      (is (= xyz (m/decode-3 (apply m/encode-3 xyz)))))
    (doseq [xy [[0 0] [1 2] [32767 0] [32767 32767]]]
      (is (= xy (m/decode-2 (apply m/encode-2 xy))))))

  (testing "round trip holds across the range"
    (let [rng (java.util.Random. 42)
          r (fn [n] (.nextInt rng n))]
      (is (every? true?
                  (repeatedly 3000
                              #(let [v [(r 1024) (r 1024) (r 1024)]]
                                 (= v (m/decode-3 (apply m/encode-3 v)))))))))

  (testing "expanding and compacting are inverse"
    (doseq [v [0 1 2 511 1023]]
      (is (= v (m/compact-3 (m/expand-3 v))))
      (is (= v (m/compact-2 (m/expand-2 v)))))))

(deftest monotone-test
  (testing "holding two axes fixed, the code rises with the third"
    ;; Not true across axes -- that is the whole nature of the curve --
    ;; but along one it has to be, or the code is not an encoding at all.
    (doseq [axis (range 3)]
      (let [code (fn [v] (apply m/encode-3 (assoc [5 5 5] axis v)))]
        (is (apply < (map code (range 0 64))))))))

(deftest quantize-test
  (testing "a value maps onto the grid it is given"
    (is (= 0 (m/quantize 0.0 0.0 1.0 4)))
    (is (= 15 (m/quantize 1.0 0.0 1.0 4)))
    (is (= 8 (m/quantize 0.5 0.0 1.0 4))))

  (testing "out of range clamps rather than wrapping"
    ;; Wrapping would put a point just outside the box at the opposite
    ;; corner of the space, which is as wrong as a code can be.
    (is (= 0 (m/quantize -100.0 0.0 1.0 10)))
    (is (= 1023 (m/quantize 100.0 0.0 1.0 10))))

  (testing "a degenerate range does not divide by zero"
    (is (= 0 (m/quantize 5.0 5.0 5.0 10)))))

(deftest bounds-test
  (testing "bounds cover every point"
    (let [pts (double-array [0.0 1.0 2.0, -1.0 5.0 0.5, 3.0 0.0 -2.0])
          [lo hi] (m/bounds pts 3)]
      (is (= [-1.0 0.0 -2.0] lo))
      (is (= [3.0 5.0 2.0] hi)))))

(deftest order-test
  (let [n 2000
        rng (java.util.Random. 7)
        pts (double-array (repeatedly (* 3 n) #(* 20.0 (- (.nextDouble rng) 0.5))))
        codes (m/codes pts n)
        order (m/order codes n)]

    (testing "the radix sort agrees with a comparison sort"
      (is (= (vec (sort-by (fn [i] [(aget ^ints codes i) i]) (range n)))
             (mapv #(aget ^ints order %) (range n)))))

    (testing "codes come out non-decreasing"
      (is (apply <= (map #(aget ^ints codes (aget ^ints order %)) (range n)))))

    (testing "every object appears exactly once"
      (is (= (set (range n)) (set (map #(aget ^ints order %) (range n))))))

    (testing "the sort is stable, so a rerun gives the same answer"
      (is (= (vec (m/order codes n)) (vec (m/order codes n))))))

  (testing "sorting puts nearby points near each other"
    ;; The only reason to take a Morton code at all.
    (let [n 3000
          rng (java.util.Random. 11)
          pts (double-array (repeatedly (* 3 n) #(* 10.0 (.nextDouble rng))))
          order (m/sorted-points pts n)
          dist (fn [i j]
                 (let [a (* 3 i) b (* 3 j)]
                   (Math/sqrt (reduce + (map (fn [k]
                                               (let [d (- (aget ^doubles pts (+ a k))
                                                          (aget ^doubles pts (+ b k)))]
                                                 (* d d)))
                                             (range 3))))))
          consecutive (/ (reduce + (map #(dist (aget ^ints order %) (aget ^ints order (inc %)))
                                        (range (dec n))))
                         (dec n))
          arbitrary (/ (reduce + (map (fn [i] (dist i (mod (* 7919 i) n))) (range n))) n)]
      (is (< consecutive (* 0.4 arbitrary))
          (str "consecutive " consecutive " should be well under arbitrary " arbitrary))))

  (testing "an empty set sorts to nothing"
    (is (zero? (alength ^ints (m/order (m/codes (double-array 0) 0) 0))))))
