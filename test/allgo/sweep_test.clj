(ns allgo.sweep-test
  (:require [allgo.spatial.hash :as hash]
            [allgo.spatial.sweep :as sweep]
            [clojure.test :refer [deftest is testing]]))

(defn- rng-doubles [seed n f]
  (let [rng (java.util.Random. seed)]
    (double-array (repeatedly n #(f (.nextDouble rng))))))

(defn- brute
  "Every overlapping pair, the way everyone writes it first."
  [^doubles mins ^doubles maxs n]
  (vec (for [i (range n) j (range (inc i) n)
             :when (every? (fn [k]
                             (and (<= (aget mins (+ (* 3 i) k)) (aget maxs (+ (* 3 j) k)))
                                  (>= (aget maxs (+ (* 3 i) k)) (aget mins (+ (* 3 j) k)))))
                           (range 3))]
         [i j])))

(defn- scene
  "`n` boxes in a 10-unit cube, all one size or wildly varying."
  [seed n {:keys [mixed?] :or {mixed? false}}]
  (let [pos (rng-doubles seed (* 3 n) #(* 10.0 %))
        half (if mixed?
               (rng-doubles (+ seed 1) (* 3 n) #(+ 0.05 (* 1.2 (Math/pow % 6))))
               (double-array (repeat (* 3 n) 0.35)))]
    (into [pos] (sweep/boxes-of pos half n))))

(deftest agreement-test
  (testing "sweep finds exactly the pairs brute force does"
    (doseq [n [0 1 2 5 50 400]]
      (let [[_ mins maxs] (scene 3 n {})
            s (sweep/sweep (max 1 n))]
        (is (= (sort (brute mins maxs n))
               (sort (sweep/overlapping-pairs s mins maxs n)))
            (str n " equal-sized boxes")))))

  (testing "including when the boxes are wildly different sizes"
    ;; The case a grid struggles with, because no one cell size suits
    ;; both a speck and something twenty times its width.
    (doseq [n [50 400]]
      (let [[_ mins maxs] (scene 5 n {:mixed? true})
            s (sweep/sweep n)]
        (is (= (sort (brute mins maxs n))
               (sort (sweep/overlapping-pairs s mins maxs n)))))))

  (testing "pairs come back ordered and once each"
    (let [[_ mins maxs] (scene 9 300 {})
          s (sweep/sweep 300)
          pairs (sweep/overlapping-pairs s mins maxs 300)]
      (is (every? (fn [[i j]] (< i j)) pairs))
      (is (= (count pairs) (count (distinct pairs))))))

  (testing "and it agrees with the spatial hash on what a sphere overlap is"
    ;; A hash pair means centres within the query distance; an AABB
    ;; overlap is a weaker condition, so every hash pair must be a sweep
    ;; pair but not the reverse.
    (let [n 400 r 0.35
          [pos mins maxs] (scene 13 n {})
          s (sweep/sweep n)
          swept (set (sweep/overlapping-pairs s mins maxs n))
          h (hash/spatial-hash (* 2 r) n)
          hashed (set (hash/overlapping-pairs h pos n (* 2 r)))]
      (is (seq hashed))
      (is (every? swept hashed)))))

(deftest axis-test
  (testing "the axis with the most spread is the one swept"
    (doseq [[axis scale] [[0 [20.0 1.0 1.0]] [1 [1.0 20.0 1.0]] [2 [1.0 1.0 20.0]]]]
      (let [n 300
            rng (java.util.Random. 17)
            pos (double-array (mapcat (fn [_] (map #(* % (.nextDouble rng)) scale))
                                      (range n)))
            [mins maxs] (sweep/spheres->boxes! pos n 0.2)]
        (is (= axis (sweep/axis-of mins maxs n))))))

  (testing "an empty scene picks something rather than failing"
    (is (= 0 (sweep/axis-of (double-array 0) (double-array 0) 0)))))

(deftest coherence-test
  ;; The reason to keep the order between frames rather than sort afresh.
  (let [n 1500 r 0.2
        pos (rng-doubles 21 (* 3 n) #(* 10.0 %))
        [mins maxs] (sweep/spheres->boxes! pos n r)
        s (sweep/sweep n)]
    (sweep/overlapping-pairs s mins maxs n)

    (testing "a sorted list has no inversions left"
      (is (zero? (sweep/inversions s mins n))))

    (testing "a small move leaves it nearly sorted"
      (let [rng (java.util.Random. 23)]
        (dotimes [i (* 3 n)]
          (aset ^doubles pos i (+ (aget ^doubles pos i) (- (* 0.02 (.nextDouble rng)) 0.01))))
        (sweep/spheres->boxes! pos n r mins maxs)
        (let [jiggled (sweep/inversions s mins n)]
          (is (< jiggled n)
              (str "after a small move: " jiggled " inversions over " n " objects"))

          (testing "and a reshuffle does not"
            (sweep/overlapping-pairs s mins maxs n)
            (let [rng (java.util.Random. 29)]
              (dotimes [i (* 3 n)] (aset ^doubles pos i (* 10.0 (.nextDouble rng)))))
            (sweep/spheres->boxes! pos n r mins maxs)
            (is (> (sweep/inversions s mins n) (* 4 jiggled))
                "a scrambled scene costs far more to re-sort")))))

    (testing "the answer is the same however scrambled the order was"
      (let [fresh (sweep/sweep n)]
        (is (= (sort (sweep/overlapping-pairs fresh mins maxs n))
               (sort (sweep/overlapping-pairs s mins maxs n))))))))

(deftest adjacency-test
  (testing "adjacency holds the same pairs, laid out for repeated reading"
    (let [n 300
          [_ mins maxs] (scene 31 n {})
          s (sweep/sweep n)
          pairs (sweep/overlapping-pairs s mins maxs n)
          {:keys [^ints starts ^ints ids] :as adj} (sweep/adjacency s mins maxs n)]
      (is (= (count pairs) (:pairs adj)))
      (is (= (count pairs) (alength ids)))
      (is (zero? (aget starts 0)))
      (is (= (count pairs) (aget starts n)))
      (testing "and every pair is recorded once, under the higher index"
        (let [from-adj (set (for [j (range n)
                                  k (range (aget starts j) (aget starts (inc j)))]
                              [(aget ids k) j]))]
          (is (= (set pairs) from-adj)))))))
