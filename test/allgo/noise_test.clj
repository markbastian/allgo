(ns allgo.noise-test
  (:require [allgo.procedural.noise :as n]
            [clojure.test :refer [deftest is testing]]))

(defn- points
  "A deterministic walk that avoids the lattice, so a test cannot pass by
  accidentally sampling only the points where noise is zero."
  [n]
  (for [i (range n)]
    [(* 0.1379 i) (* 0.2713 i) (* 0.3137 i)]))

(defn- values [f n] (mapv (fn [[x y z]] (f x y z)) (points n)))

(deftest prng-test
  (testing "the stream is in range and never falls into the fixed point"
    (loop [i 0 s (n/cell-seed 0 1 2 3)]
      (when (< i 10000)
        (is (pos? s))
        (is (< s 2147483647.0))
        (is (<= 0.0 (n/unit s)))
        (is (< (n/unit s) 1.0))
        (recur (inc i) (n/advance s)))))

  (testing "a cell's stream is a function of the cell, and of nothing else"
    (is (== (n/cell-seed 3 10 -4 7) (n/cell-seed 3 10 -4 7)))
    (is (not= (n/cell-seed 3 10 -4 7) (n/cell-seed 3 10 -4 8)))
    (is (not= (n/cell-seed 3 10 -4 7) (n/cell-seed 4 10 -4 7))))

  (testing "Poisson counts average out to the mean asked for"
    (let [mean 3.0
          counts (map #(n/poisson-count mean (/ (double %) 20000)) (range 20000))]
      (is (< (abs (- mean (/ (double (reduce + counts)) 20000))) 0.05)))))

(deftest permutation-test
  (testing "the table is a permutation of 0..255, twice"
    (let [p (n/permutation 12)]
      (is (= 512 (alength p)))
      (is (= (set (range 256)) (set (take 256 (seq p)))))
      (is (= (vec (take 256 (seq p))) (vec (drop 256 (seq p)))))))

  (testing "and a different seed shuffles it differently"
    (is (not= (vec (seq (n/permutation 1))) (vec (seq (n/permutation 2)))))))

(deftest gradient-basis-test
  (let [g (n/gradient-basis {:seed 5})]
    (testing "gradient noise is exactly zero on the lattice"
      ;; The defining property, and what keeps the grid invisible: the
      ;; value comes from the gradients, and at a lattice point the
      ;; offset the gradient is dotted with is the zero vector.
      (doseq [i (range -3 4) j (range -3 4)]
        (is (zero? (g (double i) (double j) 2.0)))))

    (testing "and lands inside [-1, 1] everywhere else"
      (let [vs (values g 20000)]
        (is (> (apply min vs) -1.0))
        (is (< (apply max vs) 1.0))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.05))))

    (testing "it is continuous: a small step gives a small change"
      (doseq [[x y z] (points 200)]
        (is (< (abs (- (g x y z) (g (+ x 1e-5) y z))) 1e-3))))

    (testing "a seed is a world"
      (is (= (values g 500) (values (n/gradient-basis {:seed 5}) 500)))
      (is (not= (values g 500) (values (n/gradient-basis {:seed 6}) 500))))))

(deftest periodic-test
  (testing "a periodic basis repeats exactly, on every axis"
    ;; What makes a field safe to bake into a tile. Not bit-for-bit: a
    ;; coordinate offset by a whole period comes back with a slightly
    ;; different fractional part, so the two agree to floating point and
    ;; not to the last bit. That is the standard a seam has to meet.
    (doseq [make [n/gradient-basis n/value-basis]
            period [4 8 16]]
      (let [b (make {:seed 5 :period period})
            p (double period)
            same? (fn [a b] (< (abs (- (double a) (double b))) 1e-12))]
        (doseq [[x y z] (points 300)]
          (is (same? (b x y z) (b (+ x p) y z)))
          (is (same? (b x y z) (b x (- y (* 2 p)) z)))
          (is (same? (b x y z) (b (+ x (* 3 p)) (+ y p) (- z p))))))))

  (testing "it is still continuous across the seam"
    (let [b (n/gradient-basis {:seed 5 :period 8})]
      (doseq [t (range 0.0 1.0 0.05)]
        (is (< (abs (- (b (- 8.0 1e-6) t 0.3) (b (+ 8.0 1e-6) t 0.3))) 1e-5)))))

  (testing "a period does not flatten the field"
    (let [b (n/gradient-basis {:seed 5 :period 8})
          vs (values b 4000)]
      (is (> (- (apply max vs) (apply min vs)) 1.0))))

  (testing "and asking for none leaves the basis exactly as it was"
    (is (= (values (n/gradient-basis {:seed 5}) 500)
           (values (n/gradient-basis {:seed 5 :period nil}) 500)))))

(deftest value-basis-test
  (let [v (n/value-basis {:seed 5})]
    (testing "value noise is not zero on the lattice -- that is the point of it"
      (is (some #(not (zero? %)) (for [i (range 6)] (v (double i) 0.0 0.0)))))

    (testing "it stays inside the table's range"
      (let [vs (values v 20000)]
        (is (>= (apply min vs) -1.0))
        (is (<= (apply max vs) 1.0))))

    (testing "and repeats the lattice value exactly at the lattice"
      ;; Two visits to the same integer point have to agree, or the
      ;; hash is not a hash.
      (is (== (v 3.0 -7.0 11.0) (v 3.0 -7.0 11.0))))))

(deftest gradient-4d-test
  (let [g4 (n/gradient-basis-4d {:seed 5})]
    (testing "zero on the 4D lattice, for the same reason as in 3D"
      (doseq [i (range -2 3) l (range -2 3)]
        (is (zero? (g4 (double i) 1.0 2.0 (double l))))))

    (testing "in range"
      (let [vs (for [i (range 8000)]
                 (g4 (* 0.137 i) (* 0.271 i) (* 0.313 i) (* 0.077 i)))]
        (is (> (apply min vs) -1.0))
        (is (< (apply max vs) 1.0))))

    (testing "a slice is a 3D basis, and moving the slice moves the world"
      (let [a (n/slice g4 0.25)
            b (n/slice g4 0.25)
            c (n/slice g4 0.75)]
        (is (= (values a 200) (values b 200)))
        (is (not= (values a 200) (values c 200)))))

    (testing "and neighbouring slices are near each other, which is what
              makes the fourth dimension usable for animation"
      (let [a (n/slice g4 0.500)
            b (n/slice g4 0.501)]
        (is (every? #(< % 0.02)
                    (map (fn [x y] (abs (- (double x) (double y))))
                         (values a 200) (values b 200))))))))

(deftest sparse-convolution-test
  (let [s (n/sparse-convolution-basis {:seed 5})]
    (testing "it is finite and roughly unit range"
      (let [vs (values s 8000)]
        (is (every? #(Double/isFinite (double %)) vs))
        (is (> (apply min vs) -2.0))
        (is (< (apply max vs) 2.0))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.1))))

    (testing "density changes the texture, not the range"
      ;; The normalisation exists so that a density knob does not double
      ;; as a brightness knob.
      (doseq [d [1.0 3.0 8.0]]
        (let [vs (values (n/sparse-convolution-basis {:seed 5 :density d}) 4000)]
          (is (< 0.5 (- (apply max vs) (apply min vs)) 3.0)))))))

(deftest vector-basis-test
  (let [vb (n/vector-basis {:seed 5})]
    (testing "three components, all in range, and none of them each other"
      (let [vs (map (fn [[x y z]] (vb x y z)) (points 2000))]
        (is (every? #(= 3 (count %)) vs))
        (is (every? (fn [v] (every? #(< (abs (double %)) 1.0) v)) vs))
        (is (not= (map first vs) (map second vs)))
        (is (not= (map second vs) (map #(nth % 2) vs)))))))

(deftest transforms-test
  (let [g (n/gradient-basis {:seed 5})]
    (testing "absolute folds the negative half up"
      (let [a (n/absolute g)]
        (doseq [[x y z] (points 500)]
          (is (== (a x y z) (abs (g x y z)))))))

    (testing "ridge is never negative and peaks where the basis crosses zero"
      (let [r (n/ridge g 1.0)]
        (doseq [[x y z] (points 500)]
          (is (>= (r x y z) 0.0)))
        ;; Exactly on the lattice the basis is zero, so the ridge is at
        ;; its offset squared.
        (is (< (abs (- 1.0 (r 2.0 3.0 4.0))) 1e-12))))

    (testing "scaled multiplies domain and range"
      (let [s (n/scaled g 2.0 3.0)]
        (doseq [[x y z] (points 200)]
          (is (< (abs (- (s x y z) (* 3.0 (g (* 2.0 x) (* 2.0 y) (* 2.0 z))))) 1e-12)))))

    (testing "translated shifts the domain"
      (let [t (n/translated g [1.5 -2.5 0.25])]
        (doseq [[x y z] (points 200)]
          (is (== (t x y z) (g (+ x 1.5) (- y 2.5) (+ z 0.25)))))))

    (testing "distortion is the identity when there is none to apply"
      (let [d (n/distorted g (n/vector-basis {:seed 9}) 0.0)]
        (doseq [[x y z] (points 200)]
          (is (== (d x y z) (g x y z))))))

    (testing "and changes things when there is"
      (let [d (n/distorted g (n/vector-basis {:seed 9}) 0.5)]
        (is (not= (values d 500) (values g 500)))))

    (testing "stepping quantises the range"
      (let [st (n/stepped g 4 0.0)
            vs (set (values st 4000))]
        (is (<= (count vs) 5))))

    (testing "sample calls a basis and agrees with calling it directly"
      (doseq [[x y z] (points 200)]
        (is (== (n/sample g x y z) (g x y z)))))))
