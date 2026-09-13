(ns allgo.fft-test
  (:require [allgo.numerics.fft :as fft]
            [clojure.test :refer [deftest is testing]]))

(defn- mag [[r i]] (Math/hypot (double r) (double i)))
(defn- close? [a b] (< (abs (- (double a) (double b))) 1e-9))

(deftest power-of-two-test
  (testing "only the powers of two, and not zero"
    (is (every? fft/power-of-two? [1 2 4 8 16 1024]))
    (is (not-any? fft/power-of-two? [0 3 5 6 7 12 1000]))
    (is (not (fft/power-of-two? -8)))))

(deftest known-transforms-test
  (testing "a delta transforms to a flat spectrum -- every frequency, equally"
    (let [s (fft/spectrum [1 0 0 0 0 0 0 0])]
      (is (every? #(close? 1.0 (mag %)) s))))

  (testing "and a constant to a single spike at zero"
    (let [s (fft/spectrum (repeat 8 3.0))]
      (is (close? 24.0 (mag (first s))))
      (is (every? #(close? 0.0 (mag %)) (rest s)))))

  (testing "a pure tone gives two spikes, at k and at n-k"
    ;; Two, because a real signal cannot tell a positive frequency from
    ;; its negative twin -- which is the whole reason the solver below
    ;; can treat wavenumber m and its conjugate as one real system.
    (let [n 32 k 5
          xs (mapv #(Math/cos (/ (* 2 Math/PI k %) n)) (range n))
          s (fft/spectrum xs)]
      (doseq [[i v] (map-indexed vector (map mag s))]
        (if (#{k (- n k)} i)
          (is (close? (/ n 2.0) v))
          (is (close? 0.0 v))))))

  (testing "the four-point transform, worked out by hand"
    (let [s (fft/spectrum [1 2 3 4])]
      (is (close? 10.0 (first (nth s 0))))
      (is (close? -2.0 (first (nth s 1))))
      (is (close? 2.0 (second (nth s 1))))
      (is (close? -2.0 (first (nth s 2))))
      (is (close? 0.0 (second (nth s 2)))))))

(deftest round-trip-test
  (testing "inverse undoes forward, and carries all the scaling"
    (doseq [n [2 8 64 256]]
      (let [orig (mapv #(+ (Math/sin (* 0.3 %)) (* 0.4 (Math/cos (* 1.7 %)))) (range n))
            re (double-array orig)
            im (double-array n)]
        (fft/forward! re im)
        (fft/inverse! re im)
        (is (every? true? (map #(close? %1 %2) orig (vec re))))
        (is (every? #(close? 0.0 %) (vec im))))))

  (testing "Parseval: the transform moves energy about but does not make any"
    (let [n 64
          xs (mapv #(Math/sin (* 0.37 %)) (range n))
          s (fft/spectrum xs)
          space (reduce + (map #(* % %) xs))
          freq (/ (reduce + (map #(let [m (mag %)] (* m m)) s)) n)]
      (is (< (abs (- space freq)) 1e-9)))))

(deftest linearity-test
  (testing "the transform of a sum is the sum of the transforms"
    (let [n 32
          a (mapv #(Math/sin (* 0.3 %)) (range n))
          b (mapv #(Math/cos (* 0.9 %)) (range n))
          sa (fft/spectrum a)
          sb (fft/spectrum b)
          sab (fft/spectrum (map + a b))]
      (is (every? true?
                  (map (fn [x y z]
                         (and (close? (+ (first x) (first y)) (first z))
                              (close? (+ (second x) (second y)) (second z))))
                       sa sb sab))))))

(deftest rejects-test
  (testing "lengths that are not powers of two, and mismatched parts"
    (is (thrown? clojure.lang.ExceptionInfo
                 (fft/forward! (double-array 6) (double-array 6))))
    (is (thrown? clojure.lang.ExceptionInfo
                 (fft/forward! (double-array 8) (double-array 4))))))
