(ns allgo.polynomial-test
  "The closed-form roots against polynomials built from their roots --
  real, repeated, complex, biquadratic, of wide range -- and, for random
  coefficients, against the polynomial itself: every root a zero of it."
  (:require [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- expand
  "The coefficients, highest first, of the monic polynomial with complex
  `roots` (conjugate pairs, so the result is real)."
  [roots]
  (let [mul (fn [cs [re im]]
              ;; multiply by (x - z), complex coefficients as [re im]
              (let [z [re im]
                    shifted (conj (vec cs) [0.0 0.0])
                    scaled (vec (cons [0.0 0.0] (map (fn [[a b]] [(- (* a (first z)) (* b (second z)))
                                                                  (+ (* a (second z)) (* b (first z)))])
                                                     cs)))]
                (mapv (fn [[a b] [c d]] [(- a c) (- b d)]) shifted scaled)))]
    (mapv first (reduce mul [[1.0 0.0]] roots))))

(defn- matched?
  "Every expected root has a found one within `tol` (relative)."
  [expected found tol]
  (every? (fn [[re im]]
            (some (fn [[r i]] (<= (math/hypot (- r re) (- i im)) (* tol (max 1.0 (math/hypot re im))))) found))
          expected))

(deftest quadratics
  (is (matched? [[2.0 0.0] [-3.0 0.0]] (poly/quadratic 1.0 1.0 -6.0) 1e-15))
  (is (matched? [[1.0 2.0] [1.0 -2.0]] (poly/quadratic 1.0 -2.0 5.0) 1e-15))
  (testing "no cancellation in the small root"
    (let [[[r1] [r2]] (poly/quadratic 1.0 -1e8 1.0)]
      (is (< (abs (- (min r1 r2) 1e-8)) 1e-22)))))

(deftest cubics
  (doseq [roots [[[1.0 0.0] [2.0 0.0] [3.0 0.0]]
                 [[-1.5 0.0] [0.25 3.0] [0.25 -3.0]]
                 [[2.0 0.0] [2.0 0.0] [-5.0 0.0]]
                 [[1e3 0.0] [1e-3 0.0] [-7.0 0.0]]]]
    (let [[a b c d] (expand roots)]
      (is (matched? roots (poly/cubic a b c d) 1e-9) (str roots))))
  (testing "a triple root, to the cube root of the precision"
    (is (matched? [[4.0 0.0]] (poly/cubic 1.0 -12.0 48.0 -64.0) 1e-5))))

(deftest quartics
  (doseq [roots [[[1.0 0.0] [2.0 0.0] [3.0 0.0] [4.0 0.0]]
                 [[-2.0 0.0] [5.0 0.0] [1.0 1.0] [1.0 -1.0]]
                 [[0.5 2.0] [0.5 -2.0] [-3.0 0.25] [-3.0 -0.25]]
                 [[2.0 0.0] [-2.0 0.0] [0.0 3.0] [0.0 -3.0]]             ; biquadratic
                 [[3.0 0.0] [3.0 0.0] [-1.0 0.0] [7.0 0.0]]
                 [[165.1 0.0] [63.3 0.0] [98.6 25.3] [98.6 -25.3]]]]     ; Roberts's profile
    (let [[a b c d e] (expand roots)]
      (is (matched? roots (poly/quartic a b c d e) 1e-9) (str roots))))
  (testing "a leading coefficient other than 1"
    (let [[a b c d e] (map #(* -3.5 %) (expand [[1.0 0.0] [-4.0 0.0] [2.0 5.0] [2.0 -5.0]]))]
      (is (matched? [[1.0 0.0] [-4.0 0.0] [2.0 5.0] [2.0 -5.0]] (poly/quartic a b c d e) 1e-9)))))

(deftest residuals
  (testing "random coefficients: each root a zero of its polynomial, relative to the terms' size"
    (let [rng (java.util.Random. 20260928)
          coef #(- (* 20.0 (.nextDouble rng)) 10.0)]
      (dotimes [_ 200]
        (let [cs (vec (repeatedly 5 coef))
              roots (apply poly/quartic cs)]
          (is (= 4 (count roots)))
          (doseq [[re im] roots]
            (let [z (math/hypot re im)
                  ;; |p(z)| against the sum of the terms' magnitudes
                  size (reduce + (map-indexed (fn [i c] (* (abs c) (math/pow z (- 4 i)))) cs))
                  [pr pi] (reduce (fn [[a b] c] [(+ (- (* a re) (* b im)) c) (+ (* a im) (* b re))]) [0.0 0.0] cs)]
              (is (< (math/hypot pr pi) (* 1e-12 size)) (str cs)))))))))

(deftest real-roots
  (is (= [-3 2] (mapv #(math/round %) (poly/real-roots (poly/quadratic 1.0 1.0 -6.0)))))
  (is (= 2 (count (poly/real-roots (poly/quartic 1.0 -2.0 -1.0 -12.0 -10.0))))))
