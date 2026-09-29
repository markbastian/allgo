(ns allgo.chapra-curve-fitting-test
  "Curve fitting (Chapra and Canale, *Numerical Methods for Engineers*,
  7th ed., chapters 17-19): regression, interpolation and splines, and
  Fourier approximation, against the book's worked examples and exact
  answers."
  (:require [allgo.array :as a]
            [allgo.numerics.fft :as fft]
            [allgo.numerics.fourier :as fourier]
            [allgo.numerics.interpolation :as interp]
            [allgo.numerics.regression :as reg]
            [allgo.numerics.special :as special]
            [allgo.numerics.splines :as splines]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(deftest regression
  (testing "Examples 17.1-17.2: the line through Table 17.1, and its statistics"
    (let [{:keys [a0 a1 sr st syx r2 r]} (reg/linear [1 2 3 4 5 6 7] [0.5 2.5 2.0 4.0 3.5 6.0 5.5])]
      (is (close? [a0 a1] [0.07142857 0.8392857] 1e-7))
      (is (close? [st sr] [22.7143 2.9911] 1e-4))
      (is (close? syx 0.7735 1e-4))
      (is (close? [r2 r] [0.868 0.932] 1e-3))))
  (testing "Example 17.5: the least-squares parabola through Table 17.4"
    (let [{:keys [coeffs syx sr st]} (reg/polynomial [0 1 2 3 4 5] [2.1 7.7 13.6 27.2 40.9 61.1] 2)]
      (is (close? coeffs [1.86071 2.35929 2.47857] 1e-5))
      (is (close? [st sr] [2513.39 3.74657] 1e-2))
      (is (close? syx 1.12 5e-3))))
  (testing "Example 17.6: multiple linear regression recovers y = 5 + 4 x1 - 3 x2"
    (is (close? (:a (reg/multiple [[0 0] [2 1] [2.5 2] [1 3] [4 6] [7 2]] [5 10 9 0 3 27])) [5.0 4.0 -3.0] 1e-10)))
  (testing "Example 17.7: standard errors and 95% confidence intervals, Table 17.2's model against measurement"
    (let [measured [10.00 16.30 23.00 27.50 31.00 35.60 39.00 41.50 42.90 45.00 46.00 45.50 46.00 49.00 50.00]
          model [8.953 16.405 22.607 27.769 32.065 35.641 38.617 41.095 43.156 44.872 46.301 47.490 48.479 49.303 49.988]
          {:keys [a syx se ci]} (reg/general [(constantly 1.0) identity] measured model)]
      (is (close? a [-0.85872 1.031592] 1e-5))
      (is (close? syx 0.863403 1e-6))
      (is (close? se [0.716372 0.018625] 1e-6))
      (testing "with Student's t for 13 degrees of freedom, 2.160368 (Excel's TINV)"
        (is (close? (special/student-t-quantile 0.975 13) 2.160368 1e-6))
        (is (close? ci [[(- -0.85872 (* 2.160368 0.716372)) (+ -0.85872 (* 2.160368 0.716372))]
                        [(- 1.031592 (* 2.160368 0.018625)) (+ 1.031592 (* 2.160368 0.018625))]]
                    2e-5)))))
  (testing "linearizations recover their exact models"
    (let [xs [1.0 2.0 3.0 4.0 5.0]]
      (is (close? ((juxt :a :b) (reg/exponential xs (map #(* 2.5 (math/exp (* 0.3 %))) xs))) [2.5 0.3] 1e-12))
      (is (close? ((juxt :a :b) (reg/power xs (map #(* 0.7 (math/pow % 1.8)) xs))) [0.7 1.8] 1e-12))
      (is (close? ((juxt :a :b) (reg/saturation-growth xs (map #(/ (* 4.0 %) (+ 2.0 %)) xs))) [4.0 2.0] 1e-12))))
  (testing "Example 17.8: Gauss-Newton for a0 (1 - e^(-a1 x)) from (1, 1)"
    (let [{:keys [p]} (reg/gauss-newton (fn [x [a0 a1]] (* a0 (- 1.0 (math/exp (- (* a1 x))))))
                                        [1.0 1.0] [0.25 0.75 1.25 1.75 2.25] [0.28 0.57 0.68 0.74 0.79])]
      (is (close? p [0.79186 1.6751] 1e-4)))))

(deftest interpolation
  ;; ln x, at the precision the book computes with
  (let [pts (mapv (fn [x] [x (math/log x)]) [1.0 4.0 6.0 5.0])]
    (testing "Example 18.3: the cubic Newton polynomial estimates ln 2 as 0.6287686"
      (is (close? (interp/divided-differences pts) [0.0 0.4620981 -0.05187311 0.007865529] 5e-8))
      (is (close? (:value (interp/newton pts 2.0)) 0.6287686 1e-7)))
    (testing "and agrees with Lagrange's form of the same polynomial"
      (is (close? (:value (interp/newton pts 3.3)) (interp/lagrange pts 3.3) 1e-12))))
  (testing "inverse interpolation: where the interpolant of 1/x through 2..6 is 0.3"
    (let [pts (map (fn [x] [x (/ 1.0 x)]) [2.0 3.0 4.0 5.0 6.0])
          x (interp/inverse pts 0.3)]
      (is (close? (interp/lagrange pts x) 0.3 1e-12))
      (is (close? x (/ 1.0 0.3) 0.02)))))

(deftest spline-fits
  (let [xs [3.0 4.5 7.0 9.0] ys [2.5 1.0 2.5 0.5]]
    (testing "Example 18.8: first-order splines, 1.3 at x = 5"
      (is (close? ((splines/linear xs ys) 5.0) 1.3 1e-12)))
    (testing "Example 18.9: quadratic splines -- the middle piece 0.64 x^2 - 6.76 x + 18.46, 0.66 at 5"
      (let [s (splines/quadratic xs ys)]
        (is (close? (s 5.0) 0.66 1e-12))
        (is (close? (s 3.7) (+ (- 3.7) 5.5) 1e-12) "the first piece a line")))
    (testing "Example 18.10: natural cubic splines -- f''(4.5) = 1.67909, f''(7) = -1.53308, 1.102886 at 5"
      (let [s (splines/cubic xs ys)]
        (is (close? (:second-derivatives (meta s)) [0.0 1.67909 -1.53308 0.0] 1e-5))
        ;; the book evaluates its six-digit coefficients
        (is (close? (s 5.0) 1.102886 1e-5))
        (is (close? (map s xs) ys 1e-12))))
    (testing "clamped and not-a-knot splines reproduce a cubic exactly"
      (let [f #(+ (* 0.5 % % %) (* -2.0 % %) % 3.0)
            df #(+ (* 1.5 % %) (* -4.0 %) 1.0)
            kx [0.0 1.0 2.5 3.0 4.2 5.0]]
        (doseq [end [[(df 0.0) (df 5.0)] :not-a-knot]]
          (let [s (splines/cubic kx (map f kx) {:end end})]
            (is (every? #(close? (s %) (f %) 1e-10) [0.3 1.7 2.9 3.6 4.9]))))))
    (testing "Example 18.11: bilinear interpolation on a heated plate, 61.2143"
      (is (close? (splines/bilinear 2.0 9.0 1.0 6.0 60.0 57.5 55.0 70.0 5.25 4.8) 61.2143 1e-4)))))

(deftest fourier-approximation
  (testing "Example 19.1: the sinusoid 1.7 + cos(4.189 t + 1.0472) fitted from ten samples"
    (let [ts (map #(* 0.15 %) (range 10))
          ys (map #(+ 1.7 (math/cos (+ (* 4.189 %) 1.0472))) ts)
          {:keys [a0 a b amplitude phase]} (fourier/sinusoids ts ys 4.189)]
      (is (close? [a0 (a 0) (b 0)] [1.7 0.5 -0.866] 1e-3))
      (is (close? [(amplitude 0) (phase 0)] [1.0 1.0472] 1e-9))))
  (testing "Example 19.2: the square wave's series, a_k = +/-4/(k pi) for odd k"
    (let [T (* 2.0 math/PI)
          square (fn [t] (if (< (abs t) (/ T 4.0)) 1.0 -1.0))
          {:keys [a0 a b]} (fourier/series square T 7 {:intervals 20000})]
      (is (close? a0 0.0 1e-3))
      (is (close? a (map (fn [k] (if (odd? k) (* (if (= 1 (mod k 4)) 1.0 -1.0) (/ 4.0 (* k math/PI))) 0.0)) (range 1 8)) 2e-3))
      (is (close? b (repeat 7 0.0) 1e-9))))
  (testing "the direct DFT, Sande-Tukey and Cooley-Tukey agree, and the inverse undoes"
    (let [xs (mapv #(+ (math/sin (* 0.7 %)) (* 0.3 (math/cos (* 2.1 %)))) (range 16))
          direct (fourier/dft xs)
          dif (fourier/sande-tukey xs)
          re (a/f64 xs) im (a/f64 16)]
      (fft/forward! re im)
      (is (close? direct dif 1e-12))
      (is (close? (map first direct) (vec re) 1e-12))
      (is (close? (map second direct) (vec im) 1e-12))
      (is (close? (map first (fourier/dft direct {:inverse? true})) xs 1e-12))))
  (testing "the power spectrum puts a pure tone's power at its frequency"
    (let [dt 0.01 n 64
          xs (map #(* 3.0 (math/sin (* 2.0 math/PI 12.5 % dt))) (range n))
          spec (fourier/power-spectrum xs dt)
          [f p] (apply max-key second spec)]
      (is (close? f 12.5 1e-12))
      ;; a sine of amplitude A carries A^2/2 of power
      (is (close? p 4.5 1e-9))
      (is (close? (reduce + (map second spec)) 4.5 1e-9)))))
