(ns allgo.numerics.fourier
  "Fourier approximation (Chapra and Canale, *Numerical Methods for
  Engineers*, chapter 19): sinusoids fitted by least squares, the
  coefficients of a continuous Fourier series, the discrete Fourier
  transform taken directly, the Sande-Tukey fast transform (decimation in
  frequency; `allgo.numerics.fft` is Cooley-Tukey's decimation in time,
  on arrays), and the power spectrum.

  Transforms follow `allgo.numerics.fft`: forward with e^(-i 2 pi j k/N)
  and unscaled, the inverse carrying 1/N. Complex numbers are `[re im]`."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.quadrature :as quad]
            [clojure.math :as math]))

(defn sinusoids
  "The least-squares fit of y = A0 + sum over k of (A_k cos k w t + B_k
  sin k w t), `harmonics` of them (default 1), to `ts` `ys` at the
  fundamental angular frequency `w`: `{:a0 :a :b :amplitude :phase}` --
  each harmonic also as C_k cos(k w t + theta_k), C_k = sqrt(A^2 + B^2)
  and theta_k = atan2(-B, A). For data spread evenly over whole periods
  the normal equations are diagonal (Chapra and Canale 19.14-19.16); any
  spacing is solved here."
  ([ts ys w] (sinusoids ts ys w 1))
  ([ts ys w harmonics]
   (let [basis (into [(constantly 1.0)]
                     (mapcat (fn [k] [#(math/cos (* k w %)) #(math/sin (* k w %))]) (range 1 (inc harmonics))))
         Z (mapv (fn [t] (mapv #(% t) basis)) ts)
         Zt (lin/transpose Z)
         c (ls/gauss (lin/mat-mul Zt Z) (lin/mat-vec Zt ys))
         a (mapv #(c (dec (* 2 %))) (range 1 (inc harmonics)))
         b (mapv #(c (* 2 %)) (range 1 (inc harmonics)))]
     {:a0 (c 0) :a a :b b
      :amplitude (mapv #(math/hypot %1 %2) a b)
      :phase (mapv #(math/atan2 (- %2) %1) a b)})))

(defn series
  "The continuous Fourier series of the function `f` of period `T` to
  `n` harmonics: `{:a0 :a :b}` with a0 = (1/T) int f, a_k = (2/T) int f
  cos k w t and b_k = (2/T) int f sin k w t over one period, w = 2 pi/T
  -- the integrals by Simpson's rule on `:intervals` (default 2000)
  intervals, which for a smooth periodic f converges very fast. `f` may
  jump, as a square wave does, at the cost of accuracy near the jumps."
  ([f T n] (series f T n {}))
  ([f T n {:keys [intervals] :or {intervals 2000}}]
   (let [w (/ (* 2.0 math/PI) T)
         half (* 0.5 T)
         avg (fn [g] (/ (quad/simpson g (- half) half intervals) T))]
     {:a0 (avg f)
      :a (mapv (fn [k] (* 2.0 (avg #(* (f %) (math/cos (* k w %)))))) (range 1 (inc n)))
      :b (mapv (fn [k] (* 2.0 (avg #(* (f %) (math/sin (* k w %)))))) (range 1 (inc n)))})))

(defn- cmul [[a b] [c d]] [(- (* a c) (* b d)) (+ (* a d) (* b c))])

(defn- ->complex [x] (if (number? x) [(double x) 0.0] x))

(defn dft
  "The discrete Fourier transform of `xs` (numbers or `[re im]`), taken
  directly -- N^2 complex products, any length N: F_k = sum_n x_n
  e^(-i 2 pi k n/N). With `:inverse? true` the inverse, carrying 1/N."
  ([xs] (dft xs {}))
  ([xs {:keys [inverse?]}]
   (let [xs (mapv ->complex xs)
         n (count xs)
         sign (if inverse? 1.0 -1.0)
         scale (if inverse? (/ 1.0 n) 1.0)]
     (mapv (fn [k]
             (let [[re im] (reduce (fn [acc j]
                                     (let [ang (/ (* sign 2.0 math/PI k j) n)]
                                       (mapv + acc (cmul (xs j) [(math/cos ang) (math/sin ang)]))))
                                   [0.0 0.0] (range n))]
               [(* scale re) (* scale im)]))
           (range n)))))

(defn sande-tukey
  "The fast Fourier transform by Sande and Tukey's decimation in
  frequency (Gentleman and Sande 1966): the sequence's first and second
  halves combined into sums, which give the even-indexed outputs, and
  twiddled differences, which give the odd -- each half a transform of
  half the length, recursively. The length must be a power of two. The
  result is in natural order."
  [xs]
  (let [xs (mapv ->complex xs)
        n (count xs)]
    (if (= n 1)
      xs
      (let [h (quot n 2)
            top (subvec xs 0 h) bottom (subvec xs h)
            sums (mapv #(mapv + %1 %2) top bottom)
            diffs (mapv (fn [j a b]
                          (let [ang (/ (* -2.0 math/PI j) n)]
                            (cmul (mapv - a b) [(math/cos ang) (math/sin ang)])))
                        (range h) top bottom)
            evens (sande-tukey sums)
            odds (sande-tukey diffs)]
        (vec (interleave evens odds))))))

(defn power-spectrum
  "The power in each frequency of the real signal `xs` sampled every
  `dt`: `[[frequency power] ...]` for frequencies 0 to the Nyquist 1/(2 dt)
  in steps of 1/(N dt), the power |F_k|^2/N^2 doubled for every frequency
  but zero and the Nyquist, which have no mirror image. Any length (by
  the direct transform when N is not a power of two)."
  [xs dt]
  (let [n (count xs)
        pow2? (and (pos? n) (zero? (bit-and n (dec n))))
        F (if pow2? (sande-tukey xs) (dft xs))]
    (vec (for [k (range (inc (quot n 2)))]
           (let [[re im] (F k)
                 p (/ (+ (* re re) (* im im)) (* n n))]
             [(/ k (* n dt)) (if (or (zero? k) (and (even? n) (= k (quot n 2)))) p (* 2.0 p))])))))
