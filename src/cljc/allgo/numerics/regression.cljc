(ns allgo.numerics.regression
  "Least-squares regression (Chapra and Canale, *Numerical Methods for
  Engineers*, chapter 17): straight lines with their error statistics,
  polynomials, several independent variables, any linear combination of
  basis functions -- with the coefficients' covariance, standard errors
  and confidence intervals -- the linearizations that fit exponential,
  power and saturation-growth models as lines, and nonlinear models by
  Gauss-Newton.

  Every fit returns its statistics alongside its coefficients: `:sr`, the
  sum of the squared residuals; `:st`, the spread about the mean; `:r2`,
  the coefficient of determination (st - sr)/st, the fraction of the
  spread the model explains; and `:syx`, the standard error of the
  estimate, sqrt(sr / (n - parameters))."
  (:require [allgo.numerics.differentiation :as d]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.polynomial :as poly]
            [allgo.numerics.special :as special]
            [clojure.math :as math]))

(defn- stats [ys predicted n-params]
  (let [n (count ys)
        mean (/ (reduce + ys) n)
        st (reduce + (map #(let [d (- % mean)] (* d d)) ys))
        sr (reduce + (map #(let [d (- %1 %2)] (* d d)) ys predicted))]
    {:sr sr :st st :r2 (if (zero? st) 1.0 (/ (- st sr) st))
     :syx (if (> n n-params) (math/sqrt (/ sr (- n n-params))) 0.0)}))

(defn general
  "The least-squares fit of y = a_0 z_0(x) + a_1 z_1(x) + ... to the data
  `xs` `ys`, `basis` the functions z_j of a data point x: the normal
  equations Z^T Z a = Z^T y solved by Gauss elimination. Returns `{:a
  :covariance :se :ci ...}` with the statistics -- the covariance
  syx^2 (Z^T Z)^-1, each coefficient's standard error its diagonal's root,
  and `:ci` the intervals a_j +/- t se_j at the `:confidence` (default
  0.95), t Student's for n - m degrees of freedom."
  ([basis xs ys] (general basis xs ys {}))
  ([basis xs ys {:keys [confidence] :or {confidence 0.95}}]
   (let [Z (mapv (fn [x] (mapv #(% x) basis)) xs)
         Zt (lin/transpose Z)
         ZtZ (lin/mat-mul Zt Z)
         a (ls/gauss ZtZ (lin/mat-vec Zt ys))
         m (count basis) n (count ys)
         {:keys [syx] :as st} (stats ys (lin/mat-vec Z a) m)
         inv (ls/inverse ZtZ)
         cov (mapv (fn [row] (mapv #(* syx syx %) row)) inv)
         se (mapv #(math/sqrt (get-in cov [% %])) (range m))
         t (when (> n m) (special/student-t-quantile (- 1.0 (* 0.5 (- 1.0 confidence))) (- n m)))]
     (merge st {:a a :covariance cov :se se
                :ci (when t (mapv (fn [ai s] [(- ai (* t s)) (+ ai (* t s))]) a se))}))))

(defn linear
  "The straight line y = a0 + a1 x through `xs` `ys`: `{:a0 :a1 :r ...}`
  with the statistics, `:r` the correlation coefficient."
  [xs ys]
  (let [{[a0 a1] :a :as fit} (general [(constantly 1.0) identity] xs ys)]
    (assoc fit :a0 a0 :a1 a1 :r (* (if (neg? a1) -1.0 1.0) (math/sqrt (max 0.0 (:r2 fit)))))))

(defn polynomial
  "The polynomial of degree `m` through `xs` `ys` by least squares:
  `{:coeffs ...}`, coefficients highest power first (as
  `allgo.numerics.polynomial` takes them), with the statistics."
  [xs ys m]
  (let [fit (general (mapv (fn [k] #(math/pow % k)) (range (inc m))) xs ys)]
    (assoc fit :coeffs (vec (reverse (:a fit))))))

(defn multiple
  "y = a0 + a1 x1 + a2 x2 + ... through the rows `X` of independent
  variables and `ys`: `{:a ...}`, a0 first, with the statistics."
  [X ys]
  (general (into [(constantly 1.0)] (map (fn [j] #(nth % j)) (range (count (first X))))) X ys))

;; -------------------------------------------------- linearizations

(defn exponential
  "y = a e^(b x) by a line through ln y against x: `{:a :b ...}`."
  [xs ys]
  (let [{:keys [a0 a1] :as fit} (linear xs (map math/log ys))]
    (assoc fit :a (math/exp a0) :b a1)))

(defn power
  "y = a x^b by a line through log y against log x: `{:a :b ...}`."
  [xs ys]
  (let [{:keys [a0 a1] :as fit} (linear (map math/log10 xs) (map math/log10 ys))]
    (assoc fit :a (math/pow 10.0 a0) :b a1)))

(defn saturation-growth
  "y = a x / (b + x) by a line through 1/y against 1/x: `{:a :b ...}`."
  [xs ys]
  (let [{:keys [a0 a1] :as fit} (linear (map #(/ 1.0 %) xs) (map #(/ 1.0 %) ys))]
    (assoc fit :a (/ 1.0 a0) :b (/ a1 a0))))

;; -------------------------------------------------------- nonlinear

(defn gauss-newton
  "The parameters `p0` of the model `(f x p)` that fit `xs` `ys` by least
  squares, by Gauss-Newton: the model linearized about the current
  parameters by its partials (central differences), the linear least
  squares for the correction solved, and repeated until the correction is
  below `:tol` relative (default 1e-12) or `:max-iter` (100). `{:p
  :iterations ...}` with the statistics. Can diverge from a poor start;
  Marquardt's damping (`allgo.numerics.optimize/marquardt`) is the
  remedy."
  ([f p0 xs ys] (gauss-newton f p0 xs ys {}))
  ([f p0 xs ys {:keys [tol max-iter] :or {tol 1e-12 max-iter 100}}]
   (let [m (count p0)]
     (loop [p (mapv double p0) k 1]
       (let [pred (mapv #(f % p) xs)
             ;; the model's partials in the parameters, one row per datum
             Z (d/jacobian (fn [p] (mapv #(f % p) xs)) p
                           {:steps (mapv #(* 1e-7 (max 1.0 (abs %))) p) :richardson? false})
             Zt (lin/transpose Z)
             dp (ls/gauss (lin/mat-mul Zt Z) (lin/mat-vec Zt (mapv - ys pred)))
             p' (mapv + p dp)]
         (if (or (>= k max-iter)
                 (every? true? (map #(<= (abs %1) (* tol (max 1.0 (abs %2)))) dp p')))
           (merge (stats ys (mapv #(f % p') xs) m) {:p p' :iterations k})
           (recur p' (inc k))))))))

(defn smoothed-derivative
  "The derivative at `x` of noisy data `xs` `ys` taken from their
  least-squares polynomial of `degree` rather than from the data
  themselves (Chapra and Canale 23.4): differentiation amplifies noise,
  and fitting first smooths it away."
  [xs ys degree x]
  (let [{:keys [coeffs]} (polynomial xs ys degree)]
    (second (poly/derivatives coeffs x 1))))
