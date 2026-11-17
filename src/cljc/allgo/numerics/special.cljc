(ns allgo.numerics.special
  "Special functions: the modified Bessel functions of the first kind and
  the Legendre functions -- associated, and the derivatives of the
  polynomials that Pines's gravity formulation is written in -- and the
  gamma and incomplete beta functions that statistics needs, with
  Student's t distribution built on them; and the error function, with
  the normal distribution built on it."
  (:require [allgo.numerics.roots :as roots]
            [allgo.search :as search]
            [clojure.math :as math])
  #?(:clj (:import [allgo.numerics.roots Interval])))

(defn bessel-i
  "The modified Bessel function of the first kind, I_n(x), for integer n >=
  0, by its power series sum (x/2)^(2k+n) / (k! (k+n)!), which converges
  for every x."
  ^double [n ^double x]
  (let [h (* 0.5 x)
        first-term (loop [t 1.0 k 1] (if (> k n) t (recur (/ (* t h) k) (inc k))))
        hh (* h h)]
    (loop [t first-term sum first-term k 1]
      (let [t (/ (* t hh) (* k (+ k n)))
            sum' (+ sum t)]
        (if (or (= sum' sum) (> k 500)) sum' (recur t sum' (inc k)))))))

(defn associated-legendre
  "The unnormalized associated Legendre functions P_nm(x), without the
  Condon-Shortley phase, to degree `deg`, as a map keyed by [n m]: `x` and
  `s` = sqrt(1 - x^2) passed apart, so a caller that has both (the sine
  and cosine of a latitude) keeps s's digits near the poles."
  [sl cl deg]
  (reduce (fn [P [n m]]
            (assoc P [n m]
                   (cond
                     (and (zero? n) (zero? m)) 1.0
                     (= n m) (* (- (* 2.0 m) 1.0) cl (P [(dec m) (dec m)]))
                     (= n (inc m)) (* (+ (* 2.0 m) 1.0) sl (P [m m]))
                     :else (/ (- (* (- (* 2.0 n) 1.0) sl (P [(dec n) m]))
                                 (* (+ n m -1.0) (P [(- n 2) m])))
                              (- n m)))))
          {}
          (for [m (range (inc deg)) n (range m (inc deg))] [n m])))

(defn derived-legendre
  "A_nm(u), the m-th derivatives of the Legendre polynomials, to degree
  `deg`: A_nn = (2n-1) A_n-1,n-1, A_n,n-1 = u A_nn, and
  (n - m) A_nm = (2n - 1) u A_n-1,m - (n + m - 1) A_n-2,m."
  [u deg]
  (reduce (fn [A [n m]]
            (assoc A [n m]
                   (cond
                     (and (zero? n) (zero? m)) 1.0
                     (= n m) (* (- (* 2.0 n) 1.0) (A [(dec n) (dec n)]))
                     (= m (dec n)) (* u (A [n n]))
                     :else (/ (- (* (- (* 2.0 n) 1.0) u (A [(dec n) m]))
                                 (* (+ n m -1.0) (A [(- n 2) m])))
                              (- n m)))))
          {}
          (for [n (range (inc deg)) m (range n -1 -1)] [n m])))

;; ------------------------------------------------- gamma and beta

(def ^:private lanczos
  "Lanczos's approximation with g = 7 and nine terms, as Godfrey computed
  the coefficients (and as widely republished): good to 15 digits."
  [0.99999999999980993 676.5203681218851 -1259.1392167224028 771.32342877765313
   -176.61502916214059 12.507343278686905 -0.13857109526572012
   9.9843695780195716e-6 1.5056327351493116e-7])

(defn log-gamma
  "ln Gamma(x) for x > 0, by Lanczos's approximation (Lanczos 1964),
  reflected for x < 1/2."
  [x]
  (if (< x 0.5)
    (- (math/log (/ math/PI (abs (math/sin (* math/PI x))))) (log-gamma (- 1.0 x)))
    (let [x (- x 1.0)
          t (+ x 7.5)
          sum (reduce + (first lanczos) (map-indexed (fn [i c] (/ c (+ x i 1.0))) (rest lanczos)))]
      (+ (* 0.5 (math/log (* 2.0 math/PI))) (* (+ x 0.5) (math/log t)) (- t) (math/log sum)))))

(defn- beta-fraction
  "The continued fraction for the incomplete beta function (DLMF 8.17.22),
  by the modified Lentz method."
  [a b x]
  (let [tiny 1e-300
        clamp #(if (< (abs %) tiny) tiny %)]
    (loop [m 1 c 1.0 d (/ 1.0 (clamp (- 1.0 (/ (* (+ a b) x) (+ a 1.0))))) h d]
      (let [m2 (* 2 m)
            ;; the even step, then the odd
            aa (/ (* m (- b m) x) (* (+ a m2 -1.0) (+ a m2)))
            d (/ 1.0 (clamp (+ 1.0 (* aa d)))) c (clamp (+ 1.0 (/ aa c)))
            h (* h d c)
            aa (/ (- (* (+ a m) (+ a b m) x)) (* (+ a m2) (+ a m2 1.0)))
            d (/ 1.0 (clamp (+ 1.0 (* aa d)))) c (clamp (+ 1.0 (/ aa c)))
            del (* d c)
            h (* h del)]
        ;; a few ulps: |del - 1| < 1e-16 would demand del = 1 exactly
        (if (or (< (abs (- del 1.0)) (* 4.0 (math/ulp 1.0))) (> m 1000)) h (recur (inc m) c d h))))))

(defn regularized-beta
  "The regularized incomplete beta function I_x(a, b), from its continued
  fraction, using the symmetry I_x(a, b) = 1 - I_(1-x)(b, a) where the
  fraction converges faster."
  [x a b]
  (cond
    (<= x 0.0) 0.0
    (>= x 1.0) 1.0
    :else
    (let [front (math/exp (+ (- (log-gamma (+ a b)) (log-gamma a) (log-gamma b))
                             (* a (math/log x)) (* b (math/log (- 1.0 x)))))]
      (if (< x (/ (+ a 1.0) (+ a b 2.0)))
        (/ (* front (beta-fraction a b x)) a)
        (- 1.0 (/ (* front (beta-fraction b a (- 1.0 x))) b))))))

(defn student-t-cdf
  "The probability that Student's t with `nu` degrees of freedom is below
  `t`."
  [t nu]
  (let [p (* 0.5 (regularized-beta (/ nu (+ nu (* t t))) (* 0.5 nu) 0.5))]
    (if (pos? t) (- 1.0 p) p)))

(defn student-t-quantile
  "The t with probability `p` below it, `nu` degrees of freedom: the
  inverse of `student-t-cdf`, by bisection (`allgo.numerics.roots`'
  crossing step) to a bracket under 1e-13."
  [p nu]
  (let [below? #(neg? (- (student-t-cdf % nu) p))
        ^Interval br (search/settle (fn [^Interval i] (< (- (.-b i) (.-a i)) 1e-13)) 202
                                    (partial roots/crossing-step below?) (roots/->Interval -1e3 1e3))]
    (* 0.5 (+ (.-a br) (.-b br)))))

(defn- erf-series-step
  "The next term of erf's series with every term positive, x e^(-x^2)
  times the sum of (2x^2)^n / (1 3 5 ... (2n+1)), and the sum so far."
  [x2 {:keys [term sum n] :as state}]
  (when (> (abs term) (* 1e-17 (abs sum)))
    (let [term (/ (* term 2.0 x2) (+ (* 2.0 n) 3.0))]
      (assoc state :term term :sum (+ sum term) :n (inc n)))))

(defn- erf-series [x]
  (let [{:keys [sum]} (last (search/step-seq (partial erf-series-step (* x x)) {:term x :sum x :n 0}))]
    (* (/ 2.0 (math/sqrt math/PI)) (math/exp (- (* x x))) sum)))

(defn- erfc-fraction-step
  "One step of Lentz's method on the continued fraction x + (1/2)/(x +
  1/(x + (3/2)/(x + ...))), the k-th partial numerator k/2."
  [x {:keys [f c d k] :as state}]
  (let [a (* 0.5 k)
        d (let [d (+ x (* a d))] (/ 1.0 (if (zero? d) 1e-300 d)))
        c (let [c (+ x (/ a c))] (if (zero? c) 1e-300 c))
        delta (* c d)]
    (when (> (abs (- delta 1.0)) 1e-16)
      (assoc state :f (* f delta) :c c :d d :k (inc k)))))

(defn- erfc-fraction
  "erfc of x >= 2.5: e^(-x^2) / (sqrt(pi) times the continued fraction)."
  [x]
  (let [{:keys [f]} (last (take 500 (search/step-seq (partial erfc-fraction-step x) {:f x :c x :d 0.0 :k 1})))]
    (/ (math/exp (- (* x x))) (* (math/sqrt math/PI) f))))

(defn erfc
  "The complementary error function, 1 - erf(x): by the continued
  fraction from x = 2.5 on, keeping its tiny values' relative precision,
  and by erf's series below."
  [x]
  (cond (>= x 2.5) (erfc-fraction x)
        (>= x 0.0) (- 1.0 (erf-series x))
        (> x -2.5) (+ 1.0 (erf-series (- x)))
        :else (- 2.0 (erfc-fraction (- x)))))

(defn erf
  "The error function, 2/sqrt(pi) times the integral of e^(-t^2) from 0
  to x: by its series with every term positive near the origin, by erfc's
  continued fraction beyond 2.5."
  [x]
  (cond (>= x 2.5) (- 1.0 (erfc-fraction x))
        (<= x -2.5) (- (erfc-fraction (- x)) 1.0)
        (neg? x) (- (erf-series (- x)))
        :else (erf-series x)))

(defn normal-pdf
  "The standard normal density at `x`."
  [x]
  (/ (math/exp (* -0.5 x x)) (math/sqrt (* 2.0 math/PI))))

(defn normal-cdf
  "The probability that a standard normal is below `x`, from erfc, so
  that the far tails keep their relative precision."
  [x]
  (* 0.5 (erfc (- (/ x (math/sqrt 2.0))))))
