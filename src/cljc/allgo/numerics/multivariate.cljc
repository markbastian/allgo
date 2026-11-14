(ns allgo.numerics.multivariate
  "Polynomials in several variables, exactly: built from variables and
  constants by sums, products and powers, evaluated, and differentiated
  term by term -- so a gradient or a Hessian is itself a polynomial, with
  no step size to choose.

  A polynomial is a map from exponent vectors to coefficients:
  `{[2 0 0] 1.0, [0 0 0] -1.0}` is x^2 - 1 in three variables. Zero
  coefficients are dropped."
  (:require [clojure.math :as math]))

(defn- clean [p] (into {} (remove (comp zero? val)) p))

(defn constant
  "The constant `c` as a polynomial in `n` variables."
  [n c]
  (clean {(vec (repeat n 0)) (double c)}))

(defn variable
  "The `i`th of `n` variables, from 0."
  [n i]
  {(assoc (vec (repeat n 0)) i 1) 1.0})

(defn add
  "The sum of polynomials."
  [& ps]
  (clean (apply merge-with + ps)))

(defn scale
  "Polynomial `p` times the number `k`."
  [p k]
  (clean (update-vals p #(* % k))))

(defn sub
  "`p` minus `q`."
  [p q]
  (add p (scale q -1.0)))

(defn mul
  "The product of polynomials."
  ([p q]
   (clean (reduce (fn [acc [e c]] (update acc e (fnil + 0.0) c))
                  {}
                  (for [[e1 c1] p [e2 c2] q] [(mapv + e1 e2) (* c1 c2)]))))
  ([p q & more] (reduce mul (mul p q) more)))

(defn power
  "`p` to the non-negative integer power `k`."
  [p k]
  (let [n (count (key (first p)))]
    (reduce mul (constant n 1.0) (repeat k p))))

(defn evaluate
  "`p` at the point `x`."
  [p x]
  (reduce-kv (fn [acc e c]
               (+ acc (* c (reduce * 1.0 (map (fn [xi k] (math/pow xi k)) x e)))))
             0.0 p))

(defn evaluator
  "A function evaluating `p` at a point, faster than `evaluate` for one
  evaluated many times: its terms laid out once, powers by
  multiplication."
  [p]
  (let [terms (vec (for [[e c] p] [(double c) (vec e)]))]
    (fn [x]
      (reduce (fn [acc [c e]]
                (+ acc (* c (reduce-kv (fn [m i k]
                                         (let [xi (double (nth x i))]
                                           (loop [m m k (long k)] (if (pos? k) (recur (* m xi) (dec k)) m))))
                                       1.0 e))))
              0.0 terms))))

(defn derivative
  "The partial derivative of `p` with respect to variable `i`."
  [p i]
  (clean (into {} (for [[e c] p :let [k (nth e i)] :when (pos? k)]
                    [(update e i dec) (* c k)]))))

(defn gradient
  "The partial derivatives of `p` in `n` variables, as polynomials."
  [p n]
  (mapv #(derivative p %) (range n)))

(defn hessian
  "The second partial derivatives of `p` in `n` variables, as a matrix of
  polynomials."
  [p n]
  (mapv #(gradient % n) (gradient p n)))
