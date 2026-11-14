(ns allgo.probability.gaussian
  "Multivariate Gaussians."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn log-density
  "log N(x; mu, C), C symmetric positive definite: through its Cholesky
  factor, whose diagonal gives the log determinant."
  [mu C x]
  (let [L (lin/cholesky C) r (lin/sub x mu)]
    (* -0.5 (+ (* (count x) (math/log (* 2.0 math/PI)))
               (* 2.0 (reduce + (map-indexed (fn [i row] (math/log (row i))) L)))
               (lin/dot r (lin/cholesky-solve C r))))))
