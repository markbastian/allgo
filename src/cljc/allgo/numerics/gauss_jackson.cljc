(ns allgo.numerics.gauss-jackson
  "The Gauss-Jackson integrator for y'' = f(t, y): Stoermer-Cowell's
  method in summed form (Jackson, \"Note on the numerical integration of
  d2x/dt2 = f(x, t)\", MNRAS 84, 1924; Berry and Healy, \"Implementation of
  Gauss-Jackson Integration for Orbit Propagation\", J. Astronautical
  Sciences 52, 2004) -- the fixed-step workhorse of orbit propagation.

  Its formulas follow from the operator calculus: with D = -ln(1 - del)/h,
  del the backward difference,

    y_n  = h^2 [ U_n-1 + sum_m g_m+2 del^m f_n ],   del^2/ln^2(1 - del) = sum g_j del^j
    y'_n = h   [ T_n   + sum_m v_m+1 del^m f_n ],   del/(-ln(1 - del)) = sum v_j del^j

  T and U the running first and second sums of f -- the \"summed form\",
  which carries the integration's constants in the sums rather than in
  differences of nearly equal positions, and so keeps round-off from
  accumulating over long arcs. Each step predicts f one step ahead by
  extrapolating its differences, corrects, evaluates and corrects again.
  The coefficients are computed exactly, as ratios, from the two series;
  the method is of order `order` + 2 in position (10 by default)."
  (:require [allgo.numerics.core :as core]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.rk :as rk]))

(defn- series-inverse
  "The power series 1/a(x), for a(0) = 1, to `n` terms."
  [a n]
  (reduce (fn [b k]
            (conj b (- (reduce + (map #(* (nth a %) (nth b (- k %))) (range 1 (inc k)))))))
          [1] (range 1 n)))

(defn- series-product [a b n]
  (mapv (fn [k] (reduce + (map #(* (nth a %) (nth b (- k %))) (range (inc k))))) (range n)))

(defn coefficients
  "`{:g :v}`, exact ratios: the series of del^2/ln^2(1 - del) and
  del/(-ln(1 - del)) to `n` terms."
  [n]
  (let [;; -ln(1 - x)/x = 1 + x/2 + x^2/3 + ...
        a (mapv #(/ 1 (inc %)) (range n))
        v (series-inverse a n)]
    {:g (series-product v v n) :v v}))

(defn- differences
  "del^0 .. del^q of the history `fs` (oldest first), at its newest point."
  [fs q]
  (loop [row (vec (reverse fs)) out [] k 0]
    (if (> k q)
      out
      (recur (mapv lin/sub (butlast row) (rest row)) (conj out (first row)) (inc k)))))

(defn integrate
  "The trajectory of y'' = (f t y) from `y0` and `v0` at `t0`, `steps` steps
  of `h`: `[[t y v] ...]`, the first `order` + 1 from a Dormand-Prince
  start at tight tolerance, the rest by Gauss-Jackson."
  ([f t0 y0 v0 h steps] (integrate f t0 y0 v0 h steps {}))
  ([f t0 y0 v0 h steps {:keys [order] :or {order 8}}]
   (let [{:keys [g v]} (coefficients (+ order 3))
         g (mapv double g) v (mapv double v)
         d (count y0)
         first-order (fn [t yv] (into (subvec yv d) (f t (subvec yv 0 d))))
         start (let [integ (rk/integrator rk/dopri54 first-order t0 (into (vec y0) v0) (* 0.1 h)
                                          {:tol-abs 1e-14 :tol-rel 1e-14})]
                 (reductions (fn [s k] (core/step-until s (+ t0 (* k h)))) integ (range 1 (inc order))))
         states (mapv (fn [s] [(:t s) (subvec (:y s) 0 d) (subvec (:y s) d)]) start)
         fs (mapv (fn [[t y]] (f t y)) states)
         combo (fn [coeffs offset diffs]
                 (reduce lin/add (map (fn [m dm] (lin/scale dm (coeffs (+ m offset)))) (range) diffs)))
         ;; the sums set so the formulas hold at the last starting point
         [tn yn vn] (peek states)
         diffs (differences fs order)
         T (lin/sub (lin/scale vn (/ 1.0 h)) (combo v 1 diffs))
         Um1 (lin/sub (lin/scale yn (/ 1.0 (* h h))) (combo g 2 diffs))]
     (loop [out states fs fs T T U (lin/add Um1 T) t tn k (count states)]
       (if (> k steps)
         out
         (let [t1 (+ t h)
               ;; predict f one step on by extrapolating its differences
               diffs (differences fs order)
               f-pred (reduce lin/add diffs)
               correct (fn [f1]
                         (let [fs1 (conj (subvec fs 1) f1)
                               d1 (differences fs1 order)
                               T1 (lin/add T f1)]
                           [(lin/scale (lin/add U (combo g 2 d1)) (* h h))
                            (lin/scale (lin/add T1 (combo v 1 d1)) h)
                            fs1 T1]))
               [y1] (correct f-pred)
               [y1 _] (correct (f t1 y1))
               f1 (f t1 y1)
               [y1 v1 fs1 T1] (correct f1)]
           (recur (conj out [t1 y1 v1]) fs1 T1 (lin/add U T1) t1 (inc k))))))))
