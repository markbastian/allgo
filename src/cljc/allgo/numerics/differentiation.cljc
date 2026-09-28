(ns allgo.numerics.differentiation
  "Numerical derivatives: the Jacobian of a vector function by central
  differences, optionally refined by Richardson extrapolation.

  A central difference (f(x + h) - f(x - h)) / 2h is wrong by a term in
  h^2; taking it at h and h/2 and combining them as (4 D(h/2) - D(h)) / 3
  cancels that term and leaves one in h^4, which at a sensible step is ten
  digits for any smooth function."
  (:require [allgo.math :as am]
            [allgo.numerics.linear :as lin]))

(defn jacobian
  "The Jacobian of `f`, a function from a vector to a vector, at `x`, one
  column per component of x. Options:

    :steps        the step h_j for each component (default 1e-5 of the
                  component, or 1e-5 where it is smaller than 1)
    :angles       the output indices to difference modulo 2 pi, so a pair
                  straddling the seam is not read as a jump
    :richardson?  refine by Richardson extrapolation (default true); false
                  gives the plain central difference"
  ([f x] (jacobian f x {}))
  ([f x {:keys [steps angles richardson?] :or {angles #{} richardson? true}}]
   (let [x (vec x)
         steps (or steps (mapv #(* 1e-5 (max 1.0 (abs %))) x))
         diff (fn [a b] (vec (map-indexed (fn [k [p q]] (let [d (- p q)] (if (angles k) (am/wrap-angle d) d)))
                                          (map vector a b))))
         column (fn [j]
                  (let [d (fn [h] (lin/scale (diff (f (update x j + h)) (f (update x j - h))) (/ 1.0 (* 2.0 h))))
                        h (steps j)]
                    (if richardson?
                      (lin/scale (lin/sub (lin/scale (d (* 0.5 h)) 4.0) (d h)) (/ 1.0 3.0))
                      (d h))))]
     (lin/transpose (mapv column (range (count x)))))))
