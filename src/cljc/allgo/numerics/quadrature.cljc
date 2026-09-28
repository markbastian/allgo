(ns allgo.numerics.quadrature
  "Definite integrals by the composite Newton-Cotes rules.

  Simpson's rule fits a parabola through each pair of intervals and is
  exact for cubics, its error falling as h^4; the trapezoid rule is exact
  only for lines -- except over a whole period of a smooth periodic
  function, where its error falls faster than any power of h, which makes
  it the rule for orbit averages.")

(defn simpson
  "The integral of `f` from `a` to `b` by Simpson's rule on `n` intervals
  (made even if odd)."
  [f a b n]
  (let [n (if (odd? n) (inc n) n)
        h (/ (- b a) n)]
    (* (/ h 3.0)
       (reduce + (for [k (range (inc n))
                       :let [w (cond (or (zero? k) (= k n)) 1.0 (odd? k) 4.0 :else 2.0)]]
                   (* w (f (+ a (* k h)))))))))

(defn trapezoid
  "The integral of `f` from `a` to `b` by the trapezoid rule on `n`
  intervals."
  [f a b n]
  (let [h (/ (- b a) n)]
    (* h (reduce + (for [k (range (inc n))
                         :let [w (if (or (zero? k) (= k n)) 0.5 1.0)]]
                     (* w (f (+ a (* k h)))))))))
