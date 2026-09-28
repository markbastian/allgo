(ns allgo.numerics.special
  "Special functions: the modified Bessel functions of the first kind and
  the Legendre functions -- associated, and the derivatives of the
  polynomials that Pines's gravity formulation is written in.")

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
