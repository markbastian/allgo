(ns allgo.numerics.differentiation
  "Numerical derivatives: the Jacobian of a vector function by central
  differences, optionally refined by Richardson extrapolation.

  A central difference (f(x + h) - f(x - h)) / 2h is wrong by a term in
  h^2; taking it at h and h/2 and combining them as (4 D(h/2) - D(h)) / 3
  cancels that term and leaves one in h^4, which at a sensible step is ten
  digits for any smooth function.

  And the rest (Chapra and Canale, *Numerical Methods for Engineers*,
  chapter 23): finite-difference formulas of any order and accuracy --
  forward, backward and centered, the book's tables generated rather than
  transcribed, by Fornberg's algorithm -- Richardson extrapolation,
  derivatives of unequally spaced data, derivatives of noisy data through
  a least-squares fit, and partial derivatives."
  (:require [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.polynomial :as poly]
            [allgo.numerics.regression :as reg]))

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

;; ------------------------------------------------ finite differences

(defn fd-weights
  "The weights that give the derivatives of orders 0 through `m` at `x0`
  from function values at the points `xs`, any spacing (Fornberg,
  \"Generation of finite difference formulas on arbitrarily spaced
  grids\", Math. Comp. 51, 1988): `(nth (fd-weights ...) k)` holds the
  weight of each point for the k-th derivative. On equal spacings these
  are the familiar formulas -- [-1/2 0 1/2]/h for a centered first
  derivative, and so on."
  [x0 xs m]
  (let [xs (vec xs) n (count xs)
        zero (vec (repeat (inc m) (vec (repeat n 0.0))))]
    (loop [i 1 c1 1.0 c4 (- (xs 0) x0) C (assoc-in zero [0 0] 1.0)]
      (if (= i n)
        C
        (let [c5 (- (xs i) x0)
              mn (min i m)
              [C c2]
              (reduce (fn [[C c2] j]
                        (let [c3 (- (xs i) (xs j))
                              c2 (* c2 c3)
                              C (if (= j (dec i))
                                  (reduce (fn [C k]
                                            (assoc-in C [k i]
                                                      (/ (* c1 (- (* k (get-in C [(dec k) (dec i)]))
                                                                  (* c4 (get-in C [k (dec i)]))))
                                                         c2)))
                                          (assoc-in C [0 i] (/ (* (- c1) c4 (get-in C [0 (dec i)])) c2))
                                          (range 1 (inc mn)))
                                  C)
                              C (reduce (fn [C k]
                                          (assoc-in C [k j] (/ (- (* c5 (get-in C [k j]))
                                                                  (* k (if (pos? k) (get-in C [(dec k) j]) 0.0)))
                                                               c3)))
                                        C (range mn -1 -1))]
                          [C c2]))
                      [C 1.0] (range i))]
          (recur (inc i) c2 c5 C))))))

(defn- math-pow [h k] (reduce * (repeat k h)))

(defn- stencil
  "The offsets, in steps, of the formula of `scheme` for the `order`-th
  derivative to error O(h^accuracy): the fewest points that give it."
  [scheme order accuracy]
  (case scheme
    :centered (let [k (+ (quot (+ order 1) 2) (quot accuracy 2) -1)]
                (range (- k) (inc k)))
    :forward (range 0 (+ order accuracy))
    :backward (range (- 1 (+ order accuracy)) 1)))

(defn derivative
  "The `:order`-th derivative (default 1) of `f` at `x` with step `h`, by
  the finite-difference formula of `:scheme` -- `:centered` (default),
  `:forward` or `:backward` -- accurate to O(h^`:accuracy`) (default 2;
  the book's high-accuracy formulas are accuracy 2 forward and backward,
  4 centered)."
  ([f x h] (derivative f x h {}))
  ([f x h {:keys [order scheme accuracy] :or {order 1 scheme :centered accuracy 2}}]
   (let [offsets (stencil scheme order accuracy)
         ws (nth (fd-weights 0.0 offsets order) order)]
     (/ (reduce + (map (fn [w k] (* w (f (+ x (* k h))))) ws offsets))
        (math-pow h order)))))

(defn richardson
  "Richardson extrapolation of an estimate `(D h)` whose error runs in
  h^p, h^(p+2), ... (p default 2): D(h/2) and D(h) combined as (2^p
  D(h/2) - D(h)) / (2^p - 1), repeated `levels` times (default 1) each
  cancelling the next term."
  ([D h] (richardson D h {}))
  ([D h {:keys [p levels] :or {p 2 levels 1}}]
   (let [base (mapv #(D (/ h (math-pow 2.0 %))) (range (inc levels)))]
     (loop [col base k 0]
       (if (= (count col) 1)
         (first col)
         (let [q (math-pow 2.0 (+ p (* 2 k)))]
           (recur (mapv #(/ (- (* q %2) %1) (- q 1.0)) col (rest col)) (inc k))))))))

(defn data-derivative
  "The derivative at `x` of the data `xs` `ys`, any spacing: the
  derivative of the parabola through the three points nearest x
  (Chapra and Canale 23.9), by `fd-weights` -- second-order accurate even
  at an end of the data. With `:points` more of them."
  ([xs ys x] (data-derivative xs ys x {}))
  ([xs ys x {:keys [points order] :or {points 3 order 1}}]
   (let [nearest (take points (sort-by #(abs (- (first %) x)) (map vector xs ys)))
         ws (nth (fd-weights x (map first nearest) order) order)]
     (reduce + (map * ws (map second nearest))))))

(defn smoothed-derivative
  "The derivative at `x` of noisy data `xs` `ys` taken from their
  least-squares polynomial of `degree` rather than from the data
  themselves (Chapra and Canale 23.4): differentiation amplifies noise,
  and fitting first smooths it away."
  [xs ys degree x]
  (let [{:keys [coeffs]} (reg/polynomial xs ys degree)]
    (second (poly/derivatives coeffs x 1))))

(defn partial-derivative
  "The partial derivative of `f`, a function of a vector, in component `i`
  at `x`, by a centered difference of step `h` (default 1e-6 of the
  component, at least 1e-6)."
  ([f x i] (partial-derivative f x i nil))
  ([f x i h]
   (let [x (mapv double x)
         h (or h (* 1e-6 (max 1.0 (abs (x i)))))]
     (/ (- (f (update x i + h)) (f (update x i - h))) (* 2.0 h)))))

(defn mixed-partial
  "The mixed partial derivative of `f` in components `i` and `j` at `x`,
  by the four-point centered formula."
  [f x i j]
  (let [x (mapv double x)
        hi (* 1e-4 (max 1.0 (abs (x i)))) hj (* 1e-4 (max 1.0 (abs (x j))))
        at (fn [a b] (f (-> x (update i + (* a hi)) (update j + (* b hj)))))]
    (/ (- (+ (at 1 1) (at -1 -1)) (at 1 -1) (at -1 1)) (* 4.0 hi hj))))
