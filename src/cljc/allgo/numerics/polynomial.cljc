(ns allgo.numerics.polynomial
  "Polynomials: evaluation, and the roots of quadratics, cubics and
  quartics in closed form (the classical solutions -- Viete's
  trigonometric one for a cubic with three real roots, Cardano's
  otherwise, Ferrari's reduction of the quartic to a resolvent cubic and
  two quadratics), each root then polished by Newton's method on the
  original polynomial.

  Coefficients are given highest power first, `[a b c]` for a x^2 + b x
  + c. Roots are complex, `[re im]`, real ones with im 0; `real-roots`
  keeps the real ones."
  (:require [clojure.math :as math]))

;; ------------------------------------------------------ complex numbers

(defn- c+ [[a b] [c d]] [(+ a c) (+ b d)])
(defn- c- [[a b] [c d]] [(- a c) (- b d)])
(defn- c* [[a b] [c d]] [(- (* a c) (* b d)) (+ (* a d) (* b c))])
(defn- c-div [[a b] [c d]] (let [m (+ (* c c) (* d d))] [(/ (+ (* a c) (* b d)) m) (/ (- (* b c) (* a d)) m)]))
(defn- c-scale [[a b] s] [(* a s) (* b s)])

(defn- c-sqrt
  "The principal square root."
  [[a b]]
  (let [m (math/hypot a b)]
    (if (zero? m)
      [0.0 0.0]
      (let [re (math/sqrt (* 0.5 (+ m (abs a))))]
        (if (>= a 0.0)
          [re (/ b (* 2.0 re))]
          [(/ (abs b) (* 2.0 re)) (if (neg? b) (- re) re)])))))

;; ------------------------------------------------------------ evaluation

(defn horner
  "The polynomial with coefficients `coeffs`, highest power first, at
  `x`."
  [coeffs x]
  (reduce (fn [acc c] (+ (* acc x) c)) 0.0 coeffs))

(defn- c-horner
  "The polynomial and its derivative at complex `z`."
  [coeffs z]
  (reduce (fn [[p dp] c] [(c+ (c* p z) [c 0.0]) (c+ (c* dp z) p)]) [[0.0 0.0] [0.0 0.0]] coeffs))

(defn- polish
  "Newton's method on the polynomial `coeffs` from root `z`, a few steps,
  keeping a step only while it lowers the residual."
  [coeffs z]
  (loop [z z k 0]
    (let [[p dp] (c-horner coeffs z)
          m (math/hypot (first dp) (second dp))]
      (if (or (> k 4) (zero? m))
        z
        (let [z' (c- z (c-div p dp))
              [p'] (c-horner coeffs z')]
          (if (< (math/hypot (first p') (second p')) (math/hypot (first p) (second p)))
            (recur z' (inc k))
            z))))))

(defn- clean
  "Snap a root's negligible imaginary part to 0."
  [[re im]]
  [re (if (<= (abs im) (* 1e-12 (max 1.0 (abs re)))) 0.0 im)])

;; ---------------------------------------------------------------- roots

(defn quadratic
  "The two roots of a x^2 + b x + c, in the form that avoids cancellation
  (q = -(b + sign(b) sqrt(b^2 - 4ac))/2, roots q/a and c/q)."
  [a b c]
  (let [disc (- (* b b) (* 4.0 a c))]
    (if (>= disc 0.0)
      (let [q (* -0.5 (+ b (* (if (neg? b) -1.0 1.0) (math/sqrt disc))))]
        (if (zero? q)
          [[0.0 0.0] [0.0 0.0]]
          [[(/ q a) 0.0] [(/ c q) 0.0]]))
      (let [re (/ (- b) (* 2.0 a)) im (/ (math/sqrt (- disc)) (* 2.0 a))]
        [[re (abs im)] [re (- (abs im))]]))))

(defn- depressed-cubic-roots
  "Roots of x^3 + a x^2 + b x + c: with Q = (a^2 - 3b)/9 and R = (2a^3 -
  9ab + 27c)/54, three real roots when R^2 < Q^3 (Viete), one otherwise
  (Cardano)."
  [a b c]
  (let [q (/ (- (* a a) (* 3.0 b)) 9.0)
        r (/ (+ (* 2.0 a a a) (* -9.0 a b) (* 27.0 c)) 54.0)
        shift (/ a 3.0)]
    (if (< (* r r) (* q q q))
      (let [theta (math/acos (/ r (math/sqrt (* q q q))))
            s (* -2.0 (math/sqrt q))]
        (mapv (fn [k] [(- (* s (math/cos (/ (+ theta (* 2.0 math/PI k)) 3.0))) shift) 0.0]) [0 1 -1]))
      (let [big-a (* (if (pos? r) -1.0 1.0) (math/cbrt (+ (abs r) (math/sqrt (- (* r r) (* q q q))))))
            big-b (if (zero? big-a) 0.0 (/ q big-a))
            re (- (* -0.5 (+ big-a big-b)) shift)
            im (* 0.5 (math/sqrt 3.0) (- big-a big-b))]
        [[(- (+ big-a big-b) shift) 0.0] [re (abs im)] [re (- (abs im))]]))))

(defn cubic
  "The three roots of a x^3 + b x^2 + c x + d."
  [a b c d]
  (let [coeffs [a b c d]]
    (mapv #(clean (polish coeffs %)) (depressed-cubic-roots (/ b a) (/ c a) (/ d a)))))

(defn quartic
  "The four roots of a x^4 + b x^3 + c x^2 + d x + e, by Ferrari: x = y -
  b/4a turns it into y^4 + p y^2 + q y + r, which with a root m of the
  resolvent cubic m^3 + p m^2 + (p^2/4 - r) m - q^2/8 splits into y^2 +/-
  sqrt(2m) y + p/2 + m -/+ q/(2 sqrt(2m)); a biquadratic, q = 0, directly."
  [a b c d e]
  (let [coeffs [a b c d e]
        [b c d e] [(/ b a) (/ c a) (/ d a) (/ e a)]
        shift (/ b 4.0)
        p (- c (* 0.375 b b))
        q (+ d (* -0.5 b c) (* 0.125 b b b))
        r (+ e (* -0.25 b d) (* 0.0625 b b c) (* (/ -3.0 256.0) b b b b))
        scale (max 1.0 (abs p) (math/sqrt (abs r)))
        ys (if (<= (abs q) (* 1e-14 scale scale scale))
             ;; y^4 + p y^2 + r = 0: y = +/- sqrt of each root of z^2 + p z + r
             (let [disc (c-sqrt [(- (* p p) (* 4.0 r)) 0.0])]
               (mapcat (fn [z] (let [s (c-sqrt z)] [s (c-scale s -1.0)]))
                       [(c-scale (c+ [(- p) 0.0] disc) 0.5) (c-scale (c- [(- p) 0.0] disc) 0.5)]))
             (let [m (->> (depressed-cubic-roots p (- (* 0.25 p p) r) (* -0.125 q q))
                          (filter #(zero? (second %)))
                          (map first)
                          (apply max))
                   s (math/sqrt (* 2.0 m))
                   quad (fn [lin const]
                          (let [disc (c-sqrt [(- (* lin lin) (* 4.0 const)) 0.0])]
                            [(c-scale (c+ [(- lin) 0.0] disc) 0.5) (c-scale (c- [(- lin) 0.0] disc) 0.5)]))]
               (concat (quad s (- (+ (* 0.5 p) m) (/ q (* 2.0 s))))
                       (quad (- s) (+ (* 0.5 p) m (/ q (* 2.0 s)))))))]
    (mapv #(clean (polish coeffs (c- % [shift 0.0]))) ys)))

(defn real-roots
  "The real ones among `roots`, ascending."
  [roots]
  (vec (sort (map first (filter #(zero? (second %)) roots)))))
