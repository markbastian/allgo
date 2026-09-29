(ns allgo.numerics.polynomial
  "Polynomials: evaluation, and the roots of quadratics, cubics and
  quartics in closed form (the classical solutions -- Viete's
  trigonometric one for a cubic with three real roots, Cardano's
  otherwise, Ferrari's reduction of the quartic to a resolvent cubic and
  two quadratics), each root then polished by Newton's method on the
  original polynomial.

  For any degree (Chapra and Canale, *Numerical Methods for Engineers*,
  chapter 7): evaluation with derivatives, deflation and division,
  Muller's method, which finds complex roots from real guesses, and
  Bairstow's, which peels off quadratic factors until every root is found.

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

(defn horner-ascending
  "The polynomial c0 + c1 x + c2 x^2 + ... with `coeffs` lowest power
  first -- the order the almanacs print their series in -- at `x`."
  [coeffs x]
  (reduce (fn [acc c] (+ (* acc x) c)) 0.0 (rseq (vec coeffs))))

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

;; --------------------------------------------------- any degree

(defn derivatives
  "The polynomial and its first `n` derivatives at `x`, `[p p' p'' ...]`,
  by repeated synthetic division (Horner's scheme carried along)."
  [coeffs x n]
  (let [deg (dec (count coeffs))
        ;; d[k] accumulates p^(k)/k!
        d (reduce (fn [d c]
                    (let [d' (reduce (fn [d k] (assoc d k (+ (* (d k) x) (d (dec k)))))
                                     d (range (min n deg) 0 -1))]
                      (assoc d' 0 (+ (* (d' 0) x) c))))
                  (vec (repeat (inc n) 0.0))
                  coeffs)]
    (vec (map-indexed (fn [k v] (* v (reduce * (range 1 (inc k))))) d))))

(defn deflate
  "The polynomial divided by (x - t): `{:quotient :remainder}`, by
  synthetic division -- the remainder the polynomial's value at t, zero
  when t is a root, and the quotient the polynomial of the other roots."
  [coeffs t]
  (let [qs (reductions (fn [r c] (+ c (* r t))) coeffs)]
    {:quotient (vec (butlast qs)) :remainder (last qs)}))

(defn long-divide
  "The polynomial `a` divided by the polynomial `d`: `{:quotient
  :remainder}`, by long division."
  [a d]
  (let [n (dec (count a)) m (dec (count d))]
    (if (< n m)
      {:quotient [0.0] :remainder (vec a)}
      (loop [r (vec a) q []]
        (if (< (dec (count r)) m)
          {:quotient q :remainder (if (empty? r) [0.0] r)}
          (let [f (/ (first r) (first d))
                r' (mapv - r (concat (map #(* f %) d) (repeat 0.0)))]
            (recur (subvec r' 1) (conj q f))))))))

(defn- c-horner-value [coeffs z] (first (c-horner coeffs z)))

(defn muller
  "A root of the polynomial (or of any function `f` of a complex number,
  as `[re im]`) by Muller's method from three guesses: the parabola
  through the last three points, its nearer root the next estimate --
  complex when the parabola misses the axis, which is how it finds complex
  roots from real starts. Until the step is below `:tol` relative
  (default 1e-14); nil after `:max-iter` (100)."
  ([coeffs x0 x1 x2] (muller coeffs x0 x1 x2 {}))
  ([coeffs x0 x1 x2 {:keys [tol max-iter] :or {tol 1e-14 max-iter 100}}]
   (let [f (if (fn? coeffs) coeffs #(c-horner-value coeffs %))
         z (fn [x] (if (number? x) [(double x) 0.0] x))]
     (loop [x0 (z x0) x1 (z x1) x2 (z x2) i 0]
       (when (< i max-iter)
         (let [h0 (c- x1 x0) h1 (c- x2 x1)
               f0 (f x0) f1 (f x1) f2 (f x2)
               d0 (c-div (c- f1 f0) h0) d1 (c-div (c- f2 f1) h1)
               a (c-div (c- d1 d0) (c+ h1 h0))
               b (c+ (c* a h1) d1)
               c f2
               rad (c-sqrt (c- (c* b b) (c-scale (c* a c) 4.0)))
               [p m] [(c+ b rad) (c- b rad)]
               den (if (> (math/hypot (first p) (second p)) (math/hypot (first m) (second m))) p m)
               dx (if (and (zero? (first den)) (zero? (second den))) [0.0 0.0] (c-div (c-scale c -2.0) den))
               x3 (c+ x2 dx)]
           (if (<= (math/hypot (first dx) (second dx)) (* tol (max 1e-300 (math/hypot (first x3) (second x3)))))
             (clean x3)
             (recur x1 x2 x3 (inc i)))))))))

(defn- bairstow-factor
  "A quadratic factor x^2 - r x - s of the polynomial `a` (highest power
  first) from the guess r s, by Bairstow's Newton iteration on the
  remainder of division by it. `[r s]`, or nil if it does not converge."
  [a r s tol max-iter]
  (let [n (dec (count a))
        ;; the book's recurrences run lowest power first
        a (vec (reverse a))]
    (loop [r r s s i 0]
      (when (< i max-iter)
        (let [b (reduce (fn [b k]
                          (assoc b k (+ (a k) (* r (get b (inc k) 0.0)) (* s (get b (+ k 2) 0.0)))))
                        (vec (repeat (+ n 3) 0.0)) (range n -1 -1))
              c (reduce (fn [c k]
                          (assoc c k (+ (b k) (* r (get c (inc k) 0.0)) (* s (get c (+ k 2) 0.0)))))
                        (vec (repeat (+ n 3) 0.0)) (range n 0 -1))
              ;; c2 dr + c3 ds = -b1, c1 dr + c2 ds = -b0
              det (- (* (c 2) (c 2)) (* (c 3) (c 1)))]
          (when-not (zero? det)
            (let [dr (/ (- (* (- (b 1)) (c 2)) (* (- (b 0)) (c 3))) det)
                  ds (/ (- (* (- (b 0)) (c 2)) (* (- (b 1)) (c 1))) det)
                  r' (+ r dr) s' (+ s ds)]
              (if (and (<= (abs dr) (* tol (max 1.0 (abs r')))) (<= (abs ds) (* tol (max 1.0 (abs s')))))
                [r' s']
                (recur r' s' (inc i))))))))))

(defn bairstow
  "Every root of the polynomial by Bairstow's method: a quadratic factor
  x^2 - r x - s found by Newton's method on the remainder of dividing by
  it, its two roots taken, the polynomial deflated by it, and on until a
  quadratic or linear factor is left -- complex roots in conjugate pairs,
  in real arithmetic throughout. The roots are then polished by Newton's
  method on the original. `:r` and `:s` are the first guess (default
  -1 and -1, as good as any)."
  ([coeffs] (bairstow coeffs {}))
  ([coeffs {:keys [r s tol max-iter] :or {r -1.0 s -1.0 tol 1e-14 max-iter 500}}]
   (let [coeffs (mapv double coeffs)
         original coeffs]
     (loop [a coeffs out []]
       (let [n (dec (count a))]
         (cond
           (= n 0) (mapv #(clean (polish original %)) out)
           (= n 1) (recur [1.0] (conj out [(/ (- (a 1)) (a 0)) 0.0]))
           (= n 2) (recur [1.0] (into out (quadratic (a 0) (a 1) (a 2))))
           :else
           (let [[r' s'] (or (some (fn [[r s]] (bairstow-factor a r s tol max-iter))
                                   [[r s] [0.5 -0.5] [1.3 0.7] [-2.1 -1.9] [0.0 1.0]])
                             (throw (ex-info "Bairstow did not converge" {:coeffs a})))
                 {:keys [quotient]} (long-divide a [1.0 (- r') (- s')])]
             (recur (mapv #(/ % (first quotient)) quotient)
                    (into out (quadratic 1.0 (- r') (- s')))))))))))
