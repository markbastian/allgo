(ns allgo.numerics.algebraic
  "Polynomials with integer coefficients, worked exactly, and the real
  numbers that are their roots -- what cylindrical algebraic
  decomposition stands on (LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, section 6.4.2): a root found only to a
  tolerance can be confused with a neighbor, or a double root with two
  close ones, and the decomposition built on it is wrong.

  Coefficients are doubles holding integers, highest power first as in
  `allgo.numerics.polynomial`; the zero polynomial is `[]`. A double holds
  every integer up to 2^53 exactly, so sums, products and the divisions
  here that are known to come out even are exact, and an operation whose
  result outgrows that throws rather than round.

  Greatest common divisors by the primitive remainder sequence, which
  keeps every coefficient an integer and divides out the content at each
  step of each division to keep them small; square-free parts from them;
  Sturm sequences, which count the distinct real roots in an interval,
  for isolating each root in an interval of its own, where bisection then
  narrows it; and the resultant of two polynomials in y whose
  coefficients are polynomials in x -- the determinant of their Sylvester
  matrix, by Bareiss's fraction-free elimination, whose divisions are
  exact -- a polynomial in x vanishing where the two share a root."
  (:require [allgo.numerics.polynomial :as poly]
            [allgo.numerics.roots :as roots]
            [allgo.search :as search]
            [clojure.math :as math]))

(def ^:private exact-limit 9007199254740992.0)

(defn- check [p]
  (when (some #(> (abs %) exact-limit) p)
    (throw (ex-info "A coefficient outgrew the integers a double holds exactly" {:p p})))
  p)

(defn trim
  "`p` without its leading zeros."
  [p]
  (vec (drop-while zero? p)))

(defn degree
  "The degree of `p`, -1 for the zero polynomial."
  [p]
  (dec (count (trim p))))

(defn integral?
  "Whether every coefficient of `p` is an integer."
  [p]
  (every? #(== % (math/floor %)) p))

(defn add
  "`p` plus `q`."
  [p q]
  (let [n (max (count p) (count q))
        pad (fn [v] (into (vec (repeat (- n (count v)) 0.0)) v))]
    (check (trim (mapv + (pad p) (pad q))))))

(defn scale
  "`p` times the number `k`."
  [p k]
  (check (trim (mapv #(* k %) p))))

(defn sub
  "`p` minus `q`."
  [p q]
  (add p (scale q -1.0)))

(defn mul
  "`p` times `q`."
  [p q]
  (if (or (empty? p) (empty? q))
    []
    (check (trim (reduce (fn [acc [i a]]
                           (reduce (fn [acc [j b]] (update acc (+ i j) + (* a b))) acc (map-indexed vector q)))
                         (vec (repeat (+ (count p) (count q) -1) 0.0))
                         (map-indexed vector p))))))

(defn derivative
  "The derivative of `p`."
  [p]
  (let [p (trim p) n (degree p)]
    (trim (vec (map-indexed (fn [i c] (* (- n i) c)) (butlast p))))))

(defn quotient
  "`p` divided by `d`, which must divide it: over the integers, when `d`
  is primitive, the quotient is integral, so long division meets only
  quotients of integers that come out even."
  [p d]
  (if (empty? (trim p))
    []
    (let [{:keys [quotient remainder]} (poly/long-divide (trim p) (trim d))]
      (when (seq (trim remainder))
        (throw (ex-info "The division does not come out even" {:p p :d d :remainder remainder})))
      (check (trim quotient)))))

(defn- prem-step
  "One step of pseudo-division by `b`: the leading term of the remainder
  cancelled after multiplying it by b's leading coefficient."
  [b {:keys [r e] :as state}]
  (when (and (seq r) (>= (degree r) (degree b)))
    (let [shift (- (degree r) (degree b))
          lead (first r)]
      (assoc state
             :r (sub (scale r (first b)) (scale (into b (repeat shift 0.0)) lead))
             :e (dec e)))))

(defn prem
  "The pseudo-remainder of `a` by `b`: the remainder of lc(b)^(deg a -
  deg b + 1) a by b, which is integral when both are."
  [a b]
  (let [a (trim a) b (trim b)
        {:keys [r e]} (last (search/step-seq (partial prem-step b) {:r a :e (inc (- (degree a) (degree b)))}))]
    (scale r (reduce * 1.0 (repeat (max 0 e) (first b))))))

(defn- euclid-step [[a b]]
  (when-not (zero? b) [b (rem a b)]))

(defn- integer-gcd [a b]
  (abs (first (last (search/step-seq euclid-step [a b])))))

(defn content
  "The greatest common divisor of `p`'s coefficients."
  [p]
  (reduce integer-gcd 0.0 p))

(defn primitive
  "`p` divided by its content, its leading coefficient made positive."
  [p]
  (let [p (trim p)]
    (if (empty? p)
      p
      (let [c (content p)]
        (scale p (/ (if (neg? (first p)) -1.0 1.0) c))))))

(defn- remainder-step
  "One step of long division by `b` in integers: the leading term of the
  remainder cancelled after multiplying it by |lc(b)|, then its content
  divided out -- positive factors, so the result stays a positive
  multiple of the remainder over the rationals, and small."
  [b r]
  (when (and (seq r) (>= (degree r) (degree b)))
    (let [shift (- (degree r) (degree b))
          s (if (neg? (first b)) -1.0 1.0)
          r (sub (scale r (abs (first b))) (scale (into b (repeat shift 0.0)) (* s (first r))))]
      (if (seq r) (scale r (/ 1.0 (content r))) r))))

(defn- remainder
  "A positive integer multiple of the remainder of `a` by `b`."
  [a b]
  (let [b (trim b)]
    (last (search/step-seq (partial remainder-step b) (trim a)))))

(defn- gcd-step [[a b]]
  (when (seq b) [b (primitive (remainder a b))]))

(defn gcd
  "The greatest common divisor of `p` and `q`, primitive, by the
  primitive remainder sequence."
  [p q]
  (let [[a b] (sort-by degree > [(primitive p) (primitive q)])]
    (primitive (first (last (search/step-seq gcd-step [a b]))))))

(defn square-free
  "The product of `p`'s distinct irreducible factors, primitive: `p`
  divided by its greatest common divisor with its derivative, so each
  root once."
  [p]
  (let [p (primitive p)]
    (if (< (degree p) 1)
      p
      (primitive (quotient p (gcd p (derivative p)))))))

(defn- sturm-step
  "The next of a Sturm sequence: minus the remainder of a by b, up to a
  positive factor, which leaves the signs alone."
  [[a b]]
  (when (seq b)
    [b (scale (remainder a b) -1.0)]))

(defn sturm
  "The Sturm sequence of the square-free `p`: p, p', and each next one
  minus the remainder of the two before, to a constant."
  [p]
  (let [p (trim p)]
    (mapv first (search/step-seq sturm-step [p (derivative p)]))))

(defn sign-changes
  "The changes of sign along `chain` at `x`, zeros skipped."
  [chain x]
  (->> chain
       (map #(poly/horner % x))
       (remove zero?)
       (partition 2 1)
       (filter (fn [[a b]] (neg? (* a b))))
       count))

(defn count-roots
  "Sturm's theorem: the distinct real roots of the polynomial whose
  Sturm sequence is `chain` in (a, b], neither a root."
  [chain a b]
  (- (sign-changes chain a) (sign-changes chain b)))

(defn root-bound
  "Cauchy's bound: every root of `p` is less than this in magnitude."
  [p]
  (let [p (trim p)]
    (inc (reduce max 0.0 (map #(abs (/ % (first p))) (rest p))))))

(defn- isolate-step
  "Split the next interval holding more than one root in two, keep one
  holding one, drop one holding none."
  [p chain {:keys [pending] :as state}]
  (when-let [[a b] (peek pending)]
    (let [state (update state :pending pop)
          n (count-roots chain a b)]
      (cond
        (zero? n) state
        (= 1 n) (update state :isolated conj [a b])
        :else (let [m (if (zero? (poly/horner p (* 0.5 (+ a b)))) (+ a (* 0.41 (- b a))) (* 0.5 (+ a b)))]
                (update state :pending conj [a m] [m b]))))))

(defn isolate
  "Intervals (a, b] each holding exactly one root of the square-free
  `p`, ascending."
  [p]
  (let [p (trim p)]
    (if (< (degree p) 1)
      []
      (let [chain (sturm p) r (root-bound p)]
        (->> (search/step-seq (partial isolate-step p chain) {:pending [[(- r) r]] :isolated []})
             last :isolated (sort-by first) vec)))))

(defn refine
  "The root of `p` in the interval (a, b], by bisection -- exactly, when
  the root is an integer, which rounding near it would otherwise blur."
  [p [a b]]
  (let [x (roots/bisect #(poly/horner p %) a b)
        k (math/rint x)]
    (if (and (< (abs (- k x)) 1e-9) (zero? (poly/horner p k))) (+ k 0.0) x)))

(defn real-roots
  "The distinct real roots of `p`, ascending: each isolated in an
  interval of its own by Sturm's sequence of p's square-free part, then
  narrowed by bisection."
  [p]
  (let [s (square-free p)]
    (mapv #(refine s %) (isolate s))))

;; --------------------------------------------- polynomials in y over x

(defn- bareiss-step
  "One step of Bareiss's elimination: below and right of pivot k, each
  entry m_ij m_kk - m_ik m_kj over the previous pivot, which divides it
  exactly."
  [{:keys [m k prev sign n] :as state}]
  (when (< k (dec n))
    (let [r (first (filter #(seq (get-in m [% k])) (range k n)))]
      (if (nil? r)
        (assoc state :k n :zero? true)
        (let [m (if (= r k) m (assoc m k (m r) r (m k)))
              pivot (get-in m [k k])
              m (reduce (fn [m [i j]]
                          (assoc-in m [i j] (quotient (sub (mul (get-in m [i j]) pivot)
                                                           (mul (get-in m [i k]) (get-in m [k j])))
                                                      prev)))
                        m
                        (for [i (range (inc k) n) j (range (inc k) n)] [i j]))]
          (assoc state :m m :k (inc k) :prev pivot :sign (if (= r k) sign (- sign))))))))

(defn determinant
  "The determinant of a square matrix of integer polynomials, by
  Bareiss's fraction-free elimination."
  [m]
  (let [n (count m)
        {:keys [m sign zero?]} (last (search/step-seq bareiss-step {:m m :k 0 :prev [1.0] :sign 1.0 :n n}))]
    (if zero? [] (scale (get-in m [(dec n) (dec n)]) sign))))

(defn resultant
  "The resultant in y of `f` and `g`, each a vector of polynomials in x,
  the coefficients of y's powers highest first: the determinant of their
  Sylvester matrix, a polynomial in x that vanishes wherever f and g
  share a root y, or both leading coefficients vanish."
  [f g]
  (let [m (dec (count f)) n (dec (count g)) size (+ m n)
        row (fn [coeffs shift] (vec (concat (repeat shift []) coeffs (repeat (- size shift (count coeffs)) []))))]
    (determinant (vec (concat (for [i (range n)] (row f i))
                              (for [i (range m)] (row g i)))))))
