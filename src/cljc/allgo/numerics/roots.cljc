(ns allgo.numerics.roots
  "Finding where a function vanishes, where a condition changes, and
  where a function peaks: bisection and fixed-point iteration (Meeus,
  *Astronomical Algorithms*, chapter 5), event location and golden-section
  search -- the slow, sure methods, for when a bracket is known."
  (:require [clojure.math :as math]))

(defn bisect
  "A root of `f` between `lo` and `hi`, where it changes sign, by halving
  the interval to the last bit.

  Slow -- one bit per evaluation -- but it cannot fail where a sign change
  is known, which is why Meeus recommends it for Kepler's equation at
  eccentricities near one, where Newton's method may not converge."
  [f lo hi]
  (let [flo (f lo)]
    (loop [lo lo hi hi flo flo i 0]
      (let [mid (* 0.5 (+ lo hi))
            fm  (f mid)]
        (cond
          (or (zero? fm) (>= i 64) (= mid lo) (= mid hi)) mid
          (= (neg? fm) (neg? flo)) (recur mid hi fm (inc i))
          :else (recur lo mid flo (inc i)))))))

(defn fixed-point
  "Iterate x <- (better x) from `start` until successive values differ by
  less than `tol`, returning the last. nil after `max-iterations` without
  settling.

  This is how most of astronomy's transcendental equations are solved:
  Kepler's equation, the time of a rising, the light time to a planet --
  each rearranged so the unknown appears once on the left, and fed back
  into itself. It converges whenever the rearranged right side changes
  more slowly than its argument."
  ([better start] (fixed-point better start 1e-15 100))
  ([better start tol max-iterations]
   (loop [x start i 0]
     (when (< i max-iterations)
       (let [n (better x)]
         (if (<= (abs (- n x)) (* tol (max 1.0 (abs n))))
           n
           (recur n (inc i))))))))

(defn crossing
  "The instant, to within `tol`, between `a` (where `(pred a)`) and `b`
  (where not) that `pred` changes."
  [pred a b tol]
  (loop [a a b b]
    (if (<= (abs (- b a)) tol)
      (* 0.5 (+ a b))
      (let [m (* 0.5 (+ a b))]
        (if (pred m) (recur m b) (recur a m))))))

(defn transitions
  "Where `pred` changes between `t0` and `t1`, sampling every `step` and
  refining each change to `tol`: `[[t from-value] ...]`. `step` must be
  shorter than the briefest interval worth finding."
  [pred t0 t1 step tol]
  (let [ts (concat (range t0 t1 step) [t1])]
    (->> (map vector ts (rest ts))
         (keep (fn [[a b]]
                 (let [pa (pred a) pb (pred b)]
                   (when (not= pa pb)
                     [(crossing #(= pa (pred %)) a b tol) pa])))))))

(def ^:private golden (/ (- (math/sqrt 5.0) 1.0) 2.0))

(defn minimize
  "The `x` between `a` and `b` where `f`, which has a single minimum there,
  is least, by golden-section search: until the bracket is narrower than
  `:tol`, or after `:max-iter` narrowings."
  ([f a b] (minimize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (loop [a a b b i 0]
     (if (or (> i max-iter) (< (abs (- b a)) tol))
       (* 0.5 (+ a b))
       (let [x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
         (if (< (f x1) (f x2)) (recur a x2 (inc i)) (recur x1 b (inc i))))))))

(defn maximize
  "The `x` between `a` and `b` where `f`, which has a single maximum there,
  is greatest, by golden-section search; options as `minimize`'s."
  ([f a b] (maximize f a b {}))
  ([f a b {:keys [tol max-iter] :or {tol 1e-10 max-iter 10000}}]
   (loop [a a b b i 0]
     (if (or (> i max-iter) (< (abs (- b a)) tol))
       (* 0.5 (+ a b))
       (let [x1 (- b (* golden (- b a))) x2 (+ a (* golden (- b a)))]
         (if (< (f x1) (f x2)) (recur x1 b (inc i)) (recur a x2 (inc i))))))))
