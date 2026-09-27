(ns allgo.numerics.fit
  "Least-squares curve fitting and the iterations that solve an equation
  by repeated improvement (Meeus, *Astronomical Algorithms*, chapters 4
  and 5).

  The fits are the closed forms: the normal equations for a line, a
  parabola, or any combination of up to three given functions, written out
  and solved by Cramer's rule. For more unknowns or weighted data the
  general machinery is `allgo.astro.estimation`; these are what an
  observer reducing a column of measurements actually reaches for.

  Points are `[[x y] ...]` throughout."
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------- chapter 4

(defn linear
  "`[a b]` of the least-squares line y = a x + b."
  [points]
  (let [n (count points)
        [sx sy sxx sxy] (reduce (fn [[sx sy sxx sxy] [x y]]
                                  [(+ sx x) (+ sy y) (+ sxx (* x x)) (+ sxy (* x y))])
                                [0.0 0.0 0.0 0.0] points)
        d (- (* n sxx) (* sx sx))]
    [(/ (- (* n sxy) (* sx sy)) d)
     (/ (- (* sy sxx) (* sx sxy)) d)]))

(defn correlation
  "The correlation coefficient: 1 or -1 for points exactly on a line, near
  zero when the line explains nothing."
  [points]
  (let [n (count points)
        [sx sy sxx syy sxy]
        (reduce (fn [[sx sy sxx syy sxy] [x y]]
                  [(+ sx x) (+ sy y) (+ sxx (* x x)) (+ syy (* y y)) (+ sxy (* x y))])
                [0.0 0.0 0.0 0.0 0.0] points)]
    (/ (- (* n sxy) (* sx sy))
       (* (math/sqrt (- (* n sxx) (* sx sx)))
          (math/sqrt (- (* n syy) (* sy sy)))))))

(defn quadratic
  "`[a b c]` of the least-squares parabola y = a x^2 + b x + c."
  [points]
  (let [N (double (count points))
        [P Q R S T U V]
        (reduce (fn [[P Q R S T U V] [x y]]
                  (let [x2 (* x x)]
                    [(+ P x) (+ Q x2) (+ R (* x x2)) (+ S (* x2 x2))
                     (+ T y) (+ U (* x y)) (+ V (* x2 y))]))
                [0.0 0.0 0.0 0.0 0.0 0.0 0.0] points)
        D (- (+ (* N Q S) (* 2.0 P Q R)) (* Q Q Q) (* P P S) (* N R R))]
    [(/ (- (+ (* N Q V) (* P R T) (* P Q U)) (* Q Q T) (* P P V) (* N R U)) D)
     (/ (- (+ (* N S U) (* P Q V) (* Q R T)) (* Q Q U) (* P S T) (* N R V)) D)
     (/ (- (+ (* Q S T) (* Q R U) (* P R V)) (* Q Q V) (* P S U) (* R R T)) D)]))

(defn functions-3
  "`[a b c]` of the least-squares fit y = a f0(x) + b f1(x) + c f2(x), for
  any three functions -- a sine and a cosine and a constant, say, for a
  periodic signal of known period."
  [points f0 f1 f2]
  (let [[M P Q R S T U V W]
        (reduce (fn [[M P Q R S T U V W] [x y]]
                  (let [y0 (f0 x) y1 (f1 x) y2 (f2 x)]
                    [(+ M (* y0 y0)) (+ P (* y0 y1)) (+ Q (* y0 y2))
                     (+ R (* y1 y1)) (+ S (* y1 y2)) (+ T (* y2 y2))
                     (+ U (* y y0)) (+ V (* y y1)) (+ W (* y y2))]))
                (repeat 9 0.0) points)
        D (- (+ (* M R T) (* 2.0 P Q S)) (* M S S) (* R Q Q) (* T P P))]
    [(/ (+ (* U (- (* R T) (* S S))) (* V (- (* Q S) (* P T))) (* W (- (* P S) (* Q R)))) D)
     (/ (+ (* U (- (* S Q) (* P T))) (* V (- (* M T) (* Q Q))) (* W (- (* P Q) (* M S)))) D)
     (/ (+ (* U (- (* P S) (* R Q))) (* V (- (* P Q) (* M S))) (* W (- (* M R) (* P P)))) D)]))

(defn function-1
  "`a` of the least-squares fit y = a f(x)."
  [points f]
  (let [[syf sff] (reduce (fn [[syf sff] [x y]]
                            (let [fx (f x)] [(+ syf (* y fx)) (+ sff (* fx fx))]))
                          [0.0 0.0] points)]
    (/ syf sff)))

;; ---------------------------------------------------------------- chapter 5

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
