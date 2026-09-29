(ns allgo.numerics.fit
  "Least-squares curve fitting as Meeus poses it (*Astronomical
  Algorithms*, chapter 4) -- a line, a parabola, or any combination of up
  to three given functions through a column of `[[x y] ...]` points,
  returning just the coefficients -- over the general machinery of
  `allgo.numerics.regression`, which also gives the fits' statistics. The
  iterations of his chapter 5 are in `allgo.numerics.roots`."
  (:require [allgo.numerics.regression :as reg]))

(defn- columns [points] [(mapv first points) (mapv second points)])

(defn linear
  "`[a b]` of the least-squares line y = a x + b."
  [points]
  (let [[xs ys] (columns points)
        {:keys [a0 a1]} (reg/linear xs ys)]
    [a1 a0]))

(defn correlation
  "The correlation coefficient: 1 or -1 for points exactly on a line, near
  zero when the line explains nothing."
  [points]
  (let [[xs ys] (columns points)] (:r (reg/linear xs ys))))

(defn quadratic
  "`[a b c]` of the least-squares parabola y = a x^2 + b x + c."
  [points]
  (let [[xs ys] (columns points)] (:coeffs (reg/polynomial xs ys 2))))

(defn functions-3
  "`[a b c]` of the least-squares fit y = a f0(x) + b f1(x) + c f2(x), for
  any three functions -- a sine and a cosine and a constant, say, for a
  periodic signal of known period."
  [points f0 f1 f2]
  (let [[xs ys] (columns points)] (:a (reg/general [f0 f1 f2] xs ys))))

(defn function-1
  "`a` of the least-squares fit y = a f(x)."
  [points f]
  (let [[xs ys] (columns points)] (first (:a (reg/general [f] xs ys)))))
