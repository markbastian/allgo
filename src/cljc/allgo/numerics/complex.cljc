(ns allgo.numerics.complex
  "Complex arithmetic on `[re im]` pairs -- the few operations the
  polynomial root finders, the Fourier transforms and the Jacobia-Roberts
  atmosphere's partial fractions need, in one place."
  (:refer-clojure :exclude [+ - * /])
  (:require [clojure.core :as core]
            [clojure.math :as math]))

(defn complex
  "`x` as a complex number: a real number becomes `[x 0]`, a pair stays."
  [x]
  (if (number? x) [(double x) 0.0] x))

(defn + [[a b] [c d]] [(core/+ a c) (core/+ b d)])
(defn - [[a b] [c d]] [(core/- a c) (core/- b d)])
(defn * [[a b] [c d]] [(core/- (core/* a c) (core/* b d)) (core/+ (core/* a d) (core/* b c))])

(defn /
  "The quotient, the denominator's modulus squared divided out."
  [[a b] [c d]]
  (let [m (core/+ (core/* c c) (core/* d d))]
    [(core// (core/+ (core/* a c) (core/* b d)) m) (core// (core/- (core/* b c) (core/* a d)) m)]))

(defn scale [[a b] s] [(core/* a s) (core/* b s)])

(defn modulus [[a b]] (math/hypot a b))

(defn sqrt
  "The principal square root."
  [[a b]]
  (let [m (math/hypot a b)]
    (if (zero? m)
      [0.0 0.0]
      (let [re (math/sqrt (core/* 0.5 (core/+ m (abs a))))]
        (if (>= a 0.0)
          [re (core// b (core/* 2.0 re))]
          [(core// (abs b) (core/* 2.0 re)) (core/* (if (neg? b) -1.0 1.0) re)])))))

(defn expi
  "e^(i theta), the unit complex number at angle `theta`."
  [theta]
  [(math/cos theta) (math/sin theta)])
