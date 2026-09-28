(ns allgo.geometry.disk
  "Overlapping disks in a plane -- or small circles on the sky, apparent
  radii seen from a point, which is what an eclipse is."
  (:require [clojure.math :as math]))

(defn overlap-fraction
  "Fraction of a disk of apparent radius `a` hidden behind one of apparent
  radius `b`, their centers `sep` apart. All angles, in radians."
  [a b sep]
  (cond
    (>= sep (+ a b))  0.0                        ; clear of each other
    (<= sep (- b a))  1.0                        ; the occulter covers it entirely
    (<= sep (- a b))  (/ (* b b) (* a a))        ; the occulter sits wholly inside
    :else
    ;; Two overlapping disks: the shared area is a pair of circular segments.
    (let [x    (/ (+ (* sep sep) (* a a) (- (* b b))) (* 2.0 sep))
          y    (math/sqrt (max 0.0 (- (* a a) (* x x))))
          area (- (+ (* a a (math/acos (/ x a)))
                     (* b b (math/acos (/ (- sep x) b))))
                  (* sep y))]
      (/ area (* math/PI a a)))))
