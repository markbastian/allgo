(ns allgo.geometry.rotation
  "Rotation matrices about the coordinate axes, in both senses.

  A frame rotation R(t) applied to v gives v's components in axes turned
  by t -- the astrodynamical convention, in which a chain of frames is a
  product of these. A vector rotation, its transpose, turns v itself by t
  in fixed axes. They differ only in sign, which is why both have names
  here rather than one being assumed: `rx` `ry` `rz` turn the frame,
  `rotate-x` `rotate-y` `rotate-z` the vector."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn rx [t] (let [s (math/sin t) k (math/cos t)]
               [[1.0 0.0 0.0] [0.0 k s] [0.0 (- s) k]]))
(defn ry [t] (let [s (math/sin t) k (math/cos t)]
               [[k 0.0 (- s)] [0.0 1.0 0.0] [s 0.0 k]]))
(defn rz [t] (let [s (math/sin t) k (math/cos t)]
               [[k s 0.0] [(- s) k 0.0] [0.0 0.0 1.0]]))

(defn rotate-x [t] (let [s (math/sin t) k (math/cos t)]
                     [[1.0 0.0 0.0] [0.0 k (- s)] [0.0 s k]]))
(defn rotate-y [t] (let [s (math/sin t) k (math/cos t)]
                     [[k 0.0 s] [0.0 1.0 0.0] [(- s) 0.0 k]]))
(defn rotate-z [t] (let [s (math/sin t) k (math/cos t)]
                     [[k (- s) 0.0] [s k 0.0] [0.0 0.0 1.0]]))

(defn chain
  "Compose rotations left to right, so `(chain a b c)` applied to a vector
  does c first."
  [& ms]
  (reduce lin/mat-mul ms))
