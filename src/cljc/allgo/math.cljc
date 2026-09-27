(ns allgo.math
  "The scalar helpers `clojure.math` does not have, in one place.

  Every one of these had been written more than once: a clamp in two
  namespaces and inline in twenty more, a lerp in three, three different
  ways of bringing an angle back onto the circle, and two `finite?`. None
  of them is hard, which is how they multiply; each copy is a chance for
  the argument order or an edge case to drift from the others.

  `clojure.math` itself stays where it is -- `[clojure.math :as math]`
  reads on both platforms and its functions compile to direct calls -- and
  this is required alongside it as `am`. Everything here is hinted to
  primitive doubles, so the physics can call it in an inner loop without
  boxing."
  (:require [clojure.math :as math]))

(defn sq
  "`x` squared."
  ^double [^double x]
  (* x x))

(defn clamp
  "`x` confined to [`lo`, `hi`].

  In the order GLSL and Peachey's C macro use: the value first. The usual
  call is guarding an inverse cosine against a dot product a rounding
  error has pushed a hair past one: `(math/acos (am/clamp d -1.0 1.0))`."
  ^double [^double x ^double lo ^double hi]
  (max lo (min hi x)))

(defn lerp
  "The linear blend from `a` at `t` = 0 to `b` at `t` = 1."
  ^double [^double a ^double b ^double t]
  (+ a (* t (- b a))))

(def ^:const two-pi
  "A full turn. The literal, so it inlines on both platforms."
  6.283185307179586)

(defn wrap-angle
  "An angle brought into (-pi, pi]: the signed turn, for differences.

  Two angles that differ by a full turn are the same direction, and the
  difference between two headings only means something once it has been
  brought here -- the short way round."
  ^double [^double theta]
  (let [r (rem theta two-pi)]
    (cond (> r math/PI) (- r two-pi)
          (<= r (- math/PI)) (+ r two-pi)
          :else r)))

(defn wrap-2pi
  "An angle brought into [0, 2 pi): the unsigned turn, for longitudes and
  anomalies."
  ^double [^double theta]
  (let [r (rem theta two-pi)]
    (if (neg? r) (+ r two-pi) r)))

(defn finite?
  "Neither NaN nor an infinity. Worth asserting after a solve."
  [^double x]
  (not (or (NaN? x) (infinite? x))))
