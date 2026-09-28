(ns allgo.geometry.vec3
  "Three-component vectors as plain `[x y z]`.

  Small enough that the representation is the interesting decision: a
  Clojure vector allocates where a flat array would not, which is the
  wrong trade for thousands of particles and the right one for tens of
  rigid bodies. Simulations here that track many things keep flat arrays
  and index them -- `allgo.simulation.boids-flat`, `allgo.physics.flip` --
  and simulations that track few use this, because a body is read far more
  often in code than in a loop."
  (:require [allgo.math :as am]
            [clojure.math :as math]))

(def zero [0.0 0.0 0.0])

(defn add
  ([[ax ay az] [bx by bz]] [(+ ax bx) (+ ay by) (+ az bz)])
  ([a b & more] (reduce add (add a b) more)))

(defn sub [[ax ay az] [bx by bz]] [(- ax bx) (- ay by) (- az bz)])

(defn negate [[x y z]] [(- x) (- y) (- z)])

(defn scale [[x y z] s]
  (let [s (double s)] [(* x s) (* y s) (* z s)]))

(defn add-scaled
  "`a + b * s`, the step that shows up in every integrator."
  [a b s]
  (add a (scale b s)))

(defn mul
  "Component by component. What multiplying by a diagonal matrix comes to,
  which is how an inertia tensor is stored when the body's own axes are
  its principal ones."
  [[ax ay az] [bx by bz]]
  [(* ax bx) (* ay by) (* az bz)])

(defn dot ^double [[ax ay az] [bx by bz]]
  (+ (* ax bx) (* ay by) (* az bz)))

(defn cross [[ax ay az] [bx by bz]]
  [(- (* ay bz) (* az by))
   (- (* az bx) (* ax bz))
   (- (* ax by) (* ay bx))])

(defn length-squared ^double [v] (dot v v))

(defn length ^double [v] (math/sqrt (dot v v)))

(defn distance ^double [a b] (length (sub a b)))

(defn normalize
  "A unit vector in the same direction, or zero if there is no direction."
  [v]
  (let [l (length v)]
    (if (zero? l) zero (scale v (/ 1.0 l)))))

(defn lerp [a b t] (add-scaled a (sub b a) t))

(defn zero-vector? [v] (zero? (length-squared v)))

(defn finite?
  "Every component a real number. Worth asserting after a solve."
  [[x y z]]
  (and (am/finite? x) (am/finite? y) (am/finite? z)))

(defn angle
  "The angle between `a` and `b`, radians, 0 to pi: atan2 of |a x b| and
  a . b, which holds its digits at every angle -- acos of the cosine
  loses half of them near 0 and pi."
  ^double [a b]
  (math/atan2 (length (cross a b)) (dot a b)))
