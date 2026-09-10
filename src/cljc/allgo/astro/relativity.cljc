(ns allgo.astro.relativity
  "The post-Newtonian correction to two-body motion (Montenbruck & Gill 3.7).

  General relativity's leading contribution to an orbit, of order 1e-8 m/s^2
  in low Earth orbit -- smaller than everything else in the chapter, and
  still not negligible. It is the same term that advances Mercury's
  perihelion, and it does the same thing here: rather than perturbing the
  orbit at random it rotates it steadily in its own plane, so a correction
  far below the noise on any single revolution accumulates into a
  measurable secular drift.

  Written in the Schwarzschild form, which is what M&G gives and what
  applies to a satellite about a non-rotating central mass. The frame
  dragging and geodetic precession terms are smaller again and left out."
  (:require [allgo.astro.constants :as c]
            [clojure.math :as math]))

(defn acceleration
  "Relativistic correction at position `r` with velocity `v`, km/s^2.

    a = GM/(c^2 r^3) * ( (4GM/r - v.v) r + 4 (r.v) v )

  The first term always points inward or outward along the radius; the
  second vanishes on a circular orbit, where r and v are perpendicular, and
  is what makes the effect depend on eccentricity."
  ([r v] (acceleration r v c/GM-earth))
  ([r v GM]
   (let [r2  (reduce + (map * r r))
         rm  (math/sqrt r2)
         v2  (reduce + (map * v v))
         rv  (reduce + (map * r v))
         k   (/ GM (* c/c-light c/c-light r2 rm))
         a   (- (/ (* 4.0 GM) rm) v2)
         b   (* 4.0 rv)]
     (mapv (fn [ri vi] (* k (+ (* a ri) (* b vi)))) r v))))

(defn perihelion-advance
  "Secular advance of perihelion per revolution, radians (M&G eq. 3.148).

  The classical result: 6 pi GM / (c^2 a (1 - e^2)). Applied to Mercury
  about the Sun it gives the 43 arcseconds a century that the correction was
  first invented to explain, which makes it a convenient check that the
  constant in front is right."
  [GM a e]
  (/ (* 6.0 math/PI GM)
     (* c/c-light c/c-light a (- 1.0 (* e e)))))
