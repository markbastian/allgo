(ns allgo.astro.srp
  "Solar radiation pressure, and the shadow that interrupts it
  (Montenbruck & Gill 3.4).

  Sunlight carries momentum, and a satellite is a sail whether or not it was
  meant to be. The acceleration is tiny -- around 1e-7 m/s^2 -- but unlike
  gravity it does not average away over an orbit: it is switched off for
  part of every revolution by the Earth's shadow, and that periodic
  interruption is what makes it accumulate.

  Which is why the shadow is modeled as a cone rather than a cylinder. The
  Sun is not a point, so its light does not stop abruptly at a boundary;
  there is a penumbra where it is partly occulted, and a satellite crossing
  it sees the force ramp rather than switch."
  (:require [allgo.astro.constants :as c]
            [clojure.math :as math]))

(defn- mag [v] (math/sqrt (reduce + (map * v v))))

(defn- circle-overlap
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

(defn shadow
  "Fraction of the Sun visible from `r`, given the Sun at `r-sun`: 1 in full
  sunlight, 0 in umbra, between the two in penumbra.

  Worked in apparent angular radii as seen from the satellite, which is the
  natural frame for an occultation and makes umbra, penumbra and the annular
  case fall out of the same comparison."
  ([r r-sun] (shadow r r-sun c/R-earth))
  ([r r-sun r-occulter]
   (let [to-sun   (mapv - r-sun r)
         d-sun    (mag to-sun)
         d-earth  (mag r)
         ;; apparent radii, and the apparent separation of their centers
         a        (math/asin (min 1.0 (/ c/R-sun d-sun)))
         b        (math/asin (min 1.0 (/ r-occulter d-earth)))
         cos-sep  (/ (reduce + (map * (mapv - r) to-sun)) (* d-earth d-sun))
         sep      (math/acos (max -1.0 (min 1.0 cos-sep)))]
     (- 1.0 (circle-overlap a b sep)))))

(defn acceleration
  "Radiation-pressure acceleration on a satellite at `r` with the Sun at
  `r-sun`, km/s^2.

  `area-to-mass` is in m^2/kg and `cr` is the reflectivity coefficient: 1 for
  a black body that absorbs everything, 2 for a flat mirror that returns the
  photons and so takes twice the momentum. Real satellites sit between.

  The 1/1000 converts the m/s^2 that pressure times area-over-mass yields
  into the km/s^2 everything else here uses."
  ([r r-sun area-to-mass cr] (acceleration r r-sun area-to-mass cr (shadow r r-sun)))
  ([r r-sun area-to-mass cr nu]
   (if (zero? nu)
     [0.0 0.0 0.0]
     (let [d  (mapv - r r-sun)
           dm (mag d)
           k  (* nu cr area-to-mass c/solar-pressure 1e-3
                 (/ (* c/AU c/AU) (* dm dm dm)))]
       (mapv #(* k %) d)))))

(defn umbra-length
  "How far behind the Earth its full shadow reaches, km.

  The Sun is larger than the Earth, so the shadow is a cone that closes to a
  point. Beyond about 1.4 million km there is no total eclipse to be had --
  which is why the Moon, at 384,000 km, can be totally eclipsed at all."
  [d-sun]
  (/ (* c/R-earth d-sun) (- c/R-sun c/R-earth)))
