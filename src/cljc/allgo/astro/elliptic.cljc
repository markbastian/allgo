(ns allgo.astro.elliptic
  "Where a planet, asteroid or comet appears from the Earth (Meeus,
  *Astronomical Algorithms*, chapter 33).

  The geometry is a subtraction -- the body's heliocentric position less
  the Earth's -- but the body is not seen where it is: its light left it
  minutes or hours ago, so the position wanted is the one it had then. The
  light time depends on the distance, which depends on the position, and
  the loop is closed by iterating, which settles in two or three passes.
  Then the same corrections as a star's: aberration, the FK5 frame, and
  nutation.

  Times are MJD (TT), angles radians, distances AU."
  (:require [allgo.astro.apparent :as apparent]
            [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.frames :as frames]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.solar :as solar]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

(defn- spherical->rect [[l b r]]
  [(* r (math/cos b) (math/cos l)) (* r (math/cos b) (math/sin l)) (* r (math/sin b))])

(defn light-time
  "Days for light to cross `distance` AU."
  [distance]
  (* c/light-time-au distance))

(defn geocentric
  "Geometric geocentric `[lon lat distance light-time]` of a VSOP87
  `planet` at `mjd-tt`, ecliptic of date: where it was when the light now
  arriving left it, corrected neither for aberration nor to FK5."
  [planet mjd-tt]
  (let [earth (vsop87/rectangular :earth mjd-tt)
        at    (fn [tau] (mapv - (vsop87/rectangular planet (- mjd-tt tau)) earth))]
    (loop [tau 0.0 i 0]
      (let [[x y z :as d] (at tau)
            delta (v3/length d)
            tau'  (light-time delta)]
        (if (or (< (abs (- tau' tau)) 1e-9) (>= i 10))
          [(am/wrap-2pi (math/atan2 y x)) (math/atan2 z (math/hypot x y)) delta tau']
          (recur tau' (inc i)))))))

(defn apparent
  "Apparent `[ra dec distance]` of a VSOP87 `planet` at `mjd-tt`: light
  time, aberration, the FK5 correction and nutation, on the true equator
  of date -- the place an almanac prints."
  [planet mjd-tt]
  (let [[l b delta] (geocentric planet mjd-tt)
        [dl db] (apparent/ecliptic-aberration [l b] mjd-tt)
        [l b] (vsop87/->fk5 [(+ l dl) (+ b db)] mjd-tt)
        [dpsi] (frames/nutation-angles mjd-tt)
        [ra dec] (coord/ecliptic->equatorial [(+ l dpsi) b] (frames/true-obliquity mjd-tt))]
    [ra dec delta]))

;; ------------------------------------------------------ from the elements

(defn orbit
  "A function of MJD (TT) giving the heliocentric position, AU, in the
  mean equatorial frame of J2000, of a body with elements `{:a :e :i :raan
  :argp :perihelion}` referred to the J2000 ecliptic -- the time of
  perihelion `:perihelion` an MJD. The orbit is a fixed ellipse: fine for
  a comet over an apparition, not for a planet over decades."
  [{:keys [a e i raan argp perihelion]}]
  (let [n   (/ c/gaussian-k (* a (math/sqrt a)))
        eps eph/obliquity-J2000
        se (math/sin eps) ce (math/cos eps)
        so (math/sin raan) co (math/cos raan)
        si (math/sin i) ci (math/cos i)
        ;; Meeus 33.7: the orbit's orientation folded into three
        ;; amplitude-phase pairs, one per axis
        F co              P (- (* so ci))
        G (* so ce)       Q (- (* co ci ce) (* si se))
        H (* so se)       R (+ (* co ci se) (* si ce))
        [A B C]  [(math/atan2 F P) (math/atan2 G Q) (math/atan2 H R)]
        [a' b' c'] [(math/hypot F P) (math/hypot G Q) (math/hypot H R)]]
    (fn [mjd-tt]
      (let [M  (* n (- mjd-tt perihelion))
            E  (kepler/kepler-equation M e)
            nu (kepler/eccentric->true E e)
            r  (* a (- 1.0 (* e (math/cos E))))]
        [(* r a' (math/sin (+ A argp nu)))
         (* r b' (math/sin (+ B argp nu)))
         (* r c' (math/sin (+ C argp nu)))]))))

(defn astrometric
  "Astrometric `[ra dec distance elongation]` at `mjd-tt` of a body whose
  heliocentric J2000 equatorial position is `(position mjd)` -- see
  `orbit`. Astrometric means referred to the mean equinox of J2000 and
  corrected for light time only, which is how comet and asteroid
  ephemerides are published: it can be plotted straight onto a J2000 star
  chart. The elongation is the body's angular distance from the Sun."
  [position mjd-tt]
  (let [sun (solar/rectangular-J2000 mjd-tt)
        at  (fn [tau] (mapv + sun (position (- mjd-tt tau))))]
    (loop [tau 0.0 i 0]
      (let [[x y z :as d] (at tau)
            delta (v3/length d)
            tau' (light-time delta)]
        (if (or (< (abs (- tau' tau)) 1e-9) (>= i 10))
          [(am/wrap-2pi (math/atan2 y x)) (math/asin (/ z delta)) delta
           (math/acos (/ (reduce + (map * d sun)) (* delta (v3/length sun))))]
          (recur tau' (inc i)))))))

(defn heliocentric->geocentric
  "Geocentric ecliptic `[lon lat distance]` of a body at heliocentric
  ecliptic `[lon lat r]` when the Earth is at `earth` (its own heliocentric
  `[lon lat r]`) -- the subtraction, in spherical coordinates."
  [body earth]
  (let [[x y z] (mapv - (spherical->rect body) (spherical->rect earth))]
    [(am/wrap-2pi (math/atan2 y x)) (math/atan2 z (math/hypot x y)) (v3/length [x y z])]))
