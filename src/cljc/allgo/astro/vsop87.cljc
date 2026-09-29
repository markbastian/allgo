(ns allgo.astro.vsop87
  "Heliocentric positions of the planets from the VSOP87 theory of
  Bretagnon and Francou (Meeus, *Astronomical Algorithms*, chapter 32 and
  appendix III).

  VSOP87 writes each coordinate as a Poisson series: sums of A cos(B + C
  tau), each sum multiplied by a power of tau, the time in Julian
  millennia from J2000. The full theory runs to some 35,000 terms and
  matches the numerical ephemeris it was fitted to within a milliarcsecond
  over four thousand years. Meeus prints a truncation of version D --
  heliocentric ecliptic longitude, latitude and radius, referred to the
  ecliptic and equinox of date -- that is good to about an arcsecond.

  The truncation here is made the same way and to about the same size,
  from the original files (`allgo.astro.vsop87-data`): every term of the
  longitude and latitude series at least 2e-7 radians, and of the radius
  series at least 2e-7 times the planet's mean distance, which keeps 2802
  terms. Checked against the full series from 2000 BC to AD 3000, the worst
  errors are 1.7 arcseconds in longitude (Jupiter) and 1.8 in latitude
  (Saturn), and under 1.5 arcseconds of the radius seen from the Sun.

  Coordinates are in VSOP87's dynamical frame, which differs from the FK5
  catalog frame by a small rotation; `->fk5` removes it, and matters only
  at the level of a tenth of an arcsecond."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.vsop87-data :as data]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def planets
  "The planets VSOP87 covers, from the Sun out. The Earth's entry is the
  Earth's center, not the Earth-Moon barycenter."
  [:mercury :venus :earth :mars :jupiter :saturn :uranus :neptune])

(defn millennia
  "Julian millennia of TT since J2000.0, VSOP87's tau."
  [mjd-tt]
  (/ (- mjd-tt c/mjd-J2000) 365250.0))

(defn- evaluate [series tau]
  (reduce (fn [acc terms]
            (+ (* acc tau)
               (reduce (fn [s [a b cc]] (+ s (* a (math/cos (+ b (* cc tau)))))) 0.0 terms)))
          0.0
          (rseq series)))

(defn heliocentric
  "Heliocentric `[lon lat r]` of `planet` at `mjd-tt`: ecliptic longitude
  and latitude in radians, referred to the mean ecliptic and equinox of
  date, and the distance from the Sun in AU."
  [planet mjd-tt]
  (let [tau (millennia mjd-tt)
        {:keys [L B R]} (data/series planet)]
    [(am/wrap-2pi (evaluate L tau)) (evaluate B tau) (evaluate R tau)]))

(defn ->fk5
  "`[lon lat]` of date in VSOP87's dynamical frame, corrected to the FK5
  frame (Meeus 32.3): a shift of 0.09 arcsecond along the ecliptic and a
  few hundredths across it."
  [[lon lat] mjd-tt]
  (let [T  (* 10.0 (millennia mjd-tt))
        l' (- lon (* c/degrees (+ (* 1.397 T) (* 0.00031 T T))))
        cl (math/cos l') sl (math/sin l')]
    [(+ lon (* c/arcsec (+ -0.09033 (* 0.03916 (+ cl sl) (math/tan lat)))))
     (+ lat (* c/arcsec 0.03916 (- cl sl)))]))

(defn rectangular
  "Heliocentric ecliptic rectangular coordinates of `planet`, AU, referred
  to the mean ecliptic and equinox of date."
  [planet mjd-tt]
  (let [[l b r] (heliocentric planet mjd-tt)]
    [(* r (math/cos b) (math/cos l))
     (* r (math/cos b) (math/sin l))
     (* r (math/sin b))]))

(defn ecliptic-of-date->J2000
  "A matrix taking rectangular coordinates on the mean ecliptic and
  equinox of `mjd-tt` to the mean equator and equinox of J2000 -- EME2000,
  the frame the rest of `allgo.astro` works in. Tilt by the obliquity of
  date, then undo the precession since J2000."
  [mjd-tt]
  (let [eps (frames/mean-obliquity mjd-tt)
        ce (math/cos eps) se (math/sin eps)]
    (lin/mat-mul (lin/transpose (frames/precession mjd-tt))
                 [[1.0 0.0 0.0] [0.0 ce (- se)] [0.0 se ce]])))

(defn equatorial-J2000
  "Heliocentric position of `planet`, AU, in EME2000 -- geometric, with no
  light time: where the planet is, which is what a model of the solar
  system draws."
  ([planet mjd-tt] (equatorial-J2000 planet mjd-tt (ecliptic-of-date->J2000 mjd-tt)))
  ([planet mjd-tt to-J2000]
   (lin/mat-vec to-J2000 (rectangular planet mjd-tt))))

(defn heliocentric-state
  "Position and velocity `[r v]` of `planet` relative to the Sun, km and
  km/s, in EME2000 -- the same shape as `allgo.astro.planets`'s, so it
  can stand in for that where arcseconds matter more than speed: in
  `allgo.astro.interplanetary/transfer` as its `:ephemeris`, say. The
  velocity is the position's central difference over two hours, whose
  error, of order the acceleration's rate times the step squared, is
  below a millimeter a second."
  [planet mjd-tt]
  (let [h (/ 1.0 24.0)
        m (ecliptic-of-date->J2000 mjd-tt)
        at #(lin/mat-vec m (rectangular planet %))
        au c/AU]
    [(mapv #(* au %) (at mjd-tt))
     (mapv #(/ (* au (- %1 %2)) (* 2.0 h 86400.0)) (at (+ mjd-tt h)) (at (- mjd-tt h)))]))
