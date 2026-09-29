(ns allgo.astro.ephemeris
  "Low-precision Sun and Moon positions, and the third-body perturbation they
  exert (Montenbruck & Gill 3.3).

  These are analytical series, not a numerical ephemeris: a handful of terms
  good to about 1' for the Sun and a few arcminutes for the Moon. That is
  far short of what a real ephemeris gives, and entirely adequate here --
  the perturbation itself is a millionth of Earth's pull, so an arcminute of
  error in the Sun's direction is far below anything else in the budget.

  The accuracy is not uniform in time, and the way it degrades is worth
  knowing. Checked against the independent planetary elements of
  `allgo.astro.planets`, this Sun agrees to 3 arcseconds at J2000 and
  then drifts by about 12.7 a year -- 8 arcminutes by 2040, past its own
  nominal accuracy. The discrepancy is almost purely in longitude, the
  latitude staying within 20 arcseconds, which identifies it as a slightly
  short mean-longitude rate rather than a frame error: 35999.02 degrees a
  century where the fitted value is 35999.37. For a perturbation that is
  still irrelevant. For pointing, or for drawing the sky far from J2000,
  prefer `planets/sun-from-earth`.

  Positions are geocentric, equatorial, referred to the mean equator and
  equinox of J2000 -- the same frame `allgo.astro.forces` integrates
  in and `allgo.astro.frames` takes as its input. Everything in this
  package shares that frame, so positions from here can be drawn alongside a
  propagated orbit without further rotation."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.time :as t]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

(def obliquity-J2000
  "Obliquity of the ecliptic at J2000, radians."
  (* 23.43929111 c/degrees))

(defn obliquity
  "Obliquity of the ecliptic at `mjd-tt`, radians -- the tilt that turns the
  Sun's annual path along the ecliptic into a seasonal swing in declination.

  Provided for work referred to the equator of date. The series below are
  not: they use `obliquity-J2000`, because their longitudes are measured
  from the J2000 equinox and rotating those by an obliquity of date would
  produce a frame that is neither one thing nor the other."
  [mjd-tt]
  (* (- 23.43929111 (* 0.0130042 (t/centuries-J2000 mjd-tt))) c/degrees))

(defn sun
  "Geocentric position of the Sun, km."
  [mjd-tt]
  (let [T (t/centuries-J2000 mjd-tt)
        M (* c/two-pi (am/frac (+ 0.9931267 (* 99.9973583 T))))
        L (* c/two-pi (am/frac (+ 0.7859444 (/ M c/two-pi)
                                  (/ (+ (* 6892.0 (math/sin M))
                                        (* 72.0 (math/sin (* 2.0 M))))
                                     1296.0e3))))
        r (- 149.619e6 (* 2.499e6 (math/cos M)) (* 0.021e6 (math/cos (* 2.0 M))))]
    ;; The mean longitude advances 35999 degrees a century, not 36001, which
    ;; identifies the series as referred to the fixed J2000 equinox rather
    ;; than the moving equinox of date. The obliquity must match: using the
    ;; one of date puts the Sun 8000 km out by 2025 and 33,000 by 2100.
    (frames/ecliptic->equatorial [(* r (math/cos L)) (* r (math/sin L)) 0.0]
                                 obliquity-J2000)))

(defn moon
  "Geocentric position of the Moon, km.

  The Moon needs many more terms than the Sun because it is genuinely
  perturbed -- by the Sun, and strongly. The arguments below are its mean
  anomaly, the Sun's, their elongation and the argument of latitude, and
  most of the series is the Sun pulling the Moon about."
  [mjd-tt]
  (let [T   (t/centuries-J2000 mjd-tt)
        L0  (am/frac (+ 0.606433 (* 1336.851344 T)))
        l   (* c/two-pi (am/frac (+ 0.374897 (* 1325.552410 T))))
        lp  (* c/two-pi (am/frac (+ 0.993133 (* 99.997361 T))))
        D   (* c/two-pi (am/frac (+ 0.827361 (* 1236.853086 T))))
        F   (* c/two-pi (am/frac (+ 0.259086 (* 1342.227825 T))))
        sin math/sin
        cos math/cos
        dL  (+ (* 22640 (sin l))       (* -4586 (sin (- l (* 2 D))))
               (* 2370 (sin (* 2 D)))  (* 769 (sin (* 2 l)))
               (* -668 (sin lp))       (* -412 (sin (* 2 F)))
               (* -212 (sin (- (* 2 l) (* 2 D))))
               (* -206 (sin (- (+ l lp) (* 2 D))))
               (* 192 (sin (+ l (* 2 D))))
               (* -165 (sin (- lp (* 2 D))))
               (* -125 (sin D))        (* -110 (sin (+ l lp)))
               (* 148 (sin (- l lp)))  (* -55 (sin (- (* 2 F) (* 2 D)))))
        L   (* c/two-pi (am/frac (+ L0 (/ dL 1296.0e3))))
        S   (+ F (/ (+ dL (* 412 (sin (* 2 F))) (* 541 (sin lp))) (/ 1.0 c/arcsec)))
        h   (- F (* 2 D))
        N   (+ (* -526 (sin h))          (* 44 (sin (+ l h)))
               (* -31 (sin (+ (- l) h))) (* -23 (sin (+ lp h)))
               (* 11 (sin (+ (- lp) h))) (* -25 (sin (+ (* -2 l) F)))
               (* 21 (sin (+ (- l) F))))
        B   (* (+ (* 18520.0 (sin S)) N) c/arcsec)
        R   (- 385000.0 (* 20905.0 (cos l))
               (* 3699.0 (cos (- (* 2 D) l)))
               (* 2956.0 (cos (* 2 D)))
               (* 570.0 (cos (* 2 l)))
               (* -246.0 (cos (- (* 2 l) (* 2 D))))
               (* 205.0 (cos (- lp (* 2 D))))
               (* 171.0 (cos (+ l (* 2 D))))
               (* 152.0 (cos (- (+ l lp) (* 2 D)))))]
    ;; L0 is likewise measured from the J2000 equinox.
    (frames/ecliptic->equatorial [(* R (cos L) (cos B)) (* R (sin L) (cos B)) (* R (sin B))]
                                 obliquity-J2000)))

;; ------------------------------------------------------------- the force

(defn third-body
  "Acceleration on a satellite at `r` from a point mass `GM` at `s`, both
  geocentric, km/s^2 (M&G eq. 3.37).

  Two terms, and the second matters more than the first: the body pulls on
  the satellite, but it also pulls on the Earth, and only the difference
  perturbs the orbit. For the Sun on a 7000 km satellite the raw pull is
  6.1e-6 km/s^2 while the perturbation is 3.1e-10 -- four orders smaller.
  Omitting the indirect term does not introduce a small error, it reports
  something else entirely."
  [GM r s]
  (let [d   (v3/sub r s)
        dm  (v3/length d)
        sm  (v3/length s)
        kd  (/ GM (* dm dm dm))
        ks  (/ GM (* sm sm sm))]
    (mapv (fn [di si] (- (- (* kd di)) (* ks si))) d s)))
