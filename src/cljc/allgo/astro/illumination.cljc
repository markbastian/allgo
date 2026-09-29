(ns allgo.astro.illumination
  "Phases, magnitudes and apparent sizes of the planets (Meeus,
  *Astronomical Algorithms*, chapters 41 and 55).

  A planet's phase is fixed by one triangle -- the Sun, the planet and the
  Earth -- and the phase angle at the planet is its only free angle. Its
  brightness follows from the same triangle: the inverse squares of both
  sides that meet at the planet, times a phase law measured for each
  planet, since how a surface scatters light at a slant depends on what
  it is made of.

  Distances are AU, angles radians."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

(defn- log10 [x] (/ (math/log x) (math/log 10.0)))

;; ------------------------------------------------------------- phase (41)

(defn phase-angle
  "The Sun-planet-Earth angle, from the planet's distance from the Sun
  `r`, from the Earth `delta`, and the Earth's from the Sun `R`."
  [r delta R]
  (math/acos (/ (- (+ (* r r) (* delta delta)) (* R R)) (* 2.0 r delta))))

(defn phase-angle-from-positions
  "The phase angle from the planet's heliocentric `[lon lat r]`, the
  Earth's `[lon0 _ R0]` and the Earth-planet distance -- the form that
  avoids a separate distance to the Sun."
  [[L B r] [L0 _ R0] delta]
  (math/acos (/ (- r (* R0 (math/cos B) (math/cos (- L L0)))) delta)))

(defn illuminated-fraction
  "Fraction of the disk lit, from the same three distances:
  (1 + cos i) / 2, written without the angle."
  [r delta R]
  (let [s (+ r delta)]
    (/ (- (* s s) (* R R)) (* 4.0 r delta))))

(defn venus-illuminated-fraction
  "Venus's illuminated fraction at `mjd-tt` from mean orbits alone -- the
  Meeus shortcut good to a few thousandths."
  [mjd-tt]
  (let [T (/ (- mjd-tt c/mjd-J2000) 36525.0)
        d c/degrees
        V (* d (+ 261.51 (* 22518.443 T)))
        M (* d (+ 177.53 (* 35999.05 T)))
        N (* d (+ 50.42 (* 58517.811 T)))
        W (+ V (* d 1.91 (math/sin M)) (* d 0.78 (math/sin N)))
        delta (math/sqrt (+ 1.52321 (* 1.44666 (math/cos W))))
        s (+ 0.72333 delta)]
    (/ (- (* s s) 1.0) (* 2.89332 delta))))

;; -------------------------------------------------------- magnitudes (41)

(defn magnitude
  "Visual magnitude of `planet` at distance `r` from the Sun and `delta`
  from the Earth, with phase angle `i` -- Müller's formulas, as Meeus
  gives them. Saturn's needs the ring as well: `opts` `{:b :du}` from
  `allgo.astro.physical/saturn-ring`.

  With `{:almanac true}`, the formulas the Astronomical Almanac used from
  1984, which also cover Pluto; Müller's date from 1893 and fit old
  observations better."
  ([planet r delta i] (magnitude planet r delta i {}))
  ([planet r delta i {:keys [b du almanac]}]
   (let [base (* 5.0 (log10 (* r delta)))
         id (/ i c/degrees)
         ring (fn [] (let [s (abs (math/sin (or b 0.0)))]
                       (+ (* 0.044 (abs (/ (or du 0.0) c/degrees))) (* -2.60 s) (* 1.25 s s))))]
     (if almanac
       (case planet
         :mercury (+ base (poly/horner-ascending [-0.42 0.038 -0.000273 0.000002] id))
         ;; Harris's phase law; some transcriptions of Meeus flip the
         ;; signs of the last two terms, which would have Venus brighten
         ;; as its phase thins
         :venus   (+ base (poly/horner-ascending [-4.40 0.0009 0.000239 -0.00000065] id))
         :mars    (+ base -1.52 (* 0.016 id))
         :jupiter (+ base -9.40 (* 0.005 id))
         :saturn  (+ base -8.88 (ring))
         :uranus  (+ base -7.19)
         :neptune (+ base -6.87)
         :pluto   (+ base -1.00))
       (case planet
         :mercury (let [s (- id 50.0)] (+ base 1.16 (* s (+ 0.02838 (* 0.0001023 s)))))
         :venus   (+ base -4.00 (* id (+ 0.01322 (* 0.0000004247 id id))))
         :mars    (+ base -1.30 (* 0.01486 id))
         :jupiter (+ base -8.93)
         :saturn  (+ base -8.68 (ring))
         :uranus  (+ base -6.85)
         :neptune (+ base -7.05))))))

;; ------------------------------------------------------ semidiameters (55)

(def semidiameters
  "Semidiameters at one AU, radians (Meeus table 55.A). Venus's is given
  to the cloud tops, which is what is seen; `:venus-surface` for the
  solid body. The giants are flattened enough for their polar values to
  differ visibly."
  (into {} (map (fn [[k v]] [k (* v c/arcsec)]))
        {:sun 959.63 :mercury 3.36 :venus 8.41 :venus-surface 8.34 :mars 4.68
         :jupiter 98.44 :jupiter-polar 92.06 :saturn 82.73 :saturn-polar 73.82
         :uranus 35.02 :neptune 33.50 :pluto 2.07}))

(defn semidiameter
  "Apparent semidiameter of `body` (a key of `semidiameters`) at `delta` AU."
  [body delta]
  (/ (semidiameters body) delta))

(defn saturn-polar-semidiameter
  "Saturn's apparent polar semidiameter at `delta` AU when its pole is
  tipped by `b` (the ring's `:b`) toward the Earth: an ellipse seen
  obliquely, between the polar and equatorial values."
  [delta b]
  (let [k (/ (semidiameters :saturn-polar) (semidiameters :saturn))
        k (- 1.0 (* k k))
        cb (math/cos b)]
    (/ (* (semidiameters :saturn) (math/sqrt (- 1.0 (* k cb cb)))) delta)))

(defn moon-semidiameter
  "The Moon's geocentric semidiameter at `distance` km. The ratio of its
  radius to the Earth's, 0.272481, is what eclipse work uses."
  [distance]
  (math/asin (* 0.272481 (/ c/R-earth distance))))

(defn moon-topocentric-semidiameter
  "The Moon's semidiameter seen from a place with parallax constants
  `[rho-sin rho-cos]` when it is at declination `dec` and hour angle `H`
  and `distance` km away: larger than the geocentric value by up to 1.8
  percent at the zenith, where the observer is an Earth radius closer."
  [distance dec H [rho-sin rho-cos]]
  (let [sp (/ c/R-earth distance)
        A (* (math/cos dec) (math/sin H))
        B (- (* (math/cos dec) (math/cos H)) (* rho-cos sp))
        C (- (math/sin dec) (* rho-sin sp))
        q (math/sqrt (+ (* A A) (* B B) (* C C)))]
    (math/asin (/ (* 0.272481 sp) q))))

(defn asteroid-diameter
  "An asteroid's diameter in km from its absolute magnitude `H` and albedo
  `albedo` -- the only way to size one that has not been visited or
  resolved."
  [H albedo]
  (math/pow 10.0 (- 3.12 (* 0.2 H) (* 0.5 (log10 albedo)))))

(defn asteroid-semidiameter
  "Apparent semidiameter of an asteroid of diameter `d` km at `delta` AU."
  [d delta]
  (* 0.0013788 c/arcsec (/ d delta)))
