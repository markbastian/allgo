(ns allgo.astro.rise
  "Rising, transit and setting (Meeus, *Astronomical Algorithms*, chapter
  15).

  A body rises when its center reaches the altitude h0 -- below the
  horizon, not on it: refraction lifts everything by some 34 arcminutes
  there, and for the Sun the event is the upper limb touching, another 16
  arcminutes. The hour angle at which that altitude is reached follows
  from the body's declination and the observer's latitude, and turns into
  a time through the sidereal time. That first answer treats the body as
  fixed for the day; for the Moon, which moves thirteen degrees a day, and
  to a lesser degree for the Sun and planets, each time is then corrected
  by interpolating the position to it.

  Results are fractions of the UT day, 0 to 1. Longitude is positive east,
  as throughout `allgo.astro`."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.frames :as frames]
            [allgo.astro.moon :as moon]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [allgo.math :as am]
            [allgo.numerics.interpolation :as interp]
            [clojure.math :as math]))

(def standard-altitude
  "h0 for a star or planet, and for the Sun's upper limb: the geometric
  altitude of the center at the instant of rising or setting."
  {:star (* (/ -34.0 60.0) c/degrees)
   :sun  (* (/ -50.0 60.0) c/degrees)})

(defn moon-standard-altitude
  "h0 for the Moon's upper limb, which depends on its parallax -- the Moon
  is close enough for the observer's place on the Earth to lower it by up
  to a degree."
  [parallax]
  (- (* 0.7275 parallax) (* (/ 34.0 60.0) c/degrees)))

(defn- frac [x] (- x (math/floor x)))

(defn approximate
  "`{:rise :transit :set}` as fractions of the UT day, for a body fixed at
  `[ra dec]`, seen from `[lat lon]`, when the apparent sidereal time at
  Greenwich at 0h UT is `theta0`; nil for a body that never crosses h0 --
  circumpolar, or never up."
  [[lat lon] h0 theta0 [ra dec]]
  (let [cH0 (/ (- (math/sin h0) (* (math/sin lat) (math/sin dec)))
               (* (math/cos lat) (math/cos dec)))]
    (when (<= -1.0 cH0 1.0)
      (let [H0 (/ (math/acos cH0) c/two-pi)
            m0 (/ (- ra lon theta0) c/two-pi)]
        {:rise (frac (- m0 H0)) :transit (frac m0) :set (frac (+ m0 H0))}))))

(defn times
  "`{:rise :transit :set}` corrected for the body's motion: `positions`
  are its apparent `[ra dec]` at 0h TT on the day before, the day, and the
  day after, and `delta-t` is TT - UT in seconds. Each time moves by
  seconds for a planet and by up to an hour for the Moon."
  [[lat lon :as place] h0 theta0 delta-t positions]
  (when-let [approx (approximate place h0 theta0 (second positions))]
    (let [unwrap (fn [as] (reductions (fn [a b] (+ a (am/wrap-angle (- b a)))) as))
          ra-t  (interp/table-3 -1.0 1.0 (unwrap (map first positions)))
          dec-t (interp/table-3 -1.0 1.0 (map second positions))
          at (fn [m]
               (let [n (+ m (/ delta-t 86400.0))
                     theta (+ theta0 (* c/two-pi 1.00273790935 m))
                     ra (interp/value-n ra-t n)]
                 [(am/wrap-angle (- (+ theta lon) ra)) (interp/value-n dec-t n)]))
          transit (let [m (:transit approx)
                        [H] (at m)]
                    (- m (/ H c/two-pi)))
          rise-set (fn [m]
                     (let [[H dec] (at m)
                           h (math/asin (+ (* (math/sin lat) (math/sin dec))
                                           (* (math/cos lat) (math/cos dec) (math/cos H))))]
                       (+ m (/ (- h h0) (* c/two-pi (math/cos dec) (math/cos lat) (math/sin H))))))]
      {:rise (rise-set (:rise approx)) :transit transit :set (rise-set (:set approx))})))

(defn- day-events
  "`times` for a body whose apparent `[ra dec]` at an MJD (TT) is
  `(position mjd)`, on the UT day starting at `mjd` (0h)."
  [position place h0 mjd]
  (let [dt (time/delta-t (time/decimal-year mjd))
        theta0 (frames/gast mjd (+ mjd (/ dt 86400.0)))]
    (times place h0 theta0 dt (map #(position (+ mjd %)) [-1.0 0.0 1.0]))))

(defn planet
  "Rising, transit and setting of `planet` on the UT day beginning at
  `mjd`, from `[lat lon]`."
  [planet place mjd]
  (day-events #(vec (take 2 (elliptic/apparent planet %))) place (:star standard-altitude) mjd))

(defn sun
  "Sunrise, the Sun's transit and sunset on the UT day beginning at `mjd`."
  [place mjd]
  (day-events #(vec (take 2 (solar/apparent-equatorial %))) place (:sun standard-altitude) mjd))

(defn moon
  "Moonrise, transit and moonset on the UT day beginning at `mjd`. The
  Moon rises some 50 minutes later each day, so on one day a month it
  does not rise at all and the time found falls outside the day."
  [place mjd]
  (let [[_ _ r] (moon/position mjd)]
    (day-events #(vec (take 2 (moon/apparent-equatorial %))) place
                (moon-standard-altitude (moon/horizontal-parallax r)) mjd)))
