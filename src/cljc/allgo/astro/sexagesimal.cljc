(ns allgo.astro.sexagesimal
  "Angles and times in base sixty: degrees, minutes and seconds of arc,
  hours, minutes and seconds of right ascension or of the day (Vallado,
  *Fundamentals of Astrodynamics and Applications*, algorithms 17-21).

  An hour of right ascension is 15 degrees, a minute of it 15 minutes of
  arc. The one trap is the sign: -0 degrees 30 minutes is half a degree
  south, and a zero of degrees cannot carry a minus sign. So going in, a
  minus on any of the parts makes the whole negative; coming out, the sign
  is kept apart, `:sign` -1 or 1, and the parts are all non-negative."
  (:require [clojure.math :as math]))

(defn- ->value
  "d + m/60 + s/3600, negative if any part is."
  [a b c]
  (let [sign (if (some neg? [a b c]) -1.0 1.0)]
    (* sign (+ (abs a) (/ (abs b) 60.0) (/ (abs c) 3600.0)))))

(defn- ->parts
  "`x` split into `{:sign :whole :minutes :seconds}`, the seconds rounded
  to `places` decimals and carried up when they round to 60."
  [x places]
  (let [sign (if (neg? x) -1 1)
        scale (math/pow 10.0 places)
        total (/ (math/round (* (abs x) 3600.0 scale)) scale)
        whole (long (math/floor (/ total 3600.0)))
        rest (- total (* 3600.0 whole))
        minutes (long (math/floor (/ rest 60.0)))
        seconds (- rest (* 60.0 minutes))]
    {:sign sign :whole whole :minutes minutes :seconds seconds}))

(defn dms->rad
  "Degrees, minutes and seconds of arc to radians (algorithm 17)."
  [d m s]
  (math/to-radians (->value d m s)))

(defn rad->dms
  "Radians to `{:sign :d :m :s}` (algorithm 18), the seconds rounded to
  `places` decimals (default 6)."
  ([a] (rad->dms a 6))
  ([a places]
   (let [{:keys [sign whole minutes seconds]} (->parts (math/to-degrees a) places)]
     {:sign sign :d whole :m minutes :s seconds})))

(defn hms->rad
  "Hours, minutes and seconds of right ascension to radians, 15 degrees
  an hour (algorithm 19)."
  [h m s]
  (math/to-radians (* 15.0 (->value h m s))))

(defn rad->hms
  "Radians to `{:sign :h :m :s}` of right ascension (algorithm 20), the
  seconds rounded to `places` decimals (default 6)."
  ([a] (rad->hms a 6))
  ([a places]
   (let [{:keys [sign whole minutes seconds]} (->parts (/ (math/to-degrees a) 15.0) places)]
     {:sign sign :h whole :m minutes :s seconds})))

(defn hms->seconds
  "A time of day, hours minutes and seconds, as seconds since midnight."
  [h m s]
  (* 3600.0 (->value h m s)))

(defn seconds->hms
  "Seconds since midnight as `{:h :m :s}` (algorithm 21), the seconds
  rounded to `places` decimals (default 6)."
  ([t] (seconds->hms t 6))
  ([t places]
   (let [{:keys [whole minutes seconds]} (->parts (/ t 3600.0) places)]
     {:h whole :m minutes :s seconds})))
