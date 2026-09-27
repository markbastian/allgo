(ns allgo.astro.sundial
  "Hour lines for a planar sundial of any orientation (Meeus,
  *Astronomical Algorithms*, chapter 58).

  The shadow of a straight stylus -- a gnomon of length `a` perpendicular
  to the dial's plane -- is cast by the Sun at hour angle H and
  declination d; where the tip's shadow falls is a central projection of
  the Sun's direction onto the plane. Each hour line is traced here
  through the Sun's declination on the first day of each zodiac sign, from
  the winter to the summer solstice, which are the curves an old sundial
  marks as its date lines.

  The plane is described by its gnomonic declination `D` -- the azimuth of
  the perpendicular to it, measured from the south toward the west, as
  dial makers do -- and its zenith distance `z`, 0 for a horizontal dial
  and 90 degrees for a vertical one. Coordinates are in the dial's plane,
  in units of the stylus: x horizontal, positive to the right as the dial
  is faced, y along its line of greatest slope, positive upward.

  Hours are local apparent solar time, 12 at noon; a dial reads the Sun,
  not a clock, and `allgo.astro.solar/equation-of-time` is the difference."
  (:require [clojure.math :as math]))

(def declinations
  "The Sun's declination when entering each zodiac sign, degrees -- the
  date lines each hour line is drawn through."
  [-23.44 -20.15 -11.47 0.0 11.47 20.15 23.44])

(defn- hour-lines [lat point]
  (let [tl (math/tan lat)]
    (for [hour (range 24)
          :let [H (* (- hour 12) 15.0 (/ math/PI 180.0))
                pts (for [d declinations
                          :let [td (math/tan (math/to-radians d))]
                          :when (<= (abs H) (math/acos (- (* tl td))))
                          :let [p (point H td)]
                          :when p]
                      p)]
          :when (seq pts)]
      {:hour hour :points (vec pts)})))

(defn general
  "Hour lines of a dial at latitude `lat` whose plane has gnomonic
  declination `D` and zenith distance `z`, for a stylus of length `a`:

    :lines   `[{:hour :points [[x y] ...]} ...]`, only where the Sun is up
             and shines on the dial's face
    :center  `[x y]` where the polar style, parallel to the Earth's axis,
             meets the plane
    :style   length of that polar style
    :psi     the angle it makes with the plane."
  [lat D a z]
  (let [sp (math/sin lat) cp (math/cos lat)
        sD (math/sin D) cD (math/cos D)
        sz (math/sin z) cz (math/cos z)
        P (- (* sp cz) (* cp sz cD))]
    {:lines (hour-lines lat
                        (fn [H td]
                          (let [sH (math/sin H) cH (math/cos H)
                                Q (+ (* sD sz sH) (* (+ (* cp cz) (* sp sz cD)) cH) (* P td))]
                            (when-not (neg? Q)
                              [(/ (* a (- (* cD sH) (* sD (- (* sp cH) (* cp td))))) Q)
                               (/ (* a (- (* cz sD sH) (* (- (* cp sz) (* sp cz cD)) cH)
                                          (* (+ (* sp sz) (* cp cz cD)) td)))
                                  Q)]))))
     :center [(* (/ a P) cp sD) (* (/ (- a) P) (+ (* sp sz) (* cp cz cD)))]
     :style (/ a (abs P))
     :psi (math/asin (abs P))}))

(defn equatorial
  "Hour lines of an equatorial dial, whose plane is parallel to the
  equator: `{:north lines :south lines}`, the face lit from spring to
  autumn and the one lit the rest of the year."
  [lat a]
  (let [lines (fn [pick]
                (hour-lines lat (fn [H td]
                                  (when (pick td)
                                    (let [x (/ (* (- a) (math/sin H)) td)
                                          y (/ (* a (math/cos H)) td)]
                                      [x (if (neg? td) y (- y))])))))]
    {:north (lines pos?) :south (lines neg?)}))

(defn horizontal
  "Hour lines of a horizontal dial at latitude `lat`, with `:center` and
  `:style` as for `general`."
  [lat a]
  (let [sp (math/sin lat) cp (math/cos lat)]
    {:lines (hour-lines lat (fn [H td]
                              (let [Q (+ (* cp (math/cos H)) (* sp td))]
                                [(/ (* a (math/sin H)) Q)
                                 (/ (* a (- (* sp (math/cos H)) (* cp td))) Q)])))
     :center [0.0 (/ (- a) (math/tan lat))]
     :style (/ a (abs sp))}))

(defn vertical
  "Hour lines of a vertical dial at latitude `lat` facing gnomonic
  declination `D`, with `:center` and `:style` as for `general`."
  [lat D a]
  (let [sp (math/sin lat) cp (math/cos lat)
        sD (math/sin D) cD (math/cos D)]
    {:lines (hour-lines lat (fn [H td]
                              (let [sH (math/sin H) cH (math/cos H)
                                    Q (- (+ (* sD sH) (* sp cD cH)) (* cp cD td))]
                                (when-not (neg? Q)
                                  [(/ (* a (+ (- (* cD sH) (* sp sD cH)) (* cp sD td))) Q)
                                   (/ (* (- a) (+ (* cp cH) (* sp td))) Q)]))))
     :center [(/ (* (- a) sD) cD) (/ (* a (math/tan lat)) cD)]
     :style (/ a (abs (* cp cD)))}))
