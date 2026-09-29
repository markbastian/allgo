(ns allgo.astro.geomagnetic
  "The planetary geomagnetic indices Kp and ap (Bartels 1949, 1957; Matzka
  et al., Space Weather 19, 2021) and the conversion between them (Vallado,
  algorithm 62).

  Kp grades each three hours' disturbance on a quasi-logarithmic scale of
  28 steps, 0 to 9 in thirds (0o, 0+, 1-, 1o, ...); ap puts the same steps
  on a linear scale, in units of 2 nT at 50 degrees geomagnetic latitude,
  from 0 to 400. The atmosphere models want one or the other -- Jacchia's
  Kp, MSIS and JB2008 ap or its daily mean Ap -- so a value between the
  steps, a daily mean say, is converted by interpolating the table
  linearly between them."
  (:require [clojure.math :as math]))

(def ap-steps
  "ap for each Kp step, Kp = 0, 1/3, 2/3, ... 9."
  [0 2 3 4 5 6 7 9 12 15 18 22 27 32 39 48 56 67 80 94 111 132 154 179 207 236 300 400])

(defn kp->ap
  "ap for `kp`, 0 to 9: the table's value on a third, interpolated
  between."
  [kp]
  (let [x (* 3.0 (max 0.0 (min 9.0 kp)))
        i (min 26 (long (math/floor x)))
        f (- x i)]
    (+ (* (- 1.0 f) (ap-steps i)) (* f (ap-steps (inc i))))))

(defn ap->kp
  "Kp for `ap`, 0 to 400: the inverse of `kp->ap`."
  [ap]
  (let [ap (max 0.0 (min 400.0 ap))
        i (or (first (filter #(<= (ap-steps %) ap (ap-steps (inc %))) (range 27))) 26)
        lo (ap-steps i) hi (ap-steps (inc i))]
    (/ (+ i (/ (- ap lo) (- hi lo))) 3.0)))
