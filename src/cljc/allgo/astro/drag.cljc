(ns allgo.astro.drag
  "Atmospheric drag and the Harris-Priester density model
  (Montenbruck & Gill 3.5).

  Drag is the one force here that is not conservative: it takes energy out
  of the orbit every revolution and never gives it back, so it is what
  eventually brings a satellite down. It is also the least predictable,
  since the density it depends on varies by an order of magnitude with solar
  activity -- which is why Harris-Priester tabulates a minimum and a maximum
  and interpolates between them rather than pretending to a single value."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as geodesy]
            [clojure.math :as math]))

(def table
  "Harris-Priester density, in g/km^3 (that is, 1e-12 kg/m^3), against
  altitude in km. The two columns are the antapex and apex of the diurnal
  bulge: the atmosphere is hotter and so puffier on the daylit side, and the
  difference between the columns grows from nothing at 100 km to a factor of
  fifteen at 1000."
  [[100 4.974e+05 4.974e+05] [120 2.490e+04 2.490e+04] [130 8.377e+03 8.710e+03]
   [140 3.899e+03 4.059e+03] [150 2.122e+03 2.215e+03] [160 1.263e+03 1.344e+03]
   [170 8.008e+02 8.758e+02] [180 5.283e+02 6.010e+02] [190 3.617e+02 4.297e+02]
   [200 2.557e+02 3.162e+02] [210 1.839e+02 2.396e+02] [220 1.341e+02 1.853e+02]
   [230 9.949e+01 1.455e+02] [240 7.488e+01 1.157e+02] [250 5.709e+01 9.308e+01]
   [260 4.403e+01 7.555e+01] [270 3.430e+01 6.182e+01] [280 2.697e+01 5.095e+01]
   [290 2.139e+01 4.226e+01] [300 1.708e+01 3.526e+01] [320 1.099e+01 2.511e+01]
   [340 7.214e+00 1.819e+01] [360 4.824e+00 1.337e+01] [380 3.274e+00 9.955e+00]
   [400 2.249e+00 7.492e+00] [420 1.558e+00 5.684e+00] [440 1.091e+00 4.355e+00]
   [460 7.701e-01 3.362e+00] [480 5.474e-01 2.612e+00] [500 3.916e-01 2.042e+00]
   [520 2.819e-01 1.605e+00] [540 2.042e-01 1.267e+00] [560 1.488e-01 1.005e+00]
   [580 1.092e-01 7.997e-01] [600 8.070e-02 6.390e-01] [620 6.012e-02 5.123e-01]
   [640 4.519e-02 4.121e-01] [660 3.430e-02 3.325e-01] [680 2.632e-02 2.691e-01]
   [700 2.043e-02 2.185e-01] [720 1.607e-02 1.779e-01] [740 1.281e-02 1.452e-01]
   [760 1.036e-02 1.190e-01] [780 8.496e-03 9.776e-02] [800 7.069e-03 8.059e-02]
   [840 4.680e-03 5.741e-02] [880 3.200e-03 4.210e-02] [920 2.210e-03 3.130e-02]
   [960 1.560e-03 2.360e-02] [1000 1.150e-03 1.810e-02]])

(def bulge-lag
  "The diurnal bulge trails the Sun by about 30 degrees of longitude: the
  atmosphere takes a couple of hours to heat, so it is thickest in the
  mid-afternoon rather than at noon."
  (* 30.0 c/degrees))

(defn geodetic-height
  "Height above the reference ellipsoid, km.

  Not r - R: the Earth is 21 km narrower through the poles, and against an
  atmospheric scale height of some 50 km that is a factor of about 1.5 in
  density. Delegated to `allgo.astro.geodesy`, which solves the same
  ellipsoid problem for coordinates."
  [r]
  (nth (geodesy/cartesian->geodetic r) 2))

(defn- interpolate
  "Exponential interpolation between table rows, which is what an atmosphere
  in hydrostatic equilibrium actually does -- density falls by a constant
  factor per scale height, so a straight line between entries would be badly
  wrong across a 20 km gap."
  [h]
  (let [rows table]
    (cond
      (<= h (first (first rows)))   [(second (first rows)) (nth (first rows) 2)]
      (>= h (first (peek rows)))    [(second (peek rows)) (nth (peek rows) 2)]
      :else
      (let [[lo hi] (first (filter (fn [[a b]] (and (<= (first a) h) (< h (first b))))
                                   (partition 2 1 rows)))
            [h0 min0 max0] lo
            [h1 min1 max1] hi
            scale (fn [r0 r1] (let [H (/ (- h0 h1) (math/log (/ r1 r0)))]
                                (* r0 (math/exp (/ (- h0 h) H)))))]
        [(scale min0 min1) (scale max0 max1)]))))

(defn density
  "Atmospheric density at `r`, kg/m^3, with the Sun at `r-sun`.

  `n` is the bulge exponent: 2 near the equator, 6 for a polar orbit. It
  controls how sharply the daylit maximum falls off toward the night side."
  ([r r-sun] (density r r-sun 4.0))
  ([r r-sun n]
   (let [h (geodetic-height r)]
     (if (> h 1000.0)
       0.0
       (let [[rho-min rho-max] (interpolate (max 100.0 h))
             ;; apex of the bulge: the Sun's direction, turned east by the lag
             [sx sy sz] r-sun
             sm  (math/sqrt (+ (* sx sx) (* sy sy) (* sz sz)))
             ra  (+ (math/atan2 sy sx) bulge-lag)
             dec (math/asin (/ sz sm))
             u   [(* (math/cos dec) (math/cos ra))
                  (* (math/cos dec) (math/sin ra))
                  (math/sin dec)]
             rm  (math/sqrt (reduce + (map * r r)))
             cos-psi (/ (reduce + (map * r u)) rm)
             ;; cos^n(psi/2), written through the half-angle identity so it
             ;; stays defined when the satellite is opposite the bulge
             f   (math/pow (max 0.0 (* 0.5 (+ 1.0 cos-psi))) (* 0.5 n))]
         (* 1e-12 (+ rho-min (* (- rho-max rho-min) f))))))))

(defn relative-velocity
  "Velocity with respect to the air, km/s.

  The atmosphere turns with the Earth, so a satellite in a prograde orbit is
  not flying through still air -- at the equator the ground moves east at
  0.46 km/s, which is 6% of orbital speed and shifts the drag direction as
  well as its size."
  [r v]
  ;; omega x r, with omega along the pole, has no z component -- the
  ;; atmosphere carries a satellite east, never up or down.
  (let [[x y] r
        [vx vy vz] v
        w c/omega-earth]
    [(+ vx (* w y)) (- vy (* w x)) vz]))

(defn acceleration
  "Drag acceleration, km/s^2.

  `area-to-mass` in m^2/kg, `cd` the drag coefficient -- about 2.2 for a
  compact satellite in free-molecular flow, where molecules bounce off
  individually rather than forming a boundary layer, and which is why it
  exceeds the value any wind tunnel would give."
  [r v r-sun area-to-mass cd & [n]]
  (let [rho (density r r-sun (or n 4.0))]
    (if (zero? rho)
      [0.0 0.0 0.0]
      (let [vr (relative-velocity r v)
            s  (math/sqrt (reduce + (map * vr vr)))
            ;; 1e3 converts the m/s^2 that the SI quantities yield into km/s^2
            ;; while velocities stay in km/s
            k  (* -0.5 cd area-to-mass rho 1e3 s)]
        (mapv #(* k %) vr)))))
