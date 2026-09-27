(ns allgo.astro.phenomena
  "Conjunctions, oppositions, elongations and stations of the planets
  from mean periods and periodic corrections (Meeus, *Astronomical
  Algorithms*, chapter 36).

  Each phenomenon recurs at the planet's synodic period on average; Meeus
  gives the mean instant as a linear function of a count k, and the
  departures from it -- up to a week or two, from the eccentricities of
  both orbits -- as a Fourier series in the mean anomaly M of the mean
  instant, with coefficients that are themselves polynomials in time. The
  result is good to a few minutes for Mercury and Venus and better than an
  hour for the outer planets, between about 1000 and 3000.

  Only the rows of Meeus's table 36.B that could be checked against an
  independent transcription are included: see `events`. The others --
  Venus's superior conjunction, elongations and stations, the remaining
  stations of Mercury and Mars, and the conjunctions of Mars, Jupiter,
  Uranus and Neptune -- are left out rather than entered unverified.

  Times are MJD (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.interpolation :refer [horner]]
            [clojure.math :as math]))

(def ^:private tables
  "Per phenomenon: `:mean` is [A B M0 M1], the mean instant A + k B and
  mean anomaly M0 + k M1 in degrees; `:terms` the constant and then
  alternating sine and cosine coefficients of M, 2M, 3M, ..., each a
  polynomial in T; `:extra` further arguments [c f] (c + f T, degrees)
  whose sine and cosine coefficients are the last rows of `:terms`. For
  Mercury's elongations, `:elongation` is the angle itself, in degrees."
  {:mercury-inferior-conjunction
   {:mean [2451612.023 115.8774771 63.5867 114.2088742]
    :terms [[0.0545 0.0002]
            [-6.2008 0.0074 0.00003]
            [-3.275 -0.0197 0.00001]
            [0.4737 -0.0052 -0.00001]
            [0.8111 0.0033 -0.00002]
            [0.0037 0.0018]
            [-0.1768 0.0 0.00001]
            [-0.0211 -0.0004]
            [0.0326 -0.0003]
            [0.0083 0.0001]
            [-0.004 0.0001]]}
   :mercury-superior-conjunction
   {:mean [2451554.084 115.8774771 6.4822 114.2088742]
    :terms [[-0.0548 -0.0002]
            [7.3894 -0.01 -0.00003]
            [3.22 0.0197 -0.00001]
            [0.8383 -0.0064 -0.00001]
            [0.9666 0.0039 -0.00003]
            [0.077 -0.0026]
            [0.2758 0.0002 -0.00002]
            [-0.0128 -0.0008]
            [0.0734 -0.0004 -0.00001]
            [-0.0122 -0.0002]
            [0.0173 -0.0002]]}
   :venus-inferior-conjunction
   {:mean [2451996.706 583.921361 82.7311 215.513058]
    :terms [[-0.0096 0.0002 -0.00001]
            [2.0009 -0.0033 -0.00001]
            [0.598 -0.0104 0.00001]
            [0.0967 -0.0018 -0.00003]
            [0.0913 0.0009 -0.00002]
            [0.0046 -0.0002]
            [0.0079 0.0001]]}
   :mars-opposition
   {:mean [2452097.382 779.936104 181.9573 48.705244]
    :terms [[-0.3088 0.0 0.00002]
            [-17.6965 0.0363 0.00005]
            [18.3131 0.0467 -0.00006]
            [-0.2162 -0.0198 -0.00001]
            [-4.5028 -0.0019 0.00007]
            [0.8987 0.0058 -0.00002]
            [0.7666 -0.005 -0.00003]
            [-0.3636 -0.0001 0.00002]
            [0.0402 0.0032]
            [0.0737 -0.0008]
            [-0.098 -0.0011]]}
   :mars-second-station
   {:mean [2452097.382 779.936104 181.9573 48.705244]
    :terms [[36.7191 0.0016 0.00003]
            [-12.6163 0.0417 -0.00001]
            [20.1218 0.0379 -0.00006]
            [-1.636 -0.019]
            [-3.9657 0.0045 0.00007]
            [1.1546 0.0029 -0.00003]
            [0.2888 -0.0073 -0.00002]
            [-0.3128 0.0017 0.00002]
            [0.2513 0.0026 -0.00002]
            [-0.0021 -0.0016]
            [-0.1497 -0.0006]]}
   :jupiter-opposition
   {:mean [2451870.628 398.884046 318.4681 33.140229]
    :terms [[-0.1029 0.0 -0.00009]
            [-1.9658 -0.0056 0.00007]
            [6.1537 0.021 -0.00006]
            [-0.2081 -0.0013]
            [-0.1116 -0.001]
            [0.0074 0.0001]
            [-0.0097 -0.0001]
            [0.0 0.0144 -0.00008]
            [0.3642 -0.0019 -0.00029]]
    :extra [[82.74 40.76]]}
   :saturn-opposition
   {:mean [2451870.17 378.091904 318.0172 12.647487]
    :terms [[-0.0209 0.0006 0.00023]
            [4.5795 -0.0312 -0.00017]
            [1.1462 -0.0351 0.00011]
            [0.0985 -0.0015]
            [0.0733 -0.0031 0.00001]
            [0.0025 -0.0001]
            [0.005 -0.0002]
            [0.0 -0.0337 0.00018]
            [-0.851 0.0044 0.00068]
            [0.0 -0.0064 0.00004]
            [0.2397 -0.0012 -0.00008]
            [0.0 -0.001]
            [0.1245 0.0006]
            [0.0 0.0024 -0.00003]
            [0.0477 -0.0005 -0.00006]]
    :extra [[82.74 40.76] [29.86 1181.36] [14.13 590.68] [220.02 1262.87]]}
   :saturn-conjunction
   {:mean [2451681.124 378.091904 131.6934 12.647487]
    :terms [[0.0172 -0.0006 0.00023]
            [-8.5885 0.0411 0.00020]
            [-1.147 0.0352 -0.00011]
            [0.3331 -0.0034 -0.00001]
            [0.1145 -0.0045 0.00002]
            [-0.0169 0.0002]
            [-0.0109 0.0004]
            [0.0 -0.0337 0.00018]
            [-0.851 0.0044 0.00068]
            [0.0 -0.0064 0.00004]
            [0.2397 -0.0012 -0.00008]
            [0.0 -0.001]
            [0.1245 0.0006]
            [0.0 0.0024 -0.00003]
            [0.0477 -0.0005 -0.00006]]
    :extra [[82.74 40.76] [29.86 1181.36] [14.13 590.68] [220.02 1262.87]]}
   :uranus-opposition
   {:mean [2451764.317 369.656035 213.6884 4.333093]
    :terms [[0.0844 -0.0006]
            [-0.1048 0.0246]
            [-5.1221 0.0104 0.00003]
            [-0.1428 0.0005]
            [-0.0148 -0.0013]
            [0.0]
            [0.0055]
            [0.0]
            [0.885]
            [0.0]
            [0.2153]]
    :extra [[207.83 8.51] [108.84 419.96]]}
   :neptune-opposition
   {:mean [2451753.122 367.486703 202.6544 2.194998]
    :terms [[-0.014 0.0 0.00001]
            [-1.3486 0.001 0.00001]
            [0.8597 0.0037]
            [-0.0082 -0.0002 0.00001]
            [0.0037 -0.0003]
            [0.0]
            [-0.5964]
            [0.0]
            [0.0728]]
    :extra [[207.83 8.51] [276.74 209.98]]}
   :mercury-greatest-east-elongation
   {:mean [2451612.023 115.8774771 63.5867 114.2088742]
    :terms [[-21.6106 0.0002]
            [-1.9803 -0.006 0.00001]
            [1.4151 -0.0072 -0.00001]
            [0.5528 -0.0005 -0.00001]
            [0.2905 0.0034 0.00001]
            [-0.1121 -0.0001 0.00001]
            [-0.0098 -0.0015]
            [0.0192]
            [0.0111 0.0004]
            [-0.0061]
            [-0.0032 -0.0001]]
    :elongation [[22.4697]
                 [-4.2666 0.0054 0.00002]
                 [-1.8537 -0.0137]
                 [0.3598 0.0008 -0.00001]
                 [-0.068 0.0026]
                 [-0.0524 -0.0003]
                 [0.0052 -0.0006]
                 [0.0107 0.0001]
                 [-0.0013 0.0001]
                 [-0.0021]
                 [0.0003]]}
   :mercury-greatest-west-elongation
   {:mean [2451612.023 115.8774771 63.5867 114.2088742]
    :terms [[21.6249 -0.0002]
            [0.1306 0.0065]
            [-2.7661 -0.0011 0.00001]
            [0.2438 -0.0024 -0.00001]
            [0.5767 0.0023]
            [0.1041]
            [-0.0184 0.0007]
            [-0.0051 -0.0001]
            [0.0048 0.0001]
            [0.0026]
            [0.0037]]
    :elongation [[22.4143 -0.0001]
                 [4.3651 -0.0048 -0.00002]
                 [2.3787 0.0121 -0.00001]
                 [0.2674 0.0022]
                 [-0.3873 0.0008 0.00001]
                 [-0.0369 -0.0001]
                 [0.0017 -0.0001]
                 [0.0059]
                 [0.0061 0.0001]
                 [0.0007]
                 [-0.0011]]}})

(def events
  "The phenomena available, the keys of `phenomenon`."
  (set (keys tables)))

(defn- series [T M rows]
  (reduce + (horner T (first rows))
          (map-indexed (fn [i [s cc]]
                         (let [a (* (inc i) M)]
                           (+ (* (math/sin a) (horner T s)) (* (math/cos a) (horner T cc)))))
                       (partition 2 (rest rows)))))

(defn phenomenon
  "MJD (TT) of `event` -- one of `events` -- nearest the decimal `year`.
  For Mercury's greatest elongations, `[mjd elongation]`, the elongation
  in radians."
  [event year]
  (let [{[A B M0 M1] :mean :keys [terms extra elongation]} (tables event)
        k  (math/floor (+ 0.5 (/ (- (+ (* 365.2425 year) 1721060.0) A) B)))
        J  (+ A (* k B))
        M  (* c/degrees (mod (+ M0 (* k M1)) 360.0))
        T  (/ (- J 2451545.0) 36525.0)
        n  (- (count terms) (* 2 (count extra)))
        dJ (+ (series T M (take n terms))
              (reduce + (map (fn [[cc f] [s co]]
                               (let [a (* c/degrees (+ cc (* f T)))]
                                 (+ (* (math/sin a) (horner T s)) (* (math/cos a) (horner T co)))))
                             extra (partition 2 (drop n terms)))))
        mjd (- (+ J dJ) c/jd-mjd-offset)]
    (if elongation
      [mjd (* c/degrees (series T M elongation))]
      mjd)))
