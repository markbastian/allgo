(ns allgo.astro.lunar-events
  "When the Moon reaches its phases, its perigee and apogee, its nodes and
  its greatest declinations (Meeus, *Astronomical Algorithms*, chapters 49
  to 52).

  All four are the same device. Each event recurs with a mean period --
  the synodic month, the anomalistic month, the draconic month, the
  tropical month -- so the mean instant of the k-th is linear in k, with a
  small secular term. The true instant departs from the mean by up to
  half a day, and those departures are Fourier series in the arguments of
  lunar theory at the mean instant. No position of the Moon is computed:
  the series are the Moon's motion, pre-solved for the moment of interest.

  Times are MJD (TT); `year` is a decimal year near the wanted event."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.interpolation :refer [horner]]
            [clojure.math :as math]))

(defn- deg [x] (* x c/degrees))

(defn- jd->mjd [jd] (- jd c/jd-mjd-offset))

(defn- series
  "Sum of `[coef args & {:cos :e}]` terms: coef (a number, or `[c0 c1]`
  for c0 + c1 T) times the sine -- or cosine -- of the linear combination
  `args` of the angles in `env`, times E per power given by `:e`."
  [rows env T E]
  (reduce (fn [s [coef args & flags]]
            (let [cos (some #{:cos} flags)
                  e   (second (drop-while #(not= :e %) flags))
                  a (reduce-kv (fn [a k m] (+ a (* m (env k)))) 0.0 args)
                  c (if (vector? coef) (+ (first coef) (* (second coef) T)) coef)]
              (+ s (* c (math/pow E (or e 0)) (if cos (math/cos a) (math/sin a))))))
          0.0 rows))

;; ------------------------------------------------------------- phases (49)

(def ^:private new-moon-terms
  [[-0.4072 {:Mp 1}]
   [0.17241 {:M 1} :e 1]
   [0.01608 {:Mp 2}]
   [0.01039 {:F 2}]
   [0.00739 {:Mp 1 :M -1} :e 1]
   [-0.00514 {:Mp 1 :M 1} :e 1]
   [0.00208 {:M 2} :e 2]
   [-0.00111 {:Mp 1 :F -2}]
   [-0.00057 {:Mp 1 :F 2}]
   [0.00056 {:Mp 2 :M 1} :e 1]
   [-0.00042 {:Mp 3}]
   [0.00042 {:M 1 :F 2} :e 1]
   [0.00038 {:M 1 :F -2} :e 1]
   [-0.00024 {:Mp 2 :M -1} :e 1]
   [-0.00017 {:Om 1}]
   [-7e-05 {:Mp 1 :M 2}]
   [4e-05 {:Mp 2 :F -2}]
   [4e-05 {:M 3}]
   [3e-05 {:Mp 1 :M 1 :F -2}]
   [3e-05 {:Mp 2 :F 2}]
   [-3e-05 {:Mp 1 :M 1 :F 2}]
   [3e-05 {:Mp 1 :M -1 :F 2}]
   [-2e-05 {:Mp 1 :M -1 :F -2}]
   [-2e-05 {:Mp 3 :M 1}]
   [2e-05 {:Mp 4}]])
(def ^:private full-moon-terms
  [[-0.40614 {:Mp 1}]
   [0.17302 {:M 1} :e 1]
   [0.01614 {:Mp 2}]
   [0.01043 {:F 2}]
   [0.00734 {:Mp 1 :M -1} :e 1]
   [-0.00515 {:Mp 1 :M 1} :e 1]
   [0.00209 {:M 2} :e 2]
   [-0.00111 {:Mp 1 :F -2}]
   [-0.00057 {:Mp 1 :F 2}]
   [0.00056 {:Mp 2 :M 1} :e 1]
   [-0.00042 {:Mp 3}]
   [0.00042 {:M 1 :F 2} :e 1]
   [0.00038 {:M 1 :F -2} :e 1]
   [-0.00024 {:Mp 2 :M -1} :e 1]
   [-0.00017 {:Om 1}]
   [-7e-05 {:Mp 1 :M 2}]
   [4e-05 {:Mp 2 :F -2}]
   [4e-05 {:M 3}]
   [3e-05 {:Mp 1 :M 1 :F -2}]
   [3e-05 {:Mp 2 :F 2}]
   [-3e-05 {:Mp 1 :M 1 :F 2}]
   [3e-05 {:Mp 1 :M -1 :F 2}]
   [-2e-05 {:Mp 1 :M -1 :F -2}]
   [-2e-05 {:Mp 3 :M 1}]
   [2e-05 {:Mp 4}]])
(def ^:private quarter-terms
  [[-0.62801 {:Mp 1}]
   [0.17172 {:M 1} :e 1]
   [-0.01183 {:Mp 1 :M 1} :e 1]
   [0.00862 {:Mp 2}]
   [0.00804 {:F 2}]
   [0.00454 {:Mp 1 :M -1} :e 1]
   [0.00204 {:M 2} :e 2]
   [-0.0018 {:Mp 1 :F -2}]
   [-0.0007 {:Mp 1 :F 2}]
   [-0.0004 {:Mp 3}]
   [-0.00034 {:Mp 2 :M -1}]
   [0.00032 {:M 1 :F 2} :e 1]
   [0.00032 {:M 1 :F -2} :e 1]
   [-0.00028 {:Mp 1 :M 2} :e 2]
   [0.00027 {:Mp 2 :M 1} :e 1]
   [-0.00017 {:Om 1}]
   [-5e-05 {:Mp 1 :M -1 :F -2}]
   [4e-05 {:Mp 2 :F 2}]
   [-4e-05 {:Mp 1 :M 1 :F 2}]
   [4e-05 {:Mp 1 :M -2}]
   [3e-05 {:Mp 1 :M 1 :F -2}]
   [3e-05 {:M 3}]
   [2e-05 {:Mp 2 :F -2}]
   [2e-05 {:Mp 1 :M -1 :F 2}]
   [-2e-05 {:Mp 3 :M 1}]])

(def ^:private planetary-coefficients
  "The fourteen planetary arguments' coefficients, days."
  [0.000325 0.000165 0.000164 0.000126 0.00011 6.2e-05 6e-05 5.6e-05 4.7e-05 4.2e-05 4e-05 3.7e-05 3.5e-05 2.3e-05])

(defn- phase-k [year q]
  (+ (math/floor (+ (- (* (- year 2000.0) 12.3685) q) 0.5)) q))

(defn- mean-phase-jd [k]
  (let [T (/ k 1236.85)]
    (+ 2451550.09766 (* 29.530588861 k)
       (* T T (horner T [0.00015437 -0.000000150 0.00000000073])))))

(defn mean-phase
  "MJD of the mean new moon (`phase` 0), first quarter (0.25), full moon
  (0.5) or last quarter (0.75) nearest `year` -- up to 14 hours from the
  true one."
  [year phase]
  (jd->mjd (mean-phase-jd (phase-k year phase))))

(defn moon-phase
  "MJD (TT) of the true new moon (`phase` 0), first quarter (0.25), full
  moon (0.5) or last quarter (0.75) nearest `year`, good to a few seconds
  over the centuries around the present."
  [year phase]
  (let [k  (phase-k year phase)
        T  (/ k 1236.85)
        E  (horner T [1.0 -0.002516 -0.0000074])
        env {:M  (deg (+ 2.5534 (* 29.10535670 k) (* T T (horner T [-0.0000014 -0.00000011]))))
             :Mp (deg (+ 201.5643 (* 385.81693528 k) (* T T (horner T [0.0107582 0.00001238 -0.000000058]))))
             :F  (deg (+ 160.7108 (* 390.67050284 k) (* T T (horner T [-0.0016118 -0.00000227 0.000000011]))))
             :Om (deg (+ 124.7746 (* -1.56375588 k) (* T T (horner T [0.0020672 0.00000215]))))}
        A  (map #(deg (+ (first %) (* (second %) k)))
                [[299.77 0.107408] [251.88 0.016321] [251.83 26.651886] [349.42 36.412478]
                 [84.66 18.206239] [141.74 53.303771] [207.17 2.453732] [154.84 7.306860]
                 [34.52 27.261239] [207.19 0.121824] [291.34 1.844379] [161.72 24.198154]
                 [239.56 25.513099] [331.55 3.592518]])
        A  (cons (- (first A) (deg (* 0.009173 T T))) (rest A))
        planetary (reduce + (map #(* %1 (math/sin %2)) planetary-coefficients A))
        {:keys [M Mp F]} env
        W  (- (+ 0.00306 (* -0.00038 E (math/cos M)) (* 0.00026 (math/cos Mp)))
              (* 0.00002 (- (math/cos (- Mp M)) (math/cos (+ Mp M)) (math/cos (* 2 F)))))
        correction (case (double phase)
                     0.0  (series new-moon-terms env T E)
                     0.5  (series full-moon-terms env T E)
                     0.25 (+ (series quarter-terms env T E) W)
                     0.75 (- (series quarter-terms env T E) W))]
    (jd->mjd (+ (mean-phase-jd k) correction planetary))))

;; -------------------------------------------------- perigee and apogee (50)

(def ^:private perigee-terms
  [[-1.6769 {:D 2}]
   [0.4589 {:D 4}]
   [-0.1856 {:D 6}]
   [0.0883 {:D 8}]
   [[-0.0773 0.00019] {:D 2 :M -1}]
   [[0.0502 -0.00013] {:M 1}]
   [-0.046 {:D 10}]
   [[0.0422 -0.00011] {:D 4 :M -1}]
   [-0.0256 {:D 6 :M -1}]
   [0.0253 {:D 12}]
   [0.0237 {:D 1}]
   [0.0162 {:D 8 :M -1}]
   [-0.0145 {:D 14}]
   [0.0129 {:F 2}]
   [-0.0112 {:D 3}]
   [-0.0104 {:D 10 :M -1}]
   [0.0086 {:D 16}]
   [0.0069 {:D 12 :M -1}]
   [0.0066 {:D 5}]
   [-0.0053 {:D 2 :F 2}]
   [-0.0052 {:D 18}]
   [-0.0046 {:D 14 :M -1}]
   [-0.0041 {:D 7}]
   [0.004 {:D 2 :M 1}]
   [0.0032 {:D 20}]
   [-0.0032 {:D 1 :M 1}]
   [0.0031 {:D 16 :M -1}]
   [-0.0029 {:D 4 :M 1}]
   [0.0027 {:D 9}]
   [0.0027 {:D 4 :F 2}]
   [-0.0027 {:D 2 :M -2}]
   [0.0024 {:D 4 :M -2}]
   [-0.0021 {:D 6 :M -2}]
   [-0.0021 {:D 22}]
   [-0.0021 {:D 18 :M -1}]
   [0.0019 {:D 6 :M 1}]
   [-0.0018 {:D 11}]
   [-0.0014 {:D 8 :M 1}]
   [-0.0014 {:D 4 :F -2}]
   [-0.0014 {:D 6 :F 2}]
   [0.0014 {:D 3 :M 1}]
   [-0.0014 {:D 5 :M 1}]
   [0.0013 {:D 13}]
   [0.0013 {:D 20 :M -1}]
   [0.0011 {:D 3 :M 2}]
   [-0.0011 {:D 4 :F 2 :M -2}]
   [-0.001 {:D 1 :M 2}]
   [-0.0009 {:D 22 :M -1}]
   [-0.0008 {:F 4}]
   [0.0008 {:D 6 :F -2}]
   [0.0008 {:D 2 :F -2 :M 1}]
   [0.0007 {:M 2}]
   [0.0007 {:F 2 :M -1}]
   [0.0007 {:D 2 :F 4}]
   [-0.0006 {:F 2 :M -2}]
   [-0.0006 {:D 2 :F -2 :M 2}]
   [0.0006 {:D 24}]
   [0.0005 {:D 4 :F -4}]
   [0.0005 {:D 2 :M 2}]
   [-0.0004 {:D 1 :M -1}]])
(def ^:private apogee-terms
  [[0.4392 {:D 2}]
   [0.0684 {:D 4}]
   [[0.0456 -0.00011] {:M 1}]
   [[0.0426 -0.00011] {:D 2 :M -1}]
   [0.0212 {:F 2}]
   [-0.0189 {:D 1}]
   [0.0144 {:D 6}]
   [0.0113 {:D 4 :M -1}]
   [0.0047 {:D 2 :F 2}]
   [0.0036 {:D 1 :M 1}]
   [0.0035 {:D 8}]
   [0.0034 {:D 6 :M -1}]
   [-0.0034 {:D 2 :F -2}]
   [0.0022 {:D 2 :M -2}]
   [-0.0017 {:D 3}]
   [0.0013 {:D 4 :F 2}]
   [0.0011 {:D 8 :M -1}]
   [0.001 {:D 4 :M -2}]
   [0.0009 {:D 10}]
   [0.0007 {:D 3 :M 1}]
   [0.0006 {:M 2}]
   [0.0005 {:D 2 :M 1}]
   [0.0005 {:D 2 :M 2}]
   [0.0004 {:D 6 :F 2}]
   [0.0004 {:D 6 :M -2}]
   [0.0004 {:D 10 :M -1}]
   [-0.0004 {:D 5}]
   [-0.0004 {:D 4 :F -2}]
   [0.0003 {:F 2 :M 1}]
   [0.0003 {:D 12}]
   [0.0003 {:D 2 :F 2 :M -1}]
   [-0.0003 {:D 1 :M -1}]])
(def ^:private perigee-parallax-terms
  [[63.224 {:D 2} :cos]
   [-6.99 {:D 4} :cos]
   [[2.834 -0.0071] {:D 2 :M -1} :cos]
   [1.927 {:D 6} :cos]
   [-1.263 {:D 1} :cos]
   [-0.702 {:D 8} :cos]
   [[0.696 -0.0017] {:M 1} :cos]
   [-0.69 {:F 2} :cos]
   [[-0.629 0.0016] {:D 4 :M -1} :cos]
   [-0.392 {:D 2 :F -2} :cos]
   [0.297 {:D 10} :cos]
   [0.26 {:D 6 :M -1} :cos]
   [0.201 {:D 3} :cos]
   [-0.161 {:D 2 :M 1} :cos]
   [0.157 {:D 1 :M 1} :cos]
   [-0.138 {:D 12} :cos]
   [-0.127 {:D 8 :M -1} :cos]
   [0.104 {:D 2 :F 2} :cos]
   [0.104 {:D 2 :M -2} :cos]
   [-0.079 {:D 5} :cos]
   [0.068 {:D 14} :cos]
   [0.067 {:D 10 :M -1} :cos]
   [0.054 {:D 4 :M 1} :cos]
   [-0.038 {:D 12 :M -1} :cos]
   [-0.038 {:D 4 :M -2} :cos]
   [0.037 {:D 7} :cos]
   [-0.037 {:D 4 :F 2} :cos]
   [-0.035 {:D 16} :cos]
   [-0.03 {:D 3 :M 1} :cos]
   [0.029 {:D 1 :M -1} :cos]
   [-0.025 {:D 6 :M 1} :cos]
   [0.023 {:M 2} :cos]
   [0.023 {:D 14 :M -1} :cos]
   [-0.023 {:D 2 :M 2} :cos]
   [0.022 {:D 6 :M -2} :cos]
   [-0.021 {:D 2 :F -2 :M -1} :cos]
   [-0.02 {:D 9} :cos]
   [0.019 {:D 18} :cos]
   [0.017 {:D 6 :F 2} :cos]
   [0.014 {:F 2 :M -1} :cos]
   [-0.014 {:D 16 :M -1} :cos]
   [0.013 {:D 4 :F -2} :cos]
   [0.012 {:D 8 :M 1} :cos]
   [0.011 {:D 11} :cos]
   [0.01 {:D 5 :M 1} :cos]
   [-0.01 {:D 20} :cos]])
(def ^:private apogee-parallax-terms
  [[-9.147 {:D 2} :cos]
   [-0.841 {:D 1} :cos]
   [0.697 {:F 2} :cos]
   [[-0.656 0.0016] {:M 1} :cos]
   [0.355 {:D 4} :cos]
   [0.159 {:D 2 :M -1} :cos]
   [0.127 {:D 1 :M 1} :cos]
   [0.065 {:D 4 :M -1} :cos]
   [0.052 {:D 6} :cos]
   [0.043 {:D 2 :M 1} :cos]
   [0.031 {:D 2 :F 2} :cos]
   [-0.023 {:D 2 :F -2} :cos]
   [0.022 {:D 2 :M -2} :cos]
   [0.019 {:D 2 :M 2} :cos]
   [-0.016 {:M 2} :cos]
   [0.014 {:D 6 :M -1} :cos]
   [0.01 {:D 8} :cos]])

(defn- apsis-args [year apogee?]
  (let [k (+ (math/floor (+ (- (* (- year 1999.97) 13.2555) (if apogee? 0.5 0.0)) 0.5))
             (if apogee? 0.5 0.0))
        T (/ k 1325.55)]
    {:k k :T T
     :jd (+ 2451534.6698 (* 27.55454989 k) (* T T (horner T [-0.0006691 -0.000001098 0.0000000052])))
     :env {:D (deg (+ 171.9179 (* 335.9106046 k) (* T T (horner T [-0.0100383 -0.00001156 0.000000055]))))
           :M (deg (+ 347.3477 (* 27.1577721 k) (* T T (horner T [-0.0008130 -0.0000010]))))
           :F (deg (+ 316.6109 (* 364.5287911 k) (* T T (horner T [-0.0125053 -0.0000148]))))}}))

(defn perigee
  "`[mjd parallax]` of the lunar perigee nearest `year`, and the Moon's
  equatorial horizontal parallax then -- from which its distance, 6378.14
  km / sin(parallax). The instant is good to half an hour; the
  departure from the mean is up to 1.5 days, from the Sun's pull on the
  orbit's line of apsides."
  [year]
  (let [{:keys [T jd env]} (apsis-args year false)]
    [(jd->mjd (+ jd (series perigee-terms env T 1.0)))
     (* c/arcsec (+ 3629.215 (series perigee-parallax-terms env T 1.0)))]))

(defn apogee
  "`[mjd parallax]` of the lunar apogee nearest `year`; see `perigee`."
  [year]
  (let [{:keys [T jd env]} (apsis-args year true)]
    [(jd->mjd (+ jd (series apogee-terms env T 1.0)))
     (* c/arcsec (+ 3245.251 (series apogee-parallax-terms env T 1.0)))]))

(defn mean-perigee-time
  "MJD of the mean perigee nearest `year`, or when `apogee?` the mean
  apogee."
  ([year] (mean-perigee-time year false))
  ([year apogee?] (jd->mjd (:jd (apsis-args year apogee?)))))

;; --------------------------------------------- passages through the nodes (51)

(def ^:private node-terms
  [[-0.4721 {:Mp 1}]
   [-0.1649 {:D 2}]
   [-0.0868 {:D 2 :Mp -1}]
   [0.0084 {:D 2 :Mp 1}]
   [-0.0083 {:D 2 :M -1} :e 1]
   [-0.0039 {:D 2 :M -1 :Mp -1} :e 1]
   [0.0034 {:Mp 2}]
   [-0.0031 {:D 2 :Mp -2}]
   [0.003 {:D 2 :M 1} :e 1]
   [0.0028 {:M 1 :Mp -1} :e 1]
   [0.0026 {:M 1} :e 1]
   [0.0025 {:D 4}]
   [0.0024 {:D 1}]
   [0.0022 {:M 1 :Mp 1} :e 1]
   [0.0017 {:Om 1}]
   [0.0014 {:D 4 :Mp -1}]
   [0.0005 {:D 2 :M 1 :Mp -1} :e 1]
   [0.0004 {:D 2 :M -1 :Mp 1} :e 1]
   [-0.0003 {:D 2 :M -2} :e 1]
   [0.0003 {:D 4 :M -1} :e 1]
   [0.0003 {:V 1}]
   [0.0003 {:P 1}]])

(defn node-passage
  "MJD (TT) of the Moon's passage through the ascending node of its orbit
  nearest `year`, or when `descending?` the descending node -- the only
  places an eclipse can happen."
  ([year] (node-passage year false))
  ([year descending?]
   (let [h (if descending? 0.5 0.0)
         k (+ (math/floor (+ (- (* (- year 2000.05) 13.4223) h) 0.5)) h)
         T (/ k 1342.23)
         E (horner T [1.0 -0.002516 -0.0000074])
         om (deg (+ 123.9767 (* -1.44098956 k) (* T T (horner T [0.0020608 0.00000214 -0.000000016]))))
         env {:D  (deg (+ 183.6380 (* 331.73735682 k) (* T T (horner T [0.0014852 0.00000209 -0.00000001]))))
              :M  (deg (+ 17.4006 (* 26.8203725 k) (* T T (horner T [0.0001186 0.00000006]))))
              :Mp (deg (+ 38.3776 (* 355.52747313 k) (* T T (horner T [0.0123499 0.000014627 -0.000000069]))))
              :Om om
              :V  (deg (+ 299.75 (* 132.85 T) (* -0.009173 T T)))
              :P  (+ om (deg (- 272.75 (* 2.3 T))))}]
     (jd->mjd (+ 2451565.1619 (* 27.212220817 k) (* T T (horner T [0.0002762 0.000000021 -0.000000000088]))
                 (series node-terms env T E))))))

;; ------------------------------------------ greatest declinations (52)

(def ^:private north-time-terms
  [[0.8975 {:F 1} :cos]
   [-0.4726 {:Mp 1}]
   [-0.103 {:F 2}]
   [-0.0976 {:D 2 :Mp -1}]
   [-0.0462 {:Mp 1 :F -1} :cos]
   [-0.0461 {:Mp 1 :F 1} :cos]
   [-0.0438 {:D 2}]
   [0.0162 {:M 1} :e 1]
   [-0.0157 {:F 3} :cos]
   [0.0145 {:Mp 1 :F 2}]
   [0.0136 {:D 2 :F -1} :cos]
   [-0.0095 {:D 2 :Mp -1 :F -1} :cos]
   [-0.0091 {:D 2 :Mp -1 :F 1} :cos]
   [-0.0089 {:D 2 :F 1} :cos]
   [0.0075 {:Mp 2}]
   [-0.0068 {:Mp 1 :F -2}]
   [0.0061 {:Mp 2 :F -1} :cos]
   [-0.0047 {:Mp 1 :F 3}]
   [-0.0043 {:D 2 :M -1 :Mp -1} :e 1]
   [-0.004 {:Mp 1 :F -2} :cos]
   [-0.0037 {:D 2 :Mp -2}]
   [0.0031 {:F 1}]
   [0.003 {:D 2 :Mp 1}]
   [-0.0029 {:Mp 1 :F 2} :cos]
   [-0.0029 {:D 2 :M -1} :e 1]
   [-0.0027 {:Mp 1 :F 1}]
   [0.0024 {:M 1 :Mp -1} :e 1]
   [-0.0021 {:Mp 1 :F -3}]
   [0.0019 {:Mp 2 :F 1}]
   [0.0018 {:D 2 :Mp -2 :F -1} :cos]
   [0.0018 {:F 3}]
   [0.0017 {:Mp 1 :F 3} :cos]
   [0.0017 {:Mp 2} :cos]
   [-0.0014 {:D 2 :Mp -1} :cos]
   [0.0013 {:D 2 :Mp 1 :F 1} :cos]
   [0.0013 {:Mp 1} :cos]
   [0.0012 {:Mp 3 :F 1}]
   [0.0011 {:D 2 :Mp -1 :F 1}]
   [-0.0011 {:D 2 :Mp -2} :cos]
   [0.001 {:D 1 :F 1} :cos]
   [0.001 {:M 1 :Mp 1} :e 1]
   [-0.0009 {:D 2 :F -2}]
   [0.0007 {:Mp 2 :F 1} :cos]
   [-0.0007 {:Mp 3 :F 1} :cos]])
(def ^:private south-time-terms
  [[-0.8975 {:F 1} :cos]
   [-0.4726 {:Mp 1}]
   [-0.103 {:F 2}]
   [-0.0976 {:D 2 :Mp -1}]
   [0.0541 {:Mp 1 :F -1} :cos]
   [0.0516 {:Mp 1 :F 1} :cos]
   [-0.0438 {:D 2}]
   [0.0112 {:M 1} :e 1]
   [0.0157 {:F 3} :cos]
   [0.0023 {:Mp 1 :F 2}]
   [-0.0136 {:D 2 :F -1} :cos]
   [0.011 {:D 2 :Mp -1 :F -1} :cos]
   [0.0091 {:D 2 :Mp -1 :F 1} :cos]
   [0.0089 {:D 2 :F 1} :cos]
   [0.0075 {:Mp 2}]
   [-0.003 {:Mp 1 :F -2}]
   [-0.0061 {:Mp 2 :F -1} :cos]
   [-0.0047 {:Mp 1 :F 3}]
   [-0.0043 {:D 2 :M -1 :Mp -1} :e 1]
   [0.004 {:Mp 1 :F -2} :cos]
   [-0.0037 {:D 2 :Mp -2}]
   [-0.0031 {:F 1}]
   [0.003 {:D 2 :Mp 1}]
   [0.0029 {:Mp 1 :F 2} :cos]
   [-0.0029 {:D 2 :M -1} :e 1]
   [-0.0027 {:Mp 1 :F 1}]
   [0.0024 {:M 1 :Mp -1} :e 1]
   [-0.0021 {:Mp 1 :F -3}]
   [-0.0019 {:Mp 2 :F 1}]
   [-0.0006 {:D 2 :Mp -2 :F -1} :cos]
   [-0.0018 {:F 3}]
   [-0.0017 {:Mp 1 :F 3} :cos]
   [0.0017 {:Mp 2} :cos]
   [0.0014 {:D 2 :Mp -1} :cos]
   [-0.0013 {:D 2 :Mp 1 :F 1} :cos]
   [-0.0013 {:Mp 1} :cos]
   [0.0012 {:Mp 3 :F 1}]
   [0.0011 {:D 2 :Mp -1 :F 1}]
   [0.0011 {:D 2 :Mp -2} :cos]
   [0.001 {:D 1 :F 1} :cos]
   [0.001 {:M 1 :Mp 1} :e 1]
   [-0.0009 {:D 2 :F -2}]
   [-0.0007 {:Mp 2 :F 1} :cos]
   [-0.0007 {:Mp 3 :F 1} :cos]])
(def ^:private north-dec-terms
  [[5.1093 {:F 1}]
   [0.2658 {:F 2} :cos]
   [0.1448 {:D 2 :F -1}]
   [-0.0322 {:F 3}]
   [0.0133 {:D 2 :F -2} :cos]
   [0.0125 {:D 2} :cos]
   [-0.0124 {:Mp 1 :F -1}]
   [-0.0101 {:Mp 1 :F 2}]
   [0.0097 {:F 1} :cos]
   [-0.0087 {:D 2 :M 1 :F -1} :e 1]
   [0.0074 {:Mp 1 :F 3}]
   [0.0067 {:D 1 :F 1}]
   [0.0063 {:Mp 1 :F -2}]
   [0.006 {:D 2 :M -1 :F -1} :e 1]
   [-0.0057 {:D 2 :Mp -1 :F -1}]
   [-0.0056 {:Mp 1 :F 1} :cos]
   [0.0052 {:Mp 1 :F 2} :cos]
   [0.0041 {:Mp 2 :F 1} :cos]
   [-0.004 {:Mp 1 :F -3} :cos]
   [0.0038 {:Mp 2 :F -1} :cos]
   [-0.0034 {:Mp 1 :F -2} :cos]
   [-0.0029 {:Mp 2}]
   [0.0029 {:Mp 3 :F 1}]
   [-0.0028 {:D 2 :M 1 :F -1} :cos :e 1]
   [-0.0028 {:Mp 1 :F -1} :cos]
   [-0.0023 {:F 3} :cos]
   [-0.0021 {:D 2 :F 1}]
   [0.0019 {:Mp 1 :F 3} :cos]
   [0.0018 {:D 1 :F 1} :cos]
   [0.0017 {:Mp 2 :F -1}]
   [0.0015 {:Mp 3 :F 1} :cos]
   [0.0014 {:D 2 :Mp 2 :F 1} :cos]
   [-0.0012 {:D 2 :Mp -2 :F -1}]
   [-0.0012 {:Mp 2} :cos]
   [-0.001 {:Mp 1} :cos]
   [-0.001 {:F 2}]
   [0.0006 {:Mp 1 :F 1}]])
(def ^:private south-dec-terms
  [[-5.1093 {:F 1}]
   [0.2658 {:F 2} :cos]
   [-0.1448 {:D 2 :F -1}]
   [0.0322 {:F 3}]
   [0.0133 {:D 2 :F -2} :cos]
   [0.0125 {:D 2} :cos]
   [-0.0015 {:Mp 1 :F -1}]
   [0.0101 {:Mp 1 :F 2}]
   [-0.0097 {:F 1} :cos]
   [0.0087 {:D 2 :M 1 :F -1} :e 1]
   [0.0074 {:Mp 1 :F 3}]
   [0.0067 {:D 1 :F 1}]
   [-0.0063 {:Mp 1 :F -2}]
   [-0.006 {:D 2 :M -1 :F -1} :e 1]
   [0.0057 {:D 2 :Mp -1 :F -1}]
   [-0.0056 {:Mp 1 :F 1} :cos]
   [-0.0052 {:Mp 1 :F 2} :cos]
   [-0.0041 {:Mp 2 :F 1} :cos]
   [-0.004 {:Mp 1 :F -3} :cos]
   [-0.0038 {:Mp 2 :F -1} :cos]
   [0.0034 {:Mp 1 :F -2} :cos]
   [-0.0029 {:Mp 2}]
   [0.0029 {:Mp 3 :F 1}]
   [0.0028 {:D 2 :M 1 :F -1} :cos :e 1]
   [-0.0028 {:Mp 1 :F -1} :cos]
   [0.0023 {:F 3} :cos]
   [0.0021 {:D 2 :F 1}]
   [0.0019 {:Mp 1 :F 3} :cos]
   [0.0018 {:D 1 :F 1} :cos]
   [-0.0017 {:Mp 2 :F -1}]
   [0.0015 {:Mp 3 :F 1} :cos]
   [0.0014 {:D 2 :Mp 2 :F 1} :cos]
   [0.0012 {:D 2 :Mp -2 :F -1}]
   [-0.0012 {:Mp 2} :cos]
   [0.001 {:Mp 1} :cos]
   [-0.001 {:F 2}]
   [0.0037 {:Mp 1 :F 1}]])

(defn greatest-declination
  "`[mjd dec]` of the Moon's greatest northern declination nearest `year`,
  or when `south?` its greatest southern one. The extreme swings between
  18.3 and 28.6 degrees over the 18.6-year cycle of the node: the major
  and minor lunar standstills that some megalithic sites are aligned on."
  ([year] (greatest-declination year false))
  ([year south?]
   (let [k (math/floor (+ (* (- year 2000.03) 13.3686) 0.5))
         T (/ k 1336.86)
         E (horner T [1.0 -0.002516 -0.0000074])
         [D0 M0 Mp0 F0 jd0] (if south?
                              [345.6676 1.3951 186.2100 145.1633 2451548.9289]
                              [152.2029 14.8591 4.6881 325.8867 2451562.5897])
         env {:D  (deg (+ D0 (* 333.0705546 k) (* T T (horner T [-0.0004214 0.00000011]))))
              :M  (deg (+ M0 (* 26.9281592 k) (* T T (horner T [-0.0000355 -0.0000001]))))
              :Mp (deg (+ Mp0 (* 356.9562794 k) (* T T (horner T [0.0103066 0.00001251]))))
              :F  (deg (+ F0 (* 1.4467807 k) (* T T (horner T [-0.0020690 -0.00000215]))))}
         jd (+ jd0 (* 27.321582247 k) (* T T (horner T [0.000119804 -0.000000141]))
               (series (if south? south-time-terms north-time-terms) env T E))
         dec (+ 23.6961 (* -0.013004 T) (series (if south? south-dec-terms north-dec-terms) env T E))]
     [(jd->mjd jd) (deg (if south? (- dec) dec))])))
