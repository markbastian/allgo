(ns allgo.astro.apparent
  "The apparent place of a star: where it is actually seen, once its
  catalog position has been carried through proper motion, precession,
  nutation and aberration (Meeus, *Astronomical Algorithms*, chapter 23).

  Aberration is the largest of the corrections a catalog does not carry --
  up to 20.5 arcseconds, the ratio of the Earth's orbital speed to the
  speed of light. Light arriving from a star is seen tilted toward the
  direction the Earth is moving, as rain falls slanted on a moving car, and
  over a year every star traces a small ellipse about its mean place.

  Positions are `[ra dec]` or `[lon lat]` in radians; times MJD (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.precession :as precession]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [allgo.math :as am]
            [allgo.numerics.interpolation :refer [horner]]
            [clojure.math :as math]))

(def aberration-constant
  "The constant of aberration, kappa: 20.49552 arcseconds."
  (* 20.49552 c/arcsec))

(defn- perihelion-longitude
  "Longitude of the Earth's perihelion, which enters because an elliptic
  orbit adds a small fixed term to the aberration, the E-terms."
  [mjd-tt]
  (* c/degrees (horner (time/centuries-J2000 mjd-tt) [102.93735 1.71946 0.00046])))

(defn nutation
  "`[dra ddec]`, nutation's shift of a star's right ascension and
  declination on the mean equator of date (Meeus 23.1)."
  [[ra dec] mjd-tt]
  (let [eps (frames/mean-obliquity mjd-tt)
        [dpsi deps] (frames/nutation-angles mjd-tt)
        se (math/sin eps) ce (math/cos eps)
        sa (math/sin ra) ca (math/cos ra)
        td (math/tan dec)]
    [(- (* (+ ce (* se sa td)) dpsi) (* ca td deps))
     (+ (* se ca dpsi) (* sa deps))]))

(defn aberration
  "`[dra ddec]`, annual aberration (Meeus 23.3), the classical formula in
  the Sun's longitude and the eccentricity of the Earth's orbit."
  [[ra dec] mjd-tt]
  (let [eps (frames/mean-obliquity mjd-tt)
        s   (solar/true-longitude mjd-tt)
        e   (solar/eccentricity mjd-tt)
        p   (perihelion-longitude mjd-tt)
        sa (math/sin ra) ca (math/cos ra)
        sd (math/sin dec) cd (math/cos dec)
        ss (math/sin s) cs (math/cos s)
        sp (math/sin p) cp (math/cos p)
        ce (math/cos eps) te (math/tan eps)
        k aberration-constant
        q1 (* ca ce)
        q2 (* ce (- (* te cd) (* sa sd)))
        q3 (* ca sd)]
    [(/ (* k (- (* e (+ (* q1 cp) (* sa sp))) (+ (* q1 cs) (* sa ss)))) cd)
     (* k (- (* e (+ (* cp q2) (* sp q3))) (+ (* cs q2) (* ss q3))))]))

(defn ecliptic-aberration
  "`[dlon dlat]`, annual aberration in ecliptic coordinates (Meeus 23.2)."
  [[lon lat] mjd-tt]
  (let [s (solar/true-longitude mjd-tt)
        e (solar/eccentricity mjd-tt)
        p (perihelion-longitude mjd-tt)
        k aberration-constant]
    [(/ (* k (- (* e (math/cos (- p lon))) (math/cos (- s lon)))) (math/cos lat))
     (* (- k) (math/sin lat) (- (math/sin (- s lon)) (* e (math/sin (- p lon)))))]))

(defn apparent-place
  "Apparent `[ra dec]` at `mjd-tt` of a star at catalog position `pos` for
  the epoch and equinox `from` with proper motion `pm` (`[pm-ra pm-dec]`,
  radians a Julian year): moved, precessed, then corrected for nutation and
  aberration. Good to about 0.01 arcsecond, the size of the terms the
  classical aberration formula leaves out."
  ([pos from mjd-tt] (apparent-place pos from mjd-tt [0.0 0.0]))
  ([pos from mjd-tt pm]
   (let [[ra dec :as mean] (precession/equatorial pos from mjd-tt pm)
         [a1 d1] (nutation mean mjd-tt)
         [a2 d2] (aberration mean mjd-tt)]
     [(am/wrap-2pi (+ ra a1 a2)) (+ dec d1 d2)])))

;; ----------------------------------------------------------- Ron-Vondrak

(def ^:private ron-vondrak
  "The 36 terms of Ron and Vondrák's expansion of the Earth's velocity
  (Meeus table 23.A). Each row is the argument, as multipliers of
  `[L2 L3 L4 L5 L6 L7 L8 L' D M' F]`, then the sine and cosine
  coefficients of X, Y and Z, each `[c0 c1]` for c0 + c1 T, in units of
  1e-8 AU a day."
  [[[0 1 0 0 0 0 0 0 0 0 0] [-1719914 -2] [-25 0] [25 -13] [1578089 156] [10 32] [684185 -358]]
   [[0 2 0 0 0 0 0 0 0 0 0] [6434 141] [28007 -107] [25697 -95] [-5904 -130] [11141 -48] [-2559 -55]]
   [[0 0 0 1 0 0 0 0 0 0 0] [715 0] [0 0] [6 0] [-657 0] [-15 0] [-282 0]]
   [[0 0 0 0 0 0 0 1 0 0 0] [715 0] [0 0] [0 0] [-656 0] [0 0] [-285 0]]
   [[0 3 0 0 0 0 0 0 0 0 0] [486 -5] [-236 -4] [-216 -4] [-446 5] [-94 0] [-193 0]]
   [[0 0 0 0 1 0 0 0 0 0 0] [159 0] [0 0] [2 0] [-147 0] [-6 0] [-61 0]]
   [[0 0 0 0 0 0 0 0 0 0 1] [0 0] [0 0] [0 0] [26 0] [0 0] [-59 0]]
   [[0 0 0 0 0 0 0 1 0 1 0] [39 0] [0 0] [0 0] [-36 0] [0 0] [-16 0]]
   [[0 0 0 2 0 0 0 0 0 0 0] [33 0] [-10 0] [-9 0] [-30 0] [-5 0] [-13 0]]
   [[0 2 0 -1 0 0 0 0 0 0 0] [31 0] [1 0] [1 0] [-28 0] [0 0] [-12 0]]
   [[0 3 -8 3 0 0 0 0 0 0 0] [8 0] [-28 0] [25 0] [8 0] [11 0] [3 0]]
   [[0 5 -8 3 0 0 0 0 0 0 0] [8 0] [-28 0] [-25 0] [-8 0] [-11 0] [-3 0]]
   [[2 -1 0 0 0 0 0 0 0 0 0] [21 0] [0 0] [0 0] [-19 0] [0 0] [-8 0]]
   [[1 0 0 0 0 0 0 0 0 0 0] [-19 0] [0 0] [0 0] [17 0] [0 0] [8 0]]
   [[0 0 0 0 0 1 0 0 0 0 0] [17 0] [0 0] [0 0] [-16 0] [0 0] [-7 0]]
   [[0 1 0 -2 0 0 0 0 0 0 0] [16 0] [0 0] [0 0] [15 0] [1 0] [7 0]]
   [[0 0 0 0 0 0 1 0 0 0 0] [16 0] [0 0] [1 0] [-15 0] [-3 0] [-6 0]]
   [[0 1 0 1 0 0 0 0 0 0 0] [11 0] [-1 0] [-1 0] [-10 0] [-1 0] [-5 0]]
   [[2 -2 0 0 0 0 0 0 0 0 0] [0 0] [-11 0] [-10 0] [0 0] [-4 0] [0 0]]
   [[0 1 0 -1 0 0 0 0 0 0 0] [-11 0] [-2 0] [-2 0] [9 0] [-1 0] [4 0]]
   [[0 4 0 0 0 0 0 0 0 0 0] [-7 0] [-8 0] [-8 0] [6 0] [-3 0] [3 0]]
   [[0 3 0 -2 0 0 0 0 0 0 0] [-10 0] [0 0] [0 0] [9 0] [0 0] [4 0]]
   [[1 -2 0 0 0 0 0 0 0 0 0] [-9 0] [0 0] [0 0] [-9 0] [0 0] [-4 0]]
   [[2 -3 0 0 0 0 0 0 0 0 0] [-9 0] [0 0] [0 0] [-8 0] [0 0] [-4 0]]
   [[0 0 0 0 2 0 0 0 0 0 0] [0 0] [-9 0] [-8 0] [0 0] [-3 0] [0 0]]
   [[2 -4 0 0 0 0 0 0 0 0 0] [0 0] [-9 0] [8 0] [0 0] [3 0] [0 0]]
   [[0 3 -2 0 0 0 0 0 0 0 0] [8 0] [0 0] [0 0] [-8 0] [0 0] [-3 0]]
   [[0 0 0 0 0 0 0 1 2 -1 0] [8 0] [0 0] [0 0] [-7 0] [0 0] [-3 0]]
   [[8 -12 0 0 0 0 0 0 0 0 0] [-4 0] [-7 0] [-6 0] [4 0] [-3 0] [2 0]]
   [[8 -14 0 0 0 0 0 0 0 0 0] [-4 0] [-7 0] [6 0] [-4 0] [3 0] [-2 0]]
   [[0 0 2 0 0 0 0 0 0 0 0] [-6 0] [-5 0] [-4 0] [5 0] [-2 0] [2 0]]
   [[3 -4 0 0 0 0 0 0 0 0 0] [-1 0] [-1 0] [-2 0] [-7 0] [1 0] [-4 0]]
   [[0 2 0 -2 0 0 0 0 0 0 0] [4 0] [-6 0] [-5 0] [-4 0] [-2 0] [-2 0]]
   [[3 -3 0 0 0 0 0 0 0 0 0] [0 0] [-7 0] [-6 0] [0 0] [-3 0] [0 0]]
   [[0 2 -2 0 0 0 0 0 0 0 0] [5 0] [-5 0] [-4 0] [-5 0] [-2 0] [-2 0]]
   [[0 0 0 0 0 0 0 1 -2 0 0] [5 0] [0 0] [0 0] [-5 0] [0 0] [-2 0]]])

(def ^:private light-speed-au-day
  "The speed of light in 1e-8 AU a day, the unit of the table above."
  17314463350.0)

(defn aberration-ron-vondrak
  "`[dra ddec]`, annual aberration from the Earth's velocity as Ron and
  Vondrák expand it -- the planetary perturbations of that velocity
  included, so good to 0.001 arcsecond rather than the 0.01 of the
  classical formula. The position is referred to the mean equinox of
  J2000, as the expansion is."
  [[ra dec] mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        args [(+ 3.1761467 (* 1021.3285546 T)) (+ 1.7534703 (* 628.3075849 T))
              (+ 6.2034809 (* 334.0612431 T)) (+ 0.5995465 (* 52.9690965 T))
              (+ 0.8740168 (* 21.3299095 T)) (+ 5.4812939 (* 7.4781599 T))
              (+ 5.3118863 (* 3.8133036 T)) (+ 3.8103444 (* 8399.6847337 T))
              (+ 5.1984667 (* 7771.3771486 T)) (+ 2.3555559 (* 8328.6914289 T))
              (+ 1.6279052 (* 8433.4661601 T))]
        coef (fn [[c0 c1]] (+ c0 (* c1 T)))
        [X Y Z] (reduce (fn [[X Y Z] [mult xs xc ys yc zs zc]]
                          (let [A  (reduce + (map * mult args))
                                sA (math/sin A) cA (math/cos A)]
                            [(+ X (* (coef xs) sA) (* (coef xc) cA))
                             (+ Y (* (coef ys) sA) (* (coef yc) cA))
                             (+ Z (* (coef zs) sA) (* (coef zc) cA))]))
                        [0.0 0.0 0.0] ron-vondrak)
        sa (math/sin ra) ca (math/cos ra)
        sd (math/sin dec) cd (math/cos dec)
        c  light-speed-au-day]
    [(/ (- (* Y ca) (* X sa)) (* c cd))
     (- (/ (- (* (+ (* X ca) (* Y sa)) sd) (* Z cd)) c))]))

(defn apparent-place-ron-vondrak
  "`apparent-place` with the Ron-Vondrák aberration, applied as it must be
  in the J2000 frame -- before precession rather than after -- for a star
  catalogued at J2000."
  ([pos mjd-tt] (apparent-place-ron-vondrak pos mjd-tt [0.0 0.0]))
  ([[ra dec] mjd-tt [pm-ra pm-dec]]
   (let [years (/ (- mjd-tt c/mjd-J2000) 365.25)
         moved [(+ ra (* years pm-ra)) (+ dec (* years pm-dec))]
         [da dd] (aberration-ron-vondrak moved mjd-tt)
         [ra dec :as mean] (precession/equatorial [(+ (first moved) da) (+ (second moved) dd)]
                                                  c/mjd-J2000 mjd-tt)
         [a1 d1] (nutation mean mjd-tt)]
     [(am/wrap-2pi (+ ra a1)) (+ dec d1)])))
