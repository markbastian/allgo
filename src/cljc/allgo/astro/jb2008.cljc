(ns allgo.astro.jb2008
  "The Jacchia-Bowman 2008 thermospheric density model, implemented from
  its published description: Bowman, Tobiska, Marcos, Huang, Lin and
  Burke, \"A new empirical thermospheric density model JB2008 using new
  solar and geomagnetic indices\", AIAA 2008-6438 (a work of the US
  government), and, for what JB2008 keeps of its predecessor, Bowman,
  Tobiska, Marcos and Valladares, \"The JB2006 empirical thermospheric
  density model\", J. Atmos. Solar-Terr. Phys. 70, 774-793 (2008). No
  code of the model's distributors was used; see `limits` below.

  JB2008 is Jacchia's 1970 model (`allgo.astro.jacchia`, `j70`) with its
  inputs replaced:

  - The nighttime minimum exospheric temperature from four solar indices
    (JB2008 eqs. 1-2) -- F10.7, S10 (EUV, 26-34 nm), M10 (Mg II core-to-
    wing, FUV) and Y10 (X-ray and Lyman-alpha) -- daily values and their
    81-day centered averages, each daily value lagged as the paper found
    best: F10 and S10 one day, M10 two, Y10 five.
  - Corrections to that temperature for local solar time and latitude
    (JB2006 eqs. 7-8, Table 4) at and above 200 km, and to the inflection
    temperature below (eq. 9).
  - A semiannual density variation whose amplitude F(z) and phase G(t)
    follow the solar indices (JB2008 eqs. 3-7, Tables 1-2).
  - A geomagnetic temperature rise, dTc, given as an input: from the Dst
    index through a storm (JB2008 eqs. 8-13, `storm-dtc`) and from ap
    otherwise (Jacchia 1970 eq. 22, `ap-dtc`).
  - A density factor above 1000 km (JB2006 eqs. 10-11, Table 5).

  The rest is Jacchia 1970's: the static profiles, the diurnal variation
  (eqs. 15-17), the seasonal-latitudinal variations of the lower
  thermosphere (24) and of helium (25).

  Heights km, latitudes and declinations degrees, times MJD (UT)."
  (:require [allgo.astro.jacchia :as j]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [clojure.math :as math]))

(def limits
  "Where the published description leaves a choice, and what was chosen:

  - The diurnal variation's R, which Jacchia 1970 lets vary from 0.27 to
    0.4: his recommended average, 0.31.
  - JB2006's local-time corrections are fitted for 700 > z >= 250 km (eq.
    7) and 200 >= z >= 140 km (eq. 9): above 700 km eq. 7 is held at its
    700 km value, and below 140 km eq. 9 at its 140 km value.
  - F in those corrections is (F10 - 100)/100: the 81-day average is
    used, the index the fits were binned by (the papers' figures).
  - JB2008 uses the 81-day centered averages in place of the July
    averages in F(z) (JB2008, after eq. 5); G(t)'s phase is omega = 2 pi
    (t - 1)/365, t the day of the year (JB2006 eq. 6).
  - The high-altitude factor from 1000 to 1500 km is the cubic with
    Table 5's value and slope at 1500 km and value 1, slope 0 at 1000 km,
    as the text specifies it (JB2006 eq. 11 as printed scales the slope
    by 500 twice).
  - Hydrogen, as in `allgo.astro.jacchia/j70`, from its density at 500 km
    fixed by the exospheric temperature.
  - The storm algorithm's event detection -- start, Dst minimum, recovery
    slope change and end -- is described but not specified (its end-time
    fit is not given), so `storm-dtc` takes the event times as inputs.")

;; ---------------------------------------------------------- temperature

(defn night-minimum
  "Tc, the global nighttime minimum exospheric temperature, K (JB2008
  eqs. 1-2), from the daily solar indices `f10` `s10` `m10` `y10` (each
  lagged as the namespace says) and their 81-day centered averages
  `f10b` `s10b` `m10b` `y10b`."
  [{:keys [f10 f10b s10 s10b m10 m10b y10 y10b]}]
  (let [w (math/pow (/ f10b 240.0) 0.25)
        fs (+ (* f10b w) (* s10b (- 1.0 w)))]
    (+ 392.4 (* 3.227 fs) (* 0.298 (- f10 f10b)) (* 2.259 (- s10 s10b))
       (* 0.312 (- m10 m10b)) (* 0.178 (- y10 y10b)))))

(def ^:private b-table
  [-0.457512297e+01 -0.512114909e+01 -0.693003609e+02 0.203716701e+03 0.703316291e+03
   -0.194349234e+04 0.110651308e+04 -0.174378996e+03 0.188594601e+04 -0.709371517e+04
   0.922454523e+04 -0.384508073e+04 -0.645841789e+01 0.409703319e+02 -0.482006560e+03
   0.181870931e+04 -0.237389204e+04 0.996703815e+03 0.361416936e+02])

(def ^:private c-table
  [-0.155986211e+02 -0.512114909e+01 -0.693003609e+02 0.203716701e+03 0.703316291e+03
   -0.194349234e+04 0.110651308e+04 -0.220835117e+03 0.143256989e+04 -0.318481844e+04
   0.328981513e+04 -0.135332119e+04 0.199956489e+02 -0.127093998e+02 0.212825156e+02
   -0.275555432e+01 0.110234982e+02 0.148881951e+03 -0.751640284e+03 0.637876542e+03
   0.127093998e+02 -0.212825156e+02 0.275555432e+01])

(def ^:private d-table
  [0.828727388e+00 0.124730376e+02 -0.318077422e+03 0.645905237e+03 -0.322577619e+03
   0.709584573e+01 0.126347673e+02 -0.744574342e+02 0.511241177e+02 -0.354195845e+01
   0.225275051e+02 -0.223771015e+02])

(defn- poly
  "c_k th^0 + c_(k+1) th^1 + ... for the coefficients `cs` (1-based) from
  k, `n` of them, starting at power `p0`."
  [cs k n p0 th]
  (reduce + (for [i (range n)] (* (cs (+ k i -1)) (math/pow th (+ p0 i))))))

(defn tc-correction
  "dTc, K, the correction to the nighttime minimum temperature for local
  solar time `lst` (hours) and latitude `lat` (deg) at height `z` km and
  81-day average flux `f10b`: JB2006 eq. 7 from 250 km up, eq. 8 from 200
  to 250 km; 0 below 200."
  [z lst lat f10b]
  (let [F (/ (- f10b 100.0) 100.0)
        th (/ lst 24.0)
        ph (math/cos (math/to-radians lat))
        B (fn [k] (b-table (dec k)))
        C (fn [k] (c-table (dec k)))]
    (cond
      (>= z 250.0)
      (let [H (/ (min z 700.0) 100.0)]
        (+ (B 1)
           (* F (poly b-table 2 6 0 th))
           (* ph (poly b-table 8 5 1 th))
           (* ph H (poly b-table 13 6 0 th))
           (* (B 19) ph)))
      (>= z 200.0)
      (let [H (/ (- z 200.0) 50.0)]
        (+ (* H (C 1))
           (* H F (poly c-table 2 6 0 th))
           (* H ph (+ (poly c-table 8 5 1 th) (C 13) (* (C 14) F) (* (C 15) F th) (* (C 16) F th th)))
           (C 17)
           (* ph (+ (poly c-table 18 3 1 th) (* (C 21) F) (* (C 22) F th) (* (C 23) F th th)))))
      :else 0.0)))

(defn tx-correction
  "dTx, K, the correction to the inflection temperature for local solar
  time `lst` (hours) at height `z` km and 81-day average flux `f10b`:
  JB2006 eq. 9, from 140 to 200 km; 0 at and above 200."
  [z lst f10b]
  (if (>= z 200.0)
    0.0
    (let [F (/ (- f10b 100.0) 100.0)
          th (/ lst 24.0)
          H (/ (max z 140.0) 100.0)
          D (fn [k] (d-table (dec k)))]
      (+ (D 1)
         (poly d-table 2 4 1 th)
         (* H (poly d-table 6 4 0 th))
         (* F (poly d-table 10 3 0 th))))))

;; --------------------------------------------------------- geomagnetic

(defn ap-dtc
  "The geomagnetic temperature rise, K, outside storms: Jacchia 1970 eq.
  22 from the 3-hour ap (lagged 6.7 hours), ap limited to 50 when no Dst
  storm is under way, as JB2008 limits it."
  [ap]
  (let [ap (min ap 50.0)]
    (+ ap (* 100.0 (- 1.0 (math/exp (* -0.08 ap)))))))

(defn main-phase-slope
  "S of JB2008 eq. 10 for a storm of minimum Dst `dst-min`, -1.40 below
  -450 nT."
  [dst-min]
  (if (< dst-min -450.0)
    -1.40
    (+ (* -1.5050e-5 dst-min dst-min) (* -1.0604e-2 dst-min) -3.20)))

(defn storm-dtc
  "The geomagnetic temperature rise, K, hour by hour through a storm, from
  hourly `dst` (nT, the first at the storm's start) and the hours of its
  minimum `i-min` and recovery slope change `i-change` counted from the
  start, beginning from `dtc0` (JB2008: Jacchia 1970's value from the ap
  at the start): the main phase by eq. 9 with eq. 10's slope, eq. 11
  wherever Dst rises, and Dst lagged 0, 1 or 2 hours for large, moderate
  and minor storms; the early recovery by eq. 12; the late by eq. 13 with
  S = -2.5; never below 0."
  [dst i-min i-change dtc0]
  (let [dst (vec dst)
        dmin (dst i-min)
        s (main-phase-slope dmin)
        lag (cond (< dmin -350.0) 0 (< dmin -250.0) 1 :else 2)
        at (fn [i] (dst (max 0 i)))]
    (reductions
     (fn [dtc i]
       (let [next (cond
                    (<= i i-min)
                    (let [d0 (at (- i 1 lag)) d1 (at (- i lag))]
                      (if (> d1 d0)
                        (- dtc (* 0.3 s (- d1 d0)))                              ; eq. 11
                        (+ (* 0.846 dtc) (* s (- d1 (* 0.870 d0))))))           ; eq. 9
                    (<= i i-change) (+ dtc (* 0.13 (dst i)))                     ; eq. 12
                    :else (+ dtc (* -2.5 (- (dst i) (dst (dec i))))))]           ; eq. 13
         (max 0.0 next)))
     dtc0
     (range 1 (count dst)))))

;; ----------------------------------------------------------- densities

(defn- day-of-year
  "The day of the year of MJD `mjd`, fractional, 1.0 at January 1 0h."
  [mjd]
  (let [[year] (time/mjd->calendar mjd)]
    (+ 1.0 (- mjd (time/calendar->mjd year 1 1)))))

(defn semiannual
  "Delta log10 rho of the semiannual variation, F(z) G(t) (JB2008 eqs.
  3-7), at height `z` km, day of year `t` and 81-day average indices
  `f10b` `s10b` `m10b`."
  [z t {:keys [f10b s10b m10b]}]
  (let [zz (/ z 1000.0)
        fsmj (- f10b (* 0.70 s10b) (* 0.04 m10b))
        fz (+ 0.269 (* -1.18e-2 fsmj) (* 2.78e-2 zz fsmj) (* -2.78e-2 zz zz fsmj) (* 3.47e-4 zz fsmj fsmj))
        fsm (- f10b (* 0.75 s10b) (* 0.37 m10b))
        w (/ (* 2.0 math/PI (- t 1.0)) 365.0)
        gt (+ -0.363 (* 8.51e-2 (math/sin w)) (* 0.240 (math/cos w)) (* -0.190 (math/sin (* 2 w)))
              (* -0.255 (math/cos (* 2 w)))
              (* fsm (+ -1.79e-2 (* 5.65e-4 (math/sin w)) (* -6.41e-4 (math/cos w))
                        (* -3.42e-3 (math/sin (* 2 w))) (* -1.25e-3 (math/cos (* 2 w))))))]
    (* fz gt)))

(defn seasonal-latitudinal
  "Jacchia 1970 eq. 24: Delta log10 rho of the lower thermosphere at height
  `z`, latitude `lat` (deg) and day of year `t`."
  [z lat t]
  (let [phi (math/to-radians lat) x (- z 90.0)]
    (* 0.02 x (math/signum phi) (math/exp (* -0.045 x)) (math/pow (math/sin phi) 2)
       (math/sin (* (/ (* 2.0 math/PI) 365.2422) (+ (- t 1.0) 100.0))))))

(defn helium-factor
  "Jacchia 1970 eq. 25: n(He)/n0(He) at latitude `lat` and the Sun's
  declination `dec` (deg) 8 days before -- A 0.5, B 2.3, p 2.5, r 4."
  [lat dec]
  (let [eps (math/to-radians 23.44) d (math/to-radians dec) phi (math/to-radians lat)
        q (/ math/PI 4.0)]
    (+ 0.5 (* 1.8 (+ (* (math/pow (/ (- eps d) (* 2.0 eps)) 2.5) (math/pow (math/sin (+ q (* 0.5 phi))) 4))
                     (* (math/pow (/ (+ eps d) (* 2.0 eps)) 2.5) (math/pow (math/sin (- q (* 0.5 phi))) 4)))))))

(defn high-altitude-factor
  "The density factor above 1000 km (JB2006 eqs. 10-11, Table 5): from
  1500 km, 0.22 - 2e-3 F + 1.15e-3 z - 2.11e-6 F z, F the 81-day average
  flux; between 1000 and 1500 km, the cubic in H = (z - 1000)/500 from
  value 1 and slope 0 at 1000 to that line's value and slope at 1500; 1
  below 1000."
  [z f10b]
  (let [line (fn [z] (+ 0.22 (* -2.0e-3 f10b) (* 1.15e-3 z) (* -2.11e-6 f10b z)))]
    (cond
      (<= z 1000.0) 1.0
      (>= z 1500.0) (line z)
      :else (let [f (line 1500.0)
                  slope (* 500.0 (+ 1.15e-3 (* -2.11e-6 f10b)))            ; dF/dH at 1500
                  h (/ (- z 1000.0) 500.0)]
              (+ 1.0 (* (- (* 3.0 (- f 1.0)) slope) h h) (* (- slope (* 2.0 (- f 1.0))) h h h))))))

(defn atmosphere
  "JB2008 at `inputs`: `:mjd`, `:alt` (km, 90 and above), `:lat` (deg),
  `:lst` (local solar time, hours), `:dec` (the Sun's declination, deg;
  by default from its low-accuracy position), the solar indices as
  `night-minimum` takes them, and `:dtc`, the geomagnetic temperature
  rise (`storm-dtc`, `ap-dtc`). Returns number densities, m^-3, `:N2 :O2
  :O :Ar :He :H`; the mass density `:rho`, kg/m^3; the temperature `:t`
  and exospheric temperature `:t-exo`, K."
  [{:keys [mjd alt lat lst dec dtc f10b] :or {dtc 0.0} :as inputs}]
  (let [dec (or dec (math/to-degrees (second (solar/equatorial-low mjd))))
        dec8 (math/to-degrees (second (solar/equatorial-low (- mjd 8.0))))
        day (day-of-year mjd)
        tc (+ (night-minimum inputs) (tc-correction alt lst lat f10b))
        tinf (+ (j/local-temperature tc (math/to-radians lat) (math/to-radians dec) lst {:m 2.5 :r 0.31}) dtc)
        params (update-in j/j70 [:tx 0] + (tx-correction alt lst f10b))
        {:keys [t n]} (j/static params tinf alt)
        n (update n :He * (helium-factor lat dec8))
        f (* (math/pow 10.0 (+ (semiannual alt day inputs) (seasonal-latitudinal alt lat day)))
             (high-altitude-factor alt f10b))
        ;; cm^-3 to m^-3, and so g/m^3 to kg/m^3
        n (update-vals n #(* % f 1e6))]
    (assoc n :rho (* 1e-3 (j/mass-density n)) :t t :t-exo tinf)))
