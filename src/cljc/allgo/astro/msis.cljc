(ns allgo.astro.msis
  "NRLMSISE-00, the Naval Research Laboratory's empirical model of the
  neutral atmosphere from the ground to the exosphere (Picone, Hedin and
  Drob, \"NRLMSISE-00 empirical model of the atmosphere: Statistical
  comparisons and scientific issues\", J. Geophys. Res. 107 (A12), 2002):
  the number densities of He, O, N2, O2, Ar, H, N and anomalous oxygen,
  the total mass density, and the temperature, as functions of day,
  time, place, solar flux and geomagnetic activity.

  Transcribed from the model's Fortran (GTD7, GTD7D, GTS7, GLOBE7, GLOB7S
  and their helpers) as NASA's CCMC ModelWeb archive distributes it, a
  work of the US government; its coefficients are in
  `allgo.astro.msis-data`, generated from the Fortran's BLOCK DATA. The
  Fortran keeps state between calls -- caches, and the magnetic-activity
  terms the last GLOBE7 call leaves for GLOB7S to use -- and this port
  does the same work afresh on every call, carrying that state through
  the calls in the Fortran's order, which is what the Fortran computes
  for any new set of inputs. It works in double precision where the
  Fortran is single, and keeps the Fortran's own constants (its
  degree-to-radian factor is 1.74533e-2, for one).

  Inputs, as `atmosphere` takes them: day of year, UT seconds, altitude
  (km), geodetic latitude and longitude (degrees), local apparent solar
  time (hours), the 81-day average of F10.7 centered on the day and the
  previous day's F10.7, and the geomagnetic index -- the daily Ap, or the
  seven-element history of 3-hour ap the model can use instead."
  (:require [allgo.astro.msis-data :as data]
            [clojure.math :as math]))

;; The Fortran's constants, as it has them
(def ^:private dgtr 1.74533e-2)                 ; degrees to radians
(def ^:private dr 1.72142e-2)                   ; days to radians of the year
(def ^:private hr 0.2618)                       ; hours to radians
(def ^:private sr 7.2722e-5)                    ; seconds to radians
(def ^:private rgas 831.4)

(defn- p
  "Element `k`, 1-based as the Fortran counts, of coefficient array `a`."
  [a k]
  (nth a (dec k)))

(defn- pd [i j] (p (nth data/pd (dec j)) i))
(defn- pdm [i j] (p (nth data/pdm (dec j)) i))
(defn- pdl [i j] (p (nth data/pdl (dec j)) i))
(defn- ptm [i] (p data/ptm i))
(defn- pavgm [i] (p data/pavgm i))
(defn- pd-col [j] (nth data/pd (dec j)))
(defn- ptl-col [j] (nth data/ptl (dec j)))
(defn- pma-col [j] (nth data/pma (dec j)))

;; ------------------------------------------------------------- switches

(defn- switches
  "TSELEC: the main-effect switches SW and cross-term switches SWC of the
  25 settings `sv` -- 0 off, 1 on, 2 main effect off but cross terms on;
  -1 for switch 9 selects the 3-hour ap history."
  [sv]
  {:sw (mapv #(rem % 2.0) sv)
   :swc (mapv #(if (#{1.0 2.0} (abs %)) 1.0 0.0) sv)})

;; ---------------------------------------------------------- gravity

(defn- glatf
  "Gravity at the surface (cm/s^2) and the effective Earth radius (km) at
  latitude `lat`."
  [lat]
  (let [c2 (math/cos (* 2.0 dgtr lat))
        gv (* 980.616 (- 1.0 (* 0.0026373 c2)))]
    [gv (* (/ (* 2.0 gv) (+ 3.085462e-6 (* 2.27e-9 c2))) 1e-5)]))

(defn- scalh [{:keys [gsurf re]} alt xm temp]
  (/ (* rgas temp) (* (/ gsurf (math/pow (+ 1.0 (/ alt re)) 2)) xm)))

;; ------------------------------------------------ horizontal variations

(defn- legendre
  "PLG(n, m), the associated Legendre functions of sin(lat) the model
  uses, as `(get-in plg [n m])` with the Fortran's indices."
  [lat]
  (let [c (math/sin (* lat dgtr)) s (math/cos (* lat dgtr))
        c2 (* c c) c4 (* c2 c2) s2 (* s s)
        m1 (let [p6 (/ (+ (* 63.0 c2 c2 c) (* -70.0 c2 c) (* 15.0 c)) 8.0)
                 p5 (/ (+ (* 35.0 c4) (* -30.0 c2) 3.0) 8.0)]
             {2 c 3 (* 0.5 (- (* 3.0 c2) 1.0)) 4 (* 0.5 (- (* 5.0 c c2) (* 3.0 c)))
              5 p5 6 p6 7 (/ (- (* 11.0 c p6) (* 5.0 p5)) 6.0)})
        m2 (let [p5 (* 2.5 (- (* 7.0 c2 c) (* 3.0 c)) s)
                 p6 (* 1.875 (+ (* 21.0 c4) (* -14.0 c2) 1.0) s)]
             {2 s 3 (* 3.0 c s) 4 (* 1.5 (- (* 5.0 c2) 1.0) s) 5 p5 6 p6
              7 (/ (- (* 11.0 c p6) (* 6.0 p5)) 5.0)})
        m3 (let [p3 (* 3.0 s2) p4 (* 15.0 s2 c) p5 (* 7.5 (- (* 7.0 c2) 1.0) s2)
                 p6 (- (* 3.0 c p5) (* 2.0 p4))
                 p7 (/ (- (* 11.0 c p6) (* 7.0 p5)) 4.0)]
             {3 p3 4 p4 5 p5 6 p6 7 p7 8 (/ (- (* 13.0 c p7) (* 8.0 p6)) 5.0)})
        m4 (let [p4 (* 15.0 s2 s) p5 (* 105.0 s2 s c)
                 p6 (/ (- (* 9.0 c p5) (* 7.0 p4)) 2.0)]
             {4 p4 5 p5 6 p6 7 (/ (- (* 11.0 c p6) (* 8.0 p5)) 3.0)})]
    (fn [n m] (get (case (int m) 1 m1 2 m2 3 m3 4 m4) (int n) 0.0))))

(defn- context
  "What every GLOBE7 and GLOB7S call shares: the switches, the Legendre
  functions, local-time harmonics, day and solar-flux terms."
  [{:keys [doy sec lat lon lst f107a f107 ap]} sv]
  (let [{:keys [sw swc]} (switches sv)]
    {:sw sw :swc swc
     :sw9 (if (neg? (sw 8)) -1.0 1.0)
     :plg (legendre lat)
     :lat lat :lon lon :sec sec :tloc lst
     :stloc (math/sin (* hr lst)) :ctloc (math/cos (* hr lst))
     :s2tloc (math/sin (* 2.0 hr lst)) :c2tloc (math/cos (* 2.0 hr lst))
     :s3tloc (math/sin (* 3.0 hr lst)) :c3tloc (math/cos (* 3.0 hr lst))
     :day (double doy) :f107a f107a :f107 f107
     :df (- f107 f107a) :dfa (- f107a 150.0)
     :ap (if (number? ap) [ap] (vec ap))}))

(defn- sw [ctx i] (nth (:sw ctx) (dec i)))
(defn- swc [ctx i] (nth (:swc ctx) (dec i)))

(defn- globe7
  "GLOBE7: the model's function G(L) of the horizontal, seasonal, solar
  and magnetic variations, for coefficient array `P`, and the state it
  leaves for GLOB7S -- the magnetic activity function of the daily Ap
  (APDF), always renewed in that mode, or of the 3-hour history (APT),
  renewed only by an array with P(52) nonzero. Returns [G state]."
  [{:keys [plg day df dfa tloc ctloc stloc c2tloc s2tloc c3tloc s3tloc sec lat lon sw9 ap] :as ctx} state P]
  (let [P (fn [k] (p P k))
        cd14 (math/cos (* dr (- day (P 14))))
        cd18 (math/cos (* 2.0 dr (- day (P 18))))
        cd32 (math/cos (* dr (- day (P 32))))
        cd39 (math/cos (* 2.0 dr (- day (P 39))))
        sw (partial sw ctx) swc (partial swc ctx)
        f1 (+ 1.0 (* (+ (* (P 48) dfa) (* (P 20) df) (* (P 21) df df)) (swc 1)))
        f2 (+ 1.0 (* (+ (* (P 50) dfa) (* (P 20) df) (* (P 21) df df)) (swc 1)))
        t1 (+ (* (P 20) df (+ 1.0 (* (P 60) dfa))) (* (P 21) df df) (* (P 22) dfa) (* (P 30) dfa dfa))
        t2 (+ (* (P 2) (plg 3 1)) (* (P 3) (plg 5 1)) (* (P 23) (plg 7 1))
              (* (P 15) (plg 3 1) dfa (swc 1))
              (* (P 27) (plg 2 1)))
        t3 (* (P 19) cd32)
        t4 (* (+ (P 16) (* (P 17) (plg 3 1))) cd18)
        t5 (* f1 (+ (* (P 10) (plg 2 1)) (* (P 11) (plg 4 1))) cd14)
        t6 (* (P 38) (plg 2 1) cd39)
        t7 (if (zero? (sw 7)) 0.0
               (let [t71 (* (P 12) (plg 3 2) cd14 (swc 5))
                     t72 (* (P 13) (plg 3 2) cd14 (swc 5))]
                 (* f2 (+ (* (+ (* (P 4) (plg 2 2)) (* (P 5) (plg 4 2)) (* (P 28) (plg 6 2)) t71) ctloc)
                          (* (+ (* (P 7) (plg 2 2)) (* (P 8) (plg 4 2)) (* (P 29) (plg 6 2)) t72) stloc)))))
        t8 (if (zero? (sw 8)) 0.0
               (let [t81 (* (+ (* (P 24) (plg 4 3)) (* (P 36) (plg 6 3))) cd14 (swc 5))
                     t82 (* (+ (* (P 34) (plg 4 3)) (* (P 37) (plg 6 3))) cd14 (swc 5))]
                 (* f2 (+ (* (+ (* (P 6) (plg 3 3)) (* (P 42) (plg 5 3)) t81) c2tloc)
                          (* (+ (* (P 9) (plg 3 3)) (* (P 43) (plg 5 3)) t82) s2tloc)))))
        t14 (if (zero? (sw 14)) 0.0
                (* f2 (+ (* (+ (* (P 40) (plg 4 4)) (* (+ (* (P 94) (plg 5 4)) (* (P 47) (plg 7 4))) cd14 (swc 5)))
                            s3tloc)
                         (* (+ (* (P 41) (plg 4 4)) (* (+ (* (P 95) (plg 5 4)) (* (P 49) (plg 7 4))) cd14 (swc 5)))
                            c3tloc))))
        ;; magnetic activity: the daily Ap, or the 3-hour history
        [t9 state]
        (if (= sw9 -1.0)
          (if (zero? (P 52))
            [0.0 state]
            (let [p25 (let [v (P 25)] (if (< v 1e-4) 1e-4 v))
                  g0 (fn [a] (+ (- a 4.0)
                                (* (- (P 26) 1.0)
                                   (+ (- a 4.0) (/ (- (math/exp (* (- (abs p25)) (- a 4.0))) 1.0) (abs p25))))))
                  ex (min 0.99999 (math/exp (/ (* -10800.0 (abs (P 52))) (+ 1.0 (* (P 139) (- 45.0 (abs lat)))))))
                  sumex (+ 1.0 (* (/ (- 1.0 (math/pow ex 19)) (- 1.0 ex)) (math/sqrt ex)))
                  apt (/ (+ (g0 (ap 1))
                            (+ (* (g0 (ap 2)) ex) (* (g0 (ap 3)) ex ex) (* (g0 (ap 4)) (math/pow ex 3))
                               (* (+ (* (g0 (ap 5)) (math/pow ex 4)) (* (g0 (ap 6)) (math/pow ex 12)))
                                  (/ (- 1.0 (math/pow ex 8)) (- 1.0 ex)))))
                         sumex)
                  state (assoc state :apt apt)]
              [(if (zero? (sw 9)) 0.0
                   (* apt (+ (P 51) (* (P 97) (plg 3 1)) (* (P 55) (plg 5 1))
                             (* (+ (* (P 126) (plg 2 1)) (* (P 127) (plg 4 1)) (* (P 128) (plg 6 1))) cd14 (swc 5))
                             (* (+ (* (P 129) (plg 2 2)) (* (P 130) (plg 4 2)) (* (P 131) (plg 6 2)))
                                (swc 7) (math/cos (* hr (- tloc (P 132))))))))
               state]))
          (let [apd (- (ap 0) 4.0)
                p44 (let [v (P 44)] (if (neg? v) 1e-5 v))
                apdf (+ apd (* (- (P 45) 1.0) (+ apd (/ (- (math/exp (* (- p44) apd)) 1.0) p44))))]
            [(if (zero? (sw 9)) 0.0
                 (* apdf (+ (P 33) (* (P 46) (plg 3 1)) (* (P 35) (plg 5 1))
                            (* (+ (* (P 101) (plg 2 1)) (* (P 102) (plg 4 1)) (* (P 103) (plg 6 1))) cd14 (swc 5))
                            (* (+ (* (P 122) (plg 2 2)) (* (P 123) (plg 4 2)) (* (P 124) (plg 6 2)))
                               (swc 7) (math/cos (* hr (- tloc (P 125))))))))
             (assoc state :apdf apdf)]))
        {:keys [apdf apt]} state
        longitude? (not (or (zero? (sw 10)) (<= lon -1000.0)))
        t11 (if (or (not longitude?) (zero? (sw 11))) 0.0
                (* (+ 1.0 (* (P 81) dfa (swc 1)))
                   (+ (* (+ (* (P 65) (plg 3 2)) (* (P 66) (plg 5 2)) (* (P 67) (plg 7 2))
                            (* (P 104) (plg 2 2)) (* (P 105) (plg 4 2)) (* (P 106) (plg 6 2))
                            (* (swc 5) (+ (* (P 110) (plg 2 2)) (* (P 111) (plg 4 2)) (* (P 112) (plg 6 2))) cd14))
                         (math/cos (* dgtr lon)))
                      (* (+ (* (P 91) (plg 3 2)) (* (P 92) (plg 5 2)) (* (P 93) (plg 7 2))
                            (* (P 107) (plg 2 2)) (* (P 108) (plg 4 2)) (* (P 109) (plg 6 2))
                            (* (swc 5) (+ (* (P 113) (plg 2 2)) (* (P 114) (plg 4 2)) (* (P 115) (plg 6 2))) cd14))
                         (math/sin (* dgtr lon))))))
        t12 (if (or (not longitude?) (zero? (sw 12))) 0.0
                (+ (* (+ 1.0 (* (P 96) (plg 2 1))) (+ 1.0 (* (P 82) dfa (swc 1)))
                      (+ 1.0 (* (P 120) (plg 2 1) (swc 5) cd14))
                      (+ (* (P 69) (plg 2 1)) (* (P 70) (plg 4 1)) (* (P 71) (plg 6 1)))
                      (math/cos (* sr (- sec (P 72)))))
                   (* (swc 11) (+ (* (P 77) (plg 4 3)) (* (P 78) (plg 6 3)) (* (P 79) (plg 8 3)))
                      (math/cos (+ (* sr (- sec (P 80))) (* 2.0 dgtr lon)))
                      (+ 1.0 (* (P 138) dfa (swc 1))))))
        t13 (cond
              (or (not longitude?) (zero? (sw 13))) 0.0
              (not= sw9 -1.0)
              (+ (* apdf (swc 11) (+ 1.0 (* (P 121) (plg 2 1)))
                    (+ (* (P 61) (plg 3 2)) (* (P 62) (plg 5 2)) (* (P 63) (plg 7 2)))
                    (math/cos (* dgtr (- lon (P 64)))))
                 (* apdf (swc 11) (swc 5)
                    (+ (* (P 116) (plg 2 2)) (* (P 117) (plg 4 2)) (* (P 118) (plg 6 2)))
                    cd14 (math/cos (* dgtr (- lon (P 119)))))
                 (* apdf (swc 12)
                    (+ (* (P 84) (plg 2 1)) (* (P 85) (plg 4 1)) (* (P 86) (plg 6 1)))
                    (math/cos (* sr (- sec (P 76))))))
              (zero? (P 52)) 0.0
              :else
              (+ (* apt (swc 11) (+ 1.0 (* (P 133) (plg 2 1)))
                    (+ (* (P 53) (plg 3 2)) (* (P 99) (plg 5 2)) (* (P 68) (plg 7 2)))
                    (math/cos (* dgtr (- lon (P 98)))))
                 (* apt (swc 11) (swc 5)
                    (+ (* (P 134) (plg 2 2)) (* (P 135) (plg 4 2)) (* (P 136) (plg 6 2)))
                    cd14 (math/cos (* dgtr (- lon (P 137)))))
                 (* apt (swc 12)
                    (+ (* (P 56) (plg 2 1)) (* (P 57) (plg 4 1)) (* (P 58) (plg 6 1)))
                    (math/cos (* sr (- sec (P 59)))))))
        t [t1 t2 t3 t4 t5 t6 t7 t8 t9 0.0 t11 t12 t13 t14]]
    [(reduce + (P 31) (map (fn [i ti] (* (abs (sw i)) ti)) (range 1 15) t)) state]))

(defn- glob7s
  "GLOB7S: the version of G(L) for the lower atmosphere, for coefficient
  array `P`, using the magnetic activity function `state` holds."
  [{:keys [plg day dfa ctloc stloc c2tloc s2tloc c3tloc s3tloc lon] :as ctx} {:keys [apdf apt]} P]
  (let [P (fn [k] (p P k))
        sw (partial sw ctx) swc (partial swc ctx)
        cd32 (math/cos (* dr (- day (P 32))))
        cd18 (math/cos (* 2.0 dr (- day (P 18))))
        cd14 (math/cos (* dr (- day (P 14))))
        cd39 (math/cos (* 2.0 dr (- day (P 39))))
        t1 (* (P 22) dfa)
        t2 (+ (* (P 2) (plg 3 1)) (* (P 3) (plg 5 1)) (* (P 23) (plg 7 1))
              (* (P 27) (plg 2 1)) (* (P 15) (plg 4 1)) (* (P 60) (plg 6 1)))
        t3 (* (+ (P 19) (* (P 48) (plg 3 1)) (* (P 30) (plg 5 1))) cd32)
        t4 (* (+ (P 16) (* (P 17) (plg 3 1)) (* (P 31) (plg 5 1))) cd18)
        t5 (* (+ (* (P 10) (plg 2 1)) (* (P 11) (plg 4 1)) (* (P 21) (plg 6 1))) cd14)
        t6 (* (P 38) (plg 2 1) cd39)
        t7 (if (zero? (sw 7)) 0.0
               (let [t71 (* (P 12) (plg 3 2) cd14 (swc 5))
                     t72 (* (P 13) (plg 3 2) cd14 (swc 5))]
                 (+ (* (+ (* (P 4) (plg 2 2)) (* (P 5) (plg 4 2)) t71) ctloc)
                    (* (+ (* (P 7) (plg 2 2)) (* (P 8) (plg 4 2)) t72) stloc))))
        t8 (if (zero? (sw 8)) 0.0
               (let [t81 (* (+ (* (P 24) (plg 4 3)) (* (P 36) (plg 6 3))) cd14 (swc 5))
                     t82 (* (+ (* (P 34) (plg 4 3)) (* (P 37) (plg 6 3))) cd14 (swc 5))]
                 (+ (* (+ (* (P 6) (plg 3 3)) (* (P 42) (plg 5 3)) t81) c2tloc)
                    (* (+ (* (P 9) (plg 3 3)) (* (P 43) (plg 5 3)) t82) s2tloc))))
        t14 (if (zero? (sw 14)) 0.0
                (+ (* (P 40) (plg 4 4) s3tloc) (* (P 41) (plg 4 4) c3tloc)))
        t9 (cond
             (== (sw 9) 1.0) (* apdf (+ (P 33) (* (P 46) (plg 3 1) (swc 2))))
             (== (sw 9) -1.0) (+ (* (P 51) apt) (* (P 97) (plg 3 1) apt (swc 2)))
             :else 0.0)
        t11 (if (or (zero? (sw 10)) (zero? (sw 11)) (<= lon -1000.0)) 0.0
                (* (+ 1.0
                      (* (plg 2 1) (+ (* (P 81) (swc 5) (math/cos (* dr (- day (P 82)))))
                                      (* (P 86) (swc 6) (math/cos (* 2.0 dr (- day (P 87)))))))
                      (* (P 84) (swc 3) (math/cos (* dr (- day (P 85)))))
                      (* (P 88) (swc 4) (math/cos (* 2.0 dr (- day (P 89))))))
                   (+ (* (+ (* (P 65) (plg 3 2)) (* (P 66) (plg 5 2)) (* (P 67) (plg 7 2))
                            (* (P 75) (plg 2 2)) (* (P 76) (plg 4 2)) (* (P 77) (plg 6 2)))
                         (math/cos (* dgtr lon)))
                      (* (+ (* (P 91) (plg 3 2)) (* (P 92) (plg 5 2)) (* (P 93) (plg 7 2))
                            (* (P 78) (plg 2 2)) (* (P 79) (plg 4 2)) (* (P 80) (plg 6 2)))
                         (math/sin (* dgtr lon))))))
        t [t1 t2 t3 t4 t5 t6 t7 t8 t9 0.0 t11 0.0 0.0 t14]]
    (reduce + 0.0 (map (fn [i ti] (* (abs (sw i)) ti)) (range 1 15) t))))

;; ------------------------------------------------------------ splines

(defn- spline
  "SPLINE: the second derivatives of the cubic spline through `x` `y`
  with end slopes `yp1` and `ypn` (Numerical Recipes' algorithm, as the
  Fortran adapts it)."
  [x y yp1 ypn]
  (let [n (count x)
        [y2 u] (loop [i 1
                      y2 [(if (> yp1 0.99e30) 0.0 -0.5)]
                      u [(if (> yp1 0.99e30) 0.0
                             (* (/ 3.0 (- (x 1) (x 0))) (- (/ (- (y 1) (y 0)) (- (x 1) (x 0))) yp1)))]]
                 (if (>= i (dec n))
                   [y2 u]
                   (let [sig (/ (- (x i) (x (dec i))) (- (x (inc i)) (x (dec i))))
                         pp (+ (* sig (y2 (dec i))) 2.0)]
                     (recur (inc i)
                            (conj y2 (/ (- sig 1.0) pp))
                            (conj u (/ (- (/ (* 6.0 (- (/ (- (y (inc i)) (y i)) (- (x (inc i)) (x i)))
                                                       (/ (- (y i) (y (dec i))) (- (x i) (x (dec i))))))
                                             (- (x (inc i)) (x (dec i))))
                                          (* sig (u (dec i))))
                                       pp))))))
        [qn un] (if (> ypn 0.99e30) [0.0 0.0]
                    [0.5 (* (/ 3.0 (- (x (dec n)) (x (- n 2))))
                            (- ypn (/ (- (y (dec n)) (y (- n 2))) (- (x (dec n)) (x (- n 2))))))])
        y2 (conj y2 (/ (- un (* qn (u (- n 2)))) (+ (* qn (y2 (- n 2))) 1.0)))]
    (reduce (fn [y2 k] (assoc y2 k (+ (* (y2 k) (y2 (inc k))) (u k)))) y2 (range (- n 2) -1 -1))))

(defn- splint
  "SPLINT: the cubic spline's value at `x`."
  [xa ya y2a x]
  (let [[klo khi] (loop [klo 0 khi (dec (count xa))]
                    (if (> (- khi klo) 1)
                      (let [k (quot (+ khi klo) 2)]
                        (if (> (xa k) x) (recur klo k) (recur k khi)))
                      [klo khi]))
        h (- (xa khi) (xa klo))
        a (/ (- (xa khi) x) h)
        b (/ (- x (xa klo)) h)]
    (+ (* a (ya klo)) (* b (ya khi))
       (/ (* (+ (* (- (* a a a) a) (y2a klo)) (* (- (* b b b) b) (y2a khi))) h h) 6.0))))

(defn- splini
  "SPLINI: the cubic spline integrated from its first node to `x`."
  [xa ya y2a x]
  (let [n (count xa)]
    (loop [yi 0.0 klo 0 khi 1]
      (if (and (> x (xa klo)) (< khi n))
        (let [xx (if (< khi (dec n)) (min x (xa khi)) x)
              h (- (xa khi) (xa klo))
              a (/ (- (xa khi) xx) h)
              b (/ (- xx (xa klo)) h)
              a2 (* a a) b2 (* b b)]
          (recur (+ yi (* (+ (/ (* (- 1.0 a2) (ya klo)) 2.0) (/ (* b2 (ya khi)) 2.0)
                             (/ (* (+ (* (+ (/ (- (+ 1.0 (* a2 a2))) 4.0) (/ a2 2.0)) (y2a klo))
                                      (* (- (/ (* b2 b2) 4.0) (/ b2 2.0)) (y2a khi)))
                                   h h)
                                6.0))
                          h))
                 (inc klo) (inc khi)))
        yi))))

;; ------------------------------------------------ temperature and density

(defn- zeta [{:keys [re]} zz zl] (/ (* (- zz zl) (+ re zl)) (+ re zz)))

(defn- densu
  "DENSU: [density temperature] at `alt` for a species of molecular weight
  `xm` (0 for the temperature alone) diffusing from density `dlb` at the
  lower boundary `zlb`: Bates's temperature profile above the first node,
  a spline in 1/T through the nodes `tn1` below it, with the first node's
  temperature and gradient from the Bates profile."
  [{:keys [gsurf re] :as g} alt dlb tinf tlb xm alpha zlb s2 zn1 tn1 tgn1]
  (let [za (zn1 0)
        z (max alt za)
        zg2 (zeta g z zlb)
        tt (- tinf (* (- tinf tlb) (math/exp (- (* s2 zg2)))))
        below? (< alt za)
        spline-part
        (when below?
          (let [ta tt
                dta (* (- tinf ta) s2 (math/pow (/ (+ re zlb) (+ re za)) 2))
                tn1 (assoc tn1 0 ta)
                tgn1 (assoc tgn1 0 dta)
                mn (count zn1)
                z (max alt (zn1 (dec mn)))
                z1 (zn1 0) z2 (zn1 (dec mn))
                t1 (tn1 0) t2 (tn1 (dec mn))
                zg (zeta g z z1)
                zgdif (zeta g z2 z1)
                xs (mapv #(/ (zeta g % z1) zgdif) zn1)
                ys (mapv #(/ 1.0 %) tn1)
                yd1 (* (/ (- (tgn1 0)) (* t1 t1)) zgdif)
                yd2 (* (/ (- (tgn1 1)) (* t2 t2)) zgdif (math/pow (/ (+ re z2) (+ re z1)) 2))
                y2out (spline xs ys yd1 yd2)
                x (/ zg zgdif)]
            {:tz (/ 1.0 (splint xs ys y2out x)) :xs xs :ys ys :y2out y2out :x x :z1 z1 :zgdif zgdif :t1 t1}))
        tz (if below? (:tz spline-part) tt)]
    (if (zero? xm)
      [tz tz]
      (let [glb (/ gsurf (math/pow (+ 1.0 (/ zlb re)) 2))
            gamma (/ (* xm glb) (* s2 rgas tinf))
            expl (let [e (math/exp (- (* s2 gamma zg2)))] (if (or (> e 50.0) (<= tt 0.0)) 50.0 e))
            densa (* dlb (math/pow (/ tlb tt) (+ 1.0 alpha gamma)) expl)]
        (if-not below?
          [densa tz]
          (let [{:keys [xs ys y2out x z1 zgdif t1]} spline-part
                glb (/ gsurf (math/pow (+ 1.0 (/ z1 re)) 2))
                gamm (/ (* xm glb zgdif) rgas)
                expl (let [e (* gamm (splini xs ys y2out x))] (if (or (> e 50.0) (<= tz 0.0)) 50.0 e))]
            [(* densa (math/pow (/ t1 tz) (+ 1.0 alpha)) (math/exp (- expl))) tz]))))))

(defn- densm-segment
  "One of DENSM's two spline segments, from its first node down to `z`:
  [density factor, temperature]."
  [{:keys [gsurf re] :as g} z xm zn tn tgn]
  (let [mn (count zn)
        z1 (zn 0) z2 (zn (dec mn))
        t1 (tn 0) t2 (tn (dec mn))
        zg (zeta g z z1)
        zgdif (zeta g z2 z1)
        xs (mapv #(/ (zeta g % z1) zgdif) zn)
        ys (mapv #(/ 1.0 %) tn)
        yd1 (* (/ (- (tgn 0)) (* t1 t1)) zgdif)
        yd2 (* (/ (- (tgn 1)) (* t2 t2)) zgdif (math/pow (/ (+ re z2) (+ re z1)) 2))
        y2out (spline xs ys yd1 yd2)
        x (/ zg zgdif)
        tz (/ 1.0 (splint xs ys y2out x))]
    (if (zero? xm)
      [1.0 tz]
      (let [glb (/ gsurf (math/pow (+ 1.0 (/ z1 re)) 2))
            gamm (/ (* xm glb zgdif) rgas)
            expl (min 50.0 (* gamm (splini xs ys y2out x)))]
        [(* (/ t1 tz) (math/exp (- expl))) tz]))))

(defn- densm
  "DENSM: [density temperature] at `alt` below the thermosphere, from
  density `d0` at the upper boundary: the stratosphere and mesosphere
  spline through nodes `zn2`, then the troposphere and lower
  stratosphere through `zn3`. With `xm` 0, the temperature alone."
  [g alt d0 xm zn3 tn3 tgn3 zn2 tn2 tgn2]
  (if (> alt (zn2 0))
    [(if (zero? xm) nil d0) nil]
    (let [[f2 tz2] (densm-segment g (max alt (peek zn2)) xm zn2 tn2 tgn2)
          [f3 tz] (if (> alt (zn3 0)) [1.0 tz2] (densm-segment g alt xm zn3 tn3 tgn3))]
      [(if (zero? xm) tz (* d0 f2 f3)) tz])))

(defn- dnet
  "DNET: the diffusive density `dd` and the fully mixed `dm` joined across
  the turbopause, over the transition scale length `zhm`."
  [dd dm zhm xmm xm]
  (let [a (/ zhm (- xmm xm))]
    (cond
      (and (zero? dm) (zero? dd)) 1.0
      (zero? dm) dd
      (zero? dd) dm
      :else (let [ylog (* a (math/log (/ dm dd)))]
              (cond (< ylog -10.0) dd
                    (> ylog 10.0) dm
                    :else (* dd (math/pow (+ 1.0 (math/exp ylog)) (/ 1.0 a))))))))

(defn- ccor
  "CCOR: the chemistry and dissociation correction, ratio `r` at the low
  end, `h1` the transition scale length, `zh` its middle."
  [alt r h1 zh]
  (let [e (/ (- alt zh) h1)]
    (math/exp (cond (> e 70.0) 0.0 (< e -70.0) r :else (/ r (+ 1.0 (math/exp e)))))))

(defn- ccor2
  "CCOR2: the O and O2 correction, with two scale lengths."
  [alt r h1 zh h2]
  (let [e1 (/ (- alt zh) h1) e2 (/ (- alt zh) h2)]
    (math/exp (cond (or (> e1 70.0) (> e2 70.0)) 0.0
                    (and (< e1 -70.0) (< e2 -70.0)) r
                    :else (/ r (+ 1.0 (* 0.5 (+ (math/exp e1) (math/exp e2)))))))))

;; ---------------------------------------------------------- thermosphere

(def ^:private alpha [-0.38 0.0 0.0 0.0 0.17 0.0 -0.38 0.0 0.0])
(def ^:private altl [200.0 300.0 160.0 250.0 240.0 450.0 320.0 450.0])

(defn- gts7
  "GTS7, the thermosphere above 72.5 km, for all species (`mass` 48) or
  N2 alone (28): number densities (cm^-3) and mass density (g/cm^3), the
  exospheric temperature and the temperature at `alt`, the lower
  thermosphere nodes, the fully mixed N2 density, and the GLOBE7 state
  the last call left."
  [{:keys [f107a day] :as ctx} g alt mass]
  (let [state0 {:apdf 0.0 :apt 0.0}
        za (pdl 16 2)
        zn1 [za 110.0 100.0 90.0 72.5]
        sw (partial sw ctx)
        pt data/pt ps data/ps
        [tinf state] (if (> alt za)
                       (let [[gv st] (globe7 ctx state0 pt)] [(* (ptm 1) (p pt 1) (+ 1.0 (* (sw 16) gv))) st])
                       [(* (ptm 1) (p pt 1)) state0])
        [g0 state] (if (> alt 72.5)
                     (let [[gv st] (globe7 ctx state ps)] [(* (ptm 4) (p ps 1) (+ 1.0 (* (sw 19) gv))) st])
                     [(* (ptm 4) (p ps 1)) state])
        [tlb state] (let [[gv st] (globe7 ctx state (pd-col 4))]
                      [(* (ptm 2) (+ 1.0 (* (sw 17) gv)) (pd 1 4)) st])
        s (/ g0 (- tinf tlb))
        [tn1 tgn1]
        (let [ptl1 (p (ptl-col 4) 1)]
          (if (< alt 300.0)
            (let [gs (fn [col] (glob7s ctx state col))
                  tn5 (/ (* (ptm 5) ptl1) (- 1.0 (* (sw 18) (sw 20) (gs (ptl-col 4)))))]
              [[0.0
                (/ (* (ptm 7) (p (ptl-col 1) 1)) (- 1.0 (* (sw 18) (gs (ptl-col 1)))))
                (/ (* (ptm 3) (p (ptl-col 2) 1)) (- 1.0 (* (sw 18) (gs (ptl-col 2)))))
                (/ (* (ptm 8) (p (ptl-col 3) 1)) (- 1.0 (* (sw 18) (gs (ptl-col 3)))))
                tn5]
               [0.0 (/ (* (ptm 9) (p (pma-col 9) 1) (+ 1.0 (* (sw 18) (sw 20) (gs (pma-col 9)))) tn5 tn5)
                       (math/pow (* (ptm 5) ptl1) 2))]])
            (let [tn5 (* (ptm 5) ptl1)]
              [[0.0 (* (ptm 7) (p (ptl-col 1) 1)) (* (ptm 3) (p (ptl-col 2) 1)) (* (ptm 8) (p (ptl-col 3) 1)) tn5]
               [0.0 (/ (* (ptm 9) (p (pma-col 9) 1) tn5 tn5) (math/pow (* (ptm 5) ptl1) 2))]])))
        zlb (ptm 6)
        du (fn [z dlb xm alpha] (densu g z dlb tinf tlb xm alpha zlb s zn1 tn1 tgn1))
        z alt
        t-at-z (second (du z 1.0 0.0 0.0))
        xmm (pdm 5 3)
        diffusive? (not (zero? (sw 15)))
        ;; N2, which fixes the turbopause and the mixing ratios
        [g28 state] (let [[gv st] (globe7 ctx state (pd-col 3))] [(* (sw 21) gv) st])
        zhf (* (pdl 25 2) (+ 1.0 (* (sw 5) (pdl 25 1) (math/sin (* dgtr (:lat ctx)))
                                    (math/cos (* dr (- day (p pt 14)))))))
        db28 (* (pdm 1 3) (math/exp g28) (pd 1 3))
        zhm28 (* (pdm 4 3) (pdl 6 2))
        b28 (first (du (* (pdm 3 3) zhf) db28 (- 28.0 xmm) (- (alpha 2) 1.0)))
        dm28 (when (and (<= z (altl 2)) diffusive?) (first (du z b28 xmm (alpha 2))))
        d3 (let [d (first (du z db28 28.0 (alpha 2)))] (if dm28 (dnet d dm28 zhm28 xmm 28.0) d))
        ;; a species: its variation, diffusive density and, below its
        ;; altitude limit, the mixing across the turbopause
        species (fn [state col-j m-j xm alpha-i altl-i]
                  (let [[gv st] (globe7 ctx state (pd-col col-j))
                        dlb (* (pdm 1 m-j) (math/exp (* (sw 21) gv)) (pd 1 col-j))
                        d (first (du z dlb xm alpha-i))]
                    {:state st :dlb dlb :d d
                     :mixed (when (and (<= z altl-i) diffusive?)
                              (let [b (first (du (pdm 3 m-j) dlb (- xm xmm) (- alpha-i 1.0)))
                                    dm (first (du z b xmm 0.0))]
                                {:b b :d (dnet d dm zhm28 xmm xm)}))}))]
    (if (= mass 28)
      {:d [0.0 0.0 d3 0.0 0.0 0.0 0.0 0.0 0.0] :tinf tinf :t t-at-z
       :tn1 tn1 :tgn1 tgn1 :dm28 dm28 :state state}
      (let [;; He
            he (species state 1 1 4.0 (alpha 0) (altl 0))
            d1 (if-let [{:keys [b d]} (:mixed he)]
                 (* d (ccor z (math/log (/ (* b28 (pdm 2 1)) b)) (* (pdm 6 1) (pdl 2 2)) (* (pdm 5 1) (pdl 1 2))))
                 (:d he))
            ;; O
            o (species (:state he) 2 2 16.0 (alpha 1) (altl 1))
            d2 (if-let [{:keys [d]} (:mixed o)]
                 (let [rl (* (pdm 2 2) (pdl 17 2) (+ 1.0 (* (sw 1) (pdl 24 1) (- f107a 150.0))))]
                   (* d
                      (ccor2 z rl (* (pdm 6 2) (pdl 4 2)) (* (pdm 5 2) (pdl 3 2)) (* (pdm 6 2) (pdl 5 2)))
                      (ccor z (* (pdm 4 2) (pdl 15 2)) (* (pdm 8 2) (pdl 14 2)) (* (pdm 7 2) (pdl 13 2)))))
                 (:d o))
            ;; O2, its departure from diffusive equilibrium corrected at every height
            o2 (species (:state o) 5 4 32.0 (alpha 3) (altl 3))
            d4 (if-not diffusive?
                 (:d o2)
                 (let [d (if-let [{:keys [b d]} (:mixed o2)]
                           (* d (ccor z (math/log (/ (* b28 (pdm 2 4)) b)) (* (pdm 6 4) (pdl 8 2)) (* (pdm 5 4) (pdl 7 2))))
                           (:d o2))
                       rc32 (* (pdm 4 4) (pdl 24 2) (+ 1.0 (* (sw 1) (pdl 24 1) (- f107a 150.0))))]
                   (* d (ccor2 z rc32 (* (pdm 8 4) (pdl 23 2)) (* (pdm 7 4) (pdl 22 2)) (* (pdm 8 4) (pdl 23 1))))))
            ;; Ar
            ar (species (:state o2) 6 5 40.0 (alpha 4) (altl 4))
            d5 (if-let [{:keys [b d]} (:mixed ar)]
                 (* d (ccor z (math/log (/ (* b28 (pdm 2 5)) b)) (* (pdm 6 5) (pdl 10 2)) (* (pdm 5 5) (pdl 9 2))))
                 (:d ar))
            ;; H
            h (species (:state ar) 7 6 1.0 (alpha 6) (altl 6))
            d7 (if-let [{:keys [b d]} (:mixed h)]
                 (* d
                    (ccor z (math/log (/ (* b28 (pdm 2 6) (abs (pdl 18 2))) b)) (* (pdm 6 6) (pdl 12 2)) (* (pdm 5 6) (pdl 11 2)))
                    (ccor z (* (pdm 4 6) (pdl 21 2)) (* (pdm 8 6) (pdl 20 2)) (* (pdm 7 6) (pdl 19 2))))
                 (:d h))
            ;; N
            n (species (:state h) 8 7 14.0 (alpha 7) (altl 7))
            d8 (if-let [{:keys [b d]} (:mixed n)]
                 (* d
                    (ccor z (math/log (/ (* b28 (pdm 2 7) (abs (pdl 3 1))) b)) (* (pdm 6 7) (pdl 2 1)) (* (pdm 5 7) (pdl 1 1)))
                    (ccor z (* (pdm 4 7) (pdl 6 1)) (* (pdm 8 7) (pdl 5 1)) (* (pdm 7 7) (pdl 4 1))))
                 (:d n))
            ;; anomalous oxygen, isothermal at its own temperature
            [g16h state] (globe7 ctx (:state n) (pd-col 9))
            db16h (* (pdm 1 8) (math/exp (* (sw 21) g16h)) (pd 1 9))
            tho (* (pdm 10 8) (pdl 7 1))
            dd (first (densu g z db16h tho tho 16.0 (alpha 8) zlb s zn1 tn1 tgn1))
            zsht (pdm 6 8) zmho (pdm 5 8)
            zsho (scalh g zmho 16.0 tho)
            d9 (* dd (math/exp (* (- (/ zsht zsho)) (- (math/exp (- (/ (- z zmho) zsht))) 1.0))))
            d6 (* 1.66e-24 (+ (* 4.0 d1) (* 16.0 d2) (* 28.0 d3) (* 32.0 d4) (* 40.0 d5) d7 (* 14.0 d8)))]
        {:d [d1 d2 d3 d4 d5 d6 d7 d8 d9] :tinf tinf :t t-at-z
         :tn1 tn1 :tgn1 tgn1 :dm28 dm28 :state state}))))

;; ---------------------------------------------------------------- GTD7

(def ^:private zn3 [32.5 20.0 15.0 10.0 0.0])
(def ^:private zn2 [72.5 55.0 45.0 32.5])
(def ^:private zmix 62.5)

(defn- gtd7
  "GTD7: number densities (cm^-3), mass density without anomalous oxygen
  (g/cm^3) and temperatures, from the ground up."
  [{:keys [alt lat] :as inputs} sv]
  (let [ctx (context inputs sv)
        sw (partial sw ctx)
        [gsurf re] (glatf (if (zero? (sw 2)) 45.0 lat))
        g {:gsurf gsurf :re re}
        xmm (pdm 5 3)
        mss (if (< alt zmix) 28 48)
        {:keys [d tinf t tn1 tgn1 dm28 state]} (gts7 ctx g (max alt (zn2 0)) mss)]
    (if (>= alt (zn2 0))
      {:d d :tinf tinf :t t}
      (let [gs (fn [j] (glob7s ctx state (pma-col j)))
            pma1 (fn [j] (p (pma-col j) 1))
            tn2-4 (/ (* (pma1 3) (pavgm 3)) (- 1.0 (* (sw 20) (sw 22) (gs 3))))
            tn2 [(tn1 4)
                 (/ (* (pma1 1) (pavgm 1)) (- 1.0 (* (sw 20) (gs 1))))
                 (/ (* (pma1 2) (pavgm 2)) (- 1.0 (* (sw 20) (gs 2))))
                 tn2-4]
            tgn2 [(tgn1 1)
                  (/ (* (pavgm 9) (pma1 10) (+ 1.0 (* (sw 20) (sw 22) (gs 10))) tn2-4 tn2-4)
                     (math/pow (* (pma1 3) (pavgm 3)) 2))]
            [tn3 tgn3] (if (> alt (zn3 0))
                         [[tn2-4 0.0 0.0 0.0 0.0] [0.0 0.0]]
                         (let [tn3-5 (/ (* (pma1 7) (pavgm 7)) (- 1.0 (* (sw 22) (gs 7))))]
                           [[tn2-4
                             (/ (* (pma1 4) (pavgm 4)) (- 1.0 (* (sw 22) (gs 4))))
                             (/ (* (pma1 5) (pavgm 5)) (- 1.0 (* (sw 22) (gs 5))))
                             (/ (* (pma1 6) (pavgm 6)) (- 1.0 (* (sw 22) (gs 6))))
                             tn3-5]
                            [(tgn2 1)
                             (/ (* (pma1 8) (pavgm 8) (+ 1.0 (* (sw 22) (gs 8))) tn3-5 tn3-5)
                                (math/pow (* (pma1 7) (pavgm 7)) 2))]]))
            ;; a linear transition to full mixing below zn2(1)
            dmc (if (> alt zmix) (- 1.0 (/ (- (zn2 0) alt) (- (zn2 0) zmix))) 0.0)
            dz28 (d 2)
            [n2 tz] (densm g alt dm28 xmm zn3 tn3 tgn3 zn2 tn2 tgn2)
            d3 (* n2 (+ 1.0 (* (- (/ (d 2) dm28) 1.0) dmc)))
            mixed (fn [i m] (* d3 (pdm 2 m) (+ 1.0 (* (- (/ (d i) (* dz28 (pdm 2 m))) 1.0) dmc))))
            d1 (mixed 0 1) d4 (mixed 3 4) d5 (mixed 4 5)
            d6 (* 1.66e-24 (+ (* 4.0 d1) (* 28.0 d3) (* 32.0 d4) (* 40.0 d5)))]
        {:d [d1 0.0 d3 d4 d5 d6 0.0 0.0 0.0] :tinf tinf :t tz}))))

;; ------------------------------------------------------------ interface

(def default-switches
  "All 25 of the model's variations on, the daily Ap used."
  (vec (repeat 25 1.0)))

(defn atmosphere
  "NRLMSISE-00 at `inputs`, a map of

    :doy    day of year, 1 to 366
    :sec    UT, seconds of the day
    :alt    altitude, km
    :lat    geodetic latitude, degrees
    :lon    longitude, degrees east
    :lst    local apparent solar time, hours -- sec/3600 + lon/15 for a
            consistent set, though the model takes the three separately
    :f107a  81-day average of F10.7, centered on the day
    :f107   the previous day's F10.7
    :ap     the daily Ap -- or a vector of seven: the daily Ap, the 3-hour
            ap now and 3, 6 and 9 hours before, and the averages of the
            eight 3-hour values 12 to 33 and 36 to 57 hours before, with
            which the model uses its 3-hour ap terms instead

  returning number densities, m^-3, `:He :O :N2 :O2 :Ar :H :N` and
  `:anomalous-O`; the mass density, kg/m^3, `:rho` of all but anomalous
  oxygen (GTD7's) and `:rho-drag` with it (GTD7D's, the one for drag above
  500 km); and the exospheric temperature `:t-exo` and temperature `:t`,
  K. Option `:switches`, 25 as the Fortran's TSELEC takes them, turns
  variations off (default all on, switch 9 -1 when a 3-hour ap history
  is given)."
  ([inputs] (atmosphere inputs {}))
  ([{:keys [ap] :as inputs} {:keys [switches]}]
   (let [sv (or switches (cond-> default-switches (not (number? ap)) (assoc 8 -1.0)))
         {[d1 d2 d3 d4 d5 d6 d7 d8 d9] :d :keys [tinf t]} (gtd7 inputs sv)
         m3 1e6]
     {:He (* d1 m3) :O (* d2 m3) :N2 (* d3 m3) :O2 (* d4 m3) :Ar (* d5 m3)
      :H (* d7 m3) :N (* d8 m3) :anomalous-O (* d9 m3)
      :rho (* d6 1e3)
      :rho-drag (* 1.66e-24 1e3 (+ (* 4.0 d1) (* 16.0 d2) (* 28.0 d3) (* 32.0 d4) (* 40.0 d5)
                                   d7 (* 14.0 d8) (* 16.0 d9)))
      :t-exo tinf :t t})))
