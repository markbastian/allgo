(ns allgo.astro.jacchia
  "Jacchia's 1971 model of the thermosphere and exosphere (Jacchia,
  \"Revised static models of the thermosphere and exosphere with
  empirical temperature profiles\", Smithsonian Astrophysical Observatory
  Special Report 332, 1971) -- J71, the model Vallado's chapter 8 builds
  its drag on, and CIRA 1972's.

  The model is static: the atmosphere at any time and place is one of a
  family of profiles fixed by a single number, the exospheric temperature
  T-inf, above constant conditions at 90 km -- 183 K and 3.46e-9 g/cm^3.
  Each profile's temperature rises from 90 km as a quartic to an
  inflection at 125 km and on toward T-inf as an arctangent; below 100
  km the air is mixed, its mean molecular mass falling with the
  dissociation of oxygen, and above it each constituent -- N2, O2, O, Ar,
  He, and H from its own level at 500 km -- diffuses in the gravity field
  alone. The equations are the report's, numbered as it numbers them;
  the integrals, numerical as the report's were, are Simpson's rule.

  The time, place, solar and geomagnetic activity enter through T-inf --
  the nighttime minimum set by the 10.7 cm flux (14), the diurnal
  variation over the globe (15-17), the geomagnetic heating (18) -- and
  through corrections to the densities themselves: the semiannual
  variation (21-23), the seasonal-latitudinal variations of the lower
  thermosphere (24) and of helium (25), and the density part of the
  geomagnetic effect in the hybrid form (20).

  Heights are geometric, km; densities g/cm^3 and number densities cm^-3,
  as the report tabulates them, except where `atmosphere` says otherwise."
  (:require [allgo.astro.solar :as solar]
            [allgo.numerics.quadrature :as quadrature]
            [clojure.math :as math]))

;; ------------------------------------------------------------ constants

(def ^:private z0 90.0)                                     ; km
(def ^:private t0 183.0)                                    ; K at z0
(def ^:private zx 125.0)                                    ; km, the inflection
(def ^:private rho0 3.46e-9)                                ; g/cm^3 at z0
(def ^:private avogadro 6.02257e23)
(def ^:private rgas 8.31432e7)                              ; erg/(K mol)
(def ^:private g0 980.665)                                  ; cm/s^2
(def ^:private re 6356.766)                                 ; km
(def ^:private m0 28.960)                                   ; sea-level mean molecular mass

(def species
  "Molecular mass and thermal diffusion coefficient of each constituent
  (section 2 and equation 6)."
  {:N2 {:m 28.0134 :alpha 0.0}
   :O2 {:m 31.9988 :alpha 0.0}
   :Ar {:m 39.948 :alpha 0.0}
   :He {:m 4.0026 :alpha -0.38}
   :O {:m 15.9994 :alpha 0.0}
   :H {:m 1.00797 :alpha 0.0}})

(defn- gravity
  "Equation (8), cm/s^2."
  [z]
  (* g0 (math/pow (+ 1.0 (/ z re)) -2)))

(defn mean-molecular-mass
  "The mean molecular mass of the mixed air from 90 to 100 km, equation (1)."
  [z]
  (let [x (- z z0)]
    (reduce (fn [acc c] (+ (* acc x) c)) 0.0
            [-6.97444e-7 1.07561e-5 -8.21895e-6 4.51103e-4 -1.19407e-2 -7.40066e-2 28.82678])))

(def j71
  "The J71 models' constants: the inflection temperature's coefficients
  (9), the homopause, the mixed air's molecular mass (1) and sea-level
  composition by volume, and hydrogen -- counted from the homopause up,
  from its density at 500 km fixed by the temperature there (7)."
  {:tx [371.6678 0.0518806 -294.3505 -0.00216222]
   :homopause 100.0
   :mass mean-molecular-mass
   :q {:N2 0.78110 :O2 0.20955 :Ar 0.0093432 :He 0.0000061471}
   :hydrogen {:from 100.0 :temperature :t500}})

(def j70
  "The constants of Jacchia's 1970 models (SAO Special Report 313), as
  Roberts (Celestial Mechanics 4, 1971) restates them: the air mixed to
  105 km, more helium, and hydrogen's density at 500 km fixed by the
  exospheric temperature -- counted, as in Jacchia's tables, from the
  homopause up."
  {:tx [444.3807 0.02385 -392.8292 -0.0021357]
   :homopause 105.0
   :mass (fn [z] (let [x (- z 100.0)]
                   (reduce (fn [acc c] (+ (* acc x) c)) 0.0
                           [9.9826e-8 1.5044e-6 -1.0210e-5 -1.0056e-5 1.2840e-4 -8.5586e-2 28.15204])))
   :q {:N2 0.78110 :O2 0.20955 :Ar 0.00934 :He 0.00001289}
   :hydrogen {:from 105.0 :temperature :tinf}})

;; ---------------------------------------------------- temperature profile

(defn inflection-temperature
  "Tx, the temperature at the inflection, 125 km, of the profile with
  exospheric temperature `tinf`: equation (9)."
  ([tinf] (inflection-temperature j71 tinf))
  ([{[a b c d] :tx} tinf] (+ a (* b tinf) (* c (math/exp (* d tinf))))))

(defn inflection-gradient
  "Gx, the gradient at the inflection, K/km: 1.90 (Tx - T0)/(zx - z0),
  equation (11)."
  ([tinf] (inflection-gradient j71 tinf))
  ([params tinf] (/ (* 1.9 (- (inflection-temperature params tinf) t0)) (- zx z0))))

(defn temperature
  "The temperature, K, at height `z` km of the profile with exospheric
  temperature `tinf`: below 125 km the quartic (10) through 183 K at 90
  km with no gradient there and the inflection's value, gradient and
  zero curvature at 125 (11); above, the arctangent (13)."
  ([tinf z] (temperature j71 tinf z))
  ([params tinf z]
   (let [tx (inflection-temperature params tinf)
         gx (inflection-gradient params tinf)]
     (if (< z zx)
       ;; with c1 = Gx and c2 = 0, T(z0) = T0 and T'(z0) = 0 fix c3 and c4
       (let [d (- z0 zx)
             ;; c3 d^3 + c4 d^4 = T0 - Tx - Gx d;  3 c3 d^2 + 4 c4 d^3 = -Gx
             r1 (- t0 tx (* gx d))
             r2 (- gx)
             c4 (/ (- (* r2 d) (* 3.0 r1)) (math/pow d 4))
             c3 (/ (- r1 (* c4 (math/pow d 4))) (math/pow d 3))
             x (- z zx)]
         (+ tx (* gx x) (* c3 x x x) (* c4 x x x x)))
       (let [a (* (/ 2.0 math/PI) (- tinf tx))
             x (- z zx)]
         (+ tx (* a (math/atan (* (/ gx a) x (+ 1.0 (* 4.5e-6 (math/pow x 2.5))))))))))))

;; ------------------------------------------------------------- densities

(defn- integral
  "The integral of `f` from `a` to `b` by Simpson's rule, split at the
  inflection, where the profile changes form, at steps of some 0.1 km."
  [f a b]
  (let [piece (fn [a b] (if (<= b a) 0.0 (quadrature/simpson f a b (max 2 (long (math/ceil (* 10.0 (- b a))))))))]
    (if (< a zx b)
      (+ (piece a zx) (piece zx b))
      (piece a b))))

(defn- mixed-density
  "The mass density, g/cm^3, of the mixed air below the homopause: the
  barometric equation (5) integrated from 90 km."
  [{:keys [mass] :as params} tinf z]
  (let [tz (temperature params tinf z)
        ;; the integral of M g/(R* T) dz, dz in cm
        path (* 1e5 (integral #(/ (* (mass %) (gravity %)) (* rgas (temperature params tinf %))) z0 z))]
    (* rho0 (/ (* (mass z) t0) (* (mass z0) tz)) (math/exp (- path)))))

(defn- mixed-number-densities
  "Number densities, cm^-3, of the mixed air of density `rho` and mean
  molecular mass `m`: equations (2) to (4)."
  [{:keys [q]} rho m]
  (let [n (/ (* avogadro rho) m)
        f (/ m m0)]
    {:N2 (* (q :N2) f n)
     :Ar (* (q :Ar) f n)
     :He (* (q :He) f n)
     :O (* 2.0 n (- 1.0 f))
     :O2 (* n (- (* f (+ 1.0 (q :O2))) 1.0))}))

(defn hydrogen-500
  "The number density of H at 500 km, cm^-3, for the temperature `t`:
  equation (7)."
  [t]
  (let [l (math/log10 t)]
    (math/pow 10.0 (+ 73.13 (* -39.40 l) (* 5.5 l l)))))

(defn mass-density
  "The mass density, g/cm^3, of number densities `n`, cm^-3."
  [n]
  (/ (reduce + (map (fn [[k v]] (* v (get-in species [k :m]))) n)) avogadro))

(defn static
  "The static model at exospheric temperature `tinf` and height `z` km,
  90 and above: `{:t :rho :n {:N2 :O2 :O :Ar :He :H} :m}` -- temperature
  (K), mass density (g/cm^3), number densities (cm^-3) and mean molecular
  mass. With `params` `j70`, Jacchia's 1970 models instead of J71. In J71
  hydrogen is counted above the homopause, as the report computed it,
  though it tabulates it only from 500 km."
  ([tinf z] (static j71 tinf z))
  ([{:keys [homopause mass hydrogen] :as params} tinf z]
   (let [temp #(temperature params tinf %)
         tz (temp z)
         g-over-t #(/ (gravity %) (temp %))
         ;; the integral of g/T from a base to z, common to every constituent
         path (fn [z-base] (if (>= z z-base) (integral g-over-t z-base z) (- (integral g-over-t z z-base))))
         diffuse (fn [n-base t-base path k]
                   (let [{:keys [m alpha]} (species k)]
                     (* n-base (math/pow (/ t-base tz) (+ 1.0 alpha))
                        (math/exp (/ (* -1e5 m path) rgas)))))
         t500 (temp 500.0)
         nh (when (> z (:from hydrogen))
              (diffuse (hydrogen-500 (if (= :tinf (:temperature hydrogen)) tinf t500)) t500 (path 500.0) :H))
         n (if (<= z homopause)
             (mixed-number-densities params (mixed-density params tinf z) (mass z))
             (let [tb (temp homopause)
                   base (mixed-number-densities params (mixed-density params tinf homopause) (mass homopause))
                   pb (path homopause)]
               (into {} (map (fn [k] [k (diffuse (base k) tb pb k)]) (keys base)))))
         n (if nh (assoc n :H nh) n)
         rho (mass-density n)]
     {:t tz :rho rho :n n :m (/ (* rho avogadro) (reduce + (vals n)))})))

;; ------------------------------------------------------------ variations

(defn night-minimum
  "Tc, the global nighttime minimum exospheric temperature, K, at Kp = 0:
  equation (14), from the previous day's 10.7 cm flux `f107` (the
  temperature lags the flux by a day) and its average over three solar
  rotations `f107a`, in 1e-22 W m^-2 Hz^-1."
  [f107 f107a]
  (+ 379.0 (* 3.24 f107a) (* 1.3 (- f107 f107a))))

(defn local-temperature
  "Tl, the uncorrected exospheric temperature at latitude `lat`, the Sun's
  declination `dec` (radians) and local solar time `lst` (hours), for
  nighttime minimum `tc`: equations (15) to (17), the maximum lagging the
  subsolar point by the terms in beta, p and gamma. J71's constants m =
  2.2 and R = 0.3 by default; `constants` `{:m :r}` for others -- J70's
  are m = 2.5 and, on average, R = 0.31."
  ([tc lat dec lst] (local-temperature tc lat dec lst {}))
  ([tc lat dec lst {:keys [m r] :or {m 2.2 r 0.3}}]
   (let [n 3.0
         beta (math/to-radians -37.0) p (math/to-radians 6.0) gamma (math/to-radians 43.0)
         h (math/to-radians (* 15.0 (- lst 12.0)))
         tau (let [t (+ h beta (* p (math/sin (+ h gamma))))]
              ;; into -pi to pi
               (- t (* 2.0 math/PI (math/floor (/ (+ t math/PI) (* 2.0 math/PI))))))
         theta (* 0.5 (abs (+ lat dec)))
         eta (* 0.5 (abs (- lat dec)))
         s (math/pow (math/sin theta) m)
         c (math/pow (math/cos eta) m)]
     (* tc (+ 1.0 (* r s)) (+ 1.0 (* (/ (* r (- c s)) (+ 1.0 (* r s))) (math/pow (math/cos (* 0.5 tau)) n)))))))

(defn geomagnetic-temperature
  "The geomagnetic rise in exospheric temperature, K, for the planetary
  index `kp` (lagged 6.7 hours): equation (18) at and above 200 km, the
  temperature part (20b) of the hybrid form below."
  [kp z]
  (if (< z 200.0)
    (+ (* 14.0 kp) (* 0.02 (math/exp kp)))
    (+ (* 28.0 kp) (* 0.03 (math/exp kp)))))

(defn geomagnetic-density
  "The density part of the hybrid geomagnetic effect below 200 km,
  Delta log10 rho: equation (20a); nothing above."
  [kp z]
  (if (< z 200.0) (+ (* 0.012 kp) (* 1.2e-5 (math/exp kp))) 0.0))

(defn- semiannual-phase
  "Phi of equation (23), from MJD `t`."
  [t]
  (/ (- t 36204.0) 365.2422))

(defn semiannual-height
  "f(z) of equation (22), the semiannual variation's amplitude at height
  `z` km."
  [z]
  (* (+ (* 5.876e-7 (math/pow z 2.331)) 0.06328) (math/exp (* -2.868e-3 z))))

(defn semiannual-time
  "g(t) of equation (22) at MJD `t`, normalized to unit amplitude."
  [t]
  (let [phi (semiannual-phase t)
        tau (+ phi (* 0.09544 (- (math/pow (+ 0.5 (* 0.5 (math/sin (+ (* 2.0 math/PI phi) 6.035)))) 1.65) 0.5)))]
    (+ 0.02835 (* 0.3817 (+ 1.0 (* 0.4671 (math/sin (+ (* 2.0 math/PI tau) 4.137))))
                  (math/sin (+ (* 4.0 math/PI tau) 4.259))))))

(defn semiannual
  "The semiannual variation, Delta log10 rho = f(z) g(t): equation (21)."
  [z t]
  (* (semiannual-height z) (semiannual-time t)))

(defn seasonal-latitudinal
  "The seasonal-latitudinal variation of the lower thermosphere, Delta
  log10 rho, at height `z` km, latitude `lat` (radians) and MJD `t`:
  equation (24)."
  [z lat t]
  (let [x (- z 90.0)]
    (* 0.014 x (math/exp (* -0.0013 x x)) (math/signum lat)
       (math/sin (+ (* 2.0 math/PI (semiannual-phase t)) 1.72))
       (math/pow (math/sin lat) 2))))

(defn helium
  "The seasonal-latitudinal variation of helium, Delta log10 n(He), at
  latitude `lat` and solar declination `dec` (radians): equation (25),
  the obliquity 23.44 degrees."
  [lat dec]
  (let [eps (math/to-radians 23.44)
        q (/ math/PI 4.0)]
    (* 0.65 (abs (/ dec eps))
       (- (math/pow (math/sin (- q (* 0.5 lat (math/signum dec)))) 3)
          (math/pow (math/sin q) 3)))))

(defn atmosphere
  "J71 at `inputs`, a map of

    :mjd    the time, MJD (UT)
    :alt    height, km, 90 and above
    :lat    geodetic latitude, degrees
    :lst    local solar time, hours
    :f107   the previous day's 10.7 cm flux
    :f107a  its average over three solar rotations
    :kp     the planetary geomagnetic index, 6.7 hours before
    :dec    the Sun's declination, degrees (by default from the Sun's
            low-accuracy position)

  returning number densities, m^-3, `:N2 :O2 :O :Ar :He :H`; the mass
  density, kg/m^3, `:rho`; the temperature `:t` and exospheric
  temperature `:t-exo`, K. The exospheric temperature carries the solar,
  diurnal and geomagnetic variations; the densities the semiannual,
  seasonal-latitudinal and (below 200 km) geomagnetic corrections, which
  scale every constituent alike, and helium its own. Option
  `:static-model`, a function of exospheric temperature and height as
  `static` is, puts the variations on other static densities -- Roberts's,
  as `allgo.astro.jacchia-roberts/atmosphere` does."
  ([inputs] (atmosphere inputs {}))
  ([{:keys [mjd alt lat lst f107 f107a kp dec]} {:keys [static-model] :or {static-model static}}]
   (let [lat (math/to-radians lat)
         dec (if dec (math/to-radians dec) (second (solar/equatorial-low mjd)))
         tinf (+ (local-temperature (night-minimum f107 f107a) lat dec lst)
                 (geomagnetic-temperature kp alt))
         {:keys [t n]} (static-model tinf alt)
         n (update n :He * (math/pow 10.0 (helium lat dec)))
         f (math/pow 10.0 (+ (semiannual alt mjd) (seasonal-latitudinal alt lat mjd) (geomagnetic-density kp alt)))
         n (update-vals n #(* % f 1e6))]
     (assoc n :rho (* 1e3 (/ (mass-density n) 1e6)) :t t :t-exo tinf))))
