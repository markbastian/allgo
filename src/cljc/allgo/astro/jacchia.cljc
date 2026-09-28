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
  (:require [allgo.numerics.quadrature :as quadrature]
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
  "Molecular mass, sea-level fraction by volume and thermal diffusion
  coefficient of each constituent (section 2 and equation 6); O's and
  H's masses are half O2's and the report's H."
  {:N2 {:m 28.0134 :q 0.78110 :alpha 0.0}
   :O2 {:m 31.9988 :q 0.20955 :alpha 0.0}
   :Ar {:m 39.948 :q 0.0093432 :alpha 0.0}
   :He {:m 4.0026 :q 0.0000061471 :alpha -0.38}
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

;; ---------------------------------------------------- temperature profile

(defn inflection-temperature
  "Tx, the temperature at the inflection, 125 km, of the profile with
  exospheric temperature `tinf`: equation (9)."
  [tinf]
  (+ 371.6678 (* 0.0518806 tinf) (* -294.3505 (math/exp (* -0.00216222 tinf)))))

(defn inflection-gradient
  "Gx, the gradient at the inflection, K/km: 1.90 (Tx - T0)/(zx - z0),
  equation (11)."
  [tinf]
  (/ (* 1.9 (- (inflection-temperature tinf) t0)) (- zx z0)))

(defn temperature
  "The temperature, K, at height `z` km of the profile with exospheric
  temperature `tinf`: below 125 km the quartic (10) through 183 K at 90
  km with no gradient there and the inflection's value, gradient and
  zero curvature at 125 (11); above, the arctangent (13)."
  [tinf z]
  (let [tx (inflection-temperature tinf)
        gx (inflection-gradient tinf)]
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
        (+ tx (* a (math/atan (* (/ gx a) x (+ 1.0 (* 4.5e-6 (math/pow x 2.5)))))))))))

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
  "The mass density, g/cm^3, from 90 to 100 km: the barometric equation
  (5) integrated from 90 km."
  [tinf z]
  (let [tz (temperature tinf z)
        mz (mean-molecular-mass z)
        m90 (mean-molecular-mass z0)
        ;; the integral of M g/(R* T) dz, dz in cm
        path (* 1e5 (integral #(/ (* (mean-molecular-mass %) (gravity %)) (* rgas (temperature tinf %))) z0 z))]
    (* rho0 (/ (* mz t0) (* m90 tz)) (math/exp (- path)))))

(defn- mixed-number-densities
  "Number densities, cm^-3, of the mixed air of density `rho` and mean
  molecular mass `m`: equations (2) to (4)."
  [rho m]
  (let [n (/ (* avogadro rho) m)
        f (/ m m0)]
    {:N2 (* (get-in species [:N2 :q]) f n)
     :Ar (* (get-in species [:Ar :q]) f n)
     :He (* (get-in species [:He :q]) f n)
     :O (* 2.0 n (- 1.0 f))
     :O2 (* n (- (* f (+ 1.0 (get-in species [:O2 :q]))) 1.0))}))

(defn- hydrogen-500
  "The number density of H at 500 km, cm^-3, for the temperature there:
  equation (7)."
  [t500]
  (let [l (math/log10 t500)]
    (math/pow 10.0 (+ 73.13 (* -39.40 l) (* 5.5 l l)))))

(defn- mass-density [n]
  (/ (reduce + (map (fn [[k v]] (* v (get-in species [k :m]))) n)) avogadro))

(defn static
  "The static model at exospheric temperature `tinf` and height `z` km,
  90 and above: `{:t :rho :n {:N2 :O2 :O :Ar :He :H} :m}` -- temperature
  (K), mass density (g/cm^3), number densities (cm^-3) and mean molecular
  mass. Hydrogen is counted above 100 km, where the diffusion begins, as
  the report computed it, though it tabulates it only from 500 km."
  [tinf z]
  (let [tz (temperature tinf z)
        g-over-t #(/ (gravity %) (temperature tinf %))
        diffuse (fn [n-base t-base z-base k]
                  (let [{:keys [m alpha]} (species k)]
                    (* n-base (math/pow (/ t-base tz) (+ 1.0 alpha))
                       (math/exp (/ (* -1e5 m (if (>= z z-base)
                                                (integral g-over-t z-base z)
                                                (- (integral g-over-t z z-base))))
                                    rgas)))))
        t500 (temperature tinf 500.0)
        nh (diffuse (hydrogen-500 t500) t500 500.0 :H)
        n (if (<= z 100.0)
            (mixed-number-densities (mixed-density tinf z) (mean-molecular-mass z))
            (let [t100 (temperature tinf 100.0)
                  base (mixed-number-densities (mixed-density tinf 100.0) (mean-molecular-mass 100.0))]
              (into {} (map (fn [k] [k (diffuse (base k) t100 100.0 k)]) (keys base)))))
        n (if (> z 100.0) (assoc n :H nh) n)
        rho (mass-density n)]
    {:t tz :rho rho :n n :m (/ (* rho avogadro) (reduce + (vals n)))}))
