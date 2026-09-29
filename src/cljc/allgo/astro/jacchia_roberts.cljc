(ns allgo.astro.jacchia-roberts
  "Roberts's analytic form of Jacchia's static models (Roberts, \"An
  analytic model for upper atmosphere densities based upon Jacchia's 1970
  models\", Celestial Mechanics 4, 1971, 368-377): the same densities as
  `allgo.astro.jacchia`'s 1970 models, but in closed form rather than by
  numerical integration -- identical to them from 90 to 125 km, and
  within some 5% above.

  From 90 to 125 km the temperature is Jacchia's quartic, which factors
  as c4 (Z - r1)(Z - r2)((Z - x)^2 + y^2), two real roots and a complex
  pair; the barometric equation below 105 km and the diffusion equation
  above, with gravity falling as (Z + Ra)^-2, then integrate by partial
  fractions into logarithms and an arctangent (Roberts's 11-13 and
  19-20). The coefficients are found here from the residues at the four
  roots -- the same partial fractions Roberts writes out, without his
  pre-multiplied constants, the roots by `allgo.numerics.polynomial`. Above 125 km Roberts replaces Jacchia's
  arctangent profile by one the diffusion equation integrates exactly,

    T = T-inf - (T-inf - Tx) exp[-((Tx - T0)/(T-inf - Tx))((Z - Zx)/(Zx - Z0))(l/(Ra + Z))],

  so that each constituent's density is a power of T and of T-inf - T
  (25-27); the length l, which Roberts fits to Jacchia's densities, is
  his Table II's values at 600, 1300 and 2000 K by Lagrange
  interpolation.

  Roberts labels O2's density at 105 km with O's formula and O's with
  O2's (his 16 and 17); they are taken here the way round Jacchia's
  mixing gives them. Heights km, densities g/cm^3 and cm^-3."
  (:require [allgo.astro.jacchia :as j]
            [allgo.numerics.complex :as cx]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

(def ^:private z0 90.0)
(def ^:private zx 125.0)
(def ^:private t0 183.0)
(def ^:private rho0 3.46e-9)
(def ^:private ra 6356.766)                                 ; km
(def ^:private g0 9.80665)                                  ; m/s^2
(def ^:private rgas 8.31432)                                ; J/(K mol)
(def ^:private m0 28.960)
(def ^:private avogadro 6.02257e23)

;; ------------------------------------------------------ complex numbers

;; ---------------------------------------------------- the quartic profile

(defn- profile
  "What Roberts's analytic model needs of exospheric temperature `tinf`:
  Tx, the quartic's leading coefficient and its roots -- r1 > r2 real, x
  +/- iy the complex pair -- in Z."
  [tinf]
  (let [tx (j/inflection-temperature j/j70 tinf)
        d (- tx t0)
        ;; T = Tx + c1 u + c3 u^3 + c4 u^4, u = Z - 125: equation (5)
        c1 (/ (* 1.9 d) 35.0) c3 (/ (* -1.7 d) (math/pow 35.0 3)) c4 (/ (* -0.8 d) (math/pow 35.0 4))
        us (poly/quartic c4 c3 0.0 c1 tx)
        real (sort > (map first (filter #(zero? (second %)) us)))
        [cx cy] (first (filter #(pos? (second %)) us))]
    {:tx tx :c4 c4 :r1 (+ zx (first real)) :r2 (+ zx (second real)) :x (+ zx cx) :y cy}))

(defn- quartic-temperature [{:keys [tx c4 r1 r2 x y]} z]
  (+ tx (* c4 (- (* (- z r1) (- z r2) (+ (* (- z x) (- z x)) (* y y)))
                 (* (- zx r1) (- zx r2) (+ (* (- zx x) (- zx x)) (* y y)))))))

(defn- fraction-integral
  "The integral from `za` to `zb` of N(Z)/((Z + Ra)^2 (Z - r1)(Z - r2)((Z - x)^2 + y^2)),
  N a polynomial of degree at most 6 given by its coefficients in powers
  of (Z - 100), highest first: by partial fractions, the residues at the
  roots and the double pole at -Ra."
  [{:keys [r1 r2 x y]} coeffs za zb]
  (let [horner (fn [z] (reduce (fn [acc c] (cx/+ (cx/* acc (cx/- z [100.0 0.0])) [c 0.0])) [0.0 0.0] coeffs))
        n (fn [v] (first (horner [v 0.0])))
        n' (fn [v] (let [k (dec (count coeffs))]
                     (first (reduce (fn [[acc i] c] [(+ (* acc (- v 100.0)) (* i c)) (dec i)])
                                    [0.0 k] (butlast coeffs)))))
        lead (if (= 7 (count coeffs)) (first coeffs) 0.0)
        q (fn [v] (+ (* (- v x) (- v x)) (* y y)))
        ;; V = (Z - r1)(Z - r2) q(Z), the denominator without the double pole
        vf (fn [v] (* (- v r1) (- v r2) (q v)))
        vf' (fn [v] (+ (* (- v r2) (q v)) (* (- v r1) (q v)) (* 2.0 (- v r1) (- v r2) (- v x))))
        p5 (/ (n (- ra)) (vf (- ra)))
        p1 (/ (- (* (n' (- ra)) (vf (- ra))) (* (n (- ra)) (vf' (- ra)))) (math/pow (vf (- ra)) 2))
        p2 (/ (n r1) (* (math/pow (+ r1 ra) 2) (- r1 r2) (q r1)))
        p3 (/ (n r2) (* (math/pow (+ r2 ra) 2) (- r2 r1) (q r2)))
        c [x y]
        [re im] (cx// (horner c)
                      (reduce cx/* [(cx/* (cx/+ c [ra 0.0]) (cx/+ c [ra 0.0]))
                                    (cx/- c [r1 0.0]) (cx/- c [r2 0.0]) [0.0 (* 2.0 y)]]))
        ln-ratio (fn [f] (math/log (/ (f zb) (f za))))]
    (+ (* lead (- zb za))
       (* p1 (ln-ratio #(+ % ra)))
       (* p5 (- (/ 1.0 (+ za ra)) (/ 1.0 (+ zb ra))))
       (* p2 (ln-ratio #(- % r1)))
       (* p3 (ln-ratio #(- % r2)))
       (* re (ln-ratio q))
       (* -2.0 im (- (math/atan (/ (- zb x) y)) (math/atan (/ (- za x) y)))))))

(def ^:private mass-coefficients
  "Jacchia's 1970 mean molecular mass, 90 to 105 km, in powers of (Z -
  100), highest first: Roberts's (6)."
  [9.9826e-8 1.5044e-6 -1.0210e-5 -1.0056e-5 1.2840e-4 -8.5586e-2 28.15204])

(defn- mean-mass [z] (reduce (fn [acc c] (+ (* acc (- z 100.0)) c)) 0.0 mass-coefficients))

;; ---------------------------------------------------------- above 125 km

(defn length-l
  "l, km, the length in the temperature profile above 125 km, as a
  function of exospheric temperature: Lagrange interpolation through
  Roberts's Table II."
  [tinf]
  (let [pts [[600.0 11825.0] [1300.0 13515.0] [2000.0 14515.0]]]
    (reduce + (for [[xi yi] pts]
                (* yi (reduce * (for [[xj] pts :when (not= xi xj)] (/ (- tinf xj) (- xi xj)))))))))

(defn- upper-temperature [tinf tx l z]
  (- tinf (* (- tinf tx) (math/exp (- (* (/ (- tx t0) (- tinf tx)) (/ (- z zx) (- zx z0)) (/ l (+ ra z))))))))

(defn temperature
  "The temperature, K, at height `z` km of Roberts's profile for
  exospheric temperature `tinf`: Jacchia's quartic to 125 km, Roberts's
  exponential above (23)."
  [tinf z]
  (let [{:keys [tx] :as pr} (profile tinf)]
    (if (<= z zx) (quartic-temperature pr z) (upper-temperature tinf tx (length-l tinf) z))))

;; ---------------------------------------------------------------- model

(defn static
  "Roberts's model at exospheric temperature `tinf` and height `z` km, 90
  and above: `{:t :rho :n :m}` as `allgo.astro.jacchia/static` gives it --
  hydrogen counted above 500 km."
  [tinf z]
  (let [{:keys [tx c4] :as pr} (profile tinf)
        ;; the barometric and diffusion integrals' common factor, g0 Ra^2/(R c4)
        k (/ (* g0 ra ra) (* rgas c4))
        rho-mixed (fn [z] (* rho0 (/ (* (mean-mass z) t0) (* (mean-mass z0) (quartic-temperature pr z)))
                             (math/exp (* (- k) (fraction-integral pr mass-coefficients z0 z)))))
        mixed-n (fn [z] (let [rho (rho-mixed z) m (mean-mass z)
                              n (/ (* avogadro rho) m) f (/ m m0)
                              q (:q j/j70)]
                          {:N2 (* (q :N2) f n) :Ar (* (q :Ar) f n) :He (* (q :He) f n)
                           :O (* 2.0 n (- 1.0 f)) :O2 (* n (- (* f (+ 1.0 (q :O2))) 1.0))}))
        alpha #(get-in j/species [% :alpha])
        mass #(get-in j/species [% :m])
        ;; 105 to 125 km, diffusion through the quartic: equation (20)
        diffuse-lower (fn [n105 z]
                        (let [t105 (quartic-temperature pr 105.0)
                              tz (quartic-temperature pr z)
                              jz (fraction-integral pr [1.0] 105.0 z)]
                          (into {} (for [[s nb] n105]
                                     [s (* nb (math/pow (/ t105 tz) (+ 1.0 (alpha s)))
                                           (math/exp (* (- k) (mass s) jz)))]))))
        l (length-l tinf)
        gamma (fn [s] (* (/ (* (mass s) g0 ra ra) (* rgas l tinf))
                         (/ (- tinf tx) (- tx t0)) (/ (- zx z0) (+ ra zx))))
        ;; above 125 km, Roberts's profile: equations (25) and (27)
        diffuse-upper (fn [s nb tb tz]
                        (* nb (math/pow (/ tb tz) (+ 1.0 (alpha s) (gamma s)))
                           (math/pow (/ (- tinf tz) (- tinf tb)) (gamma s))))
        tz (temperature tinf z)
        n (cond
            (<= z 105.0) (mixed-n z)
            (<= z zx) (diffuse-lower (mixed-n 105.0) z)
            :else (let [n125 (diffuse-lower (mixed-n 105.0) zx)]
                    (into {} (for [[s nb] n125] [s (diffuse-upper s nb tx tz)]))))
        n (if (> z 500.0)
            (let [lt (math/log10 tinf)
                  t500 (upper-temperature tinf tx l 500.0)]
              (assoc n :H (diffuse-upper :H (math/pow 10.0 (- 73.13 (* (- 39.4 (* 5.5 lt)) lt))) t500 tz)))
            n)
        rho (j/mass-density n)]
    {:t tz :rho rho :n n :m (/ (* rho avogadro) (reduce + (vals n)))}))

(defn atmosphere
  "Roberts's densities at a time and place: `allgo.astro.jacchia/atmosphere`,
  J71's solar, diurnal, geomagnetic, semiannual, seasonal-latitudinal and
  helium variations, on this model's static densities; the inputs and
  results as it has them."
  [inputs]
  (j/atmosphere inputs {:static-model static}))
