(ns procedurals.astro.tides
  "Solid Earth tides as a perturbation of the geopotential
  (Montenbruck & Gill 3.6).

  The Sun and Moon do not only pull on a satellite directly -- they also
  deform the Earth, and the deformed Earth's gravity is not what it was.
  This is a second-order effect of a small one and comes to around 1e-9
  m/s^2, but it is periodic at tidal frequencies rather than random, so it
  accumulates in the same way solar radiation pressure does.

  The response is expressed through Love numbers: the Earth's own induced
  potential is k2 times the tide-raising one, with k2 near 0.3 -- the Earth
  is not rigid, but it is far from fluid. The whole effect is expressed as a
  correction to the degree-2 harmonic coefficients, which then feed straight
  into `procedurals.astro.geopotential` without it knowing anything happened."
  (:require [procedurals.astro.constants :as c]
            [clojure.math :as math]))

(def k2
  "Nominal degree-2 Love number. The elastic response of the whole Earth:
  0 would be perfectly rigid, 1.5 perfectly fluid."
  0.30)

(defn- legendre-2
  "The three normalised degree-2 associated Legendre functions at sin(lat).
  Normalised, because published coefficients are and this correction has to
  add to them."
  [x]
  (let [s (math/sqrt (max 0.0 (- 1.0 (* x x))))]
    [(* (math/sqrt 5.0) 0.5 (- (* 3.0 x x) 1.0))            ; P20
     (* (math/sqrt (/ 5.0 3.0)) 3.0 x s)                    ; P21
     (* (math/sqrt (/ 5.0 12.0)) 3.0 s s)]))                ; P22

(defn corrections
  "Degree-2 corrections to the normalised harmonic coefficients, from a
  sequence of `[GM position]` pairs in the Earth-fixed frame.

  Returns `{:C {[2 0] .. [2 1] .. [2 2] ..} :S {[2 1] .. [2 2] ..}}`, ready to
  be merged into a field. Following the IERS form, where the tide-raising
  potential of each body is scaled by k2/5 and resolved onto the harmonics
  by its own latitude and longitude."
  ([bodies] (corrections bodies k2))
  ([bodies love]
   (reduce
    (fn [acc [GM r]]
      (let [[x y z] r
            d   (math/sqrt (+ (* x x) (* y y) (* z z)))
            lat (math/asin (/ z d))
            lon (math/atan2 y x)
            ;; the tide-raising strength: mass ratio times the cube of the
            ;; radius ratio, which is why the Moon beats the Sun three to one
            amp (* (/ love 5.0) (/ GM c/GM-earth)
                   (let [q (/ c/R-earth d)] (* q q q)))
            [p20 p21 p22] (legendre-2 (math/sin lat))]
        (-> acc
            (update-in [:C [2 0]] + (* amp p20))
            (update-in [:C [2 1]] + (* amp p21 (math/cos lon)))
            (update-in [:S [2 1]] + (* amp p21 (math/sin lon)))
            (update-in [:C [2 2]] + (* amp p22 (math/cos (* 2.0 lon))))
            (update-in [:S [2 2]] + (* amp p22 (math/sin (* 2.0 lon)))))))
    {:C {[2 0] 0.0 [2 1] 0.0 [2 2] 0.0} :S {[2 1] 0.0 [2 2] 0.0}}
    bodies)))

(defn apply-to
  "Add tidal corrections to a normalised gravity field."
  [field bodies]
  (let [{:keys [C S]} (corrections bodies)]
    (-> field
        (update :C #(merge-with + (merge {[2 0] 0.0 [2 1] 0.0 [2 2] 0.0} %) C))
        (update :S #(merge-with + (merge {[2 1] 0.0 [2 2] 0.0} %) S)))))

(defn raising-potential
  "The tide-raising potential of a body at the Earth's surface, km^2/s^2 --
  the degree-2 term of its pull, differenced across the Earth.

  Kept so that the induced correction can be checked against it: whatever
  the sign conventions, the Earth's response must come to k2 times this."
  [GM r-body surface-point]
  (let [d  (math/sqrt (reduce + (map * r-body r-body)))
        rs (math/sqrt (reduce + (map * surface-point surface-point)))
        cos-psi (/ (reduce + (map * r-body surface-point)) (* d rs))
        q  (/ rs d)]
    (* (/ GM d) q q 0.5 (- (* 3.0 cos-psi cos-psi) 1.0))))
