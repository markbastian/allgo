(ns allgo.astro.tides
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
  into `allgo.astro.geopotential` without it knowing anything happened."
  (:require [allgo.astro.constants :as c]
            [clojure.math :as math]))

(def k2
  "Nominal degree-2 Love number. The elastic response of the whole Earth:
  0 would be perfectly rigid, 1.5 perfectly fluid."
  0.30)

(defn- legendre-2
  "The three normalized degree-2 associated Legendre functions at sin(lat).
  Normalized, because published coefficients are and this correction has to
  add to them."
  [x]
  (let [s (math/sqrt (max 0.0 (- 1.0 (* x x))))]
    [(* (math/sqrt 5.0) 0.5 (- (* 3.0 x x) 1.0))            ; P20
     (* (math/sqrt (/ 5.0 3.0)) 3.0 x s)                    ; P21
     (* (math/sqrt (/ 5.0 12.0)) 3.0 s s)]))                ; P22

(defn corrections
  "Degree-2 corrections to the normalized harmonic coefficients, from a
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

(defn perturb
  "Add tidal corrections to a normalized gravity field."
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

;; ---------------------------------------------------------------- the ocean

(def k-ocean
  "Effective degree-2 response of the oceans, as a Love number comparable to
  `k2`. The oceans contribute roughly a tenth of what the solid Earth does."
  0.032)

(def ocean-lead
  "How far the ocean bulge runs *ahead* of the body raising it, radians.

  Not a lag. The Earth turns once a day while the Moon takes twenty-seven,
  so rotation drags the bulge past the sub-lunar point before friction and
  inertia can settle it. The couple that results is the reason the Moon
  recedes by 38 mm a year and the day lengthens by 2 ms a century -- the
  Earth is spinning down and handing the angular momentum to the Moon.

  A few degrees, and the sign is the whole physics: reverse it and the Moon
  would be spiraling in."
  (* 3.0 c/degrees))

(defn- lead-longitude
  "Rotate a body eastward about the pole, into the direction the Earth turns."
  [[x y z] angle]
  (let [ca (math/cos angle) sa (math/sin angle)]
    [(- (* ca x) (* sa y)) (+ (* sa x) (* ca y)) z]))

(defn ocean-corrections
  "Degree-2 harmonic corrections from ocean tides, from `[GM position]` pairs
  in the Earth-fixed frame.

  An equilibrium model: the oceans are treated as responding to the same
  tide-raising potential as the solid Earth, at about a tenth the strength
  and running a few degrees ahead. A real model carries a tabulated ocean
  tide -- FES or GOT -- resolved into Doodson constituents, because the true
  response depends on basin geometry and resonance rather than on the
  forcing alone. What is here has the right size and the right lead; it does
  not have the geography."
  ([bodies] (ocean-corrections bodies k-ocean ocean-lead))
  ([bodies love lead]
   (corrections (mapv (fn [[GM r]] [GM (lead-longitude r lead)]) bodies) love)))

(defn perturb-with-ocean
  "Add both solid and ocean tidal corrections to a normalized gravity field."
  [field bodies]
  (let [solid (corrections bodies)
        ocean (ocean-corrections bodies)
        dC    (merge-with + (:C solid) (:C ocean))
        dS    (merge-with + (:S solid) (:S ocean))]
    (-> field
        (update :C #(merge-with + (merge {[2 0] 0.0 [2 1] 0.0 [2 2] 0.0} %) dC))
        (update :S #(merge-with + (merge {[2 1] 0.0 [2 2] 0.0} %) dS)))))
