(ns allgo.astro.earth-radiation
  "Earth radiation pressure: the sunlight the Earth reflects (its albedo)
  and the heat it emits in the infrared, pushing on a satellite (Knocke,
  Ries and Tapley, \"Earth radiation pressure effects on satellites\",
  AIAA 88-4292, 1988; Vallado, chapter 8).

  Both are summed over the part of the Earth the satellite sees, divided
  into elements. Each sunlit element reflects as a diffuse (Lambertian)
  surface, its exitance its albedo times the sunlight falling on it; every
  element emits its emissivity times a quarter of the solar flux, the
  Earth's mean absorbed-and-reradiated share. An element of exitance M
  and area dA puts the irradiance M cos(its angle to the satellite) dA /
  (pi d^2) at a satellite d away, along the line from it. The albedo and
  emissivity are Knocke's zonal models, varying with latitude and season:

    a = 0.34 + a1 P1(sin lat) + 0.29 P2(sin lat),   a1 = 0.10 cos w(t - t0)
    e = 0.68 + e1 P1(sin lat) - 0.18 P2(sin lat),   e1 = -0.07 cos w(t - t0)

  w the annual frequency and t0 1981 December 22. As with direct
  sunlight, the acceleration is the irradiance over c, times the
  reflectivity coefficient and area-to-mass ratio -- here for a sphere,
  the force along each line of sight.

  Positions Earth-centered, inertial, km; the Sun's position gives the
  sunlight's direction; accelerations km/s^2."
  (:require [allgo.astro.constants :as c]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(def ^:private solar-flux
  "The solar irradiance at 1 AU, W/m^2: the radiation pressure times c."
  (* c/solar-pressure 299792458.0))

(def knocke
  "Knocke's zonal albedo and emissivity models: `(fn [lat mjd])` each."
  (let [w (/ (* 2.0 math/PI) 365.25)
        t0 44960.0                                              ; 1981 December 22, MJD
        p1 (fn [lat] (math/sin lat))
        p2 (fn [lat] (* 0.5 (- (* 3.0 (math/pow (math/sin lat) 2)) 1.0)))]
    {:albedo (fn [lat mjd] (+ 0.34 (* 0.10 (math/cos (* w (- mjd t0))) (p1 lat)) (* 0.29 (p2 lat))))
     :emissivity (fn [lat mjd] (+ 0.68 (* -0.07 (math/cos (* w (- mjd t0))) (p1 lat)) (* -0.18 (p2 lat))))}))

(defn- cap-elements
  "The elements of the Earth's surface the satellite at `r` sees: unit
  normals and areas (km^2) of `rings` rings about the sub-satellite point,
  each cut into `sectors`, out to the horizon."
  [r rings sectors]
  (let [rm (v3/length r)
        re c/R-earth
        zhat (v3/scale r (/ 1.0 rm))
        xhat (v3/normalize (v3/cross zhat (if (< (abs (nth zhat 2)) 0.9) [0.0 0.0 1.0] [1.0 0.0 0.0])))
        yhat (v3/cross zhat xhat)
        horizon (math/acos (/ re rm))
        dth (/ horizon rings)]
    (for [k (range rings)
          :let [th (* (+ k 0.5) dth)
                area (* re re (- (math/cos (* k dth)) (math/cos (* (inc k) dth))) (/ (* 2.0 math/PI) sectors))]
          j (range sectors)
          :let [ph (* (/ (+ j 0.5) sectors) 2.0 math/PI)
                n (v3/add (v3/scale zhat (math/cos th))
                          (v3/scale (v3/add (v3/scale xhat (math/cos ph)) (v3/scale yhat (math/sin ph)))
                                    (math/sin th)))]]
      {:n n :area area})))

(defn irradiance
  "The Earth's reflected and emitted irradiance at `r`, W/m^2, as vectors
  along which it pushes: `{:albedo :infrared}`. `models` gives the albedo
  and emissivity as functions of latitude (rad) and MJD, by default
  `knocke`; `:rings` and `:sectors` the grid (default 30 and 36)."
  ([r r-sun mjd] (irradiance r r-sun mjd {}))
  ([r r-sun mjd {:keys [models rings sectors] :or {models knocke rings 30 sectors 36}}]
   (let [{:keys [albedo emissivity]} models
         sun (v3/normalize r-sun)
         flux (* solar-flux (math/pow (/ c/AU (v3/length r-sun)) 2))]
     (reduce (fn [acc {:keys [n area]}]
               (let [p (v3/scale n c/R-earth)
                     d (v3/sub r p)
                     dm (v3/length d)
                     dhat (v3/scale d (/ 1.0 dm))
                     cos-out (v3/dot n dhat)]
                 (if (<= cos-out 0.0)
                   acc
                   (let [lat (math/asin (nth n 2))
                         ;; M cos / (pi d^2) dA, d and dA in the same units
                         view (/ (* cos-out area) (* math/PI dm dm))
                         cos-in (v3/dot n sun)
                         reflected (if (pos? cos-in) (* (albedo lat mjd) flux cos-in view) 0.0)
                         emitted (* (emissivity lat mjd) 0.25 flux view)]
                     (-> acc
                         (update :albedo v3/add (v3/scale dhat reflected))
                         (update :infrared v3/add (v3/scale dhat emitted)))))))
             {:albedo [0.0 0.0 0.0] :infrared [0.0 0.0 0.0]}
             (cap-elements r rings sectors)))))

(defn acceleration
  "The acceleration, km/s^2, of the Earth's albedo and infrared
  radiation on a satellite of area-to-mass ratio `area-to-mass` (m^2/kg)
  and reflectivity coefficient `cr`: the irradiance over c. Options as
  `irradiance` takes them."
  ([r r-sun mjd area-to-mass cr] (acceleration r r-sun mjd area-to-mass cr {}))
  ([r r-sun mjd area-to-mass cr opts]
   (let [{:keys [albedo infrared]} (irradiance r r-sun mjd opts)]
     (v3/scale (v3/add albedo infrared) (/ (* cr area-to-mass 1e-3) 299792458.0)))))
