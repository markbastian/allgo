(ns allgo.astro.moon
  "The Moon: its position, its phase, and how it presents its face
  (Meeus, *Astronomical Algorithms*, chapters 47, 48 and 53).

  The position is Chapront's ELP-2000/82 as Meeus abridges it: sixty
  terms each in longitude and distance and sixty in latitude, in the four
  arguments of lunar theory -- the elongation D, the Sun's and Moon's mean
  anomalies M and M', and the argument of latitude F. Good to 10
  arcseconds in longitude and 4 in latitude. The Sun's pull, which is what
  most of those terms are, weakens as the Earth's orbit grows rounder, so
  terms in M are scaled by the factor E for each power of M.

  `allgo.astro.ephemeris/moon` is a coarser series, arcminutes, for the
  force models.

  Times are MJD (TT); angles radians; the distance is in kilometers."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.frames :as frames]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

(defn- deg [x] (* x c/degrees))

(defn arguments
  "The fundamental arguments at `mjd-tt`, radians: the Moon's mean
  longitude `:L'`, the mean elongation `:D`, the Sun's and Moon's mean
  anomalies `:M` and `:M'`, the argument of latitude `:F`, and `:E`, the
  eccentricity factor."
  [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)]
    {:L' (deg (poly/horner-ascending [218.3164477 481267.88123421 -0.0015786 (/ 1.0 538841) (/ -1.0 65194000)] T))
     :D  (deg (poly/horner-ascending [297.8501921 445267.1114034 -0.0018819 (/ 1.0 545868) (/ -1.0 113065000)] T))
     :M  (deg (poly/horner-ascending [357.5291092 35999.0502909 -0.0001535 (/ 1.0 24490000)] T))
     :M' (deg (poly/horner-ascending [134.9633964 477198.8675055 0.0087414 (/ 1.0 69699) (/ -1.0 14712000)] T))
     :F  (deg (poly/horner-ascending [93.272095 483202.0175233 -0.0036539 (/ -1.0 3526000) (/ 1.0 863310000)] T))
     :E  (poly/horner-ascending [1.0 -0.002516 -0.0000074] T)}))

;; ------------------------------------------------------------ position (47)

(def ^:private longitude-distance
  "Meeus's table 47.A: multiples of D, M, M' and F, then the coefficients
  of the sine of the argument in longitude (1e-6 degree) and of its cosine
  in distance (1e-3 km)."
  [[0 0 1 0 6288774 -20905355]
   [2 0 -1 0 1274027 -3699111]
   [2 0 0 0 658314 -2955968]
   [0 0 2 0 213618 -569925]
   [0 1 0 0 -185116 48888]
   [0 0 0 2 -114332 -3149]
   [2 0 -2 0 58793 246158]
   [2 -1 -1 0 57066 -152138]
   [2 0 1 0 53322 -170733]
   [2 -1 0 0 45758 -204586]
   [0 1 -1 0 -40923 -129620]
   [1 0 0 0 -34720 108743]
   [0 1 1 0 -30383 104755]
   [2 0 0 -2 15327 10321]
   [0 0 1 2 -12528 0]
   [0 0 1 -2 10980 79661]
   [4 0 -1 0 10675 -34782]
   [0 0 3 0 10034 -23210]
   [4 0 -2 0 8548 -21636]
   [2 1 -1 0 -7888 24208]
   [2 1 0 0 -6766 30824]
   [1 0 -1 0 -5163 -8379]
   [1 1 0 0 4987 -16675]
   [2 -1 1 0 4036 -12831]
   [2 0 2 0 3994 -10445]
   [4 0 0 0 3861 -11650]
   [2 0 -3 0 3665 14403]
   [0 1 -2 0 -2689 -7003]
   [2 0 -1 2 -2602 0]
   [2 -1 -2 0 2390 10056]
   [1 0 1 0 -2348 6322]
   [2 -2 0 0 2236 -9884]
   [0 1 2 0 -2120 5751]
   [0 2 0 0 -2069 0]
   [2 -2 -1 0 2048 -4950]
   [2 0 1 -2 -1773 4130]
   [2 0 0 2 -1595 0]
   [4 -1 -1 0 1215 -3958]
   [0 0 2 2 -1110 0]
   [3 0 -1 0 -892 3258]
   [2 1 1 0 -810 2616]
   [4 -1 -2 0 759 -1897]
   [0 2 -1 0 -713 -2117]
   [2 2 -1 0 -700 2354]
   [2 1 -2 0 691 0]
   [2 -1 0 -2 596 0]
   [4 0 1 0 549 -1423]
   [0 0 4 0 537 -1117]
   [4 -1 0 0 520 -1571]
   [1 0 -2 0 -487 -1739]
   [2 1 0 -2 -399 0]
   [0 0 2 -2 -381 -4421]
   [1 1 1 0 351 0]
   [3 0 -2 0 -340 0]
   [4 0 -3 0 330 0]
   [2 -1 2 0 327 0]
   [0 2 1 0 -323 1165]
   [1 1 -1 0 299 0]
   [2 0 3 0 294 0]
   [2 0 -1 -2 0 8752]])

(def ^:private latitude-terms
  "Meeus's table 47.B: multiples of D, M, M' and F and the coefficient of
  the sine in latitude, 1e-6 degree."
  [[0 0 0 1 5128122]
   [0 0 1 1 280602]
   [0 0 1 -1 277693]
   [2 0 0 -1 173237]
   [2 0 -1 1 55413]
   [2 0 -1 -1 46271]
   [2 0 0 1 32573]
   [0 0 2 1 17198]
   [2 0 1 -1 9266]
   [0 0 2 -1 8822]
   [2 -1 0 -1 8216]
   [2 0 -2 -1 4324]
   [2 0 1 1 4200]
   [2 1 0 -1 -3359]
   [2 -1 -1 1 2463]
   [2 -1 0 1 2211]
   [2 -1 -1 -1 2065]
   [0 1 -1 -1 -1870]
   [4 0 -1 -1 1828]
   [0 1 0 1 -1794]
   [0 0 0 3 -1749]
   [0 1 -1 1 -1565]
   [1 0 0 1 -1491]
   [0 1 1 1 -1475]
   [0 1 1 -1 -1410]
   [0 1 0 -1 -1344]
   [1 0 0 -1 -1335]
   [0 0 3 1 1107]
   [4 0 0 -1 1021]
   [4 0 -1 1 833]
   [0 0 1 -3 777]
   [4 0 -2 1 671]
   [2 0 0 -3 607]
   [2 0 2 -1 596]
   [2 -1 1 -1 491]
   [2 0 -2 1 -451]
   [0 0 3 -1 439]
   [2 0 2 1 422]
   [2 0 -3 -1 421]
   [2 1 -1 1 -366]
   [2 1 0 1 -351]
   [4 0 0 1 331]
   [2 -1 1 1 315]
   [2 -2 0 -1 302]
   [0 0 1 3 -283]
   [2 1 1 -1 -229]
   [1 1 0 -1 223]
   [1 1 0 1 223]
   [0 1 -2 -1 -220]
   [2 1 -1 -1 -220]
   [1 0 1 1 -185]
   [2 -1 -2 -1 181]
   [0 1 2 1 -177]
   [4 0 -2 -1 176]
   [4 -1 -1 -1 166]
   [1 0 1 -1 -164]
   [4 0 1 -1 132]
   [1 0 -1 -1 -119]
   [4 -1 0 -1 115]
   [2 -2 0 1 107]])

(defn position
  "The Moon's geocentric `[lon lat distance]` at `mjd-tt`: ecliptic
  longitude and latitude referred to the mean equinox of date, and the
  distance between the centers in km. Add the nutation in longitude for
  the apparent longitude (`apparent`)."
  [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        {:keys [L' D M M' F E]} (arguments mjd-tt)
        A1 (deg (+ 119.75 (* 131.849 T)))
        A2 (deg (+ 53.09 (* 479264.290 T)))
        A3 (deg (+ 313.45 (* 481266.484 T)))
        e-of (fn [m] (case (long (abs m)) 0 1.0 1 E 2 (* E E)))
        [sl sr] (reduce (fn [[sl sr] [d m m' f l r]]
                          (let [a (+ (* d D) (* m M) (* m' M') (* f F))
                                e (e-of m)]
                            [(+ sl (* l e (math/sin a))) (+ sr (* r e (math/cos a)))]))
                        [(+ (* 3958 (math/sin A1)) (* 1962 (math/sin (- L' F))) (* 318 (math/sin A2))) 0.0]
                        longitude-distance)
        sb (reduce (fn [sb [d m m' f b]]
                     (+ sb (* b (e-of m) (math/sin (+ (* d D) (* m M) (* m' M') (* f F))))))
                   (+ (* -2235 (math/sin L')) (* 382 (math/sin A3)) (* 175 (math/sin (- A1 F)))
                      (* 175 (math/sin (+ A1 F))) (* 127 (math/sin (- L' M'))) (* -115 (math/sin (+ L' M'))))
                   latitude-terms)]
    [(am/wrap-2pi (+ L' (deg (* sl 1e-6)))) (deg (* sb 1e-6)) (+ 385000.56 (* sr 1e-3))]))

(defn geocentric-J2000
  "The Moon's geocentric position, km, in EME2000: `position` turned from
  the ecliptic of date into the frame the rest of `allgo.astro` uses."
  ([mjd-tt] (geocentric-J2000 mjd-tt (vsop87/ecliptic-of-date->J2000 mjd-tt)))
  ([mjd-tt to-J2000]
   (let [[l b r] (position mjd-tt)]
     (lin/mat-vec to-J2000 [(* r (math/cos b) (math/cos l))
                            (* r (math/cos b) (math/sin l))
                            (* r (math/sin b))]))))

(defn horizontal-parallax
  "The Moon's equatorial horizontal parallax at `distance` km: the Earth's
  equatorial radius as seen from it, nearly a degree."
  [distance]
  (math/asin (/ 6378.14 distance)))

(defn apparent
  "The Moon's apparent `[lon lat distance]`: `position` with the nutation
  in longitude. Aberration is negligible for the Moon -- it moves with the
  Earth -- and light time, 1.3 seconds, is within the theory's accuracy."
  [mjd-tt]
  (let [[l b r] (position mjd-tt)]
    [(am/wrap-2pi (+ l (first (frames/nutation-angles mjd-tt)))) b r]))

(defn apparent-equatorial
  "The Moon's apparent `[ra dec distance]`, geocentric, true equator of date."
  [mjd-tt]
  (let [[l b r] (apparent mjd-tt)
        [ra dec] (coord/ecliptic->equatorial [l b] (frames/true-obliquity mjd-tt))]
    [ra dec r]))

(defn mean-node
  "Longitude of the mean ascending node of the lunar orbit, which
  regresses around the ecliptic every 18.6 years."
  [mjd-tt]
  (am/wrap-2pi (deg (poly/horner-ascending [125.0445479 -1934.1362891 0.0020754 (/ 1.0 467441) (/ -1.0 60616000)] (time/centuries-J2000 mjd-tt)))))

(defn true-node
  "Longitude of the true ascending node: the mean node with its periodic
  swings, the largest 1.5 degrees at twice the Sun-node angle."
  [mjd-tt]
  (let [{:keys [D M M' F]} (arguments mjd-tt)]
    (am/wrap-2pi (+ (mean-node mjd-tt)
                    (deg (+ (* -1.4979 (math/sin (* 2 (- D F)))) (* -0.1500 (math/sin M))
                            (* -0.1226 (math/sin (* 2 D))) (* 0.1176 (math/sin (* 2 F)))
                            (* -0.0801 (math/sin (* 2 (- M' F))))))))))

(defn mean-perigee
  "Longitude of the mean perigee of the lunar orbit, which advances around
  the ecliptic every 8.85 years."
  [mjd-tt]
  (am/wrap-2pi (deg (poly/horner-ascending [83.3532465 4069.0137287 -0.0103200 (/ -1.0 80053) (/ 1.0 18999000)] (time/centuries-J2000 mjd-tt)))))

;; ------------------------------------------------ illuminated fraction (48)

(defn elongation
  "Geocentric elongation of the Moon from the Sun, from both bodies' `[ra
  dec]`."
  [moon sun]
  (coord/separation moon sun))

(defn phase-angle
  "The Sun-Moon-Earth angle, from the geocentric elongation `psi` and the
  distances of the Moon and the Sun from the Earth (same unit)."
  [psi moon-distance sun-distance]
  (math/atan2 (* sun-distance (math/sin psi)) (- moon-distance (* sun-distance (math/cos psi)))))

(defn phase-angle-approximate
  "The phase angle from the four arguments alone -- 180 degrees less the
  elongation, corrected by the six largest inequalities; a few hundredths
  of a degree."
  [mjd-tt]
  (let [{:keys [D M M']} (arguments mjd-tt)]
    (+ (- math/PI (am/wrap-2pi D))
       (deg (+ (* -6.289 (math/sin M')) (* 2.100 (math/sin M)) (* -1.274 (math/sin (- (* 2 D) M')))
               (* -0.658 (math/sin (* 2 D))) (* -0.214 (math/sin (* 2 M'))) (* -0.110 (math/sin D)))))))

(defn illuminated-fraction
  "Fraction of the disk lit at phase angle `i`."
  [i]
  (* 0.5 (+ 1.0 (math/cos i))))

(defn phase
  "The Moon's phase at `mjd-tt`: `{:i :k :chi}` -- the phase angle, the
  illuminated fraction, and the position angle of the midpoint of the
  bright limb, measured east from north, which points toward the Sun."
  [mjd-tt]
  (let [[ra dec r] (apparent-equatorial mjd-tt)
        [ra0 dec0 R] (solar/apparent-equatorial mjd-tt)
        psi (elongation [ra dec] [ra0 dec0])
        i (phase-angle psi r (* R c/AU))]
    {:i i :k (illuminated-fraction i)
     :chi (coord/position-angle [ra dec] [ra0 dec0])}))

;; ------------------------------------------- physical observations (53)

(def ^:private inclination
  "Inclination of the mean lunar equator to the ecliptic, 1.54242 degrees
  -- the Moon's axis is nearly perpendicular to its orbit's reference
  plane, which is Cassini's second law."
  (deg 1.54242))

(defn- physical-libration-terms [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        {:keys [D M M' F E]} (arguments mjd-tt)
        om (mean-node mjd-tt)
        K1 (deg (+ 119.75 (* 131.849 T)))
        K2 (deg (+ 72.56 (* 20.186 T)))
        s math/sin cs math/cos]
    {:rho (deg (+ (* -0.02752 (cs M')) (* -0.02245 (s F)) (* 0.00684 (cs (- M' (* 2 F))))
                  (* -0.00293 (cs (* 2 F))) (* -0.00085 (cs (* 2 (- F D))))
                  (* -0.00054 (cs (- M' (* 2 D)))) (* -0.00020 (s (+ M' F)))
                  (* -0.00020 (cs (+ M' (* 2 F)))) (* -0.00020 (cs (- M' F)))
                  (* 0.00014 (cs (+ M' (* 2 (- F D)))))))
     :sigma (deg (+ (* -0.02816 (s M')) (* 0.02244 (cs F)) (* -0.00682 (s (- M' (* 2 F))))
                    (* -0.00279 (s (* 2 F))) (* -0.00083 (s (* 2 (- F D))))
                    (* 0.00069 (s (- M' (* 2 D)))) (* 0.00040 (cs (+ M' F)))
                    (* -0.00025 (s (* 2 M'))) (* -0.00023 (s (+ M' (* 2 F))))
                    (* 0.00020 (cs (- M' F))) (* 0.00019 (s (- M' F)))
                    (* 0.00013 (s (+ M' (* 2 (- F D))))) (* -0.00010 (cs (- M' (* 3 F))))))
     :tau (deg (+ (* 0.02520 E (s M)) (* 0.00473 (s (* 2 (- M' F)))) (* -0.00467 (s M'))
                  (* 0.00396 (s K1)) (* 0.00276 (s (* 2 (- M' D)))) (* 0.00196 (s om))
                  (* -0.00183 (cs (- M' F))) (* 0.00115 (s (- M' (* 2 D))))
                  (* -0.00096 (s (- M' D))) (* 0.00046 (s (* 2 (- F D))))
                  (* -0.00039 (s (- M' F))) (* -0.00032 (s (- M' M D)))
                  (* 0.00027 (s (- (* 2 (- M' D)) M))) (* 0.00023 (s K2))
                  (* -0.00014 (s (* 2 D))) (* 0.00014 (cs (* 2 (- M' F))))
                  (* -0.00012 (s (- M' (* 2 F)))) (* -0.00012 (s (* 2 M')))
                  (* 0.00011 (s (* 2 (- M' M D))))))
     :F F :om om}))

(defn- librate
  "Optical plus physical libration `[l b]` of a direction at ecliptic
  `[lon lat]` (without nutation)."
  [[lon lat] {:keys [rho sigma tau F om]}]
  (let [W  (- lon om)
        sI (math/sin inclination) cI (math/cos inclination)
        A  (math/atan2 (- (* (math/sin W) (math/cos lat) cI) (* (math/sin lat) sI))
                       (* (math/cos W) (math/cos lat)))
        l' (am/wrap-2pi (- A F))
        b' (math/asin (- (* (- (math/sin W)) (math/cos lat) sI) (* (math/sin lat) cI)))
        l'' (+ (- tau) (* (+ (* rho (math/cos A)) (* sigma (math/sin A))) (math/tan b')))
        b'' (- (* sigma (math/cos A)) (* rho (math/sin A)))]
    [(am/wrap-angle (+ l' l'')) (+ b' b'') l' b']))

(defn libration
  "How the Moon turns its face to the Earth at `mjd-tt`:

    :l :b   total libration in selenographic longitude and latitude -- the
            point of the Moon's surface at the center of the disk; up to
            8 degrees east-west and 7 north-south, which is how 59 percent
            of the surface can be seen from the Earth over time
    :optical  `[l b]` the part from the geometry of the orbit
    :p      position angle of the Moon's axis of rotation."
  [mjd-tt]
  (let [[lon lat] (position mjd-tt)
        {:keys [rho sigma om] :as terms} (physical-libration-terms mjd-tt)
        [l b l' b'] (librate [lon lat] terms)
        [dpsi deps] (frames/nutation-angles mjd-tt)
        eps (+ (frames/mean-obliquity mjd-tt) deps)
        V  (+ om dpsi (/ sigma (math/sin inclination)))
        X  (* (math/sin (+ inclination rho)) (math/sin V))
        Y  (- (* (math/sin (+ inclination rho)) (math/cos V) (math/cos eps))
              (* (math/cos (+ inclination rho)) (math/sin eps)))
        w  (math/atan2 X Y)
        [ra] (coord/ecliptic->equatorial [(+ lon dpsi) lat] eps)
        P  (math/asin (/ (* (math/hypot X Y) (math/cos (- ra w))) (math/cos b)))]
    {:l l :b b :optical [l' b'] :p (am/wrap-2pi P)}))

(defn selenographic-sun
  "`[l0 b0]`, the selenographic longitude and latitude of the subsolar
  point: where on the Moon the Sun is overhead. The colongitude, 90
  degrees less l0, is the longitude of the morning terminator."
  [mjd-tt]
  (let [[lon lat dist] (position mjd-tt)
        [lam0 _ R] (solar/apparent mjd-tt)
        dR (/ dist (* R c/AU))
        lamH (+ lam0 math/PI (* dR (math/cos lat) (math/sin (- lam0 lon))))
        betH (* dR lat)
        [l0 b0] (librate [lamH betH] (physical-libration-terms mjd-tt))]
    [l0 b0]))

(defn sun-altitude
  "Altitude of the Sun above the horizon at selenographic longitude `eta`
  and latitude `theta` on the Moon, when the subsolar point is `[l0 b0]`."
  [eta theta [l0 b0]]
  (let [c0 (- (/ math/PI 2) l0)]
    (math/asin (+ (* (math/sin b0) (math/sin theta))
                  (* (math/cos b0) (math/cos theta) (math/sin (+ c0 eta)))))))

(defn- sun-correction [eta theta mjd-tt]
  (/ (/ (sun-altitude eta theta (selenographic-sun mjd-tt)) c/degrees)
     (* 12.19075 (math/cos theta))))

(defn sunrise
  "MJD of sunrise at selenographic `[eta theta]` on the Moon nearest
  `mjd-tt`: the Sun climbs about 12.19 degrees a day there, so two
  corrections from its altitude suffice."
  [eta theta mjd-tt]
  (let [t (- mjd-tt (sun-correction eta theta mjd-tt))]
    (- t (sun-correction eta theta t))))

(defn sunset
  "MJD of sunset at selenographic `[eta theta]` on the Moon nearest
  `mjd-tt`."
  [eta theta mjd-tt]
  (let [t (+ mjd-tt (sun-correction eta theta mjd-tt))]
    (+ t (sun-correction eta theta t))))
