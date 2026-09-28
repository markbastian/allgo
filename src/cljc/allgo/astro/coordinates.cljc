(ns allgo.astro.coordinates
  "Positions on the celestial sphere: the changes between its coordinate
  systems, the angles between points on it, refraction and parallax
  (Meeus, *Astronomical Algorithms*, chapters 13, 14, 16 to 20 and 40).

  Four systems cover nearly everything. Equatorial coordinates -- right
  ascension and declination -- are fixed to the Earth's axis and are what
  a telescope mount and a star catalog use. Ecliptic coordinates are fixed
  to the Earth's orbit and are what theories of the Sun, Moon and planets
  come out in. Horizontal coordinates -- azimuth and altitude -- belong to
  one observer at one instant. Galactic coordinates are fixed to the Milky
  Way. The obliquity turns the first into the second, the local sidereal
  time and latitude the first into the third.

  Positions are pairs of radians: `[ra dec]`, `[lon lat]`, `[az alt]`.
  Two conventions differ from Meeus's and follow the rest of this package
  instead: geographic longitude is positive east, as the IAU has had it
  since 1982, and azimuth is measured from north through east, as a
  surveyor or a satellite tracker measures it. Meeus counts both the other
  way -- longitude west, azimuth from the south -- so his local hour angle
  θ0 - L - α appears here as θ0 + L - α and his azimuths are 180 degrees
  from these."
  (:require [allgo.astro.constants :as c]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.interpolation :as interp]
            [clojure.math :as math]))

(defn- hav [x] (* 0.5 (- 1.0 (math/cos x))))

;; ------------------------------------------------------- transformations (13)

(defn ecliptic->equatorial
  "`[ra dec]` of ecliptic `[lon lat]`, for obliquity `eps`."
  [[lon lat] eps]
  (let [se (math/sin eps) ce (math/cos eps)
        sl (math/sin lon)
        sb (math/sin lat) cb (math/cos lat)]
    [(am/wrap-2pi (math/atan2 (- (* sl ce) (* (/ sb cb) se)) (math/cos lon)))
     (math/asin (+ (* sb ce) (* cb se sl)))]))

(defn equatorial->ecliptic
  "Ecliptic `[lon lat]` of `[ra dec]`, for obliquity `eps`."
  [[ra dec] eps]
  (let [se (math/sin eps) ce (math/cos eps)
        sa (math/sin ra)
        sd (math/sin dec) cd (math/cos dec)]
    [(am/wrap-2pi (math/atan2 (+ (* sa ce) (* (/ sd cd) se)) (math/cos ra)))
     (math/asin (- (* sd ce) (* cd se sa)))]))

(defn hour-angle
  "Local hour angle of right ascension `ra` for an observer at east
  longitude `lon` when the sidereal time at Greenwich is `theta0`: how far
  the object has turned west past the meridian."
  [ra lon theta0]
  (am/wrap-2pi (- (+ theta0 lon) ra)))

(defn equatorial->horizontal
  "`[az alt]` of `[ra dec]` for an observer at latitude `lat` and east
  longitude `lon`, when the Greenwich sidereal time is `theta0`.

  Use apparent sidereal time with apparent positions, mean with mean; and
  the altitude is geometric -- add `refraction` for where the object will
  be seen."
  [[ra dec] lat lon theta0]
  (let [H  (hour-angle ra lon theta0)
        sp (math/sin lat) cp (math/cos lat)
        sd (math/sin dec) cd (math/cos dec)
        cH (math/cos H)]
    [(am/wrap-2pi (+ math/PI (math/atan2 (math/sin H) (- (* cH sp) (* (/ sd cd) cp)))))
     (math/asin (+ (* sp sd) (* cp cd cH)))]))

(defn horizontal->equatorial
  "`[ra dec]` of `[az alt]`, the inverse of `equatorial->horizontal`."
  [[az alt] lat lon theta0]
  (let [A  (- az math/PI)                ; Meeus's azimuth, from the south
        sp (math/sin lat) cp (math/cos lat)
        sh (math/sin alt) ch (math/cos alt)
        cA (math/cos A)
        H  (math/atan2 (math/sin A) (+ (* cA sp) (* (/ sh ch) cp)))]
    [(am/wrap-2pi (- (+ theta0 lon) H))
     (math/asin (- (* sp sh) (* cp ch cA)))]))

(def galactic-north-B1950
  "The north galactic pole in B1950 equatorial coordinates, `[ra dec]`: 12h
  49m, +27.4 degrees, by the IAU definition of 1959."
  [(* 192.25 c/degrees) (* 27.4 c/degrees)])

(def galactic-node-B1950
  "Galactic longitude of the ascending node of the galactic plane on the
  B1950 equator, 33 degrees."
  (* 33.0 c/degrees))

(defn equatorial->galactic
  "Galactic `[lon lat]` of B1950 `[ra dec]`. Precess a J2000 position to
  B1950 first (`allgo.astro.precession`); the galactic system is defined
  on the older frame."
  [[ra dec]]
  (let [[a0 d0] galactic-north-B1950
        sg (math/sin d0) cg (math/cos d0)
        sd (math/sin dec) cd (math/cos dec)
        da (- a0 ra)
        x  (math/atan2 (math/sin da) (- (* (math/cos da) sg) (* (/ sd cd) cg)))]
    [(am/wrap-2pi (- (+ galactic-node-B1950 (* 1.5 math/PI)) x))
     (math/asin (+ (* sd sg) (* cd cg (math/cos da))))]))

(defn galactic->equatorial
  "B1950 `[ra dec]` of galactic `[lon lat]`."
  [[lon lat]]
  (let [[a0 d0] galactic-north-B1950
        sg (math/sin d0) cg (math/cos d0)
        sb (math/sin lat) cb (math/cos lat)
        dl (- lon galactic-node-B1950 (* 0.5 math/PI))
        y  (math/atan2 (math/sin dl) (- (* (math/cos dl) sg) (* (/ sb cb) cg)))]
    [(am/wrap-2pi (- (+ y a0) math/PI))
     (math/asin (+ (* sb sg) (* cb cg (math/cos dl))))]))

;; ------------------------------------------- parallactic angle and more (14)

(defn parallactic-angle
  "The angle at the object between the directions to the zenith and to the
  celestial pole, for hour angle `H`. It is how far a field of view
  on an altazimuth mount is turned from north-up, and so how fast a camera
  on one must be derotated."
  [lat dec H]
  (math/atan2 (math/sin H)
              (- (* (math/tan lat) (math/cos dec)) (* (math/sin dec) (math/cos H)))))

(defn parallactic-angle-at-horizon
  "The parallactic angle as the object rises or sets."
  [lat dec]
  (math/acos (/ (math/sin lat) (math/cos dec))))

(defn ecliptic-at-horizon
  "Where the ecliptic meets the horizon, for obliquity `eps`, latitude
  `lat` and local sidereal time `theta`: `[lon1 lon2 angle]`, the two
  ecliptic longitudes on the horizon -- 180 degrees apart -- and the angle
  the ecliptic makes with it. That angle is why the young Moon lies on its
  back in spring evenings and stands upright in autumn ones."
  [eps lat theta]
  (let [se (math/sin eps) ce (math/cos eps)
        sp (math/sin lat) cp (math/cos lat)
        st (math/sin theta) ct (math/cos theta)
        l  (math/atan2 (- ct) (+ (* se (/ sp cp)) (* ce st)))
        l  (if (neg? l) (+ l math/PI) l)]
    [l (+ l math/PI) (math/acos (- (* ce sp) (* se cp st)))]))

(defn ecliptic-equator-angle
  "The angle between the north pole of the ecliptic and of the equator, as
  seen from the point of the ecliptic at longitude `lon` -- the position
  angle of the ecliptic's pole there."
  [lon eps]
  (math/atan (* (- (math/cos lon)) (math/tan eps))))

(defn diurnal-path-at-horizon
  "The angle between a body's daily path and the horizon as it rises or
  sets. Small in high latitudes, which is why twilight lingers there."
  [dec lat]
  (let [tp (math/tan lat)
        b  (* (math/tan dec) tp)]
    (math/atan (/ (* (math/sqrt (- 1.0 (* b b))) (math/cos dec)) tp))))

;; ------------------------------------------------------------ refraction (16)
;;
;; The atmosphere bends light toward the vertical, so everything is seen
;; higher than it is -- by 35 arcminutes at the horizon, more than the
;; Sun's diameter, which means the Sun is already geometrically below the
;; horizon at the moment it is seen to touch it.

(defn refraction-from-apparent
  "Refraction in radians, to subtract from an apparent altitude `h` to get
  the true one. Bennett's formula, good to 0.07 arcminute down to the
  horizon, with his correction when `precise?`."
  ([h] (refraction-from-apparent h false))
  ([h precise?]
   (let [hd (/ h c/degrees)
         R  (/ 1.0 (math/tan (* c/degrees (+ hd (/ 7.31 (+ hd 4.4))))))
         R  (if precise? (- R (* 0.06 (math/sin (* c/degrees (+ 13.0 (* 14.7 R)))))) R)]
     (* R (/ c/degrees 60.0)))))

(defn refraction
  "Refraction in radians, to add to a true altitude `h` to get the
  apparent one: Saemundsson's inversion of Bennett, consistent with it to
  0.1 arcminute."
  [h]
  (let [hd (/ h c/degrees)]
    (* (/ 1.02 (math/tan (* c/degrees (+ hd (/ 10.3 (+ hd 5.11))))))
       (/ c/degrees 60.0))))

(defn refraction-high
  "Refraction for altitudes above 15 degrees from the classical series in
  the tangent of the zenith distance, to add to a true altitude or -- when
  `apparent?` -- to subtract from an apparent one. Good to an arcsecond
  there and useless near the horizon, where the series diverges."
  ([h] (refraction-high h false))
  ([h apparent?]
   (let [t (math/tan (- (* 0.5 math/PI) h))
         [a b] (if apparent? [58.276 0.0824] [58.294 0.0668])]
     (* c/arcsec (- (* a t) (* b t t t))))))

(defn refraction-scale
  "Factor for any of the refractions above at pressure `p` in millibars
  and temperature `t` in Celsius; the formulae assume 1010 mb and 10 C.
  Refraction is proportional to the air's density."
  [p t]
  (* (/ p 1010.0) (/ 283.0 (+ 273.0 t))))

;; ---------------------------------------------------- angular separation (17)

(defn separation
  "Angular distance between two points on the sphere, `[ra dec]` or
  `[lon lat]` alike.

  The textbook cos d = sin d1 sin d2 + cos d1 cos d2 cos(a1 - a2) loses all
  its precision for close pairs, where the cosine is flat -- a double star
  a few arcseconds apart is not measurable by it. This is the form that
  holds everywhere (Meeus credits Thierry Pauwels), an arctangent of the
  sine and the cosine, each computed without cancellation."
  [[a1 d1] [a2 d2]]
  (let [sd1 (math/sin d1) cd1 (math/cos d1)
        sd2 (math/sin d2) cd2 (math/cos d2)
        cda (math/cos (- a2 a1))
        x (- (* cd1 sd2) (* sd1 cd2 cda))
        y (* cd2 (math/sin (- a2 a1)))
        z (+ (* sd1 sd2) (* cd1 cd2 cda))]
    (math/atan2 (math/hypot x y) z)))

(defn separation-haversine
  "Angular distance by the haversine formula, which also keeps its
  precision for close pairs. Kept for comparison with `separation`."
  [[a1 d1] [a2 d2]]
  (* 2.0 (math/asin (math/sqrt (+ (hav (- d2 d1))
                                  (* (math/cos d1) (math/cos d2) (hav (- a2 a1))))))))

(defn position-angle
  "Position angle of the second point as seen from the first, measured
  from north through east."
  [[a1 d1] [a2 d2]]
  (let [da (- a2 a1)]
    (am/wrap-2pi (math/atan2 (math/sin da)
                             (- (* (math/cos d1) (math/tan d2))
                                (* (math/sin d1) (math/cos da)))))))

(defn minimum-separation
  "Least separation of two moving bodies, given each one's position at
  three equally spaced times: interpolates the separation itself. Poor
  when the minimum is small, where the separation has a corner rather than
  a smooth minimum; see `minimum-separation-rectangular`."
  [t1 t3 positions-1 positions-2]
  (second (interp/extremum (interp/table-3 t1 t3 (map separation positions-1 positions-2)))))

(defn minimum-separation-rectangular
  "Least separation of two moving bodies from three positions of each, by
  interpolating the second's offset from the first in rectangular
  coordinates -- smooth through the closest approach, so good even for an
  occultation. Returns `[t separation]`."
  [t1 t3 positions-1 positions-2]
  (let [uv (fn [[a1 d1] [a2 d2]]
             (let [sd1 (math/sin d1) cd1 (math/cos d1)
                   da  (- a2 a1)
                   tda (math/tan da) thda (math/tan (* 0.5 da))
                   K   (/ 1.0 (+ 1.0 (* sd1 sd1 tda thda)))
                   sdd (math/sin (- d2 d1))]
               [(* (- K) (- 1.0 (* (/ sd1 cd1) sdd)) cd1 tda)
                (* K (+ sdd (* sd1 cd1 tda thda)))]))
        [[u1 v1] [u2 v2] [u3 v3]] (map uv positions-1 positions-2)
        ut (interp/table-3 -1.0 1.0 [u1 u2 u3])
        vt (interp/table-3 -1.0 1.0 [v1 v2 v3])
        up0 (* 0.5 (- u3 u1)) vp0 (* 0.5 (- v3 v1))
        up1 (+ u1 u3 (* -2.0 u2)) vp1 (+ v1 v3 (* -2.0 v2))]
    (loop [n (/ (- (+ (* u2 up0) (* v2 vp0))) (+ (* up0 up0) (* vp0 vp0))) i 0]
      (let [u (interp/value-n ut n) v (interp/value-n vt n)
            up (+ up0 (* n up1)) vp (+ vp0 (* n vp1))
            dn (/ (- (+ (* u up) (* v vp))) (+ (* up up) (* vp vp)))]
        (if (or (< (abs dn) 1e-9) (>= i 20))
          [(+ (* 0.5 (+ t1 t3)) (* 0.5 (- t3 t1) n)) (math/hypot u v)]
          (recur (+ n dn) (inc i)))))))

;; ------------------------------------------------- planetary conjunctions (18)

(defn conjunction
  "Time and declination difference of a conjunction in right ascension,
  from five equally spaced positions of each body between `t1` and `t5`.
  Returns `[t (- dec2 dec1)]`, or nil if the bodies do not pass within the
  table. For a star, pass the same position five times."
  [t1 t5 positions-1 positions-2]
  (let [dra  (map (fn [[a1] [a2]] (am/wrap-angle (- a2 a1))) positions-1 positions-2)
        ddec (map (fn [[_ d1] [_ d2]] (- d2 d1)) positions-1 positions-2)]
    (when-let [t (interp/zero (interp/table-5 t1 t5 dra) true)]
      [t (interp/value (interp/table-5 t1 t5 ddec) t)])))

;; ------------------------------------------------ bodies in a straight line (19)

(defn- great-circle-residual
  "Zero when the three points lie on one great circle (Meeus 19.1)."
  [[a1 d1] [a2 d2] [a3 d3]]
  (+ (* (math/tan d1) (math/sin (- a2 a3)))
     (* (math/tan d2) (math/sin (- a3 a1)))
     (* (math/tan d3) (math/sin (- a1 a2)))))

(defn collinear-time
  "When a moving body, at `positions-3` five equally spaced times from `t1`
  to `t5`, crosses the great circle through two fixed ones -- the moment
  three planets, or a planet and two stars, appear in a straight line."
  [p1 p2 t1 t5 positions-3]
  (interp/zero (interp/table-5 t1 t5 (map #(great-circle-residual p1 p2 %) positions-3))))

(defn collinear-angle
  "How far three points are from lying on one great circle, as the angle
  at the middle one: pi when they do."
  [p1 p2 p3]
  (am/wrap-2pi (- (position-angle p2 p3) (position-angle p2 p1))))

(defn- unit-vector [[a d]]
  [(* (math/cos d) (math/cos a)) (* (math/cos d) (math/sin a)) (math/sin d)])

(defn- norm [v] (math/sqrt (reduce + (map * v v))))

(defn- dot [a b] (reduce + (map * a b)))

(defn great-circle-distance
  "Distance of point `p0` from the great circle through `p1` and `p2`."
  [p1 p2 p0]
  (let [n (v3/cross (unit-vector p1) (unit-vector p2))]
    (math/asin (/ (dot n (unit-vector p0)) (norm n)))))

(defn collinearity
  "`[psi omega]` for three points: the angle between the great circles
  through the first two and the last two, and the distance of the middle
  one from the great circle through the outer two -- Meeus's two measures
  of how nearly three bodies are aligned."
  [p1 p2 p3]
  (let [[u1 u2 u3] (map unit-vector [p1 p2 p3])
        n12 (v3/cross u1 u2) n23 (v3/cross u2 u3) n13 (v3/cross u1 u3)]
    [(math/acos (/ (dot n12 n23) (* (norm n12) (norm n23))))
     (math/asin (/ (dot u2 n13) (* (norm u2) (norm n13))))]))

;; ------------------------------------------------------- smallest circle (20)

(defn smallest-circle
  "Diameter of the smallest circle containing three points on the sky, and
  whether it is set by the two farthest apart (`true`) or passes through
  all three (`false`). Meeus's question: can the three be seen in one
  field of view?"
  [p1 p2 p3]
  (let [[a b c] (sort > [(separation-haversine p1 p2)
                         (separation-haversine p2 p3)
                         (separation-haversine p3 p1)])]
    (if (>= (* a a) (+ (* b b) (* c c)))
      [a true]
      [(/ (* 2.0 a b c)
          (math/sqrt (* (+ a b c) (- (+ a b) c) (- (+ b c) a) (- (+ a c) b))))
       false])))

;; -------------------------------------------------------------- parallax (40)

(defn horizontal-parallax
  "Equatorial horizontal parallax of a body `distance` AU away: the angle
  the Earth's equatorial radius subtends from it."
  [distance]
  (math/asin (/ (math/sin (* 8.794 c/arcsec)) distance)))

(defn topocentric
  "Topocentric `[ra dec]` of geocentric `[ra dec]` for a body `distance` AU
  away, seen from a place with parallax constants `[rho-sin rho-cos]`
  (`allgo.astro.geodesy/parallax-constants`) at east longitude `lon`, when
  the Greenwich sidereal time is `theta0`.

  Only the Moon moves much -- up to a degree -- but a nearby asteroid or
  the Sun's limb in a transit moves enough to matter."
  [[ra dec] distance [rho-sin rho-cos] lon theta0]
  (let [sp (math/sin (horizontal-parallax distance))
        H  (hour-angle ra lon theta0)
        sH (math/sin H) cH (math/cos H)
        sd (math/sin dec) cd (math/cos dec)
        da (math/atan2 (* (- rho-cos) sp sH) (- cd (* rho-cos sp cH)))]
    [(am/wrap-2pi (+ ra da))
     (math/atan2 (* (- sd (* rho-sin sp)) (math/cos da))
                 (- cd (* rho-cos sp cH)))]))

(defn topocentric-hour-angle
  "Topocentric `[hour-angle dec]` directly, from the geocentric hour angle
  `H` -- Meeus's alternative, which avoids the right ascension."
  [H dec distance [rho-sin rho-cos]]
  (let [sp (math/sin (horizontal-parallax distance))
        A  (* (math/cos dec) (math/sin H))
        B  (- (* (math/cos dec) (math/cos H)) (* rho-cos sp))
        C  (- (math/sin dec) (* rho-sin sp))]
    [(am/wrap-2pi (math/atan2 A B))
     (math/asin (/ C (math/sqrt (+ (* A A) (* B B) (* C C)))))]))

(defn topocentric-ecliptic
  "Topocentric ecliptic `[lon lat semidiameter]` of a body at geocentric
  ecliptic `[lon lat]` with semidiameter `s` and equatorial horizontal
  parallax `parallax`, for obliquity `eps` and local sidereal time `theta`
  -- the form eclipse and occultation work wants."
  [[lon lat] s [rho-sin rho-cos] eps theta parallax]
  (let [sl (math/sin lon) cl (math/cos lon)
        sb (math/sin lat) cb (math/cos lat)
        se (math/sin eps) ce (math/cos eps)
        st (math/sin theta) ct (math/cos theta)
        sp (math/sin parallax)
        N  (- (* cl cb) (* rho-cos sp ct))
        l' (am/wrap-2pi (math/atan2 (- (* sl cb) (* sp (+ (* rho-sin se) (* rho-cos ce st)))) N))
        cl' (math/cos l')
        b' (math/atan (/ (* cl' (- sb (* sp (- (* rho-sin ce) (* rho-cos se st))))) N))]
    [l' b' (math/asin (/ (* cl' (math/cos b') (math/sin s)) N))]))
