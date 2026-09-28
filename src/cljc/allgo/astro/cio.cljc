(ns allgo.astro.cio
  "The IAU 2006/2000A reduction between the celestial (GCRS) and
  terrestrial (ITRS) frames through the celestial intermediate origin --
  the IERS Conventions (2010), chapter 5, and the modern successor to the
  FK5 chain of `allgo.astro.reduction` (Vallado chapter 3).

    GCRS --X, Y, s--> CIRS --Earth rotation angle--> TIRS --polar motion--> ITRS

  The celestial intermediate pole's coordinates X and Y carry precession
  and nutation together, as series of some 4700 terms in the fundamental
  arguments (`allgo.astro.cio-data`); s places the CIO on the pole's
  equator, and the Earth's rotation is then one angle, linear in UT1,
  with no equation of the equinoxes to add. Earth orientation parameters
  are a map: `:xp` `:yp` polar motion and `:dx` `:dy` the IERS offsets of
  the pole from the model, radians; `:lod`, seconds. Times are MJD, the
  rotation UT1 and everything else TT. States are `[r v]`, km and km/s."
  (:require [allgo.astro.cio-data :as data]
            [allgo.astro.frames :as frames]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def ^:private as->rad (/ math/PI 648000.0))
(def ^:private twopi (* 2.0 math/PI))

(defn- centuries [mjd-tt] (/ (- mjd-tt 51544.5) 36525.0))

(defn fundamental-arguments
  "The fourteen arguments of the IERS Conventions (2003 and 2010), radians,
  at `t` Julian centuries of TT from J2000: the Moon's and Sun's mean
  anomalies l and l', F, D and the Moon's node Omega (Simon et al. 1994);
  the mean longitudes of Mercury through Neptune; and the general
  precession in longitude p_A."
  [t]
  (let [arcsec (fn [& cs] (* as->rad (am/fmod (reduce (fn [acc k] (+ (* acc t) k)) 0.0 (reverse cs)) 1296000.0)))
        lin (fn [a b] (am/fmod (+ a (* b t)) twopi))]
    [(arcsec 485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470)
     (arcsec 1287104.793048 129596581.0481 -0.5532 0.000136 -0.00001149)
     (arcsec 335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417)
     (arcsec 1072260.703692 1602961601.2090 -6.3706 0.006593 -0.00003169)
     (arcsec 450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939)
     (lin 4.402608842 2608.7903141574)
     (lin 3.176146697 1021.3285546211)
     (lin 1.753470314 628.3075849991)
     (lin 6.203480913 334.0612426700)
     (lin 0.599546497 52.9690962641)
     (lin 0.874016757 21.3299104960)
     (lin 5.481293872 7.4781598567)
     (lin 5.311886287 3.8133035638)
     (* (+ 0.02438175 (* 0.00000538691 t)) t)]))

(defn- argument [mults fa]
  (reduce + (map (fn [m a] (if (zero? m) 0.0 (* m a))) mults fa)))

(defn xy
  "`[X Y]`, radians: the celestial intermediate pole in the GCRS, IAU
  2006/2000A, at TT `mjd-tt`."
  [mjd-tt]
  (let [t (centuries mjd-tt)
        pt (vec (take 6 (iterate #(* % t) 1.0)))
        fa (fundamental-arguments t)
        poly (mapv (fn [cs] (reduce + (map * cs pt))) data/xy-polynomial)
        periodic (reduce (fn [[x y] [mults terms]]
                           (let [a (argument mults fa)
                                 sc [(math/sin a) (math/cos a)]]
                             (reduce (fn [[x y] [xy s pow amp]]
                                       (let [v (* amp (sc s) (pt pow))]
                                         (if (zero? xy) [(+ x v) y] [x (+ y v)])))
                                     [x y] terms)))
                         [0.0 0.0] data/xy-terms)]
    (mapv (fn [p q] (* as->rad (+ p (* 1e-6 q)))) poly periodic)))

(defn s
  "The CIO locator s, radians, at TT `mjd-tt` for the pole at `x` `y`:
  the series for s + XY/2, less XY/2."
  [mjd-tt x y]
  (let [t (centuries mjd-tt)
        fa (fundamental-arguments t)
        w (reduce (fn [w [k mults sn cs]]
                    (let [a (argument mults fa)]
                      (update w k + (* sn (math/sin a)) (* cs (math/cos a)))))
                  data/s-polynomial data/s-terms)
        total (reduce (fn [acc k] (+ (* acc t) (w k))) 0.0 (range 5 -1 -1))]
    (- (* as->rad total) (* 0.5 x y))))

(defn earth-rotation-angle
  "The Earth rotation angle, radians, at UT1 `mjd-ut1`: 2 pi (0.7790572732640
  + 1.00273781191135448 Du), Du days from J2000 -- the day's fraction
  taken apart from the rest to keep its digits."
  [mjd-ut1]
  (let [whole (math/floor mjd-ut1)
        f (- mjd-ut1 whole)
        du (- mjd-ut1 51544.5)]
    (am/fmod (* twopi (+ f 0.5 0.7790572732640 (* 0.00273781191135448 du))) twopi)))

(defn s-prime
  "The TIO locator s', radians: -47 microarcseconds a century."
  [mjd-tt]
  (* -47e-6 as->rad (centuries mjd-tt)))

(defn c2i-matrix
  "GCRS to CIRS: the pole's X and Y (with the IERS offsets `:dx` `:dy`) and
  s turned into a rotation, Rz(-(E + s)) Ry(d) Rz(E), E and d the pole's
  azimuth and polar distance."
  ([mjd-tt] (c2i-matrix mjd-tt {}))
  ([mjd-tt {:keys [dx dy] :or {dx 0.0 dy 0.0}}]
   (let [[x0 y0] (xy mjd-tt)
         sv (s mjd-tt x0 y0)
         x (+ x0 dx) y (+ y0 dy)
         r2 (+ (* x x) (* y y))
         e (if (pos? r2) (math/atan2 y x) 0.0)
         d (math/atan (math/sqrt (/ r2 (- 1.0 r2))))]
     (frames/chain (frames/rz (- (+ e sv))) (frames/ry d) (frames/rz e)))))

(defn polar-matrix
  "TIRS to ITRS: Rx(-yp) Ry(-xp) Rz(s')."
  [mjd-tt {:keys [xp yp] :or {xp 0.0 yp 0.0}}]
  (frames/chain (frames/rx (- yp)) (frames/ry (- xp)) (frames/rz (s-prime mjd-tt))))

(defn c2t-matrix
  "GCRS to ITRS: polar motion, the Earth rotation angle and the CIP."
  ([mjd-tt mjd-ut1] (c2t-matrix mjd-tt mjd-ut1 {}))
  ([mjd-tt mjd-ut1 eop]
   (lin/mat-mul (polar-matrix mjd-tt eop)
                (lin/mat-mul (frames/rz (earth-rotation-angle mjd-ut1)) (c2i-matrix mjd-tt eop)))))

(def ^:private omega-nominal 7.292115146706979e-5)

(defn gcrs->itrs
  "A GCRS state in the ITRS; the velocity loses the Earth's spin, at the
  rate the length of day sets."
  ([st mjd-tt mjd-ut1] (gcrs->itrs st mjd-tt mjd-ut1 {}))
  ([[r v] mjd-tt mjd-ut1 {:keys [lod] :or {lod 0.0} :as eop}]
   (let [Q (c2i-matrix mjd-tt eop)
         R (frames/rz (earth-rotation-angle mjd-ut1))
         W (polar-matrix mjd-tt eop)
         w [0.0 0.0 (* omega-nominal (- 1.0 (/ lod 86400.0)))]
         rt (lin/mat-vec R (lin/mat-vec Q r))
         vt (lin/sub (lin/mat-vec R (lin/mat-vec Q v)) (v3/cross w rt))]
     [(lin/mat-vec W rt) (lin/mat-vec W vt)])))

(defn itrs->gcrs
  "An ITRS state in the GCRS: the inverse of `gcrs->itrs`."
  ([st mjd-tt mjd-ut1] (itrs->gcrs st mjd-tt mjd-ut1 {}))
  ([[r v] mjd-tt mjd-ut1 {:keys [lod] :or {lod 0.0} :as eop}]
   (let [Qt (lin/transpose (c2i-matrix mjd-tt eop))
         Rt (lin/transpose (frames/rz (earth-rotation-angle mjd-ut1)))
         Wt (lin/transpose (polar-matrix mjd-tt eop))
         w [0.0 0.0 (* omega-nominal (- 1.0 (/ lod 86400.0)))]
         rt (lin/mat-vec Wt r)
         vt (lin/add (lin/mat-vec Wt v) (v3/cross w rt))]
     [(lin/mat-vec Qt (lin/mat-vec Rt rt)) (lin/mat-vec Qt (lin/mat-vec Rt vt))])))
