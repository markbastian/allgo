(ns allgo.astro.physical
  "How Mars, Jupiter and Saturn's rings present themselves: which of
  their latitudes faces the Earth and the Sun, which meridian is at the
  center of the disk, and which way the axis points on the sky (Meeus,
  *Astronomical Algorithms*, chapters 42, 43 and 45).

  Each is the same calculation with a different pole. The planet's
  rotation axis is fixed in space, given as the right ascension and
  declination of its north pole; the planetocentric latitude of the Earth
  is then the angle between that pole and the line of sight, and the
  position angle is the direction of the pole as projected on the sky. The
  central meridian needs the planet's rotation angle as well, at the time
  the light left it.

  Times are MJD (TT); angles radians; distances AU."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.frames :as frames]
            [allgo.astro.illumination :as illum]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.math :as am]
            [clojure.math :as math]))

(defn- geometry
  "The planet as seen from the Earth, both in FK5, light time included:
  the planet's heliocentric `[l b r]` when its light left, the Earth's
  `[l0 b0 R]`, and the geocentric rectangular `[x y z]` and distance."
  [planet mjd-tt]
  (let [[l0 b0 R] (vsop87/heliocentric :earth mjd-tt)
        [l0 b0] (vsop87/->fk5 [l0 b0] mjd-tt)
        earth [(* R (math/cos b0) (math/cos l0)) (* R (math/cos b0) (math/sin l0)) (* R (math/sin b0))]]
    (loop [delta 0.0 i 0]
      (let [tau (elliptic/light-time delta)
            [l b r] (vsop87/heliocentric planet (- mjd-tt tau))
            [l b] (vsop87/->fk5 [l b] mjd-tt)
            xyz (mapv - [(* r (math/cos b) (math/cos l)) (* r (math/cos b) (math/sin l)) (* r (math/sin b))]
                      earth)
            delta' (math/sqrt (reduce + (map * xyz xyz)))]
        (if (or (< (abs (- delta' delta)) 1e-10) (> i 10))
          {:l l :b b :r r :l0 l0 :b0 b0 :R R :xyz xyz :delta delta' :tau tau}
          (recur delta' (inc i)))))))

(defn- equatorial-of [[x y z] eps]
  (let [u (- (* y (math/cos eps)) (* z (math/sin eps)))
        v (+ (* y (math/sin eps)) (* z (math/cos eps)))]
    [(math/atan2 u x) (math/atan (/ v (math/hypot x u)))]))

(defn- position-angle-of-pole
  "Position angle of a pole at `[a0 d0]` seen at `[a d]`."
  [[a0 d0] [a d]]
  (am/wrap-2pi (math/atan2 (* (math/cos d0) (math/sin (- a0 a)))
                           (- (* (math/sin d0) (math/cos d))
                              (* (math/cos d0) (math/sin d) (math/cos (- a0 a)))))))

;; ------------------------------------------------------------------ Mars (42)

(defn mars
  "Mars's disk at `mjd-tt`:

    :de  planetocentric declination of the Earth -- how far the north pole
         is tipped toward us
    :ds  planetocentric declination of the Sun -- the Martian season
    :omega the areographic longitude of the central meridian
    :p   position angle of the north pole, east from north on the sky
    :q   position angle of the greatest defect of illumination
    :d   apparent diameter
    :k   illuminated fraction, and
    :defect the greatest defect of illumination, the width of the dark
         crescent -- Mars shows a gibbous phase near quadrature."
  [mjd-tt]
  (let [T (/ (- mjd-tt c/mjd-J2000) 36525.0)
        l-pole (+ (c/deg 352.9065) (* (c/deg 1.1733) T))
        b-pole (- (c/deg 63.2818) (* (c/deg 0.00394) T))
        {:keys [l b r l0 R xyz delta tau]} (geometry :mars mjd-tt)
        [x y z] xyz
        lam (math/atan2 y x)
        bet (math/atan (/ z (math/hypot x y)))
        sb0 (math/sin b-pole) cb0 (math/cos b-pole)
        DE (math/asin (- (* (- sb0) (math/sin bet))
                         (* cb0 (math/cos bet) (math/cos (- l-pole lam)))))
        N  (+ (c/deg 49.5581) (* (c/deg 0.7721) T))
        l' (- l (/ (c/deg 0.00697) r))
        b' (- b (/ (* (c/deg 0.000225) (math/cos (- l N))) r))
        DS (math/asin (- (* (- sb0) (math/sin b')) (* cb0 (math/cos b') (math/cos (- l-pole l')))))
        W  (+ (c/deg 11.504) (* (c/deg 350.89200025) (- (+ mjd-tt c/jd-mjd-offset) tau 2433282.5)))
        eps0 (frames/mean-obliquity mjd-tt)
        [a0 d0] (coord/ecliptic->equatorial [l-pole b-pole] eps0)
        [a d] (equatorial-of xyz eps0)
        zeta (math/atan2 (- (* (math/sin d0) (math/cos d) (math/cos (- a0 a)))
                            (* (math/sin d) (math/cos d0)))
                         (* (math/cos d) (math/sin (- a0 a))))
        ;; apparent: aberration and nutation for the pole and the planet
        [dpsi deps] (frames/nutation-angles mjd-tt)
        lam' (+ lam (/ (* (c/deg 0.005693) (math/cos (- l0 lam))) (math/cos bet)) dpsi)
        bet' (+ bet (* (c/deg 0.005693) (math/sin (- l0 lam)) (math/sin bet)))
        eps  (+ eps0 deps)
        pole (coord/ecliptic->equatorial [(+ l-pole dpsi) b-pole] eps)
        planet (coord/ecliptic->equatorial [lam' bet'] eps)
        [as ds] (coord/ecliptic->equatorial [(+ l0 math/PI) 0.0] eps)
        chi (math/atan2 (* (math/cos ds) (math/sin (- as a)))
                        (- (* (math/sin ds) (math/cos d)) (* (math/cos ds) (math/sin d) (math/cos (- as a)))))
        diameter (/ (* 9.36 c/arcsec) delta)
        k (illum/illuminated-fraction r delta R)]
    {:de DE :ds DS :omega (am/wrap-2pi (- W zeta))
     :p (position-angle-of-pole pole planet)
     :q (am/wrap-2pi (+ chi math/PI))
     :d diameter :k k :defect (* diameter (- 1.0 k))}))

;; --------------------------------------------------------------- Jupiter (43)

(defn jupiter
  "Jupiter's disk at `mjd-tt`: `:ds` and `:de`, the jovicentric
  declinations of the Sun and the Earth; `:omega1` and `:omega2`, the
  longitudes of the central meridian in System I (the equatorial belt,
  which turns in 9h50m30s) and System II (the rest of the planet, 9h55m41s)
  -- Jupiter is a fluid and does not rotate as one body -- and `:p`, the
  position angle of the north pole."
  [mjd-tt]
  (let [jd (+ mjd-tt c/jd-mjd-offset)
        d  (- jd 2433282.5)
        T1 (/ d 36525.0)
        a0 (+ (c/deg 268.0) (* (c/deg 0.1061) T1))
        d0 (- (c/deg 64.5) (* (c/deg 0.0164) T1))
        W1 (+ (c/deg 17.71) (* (c/deg 877.90003539) d))
        W2 (+ (c/deg 16.838) (* (c/deg 870.27003539) d))
        {:keys [l b r l0 R xyz delta]} (geometry :jupiter mjd-tt)
        eps0 (frames/mean-obliquity mjd-tt)
        [as ds] (coord/ecliptic->equatorial [l b] eps0)
        sd0 (math/sin d0) cd0 (math/cos d0)
        DS (math/asin (- (* (- sd0) (math/sin ds)) (* cd0 (math/cos ds) (math/cos (- a0 as)))))
        [a d] (equatorial-of xyz eps0)
        zeta (math/atan2 (- (* sd0 (math/cos d) (math/cos (- a0 a))) (* (math/sin d) cd0))
                         (* (math/cos d) (math/sin (- a0 a))))
        DE (math/asin (- (* (- sd0) (math/sin d)) (* cd0 (math/cos d) (math/cos (- a0 a)))))
        ;; phase: the central meridian is measured from the illuminated
        ;; disk's center, displaced from the geometric one
        C  (cond-> (/ (- (+ (* 2.0 r delta) (* R R)) (* r r) (* delta delta)) (* 4.0 r delta))
             (neg? (math/sin (- l l0))) -)
        w1 (am/wrap-2pi (+ W1 (- zeta) (* (c/deg -5.07033) delta) C))
        w2 (am/wrap-2pi (+ W2 (- zeta) (* (c/deg -5.02626) delta) C))
        [dpsi deps] (frames/nutation-angles mjd-tt)
        eps (+ eps0 deps)
        se (math/sin eps) ce (math/cos eps)
        sa (math/sin a) ca (math/cos a)
        sl0 (math/sin l0) cl0 (math/cos l0)
        a  (+ a (/ (* (c/deg 0.005693) (+ (* ca cl0 ce) (* sa sl0))) (math/cos d)))
        d  (+ d (* (c/deg 0.005693) (+ (* cl0 ce (- (* (/ se ce) (math/cos d)) (* sa (math/sin d))))
                                       (* ca (math/sin d) sl0))))
        nut (fn [[a d]]
              (let [td (math/tan d)]
                [(+ a (- (* (+ ce (* se (math/sin a) td)) dpsi) (* (math/cos a) td deps)))
                 (+ d (* se (math/cos a) dpsi) (* (math/sin a) deps))]))]
    {:ds DS :de DE :omega1 w1 :omega2 w2
     :p (position-angle-of-pole (nut [a0 d0]) (nut [a d]))}))

(defn jupiter-low-precision
  "Jupiter's `:ds :de :omega1 :omega2` from mean orbits alone (Meeus 43,
  second method) -- a tenth of a degree, no ephemeris needed."
  [mjd-tt]
  (let [d  (- mjd-tt c/mjd-J2000)
        V  (+ (c/deg 172.74) (* (c/deg 0.00111588) d))
        M  (+ (c/deg 357.529) (* (c/deg 0.9856003) d))
        sV (math/sin V)
        N  (+ (c/deg 20.02) (* (c/deg 0.0830853) d) (* (c/deg 0.329) sV))
        J  (- (+ (c/deg 66.115) (* (c/deg 0.9025179) d)) (* (c/deg 0.329) sV))
        A  (+ (* (c/deg 1.915) (math/sin M)) (* (c/deg 0.020) (math/sin (* 2 M))))
        B  (+ (* (c/deg 5.555) (math/sin N)) (* (c/deg 0.168) (math/sin (* 2 N))))
        K  (- (+ J A) B)
        R  (- 1.00014 (* 0.01671 (math/cos M)) (* 0.00014 (math/cos (* 2 M))))
        r  (- 5.20872 (* 0.25208 (math/cos N)) (* 0.00611 (math/cos (* 2 N))))
        delta (math/sqrt (- (+ (* r r) (* R R)) (* 2 r R (math/cos K))))
        psi (math/asin (* (/ R delta) (math/sin K)))
        dd  (- d (/ delta 173.0))
        C   (let [s (math/sin (* 0.5 psi))] (if (pos? (math/sin K)) (- (* s s)) (* s s)))
        lam (+ (c/deg 34.35) (* (c/deg 0.083091) d) (* (c/deg 0.329) sV) B)
        DS  (* (c/deg 3.12) (math/sin (+ lam (c/deg 42.8))))]
    {:ds DS
     :de (- DS (* (c/deg 2.22) (math/sin psi) (math/cos (+ lam (c/deg 22))))
            (* (c/deg 1.3) (/ (- r delta) delta) (math/sin (- lam (c/deg 100.5)))))
     :omega1 (am/wrap-2pi (+ (c/deg 210.98) (* (c/deg 877.8169088) dd) psi (- B) C))
     :omega2 (am/wrap-2pi (+ (c/deg 187.23) (* (c/deg 870.1869088) dd) psi (- B) C))}))

;; ------------------------------------------------------- Saturn's ring (45)

(def ring-edges
  "Radii of the ring's edges as fractions of the outer edge of ring A:
  the inner edge of A (the Cassini division's outer side), the outer and
  inner edges of B, and the inner edge of the dusky ring C."
  {:inner-a 0.8801 :outer-b 0.8599 :inner-b 0.6650 :inner-c 0.5486})

(defn saturn-ring
  "Saturn's ring at `mjd-tt`:

    :b      Saturnicentric latitude of the Earth, the tilt of the ring to
            the line of sight; zero at a ring-plane crossing, when the
            rings all but vanish, every fifteen years
    :b'     the same for the Sun, which lights the side it faces
    :du     difference between the Saturnicentric longitudes of the Sun
            and the Earth, for the magnitude
    :p      position angle of the northern semiminor axis of the ring
    :a :b-axis  the outer edge of ring A's apparent semi-axes."
  [mjd-tt]
  (let [T (/ (- mjd-tt c/mjd-J2000) 36525.0)
        i  (c/deg (+ 28.075216 (* -0.012998 T) (* 0.000004 T T)))
        om (c/deg (+ 169.50847 (* 1.394681 T) (* 0.000412 T T)))
        {:keys [l b r l0 xyz delta]} (geometry :saturn mjd-tt)
        [x y z] xyz
        lam (math/atan2 y x)
        bet (math/atan (/ z (math/hypot x y)))
        si (math/sin i) ci (math/cos i)
        sB (- (* si (math/cos bet) (math/sin (- lam om))) (* ci (math/sin bet)))
        N  (c/deg (+ 113.6655 (* 0.8771 T)))
        l' (- l (/ (c/deg 0.01759) r))
        b' (- b (/ (* (c/deg 0.000764) (math/cos (- l N))) r))
        U1 (math/atan2 (+ (* si (math/sin b')) (* ci (math/cos b') (math/sin (- l' om))))
                       (* (math/cos b') (math/cos (- l' om))))
        U2 (math/atan2 (+ (* si (math/sin bet)) (* ci (math/cos bet) (math/sin (- lam om))))
                       (* (math/cos bet) (math/cos (- lam om))))
        a-edge (/ (* 375.35 c/arcsec) delta)
        [dpsi deps] (frames/nutation-angles mjd-tt)
        eps (+ (frames/mean-obliquity mjd-tt) deps)
        lam' (+ lam (/ (* (c/deg 0.005693) (math/cos (- l0 lam))) (math/cos bet)) dpsi)
        bet' (+ bet (* (c/deg 0.005693) (math/sin (- l0 lam)) (math/sin bet)))
        pole (coord/ecliptic->equatorial [(+ (- om (/ math/PI 2)) dpsi) (- (/ math/PI 2) i)] eps)
        planet (coord/ecliptic->equatorial [lam' bet'] eps)]
    {:b (math/asin sB)
     :b' (math/asin (- (* si (math/cos b') (math/sin (- l' om))) (* ci (math/sin b'))))
     :du (abs (am/wrap-angle (- U1 U2)))
     :p (am/wrap-angle (position-angle-of-pole pole planet))
     :a a-edge :b-axis (* a-edge (abs sB))}))
