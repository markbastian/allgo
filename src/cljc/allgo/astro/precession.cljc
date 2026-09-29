(ns allgo.astro.precession
  "Precession between any two epochs, of equatorial and ecliptic
  coordinates and of orbital elements, with proper motion (Meeus,
  *Astronomical Algorithms*, chapters 21 and 24).

  `allgo.astro.frames/precession` carries J2000 to a date, which is what
  an orbit integrated in EME2000 needs. A star catalog or an old
  observation asks the more general question -- from 1950 to 2028, say --
  and the answer here composes that same rotation: out of the first epoch
  to J2000 and on to the second. Meeus writes the general-epoch angles out
  as polynomials in both epochs instead; the two agree to well under a
  milliarcsecond across the centuries either is fitted for.

  Epochs are MJD (TT); `allgo.astro.time/julian-epoch->mjd` converts the
  2000.0-style epochs the literature quotes. Proper motions are radians per
  Julian year."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.coordinates :as coord]
            [allgo.astro.frames :as frames]
            [allgo.astro.time :as time]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn- angles [[x y z]]
  [(am/wrap-2pi (math/atan2 y x)) (math/atan2 z (math/hypot x y))])

;; ---------------------------------------------------------------- equatorial

(defn equatorial-matrix
  "Rotation from the mean equator and equinox of `from` to those of `to`."
  [from to]
  (lin/mat-mul (frames/precession to) (lin/transpose (frames/precession from))))

(defn equatorial
  "`[ra dec]` referred to the mean equator and equinox of `from`, carried
  to those of `to`. With proper motions `[pm-ra pm-dec]`, the star is moved
  along the sky first, at the rates of the first epoch."
  ([pos from to] (equatorial pos from to [0.0 0.0]))
  ([[ra dec] from to [pm-ra pm-dec]]
   (let [years (/ (- to from) 365.25)]
     (angles (lin/mat-vec (equatorial-matrix from to)
                          (coord/unit-vector [(+ ra (* years pm-ra)) (+ dec (* years pm-dec))]))))))

(defn annual-precession
  "`[dra ddec]`, the annual precession of `[ra dec]` in radians a year
  from the rates m and n, with `from` and `to` bracketing the interval.
  Good to a few arcseconds over a few decades away from the poles, and
  the way precession was done by hand for two centuries."
  [[ra dec] from to]
  (let [T (/ (- (time/mjd->julian-epoch to) (time/mjd->julian-epoch from)) 200.0)
        T (+ T (/ (- (time/mjd->julian-epoch from) 2000.0) 100.0))
        m  (* 15.0 c/arcsec (+ 3.07496 (* 0.00186 T)))
        na (* 15.0 c/arcsec (- 1.33621 (* 0.00057 T)))
        nd (* c/arcsec (- 20.0431 (* 0.0085 T)))]
    [(+ m (* na (math/sin ra) (math/tan dec)))
     (* nd (math/cos ra))]))

(defn equatorial-approximate
  "`equatorial` by `annual-precession` and proper motion, straight-line."
  [[ra dec :as pos] from to [pm-ra pm-dec]]
  (let [[dra ddec] (annual-precession pos from to)
        years (/ (- to from) 365.25)]
    [(am/wrap-2pi (+ ra (* years (+ dra pm-ra))))
     (+ dec (* years (+ ddec pm-dec)))]))

(defn proper-motion-3d
  "`[ra dec]` at `to` of a star at `[ra dec]` at `from`, `distance` away
  and receding at `radial-velocity` (same length unit per Julian year),
  with proper motions `[pm-ra pm-dec]` -- moved in a straight line in
  space rather than along the sky. For a nearby fast star over millennia
  the two differ; Barnard's star is the classic case."
  [[ra dec] from to distance radial-velocity [pm-ra pm-dec]]
  (let [[x y z] (lin/scale (coord/unit-vector [ra dec]) distance)
        mrr (/ radial-velocity distance)
        zm  (* z pm-dec)
        mx  (- (* x mrr) (* zm (math/cos ra)) (* y pm-ra))
        my  (+ (- (* y mrr) (* zm (math/sin ra))) (* x pm-ra))
        mz  (+ (* z mrr) (* distance pm-dec (math/cos dec)))
        t   (/ (- to from) 365.25)]
    (angles [(+ x (* t mx)) (+ y (* t my)) (+ z (* t mz))])))

;; ------------------------------------------------------------------ ecliptic

(defn- ecliptic-angles
  "Meeus's eta, pi and p (21.5), in radians: the tilt between the two
  ecliptics, the longitude of its hinge, and the general precession in
  longitude between `from` and `to`."
  [from to]
  (let [T (time/centuries-J2000 from)
        t (/ (- to from) 36525.0)
        s c/arcsec
        eta (* t (+ (* s (+ 47.0029 (* -0.06603 T) (* 0.000598 T T)))
                    (* t s (+ -0.03302 (* 0.000598 T)))
                    (* t t s 0.000060)))
        pi' (+ (* c/degrees 174.876384) (* s (+ (* 3289.4789 T) (* 0.60622 T T)))
               (* t s (- (+ 869.8089 (* 0.50491 T))))
               (* t t s 0.03536))
        p   (* t (+ (* s (+ 5029.0966 (* 2.22226 T) (* -0.000042 T T)))
                    (* t s (+ 1.11113 (* -0.000042 T)))
                    (* t t s -0.000006)))]
    [eta pi' p]))

(defn ecliptic
  "Ecliptic `[lon lat]` referred to the mean ecliptic and equinox of
  `from`, carried to those of `to` -- the equinox slides along, and the
  ecliptic itself tips by some 47 arcseconds a century."
  [[lon lat] from to]
  (let [[eta pi' p] (ecliptic-angles from to)
        se (math/sin eta) ce (math/cos eta)
        sb (math/sin lat) cb (math/cos lat)
        d  (- pi' lon)
        A  (- (* ce cb (math/sin d)) (* se sb))
        B  (* cb (math/cos d))
        C  (+ (* ce sb) (* se cb (math/sin d)))]
    [(am/wrap-2pi (- (+ p pi') (math/atan2 A B)))
     (math/asin C)]))

(defn ecliptic-with-proper-motion
  "`ecliptic`, the star first moved by equatorial proper motions
  `[pm-ra pm-dec]` converted to ecliptic ones at the first epoch."
  [[lon lat :as pos] from to [pm-ra pm-dec]]
  (let [eps (frames/mean-obliquity from)
        se (math/sin eps) ce (math/cos eps)
        [ra dec] (coord/ecliptic->equatorial pos eps)
        sa (math/sin ra) ca (math/cos ra)
        sd (math/sin dec) cd (math/cos dec)
        cb (math/cos lat)
        k  (+ (* ce cd) (* se sd sa))
        ml (/ (+ (* pm-dec se ca) (* pm-ra cd k)) (* cb cb))
        mb (/ (- (* pm-dec k) (* pm-ra se ca cd)) cb)
        years (/ (- to from) 365.25)]
    (ecliptic [(+ lon (* years ml)) (+ lat (* years mb))] from to)))

;; ---------------------------------------------------- orbital elements (24)

(defn elements
  "Inclination, argument of perihelion and longitude of the node `{:i
  :argp :raan}` of an orbit, referred to the ecliptic and equinox of
  `from`, carried to those of `to`. The shape of the orbit and the time of
  perihelion are unaffected -- only the reference plane moves."
  [{:keys [i argp raan] :as el} from to]
  (let [[eta pi' p] (ecliptic-angles from to)
        se (math/sin eta) ce (math/cos eta)
        si (math/sin i) ci (math/cos i)
        d  (- raan pi')
        sd (math/sin d) cd (math/cos d)]
    (assoc el
           :i    (math/acos (+ (* ci ce) (* si se cd)))
           :raan (am/wrap-2pi (+ pi' p (math/atan2 (* si sd) (- (* ce si cd) (* se ci)))))
           :argp (am/wrap-2pi (+ argp (math/atan2 (* (- se) sd) (- (* si ce) (* ci se cd))))))))

(defn elements-B1950->J2000
  "Elements referred to the B1950 ecliptic and equinox, to J2000, by the
  fixed rotation Meeus gives for exactly that pair."
  [{:keys [i argp raan] :as el}]
  (let [S 0.0001139788 C 0.9999999935
        W  (- raan (* 174.298782 c/degrees))
        si (math/sin i) ci (math/cos i)
        sW (math/sin W) cW (math/cos W)
        A  (* si sW)
        B  (- (* C si cW) (* S ci))]
    (assoc el
           :i    (math/asin (math/hypot A B))
           :raan (am/wrap-2pi (+ (* 174.997194 c/degrees) (math/atan2 A B)))
           :argp (am/wrap-2pi (+ argp (math/atan2 (* (- S) sW) (- (* C si) (* S ci cW))))))))

(defn elements-FK4->FK5
  "Elements referred to the B1950 equinox of the FK4 catalog, to the J2000
  equinox of FK5. Differs from `elements-B1950->J2000` by the 0.525
  arcsecond equinox correction between the two catalogs; this is the one
  to use for elements published before 1984."
  [{:keys [i argp raan] :as el}]
  (let [L  (* 5.19856209 c/degrees)
        L' (* 4.50001688 c/degrees)
        J  (* 0.00651966 c/degrees)
        W  (+ L raan)
        si (math/sin i) ci (math/cos i)
        sJ (math/sin J) cJ (math/cos J)
        sW (math/sin W) cW (math/cos W)]
    (assoc el
           :i    (math/acos (- (* ci cJ) (* si sJ cW)))
           :raan (am/wrap-2pi (- (math/atan2 (* si sW) (+ (* ci sJ) (* si cJ cW))) L'))
           :argp (am/wrap-2pi (+ argp (math/atan2 (* sJ sW) (+ (* si cJ) (* ci sJ cW))))))))
