(ns allgo.astro.frames
  "Precession, nutation and the celestial reference frame
  (Montenbruck & Gill 5.2).

  The Earth's axis is not fixed. The Sun and Moon pull on the equatorial
  bulge and torque the planet like a gymbal, so the pole traces a circle of
  23.4 degrees radius every 26,000 years -- precession -- with a smaller
  wobble of 9 arcseconds superposed on it at the 18.6-year period of the
  lunar node -- nutation.

  Both matter for orbits because they move the frame the equations are
  written in. Fifty arcseconds a year is 12 meters a year at geostationary
  radius, and nutation's 17 arcseconds is 3.5 kilometers of instantaneous
  offset that no amount of careful integration will recover."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.time :as time]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

;; ------------------------------------------------------------------ matrices
;; Frame rotations, in the astrodynamical sense: R(theta) * v gives the
;; components of v in a frame turned by theta, not v turned by theta.

(defn rx [t] (let [s (math/sin t) k (math/cos t)]
               [[1.0 0.0 0.0] [0.0 k s] [0.0 (- s) k]]))
(defn ry [t] (let [s (math/sin t) k (math/cos t)]
               [[k 0.0 (- s)] [0.0 1.0 0.0] [s 0.0 k]]))
(defn rz [t] (let [s (math/sin t) k (math/cos t)]
               [[k s 0.0] [(- s) k 0.0] [0.0 0.0 1.0]]))

(defn chain
  "Compose rotations left to right, so `(chain a b c)` applied to a vector
  does c first."
  [& ms]
  (reduce lin/mat-mul ms))

;; ---------------------------------------------------------------- precession

(defn precession
  "Rotation from the mean equator and equinox of J2000 to that of `mjd-tt`
  (IAU 1976).

  Three angles, all near 2000 arcseconds a century: the equinox slides along
  the ecliptic and the pole tips, and the composition of the three carries
  one frame into the other."
  [mjd-tt]
  (let [T    (time/centuries-J2000 mjd-tt)
        T2   (* T T)
        T3   (* T2 T)
        zeta (* c/arcsec (+ (* 2306.2181 T) (* 0.30188 T2) (* 0.017998 T3)))
        z    (* c/arcsec (+ (* 2306.2181 T) (* 1.09468 T2) (* 0.018203 T3)))
        th   (* c/arcsec (- (* 2004.3109 T) (* 0.42665 T2) (* 0.041833 T3)))]
    (chain (rz (- z)) (ry th) (rz (- zeta)))))

(defn mean-obliquity
  "Obliquity of the ecliptic referred to the mean equator, radians
  (IAU 1980). It is decreasing by about 47 arcseconds a century."
  [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)]
    (* c/arcsec (+ 84381.448 (* -46.8150 T) (* -0.00059 T T) (* 0.001813 T T T)))))

(defn mean-obliquity-laskar
  "Mean obliquity from Laskar's polynomial, radians (Meeus 22.3).

  The IAU expression above is a cubic fitted over a few centuries and
  drifts by some arcseconds two thousand years out. Laskar's runs to the
  tenth power of ten-thousand-year units and holds to 0.01 arcsecond over
  a thousand years and a few arcseconds over ten thousand -- useful
  wherever the sky is wanted far from the present."
  [mjd-tt]
  (let [U (/ (time/centuries-J2000 mjd-tt) 100.0)]
    (* c/arcsec
       (reduce (fn [acc k] (+ (* acc U) k)) 0.0
               [2.45 5.79 27.87 7.12 -39.05 -249.67 -51.38 1999.25 -1.55 -4680.93
                84381.448]))))

;; ------------------------------------------------------------------ nutation

(def ^:private nutation-series
  "The IAU 1980 theory of nutation, all 106 terms: the multipliers of the
  five Delaunay arguments, then the coefficients of sin for longitude and
  cos for obliquity, in units of 0.0001 arcsecond, each with its rate per
  century. Table I of Seidelmann, '1980 IAU Theory of Nutation: The Final
  Report of the IAU Working Group on Nutation', Celestial Mechanics 27
  (1982) 79, in its order and with its term numbers, and the same terms
  the IERS distributes as nut_IAU1980.dat. Meeus's table 22.A takes its
  63 largest terms from it.

  The FK5 reduction is defined with this series, so an Earth-fixed
  position computed by it agrees with ERFA's to a tenth of a millimeter. What
  it leaves out -- the free core nutation and the errors of the theory --
  the IERS publishes as dpsi and deps corrections, which
  `nutation-angles` takes."
  ;;  l  l'  F  D Om   dpsi(sin)         deps(cos)          term
  [[0  0  0  0  1  -171996.0 -174.2  92025.0   8.9]  ; 1
   [0  0  0  0  2     2062.0    0.2   -895.0   0.5]  ; 2
   [-2  0  2  0  1       46.0    0.0    -24.0   0.0]  ; 3
   [2  0 -2  0  0       11.0    0.0      0.0   0.0]  ; 4
   [-2  0  2  0  2       -3.0    0.0      1.0   0.0]  ; 5
   [1 -1  0 -1  0       -3.0    0.0      0.0   0.0]  ; 6
   [0 -2  2 -2  1       -2.0    0.0      1.0   0.0]  ; 7
   [2  0 -2  0  1        1.0    0.0      0.0   0.0]  ; 8
   [0  0  2 -2  2   -13187.0   -1.6   5736.0  -3.1]  ; 9
   [0  1  0  0  0     1426.0   -3.4     54.0  -0.1]  ; 10
   [0  1  2 -2  2     -517.0    1.2    224.0  -0.6]  ; 11
   [0 -1  2 -2  2      217.0   -0.5    -95.0   0.3]  ; 12
   [0  0  2 -2  1      129.0    0.1    -70.0   0.0]  ; 13
   [2  0  0 -2  0       48.0    0.0      1.0   0.0]  ; 14
   [0  0  2 -2  0      -22.0    0.0      0.0   0.0]  ; 15
   [0  2  0  0  0       17.0   -0.1      0.0   0.0]  ; 16
   [0  1  0  0  1      -15.0    0.0      9.0   0.0]  ; 17
   [0  2  2 -2  2      -16.0    0.1      7.0   0.0]  ; 18
   [0 -1  0  0  1      -12.0    0.0      6.0   0.0]  ; 19
   [-2  0  0  2  1       -6.0    0.0      3.0   0.0]  ; 20
   [0 -1  2 -2  1       -5.0    0.0      3.0   0.0]  ; 21
   [2  0  0 -2  1        4.0    0.0     -2.0   0.0]  ; 22
   [0  1  2 -2  1        4.0    0.0     -2.0   0.0]  ; 23
   [1  0  0 -1  0       -4.0    0.0      0.0   0.0]  ; 24
   [2  1  0 -2  0        1.0    0.0      0.0   0.0]  ; 25
   [0  0 -2  2  1        1.0    0.0      0.0   0.0]  ; 26
   [0  1 -2  2  0       -1.0    0.0      0.0   0.0]  ; 27
   [0  1  0  0  2        1.0    0.0      0.0   0.0]  ; 28
   [-1  0  0  1  1        1.0    0.0      0.0   0.0]  ; 29
   [0  1  2 -2  0       -1.0    0.0      0.0   0.0]  ; 30
   [0  0  2  0  2    -2274.0   -0.2    977.0  -0.5]  ; 31
   [1  0  0  0  0      712.0    0.1     -7.0   0.0]  ; 32
   [0  0  2  0  1     -386.0   -0.4    200.0   0.0]  ; 33
   [1  0  2  0  2     -301.0    0.0    129.0  -0.1]  ; 34
   [1  0  0 -2  0     -158.0    0.0     -1.0   0.0]  ; 35
   [-1  0  2  0  2      123.0    0.0    -53.0   0.0]  ; 36
   [0  0  0  2  0       63.0    0.0     -2.0   0.0]  ; 37
   [1  0  0  0  1       63.0    0.1    -33.0   0.0]  ; 38
   [-1  0  0  0  1      -58.0   -0.1     32.0   0.0]  ; 39
   [-1  0  2  2  2      -59.0    0.0     26.0   0.0]  ; 40
   [1  0  2  0  1      -51.0    0.0     27.0   0.0]  ; 41
   [0  0  2  2  2      -38.0    0.0     16.0   0.0]  ; 42
   [2  0  0  0  0       29.0    0.0     -1.0   0.0]  ; 43
   [1  0  2 -2  2       29.0    0.0    -12.0   0.0]  ; 44
   [2  0  2  0  2      -31.0    0.0     13.0   0.0]  ; 45
   [0  0  2  0  0       26.0    0.0     -1.0   0.0]  ; 46
   [-1  0  2  0  1       21.0    0.0    -10.0   0.0]  ; 47
   [-1  0  0  2  1       16.0    0.0     -8.0   0.0]  ; 48
   [1  0  0 -2  1      -13.0    0.0      7.0   0.0]  ; 49
   [-1  0  2  2  1      -10.0    0.0      5.0   0.0]  ; 50
   [1  1  0 -2  0       -7.0    0.0      0.0   0.0]  ; 51
   [0  1  2  0  2        7.0    0.0     -3.0   0.0]  ; 52
   [0 -1  2  0  2       -7.0    0.0      3.0   0.0]  ; 53
   [1  0  2  2  2       -8.0    0.0      3.0   0.0]  ; 54
   [1  0  0  2  0        6.0    0.0      0.0   0.0]  ; 55
   [2  0  2 -2  2        6.0    0.0     -3.0   0.0]  ; 56
   [0  0  0  2  1       -6.0    0.0      3.0   0.0]  ; 57
   [0  0  2  2  1       -7.0    0.0      3.0   0.0]  ; 58
   [1  0  2 -2  1        6.0    0.0     -3.0   0.0]  ; 59
   [0  0  0 -2  1       -5.0    0.0      3.0   0.0]  ; 60
   [1 -1  0  0  0        5.0    0.0      0.0   0.0]  ; 61
   [2  0  2  0  1       -5.0    0.0      3.0   0.0]  ; 62
   [0  1  0 -2  0       -4.0    0.0      0.0   0.0]  ; 63
   [1  0 -2  0  0        4.0    0.0      0.0   0.0]  ; 64
   [0  0  0  1  0       -4.0    0.0      0.0   0.0]  ; 65
   [1  1  0  0  0       -3.0    0.0      0.0   0.0]  ; 66
   [1  0  2  0  0        3.0    0.0      0.0   0.0]  ; 67
   [1 -1  2  0  2       -3.0    0.0      1.0   0.0]  ; 68
   [-1 -1  2  2  2       -3.0    0.0      1.0   0.0]  ; 69
   [-2  0  0  0  1       -2.0    0.0      1.0   0.0]  ; 70
   [3  0  2  0  2       -3.0    0.0      1.0   0.0]  ; 71
   [0 -1  2  2  2       -3.0    0.0      1.0   0.0]  ; 72
   [1  1  2  0  2        2.0    0.0     -1.0   0.0]  ; 73
   [-1  0  2 -2  1       -2.0    0.0      1.0   0.0]  ; 74
   [2  0  0  0  1        2.0    0.0     -1.0   0.0]  ; 75
   [1  0  0  0  2       -2.0    0.0      1.0   0.0]  ; 76
   [3  0  0  0  0        2.0    0.0      0.0   0.0]  ; 77
   [0  0  2  1  2        2.0    0.0     -1.0   0.0]  ; 78
   [-1  0  0  0  2        1.0    0.0     -1.0   0.0]  ; 79
   [1  0  0 -4  0       -1.0    0.0      0.0   0.0]  ; 80
   [-2  0  2  2  2        1.0    0.0     -1.0   0.0]  ; 81
   [-1  0  2  4  2       -2.0    0.0      1.0   0.0]  ; 82
   [2  0  0 -4  0       -1.0    0.0      0.0   0.0]  ; 83
   [1  1  2 -2  2        1.0    0.0     -1.0   0.0]  ; 84
   [1  0  2  2  1       -1.0    0.0      1.0   0.0]  ; 85
   [-2  0  2  4  2       -1.0    0.0      1.0   0.0]  ; 86
   [-1  0  4  0  2        1.0    0.0      0.0   0.0]  ; 87
   [1 -1  0 -2  0        1.0    0.0      0.0   0.0]  ; 88
   [2  0  2 -2  1        1.0    0.0     -1.0   0.0]  ; 89
   [2  0  2  2  2       -1.0    0.0      0.0   0.0]  ; 90
   [1  0  0  2  1       -1.0    0.0      0.0   0.0]  ; 91
   [0  0  4 -2  2        1.0    0.0      0.0   0.0]  ; 92
   [3  0  2 -2  2        1.0    0.0      0.0   0.0]  ; 93
   [1  0  2 -2  0       -1.0    0.0      0.0   0.0]  ; 94
   [0  1  2  0  1        1.0    0.0      0.0   0.0]  ; 95
   [-1 -1  0  2  1        1.0    0.0      0.0   0.0]  ; 96
   [0  0 -2  0  1       -1.0    0.0      0.0   0.0]  ; 97
   [0  0  2 -1  2       -1.0    0.0      0.0   0.0]  ; 98
   [0  1  0  2  0       -1.0    0.0      0.0   0.0]  ; 99
   [1  0 -2 -2  0       -1.0    0.0      0.0   0.0]  ; 100
   [0 -1  2  0  1       -1.0    0.0      0.0   0.0]  ; 101
   [1  1  0 -2  1       -1.0    0.0      0.0   0.0]  ; 102
   [1  0 -2  2  0       -1.0    0.0      0.0   0.0]  ; 103
   [2  0  0  2  0        1.0    0.0      0.0   0.0]  ; 104
   [0  0  2  4  2       -1.0    0.0      0.0   0.0]  ; 105
   [0  1  0  1  0        1.0    0.0      0.0   0.0]]) ; 106

(defn delaunay
  "The five fundamental arguments of lunar and solar theory, radians:
  the Moon's mean anomaly, the Sun's, the Moon's argument of latitude, their
  elongation, and the longitude of the Moon's ascending node.

  Nutation is almost entirely the last of these -- the node's 18.6-year
  circuit is the 17-arcsecond term that dominates everything else."
  [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        d (fn [& coeffs] (* c/arcsec (reduce (fn [acc [i k]] (+ acc (* k (math/pow T i))))
                                             0.0 (map-indexed vector coeffs))))]
    [(+ (* c/degrees 134.96298139) (d 0.0 1717915922.6330 31.310 0.064))
     (+ (* c/degrees 357.52772333) (d 0.0 129596581.2240 -0.577 -0.012))
     (+ (* c/degrees 93.27191028)  (d 0.0 1739527263.1370 -13.257 0.011))
     (+ (* c/degrees 297.85036306) (d 0.0 1602961601.3280 -6.891 0.019))
     (+ (* c/degrees 125.04452222) (d 0.0 -6962890.5390 7.455 0.008))]))

(defn nutation-angles
  "Nutation in longitude and obliquity, radians (IAU 1980, truncated)."
  [mjd-tt]
  (let [T    (time/centuries-J2000 mjd-tt)
        args (delaunay mjd-tt)]
    (reduce (fn [[dpsi deps] [ml ms mf md mo sp sp-t ce ce-t]]
              (let [a (lin/dot [ml ms mf md mo] args)]
                [(+ dpsi (* c/arcsec 1e-4 (+ sp (* sp-t T)) (math/sin a)))
                 (+ deps (* c/arcsec 1e-4 (+ ce (* ce-t T)) (math/cos a)))]))
            [0.0 0.0]
            nutation-series)))

(defn nutation-angles-low-precision
  "Nutation in longitude and obliquity from the four largest terms, good to
  half an arcsecond and a tenth of one respectively (Meeus 22). Enough
  wherever nutation is only a correction to something coarser, as it is
  for the Sun's apparent longitude to the nearest arcsecond or two."
  [mjd-tt]
  (let [T  (time/centuries-J2000 mjd-tt)
        om (* c/degrees (- 125.04452 (* 1934.136261 T)))
        L  (* c/degrees (+ 280.4665 (* 36000.7698 T)))
        L' (* c/degrees (+ 218.3165 (* 481267.8813 T)))]
    [(* c/arcsec (+ (* -17.20 (math/sin om)) (* -1.32 (math/sin (* 2.0 L)))
                    (* -0.23 (math/sin (* 2.0 L'))) (* 0.21 (math/sin (* 2.0 om)))))
     (* c/arcsec (+ (* 9.20 (math/cos om)) (* 0.57 (math/cos (* 2.0 L)))
                    (* 0.10 (math/cos (* 2.0 L'))) (* -0.09 (math/cos (* 2.0 om)))))]))

(defn true-obliquity
  "The obliquity of the true equator of date: the mean obliquity plus the
  nutation in obliquity. This is the angle to rotate apparent ecliptic
  coordinates by."
  [mjd-tt]
  (+ (mean-obliquity mjd-tt) (second (nutation-angles mjd-tt))))

(defn nutation
  "Rotation from the mean equator and equinox of date to the true one."
  [mjd-tt]
  (let [eps (mean-obliquity mjd-tt)
        [dpsi deps] (nutation-angles mjd-tt)]
    (chain (rx (- (+ eps deps))) (rz (- dpsi)) (rx eps))))

(defn equation-of-equinoxes
  "The gap between apparent and mean sidereal time, radians. It is the
  nutation in longitude projected onto the equator -- the true equinox
  wanders, so a clock keeping time by it wanders too, by up to a second."
  [mjd-tt]
  (let [[dpsi deps] (nutation-angles mjd-tt)]
    (* dpsi (math/cos (+ (mean-obliquity mjd-tt) deps)))))

;; ------------------------------------------------------------ Earth rotation

(defn gast
  "Greenwich Apparent Sidereal Time, radians.

  Mean sidereal time reckons from the mean equinox; apparent sidereal time
  from the true one, which nutation displaces. The gap between them is the
  equation of the equinoxes, up to about a second of time.

  Note the two arguments. The rotation angle is a UT1 quantity -- it is
  where the Earth actually is -- while nutation is a TT quantity, being
  dynamics. They differ by 69 seconds today, and the Earth turns 465 m/s at
  the equator, so passing one where the other belongs puts a ground station
  30 km from where it is. That is larger than the entire precession
  correction over the first two decades of this century."
  [mjd-ut1 mjd-tt]
  (+ (time/gmst mjd-ut1) (equation-of-equinoxes mjd-tt)))

(defn earth-rotation
  "Rotation from the true equator and equinox of date to the Earth-fixed
  frame: the daily spin, and by far the largest of the four."
  [mjd-ut1 mjd-tt]
  (rz (gast mjd-ut1 mjd-tt)))

;; ------------------------------------------------------------- polar motion

(defn polar-motion
  "Rotation for the wander of the rotation pole within the Earth itself,
  `xp` and `yp` in radians.

  The pole is not fixed in the crust: it circles by some 0.3 arcseconds --
  around 9 meters at the surface -- in a 435-day Chandler wobble beating
  against an annual term. Like dUT1 this can only be measured and published,
  never predicted, and zero is the honest default when it is unknown."
  [xp yp]
  (lin/mat-mul (ry (- xp)) (rx (- yp))))

;; --------------------------------------------------- the complete transform

(defn celestial->terrestrial
  "The full rotation from the celestial frame of J2000 to the Earth-fixed
  frame (M&G eq. 5.71):

    U = polar motion * Earth rotation * nutation * precession

  read right to left. Precession and nutation carry J2000 to the true
  equator and equinox of date, Earth rotation spins that into the meridian
  of Greenwich, and polar motion accounts for the axis wandering in the
  crust. Each is nearly the identity except the third, which turns a full
  circle a day."
  ([mjd-tt mjd-ut1] (celestial->terrestrial mjd-tt mjd-ut1 0.0 0.0))
  ([mjd-tt mjd-ut1 xp yp]
   (chain (polar-motion xp yp)
          (earth-rotation mjd-ut1 mjd-tt)
          (nutation mjd-tt)
          (precession mjd-tt))))

(defn terrestrial->celestial
  "The inverse, which for a rotation is simply the transpose."
  ([mjd-tt mjd-ut1] (lin/transpose (celestial->terrestrial mjd-tt mjd-ut1)))
  ([mjd-tt mjd-ut1 xp yp] (lin/transpose (celestial->terrestrial mjd-tt mjd-ut1 xp yp))))

;; ------------------------------------------------------------------ caching

(def ^:private frame-cache (atom nil))

(def pn-bucket
  "How coarsely the slowly varying parts may be held, in days.

  A quarter of an hour, which costs at most about 4 centimeters at the
  Earth's surface -- far below the truncation of the nutation series itself,
  and below anything an orbit model cares about. The daily rotation is never
  cached and keeps its full resolution."
  0.01)

(defn- slow-parts
  "The precession-nutation product and the equation of the equinoxes, held
  together on a coarse grid.

  Both must be cached, not just the matrix. Apparent sidereal time needs the
  equation of the equinoxes, which runs the same expensive series -- cache
  one without the other and the saving is largely undone."
  [mjd-tt]
  (let [k (math/round (/ mjd-tt pn-bucket))]
    (or (when-let [[ck v] @frame-cache] (when (= ck k) v))
        (let [v [(lin/mat-mul (nutation mjd-tt) (precession mjd-tt))
                 (equation-of-equinoxes mjd-tt)]]
          (reset! frame-cache [k v])
          v))))

(defn precession-nutation
  "Nutation composed with precession, on the coarse grid.

  Worth caching because it is the expensive part and the slow part at once:
  the series behind it costs a hundred times a sidereal-time evaluation and
  changes a hundred thousand times more slowly. One entry suffices, since an
  integrator walks forward and asks for the same bucket repeatedly."
  [mjd-tt]
  (first (slow-parts mjd-tt)))

(defn celestial->terrestrial-cached
  "As `celestial->terrestrial`, but reusing the slowly varying parts across
  nearby epochs. The rotation angle itself is always recomputed."
  ([mjd-tt mjd-ut1] (celestial->terrestrial-cached mjd-tt mjd-ut1 0.0 0.0))
  ([mjd-tt mjd-ut1 xp yp]
   (let [[pn eqeq] (slow-parts mjd-tt)]
     (chain (polar-motion xp yp)
            (rz (+ (time/gmst mjd-ut1) eqeq))
            pn))))

;; ------------------------------------------------------------------- drawing

(defn y-up
  "EME2000 components into a frame with y up where EME2000 has z: the
  frame a y-up renderer like three.js draws in, with the equinox still
  along x and the celestial pole along y.

  A rotation, the quarter turn `rx(pi/2)` -- `[x z -y]`. The tempting
  `[x z y]`, swapping two axes, is a reflection, and drew every orbit
  going the wrong way round and every constellation back to front."
  [[x y z]]
  [x z (- y)])
