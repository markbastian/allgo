(ns allgo.astro.planet-orbits
  "The planets' orbits: their mean elements, when they pass perihelion,
  aphelion and their nodes, and Pluto, which VSOP87 leaves out (Meeus,
  *Astronomical Algorithms*, chapters 31 and 37 to 39).

  The mean elements are polynomials in time fitted to VSOP87 with the
  periodic terms averaged away -- the orbit a planet would follow if the
  others stopped pulling on it this century. `allgo.astro.planets` has a
  different such set, Standish's, which is fitted to the positions rather
  than derived from the theory; both are for drawing orbits and for
  first guesses, and `allgo.astro.vsop87` is for positions.

  Times are MJD (TT); angles radians; distances AU."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.precession :as precession]
            [allgo.astro.time :as time]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.interpolation :as interp]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

;; ------------------------------------------------------ mean elements (31)

(def ^:private mean-elements-of-date
  "Meeus's table 31.A: mean longitude, semi-major axis, eccentricity,
  inclination, longitude of the node and longitude of perihelion, as
  polynomials in Julian centuries from J2000, referred to the mean
  ecliptic and equinox of date. Degrees and AU."
  {:mercury {:L [252.250906 149474.0722491 0.0003035 0.000000018]
             :a [0.38709831]
             :e [0.20563175 0.000020407 -0.0000000283 -0.00000000018]
             :i [7.004986 0.0018215 -0.0000181 0.000000056]
             :node [48.330893 1.1861883 0.00017542 0.000000215]
             :peri [77.456119 1.5564776 0.00029544 0.000000009]}
   :venus   {:L [181.979801 58519.2130302 0.00031014 0.000000015]
             :a [0.72332982]
             :e [0.00677192 -0.000047765 0.0000000981 0.00000000046]
             :i [3.394662 0.0010037 -0.00000088 -0.000000007]
             :node [76.67992 0.9011206 0.00040618 -0.000000093]
             :peri [131.563703 1.4022288 -0.00107618 -0.000005678]}
   :earth   {:L [100.466457 36000.7698278 0.00030322 0.00000002]
             :a [1.000001018]
             :e [0.01670863 -0.000042037 -0.0000001267 0.00000000014]
             :i [0.0]
             :node [0.0]
             :peri [102.937348 1.7195366 0.00045688 -0.000000018]}
   :mars    {:L [355.433 19141.6964471 0.00031052 0.000000016]
             :a [1.523679342]
             :e [0.09340065 0.000090484 -0.0000000806 -0.00000000025]
             :i [1.849726 -0.0006011 0.00001276 -0.000000007]
             :node [49.558093 0.7720959 0.00001557 0.000002267]
             :peri [336.060234 1.8410449 0.00013477 0.000000536]}
   :jupiter {:L [34.351519 3036.3027748 0.0002233 0.000000037]
             :a [5.202603209 0.0000001913]
             :e [0.04849793 0.000163225 -0.0000004714 -0.00000000201]
             :i [1.303267 -0.0054965 0.00000466 -0.000000002]
             :node [100.464407 1.0209774 0.00040315 0.000000404]
             :peri [14.331207 1.6126352 0.00103042 -0.000004464]}
   :saturn  {:L [50.077444 1223.5110686 0.00051908 -0.00000003]
             :a [9.554909192 -0.0000021390 0.000000004]
             :e [0.05554814 -0.000346641 -0.0000006436 0.0000000034]
             :i [2.488879 -0.0037362 -0.00001519 0.000000087]
             :node [113.665503 0.877088 -0.00012176 -0.000002249]
             :peri [93.057237 1.9637613 0.00083753 0.000004928]}
   :uranus  {:L [314.055005 429.8640561 0.0003039 0.000000026]
             :a [19.218446062 -0.0000000372 0.00000000098]
             :e [0.04638122 -0.000027293 0.0000000789 0.00000000024]
             :i [0.773197 0.0007744 0.00003749 -0.000000092]
             :node [74.005957 0.5211278 0.00133947 0.000018484]
             :peri [173.005291 1.486379 0.00021406 0.000000434]}
   :neptune {:L [304.348665 219.8833092 0.00030882 0.000000018]
             :a [30.110386869 -0.0000001663 0.00000000069]
             :e [0.00945575 0.000006033 0.0 -0.00000000005]
             :i [1.769953 -0.0093082 -0.00000708 0.000000027]
             :node [131.784057 1.1022039 0.00025952 -0.000000637]
             :peri [48.120276 1.4262957 0.00038434 0.00000002]}})

(defn mean-elements
  "Mean elements of `planet` at `mjd-tt` referred to the ecliptic and
  equinox of date: `{:L :a :e :i :raan :argp :peri :M}` -- mean
  longitude, semi-major axis, eccentricity, inclination, node, argument
  and longitude of perihelion, and mean anomaly.

  The Earth's orbit defines the ecliptic of date, so its inclination is
  zero and its node undefined; it is reported as zero."
  [planet mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        at (fn [k] (poly/horner-ascending (get-in mean-elements-of-date [planet k]) T))
        L (c/deg (at :L)) node (c/deg (at :node)) peri (c/deg (at :peri))]
    {:L (am/wrap-2pi L) :a (at :a) :e (at :e) :i (c/deg (at :i))
     :raan (am/wrap-2pi node) :peri (am/wrap-2pi peri)
     :argp (am/wrap-2pi (- peri node)) :M (am/wrap-2pi (- L peri))}))

(defn- orbit-point
  "Unit vector, ecliptic, at argument of latitude `u` along an orbit."
  [raan i u]
  (let [so (math/sin raan) co (math/cos raan) su (math/sin u) cu (math/cos u)]
    [(- (* co cu) (* so su (math/cos i)))
     (+ (* so cu) (* co su (math/cos i)))
     (* su (math/sin i))]))

(defn- ->lon-lat [[x y z]] [(math/atan2 y x) (math/atan2 z (math/hypot x y))])

(defn- ->vector [[l b]] [(* (math/cos b) (math/cos l)) (* (math/cos b) (math/sin l)) (math/sin b)])

(defn mean-elements-J2000
  "`mean-elements` referred to the ecliptic and equinox of J2000 -- Meeus's
  table 31.B, here derived from 31.A by precessing the orbit: two points
  of it fix its plane in the new frame, and the perihelion fixes its
  orientation within that plane. Unlike the element formulae of chapter
  24 this is well defined for the Earth, whose orbit is the ecliptic of
  date and so has neither inclination nor node to start from."
  [planet mjd-tt]
  (let [{:keys [L peri raan i argp] :as el} (mean-elements planet mjd-tt)
        move (fn [v] (->vector (precession/ecliptic (->lon-lat v) mjd-tt c/mjd-J2000)))
        p1 (move (orbit-point raan i 0.0))
        p2 (move (orbit-point raan i (/ math/PI 2)))
        pp (move (orbit-point raan i argp))
        [nx ny nz :as n] (v3/cross p1 p2)
        raan' (am/wrap-2pi (math/atan2 nx (- ny)))
        node  [(math/cos raan') (math/sin raan') 0.0]
        argp' (am/wrap-2pi (math/atan2 (v3/dot (v3/cross n node) pp) (v3/dot node pp)))
        peri' (+ raan' argp')]
    (assoc el
           :i (math/acos nz) :raan raan' :argp argp'
           :peri (am/wrap-2pi peri')
           :L (am/wrap-2pi (+ L (- peri' peri))))))

;; ------------------------------------------ perihelion and aphelion (38)

(def ^:private apsis-k
  "k = factor * (year - offset): perihelion at integer k, aphelion at
  half-integer; then the JDE is a polynomial in k (Meeus table 38.A)."
  {:mercury [4.15201 2000.12 [2451590.257 87.96934963]]
   :venus   [1.62549 2000.53 [2451738.233 224.7008188 -0.0000000327]]
   :earth   [0.99997 2000.01 [2451547.507 365.2596358 0.0000000156]]
   :mars    [0.53166 2001.78 [2452195.026 686.9957857 -0.0000001187]]
   :jupiter [0.0843 2011.2 [2455636.936 4332.897065 0.0001367]]
   :saturn  [0.03393 2003.52 [2452830.12 10764.21676 0.000827]]
   :uranus  [0.0119 2051.1 [2470213.5 30694.8767 -0.00541]]
   :neptune [0.00607 2047.5 [2468895.1 60190.33 0.03429]]})

(def ^:private earth-apsis-correction
  "The Moon swings the Earth about the barycenter by up to 4700 km, which
  moves its perihelion by a day or more; five terms in k (Meeus 38.3)."
  {:args [[328.41 132.788585] [316.13 584.903153] [346.20 450.380738]
          [136.95 659.306737] [249.52 329.653368]]
   :perihelion [1.278 -0.055 -0.091 -0.056 -0.045]
   :aphelion   [-1.352 0.061 0.062 0.029 0.031]})

(defn- apsis-mean [planet year aphelion?]
  (let [[f off poly] (apsis-k (if (= planet :earth-moon) :earth planet))
        x (* f (- year off))
        k (if aphelion? (+ (math/floor x) 0.5) (math/floor (+ x 0.5)))
        jd (poly/horner-ascending poly k)
        jd (if (= planet :earth)
             (let [{:keys [args] :as corr} earth-apsis-correction]
               (+ jd (reduce + (map (fn [[a b] cc] (* cc (math/sin (c/deg (+ a (* b k))))))
                                    args (corr (if aphelion? :aphelion :perihelion))))))
             jd)]
    (- jd c/jd-mjd-offset)))

(defn perihelion
  "MJD (TT) of the perihelion of `planet` nearest the decimal `year`, from
  the mean orbit. `planet` may also be :earth-moon, the barycenter, whose
  perihelion is the unperturbed one. The mean orbit can be days out for
  the outer planets, whose perihelia the others' pulls shift: Jupiter's
  1981 aphelion by nine days. `apsis-exact` finds the true one."
  [planet year]
  (apsis-mean planet year false))

(defn aphelion
  "MJD (TT) of the aphelion of `planet` nearest the decimal `year`, from the
  mean orbit; see `perihelion`."
  [planet year]
  (apsis-mean planet year true))

(defn apsis-exact
  "`[mjd r]` of the perihelion (or when `aphelion?` the aphelion) of
  `planet` near `year`, as the extremum of its VSOP87 distance -- stepping
  `step` days from the mean-orbit estimate until the distance turns, then
  interpolating. A day suits the outer planets; the Earth, whose radius
  wobbles with the Moon, wants a step of a few minutes.

  Neptune's orbit is so nearly circular that the Sun's reflex about the
  barycenter, driven by Jupiter, makes its distance dither; its extremum
  is sought on both sides of the estimate and the deeper one kept."
  ([planet year aphelion?] (apsis-exact planet year aphelion? 1.0))
  ([planet year aphelion? step]
   (let [r (fn [t] (nth (vsop87/heliocentric planet t) 2))
         better? (if aphelion? > <)
         search (fn [t1]
                  (loop [t0 (- t1 step) t1 t1 t2 (+ t1 step)
                         r0 (r (- t1 step)) r1 (r t1) r2 (r (+ t1 step))]
                    (cond
                      (and (better? r1 r0) (better? r1 r2))
                      (interp/extremum (interp/table-3 t0 t2 [r0 r1 r2]))
                      (better? r2 r0)
                      (recur t1 t2 (+ t2 step) r1 r2 (r (+ t2 step)))
                      :else
                      (recur (- t0 step) t0 t1 (r (- t0 step)) r0 r1))))
         t (apsis-mean planet year aphelion?)]
     (if (= planet :neptune)
       (let [[_ ra :as a] (search (- t 5000.0))
             [_ rb :as b] (search (+ t 5000.0))]
         (if (better? ra rb) a b))
       (search t)))))

;; ---------------------------------------------- passages through the nodes (39)

(defn- elliptic-at [nu a e perihelion]
  (let [E (* 2.0 (math/atan (* (math/sqrt (/ (- 1.0 e) (+ 1.0 e))) (math/tan (* 0.5 nu)))))
        M (- E (* e (math/sin E)))
        n (/ c/gaussian-k (* a (math/sqrt a)))]
    [(+ perihelion (/ M n)) (* a (- 1.0 (* e (math/cos E))))]))

(defn node-passage
  "`[mjd r]` of the passage through its ascending (or when `descending?`,
  descending) node of a body on an elliptic orbit of semi-major axis `a`
  AU, eccentricity `e` and argument of perihelion `argp`, whose perihelion
  is at MJD `perihelion` -- the moments it crosses the ecliptic, and so the
  only times it can be occulted by or transit anything on it."
  ([a e argp perihelion] (node-passage a e argp perihelion false))
  ([a e argp perihelion descending?]
   (elliptic-at (if descending? (- math/PI argp) (- argp)) a e perihelion)))

(defn node-passage-parabolic
  "`[mjd r]` of a node passage on a parabolic orbit of perihelion distance
  `q`."
  ([q argp perihelion] (node-passage-parabolic q argp perihelion false))
  ([q argp perihelion descending?]
   (let [s (math/tan (* 0.5 (if descending? (- math/PI argp) (- argp))))]
     [(+ perihelion (* 27.403895 s (+ (* s s) 3.0) q (math/sqrt q)))
      (* q (+ 1.0 (* s s)))])))

;; ---------------------------------------------------------------- Pluto (37)

(def ^:private pluto-terms
  "Meeus's table 37.A, from Chapront and Francou's fit to the DE200
  ephemeris: multipliers of Jupiter's, Saturn's and Pluto's mean
  longitudes, then the sine and cosine coefficients of longitude and
  latitude (degrees) and of radius (AU)."
  [[0 0 1 -19.799805 19.850055 -5.452852 -14.974862 6.6865439 6.8951812]
   [0 0 2 0.897144 -4.954829 3.527812 1.67279 -1.1827535 -0.0332538]
   [0 0 3 0.611149 1.211027 -1.050748 0.327647 0.1593179 -0.143889]
   [0 0 4 -0.341243 -0.189585 0.17869 -0.292153 -0.0018444 0.048322]
   [0 0 5 0.129287 -0.034992 0.01865 0.10034 -0.0065977 -0.0085431]
   [0 0 6 -0.038164 0.030893 -0.030697 -0.025823 0.0031174 -0.0006032]
   [0 1 -1 0.020442 -0.009987 0.004878 0.011248 -0.0005794 0.0022161]
   [0 1 0 -0.004063 -0.005071 0.000226 -0.000064 0.0004601 0.0004032]
   [0 1 1 -0.006016 -0.003336 0.00203 -0.000836 -0.0001729 0.0000234]
   [0 1 2 -0.003956 0.003039 0.000069 -0.000604 -0.0000415 0.0000702]
   [0 1 3 -0.000667 0.003572 -0.000247 -0.000567 0.0000239 0.0000723]
   [0 2 -2 0.001276 0.000501 -0.000057 0.000001 0.0000067 -0.0000067]
   [0 2 -1 0.001152 -0.000917 -0.000122 0.000175 0.0001034 -0.0000451]
   [0 2 0 0.00063 -0.001277 -0.000049 -0.000164 -0.0000129 0.0000504]
   [1 -1 0 0.002571 -0.000459 -0.000197 0.000199 0.000048 -0.0000231]
   [1 -1 1 0.000899 -0.001449 -0.000025 0.000217 0.0000002 -0.0000441]
   [1 0 -3 -0.001016 0.001043 0.000589 -0.000248 -0.0003359 0.0000265]
   [1 0 -2 -0.002343 -0.001012 -0.000269 0.000711 0.0007856 -0.0007832]
   [1 0 -1 0.007042 0.000788 0.000185 0.000193 0.0000036 0.0045763]
   [1 0 0 0.001199 -0.000338 0.000315 0.000807 0.0008663 0.0008547]
   [1 0 1 0.000418 -0.000067 -0.00013 -0.000043 -0.0000809 -0.0000769]
   [1 0 2 0.00012 -0.000274 0.000005 0.000003 0.0000263 -0.0000144]
   [1 0 3 -0.00006 -0.000159 0.000002 0.000017 -0.0000126 0.0000032]
   [1 0 4 -0.000082 -0.000029 0.000002 0.000005 -0.0000035 -0.0000016]
   [1 1 -3 -0.000036 -0.000029 0.000002 0.000003 -0.0000019 -0.0000004]
   [1 1 -2 -0.00004 0.000007 0.000003 0.000001 -0.0000015 0.0000008]
   [1 1 -1 -0.000014 0.000022 0.000002 -0.000001 -0.0000004 0.0000012]
   [1 1 0 0.000004 0.000013 0.000001 -0.000001 0.0000005 0.0000006]
   [1 1 1 0.000005 0.000002 0 -0.000001 0.0000003 0.0000001]
   [1 1 3 -0.000001 0 0 0 0.0000006 -0.0000002]
   [2 0 -6 0.000002 0 0 -0.000002 0.0000002 0.0000002]
   [2 0 -5 -0.000004 0.000005 0.000002 0.000002 -0.0000002 -0.0000002]
   [2 0 -4 0.000004 -0.000007 -0.000007 0 0.0000014 0.0000013]
   [2 0 -3 0.000014 0.000024 0.00001 -0.000008 -0.0000063 0.0000013]
   [2 0 -2 -0.000049 -0.000034 -0.000003 0.00002 0.0000136 -0.0000236]
   [2 0 -1 0.000163 -0.000048 0.000006 0.000005 0.0000273 0.0001065]
   [2 0 0 0.000009 -0.000024 0.000014 0.000017 0.0000251 0.0000149]
   [2 0 1 -0.000004 0.000001 -0.000002 0 -0.0000025 -0.0000009]
   [2 0 2 -0.000003 0.000001 0 0 0.0000009 -0.0000002]
   [2 0 3 0.000001 0.000003 0 0 -0.0000008 0.0000007]
   [3 0 -2 -0.000003 -0.000001 0 0.000001 0.0000002 -0.000001]
   [3 0 -1 0.000005 -0.000003 0 0 0.0000019 0.0000035]
   [3 0 0 0 0 0.000001 0 0.000001 0.0000003]])

(defn pluto
  "Pluto's heliocentric `[lon lat r]` at `mjd-tt`, referred to the ecliptic
  and equinox of J2000 -- not of date, unlike VSOP87 -- in radians and AU.

  Pluto has no VSOP87 series: it was fitted separately, and only from
  1885 to 2099, the years the numerical ephemeris behind it spans. Outside
  them this is not merely less accurate but meaningless, so it returns nil
  there. Within them, a few hundredths of an arcsecond from the ephemeris
  it was fitted to; against JPL's DE441 the difference grows smoothly
  from that near 2000 to some 3 arcseconds by 2100, the older
  ephemeris's own error."
  [mjd-tt]
  (let [T (time/centuries-J2000 mjd-tt)
        y (time/mjd->julian-epoch mjd-tt)]
    (when (<= 1885.0 y 2100.0)
      (let [J (c/deg (+ 34.35 (* 3034.9057 T)))
            S (c/deg (+ 50.08 (* 1222.1138 T)))
            P (c/deg (+ 238.96 (* 144.96 T)))
            [l b r] (reduce (fn [[l b r] [i j k la lb ba bb ra rb]]
                              (let [a  (+ (* i J) (* j S) (* k P))
                                    sa (math/sin a) ca (math/cos a)]
                                [(+ l (* la sa) (* lb ca))
                                 (+ b (* ba sa) (* bb ca))
                                 (+ r (* ra sa) (* rb ca))]))
                            [0.0 0.0 0.0] pluto-terms)]
        [(am/wrap-2pi (c/deg (+ l 238.958116 (* 144.96 T))))
         (c/deg (- b 3.908239))
         (+ r 40.7241346)]))))

(defn pluto-astrometric
  "Pluto's astrometric `[ra dec distance elongation]`, J2000, as
  `allgo.astro.elliptic/astrometric` gives for any body."
  [mjd-tt]
  (when (pluto mjd-tt)
    (let [eps eph/obliquity-J2000
          se (math/sin eps) ce (math/cos eps)]
      (elliptic/astrometric
       (fn [mjd]
         (let [[l b r] (pluto mjd)
               sl (math/sin l) cl (math/cos l) sb (math/sin b) cb (math/cos b)]
           [(* r cl cb)
            (* r (- (* sl cb ce) (* sb se)))
            (* r (+ (* sl cb se) (* sb ce)))]))
       mjd-tt))))
