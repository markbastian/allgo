(ns allgo.astro.saturn-moons
  "Where Saturn's eight classical satellites appear beside the planet
  (Meeus, *Astronomical Algorithms*, chapter 46, after Dourneau's theory).

  Each satellite gets its own treatment, which is the point of the
  chapter: Mimas, Enceladus, Tethys and Dione move in near-circles
  perturbed by their resonances with one another; Rhea, Titan, Hyperion
  and Iapetus in precessing ellipses, Titan's disturbed by the Sun and
  Hyperion's -- locked in a 4:3 resonance with Titan -- by Titan. The
  orbits are referred to the B1950 equinox, as the theory was fitted, and
  then turned to the Earth's view.

  Positions are rectangular coordinates on the sky in Saturn's equatorial
  radii: X positive toward the west, Y toward the north. Times are MJD
  (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.precession :as precession]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.numerics.interpolation :refer [horner]]
            [clojure.math :as math]))

(def names [:mimas :enceladus :tethys :dione :rhea :titan :hyperion :iapetus])

(def ^:private d c/degrees)

(def ^:private s1 (math/sin (* 28.0817 d)))
(def ^:private c1 (math/cos (* 28.0817 d)))
(def ^:private s2 (math/sin (* 168.8112 d)))
(def ^:private c2 (math/cos (* 168.8112 d)))

(defn- arguments
  "Meeus's auxiliary times t1..t11 and angles W0..W8 at `jde`."
  [jde]
  (let [t1 (- jde 2411093.0)
        t2 (/ t1 365.25)
        t3 (+ (/ (- jde 2433282.423) 365.25) 1950.0)
        t4 (- jde 2411368.0)
        t5 (/ t4 365.25)
        t6 (- jde 2415020.0)
        t7 (/ t6 36525.0)
        t8 (/ t6 365.25)
        t9 (/ (- jde 2442000.5) 365.25)
        t10 (- jde 2409786.0)
        t11 (/ t10 36525.0)]
    {:t1 t1 :t2 t2 :t4 t4 :t6 t6 :t7 t7 :t8 t8 :t9 t9 :t10 t10 :t11 t11
     :W0 (* 5.095 d (- t3 1866.39))
     :W1 (+ (* 74.4 d) (* 32.39 d t2))
     :W2 (+ (* 134.3 d) (* 92.62 d t2))
     :W3 (- (* 42.0 d) (* 0.5118 d t5))
     :W4 (+ (* 276.59 d) (* 0.5118 d t5))
     :W5 (+ (* 267.2635 d) (* 1222.1136 d t7))
     :W6 (+ (* 175.4762 d) (* 1221.5515 d t7))
     :W7 (+ (* 2.4891 d) (* 0.002435 d t7))
     :W8 (- (* 113.35 d) (* 0.2597 d t7))
     :e1 (- 0.05589 (* 0.000346 t7))}))

(defn- elliptic-orbit
  "Meeus's subroutine for the four outer satellites: an ellipse of mean
  longitude `l'`, perihelion `p`, eccentricity `e` and semi-major axis `a`
  on a plane of node `om` and inclination `i` to the ecliptic, referred to
  Saturn's equator. Returns `{:lam :r :gamma :om}`."
  [l' p e a om i]
  (let [M  (- l' p)
        e2 (* e e) e3 (* e2 e) e4 (* e2 e2) e5 (* e3 e2)
        C  (+ (* (+ (* 2 e) (* -0.25 e3) (* 0.0520833333 e5)) (math/sin M))
              (* (- (* 1.25 e2) (* 0.458333333 e4)) (math/sin (* 2 M)))
              (* (- (* 1.083333333 e3) (* 0.671875 e5)) (math/sin (* 3 M)))
              (* 1.072917 e4 (math/sin (* 4 M)))
              (* 1.142708 e5 (math/sin (* 5 M))))
        g  (- om (* 168.8112 d))
        si (math/sin i) ci (math/cos i)
        sg (math/sin g) cg (math/cos g)
        a1 (* si sg)
        a2 (- (* c1 si cg) (* s1 ci))
        u  (math/atan2 a1 a2)
        psi (math/atan2 (* s1 sg) (- (* c1 si) (* s1 ci cg)))]
    {:r (/ (* a (- 1.0 e2)) (+ 1.0 (* e (math/cos (+ M C)))))
     :gamma (math/asin (math/hypot a1 a2))
     :om (+ (* 168.8112 d) u)
     :lam (- (+ l' C u) g psi)}))

(defn- mimas [{:keys [t1 t2 W0]}]
  (let [L (- (+ (* 127.64 d) (* 381.994497 d t1))
             (* 43.57 d (math/sin W0)) (* 0.72 d (math/sin (* 3 W0))) (* 0.02144 d (math/sin (* 5 W0))))
        M (- L (+ (* 106.1 d) (* 365.549 d t2)))
        C (+ (* 2.18287 d (math/sin M)) (* 0.025988 d (math/sin (* 2 M))) (* 0.00043 d (math/sin (* 3 M))))]
    {:lam (+ L C) :r (/ 3.06879 (+ 1.0 (* 0.01905 (math/cos (+ M C)))))
     :gamma (* 1.563 d) :om (- (* 54.5 d) (* 365.072 d t2))}))

(defn- enceladus [{:keys [t1 t2 W1 W2]}]
  (let [L (+ (* 200.317 d) (* 262.7319002 d t1) (* 0.25667 d (math/sin W1)) (* 0.20883 d (math/sin W2)))
        M (- L (+ (* 309.107 d) (* 123.44121 d t2)))
        C (+ (* 0.55577 d (math/sin M)) (* 0.00168 d (math/sin (* 2 M))))]
    {:lam (+ L C) :r (/ 3.94118 (+ 1.0 (* 0.00485 (math/cos (+ M C)))))
     :gamma (* 0.0262 d) :om (- (* 348.0 d) (* 151.95 d t2))}))

(defn- tethys [{:keys [t1 t2 W0]}]
  {:lam (+ (* 285.306 d) (* 190.69791226 d t1) (* 2.063 d (math/sin W0))
           (* 0.03409 d (math/sin (* 3 W0))) (* 0.001015 d (math/sin (* 5 W0))))
   :r 4.880998 :gamma (* 1.0976 d) :om (- (* 111.33 d) (* 72.2441 d t2))})

(defn- dione [{:keys [t1 t2 W1 W2]}]
  (let [L (- (+ (* 254.712 d) (* 131.53493193 d t1)) (* 0.0215 d (math/sin W1)) (* 0.01733 d (math/sin W2)))
        M (- L (+ (* 174.8 d) (* 30.82 d t2)))
        C (+ (* 0.24717 d (math/sin M)) (* 0.00033 d (math/sin (* 2 M))))]
    {:lam (+ L C) :r (/ 6.24871 (+ 1.0 (* 0.002157 (math/cos (+ M C)))))
     :gamma (* 0.0139 d) :om (- (* 232.0 d) (* 30.27 d t2))}))

(defn- rhea [{:keys [t1 t2 W3 W4]}]
  (let [p' (+ (* 342.7 d) (* 10.057 d t2))
        a1 (+ (* 0.000265 (math/sin p')) (* 0.001 (math/sin W4)))
        a2 (+ (* 0.000265 (math/cos p')) (* 0.001 (math/cos W4)))
        N  (- (* 345.0 d) (* 10.057 d t2))]
    (elliptic-orbit (+ (* 359.244 d) (* 79.6900472 d t1) (* 0.086754 d (math/sin N)))
                    (math/atan2 a1 a2) (math/hypot a1 a2) 8.725924
                    (+ (* 168.8034 d) (* 0.736936 d (math/sin N)) (* 0.041 d (math/sin W3)))
                    (+ (* 28.0362 d) (* 0.346898 d (math/cos N)) (* 0.0193 d (math/cos W3))))))

(defn- titan [{:keys [t4 W3 W4 W5 W6 W7 W8 e1]}]
  (let [L   (+ (* 261.1582 d) (* 22.57697855 d t4) (* 0.074025 d (math/sin W3)))
        i'  (+ (* 27.45141 d) (* 0.295999 d (math/cos W3)))
        om' (+ (* 168.66925 d) (* 0.628808 d (math/sin W3)))
        si' (math/sin i') ci' (math/cos i')
        a1  (* (math/sin W7) (math/sin (- om' W8)))
        a2  (- (* (math/cos W7) si') (* (math/sin W7) ci' (math/cos (- om' W8))))
        g0  (* 102.8623 d)
        psi (math/atan2 a1 a2)
        s   (math/hypot a1 a2)
        ;; the longitude of perihelion and g are coupled; three passes
        [peri g] (nth (iterate (fn [[_ g]]
                                 (let [peri (+ W4 (* 0.37515 d (- (math/sin (* 2 g)) (math/sin (* 2 g0)))))]
                                   [peri (- peri om' psi)]))
                               [nil (- W4 om' psi)])
                      3)
        e'  (+ 0.029092 (* 0.00019048 (- (math/cos (* 2 g)) (math/cos (* 2 g0)))))
        q   (* 2 (- W5 peri))
        b1  (* si' (math/sin (- om' W8)))
        b2  (- (* (math/cos W7) si' (math/cos (- om' W8))) (* (math/sin W7) ci'))
        th  (+ (math/atan2 b1 b2) W8)
        u   (+ (* 2 W5) (* -2 th) psi)
        h   (+ (* 0.9375 e' e' (math/sin q)) (* 0.1875 s s (math/sin (* 2 (- W5 th)))))]
    (elliptic-orbit (- L (* 0.254744 d (+ (* e1 (math/sin W6)) (* 0.75 e1 e1 (math/sin (* 2 W6))) h)))
                    (+ peri (* 0.159215 d (math/sin q)))
                    (+ e' (* 0.002778797 e' (math/cos q)))
                    20.216193
                    (+ om' (/ (* 0.031843 d s (math/sin u)) si'))
                    (+ i' (* 0.031843 d s (math/cos u))))))

(defn- hyperion [{:keys [t6 t8 t9 W3 W5]}]
  (let [eta  (+ (* 92.39 d) (* 0.5621071 d t6))
        zeta (- (* 148.19 d) (* 19.18 d t8))
        th   (- (* 184.8 d) (* 35.41 d t9))
        th'  (- th (* 7.5 d))
        as   (+ (* 176.0 d) (* 12.22 d t8))
        bs   (+ (* 8.0 d) (* 24.44 d t8))
        cs   (+ bs (* 5.0 d))
        peri (- (* 69.898 d) (* 18.67088 d t8))
        phi  (* 2 (- peri W5))
        chi  (- (* 94.9 d) (* 2.292 d t8))
        [sn cn] [(math/sin eta) (math/cos eta)]
        [sz cz] [(math/sin zeta) (math/cos zeta)]
        cos* math/cos sin* math/sin]
    (elliptic-orbit
     (+ (* 177.047 d) (* 16.91993829 d t6) (* 0.15648 d (sin* chi)) (* 9.142 d sn)
        (* 0.007 d (sin* (* 2 eta))) (* -0.014 d (sin* (* 3 eta))) (* 0.2275 d (sin* (+ zeta eta)))
        (* 0.2112 d (sin* (- zeta eta))) (* -0.26 d sz) (* -0.0098 d (sin* (* 2 zeta)))
        (* -0.013 d (sin* as)) (* 0.017 d (sin* bs)) (* -0.0303 d (sin* phi)))
     (+ peri (* 0.15648 d (sin* chi)) (* -0.4457 d sn) (* -0.2657 d (sin* (+ zeta eta)))
        (* -0.3573 d (sin* (- zeta eta))) (* -12.872 d sz) (* 1.668 d (sin* (* 2 zeta)))
        (* -0.2419 d (sin* (* 3 zeta))) (* -0.07 d (sin* phi)))
     (+ 0.103458 (* -0.004099 cn) (* -0.000167 (cos* (+ zeta eta))) (* 0.000235 (cos* (- zeta eta)))
        (* 0.02303 cz) (* -0.00212 (cos* (* 2 zeta))) (* 0.000151 (cos* (* 3 zeta)))
        (* 0.00013 (cos* phi)))
     (+ 24.50601 (* -0.08686 cn) (* -0.00166 (cos* (+ zeta eta))) (* 0.00175 (cos* (- zeta eta))))
     (+ (* 168.6812 d) (* 1.40136 d (cos* chi)) (* 0.68599 d (sin* W3)) (* -0.0392 d (sin* cs))
        (* 0.0366 d (sin* th')))
     (+ (* 27.3347 d) (* 0.6434886 d (cos* chi)) (* 0.315 d (cos* W3)) (* 0.018 d (cos* th))
        (* -0.018 d (cos* cs))))))

(defn- iapetus [{:keys [t4 t7 t10 t11 W4 W5]}]
  (let [L    (+ (* 261.1582 d) (* 22.57697855 d t4))
        peri' (+ (* 91.796 d) (* 0.562 d t7))
        psi  (- (* 4.367 d) (* 0.195 d t7))
        th   (- (* 146.819 d) (* 3.198 d t7))
        phi  (+ (* 60.47 d) (* 1.521 d t7))
        PHI  (- (* 205.055 d) (* 2.091 d t7))
        e'   (+ 0.028298 (* 0.001156 t11))
        peri0 (+ (* 352.91 d) (* 11.71 d t11))
        mu   (+ (* 76.3852 d) (* 4.53795125 d t10))
        i'   (* d (horner t11 [18.4602 -0.9518 -0.072 0.0054]))
        om'  (* d (horner t11 [143.198 -3.919 0.116 0.008]))
        l    (- mu peri0)
        g    (- peri0 om' psi)
        g1   (- peri0 om' phi)
        ls   (- W5 peri')
        gs   (- peri' th)
        lT   (- L W4)
        gT   (- W4 PHI)
        u1   (* 2 (- (+ l g) ls gs))
        u2   (- (+ l g1) lT gT)
        u3   (+ l (* 2 (- g ls gs)))
        u4   (- (+ lT gT) g1)
        u5   (* 2 (+ ls gs))
        [sin* cos*] [math/sin math/cos]
        w    (+ (* 0.08077 d (sin* (- g1 gT))) (* 0.02139 d (sin* (- u5 (* 2 g))))
                (* -0.00676 d (sin* u3)) (* 0.0138 d (sin* l)) (* 0.01632 d (sin* (+ l u2)))
                (* 0.03547 d (sin* u4)))
        w'   (+ (* 0.04204 d (sin* (+ u5 psi))) (* 0.00235 d (sin* (+ l g1 lT gT phi)))
                (* 0.00358 d (sin* (+ u2 phi))))]
    (elliptic-orbit
     (+ mu (* -0.04299 d (sin* u2)) (* -0.00789 d (sin* u1)) (* -0.06312 d (sin* ls))
        (* -0.00295 d (sin* (* 2 ls))) (* -0.02231 d (sin* u5)) (* 0.0065 d (sin* (+ u5 psi))))
     (+ peri0 (/ w e'))
     (+ e' (* -0.0014097 (cos* (- g1 gT))) (* 0.0003733 (cos* (- u5 (* 2 g))))
        (* 0.000118 (cos* u3)) (* 0.0002408 (cos* l)) (* 0.0002849 (cos* (+ l u2)))
        (* 0.000619 (cos* u4)))
     (+ 58.935028 (* 0.004638 (cos* u1)) (* 0.058222 (cos* u2)))
     (+ om' (/ w' (sin* i')))
     (+ i' (* 0.04204 d (cos* (+ u5 psi))) (* 0.00235 d (cos* (+ l g1 lT gT phi)))
        (* 0.0036 d (cos* (+ u2 phi)))))))

(def ^:private perspective
  [20947.0 23715.0 26382.0 29876.0 35313.0 53800.0 59222.0 91820.0])

(defn positions
  "`{satellite [x y]}` for the eight satellites at `mjd-tt`, in Saturn
  radii; see `names`."
  [mjd-tt]
  (let [[s b R] (solar/geometric mjd-tt)
        earth [(* R (math/cos b) (math/cos s)) (* R (math/cos b) (math/sin s)) (* R (math/sin b))]
        [[x y z] delta tau]
        (loop [delta 9.0 i 0]
          (let [tau (elliptic/light-time delta)
                [l b r] (vsop87/heliocentric :saturn (- mjd-tt tau))
                [l b] (vsop87/->fk5 [l b] (- mjd-tt tau))
                xyz (mapv + [(* r (math/cos b) (math/cos l)) (* r (math/cos b) (math/sin l)) (* r (math/sin b))]
                          earth)
                delta' (math/sqrt (reduce + (map * xyz xyz)))]
            (if (or (< (abs (- delta' delta)) 1e-9) (> i 10))
              [xyz delta' (elliptic/light-time delta')]
              (recur delta' (inc i)))))
        [lam0 bet0] (precession/ecliptic [(math/atan2 y x) (math/atan (/ z (math/hypot x y)))]
                                         mjd-tt (time/besselian-epoch->mjd 1950.0))
        args (arguments (+ (- mjd-tt tau) c/jd-mjd-offset))
        orbits (mapv #(% args) [mimas enceladus tethys dione rhea titan hyperion iapetus])
        xyz (conj (mapv (fn [{:keys [lam r gamma om]}]
                          (let [u (- lam om) w (- om (* 168.8112 d))]
                            [(* r (- (* (math/cos u) (math/cos w)) (* (math/sin u) (math/cos gamma) (math/sin w))))
                             (* r (+ (* (math/sin u) (math/cos w) (math/cos gamma)) (* (math/cos u) (math/sin w))))
                             (* r (math/sin u) (math/sin gamma))]))
                        orbits)
                  [0.0 0.0 1.0])
        view (mapv (fn [[X Y Z]]
                     (let [a X
                           b (- (* c1 Y) (* s1 Z))
                           c (+ (* s1 Y) (* c1 Z))
                           [a b] [(- (* c2 a) (* s2 b)) (+ (* s2 a) (* c2 b))]
                           A (- (* a (math/sin lam0)) (* b (math/cos lam0)))
                           b (+ (* a (math/cos lam0)) (* b (math/sin lam0)))]
                       [A (+ (* b (math/cos bet0)) (* c (math/sin bet0)))
                        (- (* c (math/cos bet0)) (* b (math/sin bet0)))]))
                   xyz)
        [A0 _ C0] (peek view)
        D (math/atan2 A0 C0)]
    (zipmap names
            (map (fn [[A B C] {:keys [r]} k]
                   (let [x (- (* A (math/cos D)) (* C (math/sin D)))
                         y (+ (* A (math/sin D)) (* C (math/cos D)))
                         q (/ x r)
                         x (+ x (* (/ (abs B) k) (math/sqrt (- 1.0 (* q q)))))
                         W (/ delta (+ delta (/ B 2475.0)))]
                     [(* x W) (* y W)]))
                 (pop view) orbits perspective))))
