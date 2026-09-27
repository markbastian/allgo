(ns allgo.astro.jupiter-moons
  "Where the four Galilean satellites appear beside Jupiter (Meeus,
  *Astronomical Algorithms*, chapter 44).

  Positions are rectangular coordinates on the sky in units of Jupiter's
  equatorial radius, centered on the planet: X positive toward the west
  (the direction of the sky's daily motion), Y positive toward the north
  -- Meeus's convention, and the one in which a satellite with X of 5 is
  seen five radii to the west.

  Two methods. `positions-low-precision` treats each orbit as a circle
  with its largest perturbation, and is good enough to identify which
  moon is which at the eyepiece. `positions` evaluates Lieske's E5 theory
  as Meeus abridges it: 150 periodic terms in longitude, latitude and
  radius, then a rotation of each orbit to the Earth's view, with the
  moons' light time and the perspective of the planet's disk -- a few
  hundredths of a Jovian radius, enough to time eclipses and transits to a
  minute.

  Times are MJD (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.elliptic :as elliptic]
            [allgo.astro.planet-orbits :as orbits]
            [allgo.astro.solar :as solar]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn- deg [x] (* x c/degrees))

;; --------------------------------------------------------- low precision

(defn positions-low-precision
  "`[[x y] ...]` of Io, Europa, Ganymede and Callisto, in Jupiter radii,
  from circular orbits with the main perturbation of each."
  [mjd-tt]
  (let [d  (- mjd-tt c/mjd-J2000)
        V  (deg (+ 172.74 (* 0.00111588 d)))
        M  (deg (+ 357.529 (* 0.9856003 d)))
        sV (math/sin V)
        N  (+ (deg (+ 20.02 (* 0.0830853 d))) (* (deg 0.329) sV))
        J  (- (deg (+ 66.115 (* 0.9025179 d))) (* (deg 0.329) sV))
        A  (+ (* (deg 1.915) (math/sin M)) (* (deg 0.020) (math/sin (* 2 M))))
        B  (+ (* (deg 5.555) (math/sin N)) (* (deg 0.168) (math/sin (* 2 N))))
        K  (- (+ J A) B)
        R  (- 1.00014 (* 0.01671 (math/cos M)) (* 0.00014 (math/cos (* 2 M))))
        r  (- 5.20872 (* 0.25208 (math/cos N)) (* 0.00611 (math/cos (* 2 N))))
        delta (math/sqrt (- (+ (* r r) (* R R)) (* 2 r R (math/cos K))))
        psi (math/asin (* (/ R delta) (math/sin K)))
        lam (+ (deg (+ 34.35 (* 0.083091 d))) (* (deg 0.329) sV) B)
        DS  (* (deg 3.12) (math/sin (+ lam (deg 42.8))))
        DE  (- DS (* (deg 2.22) (math/sin psi) (math/cos (+ lam (deg 22))))
               (* (deg 1.3) (/ (- r delta) delta) (math/sin (- lam (deg 100.5)))))
        dd  (- d (/ delta 173.0))
        u   (fn [a b] (+ (deg (+ a (* b dd))) psi (- B)))
        u1 (u 163.8069 203.4058646) u2 (u 358.414 101.2916335)
        u3 (u 5.7176 50.234518)     u4 (u 224.8092 21.48798)
        G  (deg (+ 331.18 (* 50.310482 dd)))
        H  (deg (+ 87.45 (* 21.569231 dd)))
        sde (math/sin DE)
        xy (fn [u r] [(* r (math/sin u)) (* (- r) (math/cos u) sde)])]
    [(xy (+ u1 (* (deg 0.473) (math/sin (* 2 (- u1 u2))))) (- 5.9057 (* 0.0244 (math/cos (* 2 (- u1 u2))))))
     (xy (+ u2 (* (deg 1.065) (math/sin (* 2 (- u2 u3))))) (- 9.3966 (* 0.0882 (math/cos (* 2 (- u2 u3))))))
     (xy (+ u3 (* (deg 0.165) (math/sin G))) (- 14.9883 (* 0.0216 (math/cos G))))
     (xy (+ u4 (* (deg 0.843) (math/sin H))) (- 26.3627 (* 0.1939 (math/cos H))))]))

;; --------------------------------------------------------------- E5 theory

(def ^:private sigma
  "The periodic terms in longitude of each satellite, degrees: a
  coefficient and its argument, as multipliers of the fundamental
  arguments (and `:deg`, a constant in degrees)."
  {:sigma1
   [[0.47259 {:l1 2 :l2 -2}]
    [-0.03478 {:pi3 1 :pi4 -1}]
    [0.01081 {:l2 1 :l3 -2 :pi3 1}]
    [0.00738 {:phi 1}]
    [0.00713 {:l2 1 :l3 -2 :pi2 1}]
    [-0.00674 {:pi1 1 :pi3 1 :PI -2 :G -2}]
    [0.00666 {:l2 1 :l3 -2 :pi4 1}]
    [0.00445 {:l1 1 :pi3 -1}]
    [-0.00354 {:l1 1 :l2 -1}]
    [-0.00317 {:psi 2 :PI -2}]
    [0.00265 {:l1 1 :pi4 -1}]
    [-0.00186 {:G 1}]
    [0.00162 {:pi2 1 :pi3 -1}]
    [0.00158 {:l1 4 :l2 -4}]
    [-0.00155 {:l1 1 :l3 -1}]
    [-0.00138 {:om3 1 :psi 1 :PI -2 :G -2}]
    [-0.00115 {:l1 2 :l2 -4 :om2 2}]
    [0.00089 {:pi2 1 :pi4 -1}]
    [0.00085 {:l1 1 :pi3 1 :PI -2 :G -2}]
    [0.00083 {:om2 1 :om3 -1}]
    [0.00053 {:om2 -1 :psi 1}]]
   :sigma2
   [[1.06476 {:l2 2 :l3 -2}]
    [0.04256 {:l1 1 :l2 -2 :pi3 1}]
    [0.03581 {:l2 1 :pi3 -1}]
    [0.02395 {:l1 1 :l2 -2 :pi4 1}]
    [0.01984 {:l2 1 :pi4 -1}]
    [-0.01778 {:phi 1}]
    [0.01654 {:l2 1 :pi2 -1}]
    [0.01334 {:l2 1 :l3 -2 :pi2 1}]
    [0.01294 {:pi3 1 :pi4 -1}]
    [-0.01142 {:l2 1 :l3 -1}]
    [-0.01057 {:G 1}]
    [-0.00775 {:psi 2 :PI -2}]
    [0.00524 {:l1 2 :l2 -2}]
    [-0.0046 {:l1 1 :l3 -1}]
    [0.00316 {:om3 1 :psi 1 :PI -2 :G -2}]
    [-0.00203 {:pi1 1 :pi3 1 :PI -2 :G -2}]
    [0.00146 {:om3 -1 :psi 1}]
    [-0.00145 {:G 2}]
    [0.00125 {:om4 -1 :psi 1}]
    [-0.00115 {:l1 1 :l3 -2 :pi3 1}]
    [-0.00094 {:l2 2 :om2 -2}]
    [0.00086 {:l1 2 :l2 -4 :om2 2}]
    [-0.00086 {:G -2 :G' 5 :deg 52.225}]
    [-0.00078 {:l2 1 :l4 -1}]
    [-0.00064 {:l3 3 :l4 -7 :pi4 4}]
    [0.00064 {:pi1 1 :pi4 -1}]
    [-0.00063 {:l1 1 :l3 -2 :pi4 1}]
    [0.00058 {:om3 1 :om4 -1}]
    [0.00056 {:psi 2 :PI -2 :G -2}]
    [0.00056 {:l2 2 :l4 -2}]
    [0.00055 {:l1 2 :l3 -2}]
    [0.00052 {:l3 3 :l4 -7 :pi3 1 :pi4 3}]
    [-0.00043 {:l1 1 :pi3 -1}]
    [0.00041 {:l2 5 :l3 -5}]
    [0.00041 {:pi4 1 :PI -1}]
    [0.00032 {:om2 1 :om3 -1}]
    [0.00032 {:l3 2 :PI -2 :G -2}]]
   :sigma3
   [[0.1649 {:l3 1 :pi3 -1}]
    [0.09081 {:l3 1 :pi4 -1}]
    [-0.06907 {:l2 1 :l3 -1}]
    [0.03784 {:pi3 1 :pi4 -1}]
    [0.01846 {:l3 2 :l4 -2}]
    [-0.0134 {:G 1}]
    [-0.01014 {:psi 2 :PI -2}]
    [0.00704 {:l2 1 :l3 -2 :pi3 1}]
    [-0.0062 {:l2 1 :l3 -2 :pi2 1}]
    [-0.00541 {:l3 1 :l4 -1}]
    [0.00381 {:l2 1 :l3 -2 :pi4 1}]
    [0.00235 {:om3 -1 :psi 1}]
    [0.00198 {:om4 -1 :psi 1}]
    [0.00176 {:phi 1}]
    [0.0013 {:l3 3 :l4 -3}]
    [0.00125 {:l1 1 :l3 -1}]
    [-0.00119 {:G -2 :G' 5 :deg 52.225}]
    [0.00109 {:l1 1 :l2 -1}]
    [-0.001 {:l3 3 :l4 -7 :pi4 4}]
    [0.00091 {:om3 1 :om4 -1}]
    [0.0008 {:l3 3 :l4 -7 :pi3 1 :pi4 3}]
    [-0.00075 {:l2 2 :l3 -3 :pi3 1}]
    [0.00072 {:pi1 1 :pi3 1 :PI -2 :G -2}]
    [0.00069 {:pi4 1 :PI -1}]
    [-0.00058 {:l3 2 :l4 -3 :pi4 1}]
    [-0.00057 {:l3 1 :l4 -2 :pi4 1}]
    [0.00056 {:l3 1 :pi3 1 :PI -2 :G -2}]
    [-0.00052 {:l2 1 :l3 -2 :pi1 1}]
    [-0.0005 {:pi2 1 :pi3 -1}]
    [0.00048 {:l3 1 :l4 -2 :pi3 1}]
    [-0.00045 {:l2 2 :l3 -3 :pi4 1}]
    [-0.00041 {:pi2 1 :pi4 -1}]
    [-0.00038 {:G 2}]
    [-0.00037 {:pi3 1 :pi4 -1 :om3 1 :om4 -1}]
    [-0.00032 {:l3 3 :l4 -7 :pi3 2 :pi4 2}]
    [0.0003 {:l3 4 :l4 -4}]
    [0.00029 {:l3 1 :pi4 1 :PI -2 :G -2}]
    [-0.00028 {:om3 1 :psi 1 :PI -2 :G -2}]
    [0.00026 {:l3 1 :PI -1 :G -1}]
    [0.00024 {:l2 1 :l3 -3 :l4 2}]
    [0.00021 {:l3 2 :PI -2 :G -2}]
    [-0.00021 {:l3 1 :pi2 -1}]
    [0.00017 {:l3 2 :pi3 -2}]]
   :sigma4
   [[0.84287 {:l4 1 :pi4 -1}]
    [0.03431 {:pi3 -1 :pi4 1}]
    [-0.03305 {:psi 2 :PI -2}]
    [-0.03211 {:G 1}]
    [-0.01862 {:l4 1 :pi3 -1}]
    [0.01186 {:om4 -1 :psi 1}]
    [0.00623 {:l4 1 :pi4 1 :PI -2 :G -2}]
    [0.00387 {:l4 2 :pi4 -2}]
    [-0.00284 {:G -2 :G' 5 :deg 52.225}]
    [-0.00234 {:pi4 -2 :psi 2}]
    [-0.00223 {:l3 1 :l4 -1}]
    [-0.00208 {:l4 1 :PI -1}]
    [0.00178 {:pi4 -2 :om4 1 :psi 1}]
    [0.00134 {:pi4 1 :PI -1}]
    [0.00125 {:l4 2 :PI -2 :G -2}]
    [-0.00117 {:G 2}]
    [-0.00112 {:l3 2 :l4 -2}]
    [0.00107 {:l3 3 :l4 -7 :pi4 4}]
    [0.00102 {:l4 1 :PI -1 :G -1}]
    [0.00096 {:l4 2 :om4 -1 :psi -1}]
    [0.00087 {:om4 -2 :psi 2}]
    [-0.00085 {:l3 3 :l4 -7 :pi3 1 :pi4 3}]
    [0.00085 {:l3 1 :l4 -2 :pi4 1}]
    [-0.00081 {:l4 2 :psi -2}]
    [0.00071 {:l4 1 :pi4 1 :PI -2 :G -3}]
    [0.00061 {:l1 1 :l4 -1}]
    [-0.00056 {:om3 -1 :psi 1}]
    [-0.00054 {:l3 1 :l4 -2 :pi3 1}]
    [0.00051 {:l2 1 :l4 -1}]
    [0.00042 {:psi 2 :PI -2 :G -2}]
    [0.00039 {:pi4 2 :om4 -2}]
    [0.00036 {:pi4 -1 :om4 -1 :psi 1 :PI 1}]
    [0.00035 {:G -1 :G' 2 :deg 188.37}]
    [-0.00035 {:l4 1 :pi4 -1 :psi -2 :PI 2}]
    [-0.00032 {:l4 1 :pi4 1 :PI -2 :G -1}]
    [0.0003 {:G -2 :G' 2 :deg 149.15}]
    [0.00029 {:l3 3 :l4 -7 :pi3 2 :pi4 2}]
    [0.00028 {:l4 1 :pi4 -1 :psi 2 :PI -2}]
    [-0.00028 {:l4 2 :om4 -2}]
    [-0.00027 {:pi3 1 :pi4 -1 :om3 1 :om4 -1}]
    [-0.00026 {:G -3 :G' 5 :deg 188.37}]
    [0.00025 {:om3 -1 :om4 1}]
    [-0.00025 {:l2 1 :l3 -3 :l4 2}]
    [-0.00023 {:l3 3 :l4 -3}]
    [0.00021 {:l4 2 :PI -2 :G -3}]
    [-0.00021 {:l3 2 :l4 -3 :pi4 1}]
    [0.00019 {:l4 1 :pi4 -1 :G -1}]
    [-0.00019 {:l4 2 :pi3 -1 :pi4 -1}]
    [-0.00018 {:l4 1 :pi4 -1 :G 1}]
    [-0.00016 {:l4 1 :pi3 1 :PI -2 :G -2}]]})

(def ^:private latitude
  "Terms of tan B, the latitude of each satellite above Jupiter's
  equator, in radians."
  [[[0.0006393 {:om1 -1 :L1 1}]
    [0.0001825 {:om2 -1 :L1 1}]
    [3.29e-5 {:om3 -1 :L1 1}]
    [-3.11e-5 {:psi -1 :L1 1}]
    [9.3e-6 {:om4 -1 :L1 1}]
    [7.5e-6 {:l2 -4 :om2 1 :L1 3 :S1 -1.9927}]
    [4.6e-6 {:psi 1 :PI -2 :G -2 :L1 1}]]
   [[0.0081004 {:om2 -1 :L2 1}]
    [0.0004512 {:om3 -1 :L2 1}]
    [-0.0003284 {:psi -1 :L2 1}]
    [0.000116 {:om4 -1 :L2 1}]
    [2.72e-5 {:l1 1 :l3 -2 :om2 1 :S2 1.0146}]
    [-1.44e-5 {:om1 -1 :L2 1}]
    [1.43e-5 {:psi 1 :PI -2 :G -2 :L2 1}]
    [3.5e-6 {:psi -1 :G 1 :L2 1}]
    [-2.8e-6 {:l1 1 :l3 -2 :om3 1 :S2 1.0146}]]
   [[0.0032402 {:om3 -1 :L3 1}]
    [-0.0016911 {:psi -1 :L3 1}]
    [0.0006847 {:om4 -1 :L3 1}]
    [-0.0002797 {:om2 -1 :L3 1}]
    [3.21e-5 {:psi 1 :PI -2 :G -2 :L3 1}]
    [5.1e-6 {:psi -1 :G 1 :L3 1}]
    [-4.5e-6 {:psi -1 :G -1 :L3 1}]
    [-4.5e-6 {:psi 1 :PI -2 :L3 1}]
    [3.7e-6 {:psi 1 :PI -2 :G -3 :L3 1}]
    [3e-6 {:l2 2 :om2 1 :L3 -3 :S3 4.03}]
    [-2.1e-6 {:l2 2 :om3 1 :L3 -3 :S3 4.03}]]
   [[-0.0076579 {:psi -1 :L4 1}]
    [0.0044134 {:om4 -1 :L4 1}]
    [-0.0005112 {:om3 -1 :L4 1}]
    [7.73e-5 {:psi 1 :PI -2 :G -2 :L4 1}]
    [1.04e-5 {:psi -1 :G 1 :L4 1}]
    [-1.02e-5 {:psi -1 :G -1 :L4 1}]
    [8.8e-6 {:psi 1 :PI -2 :G -3 :L4 1}]
    [-3.8e-6 {:psi 1 :PI -2 :G -1 :L4 1}]]])

(def ^:private radius
  "Terms of the relative departure of each radius from its mean."
  [[[-0.0041339 {:l1 2 :l2 -2}]
    [-3.87e-5 {:l1 1 :pi3 -1}]
    [-2.14e-5 {:l1 1 :pi4 -1}]
    [1.7e-5 {:l1 1 :l2 -1}]
    [-1.31e-5 {:l1 4 :l2 -4}]
    [1.06e-5 {:l1 1 :l3 -1}]
    [-6.6e-6 {:l1 1 :pi3 1 :PI -2 :G -2}]]
   [[0.0093848 {:l1 1 :l2 -1}]
    [-0.0003116 {:l2 1 :pi3 -1}]
    [-0.0001744 {:l2 1 :pi4 -1}]
    [-0.0001442 {:l2 1 :pi2 -1}]
    [5.53e-5 {:l2 1 :l3 -1}]
    [5.23e-5 {:l1 1 :l3 -1}]
    [-2.9e-5 {:l1 2 :l2 -2}]
    [1.64e-5 {:l2 2 :om2 -2}]
    [1.07e-5 {:l1 1 :l3 -2 :pi3 1}]
    [-1.02e-5 {:l2 1 :pi1 -1}]
    [-9.1e-6 {:l1 2 :l3 -2}]]
   [[-0.0014388 {:l3 1 :pi3 -1}]
    [-0.0007917 {:l3 1 :pi4 -1}]
    [0.0006342 {:l2 1 :l3 -1}]
    [-0.0001761 {:l3 2 :l4 -2}]
    [2.94e-5 {:l3 1 :l4 -1}]
    [-1.56e-5 {:l3 3 :l4 -3}]
    [1.56e-5 {:l1 1 :l3 -1}]
    [-1.53e-5 {:l1 1 :l2 -1}]
    [7e-6 {:l2 2 :l3 -3 :pi3 1}]
    [-5.1e-6 {:l3 1 :pi3 1 :PI -2 :G -2}]]
   [[-0.0073546 {:l4 1 :pi4 -1}]
    [0.0001621 {:l4 1 :pi3 -1}]
    [9.74e-5 {:l3 1 :l4 -1}]
    [-5.43e-5 {:l4 1 :pi4 1 :PI -2 :G -2}]
    [-2.71e-5 {:l4 2 :pi4 -2}]
    [1.82e-5 {:l4 1 :PI -1}]
    [1.77e-5 {:l3 2 :l4 -2}]
    [-1.67e-5 {:l4 2 :om4 -1 :psi -1}]
    [1.67e-5 {:om4 -1 :psi 1}]
    [-1.55e-5 {:l4 2 :PI -2 :G -2}]
    [1.42e-5 {:l4 2 :psi -2}]
    [1.05e-5 {:l1 1 :l4 -1}]
    [9.2e-6 {:l2 1 :l4 -1}]
    [-8.9e-6 {:l4 1 :PI -1 :G -1}]
    [-6.2e-6 {:l4 1 :pi4 1 :PI -2 :G -3}]
    [4.8e-6 {:l4 2 :om4 -2}]]])

(def ^:private mean-radius [5.90569 9.39657 14.98832 26.36273])

(defn- series [f rows args]
  (reduce (fn [s [k arg]]
            (+ s (* k (f (reduce-kv (fn [a sym m] (+ a (* m (if (= sym :deg) c/degrees (args sym)))))
                                    0.0 arg)))))
          0.0 rows))

(def ^:private perspective
  "Meeus's K: the satellites' distances in Jovian radii scaled so that a
  body in front of the disk is displaced by the planet's curvature."
  [17295.0 21819.0 27558.0 36548.0])

(defn- ecliptic-orbits
  "The E5 theory evaluated for the satellites as they were at `jd-light`
  (the Julian Day the theory is run for, less any light time), and turned
  onto the ecliptic of `mjd-tt`: `{:xyz [[x y z] x4 + the pole] :radii}`,
  in Jupiter radii about Jupiter's center. The fifth vector is the unit
  vector along Jupiter's pole, which the view needs to find the disk's
  orientation."
  [jd-light mjd-tt]
  (let [jd (+ mjd-tt c/jd-mjd-offset)
        t  (- jd-light 2443000.5)
        lin (fn [a b] (deg (+ a (* b t))))
        args0 {:l1 (lin 106.07719 203.48895579) :l2 (lin 175.73161 101.374724735)
               :l3 (lin 120.55883 50.317609207) :l4 (lin 84.44459 21.571071177)
               :pi1 (lin 97.0881 0.16138586) :pi2 (lin 154.8663 0.04726307)
               :pi3 (lin 188.184 0.00712734) :pi4 (lin 335.2868 0.00184)
               :om1 (lin 312.3346 -0.13279386) :om2 (lin 100.4411 -0.03263064)
               :om3 (lin 119.1942 -0.00717703) :om4 (lin 322.6186 -0.00175934)
               :phi (lin 199.6766 0.1737919) :psi (lin 316.5182 -0.00000208)
               :G' (lin 31.97853 0.0334597339) :PI (deg 13.469942)
               :G (+ (lin 30.23756 0.0830925701)
                     (* (deg 0.33033) (math/sin (deg (+ 163.679 (* 0.0010512 t)))))
                     (* (deg 0.03439) (math/sin (deg (- 34.486 (* 0.0161731 t))))))}
        S  (mapv #(* c/degrees (series math/sin (sigma %) args0)) [:sigma1 :sigma2 :sigma3 :sigma4])
        L  (mapv #(+ (args0 %1) %2) [:l1 :l2 :l3 :l4] S)
        args (merge args0 (zipmap [:L1 :L2 :L3 :L4] L) (zipmap [:S1 :S2 :S3 :S4] S))
        B  (mapv #(math/atan (series math/sin % args)) latitude)
        Rs (mapv #(* %1 (+ 1.0 (series math/cos %2 args))) mean-radius radius)
        ;; precession since B1950, to refer the longitudes to the equinox of date
        T0 (/ (- jd 2433282.423) 36525.0)
        P  (* (+ (deg 1.3966626) (* (deg 0.0003088) T0)) T0)
        psi (+ (:psi args) P)
        I  (deg (+ 3.120262 (* 0.0006 (/ (- jd 2415020.0) 36525.0))))
        xyz (conj (mapv (fn [l b r] [(* r (math/cos (- (+ l P) psi)) (math/cos b))
                                     (* r (math/sin (- (+ l P) psi)) (math/cos b))
                                     (* r (math/sin b))])
                        L B Rs)
                  [0.0 0.0 1.0])            ; a fictitious fifth body on the pole
        {:keys [raan i]} (orbits/mean-elements :jupiter mjd-tt)
        rot (fn [a b th] [(- (* a (math/cos th)) (* b (math/sin th)))
                          (+ (* a (math/sin th)) (* b (math/cos th)))])]
    {:radii Rs
     :xyz (mapv (fn [[X Y Z]]
                  (let [a X
                        b (- (* Y (math/cos I)) (* Z (math/sin I)))
                        c (+ (* Y (math/sin I)) (* Z (math/cos I)))
                        [a b] (rot a b (- psi raan))
                        [b c] (rot b c i)
                        [a b] (rot a b raan)]
                    [a b c]))
                xyz)}))

(defn positions-3d
  "Where Io, Europa, Ganymede and Callisto are at `mjd-tt`: `[[x y z]
  ...]` in Jupiter's equatorial radii from its center, in EME2000 --
  geometric, with no light time, for a model rather than an eyepiece.
  The same E5 evaluation as `positions`, stopped before the projection
  onto the sky."
  [mjd-tt]
  (let [to-J2000 (vsop87/ecliptic-of-date->J2000 mjd-tt)
        {:keys [xyz]} (ecliptic-orbits (+ mjd-tt c/jd-mjd-offset) mjd-tt)]
    (mapv #(lin/mat-vec to-J2000 %) (pop xyz))))

(def radii
  "The satellites' mean radii, km: the means of the triaxial radii in
  NAIF's pck00011.tpc."
  {:io 1821.5 :europa 1560.8 :ganymede 2631.2 :callisto 2410.3})

(def names [:io :europa :ganymede :callisto])

(defn positions
  "`[[x y] ...]` of Io, Europa, Ganymede and Callisto at `mjd-tt`, in
  Jupiter radii, by the E5 theory."
  [mjd-tt]
  (let [jd (+ mjd-tt c/jd-mjd-offset)
        ;; Jupiter seen from the Earth, light time included
        [s b-sun R] (solar/geometric mjd-tt)
        earth [(* R (math/cos b-sun) (math/cos s)) (* R (math/cos b-sun) (math/sin s)) (* R (math/sin b-sun))]
        [[x y z] delta tau]
        (loop [delta 5.0 i 0]
          (let [tau (elliptic/light-time delta)
                [l b r] (vsop87/heliocentric :jupiter (- mjd-tt tau))
                xyz (mapv + [(* r (math/cos b) (math/cos l)) (* r (math/cos b) (math/sin l)) (* r (math/sin b))]
                          earth)
                delta' (math/sqrt (reduce + (map * xyz xyz)))]
            (if (or (< (abs (- delta' delta)) 1e-9) (> i 10))
              [xyz delta' (elliptic/light-time delta')]
              (recur delta' (inc i)))))
        lam0 (math/atan2 y x)
        bet0 (math/atan (/ z (math/hypot x y)))
        {Rs :radii ecl :xyz} (ecliptic-orbits (- jd tau) mjd-tt)
        view (mapv (fn [[a b c]]
                     (let [a' (- (* a (math/sin lam0)) (* b (math/cos lam0)))
                           b' (+ (* a (math/cos lam0)) (* b (math/sin lam0)))]
                       [a' (+ (* c (math/sin bet0)) (* b' (math/cos bet0)))
                        (- (* c (math/cos bet0)) (* b' (math/sin bet0)))]))
                   ecl)
        [A5 _ C5] (peek view)
        D (math/atan2 A5 C5)]
    (mapv (fn [[A B C] r k]
            (let [x (- (* A (math/cos D)) (* C (math/sin D)))
                  y (+ (* A (math/sin D)) (* C (math/cos D)))
                  q (/ x r)
                  x (+ x (* (/ (abs B) k) (math/sqrt (- 1.0 (* q q)))))
                  W (/ delta (+ delta (/ B 2095.0)))]
              [(* x W) (* y W)]))
          (pop view) Rs perspective)))
