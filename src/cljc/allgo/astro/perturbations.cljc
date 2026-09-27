(ns allgo.astro.perturbations
  "General perturbations: how the orbital elements drift under a small
  force, rather than where the force carries a satellite step by step
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapters 8
  and 9).

  Gauss's form of the variation of parameters takes any perturbing
  acceleration, resolved radially (R), along-track (S) and normal to the
  orbit (W), and gives the rate of each classical element. Averaged over
  an orbit, those rates are the secular drift -- here numerically, for
  any force, and in closed form for the ones that have one: J2, drag in
  an exponential atmosphere and a steady push such as radiation pressure.

  Elements are `{:a :e :i :raan :argp :M}` as `allgo.astro.kepler` has
  them; km, seconds, radians."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kepler]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

;; ------------------------------------------------ variation of parameters

(defn rsw
  "The components of `f` along the radial, along-track and normal
  directions of the state `[r v]`."
  [[r v] f]
  (let [R (v3/normalize r)
        W (v3/normalize (v3/cross r v))
        S (v3/cross W R)]
    [(v3/dot f R) (v3/dot f S) (v3/dot f W)]))

(defn gauss-rates
  "The rates of the classical elements under a perturbing acceleration
  whose radial, along-track and normal components are `[fr fs fw]`, for an
  orbit with elements `el` -- Gauss's variational equations:

    da/dt    = 2a^2/h (e sin nu fr + p/r fs)
    de/dt    = (p sin nu fr + ((p + r) cos nu + r e) fs) / h
    di/dt    = r cos u fw / h
    dOmega/dt = r sin u fw / (h sin i)
    domega/dt = (-p cos nu fr + (p + r) sin nu fs) / (h e) - cos i dOmega/dt
    dM/dt    = n + b ((p cos nu - 2 r e) fr - (p + r) sin nu fs) / (a h e)

  p the semilatus rectum, h the angular momentum per unit mass, u the
  argument of latitude and b the semi-minor axis. `:M` here is the rate
  less the mean motion, the perturbation's part alone."
  ([el f] (gauss-rates mu el f))
  ([mu {:keys [a e i argp nu M]} [fr fs fw]]
   (let [nu (or nu (kepler/mean->true M e))
         p (* a (- 1.0 (* e e)))
         h (math/sqrt (* mu p))
         r (/ p (+ 1.0 (* e (math/cos nu))))
         u (+ argp nu)
         sn (math/sin nu) cn (math/cos nu)
         b (* a (math/sqrt (- 1.0 (* e e))))
         draan (/ (* r (math/sin u) fw) (* h (math/sin i)))]
     {:a (* (/ (* 2.0 a a) h) (+ (* e sn fr) (* (/ p r) fs)))
      :e (/ (+ (* p sn fr) (* (+ (* (+ p r) cn) (* r e)) fs)) h)
      :i (/ (* r (math/cos u) fw) h)
      :raan draan
      :argp (- (/ (+ (* (- p) cn fr) (* (+ p r) sn fs)) (* h e)) (* (math/cos i) draan))
      :M (/ (* b (- (* (- (* p cn) (* 2.0 r e)) fr) (* (+ p r) sn fs))) (* a h e))})))

(defn averaged-rates
  "The secular rates of the elements under the acceleration `(accel r v)`,
  inertial km/s^2: Gauss's rates averaged over the mean anomaly with the
  orbit held fixed -- the first-order drift, for any force. `samples`
  points are taken round the orbit (the rates are periodic in M, so an
  even spacing converges fast)."
  ([el accel] (averaged-rates mu el accel 360))
  ([mu el accel samples]
   (let [rates (for [k (range samples)
                     :let [M (* 2.0 math/PI (/ k samples))
                           el' (assoc el :M M :nu nil)
                           [r v :as s] (kepler/elements->state mu el')]]
                 (gauss-rates mu el' (rsw s (accel r v))))]
     (into {} (for [k [:a :e :i :raan :argp :M]]
                [k (/ (reduce + (map k rates)) samples)])))))

;; ------------------------------------------------------------------- J2

(defn j2-secular
  "The first-order secular rates J2 gives an orbit of semi-major axis `a`,
  eccentricity `e` and inclination `i`: `{:raan :argp :M}`, rad/s, the
  last the whole mean-anomaly rate, mean motion included.

    dOmega/dt = -3/2 n J2 (R/p)^2 cos i
    domega/dt =  3/4 n J2 (R/p)^2 (4 - 5 sin^2 i)
    dM/dt     =  n (1 + 3/4 J2 (R/p)^2 sqrt(1 - e^2) (2 - 3 sin^2 i))

  The node regresses for prograde orbits; the perigee stands still at the
  critical inclination, 63.4 degrees, where 4 = 5 sin^2 i."
  ([a e i] (j2-secular mu c/R-earth geo/J2 a e i))
  ([mu R J2 a e i]
   (let [n (kepler/mean-motion mu a)
         p (* a (- 1.0 (* e e)))
         k (* J2 (/ (* R R) (* p p)))
         s2 (math/pow (math/sin i) 2)]
     {:raan (* -1.5 n k (math/cos i))
      :argp (* 0.75 n k (- 4.0 (* 5.0 s2)))
      :M (* n (+ 1.0 (* 0.75 k (math/sqrt (- 1.0 (* e e))) (- 2.0 (* 3.0 s2)))))})))

(defn sun-synchronous-inclination
  "The inclination at which J2 turns the node at `rate` -- by default once
  a year eastward, keeping the orbit's plane fixed with respect to the Sun
  -- for semi-major axis `a` and eccentricity `e`; nil where no inclination
  turns it that fast."
  ([a e] (sun-synchronous-inclination a e (/ (* 2.0 math/PI) (* 365.2421897 86400.0))))
  ([a e rate]
   (let [{node :raan} (j2-secular a e 0.0)
         ci (/ rate node)]
     (when (<= -1.0 ci 1.0) (math/acos ci)))))
