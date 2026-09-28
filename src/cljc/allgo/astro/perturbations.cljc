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
            [allgo.astro.srp :as srp]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.quadrature :as quadrature]
            [allgo.numerics.special :as special]
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

;; ------------------------------------------------------------------ drag

(defn drag-secular
  "The secular decay drag gives an orbit of semi-major axis `a` and
  eccentricity `e` in an atmosphere whose density falls exponentially from
  `rho-p` kg/m^3 at perigee with scale height `H` km, for a ballistic
  coefficient `B` = C_D A / m, m^2/kg, and no rotation of the air:
  King-Hele's first-order results,

    da per revolution = -2 pi B rho_p a^2 exp(-c) (I0 + 2e I1)
    de per revolution = -2 pi B rho_p a exp(-c) (I1 + e/2 (I0 + I2))

  with c = ae/H and I_k the modified Bessel functions of c -- the
  density's rise at each perigee pass, integrated round the orbit.
  Returns `{:da :de}` per revolution and `:a-rate` `:e-rate` per second.
  Good while e is small; the terms dropped are e^2 against those kept."
  ([a e B rho-p H] (drag-secular mu a e B rho-p H))
  ([mu a e B rho-p H]
   (let [c (/ (* a e) H)
         k (* 2.0 math/PI B rho-p 1e3 (math/exp (- c)))
         i0 (special/bessel-i 0 c) i1 (special/bessel-i 1 c) i2 (special/bessel-i 2 c)
         da (- (* k a a (+ i0 (* 2.0 e i1))))
         de (- (* k a (+ i1 (* 0.5 e (+ i0 i2)))))
         period (kepler/period mu a)]
     {:da da :de de :a-rate (/ da period) :e-rate (/ de period)})))

(defn circular-decay-rate
  "da/dt, km/s, of a circular orbit of radius `a` in air of density `rho`
  kg/m^3, for a ballistic coefficient `B`: -B rho sqrt(mu a), the along-
  track drag put through Gauss's equation for a."
  ([a B rho] (circular-decay-rate mu a B rho))
  ([mu a B rho] (- (* B rho 1e3 (math/sqrt (* mu a))))))

(defn circular-lifetime
  "Seconds for a circular orbit to decay from radius `a0` to `a1`, the
  density at each height (km above the Earth's radius) `(density h)`
  kg/m^3: the decay rate integrated, by Simpson's rule in the radius."
  ([a0 a1 B density] (circular-lifetime mu a0 a1 B density 2000))
  ([mu a0 a1 B density n]
   (quadrature/simpson (fn [a] (/ -1.0 (circular-decay-rate mu a B (density (- a c/R-earth))))) a1 a0 n)))

;; ------------------------------------------------- a steady push: SRP

(defn eccentricity-vector
  "e, pointing to periapsis with the eccentricity for length: (v x h)/mu
  - r/|r|."
  ([s] (eccentricity-vector mu s))
  ([mu [r v]]
   (v3/sub (v3/scale (v3/cross v (v3/cross r v)) (/ 1.0 mu)) (v3/normalize r))))

(defn steady-push-secular
  "The secular rates of the angular momentum h = r x v and the eccentricity
  vector e of the orbit through state `s` under an acceleration `F` that
  holds steady over the orbit -- radiation pressure out of the Earth's
  shadow, over the day or so the Sun takes to move a degree:

    <dh/dt> = -3/2 a (e x F)
    <de/dt> =  3/2 sqrt(p/mu) (F x h/|h|)

  From dh/dt = r x F and de/dt = (F x h + v x (r x F))/mu, averaged over
  the orbit: the mean position is -3/2 a e, and the mean of r v^T, whose
  symmetric part is a derivative and averages away, leaves 1/2 F x h. The
  semi-major axis has no secular change, as the mean velocity is zero.
  Returns `{:h :e}`."
  ([s F] (steady-push-secular mu s F))
  ([mu [r v :as s] F]
   (let [h (v3/cross r v)
         e (eccentricity-vector mu s)
         a (/ 1.0 (- (/ 2.0 (v3/length r)) (/ (v3/dot v v) mu)))
         p (/ (v3/dot h h) mu)]
     {:h (v3/scale (v3/cross e F) (* -1.5 a))
      :e (v3/scale (v3/cross F (v3/normalize h)) (* 1.5 (math/sqrt (/ p mu))))})))

(defn srp-secular
  "`steady-push-secular` with the push radiation pressure gives a
  satellite of `area-to-mass` m^2/kg and reflectivity `cr` with the Sun at
  `r-sun` -- the Earth's shadow, which interrupts it, left out."
  ([s r-sun area-to-mass cr] (srp-secular mu s r-sun area-to-mass cr))
  ([mu s r-sun area-to-mass cr]
   (steady-push-secular mu s (srp/acceleration [0.0 0.0 0.0] r-sun area-to-mass cr 1.0))))

;; ------------------------------------------------ J2's short-period terms

(defn j2-osculating
  "The osculating state `[r v]` of an orbit with mean elements `el` --
  Brouwer's, as two-line elements carry them -- under J2: Brouwer's
  first-order short-period terms in Lyddane's form, as Spacetrack Report
  No. 3 (Hoots and Roehrich, 1980) gives them for SGP4, with k2 = J2 R^2/2
  and theta = cos i:

    r     = r_L (1 - 3/2 (k2/p^2) sqrt(1-e^2) (3 theta^2 - 1)) + 1/2 (k2/p) (1 - theta^2) cos 2u
    u     = u_L - 1/4 (k2/p^2) (7 theta^2 - 1) sin 2u
    Omega = Omega + 3/2 (k2/p^2) theta sin 2u
    i     = i + 3/2 (k2/p^2) theta sin i cos 2u
    r'    = r'_L - (k2/p) n (1 - theta^2) sin 2u
    r nu' = r nu'_L + (k2/p) n ((1 - theta^2) cos 2u + 3/2 (3 theta^2 - 1))

  L marking the two-body values of the mean orbit. First order in J2 and,
  as the report has them for SGP4, without the terms of order J2 e: against
  J2's motion integrated numerically, a circular low orbit drifts along
  track by J2^2, some 45 m an orbit, and one of eccentricity 0.01 by J2 e,
  some 450 m."
  ([el] (j2-osculating mu c/R-earth geo/J2 el))
  ([mu R J2 {:keys [a e i raan argp M]}]
   (let [n (math/sqrt (/ mu (* a a a)))
         p (* a (- 1.0 (* e e)))
         E (kepler/kepler-equation M e)
         rl (* a (- 1.0 (* e (math/cos E))))
         rdotl (/ (* (math/sqrt (* mu a)) e (math/sin E)) rl)
         rvdotl (/ (math/sqrt (* mu p)) rl)
         nu (kepler/eccentric->true E e)
         u (+ argp nu)
         k2 (* 0.5 J2 R R)
         kp (/ k2 p) kp2 (/ k2 (* p p))
         th (math/cos i) th2 (* th th)
         s2u (math/sin (* 2.0 u)) c2u (math/cos (* 2.0 u))
         r (+ (* rl (- 1.0 (* 1.5 kp2 (math/sqrt (- 1.0 (* e e))) (- (* 3.0 th2) 1.0))))
              (* 0.5 kp (- 1.0 th2) c2u))
         uk (- u (* 0.25 kp2 (- (* 7.0 th2) 1.0) s2u))
         node (+ raan (* 1.5 kp2 th s2u))
         inc (+ i (* 1.5 kp2 th (math/sin i) c2u))
         rdot (- rdotl (* kp n (- 1.0 th2) s2u))
         rvdot (+ rvdotl (* kp n (+ (* (- 1.0 th2) c2u) (* 1.5 (- (* 3.0 th2) 1.0)))))
         ;; the orientation vectors
         sn (math/sin node) cn (math/cos node) si (math/sin inc) ci (math/cos inc)
         su (math/sin uk) cu (math/cos uk)
         mhat [(* (- sn) ci) (* cn ci) si]
         nhat [cn sn 0.0]
         U (v3/add (v3/scale mhat su) (v3/scale nhat cu))
         V (v3/sub (v3/scale mhat cu) (v3/scale nhat su))]
     [(v3/scale U r) (v3/add (v3/scale U rdot) (v3/scale V rvdot))])))

(defn j2-propagate
  "The osculating state `dt` seconds after an epoch where the mean elements
  are `el`: the mean elements carried forward by J2's secular rates, then
  the short-period terms of `j2-osculating` put back."
  ([el dt] (j2-propagate mu c/R-earth geo/J2 el dt))
  ([mu R J2 {:keys [a e i raan argp M] :as el} dt]
   (let [rates (j2-secular mu R J2 a e i)]
     (j2-osculating mu R J2 (assoc el
                                   :raan (+ raan (* (:raan rates) dt))
                                   :argp (+ argp (* (:argp rates) dt))
                                   :M (+ M (* (:M rates) dt)))))))

(defn osculating->mean
  "The mean elements whose `j2-osculating` state is `s`, by fixed-point
  iteration on the elements: the osculating elements are the first guess,
  and each step moves the guess by how far its own osculating elements
  miss the target's. J2 is small, so a few steps settle it."
  ([s] (osculating->mean mu c/R-earth geo/J2 s))
  ([mu R J2 [r v]]
   (let [target (kepler/state->elements mu r v)
         ks [:a :e :i :raan :argp :M]
         angle? #{:i :raan :argp :M}
         d (fn [x y k] (let [dx (- x y)] (if (angle? k) (am/wrap-angle dx) dx)))]
     (loop [guess (select-keys target ks) k 0]
       (let [[r' v'] (j2-osculating mu R J2 guess)
             got (kepler/state->elements mu r' v')
             step (into {} (for [key ks] [key (d (target key) (got key) key)]))
             guess' (merge-with + guess step)]
         (if (or (> k 50) (< (abs (:a step)) 1e-9))
           guess'
           (recur guess' (inc k))))))))
