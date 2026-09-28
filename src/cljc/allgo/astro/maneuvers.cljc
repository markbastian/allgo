(ns allgo.astro.maneuvers
  "Changing orbits: transfers, plane changes, rendezvous and relative
  motion (Vallado, *Fundamentals of Astrodynamics and Applications*,
  chapter 6).

  Every impulsive maneuver here is the same arithmetic: the velocity the
  spacecraft has on one orbit, the velocity it needs on the next, at the
  one point the two share, and the difference between them. What varies
  is the choice of orbits. Hohmann's is the cheapest two-burn transfer
  between circles; the bi-elliptic beats it when the ratio of radii
  passes 11.94, by going out far beyond the target where speed is cheap;
  a one-tangent burn buys time with fuel. A plane change costs 2 v sin(di/2)
  wherever it is done, so it is done where v is least, at apoapsis, and
  best folded into a burn that happens there anyway.

  Radii in km, speeds in km/s, times in seconds, angles in radians. `mu`
  defaults to the Earth's."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

(defn speed
  "Speed at radius `r` on a conic of semi-major axis `a` (vis-viva)."
  ([r a] (speed mu r a))
  ([mu r a] (math/sqrt (- (/ (* 2.0 mu) r) (/ mu a)))))

(defn- semi-major-axis
  "a of an orbit of eccentricity `e` that passes radius `r` at true
  anomaly `nu`."
  [r e nu]
  (/ (* r (+ 1.0 (* e (math/cos nu)))) (- 1.0 (* e e))))

(defn delta-v
  "The burn that turns velocity `v1` into `v2` when the two are `theta`
  apart: the law of cosines."
  [v1 v2 theta]
  (math/sqrt (- (+ (* v1 v1) (* v2 v2)) (* 2.0 v1 v2 (math/cos theta)))))

;; ------------------------------------------------------------ transfers

(defn hohmann
  "The Hohmann transfer from radius `rinit` to `rfinal`, each on an orbit
  of eccentricity `einit` / `efinal` where the burn happens at true anomaly
  `nuinit` / `nufinal` (0 and pi for the classic circle-to-circle case):
  `{:dva :dvb :tof}`."
  ([rinit rfinal] (hohmann rinit rfinal 0.0 0.0 0.0 math/PI))
  ([rinit rfinal einit efinal nuinit nufinal]
   (let [atran (* 0.5 (+ rinit rfinal))
         ainit (semi-major-axis rinit einit nuinit)
         afinal (semi-major-axis rfinal efinal nufinal)]
     {:dva (abs (- (speed rinit atran) (speed rinit ainit)))
      :dvb (abs (- (speed rfinal afinal) (speed rfinal atran)))
      :tof (* math/PI (math/sqrt (/ (* atran atran atran) mu)))})))

(defn bi-elliptic
  "The bi-elliptic transfer from `rinit` out to `rb` and back in to
  `rfinal`: `{:dva :dvb :dvc :tof}`. Three burns, and cheaper than
  Hohmann's two when rfinal/rinit is large enough and rb larger still."
  ([rinit rb rfinal] (bi-elliptic rinit rb rfinal 0.0 0.0 0.0 math/PI))
  ([rinit rb rfinal einit efinal nuinit nufinal]
   (let [a1 (* 0.5 (+ rinit rb))
         a2 (* 0.5 (+ rb rfinal))
         ainit (semi-major-axis rinit einit nuinit)
         afinal (semi-major-axis rfinal efinal nufinal)
         t (fn [a] (* math/PI (math/sqrt (/ (* a a a) mu))))]
     {:dva (abs (- (speed rinit a1) (speed rinit ainit)))
      :dvb (abs (- (speed rb a2) (speed rb a1)))
      :dvc (abs (- (speed rfinal afinal) (speed rfinal a2)))
      :tof (+ (t a1) (t a2))})))

(defn one-tangent
  "A transfer tangent to the first orbit only, reaching `rfinal` at
  transfer true anomaly `nutran` rather than at apoapsis: faster than
  Hohmann's and dearer. `nuinit` 0 starts from periapsis of the transfer,
  pi from apoapsis. Returns `{:dva :dvb :tof :e :a :va :vb}` of the
  transfer, or nil where no conic through both points is tangent there."
  [rinit rfinal efinal nuinit nutran]
  (let [R (/ rinit rfinal)
        periapsis? (< (abs nuinit) 1e-9)
        e (if periapsis?
            (/ (- R 1.0) (- (math/cos nutran) R))
            (/ (- R 1.0) (+ (math/cos nutran) R)))
        a (if periapsis? (/ rinit (- 1.0 e)) (/ rinit (+ 1.0 e)))]
    (when (and (>= e 0.0) (< e 1.0) (pos? a))
      (let [ainit rinit
            afinal (semi-major-axis rfinal efinal nutran)
            va (speed rinit a)
            vb (speed rfinal a)
            vf (speed rfinal afinal)
            fpa-tran (math/atan (/ (* e (math/sin nutran)) (+ 1.0 (* e (math/cos nutran)))))
            fpa-final (math/atan (/ (* efinal (math/sin nutran)) (+ 1.0 (* efinal (math/cos nutran)))))
            E (math/acos (/ (+ e (math/cos nutran)) (+ 1.0 (* e (math/cos nutran)))))
            tof (* (math/sqrt (/ (* a a a) mu)) (- E (* e (math/sin E))))]
        {:dva (abs (- va (speed rinit ainit)))
         :dvb (delta-v vb vf (- fpa-tran fpa-final))
         :tof (if periapsis? tof (- (* math/PI (math/sqrt (/ (* a a a) mu))) tof))
         :e e :a a :va va :vb vb}))))

;; --------------------------------------------------------- plane changes

(defn inclination-only
  "The burn that turns the orbit plane by `di` about the radius vector at a
  point where the speed is `v` and the flight-path angle `fpa`: only the
  horizontal component turns."
  [di v fpa]
  (* 2.0 v (math/cos fpa) (math/sin (* 0.5 di))))

(defn node-only
  "Moving the ascending node by `draan` at speed `v`, leaving the
  inclination `i` alone on a circular orbit: `{:dv :u-init :u-final :i}`,
  the burn made where the two planes cross, at argument of latitude
  `:u-init` on the old orbit and `:u-final` on the new.

  On an eccentric orbit (`e` > 0) the burn is made at the orbit's highest
  latitude instead (u = 90 degrees), where turning the velocity also tilts
  the plane: the new inclination comes out larger, and is returned."
  [i e draan v fpa]
  (if (< e 1e-9)
    (let [theta (math/acos (+ (* (math/cos i) (math/cos i))
                              (* (math/sin i) (math/sin i) (math/cos draan))))
          u (math/asin (/ (* (math/sin i) (math/sin draan)) (math/sin theta)))]
      {:dv (inclination-only theta v fpa) :u-init (- math/PI u) :u-final u :i i})
    (let [i' (math/atan (/ (math/tan i) (math/cos draan)))
          theta (math/acos (+ (* (math/cos i) (math/cos i'))
                              (* (math/sin i) (math/sin i') (math/cos draan))))]
      {:dv (inclination-only theta v fpa) :u-init (/ math/PI 2)
       :u-final (math/asin (/ (math/sin i) (math/sin i')))
       :i i'})))

(defn inclination-and-node
  "Changing inclination from `iinit` to `ifinal` and the node by `draan` in
  one burn at speed `v` and flight-path angle `fpa`: `{:dv :u-init
  :u-final}`, the angle between the planes found by spherical trigonometry
  and the burn made where they cross."
  [iinit ifinal draan v fpa]
  (let [ci (math/cos iinit) si (math/sin iinit)
        cf (math/cos ifinal) sf (math/sin ifinal)
        theta (math/acos (+ (* ci cf) (* si sf (math/cos draan))))
        st (math/sin theta)]
    ;; the triangle of the two nodes and the crossing, by the cosine rule
    ;; for angles
    {:dv (inclination-only theta v fpa)
     :u-init (math/acos (/ (- (* ci (math/cos theta)) cf) (* si st)))
     :u-final (math/acos (/ (- ci (* cf (math/cos theta))) (* sf st)))}))

(defn combined
  "A Hohmann transfer from `rinit` to `rfinal` that also changes the
  inclination by `di`, split between its two burns to minimize the total:
  `{:di1 :di2 :dva :dvb :tof}`. The split is Vallado's closed form, s =
  atan(sin di / (R^1.5 + cos di)) / di with R the ratio of radii -- most
  of the turn at apoapsis, a little at periapsis, where it is nearly free
  because it just rotates a burn that is being made anyway.

  With `:optimal true` the split is instead found by minimizing the total
  numerically, which the closed form approximates."
  ([rinit rfinal di] (combined rinit rfinal di {}))
  ([rinit rfinal di {:keys [einit efinal nuinit nufinal optimal]
                     :or {einit 0.0 efinal 0.0 nuinit 0.0 nufinal math/PI}}]
   (let [atran (* 0.5 (+ rinit rfinal))
         vi (speed rinit (semi-major-axis rinit einit nuinit))
         vta (speed rinit atran)
         vf (speed rfinal (semi-major-axis rfinal efinal nufinal))
         vtb (speed rfinal atran)
         total (fn [d1] (+ (delta-v vi vta d1) (delta-v vtb vf (- di d1))))
         di1 (if optimal
               ;; golden-section search over the first burn's share
               (roots/minimize total 0.0 di {:tol 0.0 :max-iter 100})
               (let [R (/ rfinal rinit)]
                 (math/atan (/ (math/sin di) (+ (math/pow R 1.5) (math/cos di))))))]
     {:di1 di1 :di2 (- di di1)
      :dva (delta-v vi vta di1) :dvb (delta-v vtb vf (- di di1))
      :tof (* math/PI (math/sqrt (/ (* atran atran atran) mu)))})))

;; ------------------------------------------------------------ rendezvous

(defn- mean-motion [a] (math/sqrt (/ mu (* a a a))))

(defn rendezvous-same-orbit
  "Catching a target that leads by `lead` radians on the same circular
  orbit of radius `a`, by dropping into a phasing orbit for `k-int` of its
  revolutions while the target makes `k-tgt`: `{:a-phase :tof :dv}` --
  `:dv` the total of the two equal burns, signed as the first (negative
  to drop into a lower, faster orbit and catch up).

  Vallado's phase angle is the other way round, the interceptor's angle
  less the target's: his -20 degrees is a lead of 20 here."
  [a lead k-int k-tgt]
  (let [w (mean-motion a)
        tau (/ (- (* 2.0 math/PI k-tgt) lead) (* k-int w))
        a-phase (math/cbrt (* mu (math/pow (/ tau (* 2.0 math/PI)) 2.0)))]
    (when (> a-phase (* 0.5 a))
      {:a-phase a-phase :tof (* k-int tau)
       :dv (* 2.0 (- (speed a a-phase) (speed a a)))})))

(defn rendezvous-coplanar
  "Catching a target on another circular orbit in the same plane: wait
  until it leads by the angle a Hohmann transfer needs, then transfer.
  `lead` is the target's lead now; the wait is the first that works, plus
  `k` whole synodic periods. Returns `{:phase-final :wait :tof :dv}` -- the
  lead the transfer needs, the wait for it, the transfer time and the
  total burn."
  [r-int r-tgt lead k]
  (let [w-int (mean-motion r-int) w-tgt (mean-motion r-tgt)
        {:keys [dva dvb tof]} (hohmann r-int r-tgt)
        dw (- w-int w-tgt)
        needed (- math/PI (* w-tgt tof))
        ;; the lead falls at dw per second; wait until it is the one needed
        wait (/ (+ (mod (- lead needed) (* 2.0 math/PI)) (* 2.0 math/PI k)) dw)]
    {:phase-final needed :wait (if (neg? dw) (- wait) wait) :tof tof :dv (+ dva dvb)}))

(defn rendezvous-noncoplanar
  "Catching a target on a circular orbit of radius `a-tgt` from one of
  radius `a-int` inclined `di` to it. Angles are arguments of latitude
  measured from the line where the two planes cross: the interceptor at
  `u-int`, the target at `u-tgt`.

  The plan is Vallado's: coast to the node line (u = 0), enter a phasing
  orbit with its periapsis there for `k-int` revolutions while the target
  makes the phase good -- `k-tgt` extra revolutions allowed -- then a
  Hohmann transfer whose apoapsis falls on the other node, where the plane
  change is folded into the circularizing burn. Returns `{:t-coast
  :t-phase :t-trans :a-phase :dv-phase :dv-trans1 :dv-trans2}`."
  [a-int a-tgt di u-int u-tgt k-int k-tgt]
  (let [two-pi (* 2.0 math/PI)
        w-int (mean-motion a-int) w-tgt (mean-motion a-tgt)
        a-tr (* 0.5 (+ a-int a-tgt))
        t-trans (* math/PI (math/sqrt (/ (* a-tr a-tr a-tr) mu)))
        t-coast (/ (mod (- u-int) two-pi) w-int)
        u-then (+ u-tgt (* w-tgt t-coast))
        ;; where the target must be when the transfer starts, to arrive at
        ;; the far node with it
        u-needed (- math/PI (* w-tgt t-trans))
        t-phase (/ (+ (mod (- u-needed u-then) two-pi) (* two-pi k-tgt)) w-tgt)
        a-phase (math/cbrt (* mu (math/pow (/ t-phase (* two-pi k-int)) 2.0)))
        v-ph (speed a-int a-phase)]
    {:t-coast t-coast :t-phase t-phase :t-trans t-trans :a-phase a-phase
     :dv-phase (abs (- v-ph (speed a-int a-int)))
     :dv-trans1 (abs (- (speed a-int a-tr) v-ph))
     :dv-trans2 (delta-v (speed a-tgt a-tr) (speed a-tgt a-tgt) di)}))

;; ------------------------------------------------------------ low thrust

(defn low-thrust
  "Edelbaum's continuous low-thrust spiral between circular orbits of
  radii `rinit` and `rfinal` with an inclination change `di`: the total
  delta-v, sqrt(v1^2 - 2 v1 v2 cos(pi di / 2) + v2^2); and, for a thruster
  giving acceleration `accel` km/s^2 at the start while losing mass at the
  fractional rate `mdot` (negative, 1/s), the time it takes -- the rocket
  equation with that acceleration growing as the mass goes."
  [rinit rfinal di mdot accel]
  (let [v1 (math/sqrt (/ mu rinit)) v2 (math/sqrt (/ mu rfinal))
        dv (math/sqrt (+ (* v1 v1) (* v2 v2) (* -2.0 v1 v2 (math/cos (* 0.5 math/PI di)))))]
    {:dv dv :tof (/ (- (math/exp (/ (* dv mdot) accel)) 1.0) mdot)}))

;; ----------------------------------------------- Hill / Clohessy-Wiltshire

(defn hill
  "Relative motion near a target on a circular orbit of radius `a`: the
  position and velocity `[r v]` after `t` seconds of a chaser that starts
  at `r0` `v0` relative to it. Axes: x radial, y along-track, z
  cross-track. Linearized -- good while the separation is small against
  `a` -- and in closed form: the Clohessy-Wiltshire solution."
  [a [x0 y0 z0] [dx0 dy0 dz0] t]
  (let [n (mean-motion a)
        nt (* n t) s (math/sin nt) cs (math/cos nt)]
    [[(+ (* (- 4.0 (* 3.0 cs)) x0) (* (/ s n) dx0) (* (/ 2.0 n) (- 1.0 cs) dy0))
      (+ (* 6.0 (- s nt) x0) y0 (* (/ -2.0 n) (- 1.0 cs) dx0) (* (/ (- (* 4.0 s) (* 3.0 nt)) n) dy0))
      (+ (* cs z0) (* (/ s n) dz0))]
     [(+ (* 3.0 n s x0) (* cs dx0) (* 2.0 s dy0))
      (+ (* -6.0 n (- 1.0 cs) x0) (* -2.0 s dx0) (* (- (* 4.0 cs) 3.0) dy0))
      (+ (* (- n) s z0) (* cs dz0))]]))

(defn hill-intercept
  "The initial relative velocity that carries a chaser from `r0` to the
  target itself in `t` seconds, by Clohessy-Wiltshire."
  [a [x0 y0 z0] t]
  (let [n (mean-motion a)
        nt (* n t) s (math/sin nt) cs (math/cos nt)
        den (- (* 8.0 (- 1.0 cs)) (* 3.0 nt s))
        dy0 (/ (* n (+ (* x0 (- (* 6.0 nt s) (* 14.0 (- 1.0 cs)))) (* -1.0 y0 s))) den)
        dx0 (/ (* n (+ (* x0 (- 4.0 (* 3.0 cs))) (* 2.0 (- 1.0 cs) (/ dy0 n)))) (- s))
        dz0 (/ (* (- n) z0 cs) s)]
    [dx0 dy0 dz0]))
