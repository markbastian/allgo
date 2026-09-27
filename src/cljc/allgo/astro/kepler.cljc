(ns allgo.astro.kepler
  "Two-body motion and orbital elements (Montenbruck & Gill chapter 2).

  Six numbers describe a Kepler orbit, and the choice of which six is the
  whole subject. A state vector -- position and velocity -- is what an
  integrator wants, but it says nothing at a glance: two states differing in
  the last decimal may be the same orbit at different times. The classical
  elements separate the parts that stay fixed under two-body motion (size,
  shape, orientation) from the one that does not (where the satellite is),
  which is why perturbations are described as element rates.

  The price is singularities. A circular orbit has no periapsis, so the
  argument of periapsis is undefined; an equatorial orbit has no node, so
  the right ascension is undefined. Both cases are real orbits and both are
  handled here by convention rather than by returning a NaN."
  (:require [allgo.astro.constants :as c]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

;; ------------------------------------------------------------ scalar results

(defn specific-energy
  "v^2/2 - mu/r. Negative for a bound orbit, zero for escape."
  [mu r v] (- (* 0.5 (v3/dot v v)) (/ mu (v3/length r))))

(defn vis-viva
  "Speed at radius `r` on an orbit of semi-major axis `a`:
  v^2 = mu (2/r - 1/a). One equation covering every conic -- for a
  hyperbola `a` is negative and the speed stays finite at infinity."
  [mu r a] (math/sqrt (* mu (- (/ 2.0 r) (/ 1.0 a)))))

(defn mean-motion [mu a] (math/sqrt (/ mu (* a a a))))

(defn period
  "Kepler's third law. Depends only on the semi-major axis: two orbits of
  wildly different shape but equal `a` take exactly as long."
  [mu a] (/ c/two-pi (mean-motion mu a)))

;; --------------------------------------------------------------- anomalies

(defn kepler-equation
  "Solve M = E - e sin E for the eccentric anomaly, by Newton's method.

  The starting guess matters at high eccentricity, where the equation is
  nearly flat near periapsis and a poor start can wander. M + e sin M is the
  standard choice and converges in a handful of steps out to e = 0.99."
  ([M e] (kepler-equation M e 1e-13))
  ([M e tol]
   (let [M (am/wrap-2pi M)]
     (loop [E (+ M (* e (math/sin M))) n 0]
       (let [f  (- E (* e (math/sin E)) M)
             fp (- 1.0 (* e (math/cos E)))
             dE (/ f fp)]
         (if (or (< (abs dE) tol) (>= n 60))
           (- E dE)
           (recur (- E dE) (inc n))))))))

(defn eccentric->true
  "The eccentric anomaly is measured on the circumscribing circle, the true
  anomaly at the focus. They agree at periapsis and apoapsis and nowhere
  else."
  [E e]
  (am/wrap-2pi (* 2.0 (math/atan2 (* (math/sqrt (+ 1.0 e)) (math/sin (* 0.5 E)))
                                  (* (math/sqrt (- 1.0 e)) (math/cos (* 0.5 E)))))))

(defn true->eccentric [nu e]
  (am/wrap-2pi (* 2.0 (math/atan2 (* (math/sqrt (- 1.0 e)) (math/sin (* 0.5 nu)))
                                  (* (math/sqrt (+ 1.0 e)) (math/cos (* 0.5 nu)))))))

(defn eccentric->mean [E e] (am/wrap-2pi (- E (* e (math/sin E)))))
(defn mean->true [M e] (eccentric->true (kepler-equation M e) e))
(defn true->mean [nu e] (eccentric->mean (true->eccentric nu e) e))

(defn kepler-bisection
  "Solve Kepler's equation by Sinnott's binary search (Meeus 30, third
  method): fifty-odd halvings of a quarter circle, each fixing one bit of
  the eccentric anomaly.

  Newton's method can fail near e = 1 and M = 0, where the curve is almost
  flat and a first step overshoots wildly; this cannot, for any
  eccentricity below one."
  [M e]
  (let [M  (am/wrap-2pi M)
        [M sign] (if (> M math/PI) [(- c/two-pi M) -1.0] [M 1.0])]
    (loop [E (* 0.5 math/PI) d (* 0.25 math/PI) i 0]
      (if (= i 53)
        (am/wrap-2pi (* sign E))
        (recur (if (neg? (- M (- E (* e (math/sin E))))) (- E d) (+ E d))
               (* 0.5 d) (inc i))))))

(defn kepler-approximate
  "The eccentric anomaly to first order in `e` (Meeus 30.8), tan E =
  sin M / (cos M - e): a hundredth of a degree for e = 0.1, and a start for
  something better."
  [M e]
  (am/wrap-2pi (math/atan2 (math/sin M) (- (math/cos M) e))))

;; --------------------------------------------- parabolic orbits (Meeus 34)

(defn parabolic
  "`[nu r]`, true anomaly and distance in AU, `t` days after perihelion on
  a parabola of perihelion distance `q` AU about the Sun -- Barker's
  equation, which unlike Kepler's has a closed-form solution, a cube root."
  [q t]
  (let [W (/ (* 3.0 c/gaussian-k t) (* (math/sqrt 2.0) q (math/sqrt q)))
        G (* 0.5 W)
        Y (math/cbrt (+ G (math/sqrt (+ (* G G) 1.0))))
        s (- Y (/ 1.0 Y))]
    [(* 2.0 (math/atan s)) (* q (+ 1.0 (* s s)))]))

;; ---------------------------------------- near-parabolic orbits (Meeus 35)

(defn near-parabolic
  "`[nu r]` `t` days after perihelion on an orbit of perihelion distance
  `q` AU and eccentricity `e` near one, by Landgraf's method, or nil if it
  does not converge.

  A comet with e = 0.99 is poorly served by both Kepler's equation, whose
  mean anomaly then crawls, and Barker's, which ignores the difference from
  a parabola. Landgraf expands about the parabola instead, in a series
  that converges fast for e near one on either side of it."
  [q e t]
  (if (zero? t)
    [0.0 q]
    (let [q1 (/ (* c/gaussian-k (math/sqrt (/ (+ 1.0 e) q))) (* 2.0 q))
          g  (/ (- 1.0 e) (+ 1.0 e))
          q2 (* q1 t)
          s0 (/ 2.0 (* 3.0 (abs q2)))
          s0 (/ 2.0 (math/tan (* 2.0 (math/atan (math/cbrt (math/tan (* 0.5 (math/atan s0))))))))
          s0 (if (neg? t) (- s0) s0)
          tol 1e-9
          ;; the series of Landgraf's line 42 onward, summed for one s
          q3-of (fn [s]
                  (let [y (* s s)]
                    (loop [z 1.0 g1 (- (* y s)) q3 (+ q2 (/ (* 2.0 g s y) 3.0))]
                      (let [z  (inc z)
                            g1 (* (- g1) g y)
                            f  (* g1 (/ (- z (* (inc z) g)) (inc (* 2.0 z))))
                            q3 (+ q3 f)]
                        (cond
                          (or (> z 50.0) (> (abs f) 10000.0)) nil
                          (<= (abs f) tol) q3
                          :else (recur z g1 q3))))))
          solve (fn [s q3]
                  (loop [s s i 0]
                    (let [s' (/ (+ (/ (* 2.0 s s s) 3.0) q3) (+ (* s s) 1.0))]
                      (if (or (<= (abs (- s' s)) tol) (> i 1000)) s' (recur s' (inc i))))))
          s  (if (= e 1.0)
               s0
               (loop [s s0 l 0]
                 (when (<= l 50)
                   (when-let [q3 (q3-of s)]
                     (let [s' (solve s q3)]
                       (if (<= (abs (- s' s)) tol) s' (recur s' (inc l))))))))]
      (when s
        (let [nu (* 2.0 (math/atan s))]
          [(am/wrap-2pi nu) (/ (* q (+ 1.0 e)) (+ 1.0 (* e (math/cos nu))))])))))

;; ------------------------------------------------ the ellipse (Meeus 33.d)

(defn ellipse-circumference
  "Circumference of an ellipse of semi-major axis `a` and eccentricity `e`.
  Ramanujan's approximation, exact for a circle and within 0.3 percent up
  to e = 0.97, unless `exact?`, when it sums the series in
  m = (a - b)/(a + b) to convergence."
  ([a e] (ellipse-circumference a e false))
  ([a e exact?]
   (let [b (* a (math/sqrt (- 1.0 (* e e))))]
     (if-not exact?
       (* math/PI (- (* 3.0 (+ a b)) (math/sqrt (* (+ a (* 3.0 b)) (+ (* 3.0 a) b)))))
       (let [m  (/ (- a b) (+ a b))
             m2 (* m m)]
         (loop [sum 1.0 term (* 0.25 m2) nf 1.0 df 4.0]
           (let [sum' (+ sum term)]
             (if (= sum' sum)
               (/ (* c/two-pi a sum) (+ 1.0 m))
               (recur sum' (/ (* term nf nf m2) (* df df)) (+ nf 2.0) (+ df 2.0))))))))))

;; ------------------------------------------------------ state and elements

(def ^:private circular-tol 1e-11)
(def ^:private equatorial-tol 1e-11)

(defn state->elements
  "Classical elements from position and velocity.

  Returns `{:a :e :i :raan :argp :nu :M}` -- semi-major axis, eccentricity,
  inclination, right ascension of the ascending node, argument of periapsis,
  true anomaly and mean anomaly. Angles in radians.

  Where an element is undefined the convention is to set it to zero and fold
  its meaning into the next one along: a circular orbit gets argp = 0 and
  measures the anomaly from the node, an equatorial one gets raan = 0 and
  measures from the x axis."
  [mu r v]
  (let [rm   (v3/length r)
        h    (v3/cross r v)
        hm   (v3/length h)
        node (v3/cross [0.0 0.0 1.0] h)
        nm   (v3/length node)
        evec (v3/scale (v3/sub (v3/scale r (- (v3/dot v v) (/ mu rm)))
                               (v3/scale v (v3/dot r v)))
                       (/ 1.0 mu))
        e    (v3/length evec)
        en   (specific-energy mu r v)
        a    (if (< (abs en) 1e-15) ##Inf (/ (- mu) (* 2.0 en)))
        i    (math/acos (am/clamp (/ (nth h 2) hm) -1.0 1.0))
        circular?   (< e circular-tol)
        equatorial? (< nm (* equatorial-tol hm))
        raan (if equatorial? 0.0 (am/wrap-2pi (math/atan2 (nth node 1) (nth node 0))))
        argp (cond
               circular?   0.0
               equatorial? (am/wrap-2pi (math/atan2 (nth evec 1) (nth evec 0)))
               :else       (let [ang (math/acos (am/clamp (/ (v3/dot node evec) (* nm e)) -1.0 1.0))]
                             (am/wrap-2pi (if (neg? (nth evec 2)) (- c/two-pi ang) ang))))
        nu   (cond
               ;; circular and equatorial: measure from x, the only reference left
               (and circular? equatorial?)
               (am/wrap-2pi (let [ang (math/atan2 (nth r 1) (nth r 0))]
                              (if (neg? (nth h 2)) (- ang) ang)))
               ;; circular: measure from the node -- argument of latitude
               circular?
               (let [ang (math/acos (am/clamp (/ (v3/dot node r) (* nm rm)) -1.0 1.0))]
                 (am/wrap-2pi (if (neg? (nth r 2)) (- c/two-pi ang) ang)))
               :else
               (let [ang (math/acos (am/clamp (/ (v3/dot evec r) (* e rm)) -1.0 1.0))]
                 (am/wrap-2pi (if (neg? (v3/dot r v)) (- c/two-pi ang) ang))))]
    {:a a :e e :i i :raan raan :argp argp :nu nu
     :M (if (< e 1.0) (true->mean nu e) ##NaN)}))

(defn elements->state
  "Position and velocity from classical elements. The inverse of
  `state->elements`, and exact: the two round-trip to machine precision."
  [mu {:keys [a e i raan argp nu M]}]
  (let [nu (or nu (mean->true M e))
        p  (* a (- 1.0 (* e e)))
        rm (/ p (+ 1.0 (* e (math/cos nu))))
        ;; in the perifocal frame, where periapsis lies on the x axis
        rp [(* rm (math/cos nu)) (* rm (math/sin nu)) 0.0]
        k  (math/sqrt (/ mu p))
        vp [(* k (- (math/sin nu))) (* k (+ e (math/cos nu))) 0.0]
        cO (math/cos raan) sO (math/sin raan)
        cw (math/cos argp) sw (math/sin argp)
        ci (math/cos i)    si (math/sin i)
        ;; R_z(raan) R_x(i) R_z(argp), written out
        m  [[(- (* cO cw) (* sO sw ci)) (- (- (* cO sw)) (* sO cw ci)) (* sO si)]
            [(+ (* sO cw) (* cO sw ci)) (- (* cO cw ci) (* sO sw))     (- (* cO si))]
            [(* sw si)                  (* cw si)                      ci]]
        apply-m (fn [x] (mapv (fn [row] (v3/dot row x)) m))]
    [(apply-m rp) (apply-m vp)]))

(defn propagate
  "Advance a state by `dt` seconds on its Kepler orbit, analytically.

  No integration: the elements are constant under two-body motion except the
  mean anomaly, which advances uniformly. That is the whole point of the
  element set, and it is exact -- which makes it the reference an integrator
  can be judged against."
  [mu r v dt]
  (let [{:keys [a] :as el} (state->elements mu r v)
        M' (+ (:M el) (* (mean-motion mu a) dt))]
    (elements->state mu (assoc el :M (am/wrap-2pi M') :nu nil))))

;; ------------------------------------------------------ non-singular form

(defn equinoctial
  "Equinoctial elements, which have no singularity at zero eccentricity or
  zero inclination.

  Eccentricity and the argument of periapsis are folded into a vector
  (h, k) that simply goes to zero on a circular orbit rather than leaving an
  angle undefined; inclination and node likewise into (p, q). This is the
  form perturbation work and orbit determination use, because a filter
  cannot estimate an angle that does not exist."
  [{:keys [a e i raan argp nu]}]
  (let [lam (am/wrap-2pi (+ raan argp nu))]
    {:a a
     :h (* e (math/sin (+ raan argp)))
     :k (* e (math/cos (+ raan argp)))
     :p (* (math/tan (* 0.5 i)) (math/sin raan))
     :q (* (math/tan (* 0.5 i)) (math/cos raan))
     :lambda lam}))

(defn from-equinoctial
  "Classical elements back from equinoctial form."
  [{:keys [a h k p q lambda]}]
  (let [e    (math/sqrt (+ (* h h) (* k k)))
        i    (* 2.0 (math/atan (math/sqrt (+ (* p p) (* q q)))))
        raan (if (and (zero? p) (zero? q)) 0.0 (am/wrap-2pi (math/atan2 p q)))
        argp (am/wrap-2pi (- (if (and (zero? h) (zero? k)) 0.0 (math/atan2 h k)) raan))]
    {:a a :e e :i i :raan raan :argp argp
     :nu (am/wrap-2pi (- lambda raan argp))}))
