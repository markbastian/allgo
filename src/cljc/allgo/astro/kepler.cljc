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
            [clojure.math :as math]))

(defn- wrap-2pi [x] (let [r (rem x c/two-pi)] (if (neg? r) (+ r c/two-pi) r)))

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
   (let [M (wrap-2pi M)]
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
  (wrap-2pi (* 2.0 (math/atan2 (* (math/sqrt (+ 1.0 e)) (math/sin (* 0.5 E)))
                               (* (math/sqrt (- 1.0 e)) (math/cos (* 0.5 E)))))))

(defn true->eccentric [nu e]
  (wrap-2pi (* 2.0 (math/atan2 (* (math/sqrt (- 1.0 e)) (math/sin (* 0.5 nu)))
                               (* (math/sqrt (+ 1.0 e)) (math/cos (* 0.5 nu)))))))

(defn eccentric->mean [E e] (wrap-2pi (- E (* e (math/sin E)))))
(defn mean->true [M e] (eccentric->true (kepler-equation M e) e))
(defn true->mean [nu e] (eccentric->mean (true->eccentric nu e) e))

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
        evec (v3/scale (mapv - (v3/scale r (- (v3/dot v v) (/ mu rm)))
                             (v3/scale v (v3/dot r v)))
                       (/ 1.0 mu))
        e    (v3/length evec)
        en   (specific-energy mu r v)
        a    (if (< (abs en) 1e-15) ##Inf (/ (- mu) (* 2.0 en)))
        i    (math/acos (max -1.0 (min 1.0 (/ (nth h 2) hm))))
        circular?   (< e circular-tol)
        equatorial? (< nm (* equatorial-tol hm))
        raan (if equatorial? 0.0 (wrap-2pi (math/atan2 (nth node 1) (nth node 0))))
        argp (cond
               circular?   0.0
               equatorial? (wrap-2pi (math/atan2 (nth evec 1) (nth evec 0)))
               :else       (let [ang (math/acos (max -1.0 (min 1.0 (/ (v3/dot node evec) (* nm e)))))]
                             (wrap-2pi (if (neg? (nth evec 2)) (- c/two-pi ang) ang))))
        nu   (cond
               ;; circular and equatorial: measure from x, the only reference left
               (and circular? equatorial?)
               (wrap-2pi (let [ang (math/atan2 (nth r 1) (nth r 0))]
                           (if (neg? (nth h 2)) (- ang) ang)))
               ;; circular: measure from the node -- argument of latitude
               circular?
               (let [ang (math/acos (max -1.0 (min 1.0 (/ (v3/dot node r) (* nm rm)))))]
                 (wrap-2pi (if (neg? (nth r 2)) (- c/two-pi ang) ang)))
               :else
               (let [ang (math/acos (max -1.0 (min 1.0 (/ (v3/dot evec r) (* e rm)))))]
                 (wrap-2pi (if (neg? (v3/dot r v)) (- c/two-pi ang) ang))))]
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
    (elements->state mu (assoc el :M (wrap-2pi M') :nu nil))))

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
  (let [lam (wrap-2pi (+ raan argp nu))]
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
        raan (if (and (zero? p) (zero? q)) 0.0 (wrap-2pi (math/atan2 p q)))
        argp (wrap-2pi (- (if (and (zero? h) (zero? k)) 0.0 (math/atan2 h k)) raan))]
    {:a a :e e :i i :raan raan :argp argp
     :nu (wrap-2pi (- lambda raan argp))}))
