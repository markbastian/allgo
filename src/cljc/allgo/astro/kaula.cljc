(ns allgo.astro.kaula
  "Kaula's form of the geopotential: each spherical harmonic term written
  in the orbital elements (Kaula, *Theory of Satellite Geodesy*, 1966,
  chapter 3; Vallado, chapter 9),

    V_lm = (mu/a) (R/a)^l  sum_p F_lmp(i)  sum_q G_lpq(e) S_lmpq(w, M, W, theta),

    S_lmpq = C_lm cos psi + S_lm sin psi    (l - m even)
             -S_lm cos psi + C_lm sin psi   (l - m odd)
    psi = (l - 2p) w + (l - 2p + q) M + m (W - theta)

  with the inclination functions F and the eccentricity functions G --
  the latter Hansen coefficients, X^(-(l+1), l-2p)_(l-2p+q)(e). A term is
  resonant where its argument psi stands still, which is what makes the
  expansion the tool for repeat ground tracks, geosynchronous drift and
  orbital resonance: `argument-rate` gives psi's rate.

  Coefficients unnormalized, as `allgo.astro.geopotential` denormalizes
  them; angles radians."
  (:require [allgo.numerics.quadrature :as quadrature]
            [clojure.math :as math]))

(defn- fact [n] (reduce * 1.0 (range 1 (inc n))))

(defn- binomial [n k]
  (if (or (neg? k) (> k n)) 0.0 (/ (fact n) (* (fact k) (fact (- n k))))))

(defn inclination-function
  "F_lmp(i), Kaula's inclination function (his 3.61):

    sum_t (2l - 2t)! / (t! (l - t)! (l - m - 2t)! 2^(2l - 2t)) sin^(l-m-2t) i
      sum_s C(m, s) cos^s i  sum_c C(l - m - 2t + s, c) C(m - s, p - t - c) (-1)^(c - k)

  t from 0 to the lesser of p and k = floor((l - m)/2)."
  [l m p i]
  (let [k (quot (- l m) 2)
        si (math/sin i) ci (math/cos i)]
    (reduce + (for [t (range (inc (min p k)))
                    :let [a (/ (fact (- (* 2 l) (* 2 t)))
                               (* (fact t) (fact (- l t)) (fact (- l m (* 2 t))) (math/pow 2.0 (- (* 2 l) (* 2 t)))))]
                    s (range (inc m))
                    c (range (inc (+ (- l m (* 2 t)) s)))
                    :let [b (binomial (- m s) (- p t c))]
                    :when (pos? b)]
                (* a (math/pow si (- l m (* 2 t))) (binomial m s) (math/pow ci s)
                   (binomial (+ (- l m (* 2 t)) s) c) b (if (even? (- c k)) 1.0 -1.0))))))

(defn eccentricity-function
  "G_lpq(e), Kaula's eccentricity function: the Hansen coefficient
  X^(-(l+1), l-2p)_(l-2p+q)(e), from its definition --

    (1/2 pi) integral over M of (a/r)^(l+1) cos((l - 2p) f - (l - 2p + q) M) dM

  -- by the trapezoid rule, exact to rounding for a periodic integrand
  once the points resolve its harmonics."
  [l p q e]
  (let [j (- l (* 2 p)) k (+ j q)
        n (+ 64 (* 16 (abs k)) (long (/ 64.0 (max 1e-3 (- 1.0 e)))))
        integrand (fn [M]
                    (let [E (loop [E (if (< e 0.8) M math/PI) i 0]
                              (let [d (/ (- E (* e (math/sin E)) M) (- 1.0 (* e (math/cos E))))]
                                (if (or (< (abs d) 1e-15) (> i 50)) (- E d) (recur (- E d) (inc i)))))
                          r-over-a (- 1.0 (* e (math/cos E)))
                          f (* 2.0 (math/atan2 (* (math/sqrt (+ 1.0 e)) (math/sin (* 0.5 E)))
                                               (* (math/sqrt (- 1.0 e)) (math/cos (* 0.5 E)))))]
                      (* (math/pow r-over-a (- (inc l))) (math/cos (- (* j f) (* k M))))))]
    (/ (quadrature/trapezoid integrand 0.0 (* 2.0 math/PI) n) (* 2.0 math/PI))))

(defn argument
  "psi_lmpq = (l - 2p) w + (l - 2p + q) M + m (W - theta), theta the
  Greenwich sidereal angle."
  [l m p q {:keys [argp M raan]} theta]
  (+ (* (- l (* 2 p)) argp) (* (+ (- l (* 2 p)) q) M) (* m (- raan theta))))

(defn argument-rate
  "psi's rate, rad/s, from the rates of the elements and of sidereal
  time `theta-dot`: the term is resonant where it vanishes."
  [l m p q {:keys [argp M raan]} theta-dot]
  (+ (* (- l (* 2 p)) argp) (* (+ (- l (* 2 p)) q) M) (* m (- raan theta-dot))))

(defn term
  "V_lm, the potential per unit mass (km^2/s^2) of the harmonic of degree
  `l`, order `m` and unnormalized coefficients `C` `S`, at elements `el`
  `{:a :e :i :raan :argp :M}` and sidereal angle `theta`: Kaula's sum over
  p, and over q from -`q-max` to `q-max` (by default 10 + 100 e, enough
  for ten digits to e = 0.2 and degree 5, the series in q converging
  more slowly as the degree rises; for a circular orbit only q = 0
  contributes)."
  ([mu R l m C S el theta] (term mu R l m C S el theta (+ 10 (long (math/ceil (* 100.0 (:e el)))))))
  ([mu R l m C S {:keys [a e i] :as el} theta q-max]
   (let [qs (if (zero? e) [0] (range (- q-max) (inc q-max)))]
     (* (/ mu a) (math/pow (/ R a) l)
        (reduce + (for [p (range (inc l))
                        :let [f (inclination-function l m p i)]
                        q qs
                        :let [g (if (zero? e) (if (zero? q) 1.0 0.0) (eccentricity-function l p q e))
                              psi (argument l m p q el theta)]]
                    (* f g (if (even? (- l m))
                             (+ (* C (math/cos psi)) (* S (math/sin psi)))
                             (+ (* (- S) (math/cos psi)) (* C (math/sin psi)))))))))))
