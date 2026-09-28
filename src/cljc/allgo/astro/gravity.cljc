(ns allgo.astro.gravity
  "Other ways of summing the geopotential's acceleration (Vallado,
  *Fundamentals of Astrodynamics and Applications*, chapter 8): the same
  field as `allgo.astro.geopotential`, whose Montenbruck & Gill recursion
  these are checked against, reached by different routes.

  The classic one differentiates the potential in the spherical
  coordinates it is written in -- radius, latitude, longitude -- and
  turns the result Cartesian by the chain rule, which divides by the
  distance from the axis and so fails at the poles. Pines (\"Uniform
  representation of the gravitational potential and its derivatives\",
  AIAA Journal 11, 1973) avoids that by never forming latitude or
  longitude: the potential is written in the direction cosines s, t, u
  and the derivatives of the Legendre polynomials, all finite everywhere.

  Lear's and Gottlieb's algorithms, also free of singularities, are here
  in the fully normalized forms of Eckman, Brown and Adamo, *Normalization
  and Implementation of Three Gravitational Acceleration Models*
  (NASA/TP-2016-218604), transcribed from the report's MATLAB listings;
  normalized, they reach high degree. Pines's and the spherical partials
  use the coefficients unnormalized, which keeps them to moderate degree.

  Models are `allgo.astro.geopotential`'s maps. Positions and
  accelerations are Earth-fixed, km and km/s^2."
  (:require [allgo.astro.geopotential :as geo]
            [allgo.numerics.special :as special]
            [clojure.math :as math]))

(defn- coefficient [C n m] (get C [n m] (if (and (zero? n) (zero? m)) 1.0 0.0)))

;; ------------------------------------------------------- spherical partials

(defn spherical-acceleration
  "The acceleration from the field `model` to degree `degree` at `r`,
  from the potential's partial derivatives in radius, latitude and
  longitude:

    dU/dr   = -mu/r^2 sum (n+1)(R/r)^n P_nm (C cos m lon + S sin m lon)
    dU/dlat =  mu/r sum (R/r)^n (P_n,m+1 - m tan lat P_nm)(C cos + S sin)
    dU/dlon =  mu/r sum (R/r)^n m P_nm (S cos m lon - C sin m lon)

  turned Cartesian by the derivatives of r, latitude and longitude in x,
  y, z. Undefined on the polar axis, where longitude is."
  [model [x y z] degree]
  (let [{:keys [GM R C S]} (geo/denormalize model)
        rho2 (+ (* x x) (* y y)) rho (math/sqrt rho2)
        r2 (+ rho2 (* z z)) r (math/sqrt r2)
        sl (/ z r) cl (/ rho r) tl (/ z rho)
        lon (math/atan2 y x)
        P (special/associated-legendre sl cl (inc degree))
        [dr dlat dlon]
        (reduce (fn [[dr dlat dlon] [n m]]
                  (let [cnm (coefficient C n m) snm (get S [n m] 0.0)
                        cm (math/cos (* m lon)) sm (math/sin (* m lon))
                        k (/ (* GM (math/pow (/ R r) n)) r)
                        trig (+ (* cnm cm) (* snm sm))]
                    [(- dr (* (/ k r) (+ n 1.0) (P [n m]) trig))
                     (+ dlat (* k (- (get P [n (inc m)] 0.0) (* m tl (P [n m]))) trig))
                     (+ dlon (* k m (P [n m]) (- (* snm cm) (* cnm sm))))]))
                [0.0 0.0 0.0]
                (for [n (range (inc degree)) m (range (inc n))] [n m]))
        radial (- (/ dr r) (/ (* z dlat) (* r2 rho)))]
    [(- (* radial x) (* (/ dlon rho2) y))
     (+ (* radial y) (* (/ dlon rho2) x))
     (+ (* (/ dr r) z) (* (/ rho r2) dlat))]))

;; ------------------------------------------------------------------ Pines

(defn pines-acceleration
  "The acceleration from the field `model` to degree `degree` at `r`, by
  Pines's formulation. With s, t, u the direction cosines, r_m + i i_m =
  (s + i t)^m, A_nm the derived Legendre functions of u and rho_n = (mu/r^2)
  (R/r)^n:

    a1 = sum rho_n m A_nm (C r_m-1 + S i_m-1)
    a2 = sum rho_n m A_nm (S r_m-1 - C i_m-1)
    a3 = sum rho_n A_n,m+1 (C r_m + S i_m)
    a4 = -sum rho_n A_n+1,m+1 (C r_m + S i_m)

  and the acceleration is (a1 + s a4, a2 + t a4, a3 + u a4) -- finite at
  the poles and everywhere else."
  [model [x y z] degree]
  (let [{:keys [GM R C S]} (geo/denormalize model)
        r (math/sqrt (+ (* x x) (* y y) (* z z)))
        s (/ x r) t (/ y r) u (/ z r)
        A (special/derived-legendre u (+ degree 2))
        a (fn [n m] (get A [n m] 0.0))
        ;; r_m and i_m, the real and imaginary parts of (s + i t)^m
        [rm im] (loop [m 1 rs [1.0] is [0.0]]
                  (if (> m (inc degree))
                    [rs is]
                    (recur (inc m)
                           (conj rs (- (* s (peek rs)) (* t (peek is))))
                           (conj is (+ (* s (peek is)) (* t (peek rs)))))))
        [a1 a2 a3 a4]
        (reduce (fn [[a1 a2 a3 a4] [n m]]
                  (let [cnm (coefficient C n m) snm (get S [n m] 0.0)
                        rho (* (/ GM (* r r)) (math/pow (/ R r) n))
                        D (+ (* cnm (rm m)) (* snm (im m)))
                        [E F] (if (zero? m)
                                [0.0 0.0]
                                [(+ (* cnm (rm (dec m))) (* snm (im (dec m))))
                                 (- (* snm (rm (dec m))) (* cnm (im (dec m))))])]
                    [(+ a1 (* rho m (a n m) E))
                     (+ a2 (* rho m (a n m) F))
                     (+ a3 (* rho (a n (inc m)) D))
                     (- a4 (* rho (a (inc n) (inc m)) D))]))
                [0.0 0.0 0.0 0.0]
                (for [n (range (inc degree)) m (range (inc n))] [n m]))]
    [(+ a1 (* s a4)) (+ a2 (* t a4)) (+ a3 (* u a4))]))

;; --------------------------------------------- Lear and Gottlieb, normalized

(defn- normalized
  "The model's C and S, fully normalized."
  [{:keys [normalized? C S]}]
  (if normalized?
    [C S]
    (let [f (fn [m] (reduce-kv (fn [acc [n k] v] (assoc acc [n k] (/ v (geo/normalization-factor n k)))) {} m))]
      [(f C) (f S)])))

(defn lear-acceleration
  "The acceleration from the field `model` to degree `degree` at `r` by
  Lear's algorithm, normalized (NASA/TP-2016-218604, learnorm.m): the
  normalized Legendre functions and their derivatives in the sine of the
  latitude, the order-m ones carried divided by the cosine so that nothing
  is singular at the poles; the radial, east and north components summed,
  then turned Cartesian."
  [model [x y z] degree]
  (let [{:keys [GM R]} model
        [C S] (normalized model)
        cnm (fn [n m] (get C [n m] 0.0)) snm (fn [n m] (get S [n m] 0.0))
        nmax degree
        sq math/sqrt
        norm1 (fn [n] (sq (/ (+ (* 2.0 n) 1.0) (- (* 2.0 n) 1.0))))
        norm2 (fn [n] (sq (/ (+ (* 2.0 n) 1.0) (- (* 2.0 n) 3.0))))
        norm11 (fn [n] (/ (sq (/ (+ (* 2.0 n) 1.0) (* 2.0 n))) (- (* 2.0 n) 1.0)))
        norm1m (fn [n m] (sq (/ (* (- n m) (+ (* 2.0 n) 1.0)) (* (+ n m) (- (* 2.0 n) 1.0)))))
        norm2m (fn [n m] (sq (/ (* (- n m) (- n m 1.0) (+ (* 2.0 n) 1.0))
                                (* (+ n m) (+ n m -1.0) (- (* 2.0 n) 3.0)))))
        e1 (+ (* x x) (* y y))
        r2 (+ e1 (* z z))
        absr (sq r2)
        r1 (sq e1)
        sphi (/ z absr)
        cphi (/ r1 absr)
        [sm1 cm1] (if (zero? r1) [0.0 1.0] [(/ y r1) (/ x r1)])
        root3 (sq 3.0) root5 (sq 5.0)
        ;; the recursions, in maps: sin and cos of m lon, (R/r)^n, the zonal
        ;; functions and their derivatives, and the tesseral ones
        init {:rb {1 (/ R absr) 2 (math/pow (/ R absr) 2)}
              :sm {1 sm1 2 (* 2.0 cm1 sm1)}
              :cm {1 cm1 2 (- (* 2.0 cm1 cm1) 1.0)}
              :pn {1 (* root3 sphi) 2 (/ (* root5 (- (* 3.0 sphi sphi) 1.0)) 2.0)}
              :ppn {1 root3 2 (* root5 3.0 sphi)}
              :pnm {[1 1] root3 [2 2] (/ (* root5 root3 cphi) 2.0) [2 1] (* root5 root3 sphi)}
              :ppnm {[1 1] (- (* root3 sphi)) [2 2] (- (* root3 root5 sphi cphi))
                     [2 1] (* root5 root3 (- 1.0 (* 2.0 sphi sphi)))}}
        diag (reduce (fn [{:keys [rb sm cm pn ppn pnm] :as t} n]
                       (let [e (- (* 2.0 n) 1.0)
                             pnn (* e cphi (norm11 n) (pnm [(dec n) (dec n)]))]
                         (-> t
                             (assoc-in [:rb n] (* (rb (dec n)) (rb 1)))
                             (assoc-in [:sm n] (- (* 2.0 cm1 (sm (dec n))) (sm (- n 2))))
                             (assoc-in [:cm n] (- (* 2.0 cm1 (cm (dec n))) (cm (- n 2))))
                             (assoc-in [:pn n] (/ (- (* e sphi (norm1 n) (pn (dec n)))
                                                     (* (dec n) (norm2 n) (pn (- n 2))))
                                                  n))
                             (assoc-in [:ppn n] (* (norm1 n) (+ (* sphi (ppn (dec n))) (* n (pn (dec n))))))
                             (assoc-in [:pnm [n n]] pnn)
                             (assoc-in [:ppnm [n n]] (* (- n) sphi pnn)))))
                     init (range 3 (inc nmax)))
        {:keys [rb sm cm pn ppn pnm ppnm]}
        (reduce (fn [t [n m]]
                  (let [pnm (:pnm t)
                        e1 (* (- (* 2.0 n) 1.0) sphi)
                        e2 (* (- n) sphi)
                        e3 (* (norm1m n m) (get pnm [(dec n) m] 0.0))
                        e4 (+ n m)
                        e5 (/ (- (* e1 e3) (* (dec e4) (norm2m n m) (get pnm [(- n 2) m] 0.0))) (- n m))]
                    (-> t (assoc-in [:pnm [n m]] e5) (assoc-in [:ppnm [n m]] (+ (* e2 e5) (* e4 e3))))))
                diag (for [n (range 3 (inc nmax)) m (range 1 n)] [n m]))
        [a1 a3] (reduce (fn [[a1 a3] n]
                          (let [e1 (* (cnm n 0) (rb n))]
                            [(- a1 (* (inc n) e1 (pn n))) (+ a3 (* e1 (ppn n)))]))
                        [-1.0 0.0] (range 2 (inc nmax)))
        a3 (* cphi a3)
        [t1 a2 t3] (reduce (fn [[t1 a2 t3] n]
                             (let [[e1 e2 e3]
                                   (reduce (fn [[e1 e2 e3] m]
                                             (let [ts (snm n m) tc (cnm n m) tsm (sm m) tcm (cm m)
                                                   tp (get pnm [n m] 0.0)
                                                   e4 (+ (* ts tsm) (* tc tcm))]
                                               [(+ e1 (* e4 tp))
                                                (+ e2 (* m (- (* ts tcm) (* tc tsm)) tp))
                                                (+ e3 (* e4 (get ppnm [n m] 0.0)))]))
                                           [0.0 0.0 0.0] (range 1 (inc n)))]
                               [(+ t1 (* (inc n) (rb n) e1)) (+ a2 (* (rb n) e2)) (+ t3 (* (rb n) e3))]))
                           [0.0 0.0 0.0] (range 2 (inc nmax)))
        e4 (/ GM r2)
        a1 (* e4 (- a1 (* cphi t1)))
        a2 (* e4 a2)
        a3 (* e4 (+ a3 t3))
        e5 (- (* a1 cphi) (* a3 sphi))]
    [(- (* e5 cm1) (* a2 sm1)) (+ (* e5 sm1) (* a2 cm1)) (+ (* a1 sphi) (* a3 cphi))]))

(defn gottlieb-acceleration
  "The acceleration from the field `model` to degree `degree` at `r` by
  Gottlieb's algorithm, normalized (NASA/TP-2016-218604, gottliebnorm.m;
  Gottlieb, NASA CR-188243, 1993): the normalized Legendre functions of
  the direction's z cosine and the real and imaginary parts of (x + iy)^m
  over r^m, gathered into four sums from which the acceleration follows
  directly -- no angles, no singularities."
  [model [x y z] degree]
  (let [{:keys [GM R]} model
        [C S] (normalized model)
        cnm (fn [n m] (get C [n m] 0.0)) snm (fn [n m] (get S [n m] 0.0))
        nax degree mx degree
        sq math/sqrt
        norm1 (fn [n] (sq (/ (+ (* 2.0 n) 1.0) (- (* 2.0 n) 1.0))))
        norm2 (fn [n] (sq (/ (+ (* 2.0 n) 1.0) (- (* 2.0 n) 3.0))))
        norm11 (fn [n] (/ (sq (/ (+ (* 2.0 n) 1.0) (* 2.0 n))) (- (* 2.0 n) 1.0)))
        normn10 (fn [n] (sq (/ (* (+ n 1.0) n) 2.0)))
        norm1m (fn [n m] (sq (/ (* (- n m) (+ (* 2.0 n) 1.0)) (* (+ n m) (- (* 2.0 n) 1.0)))))
        norm2m (fn [n m] (sq (/ (* (- n m) (- n m 1.0) (+ (* 2.0 n) 1.0))
                                (* (+ n m) (+ n m -1.0) (- (* 2.0 n) 3.0)))))
        normn1 (fn [n m] (sq (* (+ n m 1.0) (- n m))))
        r (sq (+ (* x x) (* y y) (* z z)))
        ri (/ 1.0 r)
        xor (* x ri) yor (* y ri) zor (* z ri)
        ep zor
        reor (* R ri)
        muor2 (* GM ri ri)
        p0 (reduce (fn [p n] (assoc p [n n] (* (norm11 n) (p [(dec n) (dec n)]) (- (* 2.0 n) 1.0))))
                   {[0 0] 1.0 [1 1] (sq 3.0) [1 0] (* (sq 3.0) ep)}
                   (range 2 (inc nax)))
        pg (fn [p n m] (get p [n m] 0.0))
        {:keys [sumh sumgm sumj sumk]}
        (loop [n 2 p p0 ctil {0 1.0 1 xor} stil {0 0.0 1 yor} reorn reor
               acc {:sumh 0.0 :sumgm 1.0 :sumj 0.0 :sumk 0.0}]
          (if (> n nax)
            acc
            (let [reorn (* reorn reor)
                  n2m1 (- (* 2.0 n) 1.0) nm1 (dec n) np1 (inc n)
                  p (assoc p [n (dec n)] (* (normn1 n (dec n)) ep (pg p n n)))
                  p (assoc p [n 0] (/ (- (* n2m1 ep (norm1 n) (pg p nm1 0)) (* nm1 (norm2 n) (pg p (- n 2) 0))) n))
                  p (assoc p [n 1] (/ (- (* n2m1 ep (norm1m n 1) (pg p nm1 1)) (* n (norm2m n 1) (pg p (- n 2) 1))) nm1))
                  sumhn (* (normn10 n) (pg p n 1) (cnm n 0))
                  sumgmn (* (pg p n 0) (cnm n 0) np1)
                  p (reduce (fn [p m] (assoc p [n m] (/ (- (* n2m1 ep (norm1m n m) (pg p nm1 m))
                                                           (* (+ nm1 m) (norm2m n m) (pg p (- n 2) m)))
                                                        (- n m))))
                            p (range 2 (- n 1)))
                  ctil (assoc ctil n (- (* (ctil 1) (ctil (dec n))) (* (stil 1) (stil (dec n)))))
                  stil (assoc stil n (+ (* (stil 1) (ctil (dec n))) (* (ctil 1) (stil (dec n)))))
                  lim (min n mx)
                  [sumhn sumgmn sumjn sumkn]
                  (reduce (fn [[h g j k] m]
                            (let [mxpnm (* m (pg p n m))
                                  bnmtil (+ (* (cnm n m) (ctil m)) (* (snm n m) (stil m)))
                                  bnmtm1 (+ (* (cnm n m) (ctil (dec m))) (* (snm n m) (stil (dec m))))
                                  anmtm1 (- (* (cnm n m) (stil (dec m))) (* (snm n m) (ctil (dec m))))]
                              [(+ h (* (normn1 n m) (pg p n (inc m)) bnmtil))
                               (+ g (* (+ n m 1.0) (pg p n m) bnmtil))
                               (+ j (* mxpnm bnmtm1))
                               (- k (* mxpnm anmtm1))]))
                          [sumhn sumgmn 0.0 0.0] (range 1 (inc lim)))]
              (recur (inc n) p ctil stil reorn
                     {:sumh (+ (:sumh acc) (* reorn sumhn))
                      :sumgm (+ (:sumgm acc) (* reorn sumgmn))
                      :sumj (+ (:sumj acc) (* reorn sumjn))
                      :sumk (+ (:sumk acc) (* reorn sumkn))}))))
        lambda (+ sumgm (* ep sumh))]
    [(- (* muor2 (- (* lambda xor) sumj)))
     (- (* muor2 (- (* lambda yor) sumk)))
     (- (* muor2 (- (* lambda zor) sumh)))]))
