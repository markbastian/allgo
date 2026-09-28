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

  Models are `allgo.astro.geopotential`'s maps; coefficients are used
  unnormalized, which keeps these to moderate degree. Positions and
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
