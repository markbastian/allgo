(ns allgo.astro.geopotential
  "Earth gravity as a spherical harmonic expansion (Montenbruck & Gill 3.2).

  The Earth is not a point, and outside it the departure is described by a
  series in (R/r)^n with coefficients C_nm and S_nm. Degree 2 order 0 alone
  -- the equatorial bulge, J2 -- is a thousand times every other term put
  together, and is what makes orbit planes precess.

  Evaluation follows M&G's recursion in the quantities

    V_nm = (R/r)^(n+1) P_nm(sin(lat)) cos(m lon)
    W_nm = (R/r)^(n+1) P_nm(sin(lat)) sin(m lon)

  which are built directly from Cartesian coordinates. No latitude,
  longitude or Legendre polynomial is ever formed: the recursion is over
  x, y, z alone, which avoids both the trigonometry and the singularity at
  the poles that a spherical formulation carries."
  (:require [allgo.astro.constants :as c]
            [clojure.math :as math]))

;; ------------------------------------------------------------ normalisation

(defn- factorial-ratio
  "(n-m)! / (n+m)!, built as a product so it stays finite at high degree
  where the factorials themselves overflow."
  [n m]
  (reduce (fn [acc k] (/ acc (double k))) 1.0 (range (inc (- n m)) (inc (+ n m)))))

(defn normalisation-factor
  "The factor taking a normalised coefficient to an unnormalised one.

  Published Earth models are normalised, because raw coefficients span so
  many orders of magnitude that the series is awkward to tabulate; the
  recursion below wants them raw."
  [n m]
  (math/sqrt (* (factorial-ratio n m)
                (+ (* 2.0 n) 1.0)
                (if (zero? m) 1.0 2.0))))

(defn denormalise
  "Convert a model's normalised coefficients to the raw form the recursion
  uses. A model already raw is returned unchanged."
  [{:keys [normalised? C S] :as model}]
  (if-not normalised?
    model
    (letfn [(scale [m] (reduce-kv (fn [acc [n mm] v]
                                    (assoc acc [n mm] (* v (normalisation-factor n mm))))
                                  {} m))]
      (assoc model :C (scale C) :S (scale S) :normalised? false))))

;; ------------------------------------------------------------- the recursion

(defn- vw
  "V and W to degree `deg`, as maps keyed by [n m].

  Built column by column: each order's diagonal term comes from the one
  before it, then the column is filled downward in degree."
  [[x y z] R deg]
  (let [r2  (+ (* x x) (* y y) (* z z))
        rho (/ (* R R) r2)
        x0  (/ (* R x) r2)
        y0  (/ (* R y) r2)
        z0  (/ (* R z) r2)]
    (loop [m 0
           V {[0 0] (/ R (math/sqrt r2))}
           W {[0 0] 0.0}]
      (if (> m deg)
        [V W]
        (let [;; diagonal: V_mm from V_(m-1)(m-1)
              [V W] (if (zero? m)
                      [V W]
                      (let [k (- (* 2.0 m) 1.0)
                            vp (V [(dec m) (dec m)])
                            wp (W [(dec m) (dec m)])]
                        [(assoc V [m m] (* k (- (* x0 vp) (* y0 wp))))
                         (assoc W [m m] (* k (+ (* x0 wp) (* y0 vp))))]))
              ;; then downward in degree at fixed order
              [V W] (loop [n (inc m) V V W W]
                      (if (> n deg)
                        [V W]
                        (let [a (/ (- (* 2.0 n) 1.0) (- n m))
                              b (/ (+ n m -1.0) (- n m))
                              v1 (V [(dec n) m] 0.0)
                              w1 (W [(dec n) m] 0.0)
                              v2 (V [(- n 2) m] 0.0)
                              w2 (W [(- n 2) m] 0.0)]
                          (recur (inc n)
                                 (assoc V [n m] (- (* a z0 v1) (* b rho v2)))
                                 (assoc W [n m] (- (* a z0 w1) (* b rho w2)))))))]
          (recur (inc m) V W))))))

;; ------------------------------------------------------------ the potential

(defn potential
  "Gravitational potential U at Earth-fixed `r`, km^2/s^2.

  Kept alongside the acceleration because the two are independently derived
  -- the acceleration is not differentiated from this -- so agreement
  between it and the numerical gradient of this is a real check on both."
  [{:keys [GM R] :as model} r degree]
  (let [{:keys [C S]} (denormalise model)
        [V W] (vw r R degree)]
    (* (/ GM R)
       (reduce + (for [n (range 0 (inc degree))
                       m (range 0 (inc n))]
                   (+ (* (get C [n m] (if (and (zero? n) (zero? m)) 1.0 0.0)) (V [n m] 0.0))
                      (* (get S [n m] 0.0) (W [n m] 0.0))))))))

;; --------------------------------------------------------- the acceleration

(defn acceleration
  "Acceleration from the geopotential at Earth-fixed position `r`, km/s^2,
  in the same Earth-fixed frame (M&G eq. 3.33).

  The recursion is carried one degree past the field so that every term can
  reach the next degree up: the gradient of a term of degree n is expressed
  in V and W of degree n+1."
  [model r degree]
  (let [{:keys [GM R C S]} (denormalise model)
        [V W] (vw r R (+ degree 1))
        gr    (/ GM (* R R))]
    (reduce
     (fn [[ax ay az] [n m]]
       (let [cnm (get C [n m] (if (and (zero? n) (zero? m)) 1.0 0.0))
             snm (get S [n m] 0.0)
             v+  (fn [i j] (V [i j] 0.0))
             w+  (fn [i j] (W [i j] 0.0))]
         (if (zero? m)
           [(+ ax (* gr (- (* cnm (v+ (inc n) 1)))))
            (+ ay (* gr (- (* cnm (w+ (inc n) 1)))))
            (+ az (* gr (* (- (+ n 1.0)) (* cnm (v+ (inc n) 0)))))]
           (let [fact (* (- (+ n 2.0) m) (- (+ n 1.0) m))]
             ;; Both C and S are negated in the leading bracket. Taking the
             ;; difference instead leaves the sine term right by accident and
             ;; the cosine term wrong, which shows up only on the sectorial
             ;; and tesseral coefficients -- never on J2.
             [(+ ax (* gr 0.5 (+ (- (+ (* cnm (v+ (inc n) (inc m)))
                                       (* snm (w+ (inc n) (inc m)))))
                                 (* fact (+ (* cnm (v+ (inc n) (dec m)))
                                            (* snm (w+ (inc n) (dec m))))))))
              (+ ay (* gr 0.5 (+ (- (* snm (v+ (inc n) (inc m)))
                                    (* cnm (w+ (inc n) (inc m))))
                                 (* fact (- (* snm (v+ (inc n) (dec m)))
                                            (* cnm (w+ (inc n) (dec m))))))))
              (+ az (* gr (* (- (+ n 1.0) m)
                             (- (- (* cnm (v+ (inc n) m)))
                                (* snm (w+ (inc n) m))))))]))))
     [0.0 0.0 0.0]
     (for [n (range 0 (inc degree)) m (range 0 (inc n))] [n m]))))

;; ---------------------------------------------------------------- a model

(def J2
  "The dominant zonal coefficient, unnormalised. Everything else in the
  field is smaller by three orders of magnitude."
  1.0826359e-3)

(def earth
  "A low-degree Earth field. Normalised coefficients, as models are
  published; the zonal terms are J2, J3 and J4 converted, and the degree-2
  order-2 pair is the equatorial ellipticity."
  {:GM c/GM-earth
   :R  c/R-earth
   :normalised? true
   :C {[0 0] 1.0
       [2 0] (- (/ J2 (math/sqrt 5.0)))
       [2 2] 2.43926e-6
       [3 0] (- (/ -2.5324e-6 (math/sqrt 7.0)))
       [4 0] (- (/ -1.6193e-6 3.0))}
   :S {[2 2] -1.40027e-6}})

(defn point-mass
  "The zero-degree field: a point mass, for comparison."
  [GM R]
  {:GM GM :R R :normalised? true :C {[0 0] 1.0} :S {}})
