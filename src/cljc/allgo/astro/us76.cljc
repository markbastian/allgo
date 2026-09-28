(ns allgo.astro.us76
  "The U.S. Standard Atmosphere, 1976 (NOAA, NASA and USAF; NOAA-S/T
  76-1562, NASA-TM-X-74335; a work of the US government), and the
  piecewise-exponential atmosphere built on it for drag (Vallado,
  *Fundamentals of Astrodynamics and Applications*, chapter 8).

  Below 86 km the standard is defined analytically -- hydrostatic
  equilibrium through seven layers of constant lapse rate in geopotential
  altitude -- and computed here from its defining constants. Above, its
  densities come from integrating the diffusion of each constituent, and
  are taken here from the standard's Table I at the altitudes the
  exponential model needs.

  The exponential model gives each layer from a base altitude h0 the
  density rho0 exp(-(h - h0)/H), the scale height H chosen so the layers
  meet: continuous, and exact at every base. Altitudes are geometric, km;
  densities kg/m^3."
  (:require [clojure.math :as math]))

(def ^:private r0 6356.766)                                 ; km, the effective Earth radius
(def ^:private g0 9.80665)                                  ; m/s^2
(def ^:private R* 8314.32)                                  ; J/(kmol K)
(def ^:private M0 28.9644)                                  ; kg/kmol

(def ^:private layers
  "Base geopotential altitude (km'), base molecular-scale temperature (K)
  and lapse rate (K/km') of each layer to 86 km."
  (let [bases [[0.0 -6.5] [11.0 0.0] [20.0 1.0] [32.0 2.8] [47.0 0.0] [51.0 -2.8] [71.0 -2.0]]]
    (loop [[[h l] & more] bases T 288.15 P 101325.0 out []]
      (if-not h
        out
        (let [out (conj out {:h h :T T :P P :L l})]
          (if-let [[h' _] (first more)]
            (let [dh (- h' h)
                  T' (+ T (* l dh))
                  P' (if (zero? l)
                       (* P (math/exp (/ (* -1.0 g0 M0 dh 1000.0) (* R* T))))
                       (* P (math/pow (/ T T') (/ (* g0 M0) (* R* (/ l 1000.0))))))]
              (recur more T' P' out))
            out))))))

(defn geopotential-altitude
  "Geopotential altitude, km', of geometric altitude `z` km."
  [z]
  (/ (* r0 z) (+ r0 z)))

(defn lower
  "`{:T :P :rho}` -- molecular-scale temperature (K), pressure (Pa) and
  density (kg/m^3) -- at geometric altitude `z` km, 0 to 86, from the
  standard's defining equations."
  [z]
  (let [H (geopotential-altitude z)
        {:keys [h T P L]} (last (take-while #(<= (:h %) H) layers))
        dh (- H h)
        TM (+ T (* L dh))
        p (if (zero? L)
            (* P (math/exp (/ (* -1.0 g0 M0 dh 1000.0) (* R* T))))
            (* P (math/pow (/ T TM) (/ (* g0 M0) (* R* (/ L 1000.0))))))]
    {:T TM :P p :rho (/ (* p M0) (* R* TM))}))

(def upper-densities
  "Densities, kg/m^3, from the standard's Table I (geometric altitude,
  metric units, pages 68-73) at the altitudes above 86 km where the
  exponential model has a base."
  {90.0 3.416e-6 100.0 5.604e-7 110.0 9.708e-8 120.0 2.222e-8 130.0 8.152e-9
   140.0 3.831e-9 150.0 2.076e-9 180.0 5.194e-10 200.0 2.541e-10 250.0 6.073e-11
   300.0 1.916e-11 350.0 7.014e-12 400.0 2.803e-12 450.0 1.184e-12 500.0 5.215e-13
   600.0 1.137e-13 700.0 3.070e-14 800.0 1.136e-14 900.0 5.759e-15 1000.0 3.561e-15})

(def base-altitudes
  "The exponential model's base altitudes, km, as the standard piecewise
  model (Wertz; Vallado) places them."
  [0.0 25.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0
   180.0 200.0 250.0 300.0 350.0 400.0 450.0 500.0 600.0 700.0 800.0 900.0 1000.0])

(def exponential-table
  "`[[h0 rho0 H] ...]`: each base altitude, the standard's density there,
  and the scale height, km, that carries it to the next base's density --
  the last layer keeping the one below it."
  (let [rho (fn [h] (or (upper-densities h) (:rho (lower h))))
        rows (mapv (fn [h] [h (rho h)]) base-altitudes)
        heights (mapv (fn [[h r] [h' r']] (/ (- h' h) (math/log (/ r r')))) rows (rest rows))]
    (mapv (fn [[h r] H] [h r H]) rows (conj heights (peek heights)))))

(defn exponential-density
  "Density, kg/m^3, at geometric altitude `h` km by the piecewise-
  exponential model; above 1000 km, the last layer's scale height
  carries on."
  [h]
  (let [[h0 r0 H] (or (last (take-while #(<= (first %) h) exponential-table)) (first exponential-table))]
    (* r0 (math/exp (- (/ (- h h0) H))))))
