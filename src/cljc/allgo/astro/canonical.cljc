(ns allgo.astro.canonical
  "Canonical units (Vallado, section 3.8.1; Bate, Mueller and White,
  chapter 1): lengths in a body's radius, or for the Sun the astronomical
  unit, and times in the unit that makes its gravitational parameter one
  -- TU = sqrt(DU^3/mu). A circular orbit one DU out then has speed one
  DU/TU and period 2 pi TU, and the equations of motion lose their
  constants, which is why the older texts compute in them.

  For the Earth (JGM-3's mu and radius) a TU is 806.8 s and a DU/TU 7.905
  km/s; for the Sun a TU is 1/k days, 58.13, k the Gaussian constant."
  (:require [allgo.astro.constants :as c]
            [clojure.math :as math]))

(defn units
  "The canonical units of `body` (`:earth` or `:sun`), or of any body of
  parameter `mu` and reference length `du`: `{:du :tu :velocity :mu}`,
  km, s, km/s and km^3/s^2."
  ([body]
   (case body
     :earth (units c/GM-earth c/R-earth)
     :sun (units c/GM-sun c/AU)))
  ([mu du]
   (let [tu (math/sqrt (/ (* du du du) mu))]
     {:du du :tu tu :velocity (/ du tu) :mu mu})))

(defn ->canonical
  "A position (km), velocity (km/s) or time (s), by `kind` `:length`
  `:velocity` or `:time`, in the canonical units `u`; a vector is scaled
  component by component."
  [u kind x]
  (let [scale (/ 1.0 (case kind :length (:du u) :velocity (:velocity u) :time (:tu u)))]
    (if (number? x) (* x scale) (mapv #(* % scale) x))))

(defn <-canonical
  "The inverse of `->canonical`: back to km, km/s or s."
  [u kind x]
  (let [scale (case kind :length (:du u) :velocity (:velocity u) :time (:tu u))]
    (if (number? x) (* x scale) (mapv #(* % scale) x))))
