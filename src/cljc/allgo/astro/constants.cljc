(ns allgo.astro.constants
  "Physical and astronomical constants for the force models of Montenbruck &
  Gill chapter 3.

  Units throughout the astro package are kilometres, seconds and kilograms,
  which is what the gravitational parameters below are quoted in and what
  keeps orbital radii near 1e4 rather than 1e7. Angles are radians."
  (:require [clojure.math :as math]))

;; ------------------------------------------------------------ mathematical

(def pi math/PI)
(def two-pi (* 2.0 math/PI))
(def arcsec (/ math/PI 180.0 3600.0))
(def degrees (/ math/PI 180.0))

;; --------------------------------------------------------------- the Earth

(def GM-earth
  "Geocentric gravitational constant, km^3/s^2 (JGM-3)."
  398600.4415)

(def R-earth
  "Equatorial radius, km (WGS-84 / JGM-3)."
  6378.1363)

(def flattening
  "Earth's flattening, WGS-84."
  (/ 1.0 298.257223563))

(def omega-earth
  "Earth's rotation rate, rad/s. This is the sidereal rate, not 2*pi/86400 --
  a sidereal day is shorter than a solar one, and the difference is what
  makes the atmosphere co-rotate correctly under a satellite."
  7.2921158553e-5)

;; ------------------------------------------------------- the other bodies

(def GM-sun  1.32712440018e11)
(def GM-moon 4902.801)

(def AU
  "Astronomical unit, km (IAU 2012)."
  149597870.7)

(def R-sun  696000.0)
(def R-moon 1738.0)

;; --------------------------------------------------------------- radiation

(def solar-pressure
  "Solar radiation pressure at one astronomical unit, N/m^2 -- the solar
  constant divided by the speed of light.

  This is M&G's value, which corresponds to a solar constant of 1367 W/m^2.
  The currently accepted total solar irradiance is nearer 1361, so a model
  wanting today's number should use 4.5398e-6 instead; the difference is
  about 0.4% of the radiation-pressure acceleration."
  4.560e-6)

;; -------------------------------------------------------------- relativity

(def c-light
  "Speed of light, km/s."
  299792.458)

;; -------------------------------------------------------------------- time

(def mjd-J2000
  "Modified Julian Date of the J2000.0 epoch, 2000 January 1, 12h."
  51544.5)

(def jd-mjd-offset 2400000.5)
