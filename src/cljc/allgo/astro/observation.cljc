(ns allgo.astro.observation
  "Tracking measurements and what corrupts them
  (Montenbruck & Gill chapter 6).

  An orbit is never observed directly. What a station measures is a range, a
  Doppler shift or a pair of angles, and each is separated from the geometry
  by a chain of effects: the signal takes time to arrive, during which the
  satellite moves; the atmosphere slows it; the ionosphere slows it by an
  amount depending on frequency. Modeling those is most of the work in
  turning tracking data into an orbit.

  The corrections are not small next to the accuracy wanted. A meter matters
  for precise orbit determination, and the troposphere alone contributes 2.4
  meters straight up and twenty-five at five degrees elevation."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as geodesy]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

;; ------------------------------------------------------------------- range

(defn station-velocity
  "Inertial velocity of a ground station from the Earth's rotation, km/s.
  At the equator it is 465 m/s, which is 6% of orbital speed and cannot be
  left out of a range-rate model."
  [r-station]
  (let [[x y _] r-station]
    [(* (- c/omega-earth) y) (* c/omega-earth x) 0.0]))

(defn range-and-rate
  "Geometric range and range-rate from a station to a satellite, both given
  in the same inertial frame.

  Range-rate is the projection of relative velocity onto the line of sight,
  and it is what a Doppler measurement actually senses -- the component of
  motion along the beam, blind to everything across it."
  [r-station v-station r-sat v-sat]
  (let [d  (v3/sub r-sat r-station)
        dv (v3/sub v-sat v-station)
        rho (v3/length d)]
    {:range rho
     :range-rate (/ (v3/dot d dv) rho)}))

(defn light-time
  "Seconds for light to cross a range, and the correction that follows.

  A signal leaves the satellite before it is received, and the satellite has
  moved in between. For low Earth orbit the delay is a few milliseconds and
  the motion some tens of meters -- negligible for pointing an antenna,
  decisive for ranging at the centimeter level.

  Solved by iteration because the answer feeds its own input: how far the
  satellite was depends on when it was, which depends on how far it was."
  [r-station sat-at t-receive]
  (loop [tau 0.0 n 0]
    (let [r-sat (sat-at (- t-receive tau))
          tau'  (/ (v3/distance r-sat r-station) c/c-light)]
      (if (or (< (abs (- tau' tau)) 1e-12) (>= n 20))
        {:light-time tau' :position r-sat :range (* tau' c/c-light)}
        (recur tau' (inc n))))))

;; ------------------------------------------------------------- troposphere

(def ^:private zenith-hydrostatic 2.3)
(def ^:private zenith-wet 0.1)

(defn- chao-mapping
  "How much longer the path is at elevation `el` than straight up.

  The naive answer is 1/sin(el), which diverges at the horizon; the added
  term keeps it finite, since a real ray bends rather than grazing forever."
  [el a b]
  (/ 1.0 (+ (math/sin el) (/ a (+ (math/tan el) b)))))

(defn tropospheric-delay
  "Extra path length from the neutral atmosphere, km.

  Two components with different behavior. The hydrostatic part is 90% of it
  and is predictable from surface pressure alone. The wet part is small but
  varies with humidity along the path, which is why it is the limiting error
  in precise geodesy: nobody knows where the water vapor is."
  ([el] (tropospheric-delay el zenith-hydrostatic zenith-wet))
  ([el zhd zwd]
   (if (<= el 0.0)
     ##Inf
     (* 1e-3 (+ (* zhd (chao-mapping el 0.00143 0.0445))
                (* zwd (chao-mapping el 0.00035 0.017)))))))

;; -------------------------------------------------------------- ionosphere

(def L1 1575.42e6)
(def L2 1227.60e6)

(def ^:private iono-constant 40.3)

(defn ionospheric-delay
  "Group delay from free electrons along the path, km.

  Proportional to the total electron content and inversely to the square of
  the frequency. That dispersion is the whole reason satellite navigation
  uses two frequencies: an effect that depends on frequency can be measured
  and removed, and one that does not cannot."
  [tec-units frequency]
  (* 1e-3 (/ (* iono-constant tec-units 1e16) (* frequency frequency))))

(defn ionosphere-free
  "Combine two frequencies to cancel the ionosphere.

  The delay scales as 1/f^2, so a particular weighted difference of the two
  measurements has none of it left. The cost is noise: the combination
  amplifies whatever measurement error each carried, by about a factor of
  three."
  [obs-1 f1 obs-2 f2]
  (let [f1sq (* f1 f1) f2sq (* f2 f2)]
    (/ (- (* f1sq obs-1) (* f2sq obs-2)) (- f1sq f2sq))))

;; ------------------------------------------------------------------ angles

(defn look-angles
  "Azimuth, elevation and range from an Earth-fixed station to an
  Earth-fixed satellite position. See `allgo.astro.geodesy`."
  [r-station r-sat]
  (geodesy/look-angles r-station r-sat))

(defn right-ascension-declination
  "Equatorial angles of a target as seen from a station, radians. This is
  what an optical or radar telescope reports when it is not pointing in
  local horizon coordinates."
  [r-station r-sat]
  (let [[x y z :as d] (v3/sub r-sat r-station)
        rho (v3/length d)]
    {:right-ascension (let [a (math/atan2 y x)] (if (neg? a) (+ a c/two-pi) a))
     :declination     (math/asin (/ z rho))
     :range           rho}))

;; ------------------------------------------------------------------ Doppler

(defn doppler-shift
  "Frequency shift from line-of-sight motion, Hz.

  Negative when closing, since an approaching source crowds its waves. For a
  low orbit at L-band this runs to tens of kilohertz and reverses sign as
  the satellite passes overhead -- the signature that identifies a pass."
  [range-rate frequency]
  (* (- frequency) (/ range-rate c/c-light)))

;; ---------------------------------------------------------- the full model

(defn modeled-range
  "Range a station would measure, with the corrections chapter 6 accounts
  for: geometry, then the media, then the clocks.

  Clock offsets enter as a range because that is how a one-way system sees
  them -- a receiver clock a microsecond fast reports every satellite 300
  meters too far, which is why navigation solves for four unknowns and not
  three."
  [{:keys [geometric elevation tec frequency
           station-clock satellite-clock]
    :or   {tec 0.0 frequency L1 station-clock 0.0 satellite-clock 0.0}}]
  (let [tropo (if elevation (tropospheric-delay elevation) 0.0)
        iono  (if (pos? tec) (ionospheric-delay tec frequency) 0.0)]
    {:range (+ geometric tropo iono
               (* c/c-light (- station-clock satellite-clock)))
     :troposphere tropo
     :ionosphere iono
     :clock (* c/c-light (- station-clock satellite-clock))}))

;; ------------------------------------------------------------ carrier phase
;;
;; The carrier is a thousand times more precise than the code and useless on
;; its own, because a receiver counts cycles from an arbitrary starting point
;; and cannot tell which cycle it began on. Everything below is arranged
;; around that: an unknown integer that stays constant while lock is held,
;; and combinations designed to isolate or eliminate it.
;;
;; The other thing to keep straight is a sign. Free electrons slow the code
;; and speed the carrier by the same amount -- the group and phase velocities
;; move in opposite directions in a dispersive medium -- so a pseudorange
;; comes out too long and a phase range too short. That opposition is not a
;; nuisance; it is what makes the ionosphere measurable.

(defn wavelength [frequency] (/ c/c-light frequency))

(defn carrier-phase
  "Phase range in km: geometry, plus the media, plus an unknown whole number
  of wavelengths.

  Returned as a range rather than in cycles, which is what the estimator
  wants and keeps it comparable with a pseudorange."
  [{:keys [geometric elevation tec frequency ambiguity
           station-clock satellite-clock]
    :or   {tec 0.0 frequency L1 ambiguity 0 station-clock 0.0 satellite-clock 0.0}}]
  (let [tropo (if elevation (tropospheric-delay elevation) 0.0)
        iono  (ionospheric-delay tec frequency)]
    {:phase (+ geometric tropo (- iono)
               (* (wavelength frequency) ambiguity)
               (* c/c-light (- station-clock satellite-clock)))
     :troposphere tropo
     :ionosphere  (- iono)
     :ambiguity   (* (wavelength frequency) ambiguity)}))

(defn geometry-free
  "The difference of two phase ranges. Everything that does not depend on
  frequency -- range, clocks, troposphere -- cancels, leaving the ionosphere
  and a constant.

  It is therefore both a measurement of the ionosphere and the standard way
  to spot a cycle slip: the constant only changes when the receiver loses
  and regains lock, so a step in this quantity is a slip and nothing else."
  [phase-1 phase-2]
  (- phase-1 phase-2))

(defn widelane
  "The combination whose effective wavelength is c/(f1 - f2).

  At 86 cm against L1's 19 the ambiguity is far easier to pin to an integer,
  which is why resolution starts here and works down."
  [phase-1 f1 phase-2 f2]
  (/ (- (* f1 phase-1) (* f2 phase-2)) (- f1 f2)))

(defn narrowlane
  "The combination with wavelength c/(f1 + f2), about 11 cm. Noisier in
  ambiguity terms and quieter in every other, which is the opposite trade to
  the widelane."
  [obs-1 f1 obs-2 f2]
  (/ (+ (* f1 obs-1) (* f2 obs-2)) (+ f1 f2)))

(defn melbourne-wubbena
  "Widelane phase less narrowlane code.

  Geometry cancels because both combinations carry it identically, and the
  ionosphere cancels because it enters the two with opposite signs. What
  survives is the widelane ambiguity and noise -- an integer, observable
  directly, with no orbit or clock knowledge whatever."
  [phase-1 code-1 f1 phase-2 code-2 f2]
  (- (widelane phase-1 f1 phase-2 f2)
     (narrowlane code-1 f1 code-2 f2)))

(defn cycle-slip?
  "Whether the geometry-free combination has stepped by more than `tolerance`
  km between epochs -- the signature of the receiver having lost count.

  The default is 5 cm, comfortably under the 19 cm of a single L1 cycle,
  because the whole point is to catch one lost cycle. A threshold in meters
  would miss every slip worth detecting and report a clean series."
  ([previous current] (cycle-slip? previous current 5e-5))
  ([previous current tolerance]
   (> (abs (- current previous)) tolerance)))
