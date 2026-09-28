(ns allgo.astro.constellation
  "Constellations and their coverage (Walker, \"Satellite constellations\",
  Journal of the British Interplanetary Society 37, 1984; Vallado, chapter
  11).

  A Walker delta pattern i:t/p/f puts t satellites in p circular orbits
  of equal inclination i, their ascending nodes spaced 360/p degrees
  apart, t/p satellites evenly spaced in each, and each plane's satellites
  a phase 360 f/t degrees ahead of the plane to the west -- f, from 0 to p
  - 1, sets how the planes interleave. Coverage is geometry: a satellite
  at radius r sees the ground above elevation eps out to the Earth-
  central angle lambda = arccos(R cos eps/r) - eps, and a point is covered
  when some satellite lies within that angle of it.

  Angles radians; positions Earth-centered, km. Coverage is computed in
  the inertial frame, the Earth's rotation, where it matters, the
  caller's to supply."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(defn walker
  "The elements `{:a :e :i :raan :argp :M}` of the Walker delta pattern
  i:t/p/f -- inclination `i`, `t` satellites, `p` planes, phasing `f` --
  at semi-major axis `a`, the first satellite at the first plane's node:
  plane k's node at 2 pi k/p, and its satellite j at 2 pi j/(t/p) + 2 pi
  f k/t along the orbit."
  [a i t p f]
  (assert (zero? (mod t p)) "t must be a multiple of p")
  (let [s (quot t p)]
    (vec (for [k (range p) j (range s)]
           {:a a :e 0.0 :i i :raan (/ (* 2.0 math/PI k) p) :argp 0.0
            :M (mod (+ (/ (* 2.0 math/PI j) s) (/ (* 2.0 math/PI f k) t)) (* 2.0 math/PI))}))))

(defn coverage-angle
  "lambda, the Earth-central half-angle of the ground a satellite at
  radius `r` sees above elevation `eps`: arccos(R cos eps/r) - eps."
  ([r eps] (coverage-angle c/R-earth r eps))
  ([R r eps] (- (math/acos (/ (* R (math/cos eps)) r)) eps)))

(defn positions
  "The inertial positions of the satellites with elements `elements`,
  `dt` seconds on (two-body)."
  ([elements dt] (positions c/GM-earth elements dt))
  ([mu elements dt]
   (mapv (fn [{:keys [a M] :as el}]
           (first (kepler/elements->state mu (assoc el :M (+ M (* dt (kepler/mean-motion mu a)))))))
         elements)))

(defn in-view
  "How many of the satellites at `sats` (positions) see the ground point
  of unit direction `u` above elevation `eps`."
  ([sats u eps] (in-view c/R-earth sats u eps))
  ([R sats u eps]
   (count (filter (fn [s] (let [r (v3/length s)
                                lam (coverage-angle R r eps)]
                            (>= (/ (v3/dot s u) r) (math/cos lam))))
                  sats))))

(defn- grid
  "Ground points on an equal-area grid: `n` bands of equal area in sin
  latitude, each cut into points in proportion to its circumference."
  [n]
  (vec (for [k (range n)
             :let [z (- 1.0 (/ (* 2.0 (+ k 0.5)) n))
                   rho (math/sqrt (- 1.0 (* z z)))
                   m (max 1 (long (math/round (* 2.0 n rho))))]
             j (range m)
             :let [lon (* 2.0 math/PI (/ (+ j 0.5) m))]]
         {:u [(* rho (math/cos lon)) (* rho (math/sin lon)) z]
          :w (/ 1.0 (* n m))})))

(defn coverage
  "The coverage of the globe by the satellites at `sats` above elevation
  `eps`: `{:fraction :mean :min :max}` -- the fraction of the area seen by
  at least one, the mean number in view over the area, and the fewest and
  most in view anywhere -- on an equal-area grid of `:bands` bands
  (default 90)."
  ([sats eps] (coverage sats eps {}))
  ([sats eps {:keys [bands] :or {bands 90}}]
   (let [pts (grid bands)
         counts (mapv #(in-view sats (:u %) eps) pts)]
     {:fraction (reduce + (map (fn [{:keys [w]} n] (if (pos? n) w 0.0)) pts counts))
      :mean (reduce + (map (fn [{:keys [w]} n] (* w n)) pts counts))
      :min (apply min counts) :max (apply max counts)})))

(defn gaps
  "The coverage gaps at the ground point of unit direction `u` (inertial,
  so the Earth's rotation is left out) by the constellation `elements`
  above elevation `eps`, over `duration` seconds sampled every `step`:
  `{:covered :longest}` -- the fraction of the time some satellite is in
  view, and the longest wait, seconds."
  [elements u eps duration step]
  (let [seen (mapv #(pos? (in-view (positions elements %) u eps)) (range 0.0 duration step))
        runs (->> seen (partition-by identity) (filter (comp not first)) (map count))]
    {:covered (/ (count (filter true? seen)) (double (count seen)))
     :longest (* step (reduce max 0 runs))}))
