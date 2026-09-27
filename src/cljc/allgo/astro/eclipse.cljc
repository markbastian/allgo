(ns allgo.astro.eclipse
  "Eclipses of the Sun and Moon (Meeus, *Astronomical Algorithms*,
  chapter 54).

  An eclipse needs a new or full moon near a node of the lunar orbit, and
  the test for that is the Moon's argument of latitude F at the syzygy:
  within about 13.9 degrees of a node there may be an eclipse, beyond 21
  there cannot be. When there may be, the instant of greatest eclipse
  follows from series like those for the phases, and two quantities
  decide what kind it is: gamma, the least distance of the shadow's axis
  from the Earth's center in equatorial radii, and u, the radius of the
  umbral cone in the fundamental plane, which is negative when the cone's
  vertex falls short of the Earth -- a total rather than annular eclipse.

  These are the circumstances for the Earth as a whole: where an eclipse
  is visible needs the Besselian elements, which Meeus leaves out.

  Times are MJD (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.numerics.interpolation :refer [horner]]
            [clojure.math :as math]))

(defn- deg [x] (* x c/degrees))

(defn- circumstances
  "Greatest eclipse and gamma, u and M' for syzygy k, or nil if the Moon
  is too far from a node."
  [k c1 c2]
  (let [T (/ k 1236.85)
        F (deg (horner T [160.7108 (* 390.67050284 1236.85) -0.0016118 -0.00000227 0.000000011]))]
    (when (<= (abs (math/sin F)) 0.36)
      (let [E  (horner T [1.0 -0.002516 -0.0000074])
            M  (deg (horner T [2.5534 (* 29.10535670 1236.85) -0.0000014 -0.00000011]))
            M' (deg (horner T [201.5643 (* 385.81693528 1236.85) 0.0107582 0.00001238 -0.000000058]))
            om (deg (horner T [124.7746 (* -1.56375588 1236.85) 0.0020672 0.00000215]))
            F1 (- F (* (deg 0.02665) (math/sin om)))
            A1 (deg (horner T [299.77 (* 0.107408 1236.85) -0.009173]))
            s math/sin cs math/cos
            jd (+ 2451550.09766 (* 29.530588861 k)
                  (* T T (horner T [0.00015437 -0.000000150 0.00000000073]))
                  (* c1 (s M')) (* c2 E (s M))
                  (* 0.0161 (s (* 2 M'))) (* -0.0097 (s (* 2 F1)))
                  (* 0.0073 E (s (- M' M))) (* -0.0050 E (s (+ M' M)))
                  (* -0.0023 (s (- M' (* 2 F1)))) (* 0.0021 E (s (* 2 M)))
                  (* 0.0012 (s (+ M' (* 2 F1)))) (* 0.0006 E (s (+ (* 2 M') M)))
                  (* -0.0004 (s (* 3 M'))) (* -0.0003 E (s (+ M (* 2 F1))))
                  (* 0.0003 (s A1)) (* -0.0002 E (s (- M (* 2 F1))))
                  (* -0.0002 E (s (- (* 2 M') M))) (* -0.0002 (s om)))
            P (+ (* 0.2070 E (s M)) (* 0.0024 E (s (* 2 M))) (* -0.0392 (s M'))
                 (* 0.0116 (s (* 2 M'))) (* -0.0073 E (s (+ M' M))) (* 0.0067 E (s (- M' M)))
                 (* 0.0118 (s (* 2 F1))))
            Q (+ 5.2207 (* -0.0048 E (cs M)) (* 0.0020 E (cs (* 2 M))) (* -0.3299 (cs M'))
                 (* -0.0060 E (cs (+ M' M))) (* 0.0041 E (cs (- M' M))))
            W (abs (cs F1))]
        {:mjd (- jd c/jd-mjd-offset)
         :gamma (* (+ (* P (cs F1)) (* Q (s F1))) (- 1.0 (* 0.0048 W)))
         :u (+ 0.0059 (* 0.0046 E (cs M)) (* -0.0182 (cs M')) (* 0.0004 (cs (* 2 M')))
               (* -0.0005 (cs (+ M M'))))
         :M' M'}))))

(defn- syzygy [year q]
  (+ (math/floor (+ (- (* (- year 2000.0) 12.3685) q) 0.5)) q))

(defn solar
  "The solar eclipse at the new moon nearest `year`, or nil if there is
  none: `{:type :mjd :gamma :u :penumbra :central? :magnitude}`, where
  `:u` and `:penumbra` are the radii of the umbral and penumbral cones in
  the fundamental plane, in Earth radii, and `:type` is
  :partial, :annular, :total or :hybrid (annular along part of the track,
  total along the rest), `:mjd` is the instant of greatest eclipse, and
  `:magnitude` -- the fraction of the Sun's diameter covered -- is given
  for partial eclipses. A non-central total or annular eclipse is one
  where the umbra grazes the Earth near a pole; Meeus reports those as
  total whether the cone's vertex reaches the Earth or not."
  [year]
  (when-let [{:keys [mjd gamma u]} (circumstances (syzygy year 0.0) -0.4075 0.1721)]
    (let [g (abs gamma)]
      (when (<= g (+ 1.5433 u))
        (let [central? (< g 0.9972)
              type (cond
                     (not central?) (if (< g (+ 0.9972 (abs u))) :total :partial)
                     (neg? u) :total
                     (> u 0.0047) :annular
                     (< u (* 0.00464 (math/sqrt (- 1.0 (* gamma gamma))))) :hybrid
                     :else :annular)]
          (cond-> {:type type :mjd mjd :gamma gamma :u u :penumbra (+ u 0.5461) :central? central?}
            (= type :partial) (assoc :magnitude (/ (- (+ 1.5433 u) g) (+ 0.5461 (* 2.0 u))))))))))

(defn lunar
  "The lunar eclipse at the full moon nearest `year`, or nil: `{:type :mjd
  :gamma :magnitude :semiduration}`. `:type` is :penumbral, :partial or
  :total; `:magnitude` the fraction of the Moon's diameter in the umbra --
  or, for a penumbral eclipse, in the penumbra; `:semiduration` maps each
  phase the eclipse reaches (:total :partial :penumbral) to half its
  length, in days. `:rho` and `:sigma` are the radii of the penumbra and
  umbra at the Moon's distance, in equatorial radii of the Earth."
  [year]
  (when-let [{:keys [mjd gamma u M']} (circumstances (syzygy year 0.5) -0.4065 0.1727)]
    (let [g  (abs gamma)
          um (/ (- 1.0128 u g) 0.545)
          pm (/ (- (+ 1.5573 u) g) 0.545)
          type (cond (> um 1.0) :total (pos? um) :partial (pos? pm) :penumbral)]
      (when type
        (let [n  (+ 0.5458 (* 0.04 (math/cos M')))
              half (fn [r] (/ (math/sqrt (- (* r r) (* gamma gamma))) n 24.0))]
          {:type type :mjd mjd :gamma gamma
           :rho (+ 1.2848 u) :sigma (- 0.7403 u)
           :magnitude (if (= type :penumbral) pm um)
           :semiduration (cond-> {:penumbral (half (+ 1.5573 u))}
                           (#{:total :partial} type) (assoc :partial (half (- 1.0128 u)))
                           (= type :total) (assoc :total (half (- 0.4678 u))))})))))

(defn solar-eclipses
  "Every solar eclipse from the new moon nearest `from` to that nearest
  `to` (decimal years), as `solar` describes them."
  [from to]
  (keep #(solar (+ 2000.0 (/ % 12.3685)))
        (range (syzygy from 0.0) (inc (syzygy to 0.0)))))

(defn lunar-eclipses
  "Every lunar eclipse between the decimal years `from` and `to`."
  [from to]
  (keep #(lunar (+ 2000.0 (/ % 12.3685)))
        (range (syzygy from 0.5) (inc (syzygy to 0.5)))))
