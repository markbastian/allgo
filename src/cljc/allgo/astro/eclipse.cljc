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
  `shadow` and `obscuration` answer it at an instant instead, from the
  true positions of the Sun and Moon: where the shadow's axis meets the
  ground, and how much of the Sun is hidden from any point.

  Times are MJD (TT)."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geodesy :as geo]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.polynomial :as poly]
            [clojure.math :as math]))

(defn- circumstances
  "Greatest eclipse and gamma, u and M' for syzygy k, or nil if the Moon
  is too far from a node."
  [k c1 c2]
  (let [T (/ k 1236.85)
        F (c/deg (poly/horner-ascending [160.7108 (* 390.67050284 1236.85) -0.0016118 -0.00000227 0.000000011] T))]
    (when (<= (abs (math/sin F)) 0.36)
      (let [E  (poly/horner-ascending [1.0 -0.002516 -0.0000074] T)
            M  (c/deg (poly/horner-ascending [2.5534 (* 29.10535670 1236.85) -0.0000014 -0.00000011] T))
            M' (c/deg (poly/horner-ascending [201.5643 (* 385.81693528 1236.85) 0.0107582 0.00001238 -0.000000058] T))
            om (c/deg (poly/horner-ascending [124.7746 (* -1.56375588 1236.85) 0.0020672 0.00000215] T))
            F1 (- F (* (c/deg 0.02665) (math/sin om)))
            A1 (c/deg (poly/horner-ascending [299.77 (* 0.107408 1236.85) -0.009173] T))
            s math/sin cs math/cos
            jd (+ 2451550.09766 (* 29.530588861 k)
                  (* T T (poly/horner-ascending [0.00015437 -0.000000150 0.00000000073] T))
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
  where the umbra (or antumbra) grazes the Earth near a pole; which of
  the two it is follows from the sign of u, as for a central one --
  negative, the umbral cone's vertex beyond the Earth, total."
  [year]
  (when-let [{:keys [mjd gamma u]} (circumstances (syzygy year 0.0) -0.4075 0.1721)]
    (let [g (abs gamma)]
      (when (<= g (+ 1.5433 u))
        (let [central? (< g 0.9972)
              type (cond
                     (not central?) (if (< g (+ 0.9972 (abs u))) (if (neg? u) :total :annular) :partial)
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

;; The Moon's shadow on the Earth, from the true positions rather than the
;; series above: what the demo draws, and where on the ground an eclipse is.

(def ^:private sun-radius 695700.0)
(def ^:private moon-radius 1737.4)

(defn- disk-overlap
  "The fraction of a disk of radius `a` hidden by one of radius `b` whose
  center is `d` from its own: the lens they share over the first's area."
  [a b d]
  (cond
    (>= d (+ a b)) 0.0
    (<= d (abs (- a b))) (min 1.0 (/ (* b b) (* a a)))
    :else (let [ca (math/acos (/ (+ (* d d) (* a a) (- (* b b))) (* 2.0 d a)))
                cb (math/acos (/ (+ (* d d) (* b b) (- (* a a))) (* 2.0 d b)))
                lens (- (+ (* a a ca) (* b b cb))
                        (* 0.5 (math/sqrt (* (+ (- d) a b) (+ d a (- b)) (+ d (- a) b) (+ d a b)))))]
            (/ lens (* math/PI a a)))))

(defn obscuration
  "How much of the Sun the Moon hides from `p`: the fraction of its disk's
  area, 0 to 1, both disks taken as evenly bright. `moon`, `sun` and `p`
  are positions in km in any one frame. Whether the Sun is up at `p` is
  another question: this is the two disks alone."
  [moon sun p]
  (let [m (v3/sub moon p) s (v3/sub sun p)]
    (disk-overlap (math/asin (/ sun-radius (v3/length s)))
                  (math/asin (/ moon-radius (v3/length m)))
                  (v3/angle m s))))

(defn shadow
  "Where the Moon's shadow falls, from the geocentric Moon and Sun in km
  in the Earth-fixed frame (z through the north pole, x through
  Greenwich).

  The shadow's axis runs from the Sun through the Moon. `:gamma` is its
  least distance from the Earth's center in equatorial radii, positive
  when it passes north of it, and `:u` the umbral cone's radius there in
  the same units, negative when the cone has not yet closed -- the
  quantities `solar` gets from its series -- and `:penumbra` the
  penumbral cone's, which touches the Earth somewhere while |gamma| is
  less than 1 + `:penumbra`, the flattening aside. When the axis meets the
  ellipsoid, `:lat` and `:lon` (geodetic, radians) are where it first
  does, and `:radius` the shadow's radius there across the axis, km:
  `:kind` :total for the umbra, :annular for the antumbra beyond the
  cone's vertex. Nil when the shadow points away from the Earth -- the
  Moon on the far side of it from the Sun."
  [moon sun]
  (let [d (v3/normalize (v3/sub moon sun))
        slope (/ (- sun-radius moon-radius) (v3/distance moon sun))
        cone (fn [t] (- moon-radius (* t slope)))
        t0 (- (v3/dot moon d))
        closest (v3/add-scaled moon d t0)
        a geo/a-earth
        ;; the ellipsoid is a sphere once z is stretched by a/b
        squash (fn [[x y z]] [x y (* z (/ a geo/b-earth))])
        m' (squash moon) d' (squash d)
        A (v3/dot d' d') B (* 2.0 (v3/dot m' d')) C (- (v3/dot m' m') (* a a))
        disc (- (* B B) (* 4.0 A C))
        axis {:gamma (* (if (neg? (closest 2)) -1.0 1.0) (/ (v3/length closest) a))
              :u (/ (- (cone t0)) a)
              :penumbra (/ (+ moon-radius (* t0 (/ (+ sun-radius moon-radius) (v3/distance moon sun)))) a)}]
    (cond
      (not (pos? t0)) nil
      (neg? disc) axis
      :else
      (let [t (/ (- (- B) (math/sqrt disc)) (* 2.0 A))
            [lat lon] (geo/cartesian->geodetic (v3/add-scaled moon d t))
            r (cone t)]
        (assoc axis :lat lat :lon lon :radius (abs r) :kind (if (pos? r) :total :annular))))))
