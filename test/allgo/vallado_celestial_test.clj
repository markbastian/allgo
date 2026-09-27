(ns allgo.vallado-celestial-test
  "Vallado chapter 5's lines of sight, shadow, eclipses, twilight and
  naked-eye visibility, each checked against the geometry it rests on: a
  tangent line's closest approach, the cylinder's closed-form eclipse,
  the Sun's altitude at the twilight found, and sites placed by hand in
  daylight and in the dark."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.frames :as frames]
            [allgo.astro.rise :as rise]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [allgo.astro.visibility :as vis]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private R c/R-earth)
(defn- at-angle [r deg] (let [a (math/to-radians deg)] [(* r (math/cos a)) (* r (math/sin a)) 0.0]))

(deftest line-of-sight
  (testing "two points at 7000 km: the chord's closest approach is 7000 cos(half the angle)"
    (is (vis/sight? (at-angle 7000.0 0) (at-angle 7000.0 30)) "6761 km: clear")
    (is (not (vis/sight? (at-angle 7000.0 0) (at-angle 7000.0 60))) "6062 km: blocked")
    (is (not (vis/sight? (at-angle 7000.0 0) (at-angle 7000.0 180))) "opposite sides"))
  (testing "the segment, not the line: both ends on one side see each other"
    (is (vis/sight? [7000.0 0.0 0.0] [9000.0 100.0 0.0])))
  (testing "a point inside the Earth sees nothing"
    (is (not (vis/sight? [6000.0 0.0 0.0] [6000.0 10.0 0.0]))))
  (testing "over the pole, 6370 km up: the sphere blocks it, the ellipsoid does not"
    (let [a [-5000.0 0.0 6370.0] b [5000.0 0.0 6370.0]]
      (is (not (vis/sight? a b)))
      (is (vis/sight? a b {:oblate? true})))))

(def ^:private sun [c/AU 0.0 0.0])

(deftest shadow
  (testing "straight behind the Earth: umbra; to the side: sunlit"
    (is (= :umbra (:state (vis/shadow [-7000.0 0.0 0.0] sun))))
    (is (= :sunlit (:state (vis/shadow [0.0 7000.0 0.0] sun))))
    (is (= :sunlit (:state (vis/shadow [7000.0 0.0 0.0] sun)))))
  (testing "at the cylinder's edge the cone has it in penumbra"
    (let [edge [-7000.0 R 0.0]]
      (is (= :penumbra (:state (vis/shadow edge sun))))
      (is (< 0.0 (:fraction (vis/shadow edge sun)) 1.0))))
  (testing "the cylinder: all or nothing"
    (is (= :umbra (:state (vis/shadow [-7000.0 (- R 1.0) 0.0] sun {:model :cylinder}))))
    (is (= :sunlit (:state (vis/shadow [-7000.0 (+ R 1.0) 0.0] sun {:model :cylinder}))))
    (is (= :sunlit (:state (vis/shadow [7000.0 0.0 0.0] sun {:model :cylinder}))) "the day side")))

(deftest eclipses
  (let [r 7000.0
        n (math/sqrt (/ c/GM-earth (* r r r)))
        period (/ (* 2 math/PI) n)
        position (fn [t] (at-angle r (math/to-degrees (* n t))))
        [e & more] (vis/eclipses position (constantly sun) 0.0 period 10.0 1e-3)
        {:keys [penumbra-in umbra-in umbra-out penumbra-out]} e
        ;; the cylinder's eclipse: behind the Earth, within R of the axis
        cylinder (/ (* 2.0 (math/asin (/ R r))) n)]
    (testing "one eclipse an orbit, entered and left through the penumbra"
      (is (empty? more))
      (is (< penumbra-in umbra-in umbra-out penumbra-out)))
    (testing "centered on the anti-Sun point, half an orbit on"
      (is (< (abs (- (* 0.5 (+ umbra-in umbra-out)) (* 0.5 period))) 1e-2))
      (is (< (abs (- (* 0.5 (+ penumbra-in penumbra-out)) (* 0.5 period))) 1e-2)))
    (testing "the cone's umbra is shorter than the cylinder, its penumbra longer"
      (is (< (- umbra-out umbra-in) cylinder (- penumbra-out penumbra-in)))
      (is (< (- (- penumbra-out penumbra-in) (- umbra-out umbra-in)) 30.0)
          "a low orbit crosses the penumbra in seconds"))
    (testing "each boundary is where the shadow changes"
      (is (= 1.0 (:fraction (vis/shadow (position (- penumbra-in 0.01)) sun))))
      (is (< (:fraction (vis/shadow (position (+ penumbra-in 0.01)) sun)) 1.0))
      (is (pos? (:fraction (vis/shadow (position (- umbra-in 0.01)) sun))))
      (is (zero? (:fraction (vis/shadow (position (+ umbra-in 0.01)) sun))))))
  (testing "an interval that starts in shadow leaves the entry nil"
    (let [r 7000.0 n (math/sqrt (/ c/GM-earth (* r r r)))
          position (fn [t] (at-angle r (+ 180.0 (math/to-degrees (* n t)))))
          [e] (vis/eclipses position (constantly sun) 0.0 2000.0 10.0 1e-3)]
      (is (nil? (:penumbra-in e)))
      (is (some? (:penumbra-out e))))))

(defn- sun-altitude
  "The Sun's altitude from `[lat lon]` at UT MJD `ut`."
  [[lat lon] ut]
  (let [tt (+ ut (/ (time/delta-t (time/decimal-year ut)) 86400.0))
        [ra dec] (solar/apparent-equatorial tt)
        H (- (+ (frames/gast ut tt) lon) ra)]
    (math/asin (+ (* (math/sin lat) (math/sin dec)) (* (math/cos lat) (math/cos dec) (math/cos H))))))

(deftest twilight
  (let [place [(math/to-radians 39.7) (math/to-radians -105.0)]
        mjd 60580.0]
    (testing "at the times found the Sun's center is at the twilight's altitude"
      (doseq [kind [:civil :nautical :astronomical]
              :let [{:keys [rise set]} (rise/twilight kind place mjd)
                    h (get rise/twilight-altitude kind)]]
        (is (< (abs (- (sun-altitude place (+ mjd rise)) h)) (math/to-radians 0.01)) (str kind " dawn"))
        (is (< (abs (- (sun-altitude place (+ mjd set)) h)) (math/to-radians 0.01)) (str kind " dusk"))))
    (testing "and they bracket sunrise and sunset in order"
      (let [sun (rise/sun place mjd)
            civil (rise/twilight :civil place mjd)
            astro (rise/twilight :astronomical place mjd)
            ;; Denver's local midnight falls near 7h UT, so measure from there
            day (fn [m] (mod (- m 0.3) 1.0))]
        (is (< (day (:rise astro)) (day (:rise civil)) (day (:rise sun))))
        (is (< (day (:set sun)) (day (:set civil)) (day (:set astro)))))))
  (testing "at 60 north at midsummer the Sun never gets 18 degrees down"
    (let [place [(math/to-radians 60.0) 0.0]
          june (time/calendar->mjd 2024 6 21)]
      (is (nil? (rise/twilight :astronomical place june)))
      (is (some? (rise/twilight :civil place june)) "though it passes 6"))))

(deftest naked-eye
  ;; the Sun far along +y; sites on the equator at an angle from it
  (let [sun [0.0 c/AU 0.0]
        site (fn [deg] (let [a (math/to-radians deg)] [(* R (math/sin a)) (* R (math/cos a)) 0.0]))
        overhead (fn [deg alt] (v3/scale (v3/normalize (site deg)) (+ R alt)))]
    (testing "overhead, 1000 km up, 100 degrees round from the Sun: site dark, satellite lit"
      (is (vis/visible? (site 100) (overhead 100 1000.0) sun)))
    (testing "150 degrees round, it is in the Earth's shadow too"
      (is (not (vis/visible? (site 150) (overhead 150 1000.0) sun))))
    (testing "in daylight, nothing is visible"
      (is (not (vis/visible? (site 30) (overhead 30 1000.0) sun))))
    (testing "just after sunset the sky is still too bright"
      (is (not (vis/visible? (site 93) (overhead 93 1000.0) sun))))
    (testing "below the horizon, never"
      (is (not (vis/visible? (site 100) (overhead 160 1000.0) sun))))
    (testing "elevation: straight up is 90 degrees"
      (is (< (abs (- (vis/elevation (site 100) (overhead 100 500.0)) (/ math/PI 2))) 1e-12)))))
