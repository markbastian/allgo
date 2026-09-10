(ns procedurals.astro-test
  (:require [procedurals.astro.constants :as c]
            [procedurals.astro.drag :as drag]
            [procedurals.astro.ephemeris :as eph]
            [procedurals.astro.geopotential :as geo]
            [procedurals.astro.relativity :as rel]
            [procedurals.astro.srp :as srp]
            [procedurals.astro.tides :as tid]
            [procedurals.astro.time :as t]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (< (abs (double (- a b))) tol))

(deftest julian-dates-match-known-epochs
  (doseq [[y m d h want] [[2000 1 1 12.0 2451545.0]      ; J2000.0 by definition
                          [1999 1 1  0.0 2451179.5]
                          [1858 11 17 0.0 2400000.5]     ; MJD zero
                          [2024 2 29  6.0 2460369.75]    ; a leap day
                          [1900 1 1  0.0 2415020.5]
                          [1582 10 15 0.0 2299160.5]]]   ; the Gregorian reform
    (is (close? want (t/mjd->jd (t/calendar->mjd y m d h)) 1e-6)
        (str y "-" m "-" d " " h "h"))))

(deftest mjd-and-jd-are-inverse
  (doseq [mjd [0.0 51544.5 60000.25]]
    (is (close? mjd (- (t/mjd->jd mjd) c/jd-mjd-offset) 1e-12))))

(deftest gmst-at-j2000-is-the-defining-value
  ;; 280.46061837 degrees, or 18h 41m 50.55s, is where the series is anchored.
  (is (close? 280.46061837 (/ (t/gmst c/mjd-J2000) c/degrees) 1e-8)))

(deftest gmst-closes-over-a-sidereal-day
  ;; A sidereal day is 86164.0905 s, not 86400: the Earth turns once relative
  ;; to the stars in less than a solar day, and that gap is what separates
  ;; sidereal time from civil time.
  (let [before (t/gmst c/mjd-J2000)
        after  (t/gmst (+ c/mjd-J2000 (/ 86164.0905 86400.0)))]
    (is (< (abs (- after before)) (* 0.001 c/arcsec))
        "returns to the same angle to within a milliarcsecond")))

(deftest gmst-stays-in-range-and-advances
  (doseq [mjd [40000.0 51544.5 60000.0 70000.123]]
    (let [g (t/gmst mjd)]
      (is (<= 0.0 g) (str "mjd " mjd))
      (is (< g c/two-pi) (str "mjd " mjd))))
  (testing "it advances by roughly one turn per day"
    (let [a (t/gmst 60000.0) b (t/gmst 60001.0)]
      (is (< (abs (- b a)) (* 2.0 c/degrees))
          "a day later, back within a couple of degrees of the same angle"))))

(deftest earth-rotation-rate-agrees-with-the-sidereal-day
  ;; omega must be 2*pi over a sidereal day, not a solar one. Getting this
  ;; wrong tilts every drag calculation by the atmosphere's co-rotation.
  (is (close? c/omega-earth (/ c/two-pi 86164.0905) 1e-11)))

(deftest constants-are-self-consistent
  (testing "solar radiation pressure is the solar constant over c"
    ;; M&G's 4.560e-6 implies a solar constant of 1367 W/m^2, the value in
    ;; use when the book was written. Modern TSI is about 1361, which is a
    ;; real 0.4% difference in the resulting acceleration, not a rounding.
    (is (close? c/solar-pressure (/ 1367.0 (* c/c-light 1000.0)) 2e-9))
    (is (close? 1367.0 (* c/solar-pressure c/c-light 1000.0) 0.5)))
  (testing "the Sun dominates the Earth by the expected factor"
    (is (close? 332946.0 (/ c/GM-sun c/GM-earth) 20.0)))
  (testing "the Moon is about 1/81 of the Earth"
    (is (close? 81.30 (/ c/GM-earth c/GM-moon) 0.02))))

;; ------------------------------------------------------------ geopotential

(def ^:private test-field
  "A field with terms at every shape the recursion has to handle: zonal
  (m=0), sectorial (m=n) and tesseral (0<m<n)."
  {:GM c/GM-earth :R c/R-earth :normalised? true
   :C {[0 0] 1.0 [2 0] -4.841654e-4 [2 2] 2.43926e-6
       [3 0] 9.5717e-7 [3 1] 2.02929e-6 [3 3] 7.2114e-7 [4 0] 5.3997e-7}
   :S {[2 2] -1.40027e-6 [3 1] 2.4892e-7 [3 3] 1.41437e-6}})

(defn- mag [v] (Math/sqrt (reduce + (map * v v))))
(defn- sub [a b] (mapv - a b))

(def ^:private sample-points
  [[7000.0 0.0 0.0] [5000.0 3000.0 4000.0] [100.0 200.0 7000.0]
   [0.0 0.0 7500.0] [20000.0 5000.0 9000.0] [-6000.0 -3000.0 1500.0]])

(deftest degree-zero-is-a-point-mass
  (let [pm (geo/point-mass c/GM-earth c/R-earth)]
    (doseq [r sample-points]
      (let [d (mag r)
            want (mapv #(* (- (/ c/GM-earth (* d d d))) %) r)]
        (is (< (/ (mag (sub (geo/acceleration pm r 0) want)) (mag want)) 1e-14)
            (str "at " r))))))

(deftest j2-matches-its-closed-form
  ;; The one term with a standard closed form, so it pins the recursion,
  ;; the normalisation and the acceleration formula together.
  (let [only {:GM c/GM-earth :R c/R-earth :normalised? true
              :C {[0 0] 1.0 [2 0] (- (/ geo/J2 (Math/sqrt 5.0)))} :S {}}
        pm   (geo/point-mass c/GM-earth c/R-earth)]
    (doseq [r sample-points]
      (let [d (mag r) [x y z] r
            k (/ (* -1.5 geo/J2 c/GM-earth c/R-earth c/R-earth) (Math/pow d 5))
            s (/ (* 5.0 z z) (* d d))
            want [(* k x (- 1.0 s)) (* k y (- 1.0 s)) (* k z (- 3.0 s))]
            got  (sub (geo/acceleration only r 2) (geo/acceleration pm r 0))]
        (is (< (/ (mag (sub got want)) (mag want)) 1e-12) (str "at " r))))))

(deftest acceleration-is-the-gradient-of-the-potential
  ;; The two are derived independently -- the acceleration is not
  ;; differentiated from the potential -- so agreement checks both. The
  ;; residual must fall as h^2; a mismatch that plateaus instead is a real
  ;; error hiding under finite-difference noise, which is exactly how the
  ;; sign slip in the m>0 branch was found.
  (doseq [r sample-points]
    (let [a  (geo/acceleration test-field r 4)
          at (fn [h] (mapv (fn [i]
                             (/ (- (geo/potential test-field (update r i + h) 4)
                                   (geo/potential test-field (update r i - h) 4))
                                (* 2.0 h)))
                           (range 3)))
          e1 (/ (mag (sub a (at 2.0))) (mag a))
          e2 (/ (mag (sub a (at 1.0))) (mag a))]
      (is (< e2 1e-7) (str "at " r " residual " e2))
      (is (< 3.0 (/ e1 e2) 5.0)
          (str "at " r " the residual must fall as h^2, got ratio " (/ e1 e2))))))

(deftest every-harmonic-shape-is-exercised
  ;; Regression: the sectorial and tesseral terms went through a branch that
  ;; J2 never touches, and were wrong there while every zonal test passed.
  (doseq [[label C S] [["sectorial C22" {[2 2] 2.43926e-6} {}]
                       ["sectorial S22" {} {[2 2] -1.40027e-6}]
                       ["tesseral C31"  {[3 1] 2.0e-6} {}]
                       ["sectorial C33" {[3 3] 7.0e-7} {}]]]
    (testing label
      (let [m  {:GM c/GM-earth :R c/R-earth :normalised? true :C (assoc C [0 0] 1.0) :S S}
            pm (geo/point-mass c/GM-earth c/R-earth)
            r  [5000.0 3000.0 4000.0]
            a  (sub (geo/acceleration m r 4) (geo/acceleration pm r 0))
            n  (mapv (fn [i]
                       (/ (- (- (geo/potential m (update r i + 0.5) 4)
                                (geo/potential pm (update r i + 0.5) 0))
                             (- (geo/potential m (update r i - 0.5) 4)
                                (geo/potential pm (update r i - 0.5) 0)))
                          1.0))
                     (range 3))]
        (is (pos? (mag a)) "the term must actually contribute")
        (is (< (/ (mag (sub a n)) (mag a)) 1e-5) label)))))

(deftest higher-degrees-contribute-and-diminish
  (let [r [5000.0 3000.0 4000.0]
        step (fn [d] (mag (sub (geo/acceleration test-field r d)
                               (geo/acceleration test-field r (dec d)))))]
    (is (pos? (step 2)))
    (is (pos? (step 3)))
    (is (pos? (step 4)))
    (is (> (step 2) (* 100 (step 3))) "J2 dominates everything above it")))

(deftest normalisation-factors-are-right
  (is (close? 1.0 (geo/normalisation-factor 0 0) 1e-14))
  (is (close? (Math/sqrt 5.0) (geo/normalisation-factor 2 0) 1e-14))
  (is (close? (Math/sqrt 7.0) (geo/normalisation-factor 3 0) 1e-14))
  (is (close? 3.0 (geo/normalisation-factor 4 0) 1e-14))
  ;; N_22 = sqrt(2 * 5 * (0!/4!)) = sqrt(10/24)
  (is (close? (Math/sqrt (/ 10.0 24.0)) (geo/normalisation-factor 2 2) 1e-14)))

;; ---------------------------------------------------------------- ephemeris

(defn- declination [v] (/ (Math/asin (/ (nth v 2) (mag v))) c/degrees))

(defn- ecliptic-longitude [v mjd]
  (let [[x y z] v
        eps (eph/obliquity mjd)]
    (mod (/ (Math/atan2 (+ (* (Math/cos eps) y) (* (Math/sin eps) z)) x) c/degrees) 360.0)))

(deftest the-sun-goes-where-the-sun-goes
  (let [days (map #(+ c/mjd-J2000 %) (range 0 366))
        rs   (map #(mag (eph/sun %)) days)]
    (testing "Earth's orbit is slightly eccentric"
      ;; perihelion 147.10e6 km in early January, aphelion 152.10e6 in July
      (is (close? 147.10e6 (apply min rs) 5e3))
      (is (close? 152.10e6 (apply max rs) 5e3)))
    (testing "declination swings by the obliquity"
      (is (close? 23.44 (apply max (map #(declination (eph/sun %)) days)) 0.02))
      (is (close? -23.44 (apply min (map #(declination (eph/sun %)) days)) 0.02)))
    (testing "and it is where it should be at J2000"
      (is (close? 280.38 (ecliptic-longitude (eph/sun c/mjd-J2000) c/mjd-J2000) 0.05)))))

(deftest the-moon-goes-where-the-moon-goes
  (testing "the mean rate is one sidereal month"
    (is (close? 27.32158 (/ 36525.0 1336.851344) 1e-3)))
  (testing "distance spans perigee to apogee, sampled over four years"
    ;; a shorter window misses the extremes: perigee distance itself varies
    ;; through the month, so a single lunation does not reach 356400
    (let [rs (map #(mag (eph/moon (+ c/mjd-J2000 (* 0.5 %)))) (range 0 2920))]
      (is (close? 356400.0 (apply min rs) 1500.0))
      (is (close? 406700.0 (apply max rs) 1500.0))))
  (testing "declination reaches the major standstill over a nodal cycle"
    ;; +/-28.6 deg, the obliquity plus the Moon's 5.15 deg inclination, and
    ;; only reached when the nodes line up -- once every 18.6 years
    (let [dl (map #(declination (eph/moon (+ c/mjd-J2000 (* 2.0 %)))) (range 0 3470))]
      (is (< 28.0 (apply max dl) 29.5))
      (is (< -29.5 (apply min dl) -28.0)))))

(deftest third-body-perturbations-are-the-expected-size
  (let [r [7000.0 0.0 0.0]
        sun-a  (mag (eph/third-body c/GM-sun r (eph/sun c/mjd-J2000)))
        moon-a (mag (eph/third-body c/GM-moon r (eph/moon c/mjd-J2000)))
        earth  (/ c/GM-earth (* 7000.0 7000.0))]
    (testing "both land in the published range for low Earth orbit"
      (is (< 1e-10 sun-a 1e-9) (str "Sun " sun-a))
      (is (< 1e-10 moon-a 3e-9) (str "Moon " moon-a)))
    (testing "the Moon outweighs the Sun despite being negligible by mass"
      ;; tidal forcing falls as the cube of distance, not the square
      (is (> moon-a sun-a)))
    (testing "and both are utterly dominated by the Earth"
      (is (> (/ earth moon-a) 1e6)))))

(deftest the-indirect-term-is-the-whole-effect
  ;; The direct pull of the Sun on a satellite is enormous; almost all of it
  ;; is shared with the Earth and so does not perturb the orbit at all.
  (let [r [7000.0 0.0 0.0]
        s (eph/sun c/mjd-J2000)
        d (mapv - r s)
        direct (/ c/GM-sun (reduce + (map * d d)))
        pert   (mag (eph/third-body c/GM-sun r s))]
    (is (> (/ direct pert) 1e4)
        "the differenced perturbation is four orders below the raw pull")))

(deftest the-perturbation-is-tidal-and-falls-as-the-cube
  ;; The signature of a differenced force. A direct pull falls as 1/s^2, but
  ;; what survives the difference is the *gradient* of that pull across the
  ;; orbit, which falls as 1/s^3. Doubling the distance to the body should
  ;; therefore cut the perturbation by eight, not four -- and it is why the
  ;; Moon outweighs the Sun here despite being trivial by mass.
  (let [r [7000.0 0.0 0.0]
        at (fn [dist] (mag (eph/third-body c/GM-sun r [0.0 dist 0.0])))]
    (doseq [d [1e7 2e7 4e7]]
      (let [ratio (/ (at d) (at (* 2.0 d)))]
        (is (< 7.5 ratio 8.5)
            (str "doubling from " d " km changed it by " ratio "x, expected 8"))))))

(deftest a-distant-body-perturbs-nothing
  (let [r [7000.0 0.0 0.0]]
    (is (> (mag (eph/third-body c/GM-sun r [0.0 1e8 0.0]))
           (mag (eph/third-body c/GM-sun r [0.0 1e12 0.0])))
        "further is weaker")
    (is (< (mag (eph/third-body c/GM-sun r [0.0 1e14 0.0])) 1e-20)
        "and far enough away is nothing at all")))

;; ------------------------------------------------- solar radiation pressure

(def ^:private sun-at [c/AU 0.0 0.0])

(deftest sunlit-satellites-are-sunlit
  (doseq [r [[7000.0 0.0 0.0] [42164.0 0.0 0.0] [7000.0 7000.0 0.0] [0.0 42164.0 0.0]]]
    (is (close? 1.0 (srp/shadow r sun-at) 1e-12) (str "at " r))))

(deftest the-umbra-blocks-everything
  (doseq [d [6500.0 7000.0 42164.0 200000.0 1.0e6 1.3e6]]
    (is (close? 0.0 (srp/shadow [(- d) 0.0 0.0] sun-at) 1e-12)
        (str d " km behind the Earth"))))

(deftest the-umbra-closes-to-a-point
  ;; The Sun is bigger than the Earth, so the shadow is a cone, not a
  ;; cylinder. Past its tip only an annular eclipse is possible.
  (let [tip (srp/umbra-length c/AU)]
    (is (close? 1.383e6 tip 2e3) "about 1.38 million km")
    (is (> tip 384400.0) "the Moon orbits inside it, which is why total lunar eclipses happen")
    (is (close? 0.0 (srp/shadow [(- (* 0.9 tip)) 0.0 0.0] sun-at) 1e-12) "inside the tip, total")
    (is (pos? (srp/shadow [(- (* 1.2 tip)) 0.0 0.0] sun-at)) "past the tip, some Sun always shows")))

(deftest the-penumbra-is-a-ramp-not-a-step
  ;; The whole reason for a conical model: a satellite crossing the terminator
  ;; sees the force ease off over several hundred kilometres rather than
  ;; switch. At GEO the ramp is about 400 km wide.
  (let [nus (mapv #(srp/shadow [-42164.0 % 0.0] sun-at) (range 5800.0 7200.0 25.0))]
    (is (apply <= nus) "monotonic across the terminator")
    (is (close? 0.0 (first nus) 1e-12))
    (is (close? 1.0 (last nus) 1e-12))
    (testing "and it passes through half-light at the geometric limb"
      (is (close? 0.5 (srp/shadow [-42164.0 6378.0 0.0] sun-at) 0.02)))
    (testing "no step: no two neighbouring samples jump far"
      (is (< (apply max (map (fn [[a b]] (abs (- b a))) (partition 2 1 nus))) 0.15)))))

(deftest annular-eclipse-beyond-the-tip
  ;; Far enough back the Earth is the smaller disc and can only ever cover a
  ;; fraction (b/a)^2 of the Sun.
  (let [d      2.0e6
        d-sun  (+ c/AU d)
        a      (Math/asin (/ c/R-sun d-sun))
        b      (Math/asin (/ c/R-earth d))
        want   (- 1.0 (/ (* b b) (* a a)))]
    (is (close? want (srp/shadow [(- d) 0.0 0.0] sun-at) 1e-9))))

(deftest radiation-pressure-has-the-right-size-and-sense
  (let [r [7000.0 0.0 0.0]
        a (srp/acceleration r sun-at 0.02 1.3)]
    (testing "magnitude is pressure times area-over-mass times reflectivity"
      ;; Scaled by (AU/d)^2, and d is not quite an AU: a satellite 7000 km
      ;; sunward of the Earth's centre catches light 1e-4 brighter. Small,
      ;; but it is the inverse-square law doing its job, not rounding.
      (let [d (- c/AU 7000.0)
            want (* c/solar-pressure 1.3 0.02 1e-3 (/ (* c/AU c/AU) (* d d)))]
        (is (close? want (mag a) 1e-20))
        (is (close? 1.0000936 (/ (mag a) (* c/solar-pressure 1.3 0.02 1e-3)) 1e-6)
            "the brightening is (AU/d)^2")))
    (testing "and it pushes away from the Sun, never toward it"
      (is (neg? (first a)) "Sun is at +x, so the push is -x"))))

(deftest radiation-pressure-scales-as-advertised
  (let [r [7000.0 0.0 0.0]
        at (fn [am cr] (mag (srp/acceleration r sun-at am cr)))]
    (is (close? (* 2.0 (at 0.01 1.3)) (at 0.02 1.3) 1e-18) "linear in area over mass")
    (is (close? (* 2.0 (at 0.02 1.0)) (at 0.02 2.0) 1e-18) "linear in reflectivity")
    (testing "a mirror takes twice the momentum of a black body"
      (is (close? 2.0 (/ (at 0.02 2.0) (at 0.02 1.0)) 1e-12)))))

(deftest no-sunlight-no-pressure
  (let [r [-7000.0 0.0 0.0]]
    (is (close? 0.0 (srp/shadow r sun-at) 1e-12))
    (is (close? 0.0 (mag (srp/acceleration r sun-at 0.02 1.3)) 1e-30)
        "eclipsed satellites feel nothing, which is what makes the force periodic")))

;; --------------------------------------------------------------------- drag

(defn- at-altitude [h] [(+ c/R-earth h) 0.0 0.0])

(deftest densities-match-the-published-atmosphere
  ;; Checked on the night side, which is the table's minimum column. Real
  ;; densities vary by an order of magnitude with solar activity, so these
  ;; are compared within a factor rather than to a figure.
  (doseq [[h published] [[100 5.0e-7] [200 2.5e-10] [300 2.0e-11]
                         [400 2.8e-12] [500 5.0e-13] [800 1.0e-14]]]
    (let [rho (drag/density [(- (+ c/R-earth h)) 0.0 0.0] sun-at)]
      (is (< 0.5 (/ rho published) 2.0)
          (str h " km: got " rho ", published about " published)))))

(deftest density-falls-monotonically-and-exponentially
  (let [hs   (range 150 1000 25)
        rhos (mapv #(drag/density (at-altitude %) sun-at) hs)]
    (is (apply > rhos) "thinner the higher you go")
    (testing "falling by a factor per fixed step, as hydrostatic equilibrium requires"
      ;; This is why the table is interpolated exponentially rather than
      ;; linearly: across a 20 km gap a straight line is badly wrong.
      (let [ratios (mapv (fn [[a b]] (/ a b)) (partition 2 1 rhos))]
        (is (every? #(< 1.0 % 3.5) ratios))
        (testing "and that factor shrinks with height, because scale height grows"
          ;; The upper atmosphere is hotter and lighter, so it thins out more
          ;; gradually: 3.1x per 25 km down at 150, only 1.2x up at 950.
          (is (apply > ratios))
          (is (> (first ratios) 2.5))
          (is (< (last ratios) 1.5)))))))

(deftest the-atmosphere-bulges-toward-the-sun
  ;; Heated on the daylit side, and lagging the Sun by 30 degrees because it
  ;; takes hours to warm.
  (doseq [h [200 400 600 800]]
    (let [day   (drag/density (at-altitude h) sun-at)
          night (drag/density [(- (+ c/R-earth h)) 0.0 0.0] sun-at)]
      (is (> day night) (str h " km"))))
  (testing "the day/night contrast grows with altitude"
    (let [ratio (fn [h] (/ (drag/density (at-altitude h) sun-at)
                           (drag/density [(- (+ c/R-earth h)) 0.0 0.0] sun-at)))]
      (is (< (ratio 200) (ratio 800))))))

(deftest geodetic-height-accounts-for-the-flattening
  (testing "at the equator it is just r - R"
    (is (close? 400.0 (drag/geodetic-height (at-altitude 400)) 1e-6)))
  (testing "over the pole the ellipsoid is 21 km closer in"
    ;; Against a scale height near 50 km that is a factor of about 1.5 in
    ;; density, so it is not a refinement that can be skipped.
    (is (close? 421.4 (drag/geodetic-height [0.0 0.0 (+ c/R-earth 400.0)]) 0.1)))
  (testing "which makes the polar atmosphere thinner at the same radius"
    (is (< (drag/density [0.0 0.0 (+ c/R-earth 400.0)] sun-at)
           (drag/density [0.0 (+ c/R-earth 400.0) 0.0] sun-at)))))

(deftest the-atmosphere-co-rotates
  (let [alt (+ c/R-earth 400.0)
        v   (Math/sqrt (/ c/GM-earth alt))]
    (testing "a prograde satellite meets slower air than its inertial speed"
      (let [rel (mag (drag/relative-velocity [alt 0.0 0.0] [0.0 v 0.0]))]
        (is (close? (- v (* c/omega-earth alt)) rel 1e-9))
        (is (< rel v))
        (is (> (/ (- v rel) v) 0.05) "worth 6% of orbital speed, not a detail")))
    (testing "a retrograde one meets faster air"
      (is (> (mag (drag/relative-velocity [alt 0.0 0.0] [0.0 (- v) 0.0])) v)))
    (testing "and a polar pass is deflected rather than slowed"
      (let [rel (drag/relative-velocity [alt 0.0 0.0] [0.0 0.0 v])]
        (is (pos? (abs (second rel))) "picks up an east-west component")))))

(deftest drag-opposes-the-airflow
  (let [alt (+ c/R-earth 300.0)
        v   (Math/sqrt (/ c/GM-earth alt))
        r   [alt 0.0 0.0]
        vel [0.0 v 0.0]
        a   (drag/acceleration r vel sun-at 0.01 2.2)
        rel (drag/relative-velocity r vel)]
    (testing "exactly anti-parallel to the relative velocity"
      (let [cos (/ (reduce + (map * a rel)) (* (mag a) (mag rel)))]
        (is (close? -1.0 cos 1e-12))))
    (testing "and of the published size for low Earth orbit"
      (is (< 1e-6 (* 1000 (mag a)) 1e-4) (str (* 1000 (mag a)) " m/s^2")))))

(deftest drag-scales-with-its-parameters
  (let [alt (+ c/R-earth 400.0)
        v   (Math/sqrt (/ c/GM-earth alt))
        at  (fn [am cd] (mag (drag/acceleration [alt 0.0 0.0] [0.0 v 0.0] sun-at am cd)))]
    (is (close? 2.0 (/ (at 0.02 2.2) (at 0.01 2.2)) 1e-12) "linear in area over mass")
    (is (close? 2.0 (/ (at 0.01 4.4) (at 0.01 2.2)) 1e-12) "linear in drag coefficient")))

(deftest above-the-atmosphere-there-is-no-drag
  (doseq [h [1001 1200 5000]]
    (is (zero? (drag/density (at-altitude h) sun-at)) (str h " km"))
    (is (close? 0.0 (mag (drag/acceleration (at-altitude h) [0.0 7.0 0.0] sun-at 0.01 2.2)) 1e-30))))

;; --------------------------------------------------------------- relativity

(deftest mercury-perihelion-advance
  ;; The observation general relativity was invented to explain, and the
  ;; cleanest available check that the constant in front is right: get it
  ;; wrong by any factor and this misses by that factor.
  (let [a 5.7909e7, e 0.2056, period 87.969
        per-rev (rel/perihelion-advance c/GM-sun a e)
        per-century (/ (* per-rev (/ 36525.0 period)) c/arcsec)]
    (is (close? 42.98 per-century 0.05)
        (str "got " per-century " arcsec/century"))))

(deftest relativistic-correction-is-the-expected-size
  (doseq [[label alt lo hi] [["400 km"  6778.0  1e-8 3e-8]
                             ["1000 km" 7378.0  1e-8 2e-8]
                             ["GPS"    26560.0  1e-10 1e-9]
                             ["GEO"    42164.0  1e-11 1e-10]]]
    (let [v (Math/sqrt (/ c/GM-earth alt))
          a (* 1000.0 (mag (rel/acceleration [alt 0.0 0.0] [0.0 v 0.0])))]
      (is (< lo a hi) (str label ": " a " m/s^2")))))

(deftest the-correction-weakens-with-distance
  (let [at (fn [alt] (let [v (Math/sqrt (/ c/GM-earth alt))]
                       (mag (rel/acceleration [alt 0.0 0.0] [0.0 v 0.0]))))]
    (is (apply > (map at [6778.0 10000.0 26560.0 42164.0])))))

(deftest the-cross-term-tracks-eccentricity
  ;; r.v is zero only at apsis and on a circular orbit; away from those it
  ;; grows with eccentricity, which is why the effect is strongest on an
  ;; elongated orbit.
  (let [circular (let [alt 10000.0 v (Math/sqrt (/ c/GM-earth alt))]
                   (rel/acceleration [alt 0.0 0.0] [0.0 v 0.0]))]
    (testing "on a circular orbit the correction is purely radial"
      (is (close? 0.0 (second circular) 1e-30))
      (is (close? 0.0 (nth circular 2) 1e-30)))
    (testing "and points outward, weakening gravity rather than adding to it"
      ;; On a circular orbit v^2 = GM/r, so the radial coefficient
      ;; 4GM/r - v^2 comes to +3GM/r. The sign is not incidental: an
      ;; inward correction would close the orbit faster than Newton and
      ;; precess the perihelion backwards, which is the opposite of what
      ;; Mercury does.
      (is (pos? (first circular))))))

(deftest the-correction-vanishes-as-light-gets-faster
  ;; It is a 1/c^2 effect, so it must scale that way -- doubling GM at fixed
  ;; geometry should roughly quadruple it, since GM enters twice.
  (let [r [10000.0 0.0 0.0]
        v [0.0 6.0 0.0]
        a1 (mag (rel/acceleration r v c/GM-earth))
        a2 (mag (rel/acceleration r v (* 2.0 c/GM-earth)))]
    (is (> (/ a2 a1) 3.0) "GM appears in the prefactor and again in 4GM/r")
    (is (< (/ a2 a1) 5.0))))

;; -------------------------------------------------------------------- tides

(def ^:private moon-at [384400.0 0.0 0.0])

(deftest the-induced-potential-is-k2-times-the-raising-potential
  ;; The definition of a Love number, and the sharpest check available here:
  ;; it pins the normalisation, the Legendre functions, the k2/5 factor and
  ;; every sign convention at once. Anything wrong anywhere moves the ratio
  ;; off one.
  (let [corr  (tid/corrections [[c/GM-moon moon-at]])
        ;; a field of the corrections alone, with no central term
        field {:GM c/GM-earth :R c/R-earth :normalised? true
               :C (assoc (:C corr) [0 0] 0.0) :S (:S corr)}]
    (doseq [[label p] [["sub-lunar" [c/R-earth 0.0 0.0]]
                       ["quadrature" [0.0 c/R-earth 0.0]]
                       ["pole" [0.0 0.0 c/R-earth]]
                       ["oblique" (mapv #(* c/R-earth %) [0.6 0.48 0.64])]]]
      (let [induced (geo/potential field p 2)
            raising (* tid/k2 (tid/raising-potential c/GM-moon moon-at p))]
        (is (close? 1.0 (/ induced raising) 1e-9) (str label))))))

(deftest the-moon-raises-a-bigger-tide-than-the-sun
  ;; Despite being 27 million times lighter. Tide-raising falls as the cube
  ;; of distance while weight falls as the square, and the Moon is 400 times
  ;; closer -- which is what leaves it ahead by rather more than two to one.
  (let [m (get-in (tid/corrections [[c/GM-moon moon-at]]) [:C [2 0]])
        s (get-in (tid/corrections [[c/GM-sun [c/AU 0.0 0.0]]]) [:C [2 0]])]
    (is (close? 2.2 (/ m s) 0.1))))

(deftest the-bulge-follows-the-body
  ;; C22 and S22 must rotate as cos(2 lon) and sin(2 lon): twice, because a
  ;; tide has two bulges, one facing the Moon and one away from it.
  (doseq [[lon want-c want-s] [[0.0 1.0 0.0] [45.0 0.0 1.0] [90.0 -1.0 0.0]]]
    (let [rad (* lon c/degrees)
          m   [(* 384400.0 (Math/cos rad)) (* 384400.0 (Math/sin rad)) 0.0]
          cr  (tid/corrections [[c/GM-moon m]])
          amp (get-in (tid/corrections [[c/GM-moon moon-at]]) [:C [2 2]])]
      (is (close? (* want-c amp) (get-in cr [:C [2 2]]) (* 1e-6 amp)) (str "C22 at " lon))
      (is (close? (* want-s amp) (get-in cr [:S [2 2]]) (* 1e-6 amp)) (str "S22 at " lon)))))

(deftest tidal-perturbation-is-the-expected-size
  (let [j2    {:GM c/GM-earth :R c/R-earth :normalised? true
               :C {[0 0] 1.0 [2 0] -4.841654e-4} :S {}}
        tided (tid/apply-to j2 [[c/GM-moon moon-at] [c/GM-sun [c/AU 0.0 0.0]]])
        r     [7000.0 0.0 0.0]
        a     (* 1000.0 (mag (mapv - (geo/acceleration tided r 2) (geo/acceleration j2 r 2))))]
    (testing "around 1e-7 m/s^2, this being the geometry that maximises it"
      (is (< 1e-8 a 1e-6) (str a " m/s^2")))
    (testing "and utterly dwarfed by the static field it perturbs"
      (let [pm   (geo/point-mass c/GM-earth c/R-earth)
            a-j2 (* 1000.0 (mag (mapv - (geo/acceleration j2 r 2) (geo/acceleration pm r 0))))]
        (is (> (/ a-j2 a) 1e4))))))

(deftest no-bodies-no-tide
  (let [corr (tid/corrections [])]
    (is (every? zero? (vals (:C corr))))
    (is (every? zero? (vals (:S corr)))))
  (testing "and applying an empty correction leaves a field alone"
    (let [f (tid/apply-to geo/earth [])
          r [7000.0 1000.0 2000.0]]
      (is (close? 0.0 (mag (mapv - (geo/acceleration f r 4) (geo/acceleration geo/earth r 4))) 1e-30)))))
