(ns allgo.astro-test
  (:require [allgo.astro.constants :as c]
            [allgo.astro.drag :as drag]
            [allgo.astro.ephemeris :as eph]
            [allgo.astro.estimation :as est]
            [allgo.astro.forces :as forces]
            [allgo.astro.frames :as fr]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kep]
            [allgo.astro.observation :as obs]
            [allgo.astro.planets :as pl]
            [allgo.astro.relativity :as rel]
            [allgo.astro.srp :as srp]
            [allgo.astro.tides :as tid]
            [allgo.astro.time :as t]
            [allgo.astro.variational :as var]
            [allgo.numerics :as num]
            [allgo.numerics.linear :as lin]
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
        tided (tid/perturb j2 [[c/GM-moon moon-at] [c/GM-sun [c/AU 0.0 0.0]]])
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
    (let [f (tid/perturb geo/earth [])
          r [7000.0 1000.0 2000.0]]
      (is (close? 0.0 (mag (mapv - (geo/acceleration f r 4) (geo/acceleration geo/earth r 4))) 1e-30)))))

;; ------------------------------------------------------- the assembled model

(defn- cross [[a b cc] [d e g]]
  [(- (* b g) (* cc e)) (- (* cc d) (* a g)) (- (* a e) (* b d))])

(defn- raan
  "Right ascension of the ascending node, from the angular momentum vector.
  Undefined for an equatorial orbit, where there is no node to speak of."
  [r v]
  (let [[hx hy _] (cross r v)]
    (Math/atan2 hx (- hy))))

(def ^:private j2-only
  {:degree 2 :sun? false :moon? false
   :field {:GM c/GM-earth :R c/R-earth :normalised? true
           :C {[0 0] 1.0 [2 0] (- (/ geo/J2 (Math/sqrt 5.0)))} :S {}}})

(deftest j2-precesses-the-node-at-the-textbook-rate
  ;; The end-to-end check that chapter 3 feeds chapter 4 correctly: integrate
  ;; the force model and recover a closed-form perturbation result.
  ;;
  ;;   dOmega/dt = -3/2 J2 (R/p)^2 n cos i
  (doseq [[label alt incl] [["ISS-like" 400.0 51.6] ["polar-ish" 800.0 98.6]]]
    (testing label
      (let [a-orb   (+ c/R-earth alt)
            i       (* incl c/degrees)
            [r0 v0] (forces/circular-state a-orb i)
            n-mean  (Math/sqrt (/ c/GM-earth (* a-orb a-orb a-orb)))
            want    (* -1.5 geo/J2 (let [q (/ c/R-earth a-orb)] (* q q)) n-mean (Math/cos i))
            span    (* 2.0 86400.0)
            end     (num/step-until
                     (num/integrator (first (filter #(= "DOPRI5(4)" (:name %)) num/first-order))
                                     (forces/first-order j2-only c/mjd-J2000)
                                     0.0 (into r0 v0) 10.0 {:tol-abs 1e-10 :tol-rel 1e-10})
                     span)
            y       (:y end)
            got     (/ (- (raan (subvec y 0 3) (subvec y 3 6)) (raan r0 v0)) span)]
        (is (close? 1.0 (/ got want) 0.01)
            (str label " precessed at " (* got 86400.0 (/ 1.0 c/degrees)) " deg/day, closed form "
                 (* want 86400.0 (/ 1.0 c/degrees))))))))

(deftest a-retrograde-orbit-precesses-forward-and-can-be-sun-synchronous
  ;; cos(i) changes sign past 90 degrees, so the node drifts east instead of
  ;; west. At 98.6 degrees and 800 km the rate matches the Earth's own motion
  ;; about the Sun, which is what keeps such an orbit at a fixed local time.
  (let [a-orb (+ c/R-earth 800.0)
        i     (* 98.6 c/degrees)
        n-mean (Math/sqrt (/ c/GM-earth (* a-orb a-orb a-orb)))
        rate  (* -1.5 geo/J2 (let [q (/ c/R-earth a-orb)] (* q q)) n-mean (Math/cos i))
        deg-per-day (* rate 86400.0 (/ 1.0 c/degrees))]
    (is (pos? deg-per-day) "eastward, unlike a prograde orbit")
    (is (close? (/ 360.0 365.25) deg-per-day 0.02)
        (str "sun-synchronous needs " (/ 360.0 365.25) " deg/day, got " deg-per-day))))

(deftest the-perturbation-budget-has-the-expected-hierarchy
  (let [cfg   {:degree 8 :sun? true :moon? true :relativity? true
               :srp {:area-to-mass 0.02 :cr 1.3} :drag {:area-to-mass 0.01 :cd 2.2}}
        [r v] (forces/circular-state (+ c/R-earth 400.0) (* 51.6 c/degrees))
        b     (forces/breakdown cfg c/mjd-J2000 r v)
        m     (fn [k] (* 1000.0 (mag (get b k))))]
    (is (close? 8.7 (m :two-body) 0.2) "central term near 8.7 m/s^2 in low orbit")
    (is (> (m :two-body) (* 100 (m :harmonics))) "J2 is a small correction, not a rival")
    (is (> (m :harmonics) (* 100 (m :drag))))
    (is (> (m :moon) (m :sun)) "the Moon raises the larger tide")
    (is (> (m :srp) (m :relativity)))
    (testing "each lands in its published decade"
      (is (< 1e-3 (m :harmonics) 1e-1))
      (is (< 1e-7 (m :drag) 1e-4))
      (is (< 1e-7 (m :moon) 1e-5))
      (is (< 1e-8 (m :srp) 1e-6))
      (is (< 1e-9 (m :relativity) 1e-7)))))

(deftest earth-fixed-and-inertial-round-trip
  ;; The force model now rotates through the full chapter 5 chain rather than
  ;; sidereal time alone, so this takes a matrix rather than an angle.
  (doseq [r [[7000.0 0.0 0.0] [3000.0 -5000.0 2000.0]]
          mjd [c/mjd-J2000 (t/calendar->mjd 2025 6 1) (t/calendar->mjd 1995 2 14)]]
    (let [u (forces/earth-fixed mjd)]
      (is (< (mag (mapv - r (forces/ecef->eci (forces/eci->ecef r u) u))) 1e-10)
          (str "at mjd " mjd))
      (is (close? (mag r) (mag (forces/eci->ecef r u)) 1e-10)
          "a frame rotation cannot change an altitude"))))

(deftest confusing-tt-with-ut1-costs-thirty-kilometres
  ;; The larger of the two mistakes available here, and the easier to make:
  ;; the rotation angle is a UT1 quantity while the dynamics run on TT, and
  ;; the two differ by 64 s at J2000. The Earth turns 465 m/s at the equator.
  (let [tt  c/mjd-J2000
        r   [c/R-earth 0.0 0.0]
        turn (fn [mjd] (let [g (t/gmst mjd) ca (Math/cos g) sa (Math/sin g) [x y z] r]
                         [(+ (* ca x) (* sa y)) (+ (* (- sa) x) (* ca y)) z]))
        gap (mag (mapv - (turn tt) (turn (t/tt->utc tt))))]
    (is (< 25.0 gap 35.0) (str "got " gap " km"))))

(deftest precession-costs-less-at-the-epoch-and-more-later
  ;; With the time scale handled correctly, what remains is precession and
  ;; nutation. Nil at J2000 by construction, growing at 50 arcseconds a year.
  (let [gap (fn [mjd-tt]
              (let [ut1 (t/tt->utc mjd-tt)
                    r   [c/R-earth 0.0 0.0]
                    full (forces/eci->ecef r (forces/earth-fixed mjd-tt))
                    sidereal-only (let [g (t/gmst ut1) ca (Math/cos g) sa (Math/sin g) [x y z] r]
                                    [(+ (* ca x) (* sa y)) (+ (* (- sa) x) (* ca y)) z])]
                (mag (mapv - full sidereal-only))))]
    (is (< (gap c/mjd-J2000) 1.0) "under a kilometre at the epoch")
    (is (> (gap (t/calendar->mjd 2030 1 1)) 20.0) "tens of kilometres a generation later")
    (is (> (gap (t/calendar->mjd 2050 1 1)) (gap (t/calendar->mjd 2030 1 1))))))

(deftest nystrom-is-refused-when-the-model-needs-velocity
  ;; Drag and the relativistic correction depend on velocity, so the system
  ;; is not y'' = f(t, y) and a Nystrom method cannot represent it. Dropping
  ;; those terms to fit would give a plausible wrong answer, so it throws.
  (is (not (forces/velocity-dependent? {:degree 4 :sun? true})))
  (is (forces/velocity-dependent? {:drag {:area-to-mass 0.01 :cd 2.2}}))
  (is (forces/velocity-dependent? {:relativity? true}))
  (is (fn? (forces/second-order {:degree 4} c/mjd-J2000))
      "a conservative model is fine")
  (is (thrown? clojure.lang.ExceptionInfo
               (forces/second-order {:drag {:area-to-mass 0.01 :cd 2.2}} c/mjd-J2000)))
  (is (thrown? clojure.lang.ExceptionInfo
               (forces/second-order {:relativity? true} c/mjd-J2000))))

(deftest the-first-order-adapter-returns-velocity-then-acceleration
  (let [[r v] (forces/circular-state (+ c/R-earth 500.0) 0.5)
        f     (forces/first-order {:degree 2 :sun? false :moon? false} c/mjd-J2000)
        out   (f 0.0 (into r v))]
    (is (= 6 (count out)))
    (is (= v (subvec (vec out) 0 3)) "first three are the velocity, unchanged")
    (is (< 8.0 (* 1000.0 (mag (subvec (vec out) 3 6))) 9.5) "last three are gravity")))

;; -------------------------------------------------------------- time scales

(deftest leap-seconds-match-the-record
  (doseq [[y m d want] [[1971 6 1 0]    ; before UTC stepped at all
                        [1972 1 1 10] [1972 7 1 11] [1980 1 1 19]
                        [1999 1 1 32] [2000 1 1 32] [2005 12 31 32]
                        [2006 1 1 33] [2012 7 1 35] [2015 7 1 36]
                        [2017 1 1 37] [2026 9 10 37]]]
    (is (= want (t/tai-utc (t/calendar->mjd y m d)))
        (str y "-" m "-" d))))

(deftest leap-seconds-only-ever-increase
  ;; They are inserted to keep UTC tracking a slowing Earth; a negative leap
  ;; second is permitted in principle and has never been needed.
  (is (apply < (map second t/leap-seconds)))
  (is (apply < (map first t/leap-seconds)) "and the table is in date order"))

(deftest the-scales-differ-by-what-they-should
  (let [utc (t/calendar->mjd 2000 1 1 12.0)
        secs (fn [a b] (* 86400.0 (- a b)))]
    (is (close? 32.0 (secs (t/utc->tai utc) utc) 1e-6) "TAI - UTC is the leap count")
    (is (close? 64.184 (secs (t/utc->tt utc) utc) 1e-6) "TT - UTC adds 32.184")
    (is (close? 32.184 (secs (t/tai->tt (t/utc->tai utc)) (t/utc->tai utc)) 1e-6)
        "TT - TAI is fixed by definition, not measured")
    (is (close? 13.0 (secs (t/utc->gps utc) utc) 1e-6)
        "GPS - UTC was 13 s in 2000; it is 18 today, because GPS ignores leaps")))

(deftest gps-time-drifts-from-utc-but-never-from-tai
  (doseq [[y m d] [[1990 1 1] [2000 1 1] [2010 1 1] [2020 1 1]]]
    (let [utc (t/calendar->mjd y m d)
          gps-tai (* 86400.0 (- (t/utc->gps utc) (t/utc->tai utc)))]
      (is (close? -19.0 gps-tai 1e-6)
          (str "constant offset from atomic time at " y)))))

(deftest tt-and-utc-round-trip
  (doseq [[y m d h] [[1985 3 4 0.0] [2000 1 1 12.0] [2015 8 20 6.5] [2026 1 1 18.25]]]
    (let [utc (t/calendar->mjd y m d h)]
      (is (close? utc (t/tt->utc (t/utc->tt utc)) 1e-12) (str y "-" m "-" d)))))

(deftest tdb-stays-within-two-milliseconds-of-tt
  ;; A relativistic effect: a clock on Earth beats at a varying rate against
  ;; one at the barycentre, because Earth's distance from the Sun and its
  ;; speed both vary annually. It cannot accumulate -- it is periodic.
  (let [diffs (map (fn [d] (* 86400.0 1000.0
                              (- (t/tt->tdb (+ c/mjd-J2000 d)) (+ c/mjd-J2000 d))))
                   (range 0 400 3))]
    (is (< (apply max (map abs diffs)) 1.8) "bounded below two milliseconds")
    (is (> (apply max diffs) 1.5) "and it really does swing that far")
    (is (< (apply min diffs) -1.5))
    (testing "the swing is annual, so a year apart the offsets nearly agree"
      (is (close? (* 86400.0 1000.0 (- (t/tt->tdb c/mjd-J2000) c/mjd-J2000))
                  (* 86400.0 1000.0 (- (t/tt->tdb (+ c/mjd-J2000 365.25))
                                       (+ c/mjd-J2000 365.25)))
                  0.05)))))

(deftest a-single-double-mjd-resolves-about-a-microsecond
  ;; Worth knowing where the floor is: every conversion above adds and
  ;; subtracts seconds from a number of order 6e4 days.
  (let [mjd 59000.0
        back (fn [s] (* 86400.0 (- (+ mjd (/ s 86400.0)) mjd)))]
    (is (close? 1.0 (back 1.0) 1e-6) "a second survives")
    (is (< (abs (- 1e-3 (back 1e-3))) 1e-6) "a millisecond survives")
    (is (> (abs (- 1e-9 (back 1e-9))) 1e-11) "a nanosecond does not")))

(deftest ut1-defaults-to-utc-but-accepts-the-observed-offset
  (let [utc (t/calendar->mjd 2020 6 1)]
    (is (= utc (t/utc->ut1 utc)) "no offset given, no offset applied")
    ;; Only to a microsecond: an MJD near 59000 held in one double resolves
    ;; about 1e-11 of a day. Ample for orbit work at the metre level, and the
    ;; reason precise systems carry a two-part Julian date instead.
    (is (close? 0.3 (* 86400.0 (- (t/utc->ut1 utc 0.3) utc)) 1e-5))
    (testing "and the offset is always under a second, by construction"
      ;; leap seconds exist precisely to keep |UT1 - UTC| below 0.9 s
      (is (< (abs (* 86400.0 (- (t/utc->ut1 utc 0.9) utc))) 1.0)))))

;; ------------------------------------------------- precession and nutation

(defn- det3 [[[a b cc] [d e f] [g h i]]]
  (- (+ (* a e i) (* b f g) (* cc d h)) (+ (* cc e g) (* b d i) (* a f h))))

(defn- orthogonality-error [m]
  (let [p (fr/mul m (fr/transpose m))]
    (apply max (for [i (range 3) j (range 3)]
                 (abs (- (nth (nth p i) j) (if (= i j) 1.0 0.0)))))))

(deftest frame-rotations-are-proper-rotations
  ;; Determinant one and orthogonal: they must preserve lengths and angles
  ;; and not turn the frame inside out. A matrix that drifts off this is
  ;; stretching space.
  (doseq [yr [1900 1950 2000 2050 2200]]
    (let [mjd (t/calendar->mjd yr 1 1)]
      (doseq [[label m] [["precession" (fr/precession mjd)] ["nutation" (fr/nutation mjd)]]]
        (is (close? 1.0 (det3 m) 1e-12) (str label " " yr))
        (is (< (orthogonality-error m) 1e-12) (str label " " yr))))))

(deftest precession-is-the-identity-at-its-own-epoch
  ;; It is defined as the rotation from J2000, so at J2000 there is nothing
  ;; to rotate.
  (let [p (fr/precession c/mjd-J2000)]
    (doseq [i (range 3) j (range 3)]
      (is (close? (if (= i j) 1.0 0.0) (nth (nth p i) j) 1e-14)))))

(deftest general-precession-is-fifty-arcseconds-a-year
  ;; The equinox slides along the ecliptic, completing a circuit in about
  ;; 26,000 years. It is why the pole star changes over history.
  (let [v [1.0 0.0 0.0]
        a (fr/apply-m (fr/precession (t/calendar->mjd 2000 1 1)) v)
        b (fr/apply-m (fr/precession (t/calendar->mjd 2100 1 1)) v)
        arcsec-per-year (/ (Math/acos (max -1.0 (min 1.0 (reduce + (map * a b))))) c/arcsec 100.0)]
    (is (close? 50.29 arcsec-per-year 0.05) (str "got " arcsec-per-year))
    (testing "which comes to a full circuit in about 26 millennia"
      (is (close? 25800.0 (/ (* 360.0 3600.0) arcsec-per-year) 200.0)))))

(deftest mean-obliquity-matches-its-defining-value
  (is (close? 84381.448 (/ (fr/mean-obliquity c/mjd-J2000) c/arcsec) 1e-3)
      "23 deg 26 min 21.448 sec at J2000")
  (testing "and decreases by about 47 arcseconds a century"
    (let [drop (/ (- (fr/mean-obliquity (t/calendar->mjd 2000 1 1))
                     (fr/mean-obliquity (t/calendar->mjd 2100 1 1)))
                  c/arcsec)]
      (is (close? 46.815 drop 0.01)))))

(deftest nutation-swings-by-the-right-amount
  ;; Almost entirely the Moon's node circling in 18.6 years. The largest
  ;; single term is 17.2 arcseconds in longitude and 9.2 in obliquity; the
  ;; totals run a little beyond that once the rest of the series is added.
  (let [epochs (map #(+ c/mjd-J2000 (* 20.0 %)) (range 0 400))
        dpsi   (map #(/ (first (fr/nutation-angles %)) c/arcsec) epochs)
        deps   (map #(/ (second (fr/nutation-angles %)) c/arcsec) epochs)]
    (is (< 17.0 (apply max dpsi) 20.0))
    (is (< -20.0 (apply min dpsi) -17.0))
    (is (< 9.0 (apply max deps) 11.0))
    (is (< -11.0 (apply min deps) -9.0))
    (testing "and it is periodic, not secular -- it must not accumulate"
      (is (close? 0.0 (/ (reduce + dpsi) (count dpsi)) 3.0)
          "mean over twenty years is near zero"))))

(deftest nutation-repeats-on-the-lunar-node-cycle
  ;; 18.6 years is the period of the Moon's node, and the dominant term
  ;; tracks it exactly.
  (let [a (first (fr/nutation-angles c/mjd-J2000))
        b (first (fr/nutation-angles (+ c/mjd-J2000 (* 18.613 365.25))))]
    (is (close? a b (* 2.0 c/arcsec)) "back to nearly the same value one cycle on")))

(deftest the-equation-of-the-equinoxes-stays-within-a-second
  ;; Apparent sidereal time runs on the true equinox, which nutation moves,
  ;; so a clock keeping it wanders against the mean by up to about a second.
  (let [seconds (map #(/ (fr/equation-of-equinoxes (+ c/mjd-J2000 (* 20.0 %))) c/arcsec 15.0)
                     (range 0 400))]
    (is (< (apply max (map abs seconds)) 1.3))
    (is (> (apply max (map abs seconds)) 0.9) "and it really does reach about a second")))

(deftest rotations-compose-and-invert
  (let [mjd (t/calendar->mjd 2024 6 1)
        p   (fr/precession mjd)
        n   (fr/nutation mjd)
        both (fr/mul n p)
        v   [0.3 -0.5 0.81]]
    (is (close? 1.0 (det3 both) 1e-12))
    (testing "the transpose undoes the rotation, orthogonality being the point"
      (is (< (mag (mapv - v (fr/apply-m (fr/transpose both) (fr/apply-m both v)))) 1e-14)))
    (testing "and a rotation preserves length"
      (is (close? (mag v) (mag (fr/apply-m both v)) 1e-14)))))

;; ------------------------------------ the full celestial-terrestrial chain

(deftest the-full-transform-is-a-proper-rotation-and-inverts
  (doseq [yr [1990 2000 2025 2100]]
    (let [utc (t/calendar->mjd yr 3 15 7.25)
          tt  (t/utc->tt utc)
          u   (fr/celestial->terrestrial tt utc)
          v   [0.4 -0.6 0.6928]]
      (is (close? 1.0 (det3 u) 1e-12) (str yr))
      (is (< (orthogonality-error u) 1e-12) (str yr))
      (is (< (mag (mapv - v (fr/apply-m (fr/terrestrial->celestial tt utc) (fr/apply-m u v)))) 1e-13)
          (str yr " round trip"))
      (is (close? (mag v) (mag (fr/apply-m u v)) 1e-14) (str yr " length")))))

(deftest a-ground-station-traces-a-circle-in-inertial-space
  (let [ground [c/R-earth 0.0 0.0]
        over-a-day (map (fn [k]
                          (let [utc (+ (t/calendar->mjd 2024 1 1) (/ k 24.0))]
                            (fr/apply-m (fr/terrestrial->celestial (t/utc->tt utc) utc) ground)))
                        (range 25))]
    (testing "at constant radius, since a rotation cannot change length"
      (is (< (- (apply max (map mag over-a-day)) (apply min (map mag over-a-day))) 1e-9)))
    (testing "and after 24 hours it has overshot a full turn by a degree"
      ;; The Earth turns 360.9856 degrees per solar day: one turn relative to
      ;; the stars, plus the degree it has moved around the Sun. That excess
      ;; is exactly 360/365.25.
      (let [a (first over-a-day) b (last over-a-day)
            deg (/ (Math/acos (max -1.0 (min 1.0 (/ (reduce + (map * a b)) (* (mag a) (mag b))))))
                   c/degrees)]
        (is (close? (/ 360.0 365.25) deg 0.01) (str "got " deg " degrees"))))))

(deftest precession-makes-a-gmst-only-rotation-obsolete
  ;; Worth being explicit about: rotating by sidereal time alone ignores that
  ;; the pole and equinox have moved since J2000. The error is nil at the
  ;; epoch and grows at 50 arcseconds a year, which is tens of kilometres at
  ;; the Earth's surface within a couple of decades.
  (let [err (fn [yr]
              (let [utc (t/calendar->mjd yr 1 1)
                    tt  (t/utc->tt utc)
                    full (fr/apply-m (fr/celestial->terrestrial tt utc) [c/R-earth 0.0 0.0])
                    gmst-only (let [g (t/gmst utc)]
                                [(* c/R-earth (Math/cos g)) (* c/R-earth (- (Math/sin g))) 0.0])]
                (mag (mapv - full gmst-only))))]
    (is (< (err 2000) 1.0) "agree to under a kilometre at the epoch itself")
    (is (> (err 2025) 20.0) "but tens of kilometres apart a generation later")
    (is (> (err 2050) (err 2025)) "and it only grows")))

(deftest polar-motion-moves-the-ground-by-metres
  ;; The rotation pole wanders within the crust by about 0.3 arcseconds in a
  ;; 435-day Chandler wobble. Negligible for an orbit, decisive for geodesy.
  (let [utc (t/calendar->mjd 2024 1 1)
        tt  (t/utc->tt utc)
        p   (fn [xp yp] (fr/apply-m (fr/celestial->terrestrial tt utc xp yp) [c/R-earth 0.0 0.0]))
        shift (* 1000.0 (mag (mapv - (p 0.0 0.0) (p (* 0.3 c/arcsec) (* 0.3 c/arcsec)))))]
    (is (< 1.0 shift 20.0) (str "moved " shift " m"))
    (testing "and zero polar motion is exactly no rotation"
      (is (< (orthogonality-error (fr/polar-motion 0.0 0.0)) 1e-15))
      (is (close? 1.0 (nth (nth (fr/polar-motion 0.0 0.0) 0) 0) 1e-15)))))

(deftest apparent-sidereal-time-differs-from-mean-by-the-equinox-equation
  (doseq [yr [2000 2015 2030]]
    (let [utc (t/calendar->mjd yr 5 5 3.0)
          tt  (t/utc->tt utc)
          diff (- (fr/gast utc tt) (t/gmst utc))]
      (is (close? (fr/equation-of-equinoxes tt) diff 1e-14) (str yr))
      (is (< (abs (/ diff c/arcsec 15.0)) 1.3) "under about a second of time"))))

;; ------------------------------------------------------ geodetic coordinates

(deftest geodetic-and-cartesian-round-trip
  (doseq [[lat lon h] [[0.0 0.0 0.0] [45.0 30.0 0.0] [89.9 0.0 0.0]
                       [-33.9 151.2 0.058] [51.5 -0.1 400.0]
                       [0.0 0.0 35786.0] [-60.0 -170.0 1200.0]]]
    (let [r (gd/geodetic->cartesian (* lat c/degrees) (* lon c/degrees) h)
          [la lo hh] (gd/cartesian->geodetic r)]
      (is (close? lat (/ la c/degrees) 1e-9) (str "latitude at " lat))
      (is (close? lon (/ lo c/degrees) 1e-9) (str "longitude at " lon))
      (is (close? h hh 1e-6) (str "height at " h)))))

(deftest geodetic-latitude-differs-from-geocentric
  ;; The distinction the ellipsoid forces. Geodetic latitude is measured
  ;; from the local vertical, geocentric from the centre, and a plumb line
  ;; on an ellipsoid does not point at the middle of the Earth.
  (doseq [[lat expect-arcmin] [[0.0 0.0] [15.0 5.76] [30.0 9.98]
                               [45.0 11.55] [60.0 10.02] [90.0 0.0]]]
    (let [r  (gd/geodetic->cartesian (* lat c/degrees) 0.0 0.0)
          gc (/ (gd/geocentric-latitude r) c/degrees)]
      (is (close? expect-arcmin (* 60.0 (- lat gc)) 0.01) (str "at " lat))))
  (testing "the gap is largest near 45 degrees and vanishes at the equator and poles"
    (let [gap (fn [lat] (let [r (gd/geodetic->cartesian (* lat c/degrees) 0.0 0.0)]
                          (abs (- lat (/ (gd/geocentric-latitude r) c/degrees)))))]
      (is (> (gap 45.0) (gap 20.0)))
      (is (> (gap 45.0) (gap 70.0)))
      (is (close? 0.0 (gap 0.0) 1e-12))
      (is (close? 0.0 (gap 90.0) 1e-9)))))

(deftest the-ellipsoid-is-flattened-by-the-right-amount
  (let [equator (gd/geodetic->cartesian 0.0 0.0 0.0)
        pole    (gd/geodetic->cartesian (/ Math/PI 2.0) 0.0 0.0)]
    (is (close? c/R-earth (mag equator) 1e-9) "equatorial radius is a")
    (is (close? 6356.75 (mag pole) 0.01) "polar radius is b = a(1-f)")
    (is (close? 21.38 (- (mag equator) (mag pole)) 0.01)
        "and the difference is the 21 km that drag has to account for")))

(deftest the-local-horizon-frame-is-orthonormal
  (doseq [lat [0.0 40.0 -70.0] lon [0.0 120.0 -45.0]]
    (let [m (gd/east-north-up (* lat c/degrees) (* lon c/degrees))]
      (doseq [row m] (is (close? 1.0 (mag row) 1e-14)))
      (is (close? 0.0 (reduce + (map * (nth m 0) (nth m 1))) 1e-14))
      (is (close? 0.0 (reduce + (map * (nth m 0) (nth m 2))) 1e-14))
      (is (close? 0.0 (reduce + (map * (nth m 1) (nth m 2))) 1e-14)))))

(deftest look-angles-point-where-they-should
  (let [station (gd/geodetic->cartesian 0.0 0.0 0.0)
        alt     (+ c/R-earth 800.0)]
    (testing "straight overhead"
      (let [{:keys [elevation range]} (gd/look-angles station [alt 0.0 0.0])]
        (is (close? 90.0 (/ elevation c/degrees) 1e-9))
        (is (close? 800.0 range 1e-6))))
    (testing "north is azimuth zero, east is ninety"
      (let [north [(* alt (Math/cos 0.3)) 0.0 (* alt (Math/sin 0.3))]
            east  [(* alt (Math/cos 0.3)) (* alt (Math/sin 0.3)) 0.0]]
        (is (close? 0.0 (/ (:azimuth (gd/look-angles station north)) c/degrees) 1e-9))
        (is (close? 90.0 (/ (:azimuth (gd/look-angles station east)) c/degrees) 1e-9))
        (is (close? (:elevation (gd/look-angles station north))
                    (:elevation (gd/look-angles station east)) 1e-12)
            "symmetric directions give the same elevation")))
    (testing "the far side of the Earth is below the horizon"
      (let [{:keys [elevation]} (gd/look-angles station [(- alt) 0.0 0.0])]
        (is (close? -90.0 (/ elevation c/degrees) 1e-9))
        (is (not (gd/visible? station [(- alt) 0.0 0.0])))))
    (testing "and an elevation mask hides low passes"
      ;; 0.25 rad of geocentric angle, well inside the horizon: for a
      ;; satellite at 800 km that limit is only 27.3 degrees away, since
      ;; cos(theta) = R/(R+h). A pass is a narrow window.
      (let [low [(* alt (Math/cos 0.25)) 0.0 (* alt (Math/sin 0.25))]
            {:keys [elevation]} (gd/look-angles station low)]
        (is (< 10.0 (/ elevation c/degrees) 25.0) "a genuinely low pass")
        (is (gd/visible? station low 0.0))
        (is (not (gd/visible? station low (* 30.0 c/degrees))) "but under a 30 degree mask")))
    (testing "the horizon is closer than intuition suggests"
      ;; cos(theta) = R/(R+h) puts the geometric limit at 27.3 degrees of
      ;; geocentric angle for 800 km -- which is why tracking networks need
      ;; many stations.
      (let [limit (Math/acos (/ c/R-earth alt))]
        (is (close? 27.3 (/ limit c/degrees) 0.2))
        (is (close? 0.0 (/ (:elevation (gd/look-angles station
                                                       [(* alt (Math/cos limit)) 0.0 (* alt (Math/sin limit))]))
                           c/degrees)
                    0.01)
            "and at exactly that angle the target sits on the horizon")))))

(deftest the-rtn-frame-is-orthonormal-and-oriented
  (doseq [[r v] [[[7000.0 0.0 0.0] [0.0 7.5 0.0]]
                 [[3000.0 -5000.0 2000.0] [4.0 3.0 -1.0]]]]
    (let [[R T N] (gd/rtn-frame r v)]
      (doseq [row [R T N]] (is (close? 1.0 (mag row) 1e-14)))
      (is (close? 0.0 (reduce + (map * R T)) 1e-14))
      (is (close? 0.0 (reduce + (map * R N)) 1e-14))
      (is (close? 0.0 (reduce + (map * T N)) 1e-14))
      (testing "radial points outward and normal along the angular momentum"
        (is (pos? (reduce + (map * R r))))
        (is (pos? (reduce + (map * N (let [[a b cc] r [d e f] v]
                                       [(- (* b f) (* cc e)) (- (* cc d) (* a f)) (- (* a e) (* b d))])))))))))

(deftest orbit-errors-resolve-along-track
  ;; Why RTN exists. A small error in speed barely moves the radius but
  ;; accumulates along the track, so after a while a prediction is wrong by
  ;; far more in one direction than the others.
  (let [r  [7000.0 0.0 0.0]
        v  [0.0 7.546 0.0]
        ahead [0.0 12.0 0.0]                    ; 12 km further along the orbit
        [dr dt dn] (gd/to-rtn r v ahead)]
    (is (close? 0.0 dr 1e-12) "no radial component")
    (is (close? 12.0 dt 1e-12) "all of it along track")
    (is (close? 0.0 dn 1e-12) "and none out of plane")))

;; ------------------------------------------------------------ Kepler orbits

(def ^:private mu-e c/GM-earth)

(def ^:private orbit-cases
  ;; Deliberately including every case where an element is undefined, since
  ;; those are the ones a naive conversion returns NaN for.
  [["typical LEO"     {:a 7000.0  :e 0.01 :i 0.9 :raan 1.2 :argp 2.3 :nu 0.7}]
   ["eccentric"       {:a 26000.0 :e 0.7  :i 0.5 :raan 4.0 :argp 1.0 :nu 3.0}]
   ["circular"        {:a 7000.0  :e 0.0  :i 0.9 :raan 1.2 :argp 0.0 :nu 0.7}]
   ["equatorial"      {:a 42164.0 :e 0.01 :i 0.0 :raan 0.0 :argp 2.3 :nu 0.7}]
   ["circ+equatorial" {:a 42164.0 :e 0.0  :i 0.0 :raan 0.0 :argp 0.0 :nu 2.0}]
   ["retrograde"      {:a 7500.0  :e 0.05 :i 2.6 :raan 0.4 :argp 5.0 :nu 1.1}]
   ["polar"           {:a 8000.0  :e 0.1  :i 1.5707963 :raan 0.9 :argp 0.2 :nu 4.5}]])

(deftest elements-and-state-round-trip
  (doseq [[label el] orbit-cases]
    (testing label
      (let [[r v]   (kep/elements->state mu-e el)
            back    (kep/state->elements mu-e r v)
            [r2 v2] (kep/elements->state mu-e back)]
        (is (< (mag (mapv - r r2)) 1e-9) "position")
        (is (< (mag (mapv - v v2)) 1e-12) "velocity")
        (is (close? (:a el) (:a back) 1e-6) "semi-major axis")
        (is (close? (:e el) (:e back) 1e-12) "eccentricity")
        (is (close? (:i el) (:i back) 1e-12) "inclination")))))

(deftest degenerate-orbits-get-a-convention-not-a-nan
  ;; A circular orbit has no periapsis and an equatorial one has no node.
  ;; Both are real orbits; the elements fold the missing angle into the next
  ;; one along rather than failing.
  (doseq [[label el] orbit-cases]
    (let [[r v] (kep/elements->state mu-e el)
          back  (kep/state->elements mu-e r v)]
      (doseq [[k x] back]
        (is (or (= k :M) (and (== x x) (not= x ##Inf) (not= x ##-Inf)))
            (str label " element " k " came out " x))))))

(deftest keplers-equation-solves-to-machine-precision
  ;; M = E - e sin E, checked by residual rather than against a table.
  (doseq [e [0.0 0.1 0.5 0.9 0.99]]
    (doseq [step (range 0 32)]
      (let [M (* c/two-pi (/ step 32.0))
            E (kep/kepler-equation M e)]
        (is (< (abs (- (- E (* e (Math/sin E))) M)) 1e-12)
            (str "e=" e " M=" M))))))

(deftest anomalies-convert-both-ways
  (doseq [e [0.0 0.2 0.6 0.95]
          step (range 0 16)]
    (let [nu (* c/two-pi (/ step 16.0))
          E  (kep/true->eccentric nu e)
          M  (kep/eccentric->mean E e)]
      (is (close? nu (kep/eccentric->true E e) 1e-10) "true -> eccentric -> true")
      (is (close? nu (kep/mean->true M e) 1e-9) "true -> mean -> true")
      (testing "and they agree at the apsides, where all three coincide"
        (when (zero? step)
          (is (close? 0.0 E 1e-12))
          (is (close? 0.0 M 1e-12)))))))

(deftest keplers-third-law
  (testing "period depends only on the semi-major axis, not the shape"
    ;; Orbits of the same size but very different shape take exactly as
    ;; long. Verified by actually propagating each one and finding it back
    ;; where it started, not by comparing the formula to itself.
    (let [T (kep/period mu-e 10000.0)]
      (doseq [e [0.0 0.3 0.7 0.9]]
        (let [[r0 v0] (kep/elements->state mu-e {:a 10000.0 :e e :i 0.4 :raan 0.2 :argp 0.3 :nu 0.0})
              [r1 _]  (kep/propagate mu-e r0 v0 T)
              ;; Relative, because the round trip loses precision as the
              ;; orbit elongates: the half-angle conversions between true
              ;; and eccentric anomaly are delicate near periapsis, where a
              ;; high-eccentricity orbit sweeps fastest. A metre on a
              ;; 10,000 km orbit is 1e-7 relative, and that is the floor.
              err (/ (mag (mapv - r1 r0)) 10000.0)]
          (is (< err 1e-7)
              (str "e=" e " returned within " (* 1e5 err) " cm of its start"))))))
  (testing "and T^2 scales as a^3"
    (let [t1 (kep/period mu-e 7000.0)
          t2 (kep/period mu-e 14000.0)]
      (is (close? (Math/pow 2.0 1.5) (/ t2 t1) 1e-12))))
  (testing "geostationary radius gives exactly one sidereal day"
    ;; which is the definition of the orbit
    (is (close? 86164.09 (kep/period mu-e 42164.17) 0.5))))

(deftest vis-viva-agrees-with-the-state
  (doseq [[label el] orbit-cases]
    (let [[r v] (kep/elements->state mu-e el)]
      (is (close? (mag v) (kep/vis-viva mu-e (mag r) (:a el)) 1e-9) label))))

(deftest analytic-propagation-closes-the-orbit
  (let [el {:a 12000.0 :e 0.3 :i 0.6 :raan 1.0 :argp 2.0 :nu 0.5}
        [r0 v0] (kep/elements->state mu-e el)
        T (kep/period mu-e 12000.0)]
    (doseq [n [1 3 10]]
      (let [[r v] (kep/propagate mu-e r0 v0 (* n T))]
        (is (< (mag (mapv - r r0)) 1e-6) (str n " periods, position"))
        (is (< (mag (mapv - v v0)) 1e-9) (str n " periods, velocity"))))))

(deftest elements-are-constant-under-two-body-motion
  ;; The property that makes them elements at all.
  (let [el {:a 12000.0 :e 0.3 :i 0.6 :raan 1.0 :argp 2.0 :nu 0.5}
        [r0 v0] (kep/elements->state mu-e el)
        T (kep/period mu-e 12000.0)]
    (doseq [frac [0.0 0.25 0.5 0.75 1.0 3.7]]
      (let [[r v] (kep/propagate mu-e r0 v0 (* frac T))
            e' (kep/state->elements mu-e r v)]
        (doseq [k [:a :e :i :raan :argp]]
          (is (close? (get el k) (get e' k) 1e-8) (str k " at " frac " periods")))))))

(deftest analytic-and-numerical-propagation-agree
  ;; Two wholly independent routes to the same answer: Kepler's equation on
  ;; one side, an eighth of a million integration steps on the other.
  (doseq [[label el] [["circular" {:a 7000.0 :e 0.0 :i 0.9 :raan 1.2 :argp 0.0 :nu 0.0}]
                      ["moderate" {:a 12000.0 :e 0.3 :i 0.6 :raan 1.0 :argp 2.0 :nu 0.5}]
                      ["eccentric" {:a 26000.0 :e 0.7 :i 0.5 :raan 4.0 :argp 1.0 :nu 0.0}]]]
    (testing label
      (let [[r0 v0] (kep/elements->state mu-e el)
            T   (kep/period mu-e (:a el))
            f   (fn [_ y] (let [r (subvec (vec y) 0 3) d (mag r)]
                            (into (subvec (vec y) 3 6)
                                  (mapv #(* (- (/ mu-e (* d d d))) %) r))))
            num (:y (num/step-until
                     (num/integrator (first (filter #(= "DOPRI5(4)" (:name %)) num/first-order))
                                     f 0.0 (into r0 v0) 10.0 {:tol-abs 1e-13 :tol-rel 1e-13})
                     (* 3.0 T)))
            [ra _] (kep/propagate mu-e r0 v0 (* 3.0 T))]
        (is (< (mag (mapv - ra (subvec (vec num) 0 3))) 1e-4)
            "agreeing to well under a metre after three orbits")))))

(deftest j2-moves-the-elements-the-way-theory-says
  ;; Chapter 3 in the language of chapter 2. The whole value of elements is
  ;; that a perturbation which looks like noise in a state vector reads as a
  ;; steady drift in two angles and a wobble in the rest.
  (let [el {:a 7000.0 :e 0.01 :i 0.9 :raan 1.0 :argp 2.0 :nu 0.5}
        [r0 v0] (kep/elements->state mu-e el)
        cfg {:degree 2 :sun? false :moon? false :field geo/earth}
        after (fn [secs]
                (let [y (:y (num/step-until
                             (num/integrator (first (filter #(= "DOPRI5(4)" (:name %)) num/first-order))
                                             (forces/first-order cfg c/mjd-J2000) 0.0 (into r0 v0) 10.0
                                             {:tol-abs 1e-11 :tol-rel 1e-11}) secs))]
                  (kep/state->elements mu-e (subvec (vec y) 0 3) (subvec (vec y) 3 6))))
        d1 (after 86400.0)
        d2 (after 172800.0)]
    (testing "size, shape and tilt oscillate but do not march"
      (is (< (abs (- (:a d1) (:a el))) 20.0) "semi-major axis stays put")
      (is (< (abs (- (:e d1) (:e el))) 0.01) "so does eccentricity")
      (is (< (abs (/ (- (:i d1) (:i el)) c/degrees)) 0.1) "and inclination"))
    (testing "while the node regresses steadily, doubling in twice the time"
      (let [r1 (/ (- (:raan d1) (:raan el)) c/degrees)
            r2 (/ (- (:raan d2) (:raan el)) c/degrees)]
        (is (neg? r1) "westward for a prograde orbit")
        (is (close? 2.0 (/ r2 r1) 0.02) "linear in time, as a secular rate must be")
        (testing "at the rate the closed form predicts"
          (let [want (* -1.5 geo/J2 (let [q (/ c/R-earth 7000.0)] (* q q))
                        (Math/sqrt (/ mu-e (* 7000.0 7000.0 7000.0)))
                        (Math/cos 0.9) 86400.0 (/ 180.0 Math/PI))]
            (is (close? 1.0 (/ r1 want) 0.02))))))))

(deftest equinoctial-elements-have-no-singularity
  ;; The reason they exist: a filter cannot estimate an angle that does not
  ;; exist, and a circular orbit's argument of periapsis does not.
  (doseq [[label el] orbit-cases]
    (testing label
      (let [back (kep/from-equinoctial (kep/equinoctial el))]
        (is (close? (:a el) (:a back) 1e-9))
        (is (close? (:e el) (:e back) 1e-12))
        (is (close? (:i el) (:i back) 1e-12))
        (doseq [[k x] (kep/equinoctial el)]
          (is (and (== x x) (not= x ##Inf)) (str "equinoctial " k " is finite"))))))
  (testing "and on a circular orbit the eccentricity vector simply vanishes"
    (let [eq (kep/equinoctial {:a 7000.0 :e 0.0 :i 0.9 :raan 1.2 :argp 0.0 :nu 0.7})]
      (is (close? 0.0 (:h eq) 1e-15))
      (is (close? 0.0 (:k eq) 1e-15)))))

;; ------------------------------------------------------------ observations

(def ^:private equator-station (gd/geodetic->cartesian 0.0 0.0 0.0))

(deftest tropospheric-delay-matches-published-magnitudes
  (testing "2.4 metres straight up"
    (is (close? 2.4 (* 1000.0 (obs/tropospheric-delay (/ Math/PI 2.0))) 0.01)))
  (testing "and about twenty-five at five degrees"
    (is (< 20.0 (* 1000.0 (obs/tropospheric-delay (* 5.0 c/degrees))) 30.0)))
  (testing "growing monotonically as the ray flattens"
    (let [ds (map #(obs/tropospheric-delay (* % c/degrees)) [90 60 30 15 10 5 3])]
      (is (apply < ds))))
  (testing "and staying finite at the horizon, where 1/sin would not"
    (is (< (obs/tropospheric-delay (* 0.5 c/degrees)) 1.0)
        "a real ray bends rather than grazing forever")))

(deftest ionospheric-delay-scales-as-one-over-frequency-squared
  (testing "about 0.16 m per TEC unit at L1"
    (is (close? 0.162 (* 1000.0 (obs/ionospheric-delay 1 obs/L1)) 0.002)))
  (testing "linear in electron content"
    (is (close? 10.0 (/ (obs/ionospheric-delay 100 obs/L1) (obs/ionospheric-delay 10 obs/L1)) 1e-12)))
  (testing "and the two GPS frequencies differ by exactly (f1/f2)^2"
    (let [want (let [r (/ obs/L1 obs/L2)] (* r r))]
      (doseq [tec [1 10 100]]
        (is (close? want (/ (obs/ionospheric-delay tec obs/L2)
                            (obs/ionospheric-delay tec obs/L1))
                    1e-12))))))

(deftest the-ionosphere-free-combination-cancels-it
  ;; The reason navigation uses two frequencies. An effect that depends on
  ;; frequency can be measured and removed; one that does not cannot.
  (doseq [tec [10 50 100 200]]
    (let [truth 20000.0
          m1 (+ truth (obs/ionospheric-delay tec obs/L1))
          m2 (+ truth (obs/ionospheric-delay tec obs/L2))]
      (is (> (* 1000.0 (- m1 truth)) 1.0) "L1 alone is metres out")
      (is (< (abs (* 1000.0 (- (obs/ionosphere-free m1 obs/L1 m2 obs/L2) truth))) 1e-6)
          (str tec " TECU: the combination removes it")))))

(deftest range-rate-is-the-derivative-of-range
  ;; Checked against a numerical derivative with the station held still,
  ;; since the analytic form includes the station's own motion and a
  ;; derivative taken about a fixed point does not. That difference is the
  ;; 465 m/s the equator travels, not an error.
  (let [[r0 v0] (kep/elements->state c/GM-earth
                                     {:a (+ c/R-earth 800.0) :e 0.0 :i 0.9 :raan 0.0 :argp 0.0 :nu 0.3})]
    (doseq [t [0.0 100.0 300.0]]
      (let [[r v] (kep/propagate c/GM-earth r0 v0 t)
            analytic (:range-rate (obs/range-and-rate equator-station [0.0 0.0 0.0] r v))
            h 0.01
            numerical (/ (- (mag (mapv - (first (kep/propagate c/GM-earth r0 v0 (+ t h))) equator-station))
                            (mag (mapv - (first (kep/propagate c/GM-earth r0 v0 (- t h))) equator-station)))
                         (* 2.0 h))]
        (is (close? numerical analytic 1e-6) (str "at t=" t))))))

(deftest a-rotating-station-changes-the-range-rate
  (let [[r v] (kep/elements->state c/GM-earth
                                   {:a (+ c/R-earth 800.0) :e 0.0 :i 0.9 :raan 0.0 :argp 0.0 :nu 0.3})
        fixed   (:range-rate (obs/range-and-rate equator-station [0.0 0.0 0.0] r v))
        rotating (:range-rate (obs/range-and-rate equator-station (obs/station-velocity equator-station) r v))]
    (is (not (close? fixed rotating 0.01)) "worth hundreds of metres per second")
    (is (close? 0.4651 (mag (obs/station-velocity equator-station)) 1e-3)
        "the equator travels 465 m/s")))

(deftest light-time-shows-up-along-the-beam-not-across-it
  ;; A satellite overhead moves perpendicular to the line of sight, so the
  ;; range barely changes during the light travel time and the correction
  ;; nearly vanishes. Low on the horizon it moves along the beam and the
  ;; correction is tens of metres.
  (let [orbit (fn [nu] (kep/elements->state c/GM-earth
                                            {:a (+ c/R-earth 800.0) :e 0.0 :i 0.9
                                             :raan 0.0 :argp 0.0 :nu nu}))
        correction (fn [nu]
                     (let [[r1 v1] (orbit nu)
                           sat-at (fn [t] (first (kep/propagate c/GM-earth r1 v1 t)))
                           lt (obs/light-time equator-station sat-at 0.0)]
                       (* 1000.0 (abs (- (mag (mapv - (sat-at 0.0) equator-station)) (:range lt))))))]
    (is (< (correction 0.0) 1.0) "overhead: nothing to see")
    (is (> (correction 0.45) 10.0) "low: tens of metres")
    (is (> (correction 0.72) (correction 0.45)) "and more the lower it goes")))

(deftest light-time-is-a-few-milliseconds
  (doseq [[label alt want-ms] [["LEO" 800.0 2.67] ["GPS" 20200.0 67.4] ["GEO" 35786.0 119.4]]]
    (is (close? want-ms (* 1000.0 (/ alt c/c-light)) 0.1) label)))

(deftest doppler-reverses-sign-through-a-pass
  ;; The signature that identifies one. Approaching crowds the waves, so the
  ;; shift is positive; receding stretches them.
  (is (pos? (obs/doppler-shift -7.0 obs/L1)) "closing")
  (is (neg? (obs/doppler-shift 7.0 obs/L1)) "receding")
  (is (close? 0.0 (obs/doppler-shift 0.0 obs/L1) 1e-12) "and zero at closest approach")
  (testing "reaching tens of kilohertz at L-band for a low orbit"
    (is (< 30000.0 (abs (obs/doppler-shift -7.0 obs/L1)) 45000.0)))
  (testing "proportional to both speed and frequency"
    (is (close? 2.0 (/ (obs/doppler-shift -6.0 obs/L1) (obs/doppler-shift -3.0 obs/L1)) 1e-12))
    (is (close? (/ obs/L1 obs/L2) (/ (obs/doppler-shift -3.0 obs/L1) (obs/doppler-shift -3.0 obs/L2)) 1e-12))))

(deftest a-clock-error-looks-exactly-like-a-range-error
  ;; Which is why a navigation receiver solves for four unknowns and not
  ;; three: a microsecond of clock offset is 300 m on every satellite at once.
  (let [m (obs/modelled-range {:geometric 20000.0 :station-clock 1e-6})]
    (is (close? 0.2998 (:clock m) 1e-3) "a microsecond is about 300 metres")
    (is (close? (+ 20000.0 (:clock m)) (:range m) 1e-9)))
  (testing "and a satellite clock offset works the other way"
    (let [m (obs/modelled-range {:geometric 20000.0 :satellite-clock 1e-6})]
      (is (neg? (:clock m))))))

(deftest the-modelled-range-adds-up
  (let [m (obs/modelled-range {:geometric 20000.0
                               :elevation (* 10.0 c/degrees)
                               :tec 50 :frequency obs/L1})]
    (is (close? (+ 20000.0 (:troposphere m) (:ionosphere m) (:clock m)) (:range m) 1e-12))
    (is (close? 0.01334 (:troposphere m) 1e-4) "13 m of troposphere at ten degrees")
    (is (close? 0.00812 (:ionosphere m) 1e-4) "8 m of ionosphere at 50 TECU")))

(deftest right-ascension-and-declination-agree-with-the-geometry
  (let [station [0.0 0.0 0.0]]
    (is (close? 0.0 (:declination (obs/right-ascension-declination station [1000.0 0.0 0.0])) 1e-12))
    (is (close? 90.0 (/ (:declination (obs/right-ascension-declination station [0.0 0.0 1000.0])) c/degrees) 1e-12))
    (is (close? 90.0 (/ (:right-ascension (obs/right-ascension-declination station [0.0 1000.0 0.0])) c/degrees) 1e-12))
    (is (close? 1000.0 (:range (obs/right-ascension-declination station [0.0 1000.0 0.0])) 1e-12))))

;; --------------------------------------------------- variational equations

(defn- two-body-accel [_ r _]
  (let [d (mag r)] (mapv #(* (- (/ c/GM-earth (* d d d))) %) r)))

(defn- det-n [m]
  (let [n (count m)]
    (if (= n 1)
      (ffirst m)
      (reduce + (map-indexed
                 (fn [j x] (* (if (even? j) 1.0 -1.0) x
                              (det-n (mapv (fn [row] (vec (concat (subvec row 0 j) (subvec row (inc j)))))
                                           (rest m)))))
                 (first m))))))

(defn- propagate-variational [y0 t]
  (:y (num/step-until
       (num/integrator (first (filter #(= "DOPRI5(4)" (:name %)) num/first-order))
                       (var/rhs two-body-accel) 0.0 y0 5.0 {:tol-abs 1e-12 :tol-rel 1e-12})
       t)))

(deftest the-numerical-gradient-matches-the-closed-form
  ;; Everything but two-body has to be differentiated numerically, and a
  ;; wrong Jacobian degrades convergence rather than breaking anything, so
  ;; it would go unnoticed. Two-body is the one case with a closed form, and
  ;; therefore the only check on that machinery.
  (doseq [r [[7000.0 0.0 0.0] [5000.0 3000.0 4000.0] [42164.0 100.0 -50.0]]]
    (let [exact (var/two-body-gradient c/GM-earth r)
          got   (:d-dr (var/acceleration-gradients two-body-accel 0.0 r [0.0 7.5 0.0]))
          scale (apply max (for [i (range 3) j (range 3)] (abs (nth (nth exact i) j))))]
      (doseq [i (range 3) j (range 3)]
        (is (< (/ (abs (- (nth (nth exact i) j) (nth (nth got i) j))) scale) 1e-7)
            (str "at " r " element " i "," j))))))

(deftest two-body-gradient-is-symmetric-and-traceless
  ;; Both follow from the acceleration being the gradient of a potential.
  ;; Trace zero is Laplace's equation: gravity has no source in empty space.
  (doseq [r [[7000.0 0.0 0.0] [5000.0 3000.0 4000.0]]]
    (let [g (var/two-body-gradient c/GM-earth r)]
      (doseq [i (range 3) j (range 3)]
        (is (close? (nth (nth g i) j) (nth (nth g j) i) 1e-18) "symmetric"))
      (is (close? 0.0 (reduce + (map-indexed (fn [i row] (nth row i)) g)) 1e-15)
          "traceless, since div(a) = 0 away from the mass"))))

(deftest phi-starts-as-the-identity
  (let [{:keys [phi]} (var/unpack (var/initial [7000.0 0.0 0.0] [0.0 7.5 0.0]))]
    (doseq [i (range 6) j (range 6)]
      (is (close? (if (= i j) 1.0 0.0) (nth (nth phi i) j) 1e-15)))))

(deftest phi-is-what-it-claims-to-be
  ;; The definitive test: compare the integrated transition matrix against
  ;; finite differences of the actual propagated trajectory. Nudge each
  ;; initial component, re-propagate, see how the final state moved.
  (let [[r0 v0] (kep/elements->state c/GM-earth
                                     {:a 8000.0 :e 0.1 :i 0.6 :raan 1.0 :argp 2.0 :nu 0.5})
        T (kep/period c/GM-earth 8000.0)]
    (doseq [frac [0.1 0.5 1.0]]
      (let [t   (* frac T)
            phi (:phi (var/unpack (propagate-variational (var/initial r0 v0) t)))
            fd  (mapv (fn [j]
                        (let [h  (if (< j 3) 0.01 1e-5)
                              at (fn [sign]
                                   (propagate-variational
                                    (var/initial (if (< j 3) (update (vec r0) j + (* sign h)) r0)
                                                 (if (>= j 3) (update (vec v0) (- j 3) + (* sign h)) v0))
                                    t))]
                          (mapv (fn [a b] (/ (- a b) (* 2.0 h)))
                                (subvec (vec (at 1.0)) 0 6) (subvec (vec (at -1.0)) 0 6))))
                      (range 6))
            scale (apply max (for [i (range 6) j (range 6)] (abs (nth (nth phi i) j))))]
        (doseq [i (range 6) j (range 6)]
          (is (< (/ (abs (- (nth (nth phi i) j) (nth (nth fd j) i))) scale) 1e-6)
              (str "after " frac " orbits, element " i "," j)))))))

(deftest phase-space-volume-is-conserved
  ;; Liouville's theorem. A conservative system cannot compress phase space,
  ;; so the determinant of the transition matrix is one for all time. It is
  ;; a deep structural property and fails loudly if the variational
  ;; equations are wrong anywhere.
  (let [[r0 v0] (kep/elements->state c/GM-earth
                                     {:a 8000.0 :e 0.1 :i 0.6 :raan 1.0 :argp 2.0 :nu 0.5})
        T (kep/period c/GM-earth 8000.0)]
    (doseq [frac [0.25 1.0 3.0]]
      (is (close? 1.0 (det-n (:phi (var/unpack (propagate-variational (var/initial r0 v0) (* frac T))))) 1e-8)
          (str "after " frac " orbits")))))

(deftest transition-matrices-compose
  ;; Phi(t2, t0) = Phi(t2, t1) Phi(t1, t0). The property that lets an
  ;; estimator accumulate sensitivity across an arc rather than re-deriving
  ;; it from the epoch every time.
  (let [[r0 v0] (kep/elements->state c/GM-earth
                                     {:a 9000.0 :e 0.05 :i 0.7 :raan 0.5 :argp 1.0 :nu 0.2})
        t1 900.0 t2 2400.0
        whole (:phi (var/unpack (propagate-variational (var/initial r0 v0) t2)))
        step1 (var/unpack (propagate-variational (var/initial r0 v0) t1))
        step2 (:phi (var/unpack (propagate-variational (var/initial (:r step1) (:v step1)) (- t2 t1))))
        composed (var/mat-mul step2 (:phi step1))
        scale (apply max (for [i (range 6) j (range 6)] (abs (nth (nth whole i) j))))]
    (doseq [i (range 6) j (range 6)]
      (is (< (/ (abs (- (nth (nth whole i) j) (nth (nth composed i) j))) scale) 1e-7)
          (str "element " i "," j)))))

(deftest position-rate-is-velocity-exactly
  ;; The top half of the Jacobian is not physics and must be exact.
  (let [A (var/jacobian two-body-accel 0.0 [7000.0 100.0 -50.0] [0.1 7.5 0.3])]
    (doseq [i (range 3) j (range 6)]
      (is (close? (if (= j (+ i 3)) 1.0 0.0) (nth (nth A i) j) 1e-15)
          (str "row " i " column " j)))))

;; ------------------------------------------------------- orbit determination

(deftest cholesky-solves-and-refuses
  (testing "an exact solve on a known system"
    (let [A [[4.0 2.0 0.6] [2.0 5.0 1.0] [0.6 1.0 3.0]]
          x [1.0 -2.0 3.0]]
      (doseq [[a b] (map vector x (lin/cholesky-solve A (lin/mat-vec A x)))]
        (is (close? a b 1e-12)))))
  (testing "and an inverse that really is one"
    (let [A [[4.0 2.0 0.6] [2.0 5.0 1.0] [0.6 1.0 3.0]]
          I (lin/mat-mul A (lin/inverse A))]
      (doseq [i (range 3) j (range 3)]
        (is (close? (if (= i j) 1.0 0.0) (nth (nth I i) j) 1e-12)))))
  (testing "a system the data does not determine is refused, not fudged"
    ;; For a normal matrix this is not a numerical mishap but a statement
    ;; about the observations: some direction of the state is unobservable.
    (is (nil? (lin/cholesky-solve [[1.0 1.0] [1.0 1.0]] [1.0 1.0])) "singular")
    (is (nil? (lin/cholesky-solve [[-1.0 0.0] [0.0 1.0]] [1.0 1.0])) "not positive definite")
    (is (nil? (est/solve-batch [{:H [1.0 0.0] :residual 1.0}] 2))
        "one observation cannot fix two unknowns")))

(deftest weighted-least-squares-recovers-a-linear-fit
  ;; Before trusting it on an orbit, check it on something with an answer
  ;; that can be written down.
  (let [truth [2.0 -3.0]
        rows  (mapv (fn [t] {:H [1.0 t]
                             :residual (- (+ (* 2.0 1.0) (* -3.0 t)) 0.0)
                             :weight 1.0})
                    [0.0 1.0 2.0 3.0 4.0])
        {:keys [correction]} (est/solve-batch rows 2)]
    (doseq [[a b] (map vector truth correction)]
      (is (close? a b 1e-10)))))

(deftest weighting-does-what-weighting-should
  ;; Two contradictory observations; the answer must land nearer the one
  ;; trusted more.
  (let [heavy (:correction (est/solve-batch [{:H [1.0] :residual 10.0 :weight 100.0}
                                             {:H [1.0] :residual 0.0 :weight 1.0}] 1))
        even  (:correction (est/solve-batch [{:H [1.0] :residual 10.0 :weight 1.0}
                                             {:H [1.0] :residual 0.0 :weight 1.0}] 1))]
    (is (close? 5.0 (first even) 1e-12) "equal weights split the difference")
    (is (> (first heavy) 9.0) "a hundredfold weight nearly wins outright")))

(deftest the-kalman-update-shrinks-the-covariance
  (let [P0 [[100.0 0.0] [0.0 100.0]]
        {:keys [x P]} (est/kalman-update [0.0 0.0] P0 [1.0 0.0] 10.0 1.0)]
    (is (> (first x) 9.0) "a precise measurement moves the state most of the way")
    (is (< (nth (nth P 0) 0) (nth (nth P0 0) 0)) "and shrinks the variance it informs")
    (is (close? 100.0 (nth (nth P 1) 1) 1e-9) "leaving the unobserved one alone")))

(deftest the-covariance-stays-symmetric-and-positive
  ;; The Joseph form exists for this. A covariance that loses positive
  ;; definiteness gives a negative variance, and the filter is finished.
  (let [rng (java.util.Random. 11)]
    (loop [x [0.0 0.0] P [[10.0 1.0] [1.0 10.0]] n 0]
      (when (< n 200)
        (let [{x' :x P' :P} (est/kalman-update x P [1.0 0.3] (.nextGaussian rng) 0.01)]
          (is (close? (nth (nth P' 0) 1) (nth (nth P' 1) 0) 1e-12) "symmetric")
          (is (pos? (nth (nth P' 0) 0)) "positive variance")
          (is (pos? (nth (nth P' 1) 1)))
          (is (some? (lin/cholesky P')) "and still positive definite")
          (recur x' P' (inc n)))))))

(deftest process-noise-keeps-a-filter-listening
  ;; Without it the covariance shrinks forever and the filter stops learning,
  ;; which is the classic way to make one diverge.
  (let [phi (lin/eye 2)
        P   [[1.0 0.0] [0.0 1.0]]
        no-q  (:P (est/kalman-predict [0.0 0.0] P phi (lin/mat-scale (lin/eye 2) 0.0)))
        with-q (:P (est/kalman-predict [0.0 0.0] P phi (lin/mat-scale (lin/eye 2) 0.5)))]
    (is (close? 1.0 (nth (nth no-q 0) 0) 1e-12) "no noise, no growth")
    (is (close? 1.5 (nth (nth with-q 0) 0) 1e-12) "noise adds uncertainty back")))

(deftest an-orbit-is-recovered-from-noisy-ranges
  ;; The capstone, and the point of the whole package: chapter 2's elements
  ;; set the truth, chapter 3's J2 perturbs it, chapter 4 propagates it,
  ;; chapter 5 places the stations, chapter 6 models the measurement,
  ;; chapter 7 supplies the partials and chapter 8 solves.
  (let [mu c/GM-earth
        j2 geo/J2
        accel (fn [_ r _]
                (let [[x y z] r d (mag r) d2 (* d d)
                      k0 (- (/ mu (* d2 d)))
                      sq (/ (* 5.0 z z) d2)
                      kj (/ (* -1.5 j2 mu c/R-earth c/R-earth) (Math/pow d 5))]
                  [(+ (* k0 x) (* kj x (- 1.0 sq)))
                   (+ (* k0 y) (* kj y (- 1.0 sq)))
                   (+ (* k0 z) (* kj z (- 3.0 sq)))]))
        dopri (first (filter #(= "DOPRI5(4)" (:name %)) num/first-order))
        arc   (fn [x0 times]
                (loop [integ (num/integrator dopri (var/rhs accel) 0.0
                                             (var/initial (subvec (vec x0) 0 3) (subvec (vec x0) 3 6))
                                             10.0 {:tol-abs 1e-11 :tol-rel 1e-11})
                       ts times out []]
                  (if (empty? ts)
                    out
                    (let [s (num/step-until integ (first ts))]
                      (recur s (rest ts) (conj out (var/unpack (:y s))))))))
        stations (mapv (fn [[la lo]] (gd/geodetic->cartesian (* la c/degrees) (* lo c/degrees) 0.0))
                       [[35.0 -117.0] [-25.0 28.0] [40.0 140.0]])
        st-eci (fn [s secs]
                 (let [mjd (+ c/mjd-J2000 (/ secs 86400.0))]
                   (fr/apply-m (fr/terrestrial->celestial (t/utc->tt mjd) mjd) s)))
        truth  (let [[r v] (kep/elements->state mu {:a 7500.0 :e 0.02 :i 0.95
                                                    :raan 1.1 :argp 2.0 :nu 0.4})]
                 (vec (concat r v)))
        sigma  0.010
        ;; A longer arc, and it matters: range alone is weakly observable,
        ;; so a couple of hours from three stations leaves the geometry so
        ;; poorly conditioned that Gauss-Newton walks away from the answer
        ;; rather than toward it. Three hours and 37 observations converge.
        times  (vec (range 60 12060 60))
        t-arc  (arc truth times)
        rng    (java.util.Random. 20260910)
        obs    (vec (for [[i secs] (map-indexed vector times)
                          [idx s] (map-indexed vector stations)
                          :let [r-sat (:r (nth t-arc i))
                                se (st-eci s secs)]
                          :when (> (:elevation (gd/look-angles se r-sat)) (* 10.0 c/degrees))]
                      {:i i :t secs :station idx
                       :measured (+ (mag (mapv - r-sat se)) (* sigma (.nextGaussian rng)))}))
        guess  (mapv + truth [2.0 -1.5 1.0 0.002 0.001 -0.0015])
        step   (fn [x]
                 (let [a (arc x times)]
                   (mapv (fn [{:keys [i t station measured]}]
                           (let [{:keys [r phi]} (nth a i)
                                 d   (mapv - r (st-eci (nth stations station) t))
                                 rho (mag d)
                                 los (mapv #(/ % rho) d)]
                             {:H (mapv (fn [j] (reduce + (map-indexed
                                                          (fn [ii u] (* u (nth (nth phi ii) j))) los)))
                                       (range 6))
                              :residual (- measured rho)
                              :weight (/ 1.0 (* sigma sigma))}))
                         obs)))]
    (is (> (count obs) 30) "enough passes, well enough spread, to determine six unknowns")
    (let [rows0 (step guess)
          final (loop [x (vec guess) n 0]
                  (if (>= n 3)
                    x
                    (recur (mapv + x (:correction (est/solve-batch (step x) 6))) (inc n))))
          rows1 (step final)
          pos-err (mag (mapv - (subvec final 0 3) (subvec truth 0 3)))
          cov (:covariance (est/solve-batch rows1 6))
          formal (Math/sqrt (+ (nth (nth cov 0) 0) (nth (nth cov 1) 1) (nth (nth cov 2) 2)))]
      (is (> (est/rms rows0) 1.0) "starts kilometres out")
      (is (< (est/rms rows1) (* 2.0 sigma)) "and ends at the noise floor")
      (is (< pos-err 0.1)
          (str "recovered to " (* 1000.0 pos-err) " m from a 2.7 km initial error"))
      (testing "and the formal uncertainty is honest, not merely small"
        ;; A covariance that claims more precision than the estimate actually
        ;; has is worse than no covariance at all.
        (is (< 0.2 (/ pos-err formal) 5.0)
            (str "actual error " (* 1000 pos-err) " m against a formal "
                 (* 1000 formal) " m"))))))

(deftest a-good-fit-is-not-a-good-orbit
  ;; The failure mode that matters in practice, and the reason a covariance
  ;; is not optional. With too short an arc only two stations ever see the
  ;; satellite, giving eleven ranges for six unknowns. Gauss-Newton drives
  ;; the residuals to the noise floor -- the fit is excellent -- while the
  ;; orbit walks kilometres away from the truth. Nothing in the residuals
  ;; says so; the formal uncertainty does.
  (let [mu c/GM-earth
        accel (fn [_ r _] (let [d (mag r)] (mapv #(* (- (/ mu (* d d d))) %) r)))
        rk4   (first (filter #(= "RK4" (:name %)) num/first-order))
        arc   (fn [x0 times]
                (loop [integ (num/integrator rk4 (var/rhs accel) 0.0
                                             (var/initial (subvec (vec x0) 0 3) (subvec (vec x0) 3 6))
                                             20.0 {:adaptive? false})
                       ts times out []]
                  (if (empty? ts) out
                      (let [s (num/step-until integ (first ts))]
                        (recur s (rest ts) (conj out (var/unpack (:y s))))))))
        sites (mapv (fn [[la lo]] (gd/geodetic->cartesian (* la c/degrees) (* lo c/degrees) 0.0))
                    [[35.0 -117.0] [-25.0 28.0] [40.0 140.0] [-33.0 151.0] [51.0 0.0]])
        st-eci (fn [s secs] (let [mjd (+ c/mjd-J2000 (/ secs 86400.0))]
                              (fr/apply-m (fr/terrestrial->celestial (t/utc->tt mjd) mjd) s)))
        truth (let [[r v] (kep/elements->state mu {:a 7500.0 :e 0.02 :i 0.95
                                                   :raan 1.1 :argp 2.0 :nu 0.4})]
                (vec (concat r v)))
        sigma 0.010
        solve (fn [hours sample seed]
                (let [times (vec (range 60.0 (* 3600.0 hours) sample))
                      ta    (arc truth times)
                      rng   (java.util.Random. seed)
                      obs   (vec (for [[i secs] (map-indexed vector times)
                                       [idx s] (map-indexed vector sites)
                                       :let [rs (:r (nth ta i)) se (st-eci s secs)]
                                       :when (> (:elevation (gd/look-angles se rs)) (* 10.0 c/degrees))]
                                   {:i i :t secs :station idx
                                    :measured (+ (mag (mapv - rs se)) (* sigma (.nextGaussian rng)))}))
                      e 2.7
                      guess (mapv + truth [(* e 0.74) (* e -0.56) (* e 0.37)
                                           (* e 7.4e-4) (* e 3.7e-4) (* e -5.6e-4)])
                      rows  (fn [x] (let [a (arc x times)]
                                      (mapv (fn [{:keys [i t station measured]}]
                                              (let [{:keys [r phi]} (nth a i)
                                                    d (mapv - r (st-eci (nth sites station) t))
                                                    rho (mag d) los (mapv #(/ % rho) d)]
                                                {:H (mapv (fn [j] (reduce + (map-indexed
                                                                             (fn [ii u] (* u (nth (nth phi ii) j))) los)))
                                                          (range 6))
                                                 :residual (- measured rho)
                                                 :weight (/ 1.0 (* sigma sigma))}))
                                            obs)))
                      final (loop [x (vec guess) k 0]
                              (if (>= k 4) x
                                  (recur (mapv + x (:correction (est/solve-batch (rows x) 6))) (inc k))))
                      rs*   (rows final)
                      cov   (:covariance (est/solve-batch rs* 6))]
                  {:n (count obs)
                   :stations (count (distinct (map :station obs)))
                   :rms (est/rms rs*)
                   :error (mag (mapv - (subvec final 0 3) (subvec truth 0 3)))
                   :formal (Math/sqrt (+ (nth (nth cov 0) 0) (nth (nth cov 1) 1) (nth (nth cov 2) 2)))}))
        weak   (solve 3.0 120.0 5)
        strong (solve 6.0 150.0 5)]
    (testing "the short arc barely sees the satellite"
      (is (< (:n weak) 15))
      (is (<= (:stations weak) 2) "only two stations get a pass"))
    (testing "yet it fits the data beautifully"
      (is (< (:rms weak) (* 2.0 sigma)) "residuals at the noise floor"))
    (testing "while being kilometres wrong"
      (is (> (:error weak) 1.0)))
    (testing "and the covariance is what says so"
      (is (> (:formal weak) 1.0) "formal uncertainty is kilometres too")
      (is (> (/ (:formal weak) (:rms weak)) 20.0)
          "formal uncertainty enormous next to the residual: the tell"))
    (testing "a longer arc determines the orbit and the covariance agrees"
      (is (> (:n strong) 25))
      (is (>= (:stations strong) 3))
      (is (< (:error strong) 0.1) "tens of metres")
      (is (< (/ (:formal strong) (:rms strong)) 20.0) "and no longer flagged"))))

;; -------------------------------------------------- orthogonal least squares

(defn- lauchli
  "A matrix of full rank for any eps > 0, whose square is singular once
  eps^2 falls below machine epsilon. The standard demonstration that
  forming A^T A destroys information the problem never lacked."
  [eps]
  [[1.0 1.0 1.0] [eps 0.0 0.0] [0.0 eps 0.0] [0.0 0.0 eps]])

(defn- lauchli-rows [eps]
  (let [A (lauchli eps)
        b (lin/mat-vec A [1.0 1.0 1.0])]
    (mapv (fn [a bi] {:H a :residual bi :weight 1.0}) A b)))

(deftest normal-equations-square-the-condition-number
  ;; Not a criticism of the implementation but of the method. cond(A^T A) is
  ;; cond(A) squared, so a problem comfortably solvable in double precision
  ;; becomes singular purely from being written down that way.
  (testing "fine while the square is still conditioned"
    (is (some? (est/solve (lauchli-rows 1e-6) 3 {:method :normal}))))
  (testing "losing digits as the square approaches machine epsilon"
    (let [x (:correction (est/solve (lauchli-rows 1e-7) 3 {:method :normal}))]
      (is (some? x))
      (is (> (apply max (map (fn [a] (abs (- a 1.0))) x)) 1e-4)
          "two digits gone, with nothing to indicate it")))
  (testing "and failing outright beyond it"
    (is (nil? (est/solve (lauchli-rows 1e-8) 3 {:method :normal})))))

(deftest orthogonal-reduction-does-not
  ;; Householder never forms A^T A, so it works wherever A itself is
  ;; solvable -- which at eps = 1e-10 means cond(A) = 1e10, unremarkable.
  (doseq [eps [1e-6 1e-7 1e-8 1e-10]]
    (let [x (:correction (est/solve (lauchli-rows eps) 3 {:method :qr}))]
      (is (some? x) (str "eps " eps))
      (doseq [xi x]
        (is (close? 1.0 xi 1e-12) (str "eps " eps " gave " x))))))

(deftest both-methods-agree-where-both-work
  ;; The point of having two: on a well-conditioned problem they must be
  ;; interchangeable, in the solution and in the covariance.
  (let [g (java.util.Random. 3)
        rows (mapv (fn [_] {:H (mapv (fn [_] (.nextGaussian g)) (range 4))
                            :residual (.nextGaussian g)
                            :weight (+ 0.5 (abs (.nextGaussian g)))})
                   (range 40))
        n (est/solve rows 4 {:method :normal})
        q (est/solve rows 4 {:method :qr})]
    (doseq [[a b] (map vector (:correction n) (:correction q))]
      (is (close? a b 1e-12) "same solution"))
    (doseq [i (range 4) j (range 4)]
      (is (close? (nth (nth (:covariance n) i) j) (nth (nth (:covariance q) i) j) 1e-12)
          "and the same covariance, computed without the normal matrix"))))

(deftest qr-refuses-what-it-cannot-determine
  (testing "rank deficiency is judged against the size of the factor"
    ;; An absolute floor would pass a diagonal of 1e-15 beside a norm of 1,
    ;; and return a covariance of 1e30 with a straight face.
    (is (nil? (est/solve [{:H [1.0 1.0] :residual 1.0}
                          {:H [2.0 2.0] :residual 2.0}
                          {:H [3.0 3.0] :residual 3.0}] 2))
        "duplicate columns determine nothing"))
  (testing "and there must be at least as many observations as unknowns"
    (is (nil? (est/solve [{:H [1.0 0.0] :residual 1.0}] 2)))))

(deftest qr-honours-weights
  (let [rows [{:H [1.0] :residual 10.0 :weight 100.0}
              {:H [1.0] :residual 0.0 :weight 1.0}]]
    (is (close? (first (:correction (est/solve rows 1 {:method :normal})))
                (first (:correction (est/solve rows 1 {:method :qr})))
                1e-12)
        "weights enter as a row scaling by their square root")))

(deftest solve-defaults-to-the-safe-method
  ;; Because the cost of the cheaper one is silent: an ill-conditioned
  ;; problem does not announce itself.
  (is (some? (est/solve (lauchli-rows 1e-9) 3)) "the default copes")
  (is (nil? (est/solve (lauchli-rows 1e-9) 3 {:method :normal})) "the alternative does not"))

;; ------------------------------------------------------------------ planets

(deftest planetary-elements-give-the-right-orbits
  (doseq [[planet a-au ecc days] [[:mercury 0.3871 0.2056 87.97]
                                  [:venus   0.7233 0.0068 224.70]
                                  [:earth   1.0000 0.0167 365.26]
                                  [:mars    1.5237 0.0934 686.98]
                                  [:jupiter 5.2029 0.0484 4332.6]
                                  [:saturn  9.5367 0.0539 10759.2]
                                  [:uranus  19.1892 0.0473 30685.4]
                                  [:neptune 30.0699 0.0086 60189.0]]]
    (testing (name planet)
      (let [el (pl/elements-at planet c/mjd-J2000)]
        (is (close? a-au (/ (:a el) c/AU) 1e-4) "semi-major axis")
        (is (close? ecc (:e el) 1e-4) "eccentricity")
        ;; Kepler's third law from the fitted axis, against the observed period
        (is (< 0.999 (/ (pl/period planet c/mjd-J2000) days) 1.001) "period")))))

(deftest planets-stay-near-the-ecliptic
  ;; Every inclination is a couple of degrees or less, the solar system
  ;; being flat. A frame error would show up here first.
  (doseq [planet pl/order]
    (is (< (/ (:i (pl/elements-at planet c/mjd-J2000)) c/degrees) 7.1)
        (str (name planet) " inclination"))))

(deftest mars-comes-and-goes
  ;; Geocentric distance must swing between opposition and conjunction, which
  ;; only works if the Earth's own place comes from the same table.
  (let [ds (map #(/ (mag (pl/geocentric :mars (+ c/mjd-J2000 (* 10.0 %)))) c/AU) (range 0 1100))]
    (is (close? 0.37 (apply min ds) 0.02) "closest approach")
    (is (close? 2.68 (apply max ds) 0.03) "and the far side of the Sun")))

(deftest two-independent-ephemerides-agree-at-j2000
  ;; The M&G series and the Standish elements share no coefficients and no
  ;; derivation. Agreement is therefore a real check on both -- and on the
  ;; frame, since either being referred to a different equinox would show
  ;; immediately.
  (let [a (eph/sun c/mjd-J2000)
        b (pl/sun-from-earth c/mjd-J2000)
        sep (Math/acos (max -1.0 (min 1.0 (/ (reduce + (map * a b)) (* (mag a) (mag b))))))]
    (is (< (/ sep c/arcsec) 10.0) "within ten arcseconds at the epoch")
    (is (close? (mag a) (mag b) 5e3) "and within a few thousand km in distance")))

(deftest the-low-precision-sun-drifts-in-longitude-not-latitude
  ;; Recording a known limitation rather than leaving it to be rediscovered.
  ;; The M&G series has a mean-longitude rate 0.35 deg/century short, which
  ;; slides the Sun along the ecliptic at 12.7 arcsec a year. The frames are
  ;; fine: a rotation would move latitude too, and latitude stays put.
  (let [ecl (fn [v] (let [[x y z] v
                          ce (Math/cos eph/obliquity-J2000)
                          se (Math/sin eph/obliquity-J2000)]
                      [x (+ (* ce y) (* se z)) (- (* ce z) (* se y))]))
        gap (fn [yr]
              (let [mjd (t/calendar->mjd yr 6 15)
                    [ax ay az] (ecl (eph/sun mjd))
                    [bx by bz] (ecl (pl/sun-from-earth mjd))]
                {:lon (/ (- (Math/atan2 ay ax) (Math/atan2 by bx)) c/arcsec)
                 :lat (/ (- (Math/asin (/ az (mag [ax ay az])))
                            (Math/asin (/ bz (mag [bx by bz])))) c/arcsec)}))
        g2000 (gap 2000) g2040 (gap 2040)]
    (testing "longitude drifts secularly"
      (is (< (abs (:lon g2000)) 10.0))
      (is (> (abs (:lon g2040)) 300.0))
      (is (close? 12.7 (/ (- (abs (:lon g2040)) (abs (:lon g2000))) 40.0) 2.0)
          "about 12.7 arcsec a year"))
    (testing "latitude does not, which is what rules out a frame error"
      (is (< (abs (:lat g2040)) 30.0)))))

(deftest everything-shares-one-frame
  ;; The property that lets a satellite orbit, the Moon, the Sun and the
  ;; planets be drawn in one picture without further rotation.
  (let [mjd (t/calendar->mjd 2025 3 20)]
    (testing "the Sun's declination is bounded by the J2000 obliquity"
      (let [decl (fn [v] (/ (Math/asin (/ (nth v 2) (mag v))) c/degrees))]
        (doseq [d (range 0 365 5)]
          (is (< (abs (decl (eph/sun (+ mjd d)))) 23.5))
          (is (< (abs (decl (pl/sun-from-earth (+ mjd d)))) 23.5)))))
    (testing "and the planetary Sun agrees with the Earth row negated"
      (is (< (mag (mapv + (pl/sun-from-earth mjd) (pl/heliocentric :earth mjd))) 1e-6)))))

(deftest bodies-from-different-sources-can-simply-be-added
  ;; What the shared frame buys, and the one thing the solar system demo
  ;; depends on. The planets come from Vallado's elements and the Moon from
  ;; an unrelated analytical series; placing the Moon heliocentrically is a
  ;; vector addition only because neither needs rotating first.
  (doseq [mjd [c/mjd-J2000 (t/calendar->mjd 2025 6 1) (t/calendar->mjd 1995 2 14)]]
    (let [earth      (pl/heliocentric :earth mjd)
          moon-geo   (eph/moon mjd)
          moon-helio (mapv + earth moon-geo)]
      (is (close? (mag moon-geo) (mag (mapv - moon-helio earth)) 1e-6)
          "the sum minus the Earth gives the lunar distance back")
      (is (< 356000.0 (mag moon-geo) 407000.0) "which is a real lunar distance")
      (testing "and the Moon is never further from the Sun than Earth plus its own orbit"
        (is (< (abs (- (mag moon-helio) (mag earth))) 410000.0))))))

(deftest orbit-sampling-by-mean-anomaly-covers-the-whole-orbit
  ;; The demo traces each orbit by stepping mean anomaly rather than time, so
  ;; Neptune costs no more samples than Mercury. It has to reach both apsides.
  (doseq [planet [:mercury :earth :mars :neptune]]
    (let [el (pl/elements-at planet c/mjd-J2000)
          rs (for [i (range 129)]
               (let [nu (kep/mean->true (* c/two-pi (/ i 128.0)) (:e el))]
                 (mag (first (kep/elements->state c/GM-sun (assoc el :nu nu :M nil))))))]
      (is (close? (* (:a el) (- 1.0 (:e el))) (apply min rs) 1.0) (str (name planet) " periapsis"))
      (is (close? (* (:a el) (+ 1.0 (:e el))) (apply max rs) 1.0) (str (name planet) " apoapsis")))))

(deftest calendar-dates-round-trip
  ;; The solar system demo shows a date, which needs the inverse of the
  ;; conversion everything else uses.
  (doseq [[y m d h] [[2000 1 1 12.0] [1858 11 17 0.0] [1900 1 1 0.0]
                     [2024 2 29 6.0] [1582 10 15 0.0] [2100 12 31 23.5]]]
    (let [[y2 m2 d2 h2] (t/mjd->calendar (t/calendar->mjd y m d h))]
      (is (= [y m d] [y2 m2 d2]) (str y "-" m "-" d))
      (is (close? h h2 1e-6) "and the hour")))
  (testing "across a hundred years of consecutive days"
    (let [bad (count (for [n (range 0 36525 7)
                           :let [mjd (+ 40000.0 n)
                                 [y m d _] (t/mjd->calendar mjd)]
                           :when (> (abs (- mjd (t/calendar->mjd y m d 0.0))) 1e-9)]
                       n))]
      (is (zero? bad)))))

;; ------------------------------------------------------------- ocean tides

(deftest ocean-tides-are-a-tenth-of-the-solid-earth-tide
  (let [solid (tid/corrections [[c/GM-moon moon-at]])
        ocean (tid/ocean-corrections [[c/GM-moon moon-at]])]
    (doseq [k [[2 0] [2 2]]]
      (let [ratio (/ (get-in ocean [:C k]) (get-in solid [:C k]))]
        (is (< 0.05 ratio 0.20) (str "C" (first k) (second k) " ratio " ratio))))))

(deftest the-ocean-bulge-leads-the-body-raising-it
  ;; Not lags. The Earth turns once a day while the Moon takes twenty-seven,
  ;; so rotation drags the bulge past the sub-lunar point. The resulting
  ;; couple is why the Moon recedes 38 mm a year and the day lengthens.
  ;; Reverse this sign and the Moon would be spiralling in.
  (let [axis (fn [corr]
               (let [c22 (get-in corr [:C [2 2]])
                     s22 (get-in corr [:S [2 2]])]
                 (/ (* 0.5 (Math/atan2 s22 c22)) c/degrees)))]
    (is (close? 0.0 (axis (tid/corrections [[c/GM-moon moon-at]])) 1e-9)
        "the solid tide has no lead in this model")
    (is (close? 3.0 (axis (tid/ocean-corrections [[c/GM-moon moon-at]])) 1e-6)
        "the ocean bulge sits ahead, in the direction the Earth turns")
    (is (pos? (axis (tid/ocean-corrections [[c/GM-moon moon-at]])))
        "positive: ahead, not behind")))

(deftest the-lead-angle-is-what-it-is-set-to
  (doseq [lead [0.0 1.0 5.0 10.0]]
    (let [corr (tid/ocean-corrections [[c/GM-moon moon-at]] tid/k-ocean (* lead c/degrees))
          axis (/ (* 0.5 (Math/atan2 (get-in corr [:S [2 2]]) (get-in corr [:C [2 2]]))) c/degrees)]
      (is (close? lead axis 1e-6) (str "lead of " lead " degrees")))))

(deftest ocean-tides-add-to-the-acceleration
  (let [j2 {:GM c/GM-earth :R c/R-earth :normalised? true
            :C {[0 0] 1.0 [2 0] -4.841654e-4} :S {}}
        bodies [[c/GM-moon moon-at] [c/GM-sun [c/AU 0.0 0.0]]]
        r [7000.0 0.0 0.0]
        solid (tid/perturb j2 bodies)
        both  (tid/perturb-with-ocean j2 bodies)
        a-solid (mag (mapv - (geo/acceleration solid r 2) (geo/acceleration j2 r 2)))
        a-ocean (mag (mapv - (geo/acceleration both r 2) (geo/acceleration solid r 2)))]
    (is (pos? a-ocean) "the oceans contribute something")
    (is (< 0.05 (/ a-ocean a-solid) 0.20) "about a tenth of the solid tide")
    (is (< 1e-9 (* 1000.0 a-ocean) 1e-6) "and lands in the published decade")))

(deftest no-bodies-no-ocean-tide
  (let [corr (tid/ocean-corrections [])]
    (is (every? zero? (vals (:C corr))))
    (is (every? zero? (vals (:S corr)))))
  (let [f (tid/perturb-with-ocean geo/earth [])
        r [7000.0 1000.0 2000.0]]
    (is (close? 0.0 (mag (mapv - (geo/acceleration f r 4) (geo/acceleration geo/earth r 4))) 1e-30))))

;; ----------------------------------------------------------- carrier phase

(deftest carrier-wavelengths
  (is (close? 0.1903 (* 1000.0 (obs/wavelength obs/L1)) 1e-4) "L1 is 19.03 cm")
  (is (close? 0.2442 (* 1000.0 (obs/wavelength obs/L2)) 1e-4) "L2 is 24.42 cm")
  (is (close? 0.8619 (* 1000.0 (obs/wavelength (- obs/L1 obs/L2))) 1e-4)
      "widelane is 86 cm, which is why ambiguity resolution starts there")
  (is (close? 0.1070 (* 1000.0 (obs/wavelength (+ obs/L1 obs/L2))) 1e-4) "narrowlane is 11 cm"))

(deftest the-ionosphere-delays-code-and-advances-carrier
  ;; Group and phase velocities move opposite ways in a dispersive medium,
  ;; so a pseudorange comes out too long and a phase range too short by
  ;; exactly the same amount. That opposition is what makes the ionosphere
  ;; measurable rather than merely a nuisance.
  (doseq [tec [10 50 100]]
    (let [g 20000.0
          code (:range (obs/modelled-range {:geometric g :tec tec :frequency obs/L1}))
          ph   (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L1}))]
      (is (pos? (- code g)) "code long")
      (is (neg? (- ph g)) "carrier short")
      (is (close? 0.0 (+ (- code g) (- ph g)) 1e-15)
          "and by the same amount, to the last bit"))))

(deftest code-minus-phase-is-twice-the-ionosphere
  ;; The relation that makes cycle-slip detection work: once the ambiguity is
  ;; removed what is left is purely ionospheric and varies smoothly.
  (doseq [tec [10 50 100]]
    (let [g 20000.0 n 123456
          code (:range (obs/modelled-range {:geometric g :tec tec :frequency obs/L1}))
          ph   (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L1 :ambiguity n}))]
      ;; Tolerance set by cancellation, not by the model: this subtracts two
      ;; numbers near 20,000 km to get 0.016, so six digits go immediately
      ;; and 2e-12 km is the floor a double can offer here.
      (is (close? (* 2.0 (obs/ionospheric-delay tec obs/L1))
                  (+ (- code ph) (* (obs/wavelength obs/L1) n))
                  1e-11)
          (str tec " TECU")))))

(deftest geometry-free-keeps-only-the-ionosphere
  ;; Range, clocks and troposphere are all frequency-independent, so they
  ;; cancel in the difference. Nothing about where the satellite is survives.
  (let [tec 40
        gf (fn [g clk trop]
             (obs/geometry-free
              (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L1
                                          :station-clock clk :elevation trop}))
              (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L2
                                          :station-clock clk :elevation trop}))))
        base (gf 20000.0 0.0 nil)]
    (is (close? base (gf 25000.0 1e-6 (* 10.0 c/degrees)) 1e-12) "geometry and clock gone")
    (is (close? base (gf 30000.0 -3e-6 (* 40.0 c/degrees)) 1e-12) "troposphere gone too")
    (testing "and what remains really is the ionosphere"
      (is (close? (- (obs/ionospheric-delay tec obs/L2) (obs/ionospheric-delay tec obs/L1))
                  base 1e-12)))))

(deftest melbourne-wubbena-recovers-the-widelane-integer
  ;; Geometry cancels because both combinations carry it identically;
  ;; the ionosphere cancels because it enters them with opposite signs.
  ;; What survives is an integer, observable with no orbit or clock
  ;; knowledge whatever -- which is why ambiguity resolution can start
  ;; before an orbit is known at all.
  (let [n1 100000 n2 77000
        widelane-n (- n1 n2)
        mw (fn [g tec]
             (obs/melbourne-wubbena
              (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L1 :ambiguity n1}))
              (:range (obs/modelled-range {:geometric g :tec tec :frequency obs/L1}))
              obs/L1
              (:phase (obs/carrier-phase {:geometric g :tec tec :frequency obs/L2 :ambiguity n2}))
              (:range (obs/modelled-range {:geometric g :tec tec :frequency obs/L2}))
              obs/L2))]
    (doseq [[g tec] [[20000.0 10] [25000.0 80] [30000.0 150] [22000.0 0]]]
      (is (close? widelane-n (/ (mw g tec) (obs/wavelength (- obs/L1 obs/L2))) 1e-6)
          (str "geometry " g ", " tec " TECU recovers " widelane-n " cycles")))))

(deftest cycle-slips-show-in-the-geometry-free-combination
  ;; It changes only when the receiver loses count, so a step is a slip and
  ;; nothing else.
  (let [smooth  (obs/geometry-free 20000.0 19995.8)
        slipped (obs/geometry-free (+ 20000.0 (obs/wavelength obs/L1)) 19995.8)]
    (is (not (obs/cycle-slip? smooth (+ smooth 1e-6))) "quiet epochs are quiet")
    (is (obs/cycle-slip? smooth slipped)
        "and one lost cycle on L1 -- 19 cm -- must be visible, which is the
         whole purpose; a threshold in metres would report a clean series")
    (testing "the default threshold sits under a single wavelength"
      (is (< 5e-5 (obs/wavelength obs/L1))))))

(deftest an-ambiguity-is-a-whole-number-of-wavelengths
  (doseq [n [0 1 -5 123456]]
    (let [p (obs/carrier-phase {:geometric 20000.0 :frequency obs/L1 :ambiguity n})]
      (is (close? (* n (obs/wavelength obs/L1)) (:ambiguity p) 1e-12))
      (is (close? (+ 20000.0 (* n (obs/wavelength obs/L1))) (:phase p) 1e-9)
          "and it enters the phase range directly"))))

(deftest force-model-registry-test
  (let [mjd 58000.0
        r [7000.0 1000.0 200.0]
        v [1.0 7.0 0.5]
        mag (fn [x] (Math/sqrt (reduce + (map * x x))))]

    (testing "a breakdown sums to the acceleration it breaks down"
      ;; `acceleration` and `breakdown` used to carry separate lists of the
      ;; forces and had drifted: tides were applied by one and missing from
      ;; the other, so with tides on the parts did not add up to the whole.
      (doseq [config [{}
                      {:tides? true}
                      {:sun? false :moon? false}
                      {:drag {:area-to-mass 0.01 :cd 2.2}}
                      {:relativity? true}
                      {:srp {:area-to-mass 0.02 :cr 1.3}}
                      {:tides? true :relativity? true
                       :srp {:area-to-mass 0.02 :cr 1.3}
                       :drag {:area-to-mass 0.01 :cd 2.2}}]]
        (let [total (forces/acceleration config mjd r v)
              parts (forces/breakdown config mjd r v)
              summed (reduce (fn [acc [_ a]] (mapv + acc a)) [0.0 0.0 0.0] parts)]
          (is (< (mag (mapv - total summed)) 1e-18)
              (str "parts do not sum to the whole for " config)))))

    (testing "tides change the acceleration, and are visible in the breakdown"
      ;; Folded into the harmonics, because a tide is a perturbation of the
      ;; geopotential. What matters is that they are not lost.
      (let [plain (forces/acceleration {} mjd r v)
            tidal (forces/acceleration {:tides? true} mjd r v)]
        (is (pos? (mag (mapv - plain tidal))))
        (is (pos? (mag (mapv - (:harmonics (forces/breakdown {} mjd r v))
                             (:harmonics (forces/breakdown {:tides? true} mjd r v))))))))

    (testing "a config switches models on and off"
      (is (= #{:two-body :harmonics :sun :moon}
             (set (keys (forces/breakdown {} mjd r v)))))
      (is (= #{:two-body :harmonics}
             (set (keys (forces/breakdown {:sun? false :moon? false} mjd r v)))))
      (is (contains? (forces/breakdown {:drag {:area-to-mass 0.01 :cd 2.2}} mjd r v)
                     :drag)))

    (testing "whether velocity is needed comes from the same list"
      ;; Not from a second hand-kept list of which forces those are.
      (is (not (forces/velocity-dependent? {})))
      (is (not (forces/velocity-dependent? {:srp {:area-to-mass 0.02 :cr 1.3}})))
      (is (forces/velocity-dependent? {:drag {:area-to-mass 0.01 :cd 2.2}}))
      (is (forces/velocity-dependent? {:relativity? true}))
      (is (= (set (map :name (filter :needs-velocity? forces/force-models)))
             #{:drag :relativity})))))
