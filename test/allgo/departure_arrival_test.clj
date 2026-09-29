(ns allgo.departure-arrival-test
  "Departure and arrival in three dimensions: the hyperbolas that leave a
  parking orbit along a given asymptote and arrive into a capture orbit,
  checked against the B-plane and the elements computed independently,
  and the launch azimuth of NASA/TM-2010-216764's equation 3 against the
  plane through site and asymptote built from vectors."
  (:require [allgo.astro.bplane :as bplane]
            [allgo.astro.constants :as c]
            [allgo.astro.interplanetary :as ip]
            [allgo.astro.launch :as launch]
            [allgo.astro.rotation :as rotation]
            [allgo.astro.time :as time]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu-earth c/GM-earth)
(def ^:private mu-mars (c/GM-planet :mars))
(def ^:private z [0.0 0.0 1.0])

(defn- close? [a b tol] (every? true? (map #(< (abs (- %1 %2)) tol) a b)))
(defn- unit [v] (v3/normalize v))

(defn- radec [ra dec] [(* (math/cos dec) (math/cos ra)) (* (math/cos dec) (math/sin ra)) (math/sin dec)])

(def ^:private v-inf
  "The 2026 type II least-C3 departure: 3.02 km/s toward RLA 130.4, DLA
  23.2 degrees (EME2000)."
  (v3/scale (radec (math/to-radians 130.4) (math/to-radians 23.2)) 3.024))

(deftest the-planes
  (doseq [i-deg [28.5 51.6 90.0 150.0]]
    (testing (str "inclination " i-deg)
      (let [i (math/to-radians i-deg)
            hs (ip/planes v-inf i z)]
        (is (= 2 (count hs)))
        (doseq [h hs]
          (is (< (abs (v3/dot h (unit v-inf))) 1e-12) "the plane holds the asymptote")
          (is (< (abs (- (v3/dot h z) (math/cos i))) 1e-12) "at the inclination asked")))))
  (testing "at just the asymptote's declination the two planes are one"
    (let [hs (ip/planes v-inf (math/asin (v3/dot (unit v-inf) z)) z)]
      (is (seq hs))
      (is (every? #(close? % (first hs) 1e-6) hs))))
  (testing "no orbit less inclined than the asymptote's declination holds it"
    (is (empty? (ip/planes v-inf (math/to-radians 20.0) z)))))

(deftest leaving
  (let [rp (+ c/R-earth 407.0)                   ; the handbook's parking orbit
        orbits (ip/departure-orbits mu-earth rp v-inf (math/to-radians 28.5))]
    (is (= 2 (count orbits)))
    (doseq [{:keys [h r v dv]} orbits]
      (testing "a periapsis on the parking orbit, the burn along its motion"
        (is (< (abs (- (v3/length r) rp)) 1e-9))
        (is (< (abs (v3/dot r v)) 1e-6))
        (is (close? (unit (v3/cross r v)) h 1e-12))
        (is (< (abs (- dv (:dv (ip/departure mu-earth rp (v3/length v-inf))))) 1e-12)))
      (testing "the hyperbola's outgoing asymptote, from its elements, is the excess velocity"
        (let [hv (v3/cross r v)
              ev (v3/scale (v3/sub (v3/scale r (- (v3/dot v v) (/ mu-earth (v3/length r))))
                                   (v3/scale v (v3/dot r v)))
                           (/ 1.0 mu-earth))
              e (v3/length ev)
              nu (math/acos (/ -1.0 e))
              s (v3/add (v3/scale (unit ev) (math/cos nu)) (v3/scale (unit (v3/cross hv ev)) (math/sin nu)))
              vinf (math/sqrt (- (v3/dot v v) (/ (* 2.0 mu-earth) (v3/length r))))]
          (is (close? s (unit v-inf) 1e-12))
          (is (< (abs (- vinf (v3/length v-inf))) 1e-12))))
      (testing "and flown out for a year it is heading that way, as nearly
                as the Earth's pull, falling off as 1/r, lets it"
        (let [[_ v1] (u/propagate mu-earth r v (* 365.25 86400.0))]
          (is (close? v1 v-inf 2e-3)))))))

(deftest arriving
  (let [mjd (time/calendar->mjd 2027 8 19)
        k (rotation/pole :mars mjd)
        vinf (v3/scale (unit [0.3 -1.0 0.4]) 2.729)
        rp (+ (rotation/radius :mars) 400.0)
        ra 33000.0
        orbits (ip/arrival-orbits mu-mars rp ra vinf (math/to-radians 75.0) k)]
    (is (= 2 (count orbits)))
    (doseq [{:keys [h r v dv bt br]} orbits]
      (testing "an orbit of the inclination asked about Mars's pole"
        (is (< (abs (- (v3/dot h k) (math/cos (math/to-radians 75.0)))) 1e-12))
        (is (close? (unit (v3/cross r v)) h 1e-12)))
      (testing "the approach's B-plane, computed afresh from the periapsis state"
        (let [bp (bplane/b-plane mu-mars [r v] k)]
          (is (close? (:s bp) (unit vinf) 1e-12))
          (is (< (abs (- (:v-inf bp) (v3/length vinf))) 1e-12))
          (is (< (abs (- (bplane/periapsis-radius mu-mars (v3/length vinf) (math/hypot bt br)) rp)) 1e-6))
          (is (close? (flatten (bplane/periapsis-state mu-mars vinf bt br k)) (concat r v) 1e-9)
              "and it inverts")))
      (testing "the burn leaves the capture orbit"
        (let [v1 (v3/scale v (/ (- (v3/length v) dv) (v3/length v)))
              a (/ 1.0 (- (/ 2.0 rp) (/ (v3/dot v1 v1) mu-mars)))]
          (is (< (abs (- (- (* 2.0 a) rp) ra)) 1e-6)))))))

(deftest azimuth
  (testing "equation 3 against the plane through site and asymptote"
    (doseq [lat [28.45 -5.2 45.9 62.9] ra-site [0.0 97.0 200.0 310.0]
            rla [10.0 130.7 252.9] dla [-54.9 1.0 23.16 49.5]]
      (let [[lat ra-site rla dla] (map math/to-radians [lat ra-site rla dla])
            u (radec ra-site lat)
            s (radec rla dla)
            ;; heading from the site toward the asymptote along the circle
            d (v3/sub s (v3/scale u (v3/dot u s)))
            east (unit (v3/cross z u))
            north (v3/cross u east)
            az (launch/azimuth lat ra-site rla dla)
            expected (mod (math/atan2 (v3/dot d east) (v3/dot d north)) (* 2.0 math/PI))]
        (is (< (abs (- az expected)) 1e-12))
        (testing "and the inclination is the plane's"
          (is (< (abs (- (launch/inclination lat az)
                         (math/acos (v3/dot (unit (v3/cross u d)) z))))
                 1e-12))))))
  (testing "due east from the Cape: the latitude is the inclination"
    (is (< (abs (- (launch/inclination (math/to-radians 28.45) (/ math/PI 2.0)) (math/to-radians 28.45))) 1e-15))))

(deftest launch-windows
  (testing "a day at the Cape, the 2026 type II asymptote, the range's 40 to 115 degrees"
    (let [lat (math/to-radians 28.45) lon (math/to-radians -80.6)
          day (time/calendar->mjd 2026 10 31)
          [rla dla] (ip/asymptote v-inf day)
          limits [(math/to-radians 40.0) (math/to-radians 115.0)]
          ws (launch/windows lat lon rla dla day (inc day) limits)]
      (is (seq ws))
      (doseq [{:keys [heading open close]} ws]
        (is (< open close))
        (testing "open within the limits"
          (is (heading (launch/headings lat lon rla dla (* 0.5 (+ open close)) limits))))
        (testing "each edge where the azimuth reaches a limit, or the day's"
          (doseq [t [open close] :when (< day t (inc day))]
            (let [az (launch/azimuth lat (launch/site-right-ascension lon t) rla dla)
                  az (if (= heading :away) (mod (+ az math/PI) (* 2.0 math/PI)) az)]
              (is (some #(< (abs (- az %)) 1e-3) limits)))))))))
