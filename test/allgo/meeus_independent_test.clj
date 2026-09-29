(ns allgo.meeus-independent-test
  "Meeus's series against sources that owe nothing to Meeus
  (test/data/meeus): JPL Horizons for Pluto and the satellites of Jupiter
  and Saturn, the US Naval Observatory for the seasons and the Moon's
  phases, NASA's eclipse catalog for the eclipses -- and the full VSOP87
  Sun for the seasons' periodic terms. Each to the accuracy Meeus gives
  for his method, not to the reference's own."
  (:require [allgo.astro.eclipse :as ecl]
            [allgo.astro.jupiter-moons :as jm]
            [allgo.astro.lunar-events :as ev]
            [allgo.astro.planet-orbits :as po]
            [allgo.astro.saturn-moons :as sm]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as t]
            [allgo.geometry.vec3 :as v3]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(defn- lines [file] (str/split-lines (slurp (io/resource (str "data/meeus/" file)))))
(defn- fields [l] (str/split (str/trim l) #"\s+"))
(defn- jd->mjd [jd] (- jd 2400000.5))
(def ^:private arcsec (/ math/PI 648000.0))

(deftest pluto
  (testing "chapter 37 against DE441, 1890-2099: within 3.1 arcsec -- the theory was fitted to an older
            ephemeris, and the difference grows smoothly from a few hundredths near 2000"
    (doseq [[jd x y z] (map #(mapv parse-double (fields %)) (lines "horizons-pluto.txt"))]
      (let [[l b r] (po/pluto (jd->mjd jd))
            ours [(* r (math/cos b) (math/cos l)) (* r (math/cos b) (math/sin l)) (* r (math/sin b))]]
        (is (< (/ (v3/angle ours [x y z]) arcsec) 3.1) (str jd))))))

(defn- satellite-errors
  "For each satellite, the worst distance, km, between Horizons' position
  and ours scaled by the planet radius that fits best, and that radius."
  [file bodies positions]
  (into {}
        (for [[body k] bodies
              :let [pts (for [[b jd & xyz] (map fields (lines file)) :when (= b body)]
                          [(mapv parse-double xyz) ((positions (jd->mjd (parse-double jd))) k)])
                    radius (/ (reduce + (map (fn [[h o]] (/ (v3/length h) (v3/length o))) pts)) (count pts))]]
          [k {:radius radius :worst (apply max (map (fn [[h o]] (v3/distance h (v3/scale o radius))) pts))}])))

(deftest galilean-satellites
  (testing "chapter 44's E5 theory against Horizons, 1950-2050: a Jupiter radius near 71,420 km, and every
            satellite within 2,000 km"
    (doseq [[k {:keys [radius worst]}] (satellite-errors "horizons-jupiter-moons.txt"
                                                         (map vector ["501" "502" "503" "504"] jm/names)
                                                         (fn [mjd] (zipmap jm/names (jm/positions-3d mjd))))]
      (is (< 71400.0 radius 71450.0) (str k))
      (is (< worst 2000.0) (str k " " worst)))))

(deftest saturnian-satellites
  (testing "chapter 46's theory against Horizons, 1950-2050: a few thousand km for the inner five, more for
            Titan and Iapetus, and chaotic Hyperion some 26,000"
    (let [limits {:mimas 2500.0 :enceladus 1000.0 :tethys 1000.0 :dione 1000.0 :rhea 1000.0
                  :titan 4000.0 :hyperion 30000.0 :iapetus 12000.0}]
      (doseq [[k {:keys [worst]}] (satellite-errors "horizons-saturn-moons.txt"
                                                    (map vector (map #(str (+ 600 %)) (range 1 9)) sm/names)
                                                    sm/positions-3d)]
        (is (< worst (limits k)) (str k " " worst))))))

(defn- utc-mjd [y m d hm]
  (let [[hh mm] (map parse-long (str/split hm #":"))]
    (t/calendar->mjd (parse-long y) (parse-long m) (parse-long d) (+ hh (/ mm 60.0)))))

(deftest seasons
  (testing "chapter 27's periodic terms against the full VSOP87 Sun, 1951-2050: within Meeus's 51 seconds"
    (doseq [y (range 1951 2051) s [:march :june :september :december]]
      (is (< (* 86400.0 (abs (- (solar/season y s) (solar/season-exact y s)))) 51.0) (str y " " s))))
  (testing "and against the Naval Observatory's times, 1972-2040, to the minute they are given in, with the
            prediction of the Earth's rotation in future years: within 1.6 minutes"
    (doseq [[y m d hm] (map fields (lines "usno-seasons.txt"))]
      (let [s ({"3" :march "6" :june "9" :september "12" :december} m)]
        (is (< (* 1440.0 (abs (- (t/tt->utc (solar/season (parse-long y) s)) (utc-mjd y m d hm)))) 1.6)
            (str y " " s))))))

(deftest moon-phases
  (testing "chapter 49 against the Naval Observatory, every fourth year 1972-2040: within 1.1 minutes"
    (let [phase {"New-Moon" 0.0 "First-Quarter" 0.25 "Full-Moon" 0.5 "Last-Quarter" 0.75}]
      (doseq [[y m d hm p] (map fields (lines "usno-moon-phases.txt"))]
        (let [usno (utc-mjd y m d hm)
              ours (t/tt->utc (ev/moon-phase (+ 2000.0 (/ (- usno 51544.5) 365.25)) (phase p)))]
          (is (< (* 1440.0 (abs (- ours usno))) 1.1) (str y "-" m "-" d " " p)))))))

(def ^:private months {"Jan" 1 "Feb" 2 "Mar" 3 "Apr" 4 "May" 5 "Jun" 6 "Jul" 7 "Aug" 8 "Sep" 9 "Oct" 10 "Nov" 11 "Dec" 12})

(defn- nasa-eclipses [file]
  (for [l (lines file)
        :let [[_ y mo d hms _ _ _ ty _ g] (fields l)
              [hh mm ss] (map parse-long (str/split hms #":"))]]
    {:mjd (t/calendar->mjd (parse-long y) (months mo) (parse-long d) (+ hh (/ mm 60.0) (/ ss 3600.0)))
     :type (subs ty 0 1) :gamma (parse-double g)}))

(defn- decimal-year [mjd] (+ 2000.0 (/ (- mjd 51544.5) 365.2422)))

(deftest eclipses
  (testing "chapter 54 against NASA's catalog, 2001-2100: every solar eclipse, of the same type, its greatest
            eclipse within 1.5 minutes and gamma within 0.005"
    (doseq [{:keys [mjd type gamma]} (nasa-eclipses "nasa-solar-eclipses.txt")]
      (let [ours (ecl/solar (decimal-year mjd))]
        (is (some? ours) (str mjd))
        (is (= ({"T" :total "A" :annular "H" :hybrid "P" :partial} type) (:type ours)) (str mjd))
        (is (< (* 1440.0 (abs (- (:mjd ours) mjd))) 1.5) (str mjd))
        (is (< (abs (- (:gamma ours) gamma)) 0.005) (str mjd)))))
  (testing "and every lunar eclipse but the two faintest penumbral ones, gamma near 1.57, which Meeus's
            method is too coarse to find; 2015 April 4's, NASA's total by an umbral magnitude of 1.0008,
            partial here"
    (let [nasa (nasa-eclipses "nasa-lunar-eclipses.txt")
          found (keep (fn [e] (when-let [o (ecl/lunar (decimal-year (:mjd e)))] [e o])) nasa)]
      (is (= 226 (count found)))
      (is (every? (fn [e] (or (some #(= e (first %)) found) (> (abs (:gamma e)) 1.57))) nasa))
      (doseq [[{:keys [mjd type gamma]} ours] found]
        (when-not (< (abs (- mjd 57116.5)) 1.0)
          (is (= ({"T" :total "P" :partial "N" :penumbral} type) (:type ours)) (str mjd)))
        (is (< (* 1440.0 (abs (- (:mjd ours) mjd))) 1.5) (str mjd))
        (is (< (abs (- (:gamma ours) gamma)) 0.01) (str mjd))))))
