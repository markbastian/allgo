(ns allgo.eclipse-shadow-test
  "The Moon's shadow on the Earth (`eclipse/shadow` and `obscuration`)
  against NASA's catalog of the solar eclipses of 2001-2100
  (test/data/meeus): where each is greatest, from VSOP87 and ELP and the
  IAU's model of the Earth's rotation, at the catalog's own instant and
  delta T."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.eclipse :as ecl]
            [allgo.astro.moon :as lunar]
            [allgo.astro.rotation :as rot]
            [allgo.astro.time :as t]
            [allgo.astro.vsop87 :as vsop87]
            [allgo.geometry.vec3 :as v3]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(def ^:private months {"Jan" 1 "Feb" 2 "Mar" 3 "Apr" 4 "May" 5 "Jun" 6 "Jul" 7 "Aug" 8 "Sep" 9 "Oct" 10 "Nov" 11 "Dec" 12})

(defn- degrees
  "\"11S\" or \"3E\" as signed radians."
  [s]
  (math/to-radians (* (if (#{\S \W} (last s)) -1.0 1.0) (parse-double (subs s 0 (dec (count s)))))))

(def ^:private catalog
  (for [l (str/split-lines (slurp (io/resource "data/meeus/nasa-solar-eclipses.txt")))
        :let [[_ y mo d hms dt _ _ ty _ g mag lat lon] (str/split (str/trim l) #"\s+")
              [hh mm ss] (map parse-long (str/split hms #":"))]]
    {:mjd (t/calendar->mjd (parse-long y) (months mo) (parse-long d) (+ hh (/ mm 60.0) (/ ss 3600.0)))
     :delta-t (parse-double dt) :type (subs ty 0 1) :gamma (parse-double g)
     :magnitude (parse-double mag) :lat (degrees lat) :lon (degrees lon)}))

(defn- earth-fixed
  "The geocentric Moon and Sun, km, turned into the Earth-fixed frame."
  [mjd delta-t]
  (let [m (vsop87/ecliptic-of-date->J2000 mjd)
        r (rot/body->icrf :earth (- mjd (/ delta-t 86400.0)))
        fixed (fn [v] (mapv (fn [j] (reduce + (map #(* (get-in r [% j]) (v %)) (range 3)))) (range 3)))]
    [(fixed (lunar/geocentric-J2000 mjd m))
     (fixed (mapv #(* (- c/AU) %) (vsop87/equatorial-J2000 :earth mjd m)))]))

(defn- unit [lat lon]
  [(* (math/cos lat) (math/cos lon)) (* (math/cos lat) (math/sin lon)) (math/sin lat)])

(def ^:private eclipses
  (for [e catalog]
    (let [[moon sun] (earth-fixed (:mjd e) (:delta-t e))]
      (assoc e :moon moon :sun sun :ours (ecl/shadow moon sun)))))

(deftest gamma-test
  (testing "the axis's least distance from the center within 0.002 equatorial radii, for all 224"
    (is (= 224 (count eclipses)))
    (doseq [{:keys [mjd gamma ours]} eclipses]
      (is (< (abs (- gamma (:gamma ours))) 0.002) (str mjd))))
  (testing "and the penumbra reaching the Earth for every one, as Meeus's series has it"
    (doseq [{:keys [mjd gamma ours]} eclipses]
      (let [meeus (ecl/solar (+ 2000.0 (/ (- mjd 51544.5) 365.2422)))]
        (is (< (abs gamma) (+ 1.0 (:penumbra ours))) (str mjd))
        (is (< (abs (- (:penumbra meeus) (:penumbra ours))) 0.003) (str mjd))
        (is (< (abs (- (:u meeus) (:u ours))) 0.003) (str mjd))))))

(deftest greatest-eclipse-test
  (testing "the axis meets the Earth for every central eclipse and no partial one"
    (doseq [{:keys [mjd type gamma ours]} eclipses]
      (is (= (and (not= "P" type) (< (abs gamma) 0.9972)) (some? (:lat ours))) (str mjd))))
  (testing "and no shadow at all from a Moon beyond the Earth"
    (is (nil? (ecl/shadow [0.0 0.0 -384400.0] [0.0 0.0 1.496e8]))))
  (testing "the umbra for a total eclipse and the antumbra for an annular one"
    (doseq [{:keys [mjd type ours]} eclipses :when (and (:lat ours) (#{"T" "A"} type))]
      (is (= ({"T" :total "A" :annular} type) (:kind ours)) (str mjd))))
  (testing "where it is greatest within a degree of NASA's, which rounds to whole degrees -- but where
            the axis grazes the limb, a little error in it moves the point a long way"
    (doseq [{:keys [mjd gamma lat lon ours]} eclipses :when (and (:lat ours) (< (abs gamma) 0.95))]
      (is (< (math/to-degrees (v3/angle (unit lat lon) (unit (:lat ours) (:lon ours)))) 1.0) (str mjd)))))

(deftest obscuration-test
  (let [surface (fn [{:keys [lat lon]}] (v3/scale (unit lat lon) 6378.137))]
    (testing "at greatest eclipse the Sun is hidden entirely under the umbra, and under the antumbra as
              much as the square of NASA's magnitude, the ratio of the disks' diameters"
      (doseq [{:keys [mjd type magnitude moon sun ours]} eclipses :when (and (:lat ours) (#{"T" "A"} type))]
        (let [o (ecl/obscuration moon sun (surface ours))]
          (if (= type "T")
            (is (= 1.0 o) (str mjd))
            (is (< (abs (- o (* magnitude magnitude))) 0.004) (str mjd))))))
    (testing "partly a thousand km across the axis, and not at all beyond the penumbra's 3,500 or so"
      (doseq [{:keys [mjd moon sun ours]} eclipses :when (:lat ours)]
        (let [p (surface ours)
              across (v3/normalize (v3/cross (v3/sub moon sun) [0.0 0.0 1.0]))]
          (is (< 0.0 (ecl/obscuration moon sun (v3/add-scaled p across 1000.0)) 1.0) (str mjd))
          (is (zero? (ecl/obscuration moon sun (v3/add-scaled p across 5000.0))) (str mjd))))))
  (testing "two disks' overlap"
    (let [s [0.0 0.0 1.5e8] p [0.0 0.0 0.0]
          ;; a Moon as large as the Sun, offset by a fraction of its angular radius
          at (fn [k] (let [r (* 1737.4 (/ 1.5e8 695700.0))
                           a (math/asin (/ 695700.0 1.5e8))]
                       (v3/scale [(math/sin (* k a)) 0.0 (math/cos (* k a))] r)))]
      (is (< (abs (- 1.0 (ecl/obscuration (at 0.0) s p))) 1e-9))
      (is (zero? (ecl/obscuration (at 2.01) s p)))
      ;; offset by one radius: the lens is 2/3 - sqrt(3)/(2 pi) of a disk
      (is (< (abs (- (- (/ 2.0 3.0) (/ (math/sqrt 3.0) (* 2.0 math/PI))) (ecl/obscuration (at 1.0) s p))) 1e-6)))))
