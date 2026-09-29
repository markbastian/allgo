(ns allgo.orekit-crosscheck-test
  "Models checked otherwise only against their own source's output --
  SGP4 against its reference C++, NRLMSISE-00 against NRL's driver,
  Knocke's coefficients against memory -- held against an independent
  implementation's (Orekit 12.2, output only; test/data/orekit)."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.drag :as drag]
            [allgo.astro.earth-radiation :as er]
            [allgo.astro.msis :as msis]
            [allgo.astro.sgp4 :as sgp4]
            [allgo.astro.time :as time]
            [allgo.geometry.vec3 :as v3]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(defn- rows [file]
  (for [l (str/split-lines (slurp (io/resource (str "data/orekit/" file))))]
    (str/split (str/trim l) #"\s+")))

(deftest knocke
  (testing "Knocke's zonal albedo and emissivity, at seven latitudes on six dates"
    (let [{:keys [albedo emissivity]} er/knocke]
      (doseq [[mjd lat a e] (map #(mapv parse-double %) (rows "knocke.txt"))]
        (is (< (abs (- (albedo (math/to-radians lat) mjd) a)) 1e-9) (str mjd " " lat))
        (is (< (abs (- (emissivity (math/to-radians lat) mjd) e)) 1e-9) (str mjd " " lat))))))

(def ^:private tle-lines
  (->> (str/split-lines (slurp (io/resource "data/sgp4/SGP4-VER.TLE")))
       (remove #(or (str/starts-with? % "#") (str/blank? %)))
       vec))

(def ^:private hard-cases
  "Where Orekit parts from the reference output that allgo.astro.sgp4
  matches (see test/data/orekit/README.md)."
  #{"29:33333" "31:33335" "32:20413"})

(deftest sgp4
  (let [recs (into {} (for [i (range 0 (count tle-lines) 2)]
                        [(str (quot i 2) ":" (subs (tle-lines i) 2 7))
                         (sgp4/twoline2rv (tle-lines i) (subs (tle-lines (inc i)) 0 69) {:opsmode :i})]))
        by-case (group-by first (rows "sgp4.txt"))
        points (for [[case rs] by-case
                     :let [rec (recs case)
                           ts (map #(parse-double (second %)) rs)
                           ;; each case up to its first error, as the verification driver runs it
                           stop (or (first (filter #(pos? (:error (sgp4/sgp4 rec %))) (sort ts))) ##Inf)]
                     [_ t & xs] rs
                     :let [t (parse-double t)]
                     :when (< t stop)
                     :let [[x y z vx vy vz] (map parse-double xs)
                           {:keys [r v]} (sgp4/sgp4 rec t)]]
                 {:case case :dr (v3/distance r [x y z]) :dv (v3/distance v [vx vy vz])})]
    (testing "the verification set, every case up to its first error: millimeters, but for the three hard cases"
      (let [ordinary (remove #(hard-cases (:case %)) points)]
        (is (= 509 (count ordinary)))
        (is (every? #(< (:dr %) 5e-6) ordinary))
        (is (every? #(< (:dv %) 5e-9) ordinary))))
    (testing "where Orekit parts from the reference, it does so by more than a meter -- the reference is what we match"
      (doseq [case hard-cases]
        (is (> (apply max (map :dr (filter #(= case (:case %)) points))) 1e-3) case)))))

(deftest nrlmsise00
  (testing "945 points, 50 to 1000 km, three dates, three levels of activity: the drag density to 1e-9"
    (doseq [[mjd lat lon alt f107a f107 ap _ _ rho] (map #(mapv parse-double %) (rows "nrlmsise00.txt"))]
      (let [[y] (time/mjd->calendar mjd)
            doy (inc (long (math/floor (- mjd (time/calendar->mjd y 1 1)))))
            sec (* 86400.0 (- mjd (math/floor mjd)))
            ours (:rho-drag (msis/atmosphere {:doy doy :sec sec :alt alt :lat lat :lon lon
                                              :lst (mod (+ (/ sec 3600.0) (/ lon 15.0)) 24.0)
                                              :f107a f107a :f107 f107 :ap ap}))]
        ;; 250 km is where the model stops mixing O2 across the turbopause;
        ;; Orekit's altitude, recomputed from a position, can land a hair
        ;; either side of it and take the other branch
        (is (< (abs (- (/ ours rho) 1.0)) (if (== alt 250.0) 1e-5 1e-9)) (str mjd " " lat " " lon " " alt))))))

(deftest harris-priester
  (testing "256 points, 110 to 990 km, four directions of the Sun: to 2e-4"
    (doseq [[x y z sx sy sz rho] (map #(mapv parse-double %) (rows "harris-priester.txt"))]
      (let [ours (drag/density [x y z] (mapv #(* % c/AU) [sx sy sz]))]
        (if (zero? rho)
          ;; above the model's 1000 km ceiling, as over the pole at 990 km above the equatorial radius
          (is (zero? ours) (str [x y z]))
          (is (< (abs (- (/ ours rho) 1.0)) 2e-4) (str [x y z])))))))
