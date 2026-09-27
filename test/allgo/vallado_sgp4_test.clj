(ns allgo.vallado-sgp4-test
  "SGP4 against the verification set that comes with *Revisiting
  Spacetrack Report #3* (AIAA 2006-6753): every case of SGP4-VER.TLE, run
  as the package's C++ test driver runs it, compared with the output the
  package publishes -- Vallado's C++ and MATLAB in opsmode a, and the Java
  conversion it ships in opsmode i. test/data/sgp4/README.md has each
  file's provenance."
  (:require [allgo.astro.sgp4 :as sgp4]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(defn- data [f] (slurp (io/resource (str "data/sgp4/" f))))

(def ^:private cases
  "Each case of SGP4-VER.TLE, with the start, stop and step its line 2
  carries after column 69."
  (delay
    (->> (str/split-lines (data "SGP4-VER.TLE"))
         (remove #(or (str/blank? %) (str/starts-with? % "#")))
         (partition 2)
         (map (fn [[l1 l2]]
                (let [[start stop step] (map parse-double (str/split (str/trim (subs l2 69)) #"\s+"))]
                  {:l1 l1 :l2 l2 :satnum (parse-long (str/trim (subs l1 2 7)))
                   :start start :stop stop :step step}))))))

(defn- run
  "TestSGP4.cpp's verification loop: the state at the epoch, then from
  start to stop by step until the propagator reports an error.
  {:rows [[t r v] ...] :error [t code]}."
  [{:keys [l1 l2 start stop step]} opsmode]
  (let [rec (sgp4/twoline2rv l1 l2 {:opsmode opsmode})
        at0 (sgp4/sgp4 rec 0.0)
        failed? (pos? (:error at0))]
    (loop [t (if (> (abs start) 1.0e-8) (- start step) start)
           rows (if failed? [] [[0.0 (:r at0) (:v at0)]])
           error (when failed? [0.0 (:error at0)])]
      (if (or error (>= t stop))
        {:rows rows :error error}
        (let [t (min (+ t step) stop)
              {:keys [r v] code :error} (sgp4/sgp4 rec t)]
          (if (zero? code)
            (recur t (conj rows [t r v]) nil)
            (recur t rows [t code])))))))

(defn- numbers [line] (mapv parse-double (str/split (str/trim line) #"\s+")))

(defn- parse-ver
  "A verification .out file: [satnum rows], rows [t r v], in file order."
  [text]
  (->> (str/split-lines text)
       (remove str/blank?)
       (partition-by #(str/ends-with? (str/trim %) "xx"))
       (partition 2)
       (map (fn [[[header] lines]]
              [(parse-long (first (str/split (str/trim header) #"\s+")))
               (for [l lines :let [[t x y z vx vy vz] (numbers l)]]
                 [t [x y z] [vx vy vz]])]))))

(defn- parse-e
  "An STK .e file from the C++ driver: rows [t r v], t turned to minutes."
  [text]
  (->> (str/split-lines text)
       (drop-while #(not (str/includes? % "EphemerisTimePosVel")))
       rest
       (take-while #(not (str/includes? % "END Ephemeris")))
       (remove str/blank?)
       (map (fn [l] (let [[t x y z vx vy vz] (numbers l)] [(/ t 60.0) [x y z] [vx vy vz]])))))

(defn- worst*
  [published ours]
  (reduce (fn [[dr dv] [[_ r v] [_ r' v']]]
            [(apply max dr (map (comp abs -) r r'))
             (apply max dv (map (comp abs -) v v'))])
          [0.0 0.0]
          (map vector published ours)))

(defn- worst
  "The largest position (km) and velocity (km/s) differences between the
  published rows and a run's, having checked that the times match. Where
  a case fails at the epoch the drivers still print a row for it, of
  whatever vectors the previous case left, which is dropped."
  [published {:keys [rows error]}]
  (let [published (if (= 0.0 (first error)) (rest published) published)]
    (is (= (count published) (count rows)))
    (is (every? #(< (abs %) 1.0e-6) (map - (map first published) (map first rows))))
    (worst* published rows)))

;; Published positions are printed to 1e-8 km and velocities to 1e-9
;; km/s, and ours round to every digit of the C++ and MATLAB files: the
;; largest differences are half a unit in the last place, 5e-9 km and
;; 5e-10 km/s. (Against a build of the package's C++ at full precision
;; they agree to 3e-8 km at worst, 3.5 years out on the most eccentric
;; case, what is left being the last bits of the platform's sin, cos and
;; pow.)

(deftest cpp-verification
  (testing "every case against the C++ output in the package, opsmode a"
    (let [runs (into {} (map (juxt :satnum #(run % :a))) @cases)]
      (doseq [satnum (keys runs)
              :let [published (parse-e (data (format "cpp/%05d.e" satnum)))
                    [dr dv] (worst published (runs satnum))]]
        (is (< dr 1.0e-8) (str satnum " position " dr))
        (is (< dv 1.0e-9) (str satnum " velocity " dv))))))

(deftest matlab-verification
  (testing "every case against Vallado's MATLAB output, opsmode a, including the
            first run of 20413 that the C++'s .e files overwrite"
    ;; the file ends with a stray copy of 08195's row at 120 minutes
    (let [lines (str/split-lines (data "tmatverDec2015.out"))
          stray (take 7 (numbers (last lines)))
          row (first (filter #(= 120.0 (first %)) (second (nth (parse-ver (str/join "\n" lines)) 3))))]
      (is (= stray (flatten row))))
    (doseq [[c [satnum published]] (map vector @cases
                                        (parse-ver (str/join "\n" (butlast (str/split-lines (data "tmatverDec2015.out"))))))
            :let [[dr dv] (worst published (run c :a))]]
      (is (= satnum (:satnum c)))
      (is (< dr 1.0e-8) (str satnum " position " dr))
      (is (< dv 1.0e-9) (str satnum " velocity " dv)))))

(deftest java-verification
  (testing "every case against the package's Java conversion, opsmode i"
    (doseq [[c [satnum published]] (map vector @cases (parse-ver (data "java_sgp4_ver.out")))
            :let [[dr dv] (worst published (run c :i))]]
      (is (= satnum (:satnum c)))
      (is (< dr 1.0e-8) (str satnum " position " dr))
      (is (< dv 1.0e-9) (str satnum " velocity " dv)))))

(deftest errors
  ;; The published files stop where the propagator does. These are the
  ;; codes the package's C++ gives there; SGP4-VER.TLE's comments name
  ;; 33333's, a negative semilatus rectum, and 33334's attempt at a
  ;; negative mean motion, which the theory catches instead as the
  ;; perturbed eccentricity going out of range.
  (testing "where each failing case stops, and why"
    (is (= {22312 [494.2028672 1] 28350 [1560.0 1] 28872 [55.0 6] 29141 [440.0 6]
            33333 [25.0 4] 33334 [0.0 3] 20413 [1844345.0 6]}
           (into {} (for [c @cases :let [{:keys [error]} (run c :a)] :when error]
                      [(:satnum c) error]))))))

(deftest decay-still-has-a-state
  (testing "code 6 carries the position computed, inside the Earth"
    (let [c (first (filter #(= 28872 (:satnum %)) @cases))
          {:keys [r error]} (sgp4/sgp4 (sgp4/twoline2rv (:l1 c) (:l2 c)) 55.0)]
      (is (= 6 error))
      (is (< (math/sqrt (reduce + (map * r r))) 6378.135)))))

(deftest tle-epoch
  (testing "the epoch the C++ works out, split as it splits it"
    (let [c (first @cases)
          rec (sgp4/twoline2rv (:l1 c) (:l2 c))]
      (is (= 2451722.5 (:jdsatepoch rec)))
      (is (< (abs (- (:jdsatepochF rec) 0.78495062)) 1.0e-12))
      (is (= "00005" (:satnum rec))))))
