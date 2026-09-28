(ns allgo.msis-test
  "NRLMSISE-00 against the output of its own Fortran distribution's test
  driver (test/data/msis): all seventeen cases, fifteen to the four
  figures of the summary table and all to the three of the per-case
  blocks, the 3-hour ap cases among them. Then what the model must satisfy
  whatever its coefficients: continuity where its layers meet, and the
  mass density the sum of its species'."
  (:require [allgo.astro.msis :as msis]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(def ^:private lines (str/split-lines (slurp (io/resource "data/msis/nrlmsise00_output.txt"))))

(defn- numbers [line] (mapv parse-double (re-seq #"[-+]?\d+(?:\.\d*)?(?:E[-+]\d+)?" line)))

(def ^:private summary
  "The fifteen daily-Ap cases of the summary tables: {label [values]}."
  (let [rows (drop-while #(not (str/starts-with? % " DAY")) lines)
        row (fn [label] (->> rows
                             (filter #(str/starts-with? % (str " " label)))
                             (mapcat #(numbers (subs % 6)))
                             vec))]
    (into {} (map (fn [l] [l (row l)])
                  ["DAY" "UT" "ALT" "LAT" "LONG" "LST" "F107A" "F107 " "AP"
                   "TINF" "TG" "HE" "O " "N2" "O2" "AR" "H " "N " "ANM O" "RHO"]))))

(def ^:private blocks
  "The seventeen per-case blocks: [[He O N2 O2 Ar rho H N] [anomalous-O Tinf T]]."
  (->> lines
       (take-while #(not (str/includes? % "MSISE-00")))
       (partition-by str/blank?)
       (remove #(str/blank? (first %)))
       (mapv (fn [[a b]] [(numbers a) (numbers b)]))))

(defn- inputs [k]
  (let [v (fn [l] ((summary l) k))]
    {:doy (long (v "DAY")) :sec (v "UT") :alt (v "ALT") :lat (v "LAT") :lon (v "LONG")
     :lst (v "LST") :f107a (v "F107A") :f107 (v "F107 ") :ap (v "AP")}))

(defn- close?
  "Within half a unit in the last of `figures` significant figures; a
  printed zero allows the single-precision Fortran's underflow."
  [printed ours figures]
  (if (zero? printed)
    (< (abs ours) 1e-30)
    (<= (abs (- ours printed))
        (* 0.5 (math/pow 10.0 (- (math/floor (math/log10 (abs printed))) (dec figures)))))))

(defn- in-cgs
  "The model's outputs in the Fortran's units: cm^-3 and g/cm^3."
  [r]
  (-> (reduce #(update %1 %2 * 1e-6) r [:He :O :N2 :O2 :Ar :H :N :anomalous-O])
      (update :rho * 1e-3)))

(deftest summary-table
  (testing "the fifteen daily-Ap cases, to the table's four figures"
    (is (= 15 (count (summary "DAY"))))
    (doseq [k (range 15)]
      (let [r (in-cgs (msis/atmosphere (inputs k)))
            printed (fn [l] ((summary l) k))]
        (is (<= (abs (- (:t-exo r) (printed "TINF"))) 0.005) (str "case " (inc k) " TINF"))
        (is (<= (abs (- (:t r) (printed "TG"))) 0.005) (str "case " (inc k) " TG"))
        (doseq [[l key] [["HE" :He] ["O " :O] ["N2" :N2] ["O2" :O2] ["AR" :Ar] ["H " :H] ["N " :N]
                         ["ANM O" :anomalous-O] ["RHO" :rho]]]
          (is (close? (printed l) (key r) 4) (str "case " (inc k) " " l ": " (printed l) " " (key r))))))))

(deftest every-case
  (testing "all seventeen, the 3-hour ap cases too, to the blocks' three figures"
    (is (= 17 (count blocks)))
    (let [ap3 (vec (repeat 7 100.0))
          all (concat (map inputs (range 15))
                      [(assoc (inputs 0) :ap ap3) (assoc (inputs 0) :ap ap3 :alt 100.0)])]
      (doseq [[k in [[he o n2 o2 ar rho h n] [anom tinf t]]] (map vector (range) all blocks)]
        (let [r (in-cgs (msis/atmosphere in))]
          (doseq [[printed key] [[he :He] [o :O] [n2 :N2] [o2 :O2] [ar :Ar] [rho :rho] [h :H] [n :N]
                                 [anom :anomalous-O]]]
            (is (close? printed (key r) 3) (str "case " (inc k) " " key ": " printed " " (key r))))
          (is (close? tinf (:t-exo r) 4) (str "case " (inc k) " Tinf"))
          (is (close? t (:t r) 4) (str "case " (inc k) " T")))))))

(deftest consistency
  (let [base (assoc (inputs 0) :lst 10.0)]
    (testing "continuous where the model's layers meet: the mixing at 62.5 km and the lower nodes at 32.5"
      (doseq [z [62.5 32.5]]
        (let [below (msis/atmosphere (assoc base :alt (- z 1e-6)))
              above (msis/atmosphere (assoc base :alt (+ z 1e-6)))]
          (doseq [k [:N2 :O2 :He :Ar :rho :t]]
            (is (< (abs (- (k below) (k above))) (* 1e-5 (abs (k above)))) (str z " km " k))))))
    (testing "and at 72.5 km, where the thermosphere begins, but for the step the Fortran takes there: the
              variation of the temperature profile's shape turns on only above it, some 0.3%"
      (let [below (msis/atmosphere (assoc base :alt (- 72.5 1e-6)))
            above (msis/atmosphere (assoc base :alt (+ 72.5 1e-6)))]
        (doseq [k [:N2 :O2 :He :Ar :rho :t]]
          (is (< (abs (- (k below) (k above))) (* 5e-3 (abs (k above)))) (str "72.5 km " k)))))
    (testing "the mass densities are their species' masses summed, anomalous oxygen in the drag density alone"
      (doseq [alt [100.0 300.0 600.0 900.0]]
        (let [{:keys [He O N2 O2 Ar H N anomalous-O rho rho-drag]} (msis/atmosphere (assoc base :alt alt))
              amu 1.66e-27]
          (is (< (abs (- rho (* amu (+ (* 4 He) (* 16 O) (* 28 N2) (* 32 O2) (* 40 Ar) H (* 14 N))))) (* 1e-12 rho)))
          (is (< (abs (- rho-drag rho (* amu 16 anomalous-O))) (* 1e-12 rho))))))
    (testing "the density falls with height"
      (let [rhos (map #(:rho (msis/atmosphere (assoc base :alt %))) (range 0.0 1000.0 10.0))]
        (is (apply > rhos))))))
