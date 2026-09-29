(ns allgo.cio-iers-test
  "The IAU 2006/2000A series in `allgo.astro.cio-data`, taken from ERFA,
  against the IERS Conventions' own tables 5.2a, 5.2b and 5.2d
  (test/data/iers): every polynomial coefficient and every periodic
  amplitude the same. ERFA leaves out the zero amplitudes the tables
  print beside a nonzero partner, and only those."
  (:require [allgo.astro.cio-data :as d]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(defn- table [file] (str/split-lines (slurp (io/resource (str "data/iers/" file)))))

(defn- periodic
  "{[power multipliers :sin/:cos] microarcseconds} of an IERS table."
  [file]
  (loop [[l & more] (table file) j nil out {}]
    (cond
      (nil? l) out
      (re-find #"^\s*j = \d" l) (recur more (parse-long (second (re-find #"j = (\d)" l))) out)
      (and j (re-find #"^\s+\d+\s+-?\d" l))
      (let [[_ s c & ms] (str/split (str/trim l) #"\s+")
            ms (mapv parse-long ms)]
        (recur more j (assoc out [j ms :sin] (parse-double s) [j ms :cos] (parse-double c))))
      :else (recur more j out))))

(defn- polynomial
  "The polynomial part's coefficients, microarcseconds, T^0 up."
  [file]
  (let [line (str/replace (first (filter #(re-find #"t\^5" %) (table file))) #"t(\^\d)?" "")]
    (mapv (fn [[_ sign v]] (* (if (= sign "-") -1.0 1.0) (parse-double v)))
          (re-seq #"([+-]?)\s*(\d[\d.]*)" line))))

(defn- ours-xy [coord]
  (into {} (for [[ms terms] d/xy-terms [xy sc pow amp] terms :when (= xy coord)]
             [[pow ms (if (zero? sc) :sin :cos)] amp])))

(def ^:private ours-s
  (into {} (for [[pow ms sn cs] d/s-terms kv [[[pow ms :sin] (* 1e6 sn)] [[pow ms :cos] (* 1e6 cs)]]] kv)))

(defn- same? [iers ours]
  (and (every? (fn [[k v]] (if (contains? ours k) (< (abs (- v (ours k))) 1e-6) (zero? v))) iers)
       (every? (fn [[k v]] (or (contains? iers k) (zero? v))) ours)))

(deftest periodic-terms
  (testing "X: table 5.2a"
    (let [iers (periodic "tab5.2a.txt")]
      (is (= 1600 (count (filter (fn [[k _]] (= :sin (peek k))) iers))))
      (is (same? iers (ours-xy 0)))))
  (testing "Y: table 5.2b"
    (is (same? (periodic "tab5.2b.txt") (ours-xy 1))))
  (testing "s + XY/2: table 5.2d"
    (is (same? (periodic "tab5.2d.txt") ours-s))))

(deftest polynomials
  (testing "the polynomial parts, the tables' microarcseconds our arcseconds"
    (doseq [[file ours] [["tab5.2a.txt" (d/xy-polynomial 0)] ["tab5.2b.txt" (d/xy-polynomial 1)]
                         ["tab5.2d.txt" d/s-polynomial]]]
      (let [iers (polynomial file)]
        (is (= 6 (count iers)) file)
        (is (every? true? (map #(< (abs (- (* 1e-6 %1) %2)) 1e-12) iers ours)) file)))))
