(ns allgo.sgp4-fit-test
  "Fitting SGP4 elements to states: element sets from the verification
  set recovered from the states they generate, starting from nothing but
  the osculating elements, and a fit through noise that averages it out."
  (:require [allgo.astro.sgp4 :as sgp4]
            [allgo.astro.sgp4-fit :as fit]
            [allgo.geometry.vec3 :as v3]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(def ^:private lines
  (->> (str/split-lines (slurp (io/resource "data/sgp4/SGP4-VER.TLE")))
       (remove #(str/starts-with? % "#"))
       vec))

(defn- tle [satnum]
  (let [i (first (keep-indexed #(when (str/starts-with? %2 (str "1 " satnum)) %1) lines))]
    (sgp4/twoline2rv (lines i) (subs (lines (inc i)) 0 69))))

(defn- states [satrec times]
  (mapv (fn [t] (let [s (sgp4/sgp4 satrec t)] {:t t :r (:r s) :v (:v s)})) times))

(defn- template [truth] (select-keys truth [:epoch :ndot :nddot :bstar]))

(deftest equinoctial-round-trip
  (let [el {:no-kozai 0.06 :ecco 0.1 :inclo 1.0 :nodeo 2.0 :argpo 3.0 :mo 4.0}
        back (fit/equinoctial->elements (fit/elements->equinoctial el))]
    (doseq [k (keys el)] (is (< (abs (- (el k) (back k))) 1e-14) (str k)))))

(deftest recovers-element-sets
  (doseq [[satnum what bstar? span] [["00005" "eccentric, near-Earth" false 1440.0]
                                     ["06251" "near-circular with heavy drag, B* too" true 1440.0]
                                     ["09880" "Molniya, deep space" false 2880.0]]]
    (testing what
      (let [truth (tle satnum)
            {:keys [elements rms]} (fit/fit (cond-> (template truth) bstar? (assoc :bstar 0.0))
                                            (states truth (range 0.0 span 10.0))
                                            {:bstar? bstar?})]
        (is (< rms 1e-6))
        (doseq [k [:no-kozai :ecco :inclo :nodeo :argpo :mo]]
          (is (< (abs (- (elements k) (truth k))) 1e-12) (str satnum " " k)))
        (when bstar?
          (is (< (abs (- (:bstar elements) (:bstar truth))) 1e-12)))))))

(deftest averages-noise
  (testing "positions noisy by 100 m: the fit's residual is the noise, and its orbit far closer than that"
    (let [truth (tle "06251")
          times (range 0.0 1440.0 5.0)
          ;; a deterministic stand-in for noise, uniform on [-0.1, 0.1] km
          noise (fn [k] (* 0.2 (- (mod (* (inc k) 0.6180339887498949) 1.0) 0.5)))
          obs (map-indexed (fn [k {:keys [t r]}] {:t t :r (mapv #(+ %1 (noise (+ (* 3 k) %2))) r (range))})
                           (states truth times))
          {:keys [satrec rms]} (fit/fit (template truth) (vec obs) {:sigma-r 1.0})
          sd (/ 0.2 (math/sqrt 12.0))
          worst (apply max (map (fn [t] (v3/distance (:r (sgp4/sgp4 satrec t)) (:r (sgp4/sgp4 truth t)))) times))]
      (is (< (* 0.95 sd) rms (* 1.05 sd)) (str rms))
      (is (< worst 0.01) (str worst)))))
