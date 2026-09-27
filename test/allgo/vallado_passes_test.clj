(ns allgo.vallado-passes-test
  "Pass prediction against brute force: a day of a low satellite's
  elevation over a site, from one of the Spacetrack verification element
  sets, sampled every half second -- every rise, set and culmination the
  finder reports must be where the samples put it."
  (:require [allgo.astro.passes :as passes]
            [allgo.astro.sgp4 :as sgp4]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(def ^:private satrec
  (delay
    (let [[l1 l2] (->> (str/split-lines (slurp (io/resource "data/sgp4/SGP4-VER.TLE")))
                       (remove #(str/starts-with? % "#"))
                       (drop-while #(not (str/starts-with? % "1 06251")))
                       (take 2))]
      (sgp4/twoline2rv l1 l2))))

(def ^:private elevation
  (delay (passes/sgp4-elevation @satrec (math/to-radians 40.0) (math/to-radians -105.0) 1.6)))

(deftest against-brute-force
  (let [el @elevation
        dt (/ 0.5 60.0)
        samples (mapv (fn [k] (let [t (* k dt)] [t (el t)])) (range (int (/ 1440.0 dt))))
        found (passes/passes el 0.0 1440.0 {:step 0.5 :tol 1e-5})]
    (testing "a low satellite passes a few times a day"
      (is (< 2 (count found) 12)))
    (testing "each rise and set is where the sampled elevation crosses the horizon"
      (doseq [{:keys [rise set]} found
              [t sign] [[rise 1.0] [set -1.0]]
              :when t]
        (is (neg? (* sign (el (- t dt)))) (str t))
        (is (pos? (* sign (el (+ t dt)))) (str t))))
    (testing "and nothing is missed: every sample above the horizon lies in a pass"
      (is (every? (fn [[t e]] (or (<= e 0.0)
                                  (some (fn [{:keys [rise set]}] (<= (or rise 0.0) t (or set 1440.0))) found)))
                  samples)))
    (testing "each culmination is the highest sample of its pass"
      (doseq [{:keys [rise set max-elevation]} found
              :when (and rise set)]
        (let [highest (apply max (map second (filter #(<= rise (first %) set) samples)))]
          (is (>= (+ max-elevation 1e-9) highest))
          (is (< (- max-elevation highest) 1e-4)))))))
