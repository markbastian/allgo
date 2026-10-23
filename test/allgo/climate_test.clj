(ns allgo.climate-test
  (:require [allgo.procedural.climate :as cl]
            [allgo.procedural.tectonics :as tc]
            [clojure.test :refer [deftest is testing]]))

(def ^:private world (tc/plates {:resolution 4 :seed 2}))
(def ^:private clim (delay (cl/climate world)))

(defn- mean [xs] (/ (reduce + xs) (double (max 1 (count xs)))))

(def ^:private land
  (filter #(>= (double ((:elevation world) %)) 0.0) (range (count (:points world)))))

(defn- deg [i] (Math/toDegrees (Math/abs (cl/latitude ((:points world) i)))))

(deftest temperature-test
  (testing "warm at the equator, cold at the poles and up mountains"
    (is (= 28.0 (cl/temperature 0.0 0.0)))
    (is (= -30.0 (cl/temperature (/ Math/PI 2) 0.0)))
    (is (< (abs (- 13.5 (cl/temperature (/ Math/PI 4) 0.0))) 0.1))
    (is (= 15.0 (cl/temperature 0.0 2.0)))))

(deftest wind-test
  (let [east (fn [_] [0.0 0.0 1.0])
        at (fn [lat] [(Math/cos lat) (Math/sin lat) 0.0])
        component (fn [v u] (reduce + (map * v u)))]
    (testing "trade winds blow from the east, westerlies from the west"
      ;; at longitude 0 east is +z
      (is (neg? (component (cl/wind (at (Math/toRadians 15.0))) (east 15))))
      (is (pos? (component (cl/wind (at (Math/toRadians 45.0))) (east 45))))
      (is (neg? (component (cl/wind (at (Math/toRadians 80.0))) (east 80)))))
    (testing "the trades drift toward the equator in both hemispheres"
      (let [north (fn [lat] [(- (Math/sin lat)) (Math/cos lat) 0.0])]
        (is (neg? (component (cl/wind (at (Math/toRadians 15.0))) (north (Math/toRadians 15.0)))))
        (is (pos? (component (cl/wind (at (Math/toRadians -15.0))) (north (Math/toRadians -15.0)))))))))

(deftest rain-test
  (let [{:keys [rain wind]} @clim
        {:keys [points neighbors elevation]} world
        ;; how far the ground climbs into each cell from upwind
        climb (fn [i]
                (let [p (points i)
                      w (wind i)
                      up (filter (fn [j] (pos? (reduce + (map * (map - p (points j)) w))))
                                 (nth neighbors i))]
                  (when (seq up)
                    (- (elevation i) (mean (map #(max 0.0 (elevation %)) up))))))]
    (testing "windward slopes are wetter than lee slopes"
      (let [climbs (keep (fn [i] (when-let [c (climb i)] [c (rain i)])) land)]
        (is (> (mean (map second (filter #(> (first %) 0.05) climbs)))
               (* 2.0 (mean (map second (filter #(< (first %) -0.05) climbs))))))))
    (testing "the equator is wetter than the subtropics"
      (is (> (mean (map rain (filter #(< (deg %) 10.0) land)))
             (mean (map rain (filter #(< 20.0 (deg %) 35.0) land))))))))

(deftest rivers-test
  (let [{:keys [downstream flow river? lake?]} @clim
        n (count downstream)
        land? (set land)]
    (testing "every drop of rain finds the sea"
      (doseq [i land]
        (is (loop [c i steps 0]
              (cond (or (neg? (long c)) (not (land? c))) true
                    (> steps n) false
                    :else (recur (downstream c) (inc steps)))))))
    (testing "rivers only grow downstream"
      (doseq [i land :let [d (downstream i)] :when (and (>= d 0) (land? d))]
        (is (>= (flow d) (flow i)))))
    (testing "rivers are about the share asked for, and lakes are not rivers"
      (let [share (/ (count (filter river? land)) (double (count land)))]
        (is (< 0.02 share 0.05) (str share)))
      (is (not-any? #(and (river? %) (lake? %)) land)))))

(deftest biome-test
  (testing "Whittaker's corners"
    (is (= :tropical-rainforest (cl/biome 25.0 0.9)))
    (is (= :desert (cl/biome 25.0 0.02)))
    (is (= :taiga (cl/biome 3.0 0.5)))
    (is (= :ice (cl/biome -20.0 0.5))))
  (testing "every land cell has a biome and no sea cell does"
    (let [b (:biome @clim)]
      (is (every? some? (map b land)))
      (is (every? nil? (map b (remove (set land) (range (count b)))))))))
