(ns allgo.procedural.dungeons
  (:require [clojure.math :as math]))

; http://www.gamasutra.com/blogs/AAdonaac/20150903/252889/Procedural_Dungeon_Generation_Algorithm.php
; https://web.archive.org/web/20210823235922/https://gamasutra.com/blogs/AAdonaac/20150903/252889/Procedural_Dungeon_Generation_Algorithm.php

(defn random-point [radius]
  (let [t (* 2.0 math/PI (rand))
        u (+ (rand) (rand))
        r (if (> u 1) (- 2 u) u)]
    [(* radius r (math/cos t)) (* radius r (math/sin t))]))

(defn random-room [radius max-width max-height]
  {:center (random-point radius)
   :width (* (rand) max-width)
   :height (* (rand) max-height)})

(defn roundm [n m]
  (int (* (math/floor (/ (dec (+ n m)) m)) m)))

(defn random-int-point [radius tile-size]
  (let [t (* 2.0 math/PI (rand))
        u (+ (rand) (rand))
        r (if (> u 1) (- 2 u) u)]
    [(roundm (* radius r (math/cos t)) tile-size)
     (roundm (* radius r (math/sin t)) tile-size)]))

(defn center [s]
  (->> s
       (map :center)
       (apply map +)
       (mapv #(/ % (count s)))))

(comment
  (->> (repeatedly 20 #(random-room 20 4 4))
       center))
