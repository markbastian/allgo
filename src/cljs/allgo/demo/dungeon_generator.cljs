(ns allgo.demo.dungeon-generator
  (:require [allgo.procedural.dungeons :as pd]))

(defn rand-room []
  (let [w (rand-int 20) h (rand-int 20)
        [cx cy] (mapv + [160 160] (pd/random-int-point 160 4))]
    {:width (* w 4) :height (* h 4) :cx cx :cy cy}))

(def width 480)
(def height 480)

(defn rand-seq [n] (take n (repeatedly (comp (partial * width) rand))))
(defn rand-pts [n] (mapv vector (rand-seq n) (rand-seq n)))
(def pts (rand-pts 10))

(def rest-length 20)
(def k -2.0)

(defn point-force [p q]
  (let [v (map - p q)
        m (Math/sqrt (reduce + (map * v v)))
        l (- m rest-length)]
    (mapv #(/ (* k l %) m) v)))

(defn forces [masses]
  (loop [[i & r] (range (count masses))
         forces (vec (repeat (count masses) [0.0 0.0]))]
    (if i
      (recur r (reduce
                (fn [f j]
                  (let [x (point-force (masses i) (masses j))]
                    (-> f
                        (update i #(mapv - % x))
                        (update j #(mapv + % x)))))
                forces r))
      forces)))

#_(defn links [pts]
    (loop [[p & r] pts res []]
      (if p
        (recur r (into res (map #(mapv - p %) r)))
        res)))

(defn links [pts]
  (loop [[p & r] pts res []]
    (if p
      (recur r (into res (map (fn [l] [p l]) r)))
      res)))

#_:clj-kondo/ignore
(defn render [state]
  (let [{:keys [w h i p caves]} @state
        cw 10]
    [:div
     [:div.demo-viewport
      [:svg {:width width :height height}
       [:rect {:x 0 :y 0 :width width :height height :fill :black}]
       #_(doall (for [i (range width) j (range height)]
                  [:rect {:key i :x i :y j :width 1 :height 1 :stroke :red :fill :black}]))
       #_(doall (for [i (range 10)]
                  (let [{:keys [width height cx cy]} (rand-room)]
                    [:rect {:key i
                            :x (- cx (/ width 2))
                            :y (- cy (/ height 2))
                            :width width :height height
                            :stroke :red :fill :blue}])))
       (doall (for [[[x1 y1] [x2 y2] :as l] (links pts)]
                [:line {:x1 x1 :y1 y1 :x2 x2 :y2 y2 :stroke :green :key l}]))
       (doall (for [i (range (count pts))]
                (let [pw 2
                      ph 2
                      [cx cy] (pts i)]
                  [:rect {:key i
                          :x (- cx (/ pw 2))
                          :y (- cy (/ ph 2))
                          :width pw :height ph
                          :stroke :red :fill :blue}])))]]]))
