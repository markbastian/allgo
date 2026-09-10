(ns procedurals.delaunay-viewer
  (:require [procedurals.delaunay :as delaunay]
            [clojure.string :as str]))

(def width 320)
(def height 320)

(defn rand-points [n]
  (vec (repeatedly n #(vector (rand-int width) (rand-int height)))))

(defn init-state [] {:n 20 :points (rand-points 20)})

(defn- pts->str [pts]
  (str/join " " (map (fn [[x y]] (str x "," y)) pts)))

(defn render [state]
  (let [{:keys [points n]} @state
        triangles (delaunay/triangulate points)
        cells     (delaunay/voronoi-cells triangles)]
    [:div
     [:svg {:width width :height height :style {:background "#000"}}
      (doall (for [cell (vals cells) :when (>= (count cell) 3)]
               [:polygon {:key (hash cell) :points (pts->str cell)
                          :fill "none" :stroke "#3fa34d" :stroke-width 1}]))
      (doall (for [{:keys [points]} triangles
                   [[x1 y1] [x2 y2]] (partition 2 1 (conj points (first points)))]
               [:line {:key (hash [x1 y1 x2 y2]) :x1 x1 :y1 y1 :x2 x2 :y2 y2 :stroke "#444"}]))
      (doall (for [[x y] points]
               [:circle {:key (hash [x y]) :cx x :cy y :r 3 :fill "#f55"}]))]
     [:div
      [:input {:type :range :min 3 :max 80 :value n
               :on-change (fn [e]
                            (let [n' (-> e .-target .-value js/parseInt)]
                              (swap! state assoc :n n' :points (rand-points n'))))}]
      [:span (str " points: " n)]
      [:button {:on-click #(swap! state assoc :points (rand-points n))} "Regenerate"]]]))
