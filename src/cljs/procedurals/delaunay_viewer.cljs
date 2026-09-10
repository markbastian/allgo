(ns procedurals.delaunay-viewer
  (:require [procedurals.delaunay :as delaunay]
            [clojure.string :as str]))

(def width 480)
(def height 480)

(defn rand-points [n]
  (vec (repeatedly n #(vector (rand-int width) (rand-int height)))))

(defn init-state [] {:n 20 :points (rand-points 20)})

(defn- pts->str [pts]
  (str/join " " (map (fn [[x y]] (str x "," y)) pts)))

(defn render [state]
  (let [{:keys [points n]} @state
        triangles (delaunay/triangulate points)
        cells     (filter #(>= (count %) 3) (vals (delaunay/voronoi-cells triangles)))
        edges     (mapcat (fn [{:keys [points]}] (partition 2 1 (conj points (first points)))) triangles)]
    [:div
     [:div.demo-viewport
      [:svg {:width width :height height}
       (doall (map-indexed (fn [i cell]
                             [:polygon {:key (str "cell-" i) :points (pts->str cell)
                                        :fill "none" :stroke "#3fa34d" :stroke-width 1}])
                           cells))
       (doall (map-indexed (fn [i [[x1 y1] [x2 y2]]]
                             [:line {:key (str "edge-" i) :x1 x1 :y1 y1 :x2 x2 :y2 y2 :stroke "#444"}])
                           edges))
       (doall (map-indexed (fn [i [x y]]
                             [:circle {:key (str "pt-" i) :cx x :cy y :r 3 :fill "#f55"}])
                           points))]]
     [:div.demo-controls
      [:input {:type :range :min 3 :max 80 :value n
               :on-change (fn [e]
                            (let [n' (-> e .-target .-value js/parseInt)]
                              (swap! state assoc :n n' :points (rand-points n'))))}]
      [:span (str "points: " n)]
      [:button {:on-click #(swap! state assoc :points (rand-points n))} "Regenerate"]]]))
