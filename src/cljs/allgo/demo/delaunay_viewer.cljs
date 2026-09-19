(ns allgo.demo.delaunay-viewer
  "Reagent demo for `allgo.geometry.delaunay`: a Delaunay triangulation drawn
  over its dual Voronoi diagram.

  The one tile rendered with Reagent rather than a canvas -- it is static
  between clicks, so there is no animation loop to drive and SVG says what
  it means. It still starts and stops through the same lifecycle as the
  rest, so the shell has a single way to swap demos in and out."
  (:require [allgo.geometry.delaunay :as delaunay]
            [clojure.string :as str]
            [reagent.core :as r]
            [reagent.dom.client :as rdom]))

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
      ;; A viewBox rather than a fixed pixel size, so the drawing fills
      ;; whatever the tile has become. Points stay in the 480-unit space
      ;; they are generated in.
      [:svg {:viewBox (str "0 0 " width " " height)
             :width "100%" :height "100%"
             :preserveAspectRatio "xMidYMid meet"}
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

;; ---------------------------------------------------------------------------

(defonce ^:private state (r/atom nil))
(defonce ^:private root (atom nil))

(defn start! []
  (when-let [el (.getElementById js/document "delaunay")]
    (when-not @state (reset! state (init-state)))
    (rdom/render (or @root (reset! root (rdom/create-root el)))
                 [render state])))

(defn stop! []
  ;; Nothing is running between clicks, so there is nothing to tear down;
  ;; the mounted tree is left in place and simply hidden by the shell.
  nil)
