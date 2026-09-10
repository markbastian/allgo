(ns allgo.demo.core
  (:require [allgo.demo.delaunay-viewer :as pdv]
            [allgo.demo.dungeon-generator :as pdg]
            [allgo.procedural.cave :as c]
            [reagent.core :refer [atom]]
            [reagent.dom.client :as rdom]))

(enable-console-print!)

(println "Edits to this text should show up in your developer console.")

(defn gen-grid [{:keys [w h _i p]}]
  (vec (take 64 (c/ca-cave-iterator (c/ca-grid w h (* p 0.01))))))

(defn update-grid [m]
  (assoc m :caves (gen-grid m)))

(defonce state (atom (update-grid {:w 32 :h 32 :i 18 :p 45})))
(defonce delaunay-state (atom (pdv/init-state)))

(add-watch state :grid-watch (fn [_ _ o n]
                               (when (not= o n)
                                 (update-grid n))))

(defn- mount-root [id]
  (when-let [el (. js/document (getElementById id))]
    (rdom/create-root el)))

(defonce root (mount-root "app"))
(defonce delaunay-root (mount-root "delaunay"))

(defn render! []
  (when root
    (rdom/render root [pdg/render state]))
  (when delaunay-root
    (rdom/render delaunay-root [pdv/render delaunay-state])))

(render!)

(defn on-js-reload []
  (render!))
