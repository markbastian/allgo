(ns procedurals.core
  (:require [procedurals.cave :as c]
            [procedurals.dungeon-generator :as pdg]
            [reagent.core :refer [atom]]
            [reagent.dom.client :as rdom]))

(enable-console-print!)

(println "Edits to this text should show up in your developer console.")

(defn gen-grid [{:keys [w h _i p]}]
  (vec (take 64 (c/ca-cave-iterator (c/ca-grid w h (* p 0.01))))))

(defn update-grid [m]
  (assoc m :caves (gen-grid m)))

(defonce state (atom (update-grid {:w 32 :h 32 :i 18 :p 45})))

(add-watch state :grid-watch (fn [_ _ o n]
                               (when (not= o n)
                                 (update-grid n))))

(defonce root
  (when-let [app-context (. js/document (getElementById "app"))]
    (rdom/create-root app-context)))

(defn render! []
  (when root
    (rdom/render root [pdg/render state])))

(render!)

(defn on-js-reload []
  (render!))
