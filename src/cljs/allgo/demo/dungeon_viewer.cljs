(ns allgo.demo.dungeon-viewer
  (:require [allgo.procedural.cave :as c]
            [cljs.reader :refer [read-string]]))

(defn render [state]
  (let [{:keys [w h i _p caves]} @state
        cw 10
        grid (reduce c/mark-floor (c/init w h) (nth caves i))]
    [:div
     [:svg {:width (* w cw) :height (* h cw)}
      (doall (for [ix (range w) j (range h)]
               [:rect {:key (str ix ":" j) :x (* ix cw) :y (* j cw)
                       :width cw :height cw
                       :stroke :red :fill (case (get-in grid [ix j])
                                            :wall :blue
                                            :floor :green)}]))]
     [:div
      [:input {:type :range :min 0 :max 63 :onChange (fn [e] (swap! state assoc :i (-> e .-target .-value read-string)))} "i"]
      [:input {:type :range :min 0 :max 100 :onChange (fn [e] (swap! state assoc :p (-> e .-target .-value read-string)))} "p"]
      [:input {:type :range :min 1 :max 64 :value w :onChange #(swap! state assoc :w (-> % .-target .-value read-string))} "w"]
      [:input {:type :range :min 1 :max 64 :value h :onChange #(swap! state assoc :h (-> % .-target .-value read-string))} "h"]]
     ;[:button #(prn (with-out-str (pprint state))) "dump state"]
     ]))

(defn init-state [])
