(ns allgo.desktop.quil
  (:require [allgo.procedural.dungeons :as dungeons]
            [clojure.math :as math]
            [quil.core :as q]
            [quil.middleware :as m]))

(defn draw [{:keys [rooms width height dungeon-radius] :as _state}]
  (q/fill (q/color 255))
  (q/rect 0 0 width height)
  (q/translate (* 0.5 width) (* 0.5 height))
  (q/scale (/ (double width) (double dungeon-radius) math/PI))
  (q/stroke 0)
  (q/stroke-weight 0.1)

  (doseq [{:keys [center width height color] :as _room} rooms
          :let [[cx cy] center
                x (- cx (* 0.5 width))
                y (- cy (* 0.5 height))]]
    (q/fill color)
    (q/rect x y width height)))

(defn initial-state [{:keys [num-rooms
                             max-room-width
                             max-room-height
                             dungeon-radius] :as config}]
  (assoc config
         :rooms (->> (repeatedly
                      num-rooms
                      #(dungeons/random-room
                        dungeon-radius
                        max-room-width
                        max-room-height))
                     (mapv (fn [room] (assoc room :color
                                             (q/color
                                              (q/random 255)
                                              (q/random 255)
                                              (q/random 255))))))))

(defn setup [config]
  (q/smooth)
  (q/frame-rate 30)
  (initial-state config))

(defn sim [{:keys [rooms] :as state}]
  (loop [[room & r] rooms f {}]
    (if room
      (recur r f)
      f))
  state)

(defn launch-sketch [{:keys [width height] :as config}]
  (q/sketch
   :title "Flocking Behaviors"
   :setup #(setup config)
   :draw #'draw
   :update #'sim
    ;:mouse-clicked #'io/mouse-click
    ;:mouse-moved #'io/mouse-move
    ;:key-pressed #'io/key-pressed
    ;:key-typed #'io/key-pressed
   :middleware [m/fun-mode]
   :size [width height]))

(comment
  (launch-sketch
   {:width           500
    :height          500
    :num-rooms       40
    :dungeon-radius  20
    :max-room-width  5
    :max-room-height 5}))
