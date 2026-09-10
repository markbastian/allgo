(ns allgo.desktop.terrain3d
  (:require [allgo.desktop.terrain-shape :as terrain-shape]
            [quil.core :as q]
            [quil.middleware :as m]))

(defn setup [config]
  (q/frame-rate 30)
  (merge (terrain-shape/build-terrain config)
         {:config config :rot-x 0.55 :rot-y 0.8 :zoom 1.0}))

(defn regenerate [state]
  (merge state (terrain-shape/build-terrain (:config state))))

(defn key-pressed [state {:keys [key]}]
  (if (= key :r) (regenerate state) state))

(def max-elevation (- (/ Math/PI 2) 0.1))

(defn mouse-dragged [state {:keys [x y p-x p-y]}]
  (-> state
      (update :rot-y - (* 0.01 (- x p-x)))
      (update :rot-x (fn [rx] (-> rx (+ (* 0.01 (- y p-y)))
                                  (max (- max-elevation))
                                  (min max-elevation))))))

(defn mouse-wheel [state rotation]
  (update state :zoom #(-> % (- (* rotation 0.05)) (max 0.3) (min 4.0))))

(defn draw [{:keys [dim terrain-shape rot-x rot-y zoom]
             {:keys [cell-scale]} :config}]
  (q/background 15 15 25)
  (let [span   (* cell-scale (dec dim))
        radius (/ (* span 1.4) zoom)
        ex     (* radius (Math/cos rot-x) (Math/sin rot-y))
        ey     (* (- radius) (Math/sin rot-x))
        ez     (* radius (Math/cos rot-x) (Math/cos rot-y))]
    (q/camera ex ey ez 0 0 0 0 1 0)
    (q/perspective (/ Math/PI 3) (/ (q/width) (q/height)) 1 (* radius 10))
    (q/directional-light 220 220 200 -0.4 -1 -0.3)
    (q/ambient-light 70 70 80)
    (q/point-light 255 255 240 ex (- radius) ez)
    (q/shape terrain-shape)))

(defn launch-sketch [config]
  (q/sketch
   :title "Fractal Terrain"
   :setup #(setup config)
   :draw #'draw
   :key-pressed #'key-pressed
   :mouse-dragged #'mouse-dragged
   :mouse-wheel #'mouse-wheel
   :size [(:width-px config 900) (:height-px config 700)]
   :renderer :p3d
   :middleware [m/fun-mode]))

(comment
  (launch-sketch
   {:iterations   7
    :width        1.0
    :height-scale 200
    :cell-scale   8
    :width-px     900
    :height-px    700})

  (launch-sketch
   {:generator    :tin
    :tin-points   500
    :dim          129
    :height-scale 200
    :cell-scale   8
    :width-px     900
    :height-px    700}))
