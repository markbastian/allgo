(ns procedurals.boids-viewer
  "Canvas demo for `procedurals.boids`. Each boid is drawn as a dart pointing
  along its heading and tinted by it, so alignment emerging out of a random
  start reads as the flock converging on a single colour."
  (:require [procedurals.boids :as boids]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private background "#05070d")

(def ^:private controls
  #js {:boids          30
       :separation     1.6
       :alignment      1.0
       :cohesion       0.9
       :perception     60
       :personalSpace  24
       :speed          2.4
       :trails         true})

(defn- params []
  {:separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (.-personalSpace controls)
   :max-speed         (.-speed controls)})

(defn- resize-to-flock
  "Grow or trim the flock to match the boid-count control."
  [flock n bounds max-speed]
  (let [have (count flock)]
    (cond
      (= have n) flock
      (< n have) (subvec flock 0 n)
      :else      (into flock (repeatedly (- n have) #(boids/random-boid bounds max-speed))))))

(defn- draw-boid! [^js ctx {[x y] :pos [vx vy] :vel}]
  (let [heading (math/atan2 vy vx)
        hue     (mod (+ 180 (math/to-degrees heading)) 360)]
    (set! (.-fillStyle ctx) (str "hsl(" hue " 85% 62%)"))
    (.save ctx)
    (.translate ctx x y)
    (.rotate ctx heading)
    (.beginPath ctx)
    (.moveTo ctx 6 0)
    (.lineTo ctx -4 3)
    (.lineTo ctx -2 0)
    (.lineTo ctx -4 -3)
    (.closePath ctx)
    (.fill ctx)
    (.restore ctx)))

(defn- draw! [^js ctx flock [w h]]
  (set! (.-fillStyle ctx) (if (.-trails controls) "rgba(5, 7, 13, 0.22)" background))
  (.fillRect ctx 0 0 w h)
  (doseq [boid flock] (draw-boid! ctx boid)))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        running? (atom false)
        state  (atom {:bounds [0 0] :flock []})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(resize! []
              (let [dpr (or js/window.devicePixelRatio 1)
                    w   (.-clientWidth container)
                    h   (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) (* w dpr))
                  (set! (.-height canvas) (* h dpr))
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (.setTransform ctx dpr 0 0 dpr 0 0)
                  (set! (.-fillStyle ctx) background)
                  (.fillRect ctx 0 0 w h)
                  (swap! state assoc :bounds [w h]))))
            (reset-flock! []
              (swap! state assoc :flock (boids/flock (.-boids controls) (:bounds @state) (params))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [{:keys [bounds]} @state
                      p     (params)
                      flock (-> (:flock @state)
                                (resize-to-flock (.-boids controls) bounds (:max-speed p))
                                (boids/step bounds p))]
                  (swap! state assoc :flock flock)
                  (draw! ctx flock bounds))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (reset-flock!)
      (let [gui (GUI. #js {:container container})]
        (.add gui controls "boids" 10 400 10)
        (.add gui controls "separation" 0 3 0.1)
        (.add gui controls "alignment" 0 3 0.1)
        (.add gui controls "cohesion" 0 3 0.1)
        (.add gui controls "perception" 20 140 5)
        (.add gui controls "personalSpace" 5 60 1)
        (.add gui controls "speed" 0.5 6 0.1)
        (.add gui controls "trails")
        (.add gui #js {:reset reset-flock!} "reset"))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "boids")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
