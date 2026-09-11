(ns allgo.demo.dungeon
  "Canvas demo for `allgo.procedural.dungeons`.

  Draws the finished tile grid, so what you see is the generator's actual
  output rather than a tidied redrawing of it: hub rooms bright, the small
  rooms a corridor happened to cross dimmer behind them, and the corridor
  runs themselves dimmer still. The graph overlay shows which pairs of hubs
  the spanning tree chose to connect, and the ghost overlay shows the rooms
  separation produced that never made the final map -- usually most of
  them, which is the part of the algorithm that is hardest to picture."
  (:require [allgo.procedural.dungeons :as dungeons]
            ["lil-gui" :default GUI]))

(def ^:private background "#05070d")

(def ^:private palette
  {:corridor "#1c2740"
   :hallway  "#33598c"
   :main     "#e0a85c"
   :edge     "rgba(126, 212, 160, 0.55)"
   :ghost    "rgba(120, 140, 180, 0.16)"
   :text     "#7d8aa3"})

(def ^:private controls
  #js {:rooms         120
       :spread        150
       :sizeMean      24
       :sizeSpread    8
       :hubThreshold  1.25
       :extraEdges    0.10
       :corridorWidth 3
       :shape         "circle"
       :showGraph     true
       :showDiscarded false})

(defn- config []
  (let [ellipse? (= "ellipse" (.-shape controls))
        spread   (.-spread controls)]
    {:room-count       (.-rooms controls)
     :radius           spread
     :ellipse          (when ellipse? [(* 6 spread) (* 0.4 spread)])
     :width-mean       (.-sizeMean controls)
     :width-sd         (.-sizeSpread controls)
     :height-mean      (.-sizeMean controls)
     :height-sd        (.-sizeSpread controls)
     :main-threshold   (.-hubThreshold controls)
     :extra-edge-ratio (.-extraEdges controls)
     :corridor-width   (.-corridorWidth controls)}))

;; ---------------------------------------------------------------------------
;; Fitting the dungeon to the viewport

(defn- grid-extent
  "`[min-tx min-ty max-tx max-ty]` over every tile on the map."
  [grid]
  (reduce (fn [[x0 y0 x1 y1] [tx ty]]
            [(min x0 tx) (min y0 ty) (max x1 tx) (max y1 ty)])
          [##Inf ##Inf ##-Inf ##-Inf]
          (keys grid)))

(defn- fit
  "A tile -> pixel transform that centres the whole map in `[w h]` with a
  little air around it. Returns `[scale offset-x offset-y]`."
  [grid [w h]]
  (let [[x0 y0 x1 y1] (grid-extent grid)
        cols  (inc (- x1 x0))
        rows  (inc (- y1 y0))
        scale (* 0.94 (min (/ w cols) (/ h rows)))]
    [scale
     (+ (/ (- w (* cols scale)) 2) (* (- x0) scale))
     (+ (/ (- h (* rows scale)) 2) (* (- y0) scale))]))

;; ---------------------------------------------------------------------------

(defn- draw-tiles! [^js ctx grid [scale ox oy]]
  ;; Tiles are drawn a hair oversized; at fractional scales exact widths
  ;; leave hairline seams between neighbours that read as cracks in the
  ;; floor.
  (let [side (js/Math.ceil (+ scale 0.5))]
    (doseq [[kind tiles] (group-by val grid)]
      (set! (.-fillStyle ctx) (palette kind))
      (doseq [[[tx ty] _] tiles]
        (.fillRect ctx (+ ox (* tx scale)) (+ oy (* ty scale)) side side)))))

(defn- draw-ghosts!
  "The rooms separation produced that did not make the final map."
  [^js ctx rooms tile-size [scale ox oy]]
  (set! (.-strokeStyle ctx) (palette :ghost))
  (set! (.-lineWidth ctx) 1)
  (doseq [room rooms
          :when (= :minor (:kind room))
          :let  [[x0 y0 x1 y1] (dungeons/bounds room)
                 px (fn [v] (+ ox (* scale (/ v tile-size))))
                 py (fn [v] (+ oy (* scale (/ v tile-size))))]]
    (.strokeRect ctx (px x0) (py y0) (- (px x1) (px x0)) (- (py y1) (py y0)))))

(defn- draw-graph!
  "The edges the spanning tree kept, drawn hub centre to hub centre."
  [^js ctx rooms edges tile-size [scale ox oy]]
  (let [by-id (into {} (map (juxt :id identity)) rooms)
        px    (fn [v] (+ ox (* scale (/ v tile-size))))
        py    (fn [v] (+ oy (* scale (/ v tile-size))))]
    (set! (.-strokeStyle ctx) (palette :edge))
    (set! (.-lineWidth ctx) 1.5)
    (doseq [e edges
            :let [[u v] (vec e)
                  [ux uy] (:center (by-id u))
                  [vx vy] (:center (by-id v))]]
      (.beginPath ctx)
      (.moveTo ctx (px ux) (py uy))
      (.lineTo ctx (px vx) (py vy))
      (.stroke ctx))
    (set! (.-fillStyle ctx) (palette :edge))
    (doseq [{:keys [kind center]} rooms
            :when (= :main kind)
            :let  [[cx cy] center]]
      (.beginPath ctx)
      (.arc ctx (px cx) (py cy) 2.5 0 (* 2 js/Math.PI))
      (.fill ctx))))

(defn- draw-caption! [^js ctx {:keys [rooms edges grid]} ms [_ h]]
  (let [kinds (frequencies (map :kind rooms))]
    (set! (.-fillStyle ctx) (palette :text))
    (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
    ;; Kept short deliberately: the tile is 480px and this is monospace, so
    ;; a longer line is silently clipped rather than wrapped.
    (.fillText ctx
               (str (:main kinds 0) " hubs · "
                    (:hallway kinds 0) " halls · "
                    (:minor kinds 0) " unused · "
                    (count edges) " runs · "
                    (count grid) " tiles · "
                    ms "ms")
               10 (- h 10))))

(defn- draw! [^js ctx dungeon ms [w h :as size]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (when (seq (:grid dungeon))
    (let [tile-size (get-in dungeon [:config :tile-size])
          xf        (fit (:grid dungeon) size)]
      (when (.-showDiscarded controls)
        (draw-ghosts! ctx (:rooms dungeon) tile-size xf))
      (draw-tiles! ctx (:grid dungeon) xf)
      (when (.-showGraph controls)
        (draw-graph! ctx (:rooms dungeon) (:edges dungeon) tile-size xf))
      (draw-caption! ctx dungeon ms size))))

;; ---------------------------------------------------------------------------

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        state  (atom {:bounds [0 0] :dungeon nil :ms 0})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(render! []
              (let [{:keys [dungeon ms bounds]} @state]
                (when dungeon (draw! ctx dungeon ms bounds))))
            (resize! []
              (let [dpr (or js/window.devicePixelRatio 1)
                    w   (.-clientWidth container)
                    h   (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) (* w dpr))
                  (set! (.-height canvas) (* h dpr))
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (.setTransform ctx dpr 0 0 dpr 0 0)
                  (swap! state assoc :bounds [w h])
                  (render!))))
            (generate! []
              ;; A full generate runs a few hundred separation sweeps, so it
              ;; is done on demand rather than per frame -- this demo has
              ;; nothing to animate.
              (let [t0 (js/performance.now)
                    d  (dungeons/generate (config))]
                (swap! state assoc
                       :dungeon d
                       :ms (js/Math.round (- (js/performance.now) t0)))
                (render!)))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "rooms" 20 250 10) (.onFinishChange generate!))
          (-> (.add controls "spread" 40 300 10) (.onFinishChange generate!))
          (-> (.add controls "shape" #js ["circle" "ellipse"]) (.onFinishChange generate!))
          (-> (.add controls "sizeMean" 8 48 4) (.onFinishChange generate!))
          (-> (.add controls "sizeSpread" 0 20 1) (.onFinishChange generate!))
          (-> (.add controls "hubThreshold" 1.0 2.0 0.05) (.onFinishChange generate!))
          (-> (.add controls "extraEdges" 0 0.5 0.01) (.onFinishChange generate!))
          (-> (.add controls "corridorWidth" 1 7 2) (.onFinishChange generate!))
          (-> (.add controls "showGraph") (.onChange render!))
          (-> (.add controls "showDiscarded") (.onChange render!))
          (.add #js {:regenerate generate!} "regenerate")))
      {:start (fn [] (resize!) (render!))
       :stop  (fn [])})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "dungeon")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
