(ns allgo.demo.cave
  "Canvas demo for `allgo.procedural.cave`.

  A cellular automaton carving caverns: seed a fraction of the cells at
  random, then repeatedly keep every cell with more than four live
  neighbors. Scrubbing `iterations` walks the same seed forward one
  generation at a time, which is the point worth seeing -- noise resolves
  into rounded chambers within a handful of steps and then barely moves.

  What it does not do is guarantee those chambers reach each other, so
  `connect` joins the disconnected islands with corridors. The caption
  counts them, and turning the toggle off shows how many the automaton
  leaves stranded."
  (:require [allgo.procedural.cave :as cave]
            ["lil-gui" :default GUI]))

(def ^:private background "#05070d")
(def ^:private wall-color "#10172a")
(def ^:private floor-color "#6f9fe0")
(def ^:private text-color "#7d8aa3")

(def ^:private generations 25)

(def ^:private ^js controls
  #js {:width      64
       :height     48
       :fill       45
       :iterations 6
       :connect    true})

(defn- draw! [^js ctx grid islands [w h]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (let [rows  (count grid)
        cols  (count (first grid))
        ;; Square cells, centered: the grid's aspect ratio rarely matches
        ;; the viewport's.
        scale (min (/ w cols) (/ (- h 24) rows))
        ox    (/ (- w (* cols scale)) 2)
        oy    (/ (- (- h 24) (* rows scale)) 2)
        side  (js/Math.ceil (+ scale 0.5))]
    (doseq [r (range rows)
            c (range cols)]
      (set! (.-fillStyle ctx)
            (if (= :floor (get-in grid [r c])) floor-color wall-color))
      (.fillRect ctx (+ ox (* c scale)) (+ oy (* r scale)) side side))
    (set! (.-fillStyle ctx) text-color)
    (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
    (.fillText ctx
               (str (count (cave/floor-coords grid)) " floor cells  ·  "
                    islands (if (= 1 islands) " cavern" " caverns"))
               10 (- h 10))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        state  (atom {:bounds [0 0] :generations []})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(current []
              (let [gens (:generations @state)]
                (when (seq gens)
                  (let [grid (nth gens (min (.-iterations controls) (dec (count gens))))]
                    (if (.-connect controls) (cave/connect grid) grid)))))
            (render! []
              (when-let [grid (current)]
                (let [islands (count (cave/find-islands (cave/meadow-coords grid)))]
                  (draw! ctx grid islands (:bounds @state)))))
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
              ;; Every generation is computed once and kept, so dragging the
              ;; iterations slider scrubs an existing sequence instead of
              ;; re-running the automaton from a fresh seed each frame --
              ;; which would reseed the cave on every drag and show nothing.
              (swap! state assoc :generations
                     (vec (take generations
                                (cave/ca-cave-seq (.-width controls)
                                                  (.-height controls)
                                                  (/ (.-fill controls) 100.0)))))
              (render!))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "width" 16 128 8) (.onFinishChange generate!))
          (-> (.add controls "height" 16 96 8) (.onFinishChange generate!))
          (-> (.add controls "fill" 30 60 1) (.onFinishChange generate!))
          (-> (.add controls "iterations" 0 (dec generations) 1) (.onChange render!))
          (-> (.add controls "connect") (.onChange render!))
          (.add #js {:regenerate generate!} "regenerate")))
      {:start (fn [] (resize!) (render!))
       :stop  (fn [])})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "caves")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
