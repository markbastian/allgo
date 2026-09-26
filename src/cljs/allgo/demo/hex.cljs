(ns allgo.demo.hex
  "Canvas demo for `allgo.geometry.hex`, with `allgo.search` doing the
  pathfinding.

  The grid is drawn from `hex/corners`, and the hex under the cursor comes
  straight back from `hex/pixel->hex` -- the same inverse-matrix-and-round
  that makes hex picking a two-line problem rather than a hit-test against
  six edges. Click to move the target, paint walls, and switch overlays to
  see each algorithm drawn on the same grid.

  The path overlay shades every hex the search looked at behind the route
  it settled on, which is the difference between the algorithms made
  visible: A* reaches for the goal, Dijkstra spreads evenly in all
  directions, greedy runs at the goal and can be fooled by a wall."
  (:require [allgo.geometry.hex :as hex]
            [allgo.search :as search]
            ["lil-gui" :default GUI]
            [clojure.string :as str]))

(def ^:private background "#05070d")

(def ^:private palette
  {:grid      "#1a2236"
   :source    "#e0a85c"
   :target    "#7ed4a0"
   :wall      "#39405a"
   :overlay   "rgba(126, 212, 160, 0.30)"
   :explored  "rgba(126, 212, 160, 0.10)"
   :hover     "rgba(224, 168, 92, 0.35)"
   :text      "#7d8aa3"
   :label     "#5d6880"})

(def ^:private ^js controls
  #js {:orientation "pointy"
       :shape       "hexagon"
       :mapSize     6
       :overlay     "path"
       :algorithm   "a-star"
       :radius      3
       :labels      "none"
       :paintWalls  false
       :showSearch  true})

;; ---------------------------------------------------------------------------

(defn- board
  "The hexes of the current map shape."
  []
  (let [n (.-mapSize controls)
        o (if (= "flat" (.-orientation controls)) :flat :pointy)]
    (case (.-shape controls)
      "hexagon"       (hex/hexagon n)
      "triangle"      (hex/triangle n)
      "parallelogram" (hex/parallelogram (- n) n (- n) n)
      "rectangle"     (hex/rectangle o (inc n) (inc n))
      (hex/hexagon n))))

(defn- fit-layout
  "A layout sized so the whole board fits the viewport with a margin."
  [hexes [w h]]
  (let [o     (if (= "flat" (.-orientation controls)) :flat :pointy)
        probe (hex/layout o 1.0)
        pts   (mapcat #(hex/corners probe %) hexes)
        xs    (map first pts)
        ys    (map second pts)
        [x0 x1] [(apply min xs) (apply max xs)]
        [y0 y1] [(apply min ys) (apply max ys)]
        size  (* 0.94 (min (/ w (- x1 x0)) (/ (- h 26) (- y1 y0))))]
    (hex/layout o size
                [(- (/ w 2) (* size (/ (+ x0 x1) 2)))
                 (- (/ (- h 26) 2) (* size (/ (+ y0 y1) 2)))])))

(defn- search-for
  "Runs the selected algorithm, returning `[path explored]`."
  [source target walls on-board?]
  (let [neighbors (fn [h] (filterv #(and (on-board? %) (not (walls %))) (hex/neighbors h)))
        spec       {:start source :goal target
                    :neighbors neighbors :heuristic hex/distance}
        seq-fn     (case (.-algorithm controls)
                     "dijkstra"      search/dijkstra-seq
                     "greedy"        search/greedy-seq
                     "breadth-first" search/breadth-first-seq
                     search/a-star-seq)
        ;; Taking the goal state rather than calling the path function
        ;; directly, so the hexes the search examined can be drawn too.
        final      (search/goal-state (seq-fn spec))]
    [(or (search/path final) []) (set (keys (:came-from final)))]))

(defn- overlay-hexes
  "The hexes the current overlay highlights, and the ones a search looked
  at getting there."
  [{:keys [source target walls board-set]}]
  (let [n         (.-radius controls)
        on-board? board-set]
    (case (.-overlay controls)
      "range"     [(set (filter on-board? (hex/hexes-within source n))) #{}]
      "ring"      [(set (filter on-board? (hex/ring source n))) #{}]
      "spiral"    [(set (filter on-board? (hex/spiral source n))) #{}]
      "line"      [(set (hex/line source target)) #{}]
      "reachable" [(hex/reachable source n #(or (walls %) (not (on-board? %)))) #{}]
      "visible"   [(set (filter on-board? (hex/visible source n walls))) #{}]
      "path"      (let [[p explored] (search-for source target walls on-board?)]
                    [(set p) explored])
      [#{} #{}])))

;; ---------------------------------------------------------------------------

(defn- label-for [h]
  (case (.-labels controls)
    "axial"   (str/join "," (hex/axial h))
    "cube"    (str/join "," (hex/cube h))
    "offset"  (str/join "," (hex/->offset (if (= "flat" (.-orientation controls)) :odd-q :odd-r) h))
    "doubled" (str/join "," (hex/->doublewidth h))
    nil))

(defn- fill-hex! [^js ctx layout h color]
  (let [[[x0 y0] & more] (hex/corners layout h)]
    (set! (.-fillStyle ctx) color)
    (.beginPath ctx)
    (.moveTo ctx x0 y0)
    (doseq [[x y] more] (.lineTo ctx x y))
    (.closePath ctx)
    (.fill ctx)))

(defn- stroke-hex! [^js ctx layout h color width]
  (let [[[x0 y0] & more] (hex/corners layout h)]
    (set! (.-strokeStyle ctx) color)
    (set! (.-lineWidth ctx) width)
    (.beginPath ctx)
    (.moveTo ctx x0 y0)
    (doseq [[x y] more] (.lineTo ctx x y))
    (.closePath ctx)
    (.stroke ctx)))

(defn- draw! [^js ctx {:keys [layout hexes source target walls hover] :as state} [w h]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (let [[highlit explored] (overlay-hexes state)]
    (doseq [hx hexes]
      (fill-hex! ctx layout hx
                 (cond (walls hx)                          (palette :wall)
                       (and (.-showSearch controls)
                            (explored hx) (not (highlit hx))) (palette :explored)
                       :else                                (palette :grid)))
      (when (highlit hx) (fill-hex! ctx layout hx (palette :overlay)))
      (stroke-hex! ctx layout hx "#0b1020" 1))

    (when (and hover (some #{hover} hexes))
      (fill-hex! ctx layout hover (palette :hover)))

    (doseq [[hx color] [[source (palette :source)] [target (palette :target)]]]
      (when (some #{hx} hexes)
        (fill-hex! ctx layout hx color)
        (stroke-hex! ctx layout hx "#0b1020" 1)))

    (when-not (= "none" (.-labels controls))
      (set! (.-fillStyle ctx) (palette :label))
      (set! (.-font ctx) "9px ui-monospace, SFMono-Regular, Menlo, monospace")
      (set! (.-textAlign ctx) "center")
      (doseq [hx hexes
              :let [[x y] (hex/->pixel layout hx)]]
        (.fillText ctx (label-for hx) x (+ y 3)))
      (set! (.-textAlign ctx) "left"))

    (set! (.-fillStyle ctx) (palette :text))
    (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
    (.fillText ctx
               (str (count hexes) " hexes · "
                    (when hover (str "at " (str/join "," hover) " · "))
                    "dist " (hex/distance source target)
                    (when (= "path" (.-overlay controls))
                      (str " · path " (max 0 (dec (count highlit)))
                           " · looked at " (count explored))))
               10 (- h 10))))

;; ---------------------------------------------------------------------------

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        state  (atom {:bounds [0 0] :source [0 0] :target [3 -1]
                      :walls #{} :hover nil})]
    (set! (.-style canvas) "display:block;border-radius:8px;cursor:crosshair")
    (.appendChild container canvas)
    (letfn [(rebuild! []
              (let [hexes (board)
                    {:keys [bounds]} @state]
                (when (pos? (first bounds))
                  (swap! state assoc
                         :hexes hexes
                         :board-set (set hexes)
                         :layout (fit-layout hexes bounds)))))
            (render! []
              (let [{:keys [bounds layout]} @state]
                (when layout (draw! ctx @state bounds))))
            (regenerate! [] (rebuild!) (render!))
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
                  (rebuild!)
                  (render!))))
            (hex-at [^js e]
              (let [rect (.getBoundingClientRect canvas)]
                (hex/pixel->hex (:layout @state)
                                [(- (.-clientX e) (.-left rect))
                                 (- (.-clientY e) (.-top rect))])))]
      (.addEventListener canvas "mousemove"
                         (fn [e]
                           (let [h (hex-at e)]
                             (when (not= h (:hover @state))
                               (swap! state assoc :hover h)
                               (render!)))))
      (.addEventListener canvas "mouseleave"
                         (fn [_] (swap! state assoc :hover nil) (render!)))
      (.addEventListener canvas "click"
                         (fn [^js e]
                           (let [h (hex-at e)]
                             (when ((:board-set @state) h)
                               (cond
                                 (or (.-paintWalls controls) (.-shiftKey e))
                                 (swap! state update :walls #(if (% h) (disj % h) (conj % h)))

                                 (.-altKey e) (swap! state assoc :source h)
                                 :else        (swap! state assoc :target h))
                               (render!)))))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "orientation" #js ["pointy" "flat"]) (.onChange regenerate!))
          (-> (.add controls "shape" #js ["hexagon" "rectangle" "triangle" "parallelogram"])
              (.onChange regenerate!))
          (-> (.add controls "mapSize" 2 10 1) (.onChange regenerate!))
          (-> (.add controls "overlay" #js ["path" "range" "ring" "spiral" "line" "reachable" "visible" "none"])
              (.onChange render!))
          (-> (.add controls "algorithm" #js ["a-star" "dijkstra" "greedy" "breadth-first"])
              (.onChange render!))
          (-> (.add controls "radius" 0 8 1) (.onChange render!))
          (-> (.add controls "labels" #js ["none" "axial" "cube" "offset" "doubled"])
              (.onChange render!))
          (-> (.add controls "paintWalls") (.onChange render!))
          (-> (.add controls "showSearch") (.onChange render!))
          (.add #js {:clearWalls (fn [] (swap! state assoc :walls #{}) (render!))} "clearWalls")))
      {:start (fn [] (resize!) (render!))
       :stop  (fn [])})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "hex")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
