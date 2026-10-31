(ns allgo.demo.wumpus
  "Canvas demo for `allgo.logic.wumpus`: the hybrid agent of Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 7, exploring a wumpus world it cannot see.

  Each square is drawn as the agent knows it: green where it has proved
  there is no pit and no live wumpus, red where it has proved a pit, a W
  where it has proved the wumpus must be, and dark where it knows nothing
  -- every proof a DPLL query on its knowledge base, which the caption
  counts. It walks only on squares proved safe, by A*; when none is left
  unvisited it shoots where the wumpus might be, and when that fails it
  takes a risk on a square not proved deadly -- and sometimes dies of it.
  `reveal` shows the cave as it is."
  (:require [allgo.logic.wumpus :as wumpus]
            ["lil-gui" :default GUI]
            [clojure.math :as math]
            [clojure.string :as str]))

(def ^:private ^js controls
  #js {:seed   6
       :speed  3
       :reveal false})

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [w (wumpus/world (.-seed controls))
                    ss (wumpus/run-seq w)]
                (reset! model {:world w :rest ss :st (first ss) :clock 0 :done? false})))
            (advance! []
              ;; `speed` actions a second of real time, however often
              ;; frames come.
              (let [now (js/performance.now)
                    {:keys [clock last-time]} @model
                    clock (+ clock (* (.-speed controls) 0.001 (- now (or last-time now))))]
                (swap! model assoc :last-time now)
                (loop [clock clock]
                  (let [{:keys [rest done?]} @model]
                    (if (and (not done?) (>= clock 1.0))
                      (let [nxt (next rest)]
                        (swap! model assoc :rest (or nxt rest) :st (first (or nxt rest))
                               :done? (or (nil? nxt) (nil? (next nxt))))
                        (recur (- clock 1.0)))
                      (swap! model assoc :clock (min clock 1.0)))))))
            (draw! []
              (let [{:keys [world st]} @model
                    {:keys [state memory percept action]} st
                    {:keys [safe pits no-wumpus visited queries kb plan]} memory
                    cw (.-width canvas) ch (.-height canvas)
                    side (min (- cw 40) (- ch 170))
                    cell (/ side wumpus/size)
                    x0 (* 0.5 (- cw side)) y0 60
                    px (fn [[x y]] [(+ x0 (* cell (dec x))) (+ y0 (* cell (- wumpus/size y)))])
                    possible-w (remove no-wumpus (for [y (range 1 5) x (range 1 5)] [x y]))
                    proved-w (when (= 1 (count possible-w)) (first possible-w))
                    reveal? (.-reveal controls)
                    text (fn [s [x y] size color]
                           (set! (.-fillStyle ctx) color)
                           (set! (.-font ctx) (str "bold " (int size) "px sans-serif"))
                           (set! (.-textAlign ctx) "center") (set! (.-textBaseline ctx) "middle")
                           (.fillText ctx s x y))]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (doseq [y (range 1 5) x (range 1 5)
                        :let [sq [x y] [sx sy] (px sq)
                              cx (+ sx (* 0.5 cell)) cy (+ sy (* 0.5 cell))]]
                  (set! (.-fillStyle ctx)
                        (cond (contains? pits sq) "#5a1d1d"
                              (contains? visited sq) "#2e5a3a"
                              (contains? safe sq) "#1f3d28"
                              :else "#161b25"))
                  (.fillRect ctx (+ sx 2) (+ sy 2) (- cell 4) (- cell 4))
                  (when (contains? pits sq) (text "PIT" [cx (+ cy (* 0.3 cell))] (* 0.14 cell) "#ff8080"))
                  (when (= sq proved-w) (text "W" [cx (+ cy (* 0.3 cell))] (* 0.2 cell) "#ffb347"))
                  (when (and (not (contains? safe sq)) (not (contains? pits sq)) (not= sq proved-w))
                    (text "?" [cx cy] (* 0.3 cell) "#3a4356"))
                  ;; The cave as it is.
                  (when reveal?
                    (let [things (cond-> []
                                   (contains? (:pits world) sq) (conj "pit")
                                   (= sq (:wumpus world)) (conj (if (:wumpus-alive? state) "wumpus" "dead wumpus"))
                                   (and (= sq (:gold world)) (not (:gold? state))) (conj "gold"))]
                      (when (seq things)
                        (text (str/join " " things) [cx (+ sy (* 0.15 cell))] (* 0.1 cell) "#c9d2e0")))))
                ;; The agent: a triangle pointing the way it faces.
                (let [[sx sy] (px (:pos state))
                      cx (+ sx (* 0.5 cell)) cy (+ sy (* 0.5 cell))
                      a (case (:dir state) :east 0 :north (/ math/PI -2) :west math/PI :south (/ math/PI 2))
                      r (* 0.22 cell)
                      pt (fn [da k] [(+ cx (* k r (math/cos (+ a da)))) (+ cy (* k r (math/sin (+ a da))))])]
                  (set! (.-fillStyle ctx) (if (:alive? state) "#ffffff" "#ff4040"))
                  (.beginPath ctx)
                  (let [[ax ay] (pt 0 1.2) [bx by] (pt 2.4 1) [qx qy] (pt -2.4 1)]
                    (.moveTo ctx ax ay) (.lineTo ctx bx by) (.lineTo ctx qx qy))
                  (.closePath ctx) (.fill ctx)
                  (when (:gold? state) (text "$" [cx (- cy (* 0.3 cell))] (* 0.18 cell) "#ffd700")))
                (set! (.-textContent caption)
                      (str "felt: " (or (not-empty (str/join " " (for [[k v] percept :when v] (name k)))) "nothing")
                           "  |  then: " (if action (name action) "-")
                           (when (seq plan) (str "  |  plan: " (str/join " " (map name plan))))
                           "\nscore " (:score state)
                           "  |  " (count kb) " clauses, " queries " DPLL queries"
                           "  |  arrow " (if (:arrow? state) "unused" "spent")
                           (cond (:out? state) (if (:gold? state) "  |  climbed out with the gold" "  |  climbed out empty-handed")
                                 (not (:alive? state)) (str "  |  killed by the " (if (contains? (:pits world) (:pos state)) "pit" "wumpus"))
                                 :else "")
                           "\ngreen: proved safe  red: proved a pit  W: proved the wumpus  ?: unknown"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "seed" 1 500 1) restart!)
        (.add gui controls "speed" 0.5 20 0.5)
        (.add gui controls "reveal")
        (.add gui #js {:restart restart!} "restart")
        (.add gui #js {:next-cave (fn [] (set! (.-seed controls) (inc (.-seed controls))) (restart!))} "next-cave"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "wumpus")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
