(ns allgo.demo.city-agents
  "Canvas demo for `allgo.procedural.city.agents`: a town grown by
  agents, after Lechner, Watson, Wilensky and Felsen.

  Colored as their figures are: gray roads, red residential and blue
  commercial land, lighter where it is built denser. Nothing plans it.
  Extenders lay roads out to land the network does not serve, connectors
  close loops where the way round is too long, and developers build
  where building raises the value of the land -- commercial along the
  busy roads, housing spreading out where the roads are fewer.

  Speculators, orange, are the paper's future work, and their rules this
  code's own: each reserves a tract outlined in white at the edge of
  town, tinted while it is held, and sells it when the roads arrive --
  for a residential district, red only; a park, green, which no road
  crosses; or a campus, gold buildings on pale grounds.

  `deviation` is the paper's grid constraint: at 0 every road keeps to a
  grid line, a little more loosens the grid, and from half the spacing up
  the town grows organic. `k` is the connectors': low for a net of loops,
  high for cul-de-sacs. `speculators` at 0 is the paper's town as it
  was built. Changes take effect on `restart`."
  (:require [allgo.procedural.city.agents :as agents]
            ["lil-gui" :default GUI]))

(def ^:private ^js controls
  #js {:seed      4
       :deviation 1
       :k         2.5
       :speculators 3
       :speed     3})

(defn- color
  "A patch's color: water, road, a building by use and density, a park
  or a campus's grounds, or the land by its height -- tinted, if a
  speculator holds it, toward what it is held for."
  [world tracts k]
  (let [use (aget ^js (:use world) k)
        t (/ (dec (aget ^js (:density world) k)) 4.0)
        {:keys [kind sold?]} (get tracts (aget ^js (:tract world) k))
        c (cond
            (= 1 (aget ^js (:road world) k)) [128 128 132]
            (= 1 (aget ^js (:water world) k)) [58 98 145]
            (= use agents/residential) [(+ 175 (* 75 t)) (+ 50 (* 170 t)) (+ 50 (* 170 t))]
            (= use agents/commercial) [(+ 50 (* 170 t)) (+ 90 (* 150 t)) (+ 200 (* 50 t))]
            (= use agents/institution) [(+ 200 (* 50 t)) (+ 160 (* 60 t)) 40]
            (and sold? (= kind :park)) [70 170 70]
            (and sold? (= kind :campus)) [150 190 120]
            :else (let [e (min 1.0 (/ (aget ^js (:elev world) k) 60.0))]
                    [(+ 62 (* 50 e)) (+ 100 (* 50 e)) (+ 56 (* 38 e))]))]
    (if (and kind (not sold?) (zero? use))
      (mapv #(+ (* 0.6 %1) (* 0.4 %2)) c (case kind :park [90 220 90] :campus [230 200 60] [240 150 150]))
      c)))

(defn- outline!
  "White edges around each tract, on a canvas of `cell`-sized patches
  from `(x0, y0)`."
  [^js ctx world tracts x0 y0 cell]
  (let [n (:w world)
        tract (:tract world)
        id-at (fn [i j] (if (and (< -1 i n) (< -1 j n)) (aget ^js tract (+ i (* n j))) 0))]
    (set! (.-strokeStyle ctx) "rgba(255,255,255,0.85)")
    (set! (.-lineWidth ctx) 1)
    (.beginPath ctx)
    (doseq [{:keys [id patches]} (vals tracts)
            [i j] patches
            :when (= id (id-at i j))
            :let [x (+ x0 (* cell i)) y (+ y0 (* cell (- n 1 j)))]]
      ;; Row 0 is at the bottom, so north is up the canvas.
      (when (not= id (id-at (inc i) j)) (.moveTo ctx (+ x cell) y) (.lineTo ctx (+ x cell) (+ y cell)))
      (when (not= id (id-at (dec i) j)) (.moveTo ctx x y) (.lineTo ctx x (+ y cell)))
      (when (not= id (id-at i (inc j))) (.moveTo ctx x y) (.lineTo ctx (+ x cell) y))
      (when (not= id (id-at i (dec j))) (.moveTo ctx x (+ y cell)) (.lineTo ctx (+ x cell) (+ y cell))))
    (.stroke ctx)))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        grid (js/document.createElement "canvas")
        gctx (.getContext grid "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        world (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [w (agents/world {:seed (.-seed controls)
                                     :deviation (.-deviation controls)
                                     :k (.-k controls)
                                     :speculators (.-speculators controls)})]
                (set! (.-width grid) (:w w))
                (set! (.-height grid) (:w w))
                (reset! world w)))
            (draw! []
              (let [w @world
                    n (:w w)
                    tracts @(:tracts w)
                    img (.createImageData gctx n n)
                    data (.-data img)]
                (dotimes [j n]
                  (dotimes [i n]
                    (let [[r g b] (color w tracts (+ i (* n j)))
                          ;; Row 0 is the south edge: drawn at the bottom.
                          o (* 4 (+ i (* n (- n 1 j))))]
                      (aset data o r) (aset data (+ o 1) g) (aset data (+ o 2) b) (aset data (+ o 3) 255))))
                (.putImageData gctx img 0 0)
                (let [cw (.-width canvas) ch (.-height canvas)
                      side (min cw ch)
                      x0 (* 0.5 (- cw side)) y0 (* 0.5 (- ch side))]
                  (set! (.-fillStyle ctx) "#1a1f28")
                  (.fillRect ctx 0 0 cw ch)
                  (set! (.-imageSmoothingEnabled ctx) false)
                  (.drawImage ctx grid x0 y0 side side)
                  (outline! ctx w tracts x0 y0 (/ side n))
                  ;; The agents themselves, as dots.
                  (let [cell (/ side n)]
                    (doseq [{:keys [kind] [i j] :at} @(:agents w)]
                      (set! (.-fillStyle ctx) (case kind :extender "#ffe066" :connector "#ff7bd5"
                                                    :speculator "#ff9a3c" :residential "#ffffff" "#9fe0ff"))
                      (.fillRect ctx (+ x0 (* cell i) (* 0.25 cell)) (+ y0 (* cell (- n 1 j)) (* 0.25 cell))
                                 (* 0.5 cell) (* 0.5 cell)))))
                (let [u (frequencies (array-seq (:use w)))
                      sold (frequencies (map :kind (filter :sold? (vals tracts))))]
                  (set! (.-textContent caption)
                        (str "tick " @(:ticks w) "  |  " (count @(:roads w)) " roads, "
                             (get u agents/residential 0) " residential and "
                             (get u agents/commercial 0) " commercial patches  |  sold "
                             (get sold :residential 0) " districts, " (get sold :park 0) " parks, "
                             (get sold :campus 0) " campus  |  "
                             "agents: extenders yellow, connectors pink, developers white and blue, speculators orange")))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (dotimes [_ (.-speed controls)] (agents/tick! @world))
                (draw!)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h))))]
      (let [gui (GUI. #js {:container container})]
        (.add gui controls "seed" 1 99 1)
        (.add gui controls "deviation" 0 6 1)
        (.add gui controls "k" 1.2 6 0.1)
        (.add gui controls "speculators" 0 8 1)
        (.add gui controls "speed" 1 20 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "city-agents")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
