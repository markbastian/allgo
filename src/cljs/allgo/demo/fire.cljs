(ns allgo.demo.fire
  "Canvas demo for `allgo.physics.fire`, after Ten Minute Physics 21.

  The solver is the wind tunnel from demo 17 with nothing added to it.
  What makes this fire rather than smoke is three passes over the scalar
  field the tunnel was already carrying as dye: hot cells are pulled
  upward, hot cells cool, and vortices are spawned at the source.

  `lift` against `buoyancy` is the pair worth playing with -- the first is
  how fast a cell at full heat wants to rise, the second how quickly it
  gets there. Turn `lift` to zero and the fire stops rising entirely,
  which is the point: there is no gravity in this scene, and buoyancy is
  the only thing moving anything.

  `fireCooling` against `smokeCooling` is what separates flame from
  smoke. Set them equal and the plume dims uniformly instead of turning
  grey at a height.

  `vortexRate` is the honest one. Turn it to zero and watch the plume go
  flat and sheet-like: a grid at this resolution smooths its own eddies
  away within a few cells, and the curl has to be put back by hand.

  Drag to move the burning disc. `show` swaps the fire palette for the
  raw temperature, the flow speed, or the vortices themselves."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.fire :as fire]
            [allgo.physics.vortex :as vortex]
            ["lil-gui" :default GUI]))

(def ^:private ^js controls
  #js {:resolution   140
       :source       "disc"
       :iterations   10
       :lift         3.0
       :buoyancy     6.0
       :fireCooling  1.2
       :smokeCooling 0.3
       :vortexRate   2000
       :vortexRadius 0.05
       :vortexOmega  20.0
       :show         "fire"})

(defn- emitters [{:keys [ny h]} [ox oy]]
  (let [disc {:kind :disc :x ox :y oy :radius (* 0.09 (* ny h))}
        floor {:kind :floor :rows (max 2 (quot ny 40))}]
    (case (.-source controls)
      "disc"  [disc]
      "floor" [floor]
      "both"  [disc floor]
      [disc])))

(defn- build [[w h]]
  (let [ny   (.-resolution controls)
        nx   (max 8 (js/Math.round (* ny (/ w (max 1 h)))))
        cell (/ 1.0 ny)
        f    (fire/fire {:nx nx :ny ny :h cell :vortices 400})]
    {:fire f :cell cell
     :obstacle [(* 0.5 (:nx f) cell) (* 0.25 (:ny f) cell)]}))

;; ---------------------------------------------------------------------------
;; Colour

(defn- fire-colour
  "Black through dark smoke to red, orange and yellow.

  The three bands are not decoration: the first covers everything below
  the smoke threshold, so the palette changes character at exactly the
  temperature the cooling rate does."
  [t]
  (let [t (max 0.0 (min 1.0 t))]
    (cond
      (< t 0.3) (let [s (/ t 0.3) g (js/Math.round (* 255 0.2 s))] [g g g])
      (< t 0.5) (let [s (/ (- t 0.3) 0.2)]
                  [(js/Math.round (* 255 (+ 0.2 (* 0.8 s))))
                   (js/Math.round (* 255 0.1))
                   (js/Math.round (* 255 0.1))])
      :else     (let [s (/ (- t 0.5) 0.48)]
                  [255 (js/Math.round (* 255 (min 1.0 s))) 0]))))

(defn- sci-colour [v lo hi]
  (let [t (max 0.0 (min 0.999 (/ (- v lo) (max 1e-9 (- hi lo)))))
        n (js/Math.floor (/ t 0.25))
        s (/ (- t (* n 0.25)) 0.25)
        [r g b] (case n
                  0 [0.0 s 1.0]
                  1 [0.0 1.0 (- 1.0 s)]
                  2 [s 1.0 0.0]
                  [1.0 (- 1.0 s) 0.0])]
    [(js/Math.round (* 255 r)) (js/Math.round (* 255 g)) (js/Math.round (* 255 b))]))

(defn- draw! [^js ctx {:keys [fire cell obstacle]} [w h]]
  (let [{:keys [nx ny ^floats smoke ^floats u ^floats v vortices]} fire
        img  (.createImageData ctx nx ny)
        data (.-data img)
        mode (.-show controls)]
    (dotimes [j ny]
      (dotimes [i nx]
        (let [k (+ (* i ny) j)
              [r g b]
              (case mode
                "temperature" (let [c (js/Math.round (* 255 (max 0.0 (min 1.0 (aget smoke k)))))]
                                [c c c])
                "speed" (sci-colour (js/Math.hypot (aget u k) (aget v k)) 0.0 3.0)
                (fire-colour (aget smoke k)))
              ;; ImageData runs top-down; the grid runs bottom-up.
              px (* 4 (+ i (* (- ny 1 j) nx)))]
          (aset data px r)
          (aset data (+ px 1) g)
          (aset data (+ px 2) b)
          (aset data (+ px 3) 255))))
    (.putImageData ctx img 0 0)
    (.drawImage ctx (.-canvas ctx) 0 0 nx ny 0 0 w h)
    (when (= "vortices" mode)
      (let [scale (/ w (* nx cell))]
        (set! (.-strokeStyle ctx) "#6fd0ff")
        (dotimes [k (vortex/vortex-count vortices)]
          (let [{:keys [x y omega]} (vortex/vortex vortices k)]
            (.beginPath ctx)
            (.arc ctx (* x scale) (- h (* y scale))
                  (* (.-vortexRadius controls) scale) 0 (* 2 js/Math.PI))
            (set! (.-lineWidth ctx) (if (pos? omega) 2 1))
            (.stroke ctx)))))
    (when (not= "floor" (.-source controls))
      (let [scale (/ w (* nx cell))
            [ox oy] obstacle]
        (set! (.-strokeStyle ctx) "rgba(255,255,255,0.25)")
        (set! (.-lineWidth ctx) 1)
        (.beginPath ctx)
        (.arc ctx (* ox scale) (- h (* oy scale))
              (* 0.09 ny cell scale) 0 (* 2 js/Math.PI))
        (.stroke ctx)))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state  (atom {:bounds [0 0]})]
    (set! (.-style canvas) "display:block;border-radius:8px;cursor:crosshair")
    (.appendChild container canvas)
    (letfn [(rebuild! []
              (let [{:keys [bounds]} @state]
                (when (pos? (first bounds))
                  (swap! state merge (build bounds)))))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (set! (.-imageSmoothingEnabled ctx) false)
                  (swap! state assoc :bounds [w h])
                  (rebuild!))))
            (move-source! [^js e]
              (let [{:keys [fire cell bounds]} @state
                    rect (.getBoundingClientRect canvas)
                    [w h] bounds
                    scale (/ w (* (:nx fire) cell))]
                (swap! state assoc :obstacle
                       [(/ (- (.-clientX e) (.-left rect)) scale)
                        (/ (- h (- (.-clientY e) (.-top rect))) scale)])))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      {:keys [fire obstacle] :as st} @state]
                  (fire/step! fire
                              {:iterations    (.-iterations controls)
                               :lift          (.-lift controls)
                               :buoyancy      (.-buoyancy controls)
                               :fire-cooling  (.-fireCooling controls)
                               :smoke-cooling (.-smokeCooling controls)
                               :vortex-rate   (.-vortexRate controls)
                               :vortex-radius (.-vortexRadius controls)
                               :vortex-omega  (.-vortexOmega controls)
                               :emitters      (emitters fire obstacle)})
                  (draw! ctx st (:bounds st))
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.addEventListener canvas "mousemove"
                         (fn [^js e] (when (pos? (.-buttons e)) (move-source! e))))
      (.addEventListener canvas "mousedown" move-source!)
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "resolution" 60 220 10) (.onFinishChange rebuild!))
          (-> (.add controls "source" #js ["disc" "floor" "both"]) (.onChange rebuild!))
          (.add controls "iterations" 1 40 1)
          (.add controls "lift" 0 8 0.1)
          (.add controls "buoyancy" 0 20 0.5)
          (.add controls "fireCooling" 0.1 4 0.1)
          (.add controls "smokeCooling" 0 2 0.05)
          (.add controls "vortexRate" 0 8000 100)
          (.add controls "vortexRadius" 0.01 0.15 0.005)
          (.add controls "vortexOmega" 0 60 1)
          (.add controls "show" #js ["fire" "temperature" "speed" "vortices"])
          (.add #js {:reset rebuild!} "reset")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "fire")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
