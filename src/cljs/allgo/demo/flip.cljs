(ns allgo.demo.flip
  "Canvas demo for `allgo.physics.flip`, after Ten Minute Physics 18.

  A dam break. The grid from demo 17 is still here doing the pressure
  solve, but nothing is advected on it any more -- the fluid is the
  particles, and they simply move. That is the whole point: advection
  stops being an interpolation and starts being transport, so the detail
  survives.

  `flipRatio` is the parameter worth playing with. At 0 the particles take
  the grid velocity outright, every particle in a cell is averaged with its
  neighbours, and the water turns to treacle. At 1 they take only the
  grid's change, nothing is averaged away, and the splash keeps every eddy
  along with enough noise to look like it is boiling. 0.9 is the usual
  compromise.

  `compensateDrift` is easy to miss until it is off: turn it off and watch
  the body of water slowly lose volume as the particles crowd, because
  nothing else notices a cell holding more than its share.

  `separate` is what keeps the particles evenly spread; without it they
  clump into the middle of cells and the pressure solve reads the gaps as
  empty space."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.flip :as flip]
            ["lil-gui" :default GUI]))

(def ^:private ^js controls
  #js {:resolution      60
       :flipRatio       0.9
       :iterations      50
       :separate        true
       :compensateDrift true
       :gravity         9.81
       :show            "particles"
       :obstacle        true})

(def ^:private obstacle-radius 0.12)

(defn- build
  "A tank sized so its cells are square in the viewport."
  [[w h]]
  (let [ny   (.-resolution controls)
        nx   (max 8 (js/Math.round (* ny (/ w (max 1 h)))))
        cell (/ 1.0 ny)
        ;; The reference's ratio: a particle is 0.3 of a cell, which is
        ;; what makes a full cell hold about three of them.
        r    (* 0.3 cell)
        ;; A staggered lattice at this radius; leave room to spare.
        est  (js/Math.ceil (* 1.6 (/ (* 0.5 nx cell 0.9 ny cell) (* 2 r r))))
        f    (-> (flip/flip-fluid {:nx nx :ny ny :h cell :particle-radius r
                                   :max-particles est})
                 (flip/close-border!))]
    (flip/fill-block! f (* 2 cell) (* 2 cell) (* 0.4 nx cell) (* 0.9 ny cell))
    {:fluid f :cell cell
     :obstacle [(* 0.65 nx cell) (* 0.4 ny cell) (* obstacle-radius ny cell)]}))

;; ---------------------------------------------------------------------------
;; Colour

(defn- sci-colour
  "Blue through green to red, the same ramp as the Eulerian demo."
  [v lo hi]
  (let [t (max 0.0 (min 0.999 (/ (- v lo) (max 1e-9 (- hi lo)))))
        m 0.25
        n (js/Math.floor (/ t m))
        s (/ (- t (* n m)) m)
        [r g b] (case n
                  0 [0.0 s 1.0]
                  1 [0.0 1.0 (- 1.0 s)]
                  2 [s 1.0 0.0]
                  [1.0 (- 1.0 s) 0.0])]
    (str "rgb(" (js/Math.round (* 255 r)) ","
         (js/Math.round (* 255 g)) "," (js/Math.round (* 255 b)) ")")))

(defn- draw!
  [^js ctx {:keys [fluid cell obstacle]} [w h]]
  (let [{:keys [nx ny ^floats pos ^floats vel]} fluid
        n     (flip/particle-count fluid)
        scale (/ w (* nx cell))
        ;; The grid runs bottom-up; the canvas runs top-down.
        sy    (fn [y] (- h (* y scale)))
        mode  (.-show controls)
        r     (* (:particle-radius fluid) scale)]
    (set! (.-fillStyle ctx) "#05070d")
    (.fillRect ctx 0 0 w h)
    (when (= "grid" mode)
      (let [{:keys [^ints cell-type]} fluid]
        (dotimes [i nx]
          (dotimes [j ny]
            (let [t (aget cell-type (+ (* i ny) j))]
              (when (not= 1 t)                          ; skip air
                (set! (.-fillStyle ctx) (if (= 2 t) "#28304a" "#123055"))
                (.fillRect ctx (* i cell scale) (- h (* (inc j) cell scale))
                           (js/Math.ceil (* cell scale))
                           (js/Math.ceil (* cell scale)))))))))
    (when (.-obstacle controls)
      (let [[ox oy orad] obstacle]
        (set! (.-fillStyle ctx) "#4f7ac0")
        (.beginPath ctx)
        (.arc ctx (* ox scale) (sy oy) (* orad scale) 0 (* 2 js/Math.PI))
        (.fill ctx)))
    (dotimes [i n]
      (let [b (* 2 i)]
        (set! (.-fillStyle ctx)
              (if (= "speed" mode)
                (sci-colour (js/Math.hypot (aget vel b) (aget vel (inc b))) 0.0 3.0)
                "#5fa8e8"))
        (.beginPath ctx)
        (.arc ctx (* (aget pos b) scale) (sy (aget pos (inc b)))
              (max 1.0 r) 0 (* 2 js/Math.PI))
        (.fill ctx)))))

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
                  (swap! state assoc :bounds [w h])
                  (rebuild!))))
            (move-obstacle! [^js e]
              (let [{:keys [fluid cell bounds obstacle]} @state
                    rect (.getBoundingClientRect canvas)
                    [w h] bounds
                    scale (/ w (* (:nx fluid) cell))
                    fx (/ (- (.-clientX e) (.-left rect)) scale)
                    fy (/ (- h (- (.-clientY e) (.-top rect))) scale)]
                (when (.-obstacle controls)
                  (swap! state assoc :obstacle [fx fy (nth obstacle 2)]))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      {:keys [fluid obstacle] :as st} @state]
                  (flip/step! fluid
                              {:gravity    (- (.-gravity controls))
                               :iterations (.-iterations controls)
                               :flip-ratio (.-flipRatio controls)
                               :separation (if (.-separate controls) 2 0)
                               :drift?     (.-compensateDrift controls)
                               :obstacle   (when (.-obstacle controls) obstacle)})
                  (draw! ctx st (:bounds st))
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.addEventListener canvas "mousemove"
                         (fn [^js e] (when (pos? (.-buttons e)) (move-obstacle! e))))
      (.addEventListener canvas "mousedown" move-obstacle!)
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "resolution" 30 110 5) (.onFinishChange rebuild!))
          (.add controls "flipRatio" 0 1 0.05)
          (.add controls "iterations" 1 100 1)
          (.add controls "separate")
          (.add controls "compensateDrift")
          (.add controls "gravity" 0 20 0.5)
          (.add controls "show" #js ["particles" "speed" "grid"])
          (.add controls "obstacle")
          (.add #js {:reset rebuild!} "reset")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "flip")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
