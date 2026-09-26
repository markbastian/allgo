(ns allgo.demo.fluid
  "Canvas demo for `allgo.physics.fluid`, after Ten Minute Physics 17.

  A wind tunnel with an obstacle you can drag. What is actually being
  simulated is a grid of velocities; the dye is only paint, carried along
  to make the flow visible.

  Worth doing: drop `iterations` to 2 and watch the flow go slack and
  compressible, because the projection no longer has time to drive the
  divergence out. Then turn `overRelax` off at the same iteration count
  and watch it get worse again -- the overshoot is worth several sweeps.

  Then leave it low and switch `solver` to multigrid, where an iteration
  means a V-cycle rather than a sweep. Eight of them beat two thousand
  sweeps on a 128 grid, because a sweep moves information one cell and a
  cycle moves it across the whole grid. `overRelax` does nothing except
  for Gauss-Seidel; it is that solver's own trick.
  `show` switches between the dye, the speed of the flow, and the pressure
  the projection worked out."
  (:require [allgo.demo.fps :as fps]
            [allgo.physics.fluid :as fluid]
            ["lil-gui" :default GUI]))

(def ^:private ^js controls
  #js {:resolution 96
       :speed      2.0
       :solver     "gauss-seidel"
       :iterations 40
       :overRelax  true
       :gravity    0.0
       :show       "dye"
       :obstacle   true})

(defn- build
  "A tunnel sized so its cells are square in the viewport."
  [[w h]]
  (let [ny  (.-resolution controls)
        nx  (max 8 (js/Math.round (* ny (/ w (max 1 h)))))
        cell (/ 1.0 ny)
        f   (-> (fluid/fluid nx ny cell 1000.0)
                (fluid/close-border! true))]
    (fluid/wind-tunnel! f (.-speed controls) 0.1)
    (when (.-obstacle controls)
      (fluid/disk! f (* 0.35 nx cell) (* 0.5 ny cell) (* 0.12 ny cell)))
    {:fluid f :cell cell
     :obstacle [(* 0.35 nx cell) (* 0.5 ny cell) (* 0.12 ny cell)]}))

;; ---------------------------------------------------------------------------
;; Color

(defn- sci-color
  "Blue through green to red, the usual scientific ramp."
  [v lo hi]
  (let [t (max 0.0 (min 0.999 (/ (- v lo) (max 1e-9 (- hi lo)))))
        m (/ 0.25 1.0)
        n (js/Math.floor (/ t m))
        s (/ (- t (* n m)) m)
        [r g b] (case n
                  0 [0.0 s 1.0]
                  1 [0.0 1.0 (- 1.0 s)]
                  2 [s 1.0 0.0]
                  [1.0 (- 1.0 s) 0.0])]
    [(js/Math.round (* 255 r)) (js/Math.round (* 255 g)) (js/Math.round (* 255 b))]))

(defn- draw!
  "One pixel per cell, written straight into an ImageData buffer.

  Per-cell fillRect at this resolution costs more than the simulation
  does; writing the bytes and blitting once does not."
  [^js ctx {:keys [fluid]} [w h]]
  (let [{:keys [nx ny ^js smoke ^js u ^js v ^js p ^js s]} fluid
        img   (.createImageData ctx nx ny)
        data  (.-data img)
        mode  (.-show controls)
        ;; Pressure needs its own range each frame; it has no fixed scale.
        [lo hi] (if (= "pressure" mode)
                  (loop [i 0 lo ##Inf hi ##-Inf]
                    (if (= i (* nx ny))
                      [lo hi]
                      (recur (inc i) (min lo (aget p i)) (max hi (aget p i)))))
                  [0.0 1.0])]
    (dotimes [j ny]
      (dotimes [i nx]
        (let [k (+ (* i ny) j)
              [r g b]
              (cond
                (zero? (aget s k)) [40 48 70]
                (= "pressure" mode) (sci-color (aget p k) lo hi)
                (= "speed" mode)
                (let [sp (js/Math.hypot (aget u k) (aget v k))]
                  (sci-color sp 0.0 (* 2.2 (.-speed controls))))
                :else (let [m (aget smoke k)
                            c (js/Math.round (* 255 (max 0.0 (min 1.0 m))))]
                        [c c c]))
              ;; ImageData runs top-down; the grid runs bottom-up.
              px (* 4 (+ i (* (- ny 1 j) nx)))]
          (aset data px r)
          (aset data (+ px 1) g)
          (aset data (+ px 2) b)
          (aset data (+ px 3) 255))))
    (.putImageData ctx img 0 0)
    ;; Blow the grid up to the viewport.
    (.drawImage ctx (.-canvas ctx) 0 0 nx ny 0 0 w h)))

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
            (move-obstacle! [^js e]
              (let [{:keys [fluid cell bounds obstacle]} @state
                    rect (.getBoundingClientRect canvas)
                    [w h] bounds
                    fx (* (/ (- (.-clientX e) (.-left rect)) w) (:nx fluid) cell)
                    fy (* (- 1.0 (/ (- (.-clientY e) (.-top rect)) h)) (:ny fluid) cell)]
                (when (.-obstacle controls)
                  (fluid/disk! fluid fx fy (nth obstacle 2))
                  (swap! state assoc :obstacle [fx fy (nth obstacle 2)]))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      {:keys [fluid] :as st} @state]
                  ;; The inflow is driven every frame: the projection would
                  ;; otherwise bleed it away within a few steps.
                  (fluid/wind-tunnel! fluid (.-speed controls) 0.1)
                  (fluid/step! fluid {:gravity (.-gravity controls)
                                      :solver (keyword (.-solver controls))
                                      :iterations (.-iterations controls)
                                      :over-relaxation (if (.-overRelax controls) 1.9 1.0)})
                  (draw! ctx st (:bounds st))
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.addEventListener canvas "mousemove"
                         (fn [^js e] (when (pos? (.-buttons e)) (move-obstacle! e))))
      (.addEventListener canvas "mousedown" move-obstacle!)
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "resolution" 32 160 8) (.onFinishChange rebuild!))
          (.add controls "speed" 0 6 0.1)
          (.add controls "solver" #js ["gauss-seidel" "conjugate-gradient" "multigrid"])
          (.add controls "iterations" 1 100 1)
          (.add controls "overRelax")
          (.add controls "gravity" -20 20 0.5)
          (.add controls "show" #js ["dye" "speed" "pressure"])
          (-> (.add controls "obstacle") (.onChange rebuild!))
          (.add #js {:reset rebuild!} "reset")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "fluid")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
