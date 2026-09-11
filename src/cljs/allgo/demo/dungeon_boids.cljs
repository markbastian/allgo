(ns allgo.demo.dungeon-boids
  "`allgo.simulation.boids` loose inside `allgo.procedural.dungeons`.

  The dungeon is generated once and its walls become obstacles, so the
  flock has to find its way through the corridors rather than across open
  canvas. It is a fair test of both halves: the boids' avoidance is the
  same GJK query used against a handful of circles in the 2D flocking
  demo, now answering against a few hundred wall segments, and the
  dungeon's corridors have to actually be navigable for the flock to
  spread beyond the room it started in.

  Watching it, the flock behaves rather like water: it pools in the hub
  rooms, files single-column down the corridors, and splits at junctions."
  (:require [allgo.demo.dungeon :as dungeon-demo]
            [allgo.demo.fps :as fps]
            [allgo.procedural.dungeons :as dungeons]
            [allgo.simulation.boids :as boids]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private background "#05070d")

(def ^:private controls
  #js {:boids       50
       :rooms       45
       :separation  1.4
       :alignment   1.2
       :cohesion    0.9
       :perception  34
       :speed       1.7
       :avoidance   2.6
       :wallMargin  9
       :trails      true
       :showDungeon true})

(defn- dungeon-config []
  ;; Wider corridors and fewer rooms than the standalone dungeon demo: the
  ;; flock needs passages several boid-widths across, and a 150-room
  ;; dungeon scaled into a 480px tile leaves corridors about as wide as a
  ;; single boid.
  {:room-count     (.-rooms controls)
   :radius         70.0
   :width-mean     22.0
   :width-sd       7.0
   :height-mean    22.0
   :height-sd      7.0
   :corridor-width 5
   :min-hubs       4})

;; ---------------------------------------------------------------------------
;; Turning the dungeon's walls into obstacles

(defn- wall-shell
  "The ring of solid tiles touching the floor.

  Only this shell can be hit, so the rest of the void is not worth
  modelling -- and there is no way to enumerate it anyway, the void being
  unbounded."
  [grid]
  (into #{}
        (comp (mapcat (fn [[tx ty]]
                        (for [dx [-1 0 1] dy [-1 0 1]
                              :when (not (and (zero? dx) (zero? dy)))]
                          [(+ tx dx) (+ ty dy)])))
              (remove #(contains? grid %)))
        (keys grid)))

(defn- wall-obstacles
  "One box per shell tile, in pixel space.

  Merging each row into long runs is the obvious optimisation and it makes
  things slower. Obstacles are culled against a bounding sphere before the
  GJK query, and a twenty-tile run has a sphere wide enough to catch most
  of the flock, so the cull stops doing anything and GJK runs against the
  whole wall. Measured over 20 steps with 60 boids: 65 merged runs cost
  13.5ms a step, 138 single tiles cost 11.8ms."
  [grid [scale ox oy]]
  (mapv (fn [[tx ty]]
          (boids/box-obstacle
           [(+ ox (* tx scale)) (+ oy (* ty scale))]
           [(+ ox (* (inc tx) scale)) (+ oy (* (inc ty) scale))]))
        (wall-shell grid)))

(defn- spawn
  "`n` boids dropped on random floor tiles, so none start inside a wall."
  [grid n [scale ox oy] max-speed]
  (let [floor (vec (keys grid))]
    (vec (repeatedly
          n
          #(let [[tx ty] (rand-nth floor)
                 heading (* 2 math/PI (rand))
                 speed   (* max-speed (+ 0.5 (rand 0.5)))]
             {:pos [(+ ox (* (+ 0.5 tx) scale)) (+ oy (* (+ 0.5 ty) scale))]
              :vel [(* speed (math/cos heading)) (* speed (math/sin heading))]})))))

(defn- params [obstacles]
  {:edges             :bounce
   :obstacles         obstacles
   :boid-radius       1.5
   :avoid-radius      (.-wallMargin controls)
   :avoid-weight      (.-avoidance controls)
   :separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (* 0.4 (.-perception controls))
   :max-speed         (.-speed controls)})

;; ---------------------------------------------------------------------------

(defn- draw-dungeon! [^js ctx grid [scale ox oy]]
  (let [side (js/Math.ceil (+ scale 0.5))]
    (doseq [[kind tiles] (group-by val grid)]
      (set! (.-fillStyle ctx) (case kind
                                :main "#1b2338"
                                :hallway "#161d30"
                                "#111726"))
      (doseq [[[tx ty] _] tiles]
        (.fillRect ctx (+ ox (* tx scale)) (+ oy (* ty scale)) side side)))))

(defn- draw-boid! [^js ctx {[x y] :pos [vx vy] :vel}]
  (let [heading (math/atan2 vy vx)
        hue     (mod (+ 180 (math/to-degrees heading)) 360)]
    (set! (.-fillStyle ctx) (str "hsl(" hue " 85% 64%)"))
    (.save ctx)
    (.translate ctx x y)
    (.rotate ctx heading)
    (.beginPath ctx)
    (.moveTo ctx 5 0)
    (.lineTo ctx -3.5 2.5)
    (.lineTo ctx -1.5 0)
    (.lineTo ctx -3.5 -2.5)
    (.closePath ctx)
    (.fill ctx)
    (.restore ctx)))

(defn- draw! [^js ctx flock grid xf [w h]]
  ;; Trails are painted as a translucent wash rather than a clear, which
  ;; is what makes the flow through a corridor legible as a stream.
  (set! (.-fillStyle ctx) (if (.-trails controls) "rgba(5, 7, 13, 0.3)" background))
  (.fillRect ctx 0 0 w h)
  (when (.-showDungeon controls)
    (draw-dungeon! ctx grid xf))
  (doseq [boid flock] (draw-boid! ctx boid)))

(defn init! [^js container]
  (let [canvas    (js/document.createElement "canvas")
        ctx       (.getContext canvas "2d")
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom {:bounds [0 0] :flock [] :grid {} :xf [1 0 0] :obstacles []})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(rebuild! []
              ;; The fit transform depends on the viewport, so the walls are
              ;; rebuilt on resize as well as on regenerate -- obstacles live
              ;; in pixel space, not tile space.
              (let [{:keys [bounds dungeon]} @state]
                (when (and dungeon (pos? (first bounds)))
                  (let [xf (dungeon-demo/fit (:grid dungeon) bounds)]
                    (swap! state assoc
                           :grid (:grid dungeon)
                           :xf xf
                           :obstacles (wall-obstacles (:grid dungeon) xf))))))
            (reseed! []
              (let [{:keys [grid xf]} @state]
                (when (seq grid)
                  (swap! state assoc :flock
                         (spawn grid (.-boids controls) xf (.-speed controls))))))
            (generate! []
              (swap! state assoc :dungeon (dungeons/generate (dungeon-config)))
              (rebuild!)
              (reseed!))
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
                  (set! (.-fillStyle ctx) background)
                  (.fillRect ctx 0 0 w h)
                  (swap! state assoc :bounds [w h])
                  (rebuild!)
                  (reseed!))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      {:keys [bounds flock grid xf obstacles]} @state
                      flock (boids/step flock bounds (params obstacles))]
                  (swap! state assoc :flock flock)
                  (draw! ctx flock grid xf bounds)
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (.add controls "boids" 10 200 10)
          (-> (.add controls "rooms" 20 90 5) (.onFinishChange generate!))
          (.add controls "separation" 0 3 0.1)
          (.add controls "alignment" 0 3 0.1)
          (.add controls "cohesion" 0 3 0.1)
          (.add controls "perception" 10 80 2)
          (.add controls "speed" 0.4 4 0.1)
          (.add controls "avoidance" 0 6 0.1)
          (.add controls "wallMargin" 3 24 1)
          (.add controls "trails")
          (.add controls "showDungeon")
          (.add #js {:regenerate generate!} "regenerate")
          (.add #js {:reseed reseed!} "reseed")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "dungeon-boids")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
