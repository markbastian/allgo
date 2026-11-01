(ns allgo.demo.localization
  "Canvas demo for `allgo.robotics.localization`: Monte Carlo localization,
  after Russell and Norvig, *Artificial Intelligence: A Modern Approach*,
  4th edition, chapter 26.

  The robot (white, its scan in faint rays) drives a loop round a room it
  has a map of, but it starts not knowing where it is: the particles (in
  orange) begin spread over the whole room. Each step they move as the
  robot says it moved, each with its own noise, are weighed by how well
  the scan they would see matches the one the robot saw, and are
  resampled. Within a few steps they gather into a few clusters -- places
  that look alike -- and then, as the robot moves on and the look-alikes
  stop agreeing, into one, around the robot. A few percent are thrown
  anywhere at random each step, so a filter that has settled on the
  wrong place can still recover."
  (:require [allgo.random :as random]
            [allgo.robotics.localization :as mcl]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:particles 1500 :seed 6 :speed 8})

(def ^:private route
  ;; A closed loop across the open top of the room: east, a right turn of
  ;; exactly a quarter, down, west, up. The turns are timed to make
  ;; quarter turns, so that without noise the loop would close.
  (let [w (- (/ (/ math/PI 2) 1.6))]
    (cycle (concat (repeat 40 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 5 [1.0 0.0]) (repeat 16 [0.5 w])
                   (repeat 40 [1.0 0.0]) (repeat 16 [0.5 w]) (repeat 5 [1.0 0.0]) (repeat 16 [0.5 w])))))

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
              (let [world (assoc (mcl/room) :rng (random/rng (.-seed controls)))
                    states (mcl/mcl-seq world [1.5 4.5 0.0] route (.-particles controls))]
                (reset! model {:world world :rest states :st (first states) :step 0 :clock 0.0 :last-time nil})))
            (advance! []
              (let [now (js/performance.now)
                    {:keys [clock last-time]} @model
                    clock (+ clock (* (.-speed controls) 0.001 (- now (or last-time now))))]
                (swap! model assoc :last-time now)
                ;; One update a frame at most: each costs a scan per particle.
                (if (>= clock 1.0)
                  (swap! model (fn [{:keys [rest step] :as m}]
                                 (assoc m :rest (next rest) :st (fnext rest) :step (inc step) :clock (min 1.0 (- clock 1.0)))))
                  (swap! model assoc :clock clock))))
            (draw! []
              (let [{:keys [world st step]} @model
                    {:keys [truth particles]} st
                    {:keys [walls angles max-range]} world
                    cw (.-width canvas) ch (.-height canvas)
                    scale (min (/ (- cw 40) 10.0) (/ (- ch 120) 6.0))
                    x0 (* 0.5 (- cw (* 10 scale))) y0 50
                    px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 6 y)))])
                    [[ex ey] spread] (mcl/estimate particles)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (set! (.-fillStyle ctx) "#161b25")
                (let [[ax ay] (px [0 6])] (.fillRect ctx ax ay (* 10 scale) (* 6 scale)))
                (set! (.-strokeStyle ctx) "#c9d2e0") (set! (.-lineWidth ctx) 3)
                (doseq [[a b] walls :let [[ax ay] (px a) [bx by] (px b)]]
                  (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
                (set! (.-fillStyle ctx) "rgba(255,160,60,0.5)")
                (doseq [p particles :let [[x y] (px p)]] (.fillRect ctx (- x 1.5) (- y 1.5) 3 3))
                ;; The scan, from where the robot really is.
                (let [[tx ty heading] truth [rx ry] (px truth)]
                  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.25)") (set! (.-lineWidth ctx) 1)
                  (doseq [a angles :let [d (mcl/raycast walls [tx ty] (+ heading a) max-range)
                                         [hx hy] (px [(+ tx (* d (math/cos (+ heading a)))) (+ ty (* d (math/sin (+ heading a))))])]]
                    (.beginPath ctx) (.moveTo ctx rx ry) (.lineTo ctx hx hy) (.stroke ctx))
                  (set! (.-fillStyle ctx) "#ffffff")
                  (.beginPath ctx) (.arc ctx rx ry 7 0 (* 2 math/PI)) (.fill ctx)
                  (set! (.-strokeStyle ctx) "#0d1118") (set! (.-lineWidth ctx) 2)
                  (.beginPath ctx) (.moveTo ctx rx ry) (.lineTo ctx (+ rx (* 10 (math/cos heading))) (- ry (* 10 (math/sin heading)))) (.stroke ctx))
                (let [[sx sy] (px [ex ey])]
                  (set! (.-strokeStyle ctx) "#7dffa0") (set! (.-lineWidth ctx) 2)
                  (.beginPath ctx) (.arc ctx sx sy 10 0 (* 2 math/PI)) (.stroke ctx))
                (set! (.-textContent caption)
                      (str "step " step "  |  estimate off by "
                           (.toFixed (math/sqrt (+ (math/pow (- ex (truth 0)) 2) (math/pow (- ey (truth 1)) 2))) 2)
                           " m, cluster spread " (.toFixed spread 2) " m"
                           "\nwhite: the robot and its scan  orange: particles  green ring: their estimate"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "particles" 100 4000 100) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 30 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "localization")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
