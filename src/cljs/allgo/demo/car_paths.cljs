(ns allgo.demo.car-paths
  "Canvas demo for `allgo.robotics.dubins`: the shortest paths for a car,
  after LaValle, *Planning Algorithms*, Cambridge University Press, 2006,
  section 15.3.

  The car starts at the green pose; press anywhere to put the goal there
  and drag to aim it. Dubins' car only drives forward and turns no
  tighter than its radius, and its shortest path is always one of six
  words of turns and straights -- every word that joins the two poses is
  drawn faintly, the shortest in gold. Reeds and Shepp's car may also
  reverse: backward stretches are drawn in coral, and a goal just behind
  the car is a short reverse away rather than a loop."
  (:require [allgo.robotics.dubins :as dubins]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:car "Reeds-Shepp" :radius 1.0 :speed 4 :candidates true})

(def ^:private start [4.0 4.0 0.0])

(defn- paths [car goal rho]
  (if (= car "Dubins")
    (dubins/dubins-paths start goal rho)
    (dubins/reeds-shepp-paths start goal rho)))

(defn- pieces
  "The path's segments as `[backward? poses]`, eight poses a segment."
  [{:keys [word lengths]} rho]
  (:out (reduce (fn [{:keys [pose out]} [kind len]]
                  (let [poses (mapv #(dubins/segment pose kind (* rho len (/ % 8.0)) rho) (range 9))]
                    {:pose (peek poses) :out (conj out [(neg? len) poses])}))
                {:pose start :out []}
                (map vector word lengths))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        state (atom {:goal [8.0 5.5 (/ math/PI 2)] :t 0.0 :dragging false})]
    (set! (.-style canvas) "display:block;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(layout []
              (let [cw (.-width canvas) ch (.-height canvas)
                    scale (min (/ (- cw 40) 12.0) (/ (- ch 110) 8.0))
                    x0 (* 0.5 (- cw (* 12 scale))) y0 40]
                {:scale scale :x0 x0 :y0 y0
                 :px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 8 y)))])
                 :world (fn [cx cy] [(/ (- cx x0) scale) (- 8 (/ (- cy y0) scale))])}))
            (polyline! [px pts color width]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
              (.beginPath ctx)
              (let [[x y] (px (first pts))] (.moveTo ctx x y))
              (doseq [p (rest pts) :let [[x y] (px p)]] (.lineTo ctx x y))
              (.stroke ctx))
            (car! [px scale [x y th] color]
              ;; a car a metre long, its nose a brighter triangle
              (let [corner (fn [f s] [(+ x (* f (math/cos th)) (* s (- (math/sin th)))) (+ y (* f (math/sin th)) (* s (math/cos th)))])
                    body [(corner 0.6 0.3) (corner -0.4 0.3) (corner -0.4 -0.3) (corner 0.6 -0.3)]]
                (set! (.-fillStyle ctx) color)
                (.beginPath ctx)
                (let [[a b] (px (first body))] (.moveTo ctx a b))
                (doseq [p (rest body) :let [[a b] (px p)]] (.lineTo ctx a b))
                (.closePath ctx) (.fill ctx)
                (let [[a b] (px (corner 0.6 0.0))]
                  (set! (.-fillStyle ctx) "#0d1118")
                  (.beginPath ctx) (.arc ctx a b (* 0.08 scale) 0 (* 2 math/PI)) (.fill ctx))))
            (draw! []
              (let [{:keys [goal t]} @state
                    {:keys [px scale]} (layout)
                    rho (.-radius controls)
                    car (.-car controls)
                    ps (paths car goal rho)
                    best (first ps)
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (let [[ax ay] (px [0 8]) [bx by] (px [12 0])]
                  (set! (.-fillStyle ctx) "#161b25")
                  (.fillRect ctx ax ay (- bx ax) (- by ay)))
                (when (.-candidates controls)
                  (doseq [p (rest ps) [_ poses] (pieces p rho)]
                    (polyline! px poses "rgba(140,170,220,0.3)" 1)))
                (when best
                  (let [segs (pieces best rho)]
                    (doseq [[back? poses] segs]
                      (polyline! px poses (if back? "#ff8f70" "#ffd479") 3))
                    ;; the car driven along it, round and round
                    (let [all (vec (mapcat (fn [[_ poses]] (rest poses)) segs))
                          k (mod (long t) (max 1 (count all)))]
                      (car! px scale (if (seq all) (all k) start) "#ffffff"))))
                (car! px scale start "rgba(125,255,160,0.8)")
                (car! px scale goal "rgba(255,143,112,0.8)")
                (set! (.-textContent caption)
                      (str car (if (= car "Dubins") " car, forward only" " car, forward and back")
                           "  |  "
                           (if best
                             (str "shortest " (:word best) ", " (.toFixed (:length best) 2) " m"
                                  (when (next ps)
                                    (str "; then " (apply str (interpose ", " (map #(str (:word %) " " (.toFixed (:length %) 2)) (take 4 (rest ps))))))))
                             "no path")
                           "\ngreen: start  red: goal  gold: forward  coral: reverse  faint: the other words  press and drag: place and aim the goal"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (swap! state update :t + (* 0.25 (.-speed controls)))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))
            (pointer [^js e]
              (let [r (.getBoundingClientRect canvas)]
                ((:world (layout)) (- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r)))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (let [[x y] (pointer e)]
                             (swap! state (fn [s] (assoc s :dragging true :t 0.0 :goal [x y (get-in s [:goal 2])]))))))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (when (:dragging @state)
                             (let [[x y] (pointer e)
                                   [gx gy] (:goal @state)]
                               (when (> (math/hypot (- x gx) (- y gy)) 0.2)
                                 (swap! state assoc-in [:goal 2] (math/atan2 (- y gy) (- x gx))))))))
      (doseq [ev ["pointerup" "pointercancel" "pointerleave"]]
        (.addEventListener canvas ev (fn [_] (swap! state assoc :dragging false))))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "car" #js ["Dubins" "Reeds-Shepp"]) #(swap! state assoc :t 0.0))
        (.add gui controls "radius" 0.5 2.0 0.1)
        (.add gui controls "candidates")
        (.add gui controls "speed" 1 20 1))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "car-paths")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
