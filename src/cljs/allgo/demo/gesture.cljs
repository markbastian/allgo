(ns allgo.demo.gesture
  "Canvas demo for `allgo.procedural.gesture`: a figure with no
  keyframes, moved by sines and noise.

  Pick an action for each layer -- the body, the arms over it, the head
  over those -- and the figure fades into it over `transition` seconds.
  The strip along the bottom traces one angle from each layer through
  time, with a mark at every change: the curves bend into their new
  rhythm rather than jumping to it, even when an action is changed
  again halfway through a transition. Leave it at `rest` or `fidget`
  for a while to see the point of noise: the figure never does quite
  the same thing twice, and never stops.

  Drag to walk around it."
  (:require [allgo.procedural.gesture :as gesture]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:body       "fidget"
       :arms       "none"
       :head       "none"
       :transition 1.0
       :speed      1.0
       :seed       0
       :yaw        25})

(def ^:private traces
  "One angle from each layer, for the strip: `[label joint axis color]`."
  [["chest pitch (body)" :chest 0 "#7fb2e5"]
   ["right shoulder (arms)" :r-shoulder 2 "#f2c14e"]
   ["head turn (head)" :head 1 "#e57f9a"]])

(def ^:private history-seconds 8.0)

(defn- projector
  "World to screen: turned by `yaw` about the vertical, seen a little
  from above, in perspective."
  [w h yaw]
  (let [cy (math/cos yaw) sy (math/sin yaw)
        pitch 0.18 cp (math/cos pitch) sp (math/sin pitch)
        scale (* 0.3 (min h (* 1.6 w)))
        cx (* 0.5 w) base (* 0.78 h)]
    (fn [[x y z]]
      (let [x1 (+ (* cy x) (* sy z)) z1 (- (* cy z) (* sy x))
            y (- y 0.9)
            y2 (- (* cp y) (* sp z1)) z2 (+ (* sp y) (* cp z1))
            k (/ 3.2 (- 3.2 z2))]
        [(+ cx (* scale k x1)) (- base (* scale k (+ y2 0.9))) z2 k]))))

(def ^:private widths
  {:pelvis 0.13 :chest 0.12 :neck 0.06 :l-hip 0.07 :r-hip 0.07})

(defn- draw-figure! [^js ctx project points]
  (let [scale-of (fn [[_ _ _ k]] k)
        segs (for [[a b] gesture/bones]
               (let [pa (project (points a)) pb (project (points b))]
                 {:a pa :b pb :depth (* 0.5 (+ (nth pa 2) (nth pb 2)))
                  :width (get widths a 0.055)
                  :color (cond (#{:pelvis :chest :neck} a) "#d9a46c"
                               (= "l" (subs (name a) 0 1)) "#e9c49a"
                               :else "#b98552")}))
        [hx hy hz hk] (project (mapv #(* 0.5 (+ %1 %2)) (points :head) (points :crown)))
        head {:head true :depth hz}
        unit (/ (- (first (project [1.0 0.9 0.0])) (first (project [0.0 0.9 0.0]))) 1.0)]
    ;; Shadows first: every bone flattened onto the floor.
    (set! (.-strokeStyle ctx) "rgba(0,0,0,0.22)")
    (set! (.-lineCap ctx) "round")
    (doseq [[a b] gesture/bones]
      (let [[ax _ az] (points a) [bx _ bz] (points b)
            [sx sy _ k] (project [ax 0.0 az]) [ex ey] (project [bx 0.0 bz])]
        (set! (.-lineWidth ctx) (* unit k 0.08))
        (.beginPath ctx) (.moveTo ctx sx sy) (.lineTo ctx ex ey) (.stroke ctx)))
    ;; Back to front.
    (doseq [s (sort-by :depth (cons head segs))]
      (if (:head s)
        (do (set! (.-fillStyle ctx) "#e2b07a")
            (.beginPath ctx) (.arc ctx hx hy (* unit hk 0.11) 0 (* 2 math/PI)) (.fill ctx))
        (let [{[ax ay] :a [bx by] :b :keys [width color]} s]
          (set! (.-strokeStyle ctx) color)
          (set! (.-lineWidth ctx) (* unit width (* 0.5 (+ (scale-of (:a s)) (scale-of (:b s))))))
          (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))))))

(defn- draw-floor! [^js ctx project]
  (set! (.-strokeStyle ctx) "rgba(200,210,225,0.12)")
  (set! (.-lineWidth ctx) 1)
  (doseq [i (range -6 7)]
    (let [a (* 0.25 i)
          [x0 y0] (project [a 0.0 -1.5]) [x1 y1] (project [a 0.0 1.5])
          [x2 y2] (project [-1.5 0.0 a]) [x3 y3] (project [1.5 0.0 a])]
      (.beginPath ctx) (.moveTo ctx x0 y0) (.lineTo ctx x1 y1)
      (.moveTo ctx x2 y2) (.lineTo ctx x3 y3) (.stroke ctx))))

(defn- draw-strip! [^js ctx w h history changes now]
  (let [top (- h 96) height 70 left 16 right (- w 16)
        x-of (fn [t] (+ left (* (- right left) (/ (- t (- now history-seconds)) history-seconds))))]
    (set! (.-fillStyle ctx) "rgba(10,14,22,0.55)")
    (.fillRect ctx (- left 6) (- top 18) (+ 12 (- right left)) (+ height 30))
    (set! (.-strokeStyle ctx) "rgba(255,255,255,0.35)")
    (set! (.-lineWidth ctx) 1)
    (doseq [[t _] changes :when (> t (- now history-seconds))]
      (.beginPath ctx) (.moveTo ctx (x-of t) top) (.lineTo ctx (x-of t) (+ top height)) (.stroke ctx))
    (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
    (doseq [[n [label joint axis color]] (map-indexed vector traces)]
      (set! (.-fillStyle ctx) color)
      (.fillText ctx label (+ left (* n 190)) (- top 5))
      (set! (.-strokeStyle ctx) color)
      (set! (.-lineWidth ctx) 1.5)
      (.beginPath ctx)
      (doseq [[i [t pose]] (map-indexed vector history)]
        (let [v (max (/ math/PI -3.0) (min (/ math/PI 3.0) (nth (get pose joint [0.0 0.0 0.0]) axis)))
              ;; Plus or minus 60 degrees across the strip, and held to it.
              y (+ top (* 0.5 height) (* -0.5 height (/ v (/ math/PI 3.0))))]
          (if (zero? i) (.moveTo ctx (x-of t) y) (.lineTo ctx (x-of t) y))))
      (.stroke ctx))
    (when-let [[_ text] (last changes)]
      (set! (.-fillStyle ctx) "rgba(255,255,255,0.7)")
      (.fillText ctx text (- right 220) (- top 5)))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        running? (atom false)
        clock (atom 0.0)
        last-frame (atom nil)
        performance (atom (gesture/start 0.0 {:body (.-body controls)}))
        history (atom [])
        changes (atom [])]
    (set! (.-style canvas) "display:block;cursor:grab;touch-action:none")
    (.appendChild container canvas)
    (letfn [(resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h))))
            (change! [layer]
              (let [action (aget controls (name layer))]
                (swap! performance gesture/change layer action @clock (.-transition controls))
                (swap! changes conj [@clock (str (name layer) " -> " action)])))
            (frame [now]
              (when @running?
                (js/requestAnimationFrame frame)
                (let [dt (if-let [t0 @last-frame] (min 0.1 (/ (- now t0) 1000.0)) 0.0)
                      t (swap! clock + (* dt (.-speed controls)))
                      pose (gesture/pose @performance t {:transition (.-transition controls)
                                                         :seed (.-seed controls)})
                      points (gesture/figure pose)
                      w (.-width canvas) h (.-height canvas)
                      project (projector w h (* (.-yaw controls) (/ math/PI 180)))]
                  (reset! last-frame now)
                  (swap! history (fn [hs] (into [] (drop-while #(< (first %) (- t history-seconds 0.5)))
                                                (conj hs [t pose]))))
                  (let [g (.createLinearGradient ctx 0 0 0 h)]
                    (.addColorStop g 0 "#27324a") (.addColorStop g 1 "#11151f")
                    (set! (.-fillStyle ctx) g)
                    (.fillRect ctx 0 0 w h))
                  (draw-floor! ctx project)
                  (draw-figure! ctx project points)
                  (draw-strip! ctx w h @history @changes t))))]
      (let [gui (GUI. #js {:container container})
            names (fn [layer] (clj->js (sort (keys (gesture/actions layer)))))]
        (-> (.add gui controls "body" (names :body)) (.onChange #(change! :body)))
        (-> (.add gui controls "arms" (names :arms)) (.onChange #(change! :arms)))
        (-> (.add gui controls "head" (names :head)) (.onChange #(change! :head)))
        (.add gui controls "transition" 0 3 0.1)
        (.add gui controls "speed" 0.25 2 0.05)
        (.add gui controls "seed" 0 99 1)
        (.listen (.add gui controls "yaw" -180 180 1)))
      (let [drag (atom nil)]
        (.addEventListener canvas "pointerdown"
                           (fn [^js e]
                             (try (.setPointerCapture canvas (.-pointerId e)) (catch :default _ nil))
                             (reset! drag (.-clientX e))))
        (.addEventListener canvas "pointermove"
                           (fn [^js e]
                             (when-let [x @drag]
                               (set! (.-yaw controls)
                                     (-> (.-yaw controls) (+ (* 0.4 (- (.-clientX e) x)))
                                         (+ 540) (mod 360) (- 180)))
                               (reset! drag (.-clientX e)))))
        (doseq [ev ["pointerup" "pointercancel"]]
          (.addEventListener canvas ev (fn [_] (reset! drag nil)))))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running?
                       (reset! running? true)
                       (reset! last-frame nil)
                       (resize!)
                       (js/requestAnimationFrame frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "gesture")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
