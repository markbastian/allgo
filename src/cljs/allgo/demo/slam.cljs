(ns allgo.demo.slam
  "Canvas demo for `allgo.robotics.localization`'s EKF and
  `allgo.robotics.slam`: the extended Kalman filter of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 26,
  locating a robot among landmarks and then mapping them too.

  The robot (white) drives laps of a field of landmarks (gray crosses),
  steering for the corners of its route by where it believes it is,
  measuring the range and bearing of those within reach (the faint
  lines). Dead reckoning -- the controls alone -- drifts away. EKF
  localization knows the map and holds the robot to within centimeters,
  its uncertainty (the orange ellipse) swelling between sightings and
  shrinking at each. EKF SLAM is not given the map: each landmark is
  placed where it is first seen, with an ellipse of its own, and all of
  them sharpen together when the robot comes round and sees the first
  ones again."
  (:require [allgo.random :as random]
            [allgo.robotics.localization :as loc]
            [allgo.robotics.slam :as slam]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "EKF SLAM" :seed 1 :speed 3})

(def ^:private world (loc/landmark-field))
(def ^:private start [3.0 3.5 0.0])

(def ^:private route
  ;; The corners of the lap, driven to in turn.
  [[15.5 3.5] [15.5 10.0] [3.0 10.0] [3.0 3.5]])

(defn- autopilot
  "The control toward the next corner, steering by where the robot
  believes it is -- so a wrong belief drives it wrong: `[[v w] corner']`."
  [[x y h] corner]
  (let [[tx ty] (route corner)
        corner (if (< (math/hypot (- tx x) (- ty y)) 0.6) (mod (inc corner) (count route)) corner)
        [tx ty] (route corner)
        err (loc/wrap-angle (- (math/atan2 (- ty y) (- tx x)) h))]
    [[1.0 (max -1.5 (min 1.5 (* 2.0 err)))] corner]))

(defn- initial [method]
  {:truth start :k 0 :corner 0 :trail [] :est-trail [] :obs []
   :belief (case method
             "EKF SLAM" (slam/initial-belief start (count (:landmarks world)))
             "EKF localization" [start [[0.01 0.0 0.0] [0.0 0.01 0.0] [0.0 0.0 0.001]]]
             "dead reckoning" start)})

(defn- estimate [method belief]
  (case method
    "EKF SLAM" (slam/pose belief)
    "EKF localization" (first belief)
    "dead reckoning" belief))

(defn- step [method rng {:keys [truth k belief trail est-trail corner]}]
  (let [{:keys [dt sigma-v sigma-w]} world
        [[v w :as u] corner] (autopilot (estimate method belief) corner)
        truth' (loc/move truth v w dt rng sigma-v sigma-w)
        z (loc/observe-landmarks world truth' rng)
        belief' (case method
                  "EKF SLAM" (slam/slam-step world belief u z)
                  "EKF localization" (loc/ekf-localization-step world belief u z)
                  "dead reckoning" (loc/drive belief v w dt))
        est (estimate method belief')
        keep-last (fn [v x] (let [v (conj v x)] (if (> (count v) 1200) (subvec v 1) v)))]
    {:truth truth' :k (inc k) :belief belief' :obs z :corner corner
     :trail (keep-last trail truth') :est-trail (keep-last est-trail est)}))

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
              (reset! model {:rng (random/rng (.-seed controls)) :st (initial (.-method controls))}))
            (advance! []
              (swap! model (fn [{:keys [rng st] :as m}]
                             (assoc m :st (nth (iterate #(step (.-method controls) rng %) st) (.-speed controls))))))
            (ellipse! [px scale [x y] [[a b] [_ d]] color]
              ;; two standard deviations along the covariance's axes
              (let [tr (* 0.5 (+ a d)) disc (math/sqrt (+ (math/pow (* 0.5 (- a d)) 2) (* b b)))
                    l1 (max 0.0 (+ tr disc)) l2 (max 0.0 (- tr disc))
                    angle (* 0.5 (math/atan2 (* 2 b) (- a d)))
                    [cx cy] (px [x y])]
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 1.5)
                (.beginPath ctx)
                (.ellipse ctx cx cy (max 1 (* scale 2 (math/sqrt l1))) (max 1 (* scale 2 (math/sqrt l2))) (- angle) 0 (* 2 math/PI))
                (.stroke ctx)))
            (draw! []
              (let [{{:keys [truth belief trail est-trail obs k]} :st} @model
                    method (.-method controls)
                    cw (.-width canvas) ch (.-height canvas)
                    scale (min (/ (- cw 40) 20.0) (/ (- ch 110) 14.0))
                    x0 (* 0.5 (- cw (* 20 scale))) y0 40
                    px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 14 y)))])
                    est (peek est-trail)
                    err (when est (math/sqrt (+ (math/pow (- (est 0) (truth 0)) 2) (math/pow (- (est 1) (truth 1)) 2))))]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-fillStyle ctx) "#161b25")
                (let [[ax ay] (px [0 14])] (.fillRect ctx ax ay (* 20 scale) (* 14 scale)))
                ;; the true landmarks
                (set! (.-strokeStyle ctx) "#6b7894") (set! (.-lineWidth ctx) 2)
                (doseq [m (:landmarks world) :let [[x y] (px m)]]
                  (.beginPath ctx) (.moveTo ctx (- x 5) (- y 5)) (.lineTo ctx (+ x 5) (+ y 5))
                  (.moveTo ctx (- x 5) (+ y 5)) (.lineTo ctx (+ x 5) (- y 5)) (.stroke ctx))
                ;; what the sensor sees
                (let [[tx ty] (px truth)]
                  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.18)") (set! (.-lineWidth ctx) 1)
                  (doseq [[id] obs :let [[mx my] (px ((:landmarks world) id))]]
                    (.beginPath ctx) (.moveTo ctx tx ty) (.lineTo ctx mx my) (.stroke ctx)))
                ;; trails
                (doseq [[pts color] [[trail "rgba(255,255,255,0.5)"] [est-trail "rgba(255,170,80,0.9)"]]
                        :when (> (count pts) 1)]
                  (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 1.5)
                  (.beginPath ctx)
                  (let [[x y] (px (first pts))] (.moveTo ctx x y))
                  (doseq [p (rest pts) :let [[x y] (px p)]] (.lineTo ctx x y))
                  (.stroke ctx))
                ;; the estimate's uncertainty
                (case method
                  "EKF SLAM" (do (ellipse! px scale (slam/pose belief) (mapv #(subvec % 0 2) (subvec (:cov belief) 0 2)) "#ffaa50")
                                 (doseq [[id m] (slam/landmarks belief)]
                                   (ellipse! px scale m (slam/landmark-cov belief id) "#7dc4ff")
                                   (let [[x y] (px m)] (set! (.-fillStyle ctx) "#7dc4ff") (.fillRect ctx (- x 2) (- y 2) 4 4))))
                  "EKF localization" (let [[m P] belief] (ellipse! px scale m (mapv #(subvec % 0 2) (subvec P 0 2)) "#ffaa50"))
                  nil)
                ;; the robot
                (let [[x y] (px truth) h (truth 2)]
                  (set! (.-fillStyle ctx) "#ffffff")
                  (.beginPath ctx) (.arc ctx x y 6 0 (* 2 math/PI)) (.fill ctx)
                  (set! (.-strokeStyle ctx) "#0d1118") (set! (.-lineWidth ctx) 2)
                  (.beginPath ctx) (.moveTo ctx x y) (.lineTo ctx (+ x (* 9 (math/cos h))) (- y (* 9 (math/sin h)))) (.stroke ctx))
                (set! (.-textContent caption)
                      (str method "  |  step " k (when err (str "  |  estimate off by " (.toFixed err 2) " m"))
                           (when (= method "EKF SLAM") (str "  |  " (count (:seen belief)) " of " (count (:landmarks world)) " landmarks mapped"))
                           "\nwhite: the robot and its true path  orange: the estimate and its 2-sigma ellipse"
                           (when (= method "EKF SLAM") "  blue: the mapped landmarks")))))
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
        (.onChange (.add gui controls "method" #js ["EKF SLAM" "EKF localization" "dead reckoning"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 10 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "slam")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
