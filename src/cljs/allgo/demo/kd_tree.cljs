(ns allgo.demo.kd-tree
  "Canvas demo for `allgo.spatial.kdtree`: finding the nearest of many
  points without measuring them all, after LaValle, *Planning
  Algorithms*, Cambridge University Press, 2006, section 5.5.2.

  The points are split at the median of x, each half at the median of y,
  and so on -- the lines drawn are the splits, fainter as they go deeper.
  Move the pointer: the query descends to the cell it lies in, then backs
  out, crossing a split only when the splitting line is nearer than the
  best point found so far. The nodes it visited are lit; the count says
  how few of all the points it had to measure. Nearest, the k nearest,
  or every point within a radius."
  (:require [allgo.random :as random]
            [allgo.spatial.kdtree :as kdtree]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:query "nearest" :points 2000 :k 8 :radius 0.08 :seed 1 :layout "clusters"})

(defn- make-points [n layout rng]
  (vec (case layout
         "uniform" (repeatedly n (fn [] [(rng) (rng)]))
         "clusters" (let [centers (vec (repeatedly 7 (fn [] [(+ 0.1 (* 0.8 (rng))) (+ 0.1 (* 0.8 (rng)))])))]
                      (repeatedly n (fn [] (let [[cx cy] (random/pick rng centers)]
                                             [(min 1.0 (max 0.0 (random/gaussian rng cx 0.07)))
                                              (min 1.0 (max 0.0 (random/gaussian rng cy 0.07)))])))))))

(defn- splits
  "Every split of `tree` as `[axis value lo hi depth]` over the box it
  divides, `[[x0 y0] [x1 y1]]`."
  [tree box depth]
  (when-let [{:keys [point axis left right]} tree]
    (let [[[x0 y0] [x1 y1]] box v (nth point axis)]
      (concat [[axis v box depth]]
              (splits left (if (zero? axis) [[x0 y0] [v y1]] [[x0 y0] [x1 v]]) (inc depth))
              (splits right (if (zero? axis) [[v y0] [x1 y1]] [[x0 v] [x1 y1]]) (inc depth))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        pointer (atom [0.5 0.5])]
    (set! (.-style canvas) "display:block;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [pts (make-points (.-points controls) (.-layout controls) (random/rng (.-seed controls)))
                    tree (kdtree/build pts)]
                (reset! model {:points pts :tree tree :splits (vec (splits tree [[0.0 0.0] [1.0 1.0]] 0))})))
            (layout []
              (let [cw (.-width canvas) ch (.-height canvas)
                    side (min (- cw 80) (- ch 120))
                    x0 (* 0.5 (- cw side)) y0 40]
                {:side side :x0 x0 :y0 y0
                 :px (fn [[u v]] [(+ x0 (* side u)) (+ y0 (* side (- 1 v)))])
                 :uv (fn [x y] [(/ (- x x0) side) (- 1 (/ (- y y0) side))])}))
            (draw! []
              (let [{:keys [points tree splits]} @model
                    {:keys [side x0 y0 px]} (layout)
                    cw (.-width canvas) ch (.-height canvas)
                    q @pointer
                    query (.-query controls)
                    result (case query
                             "nearest" (kdtree/nearest tree q)
                             "k nearest" (kdtree/k-nearest tree q (.-k controls))
                             "within radius" (kdtree/within tree q (.-radius controls)))
                    found (if (= query "nearest") [(first result)] result)
                    seen (set (:measured (meta result)))]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x0 y0 side side)
                (doseq [[axis v [[bx0 by0] [bx1 by1]] depth] splits :when (< depth 9)]
                  (let [[a b] (if (zero? axis) [[v by0] [v by1]] [[bx0 v] [bx1 v]])
                        [ax ay] (px a) [ex ey] (px b)]
                    (set! (.-strokeStyle ctx) (str "rgba(140,170,220," (max 0.08 (- 0.7 (* 0.07 depth))) ")"))
                    (set! (.-lineWidth ctx) (max 0.5 (- 2.5 (* 0.25 depth))))
                    (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx ex ey) (.stroke ctx)))
                (doseq [p points :let [[x y] (px p)]]
                  (set! (.-fillStyle ctx) (if (seen p) "#ffd479" "rgba(201,210,224,0.45)"))
                  (.fillRect ctx (- x 1.5) (- y 1.5) 3 3))
                (let [[qx qy] (px q)]
                  (doseq [p found :let [[x y] (px p)]]
                    (set! (.-strokeStyle ctx) "#7dffa0") (set! (.-lineWidth ctx) 1.5)
                    (.beginPath ctx) (.moveTo ctx qx qy) (.lineTo ctx x y) (.stroke ctx)
                    (set! (.-fillStyle ctx) "#7dffa0") (.beginPath ctx) (.arc ctx x y 4 0 (* 2 math/PI)) (.fill ctx))
                  (when (= query "within radius")
                    (set! (.-strokeStyle ctx) "rgba(125,255,160,0.6)") (set! (.-lineWidth ctx) 1)
                    (.beginPath ctx) (.arc ctx qx qy (* side (.-radius controls)) 0 (* 2 math/PI)) (.stroke ctx))
                  (set! (.-fillStyle ctx) "#ffffff") (.beginPath ctx) (.arc ctx qx qy 5 0 (* 2 math/PI)) (.fill ctx))
                (set! (.-textContent caption)
                      (str query "  |  " (count found) " found, " (count seen) " of " (count points) " points measured ("
                           (.toFixed (* 100 (/ (count seen) (count points))) 1) "%)"
                           "\nlines: the tree's splits, fainter deeper  gold: the points the search measured  green: what it found  move the pointer: the query"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (let [r (.getBoundingClientRect canvas)
                                 [u v] ((:uv (layout)) (- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r)))]
                             (reset! pointer [(max 0.0 (min 1.0 u)) (max 0.0 (min 1.0 v))]))))
      (let [gui (GUI. #js {:container container})]
        (.add gui controls "query" #js ["nearest" "k nearest" "within radius"])
        (.add gui controls "k" 1 40 1)
        (.add gui controls "radius" 0.01 0.3 0.01)
        (.onChange (.add gui controls "points" 100 20000 100) restart!)
        (.onChange (.add gui controls "layout" #js ["clusters" "uniform"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "kd-tree")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
