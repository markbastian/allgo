(ns allgo.demo.algebraic-decomposition
  "Canvas demo for `allgo.robotics.cylindrical`: the cylindrical algebraic
  decomposition of the plane by a few polynomials, after LaValle,
  *Planning Algorithms*, Cambridge University Press, 2006, section 6.4.2.

  First the curves, the polynomials' zeros. Then the projection: on the
  axis below, the critical values, each where a curve folds over in y
  (gold), two curves cross (red) or a curve runs off to infinity or is a
  vertical line (blue), with a vertical line up through each. Then the
  cells, a cylinder at a time from the left: over each open interval the
  curves cut the cylinder into sectors, shaded by whether the planning
  problem's free space holds them; then a sample point in each, and the
  roadmap joining the free 2-cells through the vertical 1-cells beside
  them; last a route between two points, kept inside the cells it
  crosses. The gingerbread face is Example 3.1's, decomposed as in Figure
  6.35 into 37 2-cells, 64 1-cells and 28 0-cells. Hover a cell for its
  signs; click to move an end of the route."
  (:require [allgo.numerics.multivariate :as mv]
            [allgo.robotics.cylindrical :as cad]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:scene "gingerbread face" :free "outside the obstacles" :numbers true :speed 3})

(def ^:private x (mv/variable 2 0))
(def ^:private y (mv/variable 2 1))
(defn- c [k] (mv/constant 2 k))

(defn- circle [cx cy r2]
  (mv/add (mv/power (mv/sub x (c cx)) 2) (mv/power (mv/sub y (c cy)) 2) (c (- r2))))

(def ^:private scenes
  {"gingerbread face"
   {:polys [(circle 0 0 25)
            (mv/scale (circle -2 2 1) -1.0)
            (mv/scale (circle 2 2 1) -1.0)
            (mv/scale (mv/add (mv/power x 2) (mv/scale (mv/power (mv/add y (c 2)) 2) 4.0) (c -4)) -1.0)]
    :names ["head" "left eye" "right eye" "mouth"]
    ;; Example 3.1: the face is where every primitive is at most zero
    :obstacle (fn [signs] (every? #(<= % 0) signs))
    :interior (fn [signs] (every? neg? signs))
    :view [[-7.0 -6.5] [7.0 6.5]]
    :ends [[-6.2 -5.2] [5.8 4.8]]
    :inner-ends [[-3.5 0.0] [3.5 -0.5]]
    :blurb "Example 3.1 and Figure 6.35"}
   "crossing discs"
   {:polys [(circle -2 0 9) (circle 2 0 9) (circle 0 3 9)]
    :names ["left" "right" "top"]
    :obstacle (fn [signs] (some #(<= % 0) signs))
    :interior (fn [signs] (some neg? signs))
    :view [[-7.0 -4.5] [7.0 7.5]]
    :ends [[-6.0 -3.8] [6.0 6.8]]
    :inner-ends [[-4.0 0.0] [3.5 0.5]]
    :blurb "three discs: folds, and crossings where their circles meet"}
   "polygons"
   {:polys [(mv/sub x (c 1)) (mv/add x (c 1)) (mv/sub y (c 1)) (mv/add y (c 1))
            (mv/sub (mv/add x y) (c 4)) (mv/sub x (c 2)) y]
    :names ["x - 1" "x + 1" "y - 1" "y + 1" "x + y - 4" "x - 2" "y"]
    :obstacle (fn [[a b cc dd e f g]] (or (and (<= a 0) (>= b 0) (<= cc 0) (>= dd 0))
                                          (and (<= e 0) (>= f 0) (>= g 0))))
    :interior (fn [[a b cc dd e f g]] (or (and (neg? a) (pos? b) (neg? cc) (pos? dd))
                                          (and (neg? e) (pos? f) (pos? g))))
    :view [[-3.5 -2.5] [6.5 4.5]]
    :ends [[-3.0 1.8] [5.5 0.6]]
    :inner-ends [[-0.5 0.5] [0.5 -0.5]]
    :blurb "a square and a triangle as half-planes: Figure 6.19's decomposition, made algebraic"}})

(def ^:private why-colors {:fold "#ffd479" :crossing "#ff7a7a" :leading "#7fb2ff"})

(defn- spaced
  "n + 1 values from a to b, bunched toward the ends, where the curves fold."
  [a b n]
  (mapv #(+ a (* (- b a) 0.5 (- 1 (math/cos (* math/PI (/ % n)))))) (range (inc n))))

(defn- geometry
  "What to draw of decomposition `d` inside `view`: each interval
  column's sections as polylines, and each of its sectors as a polygon."
  [d [[vx0 vy0] [vx1 vy1]]]
  (into {}
        (for [{:keys [index kind lo hi roots]} (:columns d)
              :when (= kind :interval)
              :let [a (max (or lo vx0) vx0) b (min (or hi vx1) vx1)]
              :when (< a b)
              :let [w (- b a)
                    xs (spaced (if (= a lo) (+ a (* 1e-7 w)) a) (if (= b hi) (- b (* 1e-7 w)) b) 48)
                    n (count roots)
                    columns (keep (fn [x] (let [ys (cad/sections d x)] (when (= n (count ys)) [x ys]))) xs)
                    clamp #(max (- vy0 1.0) (min (+ vy1 1.0) %))
                    curve (fn [k] (mapv (fn [[x ys]] [x (clamp (ys k))]) columns))
                    edge (fn [k fallback] (if (< -1 k n) (curve k) (mapv (fn [[x]] [x fallback]) columns)))]]
          [index {:curves (mapv curve (range n))
                  :sectors (vec (for [s (range (inc n))]
                                  (vec (concat (edge (dec s) (- vy0 1.0)) (rseq (edge s (+ vy1 1.0)))))))}])))

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
    (letfn [(scene [] (scenes (.-scene controls)))
            (free? []
              (let [{:keys [obstacle interior]} (scene)]
                (if (= "outside the obstacles" (.-free controls)) (complement obstacle) interior)))
            (plan! [ends]
              (let [{:keys [rm]} @model]
                (swap! model assoc :ends ends :route (cad/route rm (ends 0) (ends 1)))))
            (restart! []
              (let [{:keys [polys view ends inner-ends]} (scene)
                    d (cad/decompose polys)
                    rm (cad/roadmap d (free?))
                    comps (cad/components rm)]
                (reset! model {:d d :rm rm :geometry (geometry d view)
                               :component (into {} (for [[k comp] (map-indexed vector comps) id comp] [id k]))
                               :components (count comps)
                               :links (vec (for [[m ns] (:graph rm) n (keys ns) :when (neg? (compare m n))]
                                             (cad/link-path rm m n)))
                               :t 0.0 :hover nil :next 0})
                (plan! (if (= "outside the obstacles" (.-free controls)) ends inner-ends))))
            (layout []
              (let [cw (.-width canvas) ch (.-height canvas)
                    [[vx0 vy0] [vx1 vy1]] (:view (scene))
                    k (min (/ (- cw 60) (- vx1 vx0)) (/ (- ch 190) (- vy1 vy0)))
                    ox (* 0.5 (- cw (* k (- vx1 vx0))))
                    oy 50]
                {:k k
                 :px (fn [[px py]] [(+ ox (* k (- px vx0))) (+ oy (* k (- vy1 py)))])
                 :world (fn [sx sy] [(+ vx0 (/ (- sx ox) k)) (- vy1 (/ (- sy oy) k))])
                 :axis-y (+ oy (* k (- vy1 vy0)) 40)}))
            (polyline! [px pts color width]
              (when (seq pts)
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                (.beginPath ctx)
                (doseq [[i p] (map-indexed vector pts) :let [[sx sy] (px p)]]
                  (if (zero? i) (.moveTo ctx sx sy) (.lineTo ctx sx sy)))
                (.stroke ctx)))
            (polygon! [px pts color]
              (set! (.-fillStyle ctx) color)
              (.beginPath ctx)
              (doseq [[i p] (map-indexed vector pts) :let [[sx sy] (px p)]]
                (if (zero? i) (.moveTo ctx sx sy) (.lineTo ctx sx sy)))
              (.closePath ctx)
              (.fill ctx))
            (dot! [px p r color]
              (let [[sx sy] (px p)]
                (set! (.-fillStyle ctx) color) (.beginPath ctx) (.arc ctx sx sy r 0 (* 2 math/PI)) (.fill ctx)))
            (hue [k] (str "hsla(" (mod (+ 150 (* 67 k)) 360) ",55%,45%,0.40)"))
            (draw! []
              (let [{:keys [d rm geometry component components links route ends t hover]} @model
                    {:keys [px axis-y]} (layout)
                    {[[vx0 vy0] [vx1 vy1]] :view :keys [blurb names]} (scene)
                    cw (.-width canvas) ch (.-height canvas)
                    columns (:columns d)
                    ncol (count columns)
                    stage (fn [s] (max 0.0 (min 1.0 (- t s))))
                    shown (long (math/floor (* (stage 2) ncol)))
                    free (:free rm)
                    in-view (fn [[sx sy]] (and (<= vx0 sx vx1) (<= vy0 sy vy1)))]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (.save ctx)
                (let [[ax ay] (px [vx0 vy1]) [bx by] (px [vx1 vy0])]
                  (.beginPath ctx) (.rect ctx ax ay (- bx ax) (- by ay)) (.clip ctx)
                  (set! (.-fillStyle ctx) "#141a24") (.fillRect ctx ax ay (- bx ax) (- by ay)))
                ;; the cells, a cylinder at a time
                (doseq [{:keys [index]} (take shown columns)
                        :when (even? index)
                        :let [{:keys [sectors]} (geometry index)]
                        [s poly] (map-indexed vector sectors)
                        :let [id [index (* 2 s)]]]
                  (polygon! px poly (cond (= id hover) "rgba(255,255,255,0.22)"
                                          (free id) (hue (component id 0))
                                          :else "rgba(70,80,100,0.35)")))
                ;; the vertical lines over the critical values
                (doseq [[i {:keys [x why]}] (map-indexed vector (:critical d))
                        :let [col (inc (* 2 i))
                              f (stage 1)]
                        :when (pos? f)]
                  (let [{:keys [roots]} (columns col)
                        ys (concat [(- vy0 1.0)] (map :y roots) [(+ vy1 1.0)])]
                    (doseq [[j [y0 y1]] (map-indexed vector (partition 2 1 ys))
                            :let [id [col (* 2 j)]]]
                      (polyline! px [[x y0] [x (+ y0 (* f (- y1 y0)))]]
                                 (cond (= id hover) "#ffffff"
                                       (and (>= shown col) (free id)) "rgba(220,230,255,0.75)"
                                       :else (str (why-colors (:kind (first why))) "66"))
                                 (if (= id hover) 2.5 1)))))
                ;; the curves
                (doseq [[_ {:keys [curves]}] geometry c curves]
                  (polyline! px c "#e6ecf5" 1.6))
                ;; samples and the roadmap
                (when (pos? (stage 3))
                  (doseq [{:keys [cells]} columns {:keys [sample dim]} cells
                          :when (in-view sample)]
                    (dot! px sample (case dim 2 2.6 1 1.8 1.2) (case dim 2 "#ffffff" 1 "#9fb4d8" "#56627a"))))
                (when (pos? (stage 4))
                  (doseq [l (take (long (* (stage 4) (count links))) links)]
                    (polyline! px l "rgba(255,255,255,0.45)" 1)))
                (when (and (.-numbers controls) (pos? (stage 3)))
                  (set! (.-font ctx) "10px ui-monospace,Menlo,monospace")
                  (set! (.-fillStyle ctx) "#c9d2e0")
                  (doseq [[n {:keys [sample]}] (map-indexed vector (filter #(= 2 (:dim %)) (cad/cells d)))
                          :when (in-view sample)
                          :let [[sx sy] (px sample)]]
                    (.fillText ctx (str (inc n)) (+ sx 4) (- sy 4))))
                (when (pos? (stage 5))
                  (when route (polyline! px route "#7dffa0" 3))
                  (dot! px (ends 0) 5 "#7dffa0")
                  (dot! px (ends 1) 5 "#ff8f70"))
                (.restore ctx)
                ;; the line below: its decomposition by the projection
                (let [[ax] (px [vx0 0]) [bx] (px [vx1 0])]
                  (set! (.-strokeStyle ctx) "#5b6884") (set! (.-lineWidth ctx) 1)
                  (.beginPath ctx) (.moveTo ctx ax axis-y) (.lineTo ctx bx axis-y) (.stroke ctx)
                  (set! (.-font ctx) "10px ui-monospace,Menlo,monospace")
                  (set! (.-fillStyle ctx) "#8a96ab")
                  (.fillText ctx "R" (+ bx 6) (+ axis-y 3))
                  (when (pos? (stage 1))
                    (doseq [{:keys [x why]} (:critical d)
                            :when (<= vx0 x vx1)
                            :let [[sx] (px [x 0])]]
                      (set! (.-fillStyle ctx) (why-colors (:kind (first why))))
                      (.beginPath ctx) (.arc ctx sx axis-y 4 0 (* 2 math/PI)) (.fill ctx)
                      (.fillText ctx (.toFixed x 2) (- sx 12) (+ axis-y 16)))))
                (let [cells (cad/cells d)
                      {v 0 e 1 f 2} (frequencies (map :dim cells))
                      cell (when hover (cad/cell d hover))]
                  (set! (.-textContent caption)
                        (str (.-scene controls) ", " blurb
                             "\n" (count (:critical d)) " critical values; " f " 2-cells, " e " 1-cells, " v " 0-cells"
                             "; the free space in " components (if (= 1 components) " piece" " pieces")
                             (when route (str "; route " (count route) " points"))
                             (when (and (pos? (stage 5)) (nil? route)) "; the two ends are not joined")
                             (when cell
                               (str "\ncell " (:id cell) ", dimension " (:dim cell) ": "
                                    (apply str (interpose ", " (map (fn [nm s] (str nm " " (case s 1 "+" -1 "-" "0"))) names (:signs cell))))
                                    (if (free (:id cell)) "  free" "")))
                             "\ngold: fold  red: crossing  blue: leading coefficient or vertical line   hover: a cell's signs  click: move an end")))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (swap! model update :t #(min 6.0 (+ % (* 0.004 (.-speed controls)))))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))
            (at-pointer [^js e]
              (let [r (.getBoundingClientRect canvas)
                    {:keys [world]} (layout)
                    [wx wy :as p] (world (- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r)))
                    [[vx0 vy0] [vx1 vy1]] (:view (scene))]
                (when (and (<= vx0 wx vx1) (<= vy0 wy vy1)) p)))]
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (swap! model assoc :hover (some->> (at-pointer e) (cad/locate (:d @model))))))
      (.addEventListener canvas "pointerleave" (fn [_] (swap! model assoc :hover nil)))
      (.addEventListener canvas "click"
                         (fn [^js e]
                           (when-let [p (at-pointer e)]
                             (let [{:keys [rm ends next]} @model]
                               (when ((:graph rm) (cad/locate rm p))
                                 (swap! model assoc :next (- 1 next) :t (max 5.0 (:t @model)))
                                 (plan! (assoc ends next p)))))))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "scene" (clj->js (keys scenes))) restart!)
        (.onChange (.add gui controls "free" #js ["outside the obstacles" "inside them"]) restart!)
        (.add gui controls "numbers")
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
                   (when-let [container (js/document.getElementById "algebraic-decomposition")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
