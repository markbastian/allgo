(ns allgo.demo.canny-roadmap
  "Canvas demo for `allgo.robotics.canny`: a roadmap of a surface made of
  its critical curves, after LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, section 6.4.3.

  The surface is drawn as cross-sections, faint. Then the roadmap is
  built a stage at a time: the silhouette (gold) -- where the surface,
  seen along the view direction, turns away; the critical points (white)
  of the sweep coordinate, where the silhouette turns back in the sweep;
  the slices through them (blue), which on the torus join the inner ring
  to the outer, as in Figures 6.37 to 6.39 -- turn slices off to see the
  roadmap fall apart; and the junctions where slices cross silhouettes.
  Last, a route between two query points (green), declared critical so
  slices through them join them on. Drag to turn the view; click the
  surface to move a query point."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.numerics.multivariate :as mv]
            [allgo.robotics.canny :as canny]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:surface "torus" :slices true :speed 3})

(def ^:private x (mv/variable 3 0))
(def ^:private y (mv/variable 3 1))
(def ^:private z (mv/variable 3 2))

(defn- torus [R r]
  (let [s (mv/add (mv/power x 2) (mv/power y 2) (mv/power z 2) (mv/constant 3 (- (* R R) (* r r))))]
    (mv/sub (mv/power s 2) (mv/scale (mv/add (mv/power x 2) (mv/power y 2)) (* 4 R R)))))

(def ^:private scenes
  {"sphere" {:f (mv/add (mv/power x 2) (mv/power y 2) (mv/power z 2) (mv/constant 3 -1))
             :bounds [[-1.5 -1.5 -1.5] [1.5 1.5 1.5]] :frame canny/axes
             :queries [[0.3 0.5 0.8] [-0.2 -0.7 -0.6]] :blurb "Example 6.7"}
   "torus" {:f (torus 1.0 0.4) :bounds [[-1.6 -1.6 -0.6] [1.6 1.6 0.6]] :frame canny/axes
            :queries [[0.1 1.0 0.5] [-0.6 -0.2 -0.3]] :blurb "Example 6.8"}
   "double torus" {:f (mv/add (mv/power (mv/sub (mv/mul (mv/power x 2) (mv/sub (mv/constant 3 1) (mv/power x 2))) (mv/power y 2)) 2)
                              (mv/power z 2) (mv/constant 3 -0.01))
                   :bounds [[-1.3 -0.8 -0.3] [1.3 0.8 0.3]] :frame (canny/turned 0.05 0.03 0.07)
                   :queries [[-0.7 0.5 0.08] [0.7 -0.5 -0.08]] :blurb "a tube around a figure eight, genus two"}})

(defn- wire
  "Cross-sections of the surface, for drawing it: planes across x and y."
  [s bounds]
  (let [[lo hi] bounds
        step (/ (v/distance lo hi) 60)
        inside? (fn [p] (every? true? (map <= lo p hi)))
        opts {:step step :stops [] :inside? inside? :limit 2000 :reach (* 0.25 (v/distance lo hi))}]
    (vec (for [[k frame] [[0 canny/axes] [1 {:sweep [0.0 1.0 0.0] :across [0.0 0.0 1.0] :view [1.0 0.0 0.0]}]]
               i (range 1 12)
               :let [c (+ (nth lo k) (* (- (nth hi k) (nth lo k)) (/ i 12)))]
               curve (canny/curves (canny/slice-system s frame c) (canny/plane-seeds frame bounds c 10) opts)]
           curve))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        view (atom {:yaw 0.6 :pitch 1.0 :drag nil :moved? false})]
    (set! (.-style canvas) "display:block;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(scene [] (scenes (.-surface controls)))
            (build! [queries]
              (let [{:keys [f bounds frame]} (scene)
                    rm (canny/roadmap f bounds {:frame frame :queries queries :slices? (.-slices controls)})]
                (swap! model assoc :rm rm :queries queries :route (canny/route rm 0 1)
                       :components (count (canny/components rm))
                       ;; a query no slice joins on is a component alone; the roadmap's own pieces are the rest
                       :pieces (count (remove (fn [c] (every? #(= :query (first %)) c)) (canny/components rm))))))
            (restart! []
              (let [{:keys [f bounds queries]} (scene)
                    s (canny/surface f)]
                (reset! model {:surface s :wire (wire s bounds) :t 0.0 :next-query 0})
                (build! (mapv #(canny/onto-surface s %) queries))))
            (project [p]
              (let [{:keys [yaw pitch]} @view
                    [px py pz] p
                    ;; turn about z by the yaw, then tip by the pitch
                    x1 (- (* px (math/cos yaw)) (* py (math/sin yaw)))
                    y1 (+ (* px (math/sin yaw)) (* py (math/cos yaw)))
                    y2 (- (* y1 (math/cos pitch)) (* pz (math/sin pitch)))
                    z2 (+ (* y1 (math/sin pitch)) (* pz (math/cos pitch)))
                    cw (.-width canvas) ch (.-height canvas)
                    {[lo hi] :bounds} (scene)
                    k (/ (* 0.85 (min (- cw 260) (- ch 120))) (reduce max (map - hi lo)))]
                [(+ (* 0.5 (- cw 230)) (* k x1)) (- (+ 40 (* 0.5 (- ch 100))) (* k z2)) y2]))
            (polyline! [pts closed? color width]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
              (.beginPath ctx)
              (doseq [[i p] (map-indexed vector pts) :let [[sx sy] (project p)]]
                (if (zero? i) (.moveTo ctx sx sy) (.lineTo ctx sx sy)))
              (when closed? (.closePath ctx))
              (.stroke ctx))
            (revealed [curves frac]
              ;; the first `frac` of each curve
              (for [{:keys [points closed?]} curves
                    :let [n (max 1 (long (math/ceil (* frac (count points)))))]]
                [(subvec points 0 (min n (count points))) (and closed? (>= frac 1))]))
            (dot! [p r color]
              (let [[sx sy] (project p)]
                (set! (.-fillStyle ctx) color) (.beginPath ctx) (.arc ctx sx sy r 0 (* 2 math/PI)) (.fill ctx)))
            (draw! []
              (let [{:keys [wire rm t route components pieces queries]} @model
                    {:keys [silhouettes critical slices junctions frame]} rm
                    cw (.-width canvas) ch (.-height canvas)
                    stage (long (math/floor t))
                    frac #(cond (< stage %) nil (> stage %) 1.0 :else (- t stage))]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (doseq [{:keys [points closed?]} wire]
                  (polyline! points closed? "rgba(140,170,220,0.22)" 1))
                (when-let [f (frac 0)]
                  (doseq [[pts closed?] (revealed silhouettes f)] (polyline! pts closed? "#ffd479" 2.5)))
                (when-let [f (frac 1)]
                  (doseq [{:keys [point]} critical] (dot! point (* 6 (min 1 (* 2 f))) "#ffffff")))
                (when-let [f (frac 2)]
                  (doseq [[pts closed?] (revealed (mapcat :curves slices) f)] (polyline! pts closed? "#7fd1ff" 2.5)))
                (when (frac 3)
                  (doseq [{:keys [point]} junctions] (dot! point 4 "#ff8f70")))
                (when-let [f (frac 4)]
                  (when route (polyline! (subvec route 0 (max 1 (long (* f (count route))))) false "#7dffa0" 4)))
                (doseq [q queries] (dot! q 6 "#7dffa0"))
                ;; the sweep direction, from the middle of the box
                (let [{[lo hi] :bounds} (scene)
                      o (v/scale (v/add lo hi) 0.5)
                      tip (v/add o (v/scale (:sweep frame) (* 0.55 (- (nth hi 0) (nth lo 0)))))]
                  (polyline! [(v/add o (v/scale (:sweep frame) (* -0.55 (- (nth hi 0) (nth lo 0))))) tip] false "rgba(255,255,255,0.25)" 1)
                  (let [[sx sy] (project tip)]
                    (set! (.-fillStyle ctx) "rgba(255,255,255,0.5)") (set! (.-font ctx) "12px ui-monospace,Menlo,monospace")
                    (.fillText ctx "g1" (+ sx 6) sy)))
                (set! (.-textContent caption)
                      (str (.-surface controls) ", " (:blurb (scene)) "  |  "
                           (case (min stage 4)
                             0 (str "the silhouette: " (count silhouettes) " curve" (when (> (count silhouettes) 1) "s"))
                             1 (str (count critical) " critical points of g1, where the silhouette turns back in the sweep")
                             2 (if (.-slices controls)
                                 (str "slices through the critical points and the queries: " (count (mapcat :curves slices)) " curves")
                                 "no slices: the silhouettes and critical points alone")
                             3 (str (count junctions) " junctions where slices cross silhouettes")
                             4 (if route
                                 (str "a route along the roadmap, " (count route) " points")
                                 "no route: the roadmap is in pieces"))
                           "  |  the roadmap is " (if (= 1 pieces) "connected" (str "in " pieces " pieces"))
                           (when (> components pieces) ", the queries apart from it")
                           "\ngold: silhouette  white: critical points  blue: slices  red: junctions  green: queries and route"
                           "  drag: turn  click: move a query"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (when-not (:drag @view)
                  (swap! view update :yaw + 0.002))
                (swap! model update :t #(min 5.0 (+ % (* 0.004 (.-speed controls)))))
                (draw!)))
            (pick [sx sy]
              ;; the nearest cross-section point under the pointer, nearest the viewer
              (let [{:keys [wire surface]} @model
                    hits (for [{:keys [points]} wire p points
                               :let [[px py depth] (project p)]
                               :when (< (+ (math/pow (- px sx) 2) (math/pow (- py sy) 2)) 144)]
                           [p depth])]
                (when (seq hits)
                  (canny/onto-surface surface (first (apply min-key second hits))))))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e] (swap! view assoc :drag [(.-clientX e) (.-clientY e)] :moved? false)))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (when-let [[x0 y0] (:drag @view)]
                             (let [dx (- (.-clientX e) x0) dy (- (.-clientY e) y0)]
                               (when (> (+ (abs dx) (abs dy)) 2)
                                 (swap! view (fn [w] (-> w
                                                         (update :yaw + (* 0.01 dx))
                                                         (update :pitch #(max 0.0 (min math/PI (+ % (* 0.01 dy)))))
                                                         (assoc :drag [(.-clientX e) (.-clientY e)] :moved? true)))))))))
      (.addEventListener canvas "pointerup"
                         (fn [^js e]
                           (let [moved? (:moved? @view)]
                             (swap! view assoc :drag nil)
                             (when-not moved?
                               (let [r (.getBoundingClientRect canvas)]
                                 (when-let [p (pick (- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r)))]
                                   (let [{:keys [queries next-query]} @model]
                                     (swap! model assoc :next-query (- 1 next-query) :t 4.0)
                                     (build! (assoc queries next-query p)))))))))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "surface" (clj->js (keys scenes))) restart!)
        (.onChange (.add gui controls "slices") restart!)
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
                   (when-let [container (js/document.getElementById "canny-roadmap")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
