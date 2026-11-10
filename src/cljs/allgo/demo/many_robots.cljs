(ns allgo.demo.many-robots
  "Canvas demo for `allgo.robotics.coordination`: planning with time and
  with other robots, after LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, sections 7.1 and 7.2.

  Velocity tuning: the robot's path across the road is fixed, and cars
  cross it; it chooses only when to go, when to wait and how fast. The
  inset is the plan's space, how far along the path (up) by time
  (across): dark where a car would be too close, and the schedule a way
  through the light that never goes back in time. The coordination
  diagram does the same for two robots on fixed paths that cross twice:
  the square of their two progresses, dark where they would collide,
  crossed from corner to corner -- each may wait, neither go back.
  Prioritized planning moves a fleet through a warehouse, each robot
  planned in space and time around the ones planned before it: fast,
  but a robot planned late can find its way blocked for good."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.random :as random]
            [allgo.robotics.coordination :as coordination]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:scene "velocity tuning" :seed 1 :speed 4})

;; ---------------------------------------------------------------------------
;; Velocity tuning

(def ^:private road-path [[0.5 4.0] [11.5 4.0]])
(def ^:private tuning {:ds 0.1 :dt 0.1 :max-step 1 :clearance 0.8 :horizon 600})

(defn- car
  "A car driving along x = `x` at `v`, looping over the room's height."
  [x v y0]
  (fn [t] [x (- (mod (+ y0 (* v t)) 10.0) 1.0)]))

(defn- tuning-scene [rng]
  (let [cars (mapv (fn [x] (car x (* (if (< (rng) 0.5) -1 1) (+ 1.0 (* 1.5 (rng)))) (* 10 (rng)))) [2.5 4.5 6.5 8.5 10.0])
        pts (coordination/samples road-path (:ds tuning))
        schedule (coordination/velocity-tuning road-path cars tuning)
        ticks (+ 20 (long (/ (or (second (peek schedule)) 0.0) (:dt tuning))))]
    {:cars cars :pts pts :schedule schedule
     :blocked (vec (for [k (range ticks)]
                     (let [t (* k (:dt tuning))]
                       (mapv (fn [p] (some #(<= (polygon/distance p (% t)) (:clearance tuning)) cars)) pts))))}))

;; ---------------------------------------------------------------------------
;; Two robots

(def ^:private path-a [[1.0 2.0] [6.0 6.5] [11.0 2.0]])
(def ^:private path-b [[1.0 6.0] [6.0 1.5] [11.0 6.0]])

(defn- diagram-scene []
  (let [{[na nb] :size blocked :blocked :as d} (coordination/coordination-diagram path-a path-b 0.1 0.8)]
    {:diagram d :schedule (coordination/coordinate d)
     ;; grid[i][j], for the bitmap
     :grid (vec (for [i (range na)] (vec (for [j (range nb)] (contains? blocked [i j])))))
     :a (coordination/samples path-a 0.1) :b (coordination/samples path-b 0.1)}))

;; ---------------------------------------------------------------------------
;; A warehouse

(def ^:private gw 24)
(def ^:private gh 14)
(def ^:private shelves
  (set (for [x0 [3 8 13 18] y0 [3 9] x (range x0 (+ x0 3)) y [y0 (inc y0)]] [x y])))
(defn- floor? [[x y]] (and (< -1 x gw) (< -1 y gh) (not (shelves [x y]))))

(defn- fleet-scene [rng]
  (let [side (fn [x] (vec (map #(vector x %) (random/shuffle rng (vec (range gh))))))
        left (side 0) right (side (dec gw))
        robots (vec (concat (map vector (take 4 left) (take 4 right)) (map vector (take 4 (drop 4 right)) (take 4 (drop 4 left)))))
        robots (random/shuffle rng robots)
        paths (coordination/prioritized floor? robots 120)]
    {:robots robots :paths paths}))

(def ^:private colors ["#ffd479" "#7fd1ff" "#ff8f70" "#7dffa0" "#d8a6ff" "#ffb3d9" "#9fe8d8" "#f0f0a0"])

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        st-img (js/document.createElement "canvas")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls))]
                (reset! model {:scene (case (.-scene controls)
                                        "velocity tuning" (tuning-scene rng)
                                        "coordination diagram" (diagram-scene)
                                        "prioritized planning" (fleet-scene rng))
                               :t 0.0 :drawn nil})))
            (layout [frac]
              ;; the room, 12 by 8, in the left `frac` of the canvas
              (let [cw (.-width canvas) ch (.-height canvas)
                    avail (* frac cw)
                    scale (min (/ (- avail 40) 12.0) (/ (- ch 110) 8.0))
                    x0 (* 0.5 (- avail (* 12 scale))) y0 40]
                {:scale scale :px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 8 y)))])}))
            (polyline! [px pts color width]
              (when (next pts)
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                (.beginPath ctx)
                (let [[x y] (px (first pts))] (.moveTo ctx x y))
                (doseq [p (rest pts) :let [[x y] (px p)]] (.lineTo ctx x y))
                (.stroke ctx)))
            (dot! [px p r color]
              (let [[x y] (px p)]
                (set! (.-fillStyle ctx) color)
                (.beginPath ctx) (.arc ctx x y r 0 (* 2 math/PI)) (.fill ctx)))
            (room! [px]
              (let [[ax ay] (px [0 8]) [bx by] (px [12 0])]
                (set! (.-fillStyle ctx) "#161b25")
                (.fillRect ctx ax ay (- bx ax) (- by ay))))
            (inset [cw ch]
              ;; a square below the controls, in the right two fifths
              (let [side (min (- (* 0.4 cw) 40) (- ch 200))]
                [(+ (* 0.6 cw) (* 0.5 (- (* 0.4 cw) side))) 130 side]))
            (bitmap! [grid]
              ;; grid[x][y] true where blocked, cached
              (when-not (identical? grid (:drawn @model))
                (let [w (count grid) h (count (first grid))
                      c2 (.getContext st-img "2d")]
                  (set! (.-width st-img) w) (set! (.-height st-img) h)
                  (let [img (.createImageData c2 w h) data (.-data img)]
                    (dotimes [x w]
                      (dotimes [y h]
                        (let [o (* 4 (+ x (* w (- h 1 y)))) hit (get-in grid [x y])]
                          (aset data o (if hit 60 30)) (aset data (+ o 1) (if hit 70 45))
                          (aset data (+ o 2) (if hit 95 60)) (aset data (+ o 3) 255))))
                    (.putImageData c2 img 0 0))
                  (swap! model assoc :drawn grid))))
            (draw-tuning! [{:keys [cars pts schedule blocked]} t]
              (let [{:keys [px scale]} (layout 0.6)
                    cw (.-width canvas) ch (.-height canvas)
                    dt (:dt tuning)
                    end (if (seq schedule) (second (peek schedule)) 0.0)
                    t (mod t (+ end 2.0))
                    [i] (or (last (take-while #(<= (second %) t) schedule)) [0])]
                (room! px)
                (polyline! px road-path "rgba(255,212,121,0.4)" 2)
                (doseq [c cars]
                  (let [p (c t)]
                    (dot! px p (* scale (:clearance tuning)) "rgba(255,143,112,0.15)")
                    (dot! px p (* scale 0.3) "#ff8f70")))
                (dot! px (pts i) (* scale 0.25) "#ffffff")
                (bitmap! blocked)
                (let [[x0 y0 side] (inset cw ch)
                      k-max (count blocked) n (count pts)
                      ipx (fn [[ii tt]] [(+ x0 (* side (/ tt (* dt k-max)))) (+ y0 (* side (- 1 (/ ii (dec n)))))])]
                  (set! (.-imageSmoothingEnabled ctx) false)
                  (.drawImage ctx st-img x0 y0 side side)
                  (polyline! ipx schedule "#ffd479" 2)
                  (dot! ipx [i (min t end)] 4 "#ffffff")
                  (set! (.-fillStyle ctx) "#9fb4d8") (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
                  (.fillText ctx "distance along the path (up) by time (across)" x0 (- y0 8)))
                (if (seq schedule)
                  (str "across in " (.toFixed end 1) " s, " (count (filter (fn [[[a] [b]]] (= a b)) (partition 2 1 schedule)))
                       " ticks of waiting; the cars' discs are the clearance it keeps")
                  "no schedule gets across")))
            (draw-diagram! [{:keys [diagram schedule a b grid]} t]
              (let [{:keys [px scale]} (layout 0.6)
                    cw (.-width canvas) ch (.-height canvas)
                    [na nb] (:size diagram)
                    k (mod (long t) (+ 30 (count schedule)))
                    [i j] (if (seq schedule) (schedule (min k (dec (count schedule)))) [0 0])]
                (room! px)
                (polyline! px path-a "rgba(127,209,255,0.5)" 2)
                (polyline! px path-b "rgba(255,143,112,0.5)" 2)
                (dot! px (a i) (* scale 0.4) "#7fd1ff")
                (dot! px (b j) (* scale 0.4) "#ff8f70")
                (bitmap! grid)
                (let [[x0 y0 side] (inset cw ch)
                      dpx (fn [[ii jj]] [(+ x0 (* side (/ ii (dec na)))) (+ y0 (* side (- 1 (/ jj (dec nb)))))])]
                  (set! (.-imageSmoothingEnabled ctx) false)
                  (.drawImage ctx st-img x0 y0 side side)
                  (polyline! dpx schedule "#ffd479" 2)
                  (dot! dpx [i j] 4 "#ffffff")
                  (set! (.-fillStyle ctx) "#9fb4d8") (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
                  (.fillText ctx "blue's progress (across) by coral's (up)" x0 (- y0 8)))
                (if (seq schedule)
                  (str (count (:blocked diagram)) " of " (* na nb) " pairs of positions collide; the schedule takes "
                       (dec (count schedule)) " ticks against " (dec (max na nb)) " for either alone")
                  "each blocks the other")))
            (draw-fleet! [{:keys [robots paths]} t]
              (let [cw (.-width canvas) ch (.-height canvas)
                    cell (min (/ (- cw 40) gw) (/ (- ch 110) gh))
                    x0 (* 0.5 (- cw (* gw cell))) y0 40
                    center (fn [[x y]] [(+ x0 (* cell (+ x 0.5))) (+ y0 (* cell (- gh y 0.5)))])
                    longest (reduce max 1 (map count (filter some? paths)))
                    t (mod t (+ longest 10))
                    at (fn [p] (let [k (long t) f (- t k)
                                     a (p (min k (dec (count p)))) b (p (min (inc k) (dec (count p))))]
                                 [(+ (a 0) (* f (- (b 0) (a 0)))) (+ (a 1) (* f (- (b 1) (a 1))))]))]
                (doseq [x (range gw) y (range gh)
                        :let [[cx cy] (center [x y])]]
                  (set! (.-fillStyle ctx) (if (shelves [x y]) "#3a4458" "#161b25"))
                  (.fillRect ctx (+ (- cx (* 0.5 cell)) 0.5) (+ (- cy (* 0.5 cell)) 0.5) (- cell 1) (- cell 1)))
                (doseq [[k [[_ goal] p]] (map-indexed vector (map vector robots paths))
                        :let [color (colors k) [gx gy] (center goal)]]
                  (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 2)
                  (.strokeRect ctx (- gx (* 0.35 cell)) (- gy (* 0.35 cell)) (* 0.7 cell) (* 0.7 cell))
                  (if p
                    (do (polyline! center p (str color "55") 1.5)
                        (dot! center (at p) (* 0.38 cell) color)
                        (let [[x y] (center (at p))]
                          (set! (.-fillStyle ctx) "#0d1118") (set! (.-font ctx) "bold 10px ui-monospace,Menlo,monospace")
                          (.fillText ctx (str (inc k)) (- x 3) (+ y 3))))
                    (let [[sx sy] (center (first (robots k)))]
                      (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 2)
                      (.beginPath ctx) (.moveTo ctx (- sx 6) (- sy 6)) (.lineTo ctx (+ sx 6) (+ sy 6))
                      (.moveTo ctx (+ sx 6) (- sy 6)) (.lineTo ctx (- sx 6) (+ sy 6)) (.stroke ctx))))
                (str (count (filter some? paths)) " of " (count robots) " robots planned, numbered by priority, "
                     (count (coordination/conflicts (vec paths))) " collisions"
                     (when (some nil? paths) "; a cross marks one planned too late to find a way"))))
            (draw! []
              (let [{:keys [scene t]} @model
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (let [label (case (.-scene controls)
                              "velocity tuning" (draw-tuning! scene t)
                              "coordination diagram" (draw-diagram! scene t)
                              "prioritized planning" (draw-fleet! scene t))]
                  (set! (.-textContent caption)
                        (str (.-scene controls) "  |  " label
                             "\n" (case (.-scene controls)
                                    "velocity tuning" "white: the robot  coral: the cars  inset: dark where a car is too close, gold the schedule"
                                    "coordination diagram" "inset: dark where the two collide; gold the schedule, never back in either"
                                    "prioritized planning" "squares: the goals  faint: each robot's planned path"))))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (swap! model update :t + (* (.-speed controls)
                                            (case (.-scene controls)
                                              "velocity tuning" 0.01
                                              "coordination diagram" 0.25
                                              "prioritized planning" 0.02)))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "scene" #js ["velocity tuning" "coordination diagram" "prioritized planning"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 20 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "many-robots")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
