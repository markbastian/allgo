(ns allgo.demo.grid-navigation
  "Canvas demo for `allgo.search.incremental` and
  `allgo.robotics.navigation`: planning on a grid when the map is wrong,
  and plans that are feedback rather than paths -- after LaValle,
  *Planning Algorithms*, Cambridge University Press, 2006, sections
  12.3.2 and 8.2.

  D* Lite's robot believes the unseen grid is clear. It sees walls only
  within two cells of itself (outlined faintly until then), and each
  time it finds one in its way it repairs its plan (gold) rather than
  planning again from nothing -- the count of cells it has expanded
  stays far below what a fresh search at every move would cost. The
  wavefront (NF1) labels every cell with its distance from the goal, a
  ring at a time; the optimal navigation function does the same with
  diagonal moves weighed root two. Either way every cell knows which way
  is downhill, with no local minima to stall in: the robots dropped
  around the map all roll down to the goal. Click to move the goal."
  (:require [allgo.random :as random]
            [allgo.robotics.navigation :as navigation]
            [allgo.search.incremental :as incremental]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "D* Lite" :seed 1 :speed 2})

(def ^:private cols 40)
(def ^:private rows 26)
(def ^:private start [1 1])

(defn- walls-for
  "Walls for `rng`: four vertical walls, each with a gap or two, and a few
  horizontal runs, kept off the start."
  [rng]
  (let [vertical (for [x [8 16 24 32]
                       :let [gaps (set (mapcat (fn [_] (let [g (random/below rng (- rows 3))] [g (inc g) (+ g 2)])) (range (inc (random/below rng 2)))))]
                       y (range rows) :when (not (gaps y))]
                   [x y])
        horizontal (for [_ (range 7)
                         :let [y (+ 2 (random/below rng (- rows 4))) x0 (random/below rng (- cols 10)) len (+ 4 (random/below rng 8))]
                         x (range x0 (+ x0 len))]
                     [x y])]
    (disj (set (concat vertical horizontal)) start [2 1] [1 2])))

(defn- inside? [[x y]] (and (< -1 x cols) (< -1 y rows)))

(defn- maze
  "The first walls `rng` gives that leave a way from the start to `goal`
  -- the wavefront from the goal reaching the start."
  [rng goal]
  (first (filter (fn [walls]
                   (contains? (navigation/wavefront #(and (inside? %) (not (walls %))) goal) start))
                 (map #(disj % goal) (repeatedly #(walls-for rng))))))

(defn- nbrs4 [[x y]]
  (filter inside? [[(inc x) y] [(dec x) y] [x (inc y)] [x (dec y)]]))

(defn- manhattan [[a b] [c d]] (+ (abs (- a c)) (abs (- b d))))

(defn- sense
  "The walls within two cells of `p`, as the edges into them that turn out
  impassable."
  [walls [x y]]
  (for [dx (range -2 3) dy (range -2 3)
        :let [q [(+ x dx) (+ y dy)]]
        :when (walls q)
        n (nbrs4 q)]
    [q n ##Inf]))

(defn- with-trail
  "The next D* Lite state, with the cells the robot has been in."
  [{:keys [trail]} st]
  {:st st :trail (conj trail (:start st))})

(defn- states [method walls goal rng]
  (let [free? (fn [c] (and (inside? c) (not (walls c))))
        drops (vec (take 10 (filter free? (repeatedly #(vector (random/below rng cols) (random/below rng rows))))))]
    (case method
      "D* Lite"
      (let [problem {:start start :goal goal :neighbors nbrs4 :cost (fn [_ _] 1.0) :heuristic manhattan}
            [s0 & more] (incremental/navigate-seq problem (partial sense walls))]
        (reductions with-trail {:st s0 :trail [(:start s0)]} more))
      "NF1 wavefront"
      (let [waves (navigation/wavefront-seq free? goal)
            nf (:values (last waves))]
        (concat waves (repeat {:values nf :paths (keep #(navigation/descend nf free? %) drops)})))
      "optimal navigation function"
      (let [nf (navigation/optimal-navigation free? goal)]
        (repeat {:values nf :paths (keep #(navigation/descend nf free? %) drops)})))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        goal (atom [38 24])]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls))
                    walls (maze rng @goal)
                    ss (states (.-method controls) walls @goal rng)]
                (reset! model {:walls walls :rest ss :st (first ss) :t 0.0})))
            (advance! []
              (swap! model (fn [{:keys [rest] :as m}]
                             (let [ahead (nthrest rest (.-speed controls))]
                               (cond (seq ahead) (assoc m :rest ahead :st (first ahead))
                                     (next rest) (let [l (last rest)] (assoc m :rest (list l) :st l))
                                     :else m)))))
            (layout []
              (let [cw (.-width canvas) ch (.-height canvas)
                    cell (min (/ (- cw 40) cols) (/ (- ch 110) rows))
                    x0 (* 0.5 (- cw (* cols cell))) y0 40]
                {:cell cell :x0 x0 :y0 y0
                 ;; the top-left corner of cell [x y], y up
                 :corner (fn [[x y]] [(+ x0 (* cell x)) (+ y0 (* cell (- rows 1 y)))])
                 :center (fn [[x y]] [(+ x0 (* cell (+ x 0.5))) (+ y0 (* cell (- rows y 0.5)))])}))
            (fill-cell! [corner cell c color]
              (let [[x y] (corner c)]
                (set! (.-fillStyle ctx) color)
                (.fillRect ctx (+ x 0.5) (+ y 0.5) (- cell 1) (- cell 1))))
            (polyline! [center cells color width]
              (when (next cells)
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                (.beginPath ctx)
                (let [[x y] (center (first cells))] (.moveTo ctx x y))
                (doseq [c (rest cells) :let [[x y] (center c)]] (.lineTo ctx x y))
                (.stroke ctx)))
            (dot! [center c r color]
              (let [[x y] (center c)]
                (set! (.-fillStyle ctx) color)
                (.beginPath ctx) (.arc ctx x y r 0 (* 2 math/PI)) (.fill ctx)))
            (along [path t]
              ;; down the path, a pause at the goal, and again
              (path (min (dec (count path)) (mod (long t) (+ 15 (count path))))))
            (draw! []
              (let [{:keys [walls st t]} @model
                    method (.-method controls)
                    {:keys [cell corner center]} (layout)
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (doseq [x (range cols) y (range rows)] (fill-cell! corner cell [x y] "#161b25"))
                (let [label
                      (case method
                        "D* Lite"
                        (let [{:keys [st trail]} st
                              known (set (map first (keys (:changed st))))
                              plan (incremental/planned-path st)]
                          (doseq [w walls]
                            (if (known w)
                              (fill-cell! corner cell w "#5b6884")
                              (let [[x y] (corner w)]
                                (set! (.-strokeStyle ctx) "rgba(91,104,132,0.5)") (set! (.-lineWidth ctx) 1)
                                (.strokeRect ctx (+ x 1.5) (+ y 1.5) (- cell 3) (- cell 3)))))
                          (polyline! center plan "#ffd479" 2)
                          (polyline! center trail "rgba(125,255,160,0.8)" 2)
                          (dot! center (:start st) (* 0.4 cell) "#ffffff")
                          (str (dec (count trail)) " moves, " (count known) " walls seen, " (:expanded st) " cells expanded in all"
                               (cond (= (:start st) (:goal st)) " -- arrived"
                                     (seq plan) (str ", plan " (dec (count plan)) " moves as the map now stands")
                                     :else " -- no way through")))

                        ("NF1 wavefront" "optimal navigation function")
                        (let [{:keys [values wave paths]} st
                              top (max 1.0 (reduce max 0.0 (vals values)))]
                          (doseq [[c v] values]
                            (fill-cell! corner cell c (str "hsla(" (- 210 (* 200 (/ v top))) ",70%,55%,0.55)")))
                          (doseq [c wave] (fill-cell! corner cell c "rgba(255,255,255,0.6)"))
                          (doseq [p paths]
                            (polyline! center p "rgba(255,212,121,0.6)" 1.5)
                            (dot! center (along p t) (* 0.4 cell) "#ffffff"))
                          (if paths
                            (str (count values) " cells labeled, the farthest " (.toFixed top 1)
                                 (if (= method "NF1 wavefront") " steps" " m")
                                 " from the goal; each robot steps to its lowest neighbor and every one arrives")
                            (str "ring " (long (get values (first wave) 0)) ", " (count values) " cells labeled"))))]
                  (doseq [w walls :when (not= method "D* Lite")] (fill-cell! corner cell w "#5b6884"))
                  (dot! center start (* 0.3 cell) "#7dffa0")
                  (dot! center @goal (* 0.3 cell) "#ff8f70")
                  (set! (.-textContent caption)
                        (str method "  |  " label
                             "\ngreen: start  red: goal  gold: the plan  click: move the goal")))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (swap! model update :t + (* 0.1 (.-speed controls)))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (let [{:keys [cell x0 y0]} (layout)
                                 r (.getBoundingClientRect canvas)
                                 x (long (math/floor (/ (- (- (.-clientX e) (.-left r)) x0) cell)))
                                 y (- rows 1 (long (math/floor (/ (- (- (.-clientY e) (.-top r)) y0) cell))))]
                             (when (and (inside? [x y]) (not ((:walls @model) [x y])))
                               (reset! goal [x y])
                               (restart!)))))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "method" #js ["D* Lite" "NF1 wavefront" "optimal navigation function"]) restart!)
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
                   (when-let [container (js/document.getElementById "grid-navigation")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
