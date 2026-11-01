(ns allgo.demo.motion-planning
  "Canvas demo for `allgo.robotics.motion`: the motion planners of Russell
  and Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 26, in one room.

  The robot (white) is a disk; the obstacles are drawn as they are and,
  faintly, grown by its radius -- the configuration space, where it is a
  point. The visibility graph joins the grown corners that see each
  other and gives the shortest path; the Voronoi graph runs between the
  obstacles, equally far from each, and gives the safest. Cell decomposition searches the
  wholly free cells, shown as the A* search spreads through them. PRM
  scatters a roadmap and searches it; RRT grows a tree toward random
  samples until it touches the goal; RRT* keeps growing and rewiring,
  and its path keeps shortening. The potential field slides downhill
  and can stall. The two-link arm plans in its own configuration space,
  drawn beside it: dark where the arm would hit something. Click to move
  the goal."
  (:require [allgo.random :as random]
            [allgo.robotics.motion :as motion]
            [allgo.search :as search]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "RRT*" :seed 1 :speed 4})

(def ^:private world
  {:bounds [[0.0 0.0] [12.0 8.0]] :radius 0.25
   :obstacles [[[2.0 0.0] [2.6 0.0] [2.6 5.5] [2.0 5.5]]
               [[4.5 2.5] [6.5 2.0] [7.0 4.0] [5.0 4.8]]
               [[5.0 6.2] [9.0 6.2] [9.0 6.8] [5.0 6.8]]
               [[8.5 0.8] [10.0 0.8] [10.0 4.8] [9.4 4.8] [9.4 1.4] [8.5 1.4]]
               [[3.6 6.0] [4.2 6.0] [4.2 8.0] [3.6 8.0]]]})

(def ^:private start [1.0 1.0])

(def ^:private arm (motion/arm [6.0 4.0] [2.2 1.6]))
(def ^:private arm-world
  {:obstacles [[[8.0 5.0] [9.0 5.0] [9.0 6.0] [8.0 6.0]]
               [[3.0 1.5] [4.0 1.5] [4.0 2.6] [3.0 2.6]]
               [[6.6 0.6] [7.6 0.6] [7.6 1.4] [6.6 1.4]]]})
(def ^:private arm-n 90)

(defn- states [method goal]
  (let [rng (random/rng (.-seed controls))]
    (case method
      "visibility graph" (let [g (motion/visibility-graph world start goal)]
                           (repeat {:graph g :path (motion/visibility-path world start goal)}))
      "Voronoi graph" (let [g (motion/voronoi-graph world 0.15)]
                        (repeat {:graph g :path (motion/voronoi-path world g start goal)}))
      "cell decomposition" (let [grid (motion/cell-grid world 0.25)
                                 ;; up to and including the step that reaches the goal
                                 [before after] (split-with (complement :solved?) (motion/cell-seq grid start goal))]
                             (map (fn [st] {:grid grid :search st}) (concat before (take 1 after))))
      "PRM" (map-indexed (fn [k rm] {:roadmap rm :path (when (and (pos? k) (zero? (mod k 10)))
                                                         (motion/prm-path world rm start goal))})
                         (motion/prm-seq world rng))
      "RRT" (motion/rrt-seq world start goal rng)
      "RRT*" (motion/rrt-seq world start goal rng :star? true)
      "potential field" (reductions conj [] (take 2000 (motion/potential-seq world start goal :influence 0.8 :repulse 0.1)))
      "two-link arm" (let [grid (motion/arm-cspace arm arm-world arm-n)]
                       (repeat {:cspace grid :plan (motion/arm-plan grid [-3.0 0.0] [0.0 0.6])})))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        cspace-img (js/document.createElement "canvas")
        running? (atom false)
        model (atom nil)
        goal (atom [11.0 7.2])]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [ss (states (.-method controls) @goal)]
                (reset! model {:rest ss :st (first ss) :step 0 :path-t 0.0 :drawn-cspace nil})))
            (advance! []
              (swap! model (fn [{:keys [rest step] :as m}]
                             (let [n (.-speed controls)
                                   ahead (nthrest rest n)]
                               (cond
                                 (seq ahead) (assoc m :rest ahead :st (first ahead) :step (+ step n))
                                 ;; Fewer than n steps left: land on the last,
                                 ;; which for a search is the one that solved it.
                                 (next rest) (let [l (last rest)] (assoc m :rest (list l) :st l :step (+ step (dec (count rest)))))
                                 :else m)))))
            (layout []
              (let [cw (.-width canvas) ch (.-height canvas)
                    arm? (= "two-link arm" (.-method controls))
                    avail (if arm? (* 0.6 cw) cw)
                    scale (min (/ (- avail 40) 12.0) (/ (- ch 110) 8.0))
                    x0 (* 0.5 (- avail (* 12 scale))) y0 40]
                {:scale scale :x0 x0 :y0 y0
                 :px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 8 y)))])}))
            (poly! [px pts fill stroke]
              (.beginPath ctx)
              (let [[x y] (px (first pts))] (.moveTo ctx x y))
              (doseq [p (rest pts) :let [[x y] (px p)]] (.lineTo ctx x y))
              (.closePath ctx)
              (when fill (set! (.-fillStyle ctx) fill) (.fill ctx))
              (when stroke (set! (.-strokeStyle ctx) stroke) (.stroke ctx)))
            (line! [px a b color width]
              (let [[ax ay] (px a) [bx by] (px b)]
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx)))
            (dot! [px p r color]
              (let [[x y] (px p)]
                (set! (.-fillStyle ctx) color)
                (.beginPath ctx) (.arc ctx x y r 0 (* 2 math/PI)) (.fill ctx)))
            (path! [px path color]
              (doseq [[a b] (partition 2 1 path)] (line! px a b color 3)))
            (robot-on! [px scale path t]
              ;; The robot driven along the path, round and round.
              (when (and (seq path) (> (count path) 1))
                (let [len (motion/path-length path)
                      s (mod t len)
                      p (loop [[[a b] & more] (partition 2 1 path) s s]
                          (let [d (math/sqrt (+ (math/pow (- (b 0) (a 0)) 2) (math/pow (- (b 1) (a 1)) 2)))]
                            (if (or (nil? more) (<= s d))
                              [(+ (a 0) (* (/ (min s d) (max d 1e-9)) (- (b 0) (a 0)))) (+ (a 1) (* (/ (min s d) (max d 1e-9)) (- (b 1) (a 1))))]
                              (recur more (- s d)))))]
                  (dot! px p (* scale (:radius world)) "#ffffff"))))
            (draw-room! [px]
              (let [[[bx0 by0] [bx1 by1]] (:bounds world)]
                (poly! px [[bx0 by0] [bx1 by0] [bx1 by1] [bx0 by1]] "#161b25" nil))
              (doseq [g (motion/grown world)] (poly! px g "rgba(120,140,170,0.18)" nil))
              (set! (.-lineWidth ctx) 1)
              (doseq [o (:obstacles world)] (poly! px o "#3a4458" "#5b6884"))
              (dot! px start 5 "#7dffa0")
              (dot! px @goal 5 "#ff8f70")
              (set! (.-lineWidth ctx) 1))
            (draw! []
              (let [{:keys [st path-t]} @model
                    method (.-method controls)
                    {:keys [px scale]} (layout)
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (if (= method "two-link arm")
                  (draw-arm! st px)
                  (do (draw-room! px)
                      (draw-method! method st px scale path-t)))))
            (draw-method! [method st px scale path-t]
              (let [label
                    (case method
                      "visibility graph"
                      (let [{{:keys [nodes edges]} :graph path :path} st]
                        (doseq [[i j] edges] (line! px (nodes i) (nodes j) "rgba(140,170,220,0.25)" 1))
                        (doseq [n nodes] (dot! px n 2.5 "#9fb4d8"))
                        (path! px path "#ffd479")
                        (robot-on! px scale path path-t)
                        (str (count edges) " sight lines, shortest path " (.toFixed (motion/path-length path) 2) " m"))

                      "Voronoi graph"
                      (let [{{:keys [nodes edges]} :graph path :path} st]
                        (doseq [[i j] edges] (line! px (nodes i) (nodes j) "rgba(140,170,220,0.55)" 1.5))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count edges) " graph edges, the points equidistant from two obstacles"
                             (if path (str "; path " (.toFixed (motion/path-length path) 2) " m, at least "
                                           (.toFixed (reduce min (map #(motion/clearance world %) (subvec path 1 (dec (count path))))) 2)
                                           " m from every obstacle along the graph")
                                 "; the goal cannot reach the graph")))

                      "cell decomposition"
                      (let [{:keys [grid search]} st
                            {:keys [h origin free]} grid
                            [ox oy] origin
                            cell-px (fn [[i j]] (px [(+ ox (* h i)) (+ oy (* h (inc j)))]))
                            seen (keys (:came-from search))
                            p (when (:solved? search)
                                (let [cells (search/path search)]
                                  (vec (concat [start] (map (fn [[i j]] [(+ ox (* h (+ i 0.5))) (+ oy (* h (+ j 0.5)))]) (rest (butlast cells))) [@goal]))))]
                        (set! (.-fillStyle ctx) "rgba(120,200,140,0.10)")
                        (doseq [c free :let [[x y] (cell-px c)]] (.fillRect ctx (+ x 0.5) (+ y 0.5) (- (* scale h) 1) (- (* scale h) 1)))
                        (set! (.-fillStyle ctx) "rgba(120,170,255,0.35)")
                        (doseq [c seen :let [[x y] (cell-px c)]] (.fillRect ctx (+ x 0.5) (+ y 0.5) (- (* scale h) 1) (- (* scale h) 1)))
                        (when p (path! px p "#ffd479") (robot-on! px scale p path-t))
                        (str (count free) " free cells, " (count seen) " reached"
                             (if p (str ", path " (.toFixed (motion/path-length p) 2) " m") ", searching")))

                      "PRM"
                      (let [{{:keys [nodes edges]} :roadmap} st
                            path (:last-path @model)]
                        (doseq [[i j] edges] (line! px (nodes i) (nodes j) "rgba(140,170,220,0.3)" 1))
                        (doseq [n nodes] (dot! px n 2 "#9fb4d8"))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count nodes) " samples, " (count edges) " links"
                             (if path (str ", path " (.toFixed (motion/path-length path) 2) " m") ", not yet connected")))

                      ("RRT" "RRT*")
                      (let [{:keys [nodes parent]} st
                            path (motion/tree-path st)]
                        (doseq [[i p] (map-indexed vector parent) :when p]
                          (line! px (nodes p) (nodes i) "rgba(140,170,220,0.4)" 1))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count nodes) " nodes"
                             (if path (str ", path " (.toFixed (motion/path-length path) 2) " m") ", growing")))

                      "potential field"
                      (let [trail st]
                        (doseq [[a b] (partition 2 1 trail)] (line! px a b "#ffd479" 2))
                        (when (seq trail) (dot! px (peek trail) (* scale (:radius world)) "#ffffff"))
                        (let [p (peek trail)]
                          (if (= p @goal) "arrived"
                              (str "step " (count trail) (when (and p (> (count trail) 50)
                                                                    (= (peek trail) (nth trail (- (count trail) 3))))
                                                           "  -- stuck in a local minimum")))))
                      "")]
                (set! (.-textContent caption)
                      (str method "  |  " label
                           "\nwhite: the robot  faint: obstacles grown by its radius  gold: the path  click: move the goal"))))
            (draw-arm! [{:keys [cspace plan]} px]
              (let [cw (.-width canvas) ch (.-height canvas)
                    t (:path-t @model)
                    k (when (seq plan) (mod (long (* 6 t)) (* 2 (count plan))))
                    angles (when k (plan (if (< k (count plan)) k (- (* 2 (count plan)) k 1))))]
                (doseq [o (:obstacles arm-world)] (set! (.-lineWidth ctx) 1) (poly! px o "#3a4458" "#5b6884"))
                (when angles
                  (doseq [[a b] (motion/arm-links arm angles)] (line! px a b "#ffffff" 6))
                  (dot! px (:base arm) 6 "#7dffa0"))
                ;; configuration space, cached as an image
                (when-not (identical? cspace (:drawn-cspace @model))
                  (set! (.-width cspace-img) arm-n) (set! (.-height cspace-img) arm-n)
                  (let [c2 (.getContext cspace-img "2d") img (.createImageData c2 arm-n arm-n) data (.-data img)]
                    (dotimes [i arm-n]
                      (dotimes [j arm-n]
                        (let [o (* 4 (+ i (* arm-n (- arm-n 1 j)))) hit (get-in cspace [i j])]
                          (aset data o (if hit 60 30)) (aset data (+ o 1) (if hit 70 45))
                          (aset data (+ o 2) (if hit 95 60)) (aset data (+ o 3) 255))))
                    (.putImageData c2 img 0 0))
                  (swap! model assoc :drawn-cspace cspace))
                ;; below the controls, in the right two fifths
                (let [side (min (- (* 0.4 cw) 40) (- ch 200))
                      x0 (+ (* 0.6 cw) (* 0.5 (- (* 0.4 cw) side))) y0 130
                      cpx (fn [[t1 t2]] [(+ x0 (* side (/ (+ t1 math/PI) (* 2 math/PI)))) (+ y0 (* side (- 1 (/ (+ t2 math/PI) (* 2 math/PI)))))])]
                  (set! (.-imageSmoothingEnabled ctx) false)
                  (.drawImage ctx cspace-img x0 y0 side side)
                  (doseq [[a b] (partition 2 1 plan)
                          :when (and (< (abs (- (a 0) (b 0))) 1) (< (abs (- (a 1) (b 1))) 1))]
                    (line! cpx a b "#ffd479" 2))
                  (when angles (dot! cpx angles 4 "#ffffff"))
                  (set! (.-fillStyle ctx) "#9fb4d8") (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
                  (.fillText ctx "configuration space: shoulder (across) by elbow (up)" x0 (- y0 8)))
                (set! (.-textContent caption)
                      (str "two-link arm  |  " (if (seq plan) (str (count plan) " cells on the plan") "no plan")
                           "\nright: every pair of joint angles, dark where the arm hits an obstacle; the plan is a path through the light"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (swap! model update :path-t + 0.05)
                (when (= "PRM" (.-method controls))
                  (when-let [p (:path (:st @model))] (swap! model assoc :last-path p)))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (when-not (= "two-link arm" (.-method controls))
                             (let [{:keys [scale x0 y0]} (layout)
                                   r (.getBoundingClientRect canvas)
                                   x (/ (- (- (.-clientX e) (.-left r)) x0) scale)
                                   y (- 8 (/ (- (- (.-clientY e) (.-top r)) y0) scale))]
                               (when (motion/free? world [x y])
                                 (reset! goal [x y])
                                 (restart!))))))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "method" #js ["visibility graph" "Voronoi graph" "cell decomposition" "PRM" "RRT" "RRT*"
                                                    "potential field" "two-link arm"])
                   restart!)
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
                   (when-let [container (js/document.getElementById "motion-planning")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
