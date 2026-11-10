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
  and can stall. Subdivision splits the mixed cells until a way opens;
  the bidirectional RRT grows a tree from each end until a milestone
  joins them; hybrid A* searches cells but drives a car's own arcs
  between them, so its path can be driven; trajectory optimization bends
  a straight line out of an obstacle by gradient descent. The two-link
  arm plans in its own configuration space,
  drawn beside it: dark where the arm would hit something.

  From LaValle's *Planning Algorithms* (Cambridge University Press,
  2006): the roadmap grown from the Halton sequence, which leaves no
  large gap for long; the visibility roadmap, a few guards and the
  connectors that join them; randomized potential fields, which walk at
  random when descent stalls; the exact vertical decomposition into
  trapezoids, and the boustrophedon sweep that covers every one; and a
  robot with momentum -- a double integrator -- planned on a lattice of
  its exact motions by A*, or by a kinodynamic RRT. Click to move the
  goal."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.random :as random]
            [allgo.robotics.decomposition :as decomposition]
            [allgo.robotics.kinodynamic :as kinodynamic]
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

;; The exact decomposition is of a point robot's free space.
(def ^:private point-world (assoc world :radius 0.0))

(def ^:private double-integrator {:amax 0.5 :vmax 1.0 :dt 1.0 :tolerance 0.4 :slow 0.3})

(defn- until-solved
  "The search states up to and including the one that solves it."
  [states]
  (let [[before after] (split-with (complement :solved?) states)]
    (concat before (take 1 after))))

(defn- arc-points
  "Points along the arcs of a double integrator's `states` under
  `controls`, four to an arc."
  [states controls dt]
  (into [(subvec (first states) 0 2)]
        (mapcat (fn [s u] (map #(subvec (kinodynamic/propagate s u (* dt (/ % 4))) 0 2) (range 1 5))) states controls)))

;; trajectory optimization starts from a straight line through an obstacle
(def ^:private traj-start [3.2 2.2])
(def ^:private traj-goal [8.0 2.8])

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
      "Halton PRM" (map-indexed (fn [k rm] {:roadmap rm :path (when (and (pos? k) (zero? (mod k 10)))
                                                                (motion/prm-path world rm start goal))})
                                (motion/halton-prm-seq world))
      "visibility PRM" (map-indexed (fn [k rm] (assoc rm :path (when (zero? (mod k 10)) (motion/visibility-roadmap-path world rm start goal))))
                                    (take 3000 (motion/visibility-prm-seq world rng)))
      "randomized potential field" (let [[before after] (split-with #(not= goal (:p %))
                                                                    (take 30000 (motion/randomized-potential-seq world start goal rng :influence 0.8 :repulse 0.1)))]
                                     (concat before (take 1 after)))
      "vertical decomposition" (let [cells (decomposition/vertical-decomposition point-world)
                                     doors (decomposition/adjacency cells)]
                                 (repeat {:cells cells :doors doors :path (decomposition/path cells doors start goal)}))
      "boustrophedon coverage" (let [cells (decomposition/vertical-decomposition point-world)
                                     {:keys [path]} (decomposition/boustrophedon point-world start 0.3)]
                                 (map (fn [k] {:cells cells :path (subvec path 0 k)}) (range 2 (inc (count path)))))
      "lattice planner" (map (fn [st] {:search st :plan (kinodynamic/lattice-result start double-integrator st)})
                             (until-solved (kinodynamic/lattice-seq world start goal double-integrator)))
      "kinodynamic RRT" (kinodynamic/rrt-seq world start goal rng double-integrator)
      "PRM" (map-indexed (fn [k rm] {:roadmap rm :path (when (and (pos? k) (zero? (mod k 10)))
                                                         (motion/prm-path world rm start goal))})
                         (motion/prm-seq world rng))
      "RRT" (motion/rrt-seq world start goal rng)
      "RRT*" (motion/rrt-seq world start goal rng :star? true)
      "potential field" (reductions conj [] (take 2000 (motion/potential-seq world start goal :influence 0.8 :repulse 0.1)))
      "subdivision" (motion/subdivision-seq world start goal 1.0 0.1)
      "bidirectional RRT" (motion/bidirectional-rrt-seq world start goal rng)
      "hybrid A*" (let [ss (motion/hybrid-a-star-seq world (conj start (/ math/PI 2)) goal {:min-radius 0.8 :step 0.6 :cell 0.4})
                        [before after] (split-with (comp nil? :path) ss)]
                    (concat (take 30000 before) (take 1 after)))
      "trajectory optimization" (take 400 (motion/trajectory-optimization-seq
                                           world (vec (for [i (range 25)] (polygon/lerp traj-start traj-goal (/ i 24.0))))
                                           {:epsilon 0.6 :lambda 0.02 :rate 0.05}))
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

                      "Halton PRM"
                      (let [{{:keys [nodes edges]} :roadmap} st
                            path (:last-path @model)]
                        (doseq [[i j] edges] (line! px (nodes i) (nodes j) "rgba(140,170,220,0.3)" 1))
                        (doseq [n nodes] (dot! px n 2 "#9fb4d8"))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count nodes) " Halton samples -- the same every time, each in the largest gap left -- " (count edges) " links"
                             (if path (str ", path " (.toFixed (motion/path-length path) 2) " m") ", not yet connected")))

                      "visibility PRM"
                      (let [{:keys [nodes edges guards]} st
                            path (:last-path @model)
                            guard? (set guards)]
                        (doseq [[i j] edges] (line! px (nodes i) (nodes j) "rgba(140,170,220,0.6)" 1.5))
                        (doseq [[i n] (map-indexed vector nodes)] (dot! px n (if (guard? i) 5 3) (if (guard? i) "#7fd1ff" "#ffd479")))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count guards) " guards (blue), each seeing what no other does, and " (- (count nodes) (count guards))
                             " connectors (gold) joining guards that see each other's samples"
                             (if path (str "; path " (.toFixed (motion/path-length path) 2) " m") "; not yet connected")))

                      "randomized potential field"
                      (let [{:keys [p walking trail]} st]
                        (doseq [[a b] (partition 2 1 trail)] (line! px a b "rgba(255,212,121,0.7)" 1.5))
                        (dot! px p (* scale (:radius world)) "#ffffff")
                        (cond (= p @goal) (str "arrived after " (count trail) " steps")
                              (pos? walking) (str "step " (count trail) " -- stalled, walking at random for " walking " more")
                              :else (str "step " (count trail) " -- descending the potential")))

                      "vertical decomposition"
                      (let [{:keys [cells doors path]} st]
                        (doseq [[k c] (map-indexed vector cells)]
                          (set! (.-lineWidth ctx) 1)
                          (poly! px (decomposition/corners c) (if (even? k) "rgba(120,200,140,0.26)" "rgba(120,170,255,0.26)") "rgba(180,200,230,0.7)"))
                        (doseq [d (vals doors)] (dot! px d 2.5 "#9fb4d8"))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count cells) " trapezoids from a vertical up and down from every corner, " (quot (count doors) 2) " doors between them (a point robot)"
                             (if path (str "; path " (.toFixed (motion/path-length path) 2) " m through the doors") "; the goal is unreachable")))

                      "boustrophedon coverage"
                      (let [{:keys [cells path]} st]
                        (doseq [[k c] (map-indexed vector cells)]
                          (set! (.-lineWidth ctx) 1)
                          (poly! px (decomposition/corners c) (if (even? k) "rgba(120,200,140,0.10)" "rgba(120,170,255,0.10)") "rgba(160,180,210,0.35)"))
                        (doseq [[a b] (partition 2 1 path)] (line! px a b "#ffd479" 2))
                        (dot! px (peek path) 5 "#ffffff")
                        (str "each trapezoid swept up and down in strips 0.3 m apart, the cells in depth-first order: "
                             (.toFixed (motion/path-length path) 1) " m so far"))

                      "lattice planner"
                      (let [{:keys [search plan]} st
                            seen (keys (:came-from search))]
                        (doseq [c seen] (dot! px (subvec (kinodynamic/lattice-state start double-integrator c) 0 2) 1.5 "rgba(140,170,220,0.6)"))
                        (when plan
                          (let [pts (arc-points (:states plan) (:controls plan) (:dt double-integrator))]
                            (path! px pts "#ffd479")
                            (robot-on! px scale pts path-t)))
                        (str (count seen) " lattice states reached -- position and velocity, each control full or no acceleration on each axis for a second"
                             (if plan (str "; at rest at the goal in " (:time plan) " s, the quickest the lattice holds") "; searching")))

                      "kinodynamic RRT"
                      (let [{:keys [nodes parent control] :as tree} st
                            plan (kinodynamic/rrt-plan tree)]
                        (doseq [[i p] (map-indexed vector parent) :when p]
                          (let [pts (arc-points [(nodes p)] [(control i)] (:dt double-integrator))]
                            (doseq [[a b] (partition 2 1 pts)] (line! px a b "rgba(140,170,220,0.4)" 1))))
                        (when plan
                          (let [pts (arc-points (:states plan) (:controls plan) (:dt double-integrator))]
                            (path! px pts "#ffd479")
                            (robot-on! px scale pts path-t)))
                        (str (count nodes) " states, each the best of nine controls from the nearest toward a random state"
                             (if plan (str "; slow at the goal after " (* (:dt double-integrator) (count (:controls plan))) " s") "; growing")))

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

                      "subdivision"
                      (let [{:keys [cells path]} st]
                        (doseq [[x y s kind] cells :let [[px0 py0] (px [x (+ y s)])]]
                          (set! (.-fillStyle ctx) (case kind :free "rgba(120,200,140,0.16)" :blocked "rgba(255,143,112,0.10)" "rgba(255,212,121,0.22)"))
                          (.fillRect ctx (+ px0 0.5) (+ py0 0.5) (- (* scale s) 1) (- (* scale s) 1)))
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count cells) " cells: green free, red blocked, yellow mixed, split in four each step"
                             (if path (str "; path " (.toFixed (motion/path-length path) 2) " m") "; no way through the free cells yet")))

                      "bidirectional RRT"
                      (let [{:keys [a b path]} st
                            tree! (fn [{:keys [nodes parent]} color] (doseq [[i p] (map-indexed vector parent) :when p] (line! px (nodes p) (nodes i) color 1)))]
                        (tree! a "rgba(125,255,160,0.5)") (tree! b "rgba(255,143,112,0.5)")
                        (when path (path! px path "#ffd479") (robot-on! px scale path path-t))
                        (str (count (:nodes a)) " + " (count (:nodes b)) " nodes, one tree from the start (green), one from the goal (red)"
                             (if path (str "; joined, path " (.toFixed (motion/path-length path) 2) " m") "; growing toward each milestone")))

                      "hybrid A*"
                      (let [{:keys [best path]} st]
                        (doseq [{:keys [pose]} (vals best)] (dot! px (subvec pose 0 2) 1.5 "rgba(140,170,220,0.5)"))
                        (when path
                          (doseq [[a b] (partition 2 1 path)
                                  :let [k (some (fn [k] (when (< (polygon/distance (subvec (motion/arc a k 0.6) 0 2) (subvec b 0 2)) 1e-9) k)) [1.25 0.0 -1.25])
                                        pts (map #(subvec (motion/arc a (or k 0.0) (* 0.6 (/ % 8.0))) 0 2) (range 9))]]
                            (doseq [[u v] (partition 2 1 pts)] (line! px u v "#ffd479" 3))))
                        (str (count best) " cells reached" (if path (str "; a path of " (dec (count path)) " arcs, none tighter than the car can turn") "; searching")))

                      "trajectory optimization"
                      (let [path st]
                        (path! px path "#ffd479")
                        (doseq [p path] (dot! px p 2.5 "#ffffff"))
                        (dot! px traj-start 5 "#7dffa0") (dot! px traj-goal 5 "#ff8f70")
                        (str "from the straight line, a gradient step at a time: " (if (motion/path-free? world path) "clear of the obstacle" "still through it")
                             "; the efficiency term pulls it straight, the obstacles' cost field pushes it out"))

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
                (when (#{"PRM" "Halton PRM" "visibility PRM"} (.-method controls))
                  (when-let [p (or (:path (:st @model)) (:path (:roadmap (:st @model))))] (swap! model assoc :last-path p)))
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
        (.onChange (.add gui controls "method" #js ["visibility graph" "Voronoi graph" "cell decomposition" "vertical decomposition"
                                                    "boustrophedon coverage" "PRM" "Halton PRM" "visibility PRM" "RRT" "RRT*"
                                                    "potential field" "randomized potential field" "subdivision" "bidirectional RRT"
                                                    "hybrid A*" "lattice planner" "kinodynamic RRT"
                                                    "trajectory optimization" "two-link arm"])
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
