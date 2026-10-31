(ns allgo.demo.search
  "Canvas demo for `allgo.search` and `allgo.search.online`: Russell and
  Norvig's chapter 3 searches, and section 4.5's online agents, on one
  maze, a step at a time.

  Every algorithm is a lazy seq of search states, and this walks it --
  `speed` steps a frame -- drawing each cell by how often the search has
  expanded it, the frontier waiting, and the path it is on or has found.
  Mud costs `mud` to enter, so breadth-first's fewest moves are not the
  cheapest; uniform cost spreads evenly, A* leans toward the goal, and
  weighted A* and greedy lean harder, for worse routes. The deepening
  searches and RBFS are tree searches in linear memory, and the shading
  shows their price: the same cells expanded again and again, which on a
  maze full of loops can go on longer than anyone will watch -- the
  caption counts it, and gives up at a million steps. Bidirectional
  search grows from both ends, green and blue. The online agents know
  nothing of the maze and must walk it: online DFS backtracks on foot;
  LRTA* learns which cells are dead ends and stops trying them.

  The caption keeps the memory each holds -- reached table and frontier,
  or stack -- beside the steps it took."
  (:require [allgo.search :as search]
            [allgo.search.online :as online]
            [allgo.search.problems :as problems]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:algorithm "A*"
       :size      "medium"
       :seed      3
       :mud       5
       :weight    2
       :beam      16
       :speed     8})

(def ^:private algorithms
  {"breadth-first"       search/breadth-first-seq
   "depth-first"         search/depth-first-seq
   "uniform cost"        search/dijkstra-seq
   "greedy best-first"   search/greedy-seq
   "A*"                  search/a-star-seq
   "weighted A*"         search/weighted-a-star-seq
   "beam"                search/beam-seq
   "iterative deepening" search/iterative-deepening-seq
   "IDA*"                search/ida-star-seq
   "recursive best-first" search/recursive-best-first-seq
   "bidirectional"       (fn [g] (search/bidirectional-seq (dissoc g :heuristic)))
   "bidirectional A*"    search/bidirectional-seq
   "online DFS agent"    :online-dfs
   "LRTA* agent"         :lrta})

(def ^:private order
  "The algorithms in the book's order, for the menu."
  ["breadth-first" "depth-first" "uniform cost" "greedy best-first" "A*" "weighted A*" "beam"
   "iterative deepening" "IDA*" "recursive best-first" "bidirectional" "bidirectional A*"
   "online DFS agent" "LRTA* agent"])

(def ^:private sizes {"small" [16 10] "medium" [32 20] "large" [56 34]})

(def ^:private give-up 1000000)

(def ^:private directions {:e [1 0] :n [0 1] :w [-1 0] :s [0 -1]})

(defn- maze-problem
  "The grid as an online agent meets it: four actions everywhere, a move
  into a wall or off the grid leaving it where it was."
  [{:keys [w h walls start goal cost heuristic]}]
  (let [open? (fn [[x y]] (and (< -1 x w) (< -1 y h) (not (walls [x y]))))]
    {:start start
     :actions (constantly [:e :n :w :s])
     :result (fn [[x y] a] (let [[dx dy] (directions a) c [(+ x dx) (+ y dy)]] (if (open? c) c [x y])))
     :goal? #(= % goal)
     :cost (fn [s _ s'] (if (= s s') 1 (cost s s')))
     :h (fn [c] (heuristic c goal))}))

(defn- node-path [n] (reverse (map :state (take-while some? (iterate :parent n)))))

(defn- came-from-path [table s]
  (when (contains? table s) (take-while some? (iterate table s))))

(defn- frontier-states [frontier]
  (cond
    (nil? frontier) nil
    (instance? search/PriorityFrontier frontier) (keys (:priority-of frontier))
    :else (map first frontier)))

(defn- view
  "What to draw of a search state: `{:frontier :frontier2 :path :memory}`."
  [st online?]
  (cond
    online? {:memory (count (or (:H (:memory st)) (:result (:memory st))))}
    (:fwd st) {:frontier (frontier-states (:frontier (:fwd st)))
               :frontier2 (frontier-states (:frontier (:bwd st)))
               :path (or (search/path st)
                         (when-let [c (:current st)]
                           (came-from-path (:came-from ((:direction st) st)) c)))
               :memory (+ (count (:came-from (:fwd st))) (count (:came-from (:bwd st))))}
    (contains? st :stack) (if (:node st)
                            {:frontier (map :state (:stack st))
                             :path (or (search/path st) (node-path (:node st)))
                             :memory (count (:stack st))}
                            ;; RBFS: the frames are the path; each holds
                            ;; its successors.
                            {:frontier (mapcat (fn [fr] (map :state (:succ fr))) (:stack st))
                             :path (or (search/path st) (map (comp :state :node) (:stack st)))
                             :memory (reduce + (map (fn [fr] (inc (count (:succ fr)))) (:stack st)))})
    :else {:frontier (frontier-states (:frontier st))
           :path (or (search/path st) (some->> (:current st) (came-from-path (:came-from st))))
           :memory (+ (count (:came-from st)) (count (frontier-states (:frontier st))))}))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        ;; The search: the grid, the rest of the seq, the state now, how
        ;; often each cell was expanded, and whether it is over.
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [[w h] (sizes (.-size controls))
                    g (problems/grid {:w w :h h :seed (.-seed controls) :mud (.-mud controls)})
                    g (assoc g :weight (.-weight controls) :beam-width (.-beam controls))
                    alg (algorithms (.-algorithm controls))
                    online? (keyword? alg)
                    states (if online?
                             (online/run-seq (if (= alg :lrta) online/lrta-star-agent online/online-dfs-agent)
                                             (maze-problem g) give-up)
                             (alg g))]
                (reset! model {:grid g :online? online? :rest states :st (first states)
                               :visits {} :steps 0 :done? false
                               :optimal (search/path-cost g (search/a-star g))})))
            (advance! []
              (let [{:keys [rest online? done?]} @model]
                (when-not done?
                  (loop [k 0 rest rest visits (:visits @model) steps (:steps @model)]
                    (let [st (first rest)
                          cell (if online? (:state st) (:current st))
                          visits (if cell (update visits cell (fnil inc 0)) visits)
                          finished? (or (nil? (next rest)) (:solved? st)
                                        (and online? (= :stop (:action st)))
                                        (>= steps give-up))]
                      (if (or finished? (>= k (dec (.-speed controls))))
                        (swap! model assoc :rest (if finished? rest (next rest)) :st st :visits visits
                               :steps (inc steps) :done? finished?)
                        (recur (inc k) (next rest) visits (inc steps))))))))
            (draw! []
              (let [{:keys [grid st visits steps done? online? optimal]} @model
                    {:keys [w h walls mud start goal]} grid
                    cw (.-width canvas) ch (- (.-height canvas) 40)
                    cell (min (/ cw w) (/ ch h))
                    x0 (* 0.5 (- cw (* cell w))) y0 (* 0.5 (- ch (* cell h)))
                    px (fn [[x y]] [(+ x0 (* cell x)) (+ y0 (* cell (- h 1 y)))])
                    {:keys [frontier frontier2 path memory]} (view st online?)
                    ;; At least 4, so that a cell expanded once is not
                    ;; drawn as hot as one expanded over and over.
                    most (max 4 (reduce max 0 (vals visits)))]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 (.-width canvas) (.-height canvas))
                (dotimes [x w]
                  (dotimes [y h]
                    (let [c [x y] [sx sy] (px c)
                          v (get visits c 0)]
                      (set! (.-fillStyle ctx)
                            (cond
                              (walls c) "#5a6274"
                              (pos? v) (let [t (/ (math/log (inc v)) (math/log (inc most)))]
                                         (str "rgb(" (int (+ 60 (* 180 t))) "," (int (+ 70 (* 60 t))) "," (int (- 120 (* 70 t))) ")"))
                              (mud c) "#4a3a26"
                              :else "#141922"))
                      (.fillRect ctx sx sy (- cell 1) (- cell 1)))))
                (doseq [[cells color] [[frontier "rgba(255,214,90,0.85)"] [frontier2 "rgba(110,200,255,0.85)"]]
                        c cells
                        :let [[sx sy] (px c)]]
                  (set! (.-fillStyle ctx) color)
                  (.fillRect ctx (+ sx (* 0.3 cell)) (+ sy (* 0.3 cell)) (* 0.4 cell) (* 0.4 cell)))
                (when (and (seq path) (not online?))
                  (set! (.-strokeStyle ctx) (if (:solved? st) "#7dffa0" "rgba(255,255,255,0.8)"))
                  (set! (.-lineWidth ctx) (max 2 (* 0.18 cell)))
                  (.beginPath ctx)
                  (doseq [[k c] (map-indexed vector path)
                          :let [[sx sy] (px c)]]
                    ((if (zero? k) #(.moveTo ctx %1 %2) #(.lineTo ctx %1 %2)) (+ sx (* 0.5 cell)) (+ sy (* 0.5 cell))))
                  (.stroke ctx))
                (doseq [[c color] [[start "#ffb050"] [goal "#7dffa0"] [(if online? (:state st) (:current st)) "#ffffff"]]
                        :when c
                        :let [[sx sy] (px c)]]
                  (set! (.-fillStyle ctx) color)
                  (.beginPath ctx)
                  (.arc ctx (+ sx (* 0.5 cell)) (+ sy (* 0.5 cell)) (* 0.32 cell) 0 (* 2 math/PI))
                  (.fill ctx))
                (set! (.-textContent caption)
                      (str (.-algorithm controls) "  |  step " steps
                           (when-let [l (:limit st)] (str "  |  bound " l))
                           "  |  in memory " memory
                           "  |  cells expanded " (count visits)
                           (cond
                             (and online? done? (= :stop (:action st)))
                             (str "  |  arrived, " (:steps st) " moves costing " (:cost st)
                                  " (best route " optimal ")")
                             online? (str "  |  " (:steps st) " moves so far")
                             (:solved? st) (str "  |  path cost " (search/path-cost grid (search/path st))
                                                " (best " optimal ")")
                             (>= steps give-up) "  |  gave up"
                             done? "  |  no route found"
                             :else "")
                           "\nshade: times expanded   yellow: frontier" (when frontier2 "   blue: backward frontier")
                           "   brown: mud"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "algorithm" (clj->js order)) restart!)
        (.onChange (.add gui controls "size" #js ["small" "medium" "large"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.onChange (.add gui controls "mud" 1 9 1) restart!)
        (.onChange (.add gui controls "weight" 1 5 0.1) restart!)
        (.onChange (.add gui controls "beam" 1 64 1) restart!)
        (.add gui controls "speed" 1 400 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "search")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
