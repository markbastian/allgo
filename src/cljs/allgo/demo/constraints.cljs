(ns allgo.demo.constraints
  "Canvas demo for `allgo.search.csp`: Russell and Norvig's constraint
  satisfaction (*Artificial Intelligence: A Modern Approach*, 4th
  edition, chapter 6) on a sudoku, a map to color and n-queens, a step at
  a time.

  Backtracking assigns one value a step; the small digits in a sudoku
  cell, and the faint colors on an unassigned region, are the values its
  domain still holds, so forward checking and MAC can be seen pruning
  them -- and the caption counts the assignments each saves. MRV picks the
  variable with fewest values left; LCV tries the value that rules out
  least first. AC-3 alone revises one arc a step and stops when every arc
  is consistent: on an easy sudoku it fills much of the grid without a
  single guess. Min-conflicts starts from a complete assignment and moves
  one conflicted variable at a time to its least conflicted value -- the
  conflicted ones drawn red. Cutset conditioning, on the map, fixes the
  outlined regions whose removal leaves a forest, and solves the forest
  with the tree solver, which never backtracks.

  The four sudokus were made by `allgo.search.problems/sudoku-puzzle`,
  seeds 11 to 14, each with only one solution."
  (:require [allgo.random :as random]
            [allgo.search.csp :as csp]
            [allgo.search.problems :as problems]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private puzzles
  {"sudoku, 45 blanks" "030470026000000000074056010053201000802000050917000000008605730040007290720304681"
   "sudoku, 52 blanks" "000050001010000050070149000040008170061094800002070000050086000000010900030900768"
   "sudoku, 54 blanks" "107508000020030007003007000000050004000004600000700513001005000030200400605001920"
   "sudoku, 56 blanks" "038760000100040000090000020080900105000000080610430700000000579000000000309002060"})

(def ^:private ^js controls
  #js {:problem   "sudoku, 52 blanks"
       :algorithm "backtracking"
       :select    "MRV"
       :order     "as given"
       :inference "forward checking"
       :regions   30
       :colors    4
       :queens    12
       :seed      2
       :speed     4})

(def ^:private problem-names
  (into (vec (sort (keys puzzles))) ["map coloring" "n-queens"]))

(def ^:private palette ["#e8453c" "#ffcc33" "#4aa3df" "#6cc070" "#b37be0" "#f08a3c"])

(defn- problem []
  (let [name (.-problem controls)]
    (cond
      (puzzles name) {:kind :sudoku
                      :grid (mapv #(- (.charCodeAt % 0) 48) (puzzles name))
                      :csp (problems/sudoku (mapv #(- (.charCodeAt % 0) 48) (puzzles name)))}
      (= name "map coloring") (let [{:keys [mesh csp]} (problems/region-map (.-regions controls) (.-colors controls) (.-seed controls))]
                                {:kind :map :mesh mesh :csp csp})
      :else {:kind :queens :csp (problems/queens-csp (.-queens controls))})))

(defn- states [p]
  (case (.-algorithm controls)
    "backtracking" (csp/backtracking-seq p {:select (if (= "MRV" (.-select controls)) :mrv :first)
                                            :order (if (= "LCV" (.-order controls)) :lcv :given)
                                            :inference (case (.-inference controls)
                                                         "forward checking" :forward-checking
                                                         "MAC" :mac
                                                         :none)})
    "AC-3 alone" (csp/ac-3-seq p)
    "min-conflicts" (csp/min-conflicts-seq p {:rng (random/rng (.-seed controls))})
    "cutset conditioning" (csp/cutset-conditioning-seq p)))

(defn- view
  "What a state shows: `{:assignment :domains :done? :solved? :cutset
  :focus}`."
  [st]
  (cond
    (contains? st :stack)                       ; backtracking
    (let [top (peek (:stack st))
          v (:value st)]
      {:assignment (if (map? v) v (or (:assignment top) {}))
       :domains (or (:domains top) {})
       :done? (empty? (:stack st))
       :solved? (map? v)
       :focus (:var top)})
    (contains? st :queue)                       ; AC-3
    {:assignment (into {} (for [[k d] (:domains st) :when (= 1 (count d))] [k (first d)]))
     :domains (:domains st)
     :done? (or (:failed? st) (empty? (:queue st)))
     :focus (first (:arc st))}
    (contains? st :current)                     ; min-conflicts
    {:assignment (:current st) :focus (:var st) :min-conflicts? true}
    :else                                       ; cutset conditioning
    {:assignment (:assignment st) :domains (:domains st)
     :done? (or (:solved? st) (:failed? st)) :solved? (:solved? st)
     :cutset (set (keys (:cutset-assignment st))) :focus (:var st)}))

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
              ;; A sudoku or a board of queens has a cutset of most of its
              ;; variables, and finding even one assignment to it can take
              ;; longer than a frame: cutset conditioning is for the map.
              (when (and (= "cutset conditioning" (.-algorithm controls))
                         (not= "map coloring" (.-problem controls)))
                (set! (.-algorithm controls) "backtracking"))
              (let [p (problem)
                    ss (states (:csp p))]
                (reset! model (assoc p :rest ss :st (first ss) :steps 0 :done? false))))
            (advance! []
              (let [{:keys [rest done?]} @model]
                (when-not done?
                  (loop [k 0 rest rest steps (:steps @model)]
                    (let [st (first rest)
                          finished? (or (nil? (next rest)) (:done? (view st)))]
                      (if (or finished? (>= k (dec (.-speed controls))))
                        (swap! model assoc :rest (if finished? rest (next rest)) :st st :steps steps :done? finished?)
                        (recur (inc k) (next rest) (inc steps))))))))
            (conflicted [p assignment]
              (set (filter #(and (contains? assignment %) (pos? (csp/conflicts p assignment % (assignment %))))
                           (:variables p))))
            (draw-sudoku! [{:keys [grid]} {:keys [assignment domains focus]} bad x0 y0 side]
              (let [cell (/ side 9)]
                (dotimes [k 81]
                  (let [r (quot k 9) c (mod k 9) x (+ x0 (* c cell)) y (+ y0 (* r cell))
                        v (assignment k) d (domains k)]
                    (set! (.-fillStyle ctx) (cond (= k focus) "#3a3420" (bad k) "#4a2020" :else "#161b25"))
                    (.fillRect ctx x y (- cell 1) (- cell 1))
                    (set! (.-textAlign ctx) "center") (set! (.-textBaseline ctx) "middle")
                    (cond
                      (pos? (grid k))
                      (do (set! (.-fillStyle ctx) "#ffffff")
                          (set! (.-font ctx) (str "bold " (int (* 0.55 cell)) "px sans-serif"))
                          (.fillText ctx (str (grid k)) (+ x (* 0.5 cell)) (+ y (* 0.52 cell))))
                      v
                      (do (set! (.-fillStyle ctx) "#7dd3ff")
                          (set! (.-font ctx) (str (int (* 0.55 cell)) "px sans-serif"))
                          (.fillText ctx (str v) (+ x (* 0.5 cell)) (+ y (* 0.52 cell))))
                      d
                      (do (set! (.-fillStyle ctx) "#6b7488")
                          (set! (.-font ctx) (str (int (* 0.22 cell)) "px sans-serif"))
                          (doseq [n d]
                            (.fillText ctx (str n) (+ x (* cell (+ 0.2 (* 0.3 (mod (dec n) 3)))))
                                       (+ y (* cell (+ 0.2 (* 0.3 (quot (dec n) 3)))))))))))
                (set! (.-strokeStyle ctx) "#c9d2e0")
                (set! (.-lineWidth ctx) 2)
                (doseq [i (range 4)]
                  (.beginPath ctx) (.moveTo ctx (+ x0 (* 3 i cell)) y0) (.lineTo ctx (+ x0 (* 3 i cell)) (+ y0 side)) (.stroke ctx)
                  (.beginPath ctx) (.moveTo ctx x0 (+ y0 (* 3 i cell))) (.lineTo ctx (+ x0 side) (+ y0 (* 3 i cell))) (.stroke ctx))))
            (draw-map! [{:keys [mesh]} {:keys [assignment domains cutset focus]} bad x0 y0 side]
              ;; Cells on the rim reach past the map: cut them at its edge.
              (.save ctx)
              (.beginPath ctx) (.rect ctx x0 y0 side side) (.clip ctx)
              (doseq [c (:centers mesh)
                      :let [ring (mapv #(:point ((:corners mesh) %)) (:corners c))
                            id (:id c)
                            v (assignment id)]
                      :when (>= (count ring) 3)]
                (.beginPath ctx)
                (doseq [[k [x y]] (map-indexed vector ring)]
                  ((if (zero? k) #(.moveTo ctx %1 %2) #(.lineTo ctx %1 %2)) (+ x0 (* side x)) (+ y0 (* side (- 1 y)))))
                (.closePath ctx)
                (set! (.-fillStyle ctx) (if v (palette v) "#1d2330"))
                (.fill ctx)
                ;; Unassigned: a dot of each color still possible.
                (when (and (nil? v) domains (domains id))
                  (let [[px py] (:point c) ds (domains id)]
                    (doseq [[k col] (map-indexed vector ds)]
                      (set! (.-fillStyle ctx) (palette col))
                      (.beginPath ctx)
                      (.arc ctx (+ x0 (* side px) (* 7 (- k (* 0.5 (dec (count ds)))))) (+ y0 (* side (- 1 py))) 2.5 0 (* 2 math/PI))
                      (.fill ctx))))
                (set! (.-strokeStyle ctx) (cond (bad id) "#ff3030" (= id focus) "#ffffff" (and cutset (cutset id)) "#ffffff" :else "#0d1118"))
                (set! (.-lineWidth ctx) (if (or (bad id) (= id focus) (and cutset (cutset id))) 3 1))
                (.stroke ctx))
              (.restore ctx))
            (draw-queens! [{:keys [csp]} {:keys [assignment focus]} bad x0 y0 side]
              (let [n (count (:variables csp)) cell (/ side n)]
                (dotimes [c n]
                  (dotimes [r n]
                    (set! (.-fillStyle ctx) (if (even? (+ c r)) "#2a3140" "#202633"))
                    (.fillRect ctx (+ x0 (* c cell)) (+ y0 (* cell (- n 1 r))) cell cell)))
                (doseq [[c r] assignment]
                  (set! (.-fillStyle ctx) (cond (bad c) "#ff5050" (= c focus) "#ffffff" :else "#ffd479"))
                  (.beginPath ctx)
                  (.arc ctx (+ x0 (* cell (+ c 0.5))) (+ y0 (* cell (- n r 0.5))) (* 0.35 cell) 0 (* 2 math/PI))
                  (.fill ctx))))
            (draw! []
              (let [{:keys [kind csp st steps done?] :as m} @model
                    v (view st)
                    bad (conflicted csp (:assignment v))
                    cw (.-width canvas) ch (.-height canvas)
                    side (min (- cw 40) (- ch 134))
                    x0 (* 0.5 (- cw side)) y0 44]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (case kind
                  :sudoku (draw-sudoku! m v bad x0 y0 side)
                  :map (draw-map! m v bad x0 y0 side)
                  :queens (draw-queens! m v bad x0 y0 side))
                (set! (.-textContent caption)
                      (let [c (:context st)]
                        (str (.-algorithm controls)
                             (when (= "backtracking" (.-algorithm controls))
                               (str " (" (.-select controls) ", " (.-order controls) ", " (.-inference controls) ")"))
                             "  |  step " steps
                             (when c (str "  |  " (:assignments c) " assignments, " (:backtracks c) " backtracks"))
                             (when-let [k (:checks st)] (str "  |  " k " constraint checks"))
                             (when (:cutset v) (str "  |  cutset of " (count (:cutset v))))
                             "  |  " (count (:assignment v)) " of " (count (:variables csp)) " assigned"
                             (cond
                               (and done? (or (:solved? v) (and (:min-conflicts? v) (empty? bad)
                                                                (= (count (:assignment v)) (count (:variables csp))))))
                               "  |  solved"
                               (and done? (= "AC-3 alone" (.-algorithm controls)))
                               (if (:failed? st) "  |  a domain emptied: no solution" "  |  arc consistent")
                               done? "  |  no solution"
                               :else "")
                             "\n" (case kind
                                    :sudoku "white: given  blue: assigned  small digits: values still in the domain"
                                    :map "dots: colors still in the domain  outlined white: the cutset, or the variable just set"
                                    :queens "")
                             "  red: in conflict")))))
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
        (.onChange (.add gui controls "problem" (clj->js problem-names)) restart!)
        (.onChange (.listen (.add gui controls "algorithm" #js ["backtracking" "AC-3 alone" "min-conflicts" "cutset conditioning"])) restart!)
        (.onChange (.add gui controls "select" #js ["first" "MRV"]) restart!)
        (.onChange (.add gui controls "order" #js ["as given" "LCV"]) restart!)
        (.onChange (.add gui controls "inference" #js ["none" "forward checking" "MAC"]) restart!)
        (.onChange (.add gui controls "regions" 10 200 1) restart!)
        (.onChange (.add gui controls "colors" 2 6 1) restart!)
        (.onChange (.add gui controls "queens" 4 60 1) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 500 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "constraints")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
