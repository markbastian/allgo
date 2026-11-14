(ns allgo.robotics.coastal
  "Coastal navigation: planning so as to stay localized. Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  section 26.6 -- \"the coastal navigation heuristic which requires the
  robot to stay near known landmarks\" -- and the information gain in the
  cost that it generalizes to.

  The robot moves on a grid by dead reckoning, so its uncertainty about
  where it is grows a step at a time, until it comes within sensing range
  of a landmark and relocalizes. The planner searches over the robot's
  belief, not just its cell: a state is a cell and how uncertain the
  robot is there, `[cell spread]`, the spread growing by one a step to a
  cap and falling to zero near a landmark. Each step costs one, plus
  `:weight` times the spread it ends with -- uncertainty as a cost, the
  coastal heuristic's preference for knowing where one is. And a goal
  may demand precision: docking needs the spread at most `:dock`. The
  shortest path across open floor arrives lost; the planner instead
  hugs the coast of landmarks, or detours to one before docking.

  A map is `{:size [w h] :walls #{[x y]} :landmarks #{[x y]} :range r
  :cap c}`; moves are the four neighbors."
  (:require [allgo.search :as search]))

(defn- in-range? [{:keys [landmarks range]} [x y]]
  (some (fn [[lx ly]] (<= (+ (abs (- lx x)) (abs (- ly y))) range)) landmarks))

(defn spread-after
  "The robot's spread after a step into `cell` with spread `s`: zero in
  range of a landmark, else one more, to the cap."
  [{:keys [cap] :as m} cell s]
  (if (in-range? m cell) 0 (min cap (inc s))))

(defn- moves [{[w h] :size :keys [walls]} [x y]]
  (for [[dx dy] [[1 0] [-1 0] [0 1] [0 -1]]
        :let [c [(+ x dx) (+ y dy)]]
        :when (and (< -1 (first c) w) (< -1 (second c) h) (not (walls c)))]
    c))

(defn plan
  "The cheapest way from `start` to `goal` over beliefs, `[cell spread]`
  states from the start's: each step costing 1 plus `weight` times the
  spread it ends with, the goal reached only with spread at most `dock`
  (default the cap). nil if there is none."
  [m start goal & {:keys [weight dock] :or {weight 0.0}}]
  (let [dock (or dock (:cap m))
        s0 [start (spread-after m start (:cap m))]]
    (some->> (search/dijkstra-seq {:start s0
                                   :goal? (fn [[c s]] (and (= c goal) (<= s dock)))
                                   :neighbors (fn [[c s]] (for [c' (moves m c)] [c' (spread-after m c' s)]))
                                   :cost (fn [_ [_ s']] (+ 1.0 (* weight s')))})
             search/goal-state search/path)))

(defn shortest
  "The shortest way from `start` to `goal`, ignoring uncertainty, as
  `[cell spread]` states: what the robot believes along it."
  [m start goal]
  (plan m start goal :weight 0.0 :dock (:cap m)))
