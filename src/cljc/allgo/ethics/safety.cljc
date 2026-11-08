(ns allgo.ethics.safety
  "AI safety: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, section 27.3.7.

    fault tree     safety engineering's AND/OR tree of how a system can
                   fail, a probability for each root cause: the chance of
                   the top failure, the causes independent, and the
                   minimal cut sets -- the smallest sets of causes that
                   bring it about
    FMEA           failure modes and effects analysis, each mode ranked by
                   severity times occurrence times detection
    low impact     an agent maximizing its utility less a weighted count
                   of the changes it makes to the world (after Armstrong
                   and Levinstein, \"Low Impact Artificial Intelligences\",
                   2017): told only to fetch the coffee, it goes round the
                   vase rather than through it -- not because it knows
                   vases matter, but because it prefers to change nothing
                   it need not

  A fault tree is `{:and [...]}`, `{:or [...]}` or `{:event name}`."
  (:require [allgo.decision.mdp :as mdp]))

;; ---------------------------------------------------------------------------
;; Fault trees

(defn failure-probability
  "The probability of the tree's top event, its root causes independent
  with the probabilities `p` (event to probability): an AND the product of
  its children's, an OR one less the product of their complements."
  [tree p]
  (cond (:event tree) (p (:event tree))
        (:and tree) (reduce * (map #(failure-probability % p) (:and tree)))
        (:or tree) (- 1.0 (reduce * (map #(- 1.0 (failure-probability % p)) (:or tree))))))

(defn cut-sets
  "The minimal cut sets of the tree: an OR's are its children's together,
  an AND's every union of one from each child, and any set holding
  another left out. A set of sets of events."
  [tree]
  (let [sets (cond (:event tree) #{#{(:event tree)}}
                   (:or tree) (set (mapcat cut-sets (:or tree)))
                   (:and tree) (reduce (fn [acc child] (set (for [a acc b (cut-sets child)] (into a b))))
                                       #{#{}} (:and tree)))]
    (set (remove (fn [s] (some #(and (not= % s) (every? s %)) sets)) sets))))

(defn fmea
  "Failure modes, each `{:mode :severity :occurrence :detection}` rated 1
  to 10 (detection 10 the hardest to catch), ranked by their risk
  priority number, the product, worst first."
  [modes]
  (sort-by (comp - :rpn) (map #(assoc % :rpn (* (:severity %) (:occurrence %) (:detection %))) modes)))

;; ---------------------------------------------------------------------------
;; A low-impact agent

(defn coffee-world
  "A corridor world for the coffee fetcher: `rows` strings, `.` floor, `#`
  wall, `V` a vase, `S` the start and `C` the coffee. The state is the
  agent's cell and the set of vases broken; walking onto a vase breaks
  it. The designer's reward: `step` each move and 1 for the coffee --
  nothing about vases. With `impact` above 0, each broken vase also costs
  that much a step, the low-impact agent's penalty for changing the
  world. Deterministic moves; an MDP as `allgo.decision.mdp` takes it."
  [rows step impact gamma]
  (let [h (count rows)
        cell (fn [ch] (for [[r row] (map-indexed vector rows) [x c] (map-indexed vector row) :when (= c ch)] [x (- h 1 r)]))
        floor (set (for [[r row] (map-indexed vector rows) [x c] (map-indexed vector row) :when (not= c \#)] [x (- h 1 r)]))
        vases (set (cell \V))
        start (first (cell \S)) coffee (first (cell \C))
        subsets (reduce (fn [acc v] (into acc (map #(conj % v) acc))) #{#{}} vases)
        states (vec (for [c floor b subsets] [c b]))
        moves [[1 0] [-1 0] [0 1] [0 -1]]]
    {:states states :start [start #{}] :gamma gamma
     :terminals (set (filter #(= coffee (first %)) states))
     :actions (fn [[c]] (if (= c coffee) [] moves))
     :transitions (fn [[[x y] broken :as s] a]
                    (if (nil? a)
                      [[1.0 s]]
                      (let [[dx dy] a c [(+ x dx) (+ y dy)]]
                        (if (floor c)
                          [[1.0 [c (if (vases c) (conj broken c) broken)]]]
                          [[1.0 s]]))))
     :reward (into {} (for [[c broken :as s] states]
                        [s (- (if (= c coffee) 1.0 step) (* impact (count broken)))]))}))

(defn route
  "The cells the agent visits acting on the optimal policy of `world`
  from its start, and the vases broken: `{:cells :broken}`."
  [world]
  (let [U (mdp/value-iteration world 1e-6)
        pi (mdp/best-policy world U)]
    (loop [s (:start world) cells [(first (:start world))] n 0]
      (if-let [a (and (< n 100) (pi s))]
        (let [[[_ s']] ((:transitions world) s a)]
          (recur s' (conj cells (first s')) (inc n)))
        {:cells cells :broken (second s)}))))
