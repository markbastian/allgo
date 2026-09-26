(ns allgo.procedural.cave
  (:require [clojure.set :refer [difference intersection]]))

(defn grid-data->ascii-lines [grid]
  (mapv
   (fn [row] (apply str (map {:wall "#" :floor " "} row)))
   grid))

(defn prngrid [grid]
  (doseq [row (grid-data->ascii-lines grid)] (prn row)))

(defn init [w h] (vec (repeat h (vec (repeat w :wall)))))

(defn mark-floor [grid cell]
  (cond-> grid (get-in grid cell) (assoc-in cell :floor)))

(defn ortho-neighbors [[x y]]
  [[(inc x) y] [(dec x) y] [x (inc y)] [x (dec y)]])

(defn all-neighbors [[x y]]
  (map vector
       ((juxt identity inc inc identity dec dec dec identity inc) x)
       ((juxt identity identity inc inc inc identity dec dec dec) y)))

(defn random-walk-step
  ([start] (mapv + start (rand-nth [[0 1] [0 -1] [1 0] [-1 0]])))
  ([grid start]
   (first (filter (partial get-in grid) (repeatedly #(random-walk-step start))))))

(defn random-walk-cave-step [{:keys [current-location grid] :as m}]
  (let [next-location (random-walk-step grid current-location)]
    (-> m
        (assoc :current-location next-location)
        (assoc-in (into [:grid] next-location) :floor))))

(defn random-walk-cave-seq [start grid]
  (->> {:current-location start :grid (assoc-in grid start :floor)}
       (iterate random-walk-cave-step)
       (map :grid)
       distinct))

(defn wander-caverns [start grid iterations]
  (nth (random-walk-cave-seq start grid) iterations))

;The frontier cave strategy randomly marks any "frontier" location of a grid as floor and
;then adds all wall neighbors of the selected location to the frontier. This tends to create
;larger more cavernous spaces than random walk as the space can expand from any location,
;not the location of the "walker"
(defn frontier-cave-step [{:keys [frontier grid] :as m}]
  (when-some [n (some-> frontier seq rand-nth)]
    (-> m
        (assoc-in (into [:grid] n) :floor)
        (update :frontier disj n)
        (update :frontier into (filter #(= :wall (get-in grid %)) (ortho-neighbors n))))))

(defn frontier-cave-seq [start grid]
  (->> {:frontier #{start} :grid grid}
       (iterate frontier-cave-step)
       (map :grid)
       (take-while identity)))

(defn frontier-caverns [start grid iterations]
  (nth (frontier-cave-seq start grid) iterations))

;Cellular automata caves
(defn ca-grid [w h n]
  (random-sample n (for [row (range h) col (range w)] [row col])))

(defn ca-cave-step [grid]
  (->> grid
       (mapcat all-neighbors)
       frequencies
       (filter (fn [[_ c]] (> c 4)))
       (map first)))

(def ca-cave-iterator #(iterate ca-cave-step %))

(defn ca-cave-seq [w h pct]
  (let [grid (init w h)]
    (->> (ca-grid w h pct)
         ca-cave-iterator
         (map (partial reduce mark-floor grid)))))

(defn ca-caverns [w h pct iterations]
  (nth (ca-cave-seq w h pct) iterations))

(defn floor-coords [grid]
  (set (for [i (range (count grid)) j (range (count (grid i)))
             :when (= :floor (get-in grid [i j]))]
         [i j])))

(defn advance [{:keys [frontier unvisited] :as m}]
  (let [u (difference unvisited frontier)
        f (intersection u (set (mapcat ortho-neighbors frontier)))]
    (-> m
        (update :visited conj frontier)
        (assoc :unvisited u)
        (assoc :frontier f))))

(defn meadow-coords [grid]
  (for [i (range (count grid))
        j (range (count (grid i)))
        :let [cell (get-in grid [i j])]
        :when (or (= :floor cell) (= " " (str cell)))]
    [i j]))

(defn find-islands [[f & r]]
  (loop [frontier [f] unvisited (set r) visited #{} islands []]
    (if (first frontier)
      (let [front (filter unvisited (distinct (mapcat ortho-neighbors frontier)))
            u (difference unvisited (set frontier))
            v (into visited frontier)]
        (if (seq front)
          (recur front u v islands)
          (recur [(first u)] (disj u (first u)) (empty visited) (conj islands v))))
      islands)))

(defn center [island]
  (mapv
   (fn [v] (Math/round (double (/ v (count island)))))
   (apply mapv + island)))

(defn step-toward [a b]
  (letfn [(signum [x] (cond (pos? x) 1 (neg? x) -1 :else 0))]
    (let [[dx dy] (map - b a)
          delta (if (> (abs dx) (abs dy)) [(signum dx) 0] [0 (signum dy)])]
      (mapv + a delta))))

(defn path-to
  ([start finish]
   (->> (iterate #(step-toward % finish) start)
        (take-while (complement #{finish}))))
  ([[start finish]] (path-to start finish)))

(defn shuffle-path-to
  "Generate a path from start to finish (inclusive of both ends) that randomly
  shuffles steps, producing equivalent manhattan distances along the path."
  ([start finish]
   (letfn [(signum [x] (cond (pos? x) 1 (neg? x) -1 :else 0))]
     (let [[dx dy] (map - finish start)
           x-steps (repeat (abs dx) [(signum dx) 0])
           y-steps (repeat (abs dy) [0 (signum dy)])]
       (->> (into x-steps y-steps)
            shuffle
            (reductions (partial mapv +) start)))))
  ([[start finish]] (shuffle-path-to start finish)))

(comment
  (shuffle-path-to [0 0] [10 10])
  (shuffle-path-to [[0 0] [10 10]]))

(defn connect
  "Add connective :floor cells to each cavern to ensure all are connected."
  [cavern-data]
  (let [islands (-> cavern-data meadow-coords find-islands)
        centers (map center islands)
        links (take (dec (count islands)) (partition 2 1 (shuffle centers)))]
    (reduce
     (fn [acc coord] (assoc-in acc coord :floor))
     cavern-data
     (mapcat shuffle-path-to links))))

(defn connect-caverns [caverns]
  (->> caverns connect grid-data->ascii-lines))

(comment
  (let [cavern-data (ca-caverns 32 64 0.45 18)]
    (->> cavern-data
         connect
         grid-data->ascii-lines))

  (connect (ca-caverns 32 64 0.45 18))
  (connect-caverns (ca-caverns 32 64 0.45 18))
  (connect-caverns (wander-caverns [16 16] (init 32 32) 200))
  (connect-caverns (frontier-caverns [16 16] (init 32 32) 500)))
