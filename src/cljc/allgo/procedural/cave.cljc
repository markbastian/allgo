(ns allgo.procedural.cave
  (:require [allgo.random :as random]
            [clojure.math :as math]
            [clojure.set :refer [difference intersection]]))

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
  "One step in a random direction; with a grid, one that stays on it.
  `rng` is a function of no arguments returning [0, 1): `rand`, or an
  `allgo.random/rng` to get the same cave again."
  ([rng start] (mapv + start (random/pick rng [[0 1] [0 -1] [1 0] [-1 0]])))
  ([rng grid start]
   (first (filter (partial get-in grid) (repeatedly #(random-walk-step rng start))))))

(defn random-walk-cave-step [{:keys [current-location grid rng] :or {rng rand} :as m}]
  (let [next-location (random-walk-step rng grid current-location)]
    (-> m
        (assoc :current-location next-location)
        (assoc-in (into [:grid] next-location) :floor))))

(defn random-walk-cave-seq
  ([start grid] (random-walk-cave-seq start grid rand))
  ([start grid rng]
   (->> {:current-location start :grid (assoc-in grid start :floor) :rng rng}
        (iterate random-walk-cave-step)
        (map :grid)
        distinct)))

(defn wander-caverns
  ([start grid iterations] (wander-caverns start grid iterations rand))
  ([start grid iterations rng]
   (nth (random-walk-cave-seq start grid rng) iterations)))

;The frontier cave strategy randomly marks any "frontier" location of a grid as floor and
;then adds all wall neighbors of the selected location to the frontier. This tends to create
;larger more cavernous spaces than random walk as the space can expand from any location,
;not the location of the "walker"
(defn frontier-cave-step [{:keys [frontier grid rng] :or {rng rand} :as m}]
  (when-some [n (some->> frontier seq (random/pick rng))]
    (-> m
        (assoc-in (into [:grid] n) :floor)
        (update :frontier disj n)
        (update :frontier into (filter #(= :wall (get-in grid %)) (ortho-neighbors n))))))

(defn frontier-cave-seq
  ([start grid] (frontier-cave-seq start grid rand))
  ([start grid rng]
   (->> {:frontier #{start} :grid grid :rng rng}
        (iterate frontier-cave-step)
        (map :grid)
        (take-while identity))))

(defn frontier-caverns
  ([start grid iterations] (frontier-caverns start grid iterations rand))
  ([start grid iterations rng]
   (nth (frontier-cave-seq start grid rng) iterations)))

;Cellular automata caves
(defn ca-grid
  ([w h n] (ca-grid w h n rand))
  ([w h n rng]
   (random/sample rng n (for [row (range h) col (range w)] [row col]))))

(defn ca-cave-step [grid]
  (->> grid
       (mapcat all-neighbors)
       frequencies
       (filter (fn [[_ c]] (> c 4)))
       (map first)))

(def ca-cave-iterator #(iterate ca-cave-step %))

(defn ca-cave-seq
  ([w h pct] (ca-cave-seq w h pct rand))
  ([w h pct rng]
   (let [grid (init w h)]
     (->> (ca-grid w h pct rng)
          ca-cave-iterator
          (map (partial reduce mark-floor grid))))))

(defn ca-caverns
  ([w h pct iterations] (ca-caverns w h pct iterations rand))
  ([w h pct iterations rng]
   (nth (ca-cave-seq w h pct rng) iterations)))

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
   (fn [v] (math/round (double (/ v (count island)))))
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
  ([start finish] (shuffle-path-to rand start finish))
  ([rng start finish]
   (letfn [(signum [x] (cond (pos? x) 1 (neg? x) -1 :else 0))]
     (let [[dx dy] (map - finish start)
           x-steps (repeat (abs dx) [(signum dx) 0])
           y-steps (repeat (abs dy) [0 (signum dy)])]
       (->> (into x-steps y-steps)
            (random/shuffle rng)
            (reductions (partial mapv +) start))))))

(comment
  (shuffle-path-to [0 0] [10 10]))

(defn connect
  "Add connective :floor cells to each cavern to ensure all are connected."
  ([cavern-data] (connect cavern-data rand))
  ([cavern-data rng]
   (let [islands (-> cavern-data meadow-coords find-islands)
         centers (map center islands)
         links (take (dec (count islands)) (partition 2 1 (random/shuffle rng centers)))]
     (reduce
      (fn [acc coord] (assoc-in acc coord :floor))
      cavern-data
      (mapcat (fn [[a b]] (shuffle-path-to rng a b)) links)))))

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
