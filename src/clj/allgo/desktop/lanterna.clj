(ns allgo.desktop.lanterna
  (:require [allgo.procedural.cave :as cave]
            [clojure.core.async :as ca]
            [lanterna.screen :as s]))

(defonce scr (s/get-screen))

;; TODO - Connect before other steps
;(def cave-data (cave/connect (cave/ca-caverns 64 32 0.45 50)))
;(def floor-cells (cave/floor-coords cave-data))
;(def start (rand-nth (vec floor-cells)))
;(def cave (cave/grid-data->ascii-lines cave-data))

(defn move-to [{:keys [cave-data] :as state} coord]
  (cond-> state
    (= :floor (get-in cave-data coord))
    (assoc :coord coord)))

(defn move-up [{[row col] :coord :as state}]
  (move-to state [(dec row) col]))

(defn move-down [{[row col] :coord :as state}]
  (move-to state [(inc row) col]))

(defn move-left [{[row col] :coord :as state}]
  (move-to state [row (dec col)]))

(defn move-right [{[row col] :coord :as state}]
  (move-to state [row (inc col)]))

(defn move-agent [{:keys [coord cave-data]}
                  {:keys [location tick tick-offset] :as agent}
                  global-tick]
  (if (zero? (rem global-tick (+ tick tick-offset)))
    (let [dirs [[0 0] [0 1] [0 -1] [1 0] [-1 0]]
          dir  (rand-nth dirs)
          nloc (mapv + location dir)]
      (cond-> agent
        (and
         (not= coord nloc)
         (= :floor (get-in cave-data nloc)))
        (assoc :location nloc)))
    agent))

(defn move-agents [state global-tick]
  (update state :agents (fn [agents] (set (map #(move-agent state % global-tick) agents)))))

;; Do this instead....
;; - Put everything in the cave-data grid
;; - Note that we can place items on the floor
;; - When they move, we replace the space with :floor
;; - Maybe someday: Each cell is a stack.
;;   - You move out and we just pop
;;   - Moving in is a push (if the peek is valid)
;; - coord, rat, etc. are all agents. we track their cell location along with
;;   anything else that is needed.
(defn initial-state []
  (let [cave-data   (cave/connect (cave/ca-caverns 64 32 0.45 50))
        floor-cells (cave/floor-coords cave-data)
        start       (rand-nth (vec floor-cells))
        cave        (cave/grid-data->ascii-lines cave-data)]
    {:cave      cave
     :cave-data cave-data
     :coord     start
     :agents    (set (repeatedly 5 (fn [] {:id          (random-uuid)
                                           :location    (rand-nth (vec (disj floor-cells start)))
                                           :tick        5
                                           :tick-offset (rand-int 5)
                                           :hp          2})))}))

(defonce state (atom (initial-state)))

(comment
  (reset! state (initial-state)))

(defn init [{:keys [coord agents cave]}]
  (let [[row col] coord]
    (s/clear scr)
    (doseq [row-index (range (count cave))]
      (s/put-string scr 0 row-index (cave row-index)))
    (doseq [{[col row] :location} agents]
      (s/put-string scr col row "R"))
    (s/move-cursor scr col row)
    (s/redraw scr)))

(defn add-move-watch [state]
  (add-watch
   state
   :move-watch
   (fn [_context _key {o :coord} {[row col :as n] :coord}]
     (when-not (= o n)
       (s/move-cursor scr col row)
       (s/redraw scr)))))

(defn add-agents-watch [state]
  (add-watch
   state
   :agents
   (fn [_context _key {agents :agents} {agents' :agents}]
     (when-not (= agents agents')
       (doseq [{:keys [location]} agents
               :let [[row col] location]]
         (s/put-string scr col row " "))
       (doseq [{:keys [location]} agents'
               :let [[row col] location]]
         (s/put-string scr col row "R"))
       (s/redraw scr)))))

(defn handle-input! [state]
  (ca/go-loop [key-pressed (s/get-key scr) tick 0]
    (swap! state move-agents tick)
    (if key-pressed
      (let [res (case key-pressed
                  :up (swap! state move-up)
                  :down (swap! state move-down)
                  :right (swap! state move-right)
                  :left (swap! state move-left)
                  :escape :kill
                  (println key-pressed))]
        (if-not (= :kill res)
          (do
            (ca/<! (ca/timeout 100))
            (recur (s/get-key scr) (inc tick)))
          (s/stop scr)))
      (do
        (ca/<! (ca/timeout 100))
        (recur (s/get-key scr) (inc tick))))))

(comment
  (do
    (remove-watch state :move-watch)
    (remove-watch state :agents)
    (s/start scr)
    (reset! state (initial-state))
    (init @state)
    (add-move-watch state)
    (add-agents-watch state)
    (handle-input! state)))

(comment
  (s/stop scr))
