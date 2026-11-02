(ns allgo.demo.contract-net
  "Canvas demo for `allgo.decision.contract-net`: the contract net of
  Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
  edition, chapter 18, a message at a time.

  A dispatcher (the square) has jobs to place with four vehicles, each
  starting from home. It announces each job to every vehicle; each
  replies with a bid -- how much the job would lengthen its route -- or
  declines if it lacks the skill (only some have cranes, marked C) or is
  full; the cheapest bid wins the award, and the rest are told they lost.
  A job with a part (the ringed ones) makes its winner a manager in turn,
  announcing the part itself. When every job is placed, recontracting:
  each vehicle replans its route, the best through its jobs, and then
  they offer jobs to one another -- one at a time, a job in exchange for
  one of the taker's (a swap), several together (a cluster), a job each
  among three or more vehicles (an M-contract), or any reallocation of a
  few jobs (an OCSM-contract) -- and a contract is made when it shortens
  the routes of those party to it together. Cost only ever falls. The
  offers are weighed a few milliseconds a frame, M- and OCSM-contracts
  being many. Then a new day."
  (:require [allgo.decision.contract-net :as cn]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]
            [clojure.string :as str]))

(def ^:private ^js controls #js {:jobs 12 :seed 1 :speed 4 :contracts "one, swap, cluster"})

(def ^:private contract-types
  {"one, swap, cluster" [:contracts #{:o :s :c}]
   "one task only" [:contracts #{:o}]
   "one, swap, cluster, M" [:contracts #{:o :s :c :m}]
   "OCSM, up to 3 jobs" [:contracts #{:ocsm} :max-ocsm 3]})

(def ^:private kind {:o "" :s ", a swap" :c ", a cluster" :m ", an M-contract" :ocsm ", an OCSM-contract"})

(def ^:private colors {:a "#7dc4ff" :b "#ffaa50" :c "#7dffa0" :d "#ff8fc8"})

(def ^:private fleet
  [{:id :a :home [1.0 1.0] :skills #{:haul} :capacity 6}
   {:id :b :home [11.0 1.0] :skills #{:haul :crane} :capacity 6}
   {:id :c :home [11.0 7.0] :skills #{:haul} :capacity 6}
   {:id :d :home [1.0 7.0] :skills #{:haul :crane} :capacity 6}])

(def ^:private dispatcher-at [6.0 4.0])

(defn- day [seed n]
  (let [rng (random/rng seed)
        job (fn [id] {:id id :at [(+ 0.5 (* 11 (rng))) (+ 0.5 (* 7 (rng)))]
                      :needs (if (< (rng) 0.25) :crane :haul)})]
    (vec (for [i (range 1 (inc n))]
           (let [j (job i)]
             (if (< (rng) 0.15)
               (assoc j :subtasks [(assoc (job (+ 100 i)) :needs :crane)])
               j))))))

(defn- run
  "The day's protocol, every message, and its recontracting, still to be
  walked: of that, only the contracts made are shown, appended to
  `:states` as they are found."
  [seed n]
  (let [protocol (vec (cn/contract-net-seq fleet (day seed n) :dispatcher))]
    {:states (mapv #(assoc % :phase :protocol) protocol)
     :re (rest (apply cn/recontract-seq (peek protocol) (contract-types (.-contracts controls))))
     :offers 0
     :first-cost (cn/total-cost (peek protocol))}))

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
    (letfn [(restart! [] (reset! model (assoc (run (.-seed controls) (.-jobs controls)) :i 0 :t 0.0 :last nil :day 0)))
            (advance! []
              (let [now (js/performance.now)
                    {:keys [last t i states]} @model
                    t (+ t (* (.-speed controls) 0.001 (- now (or last now))))]
                (swap! model assoc :last now)
                (cond (< t 1.0) (swap! model assoc :t t)
                      (< (inc i) (count states)) (swap! model assoc :t 0.0 :i (inc i))
                      (:re @model) (seek!)
                      (< i (+ (count states) (* 3 (.-speed controls)))) (swap! model assoc :t 0.0 :i (inc i))
                      ;; a pause at the end, then a new day
                      :else (let [d (inc (:day @model))]
                              (reset! model (assoc (run (+ (.-seed controls) d) (.-jobs controls)) :i 0 :t 0.0 :last now :day d))))))
            (seek! []
              ;; Weigh offers until a contract is made, or for 8 ms: a step weighs
              ;; a few hundred contracts at most.
              (let [deadline (+ (js/performance.now) 8)]
                (loop [re (:re @model) k (:offers @model) seen (:seen @model)]
                  (if-let [s (seq re)]
                    (let [st (first s) k (if (:concluded st) (inc k) k)]
                      (cond (or (:replanned st) (:moved st))
                            (swap! model #(-> % (assoc :re (rest s) :offers k :seen st)
                                              (update :states conj (assoc st :phase :recontract :offers k))))
                            (> (js/performance.now) deadline) (swap! model assoc :re (rest s) :offers k :seen st)
                            :else (recur (rest s) k st)))
                    (swap! model #(-> % (assoc :re nil :offers k)
                                      (update :states conj (assoc seen :phase :recontract :offers k :moved nil :replanned false))))))))
            (contract [{:keys [type moves]}]
              (str (str/join "; " (for [[[from to] ms] (sort (group-by val moves))
                                        :let [ts (sort (map key ms))]]
                                    (str (name from) " -> " (name to) ": " (if (next ts) "jobs " "job ") (str/join ", " ts))))
                   (kind type)))
            (text! [s x y color & [size align]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (set! (.-textAlign ctx) (or align "left")) (.fillText ctx s x y) (set! (.-textAlign ctx) "left"))
            (draw! []
              (let [{:keys [states i t first-cost re offers]} @model
                    st (nth states (min i (dec (count states))))
                    cw (.-width canvas) ch (.-height canvas)
                    scale (min (/ (- cw 60) 12.0) (/ (- ch 130) 8.0))
                    x0 (* 0.5 (- cw (* 12 scale))) y0 50
                    px (fn [[x y]] [(+ x0 (* scale x)) (+ y0 (* scale (- 8 y)))])
                    where (fn [who] (if (= who :dispatcher) dispatcher-at (:home ((:agents st) who))))
                    tasks (:tasks st)
                    msg (:delivered st)
                    moved (:moved st)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-fillStyle ctx) "#161b25")
                (let [[ax ay] (px [0 8])] (.fillRect ctx ax ay (* 12 scale) (* 8 scale)))
                ;; routes
                (doseq [[id route] (:routes st) :when (seq route)]
                  (set! (.-strokeStyle ctx) (colors id)) (set! (.-lineWidth ctx) 2)
                  (.beginPath ctx)
                  (let [[x y] (px (:home ((:agents st) id)))] (.moveTo ctx x y))
                  (doseq [tid route :let [[x y] (px (:at (tasks tid)))]] (.lineTo ctx x y))
                  (.stroke ctx))
                ;; jobs
                (doseq [[tid {:keys [at needs subtasks]}] tasks]
                  (let [[x y] (px at)
                        owner ((:awarded st) tid)
                        hot? (or (= tid (:task msg)) (contains? (:moves moved) tid))]
                    (set! (.-fillStyle ctx) (if owner (colors owner) "#c9d2e0"))
                    (.beginPath ctx) (.arc ctx x y (if hot? 7 5) 0 (* 2 math/PI)) (.fill ctx)
                    (when (seq subtasks)
                      (set! (.-strokeStyle ctx) "#c9d2e0") (set! (.-lineWidth ctx) 1.5)
                      (.beginPath ctx) (.arc ctx x y 10 0 (* 2 math/PI)) (.stroke ctx))
                    (when (= needs :crane) (text! "C" (+ x 8) (- y 6) "#ffd479" 11))))
                ;; dispatcher and vehicles
                (let [[x y] (px dispatcher-at)]
                  (set! (.-fillStyle ctx) "#e6edf7") (.fillRect ctx (- x 9) (- y 9) 18 18)
                  (text! "dispatcher" x (+ y 24) "#9fb4d8" 11 "center"))
                (doseq [{:keys [id home skills]} fleet]
                  (let [[x y] (px home)]
                    (set! (.-fillStyle ctx) (colors id))
                    (.beginPath ctx) (.moveTo ctx x (- y 11)) (.lineTo ctx (+ x 10) (+ y 8)) (.lineTo ctx (- x 10) (+ y 8)) (.fill ctx)
                    (text! (str (name id) (when (skills :crane) " C") "  " (count ((:routes st) id)) " jobs") x (+ y 24) (colors id) 11 "center")))
                ;; the message in flight
                (when (and msg (= :protocol (:phase st)))
                  (let [[ax ay] (px (where (:from msg))) [bx by] (px (where (:to msg)))
                        f (min 1.0 t)
                        [x y] [(+ ax (* f (- bx ax))) (+ ay (* f (- by ay)))]
                        color (case (:type msg) :announce "#e6edf7" :bid "#ffd479" :decline "#6b7894" :award "#7dffa0" :reject "#ff8f70")]
                    (set! (.-strokeStyle ctx) "rgba(255,255,255,0.15)") (set! (.-lineWidth ctx) 1)
                    (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx)
                    (set! (.-fillStyle ctx) color)
                    (.beginPath ctx) (.arc ctx x y 5 0 (* 2 math/PI)) (.fill ctx)
                    (text! (str (name (:type msg)) (when (:cost msg) (str " " (.toFixed (:cost msg) 2)))) (+ x 8) (- y 8) color 12)))
                (when moved
                  (let [{:keys [moves gain] :as m} moved [x y] (px (where (first (val (first moves)))))]
                    (text! (str (contract m) ", saves " (.toFixed gain 2) " m") x (- y 18) "#7dffa0" 12 "center")))
                (set! (.-textContent caption)
                      (str (if (= :recontract (:phase st))
                             (str "recontracting  |  " (cond (:replanned st) "each vehicle replans its route"
                                                             (:moved st) (contract (:moved st))
                                                             :else (str "nothing more to gain"))
                                  "  |  " (:offers st) " offers"
                                  (when (and re (>= (inc i) (count states))) (str ", weighing offer " (inc offers)))
                                  "  |  total route " (.toFixed (cn/total-cost st) 2) " m, from " (.toFixed first-cost 2) " after the first round")
                             (str "contract net  |  " (:messages st) " messages  |  " (count (:awarded st)) " of " (count tasks) " jobs placed"
                                  (when (seq (:unplaced st)) (str ", " (count (:unplaced st)) " unplaced"))
                                  "  |  total route " (.toFixed (cn/total-cost st) 2) " m"))
                           "\nwhite: announce  gold: bid  grey: decline  green: award  red: reject  |  C: needs, or has, a crane; ringed jobs have a crane part their winner subcontracts"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "jobs" 4 20 1) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 30 1)
        (.onChange (.add gui controls "contracts" (clj->js (keys contract-types))) restart!)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "contract-net")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
