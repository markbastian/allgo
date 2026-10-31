(ns allgo.demo.local-search
  "Canvas demo for `allgo.search.local`: Russell and Norvig's local search (chapter 4)
  on n-queens and on a traveling salesperson's tour, a move a step.

  On the queens' board, red lines join the queens that attack each other;
  the value is minus their number, and zero is solved. Steepest ascent
  stops on the first peak -- on eight queens about six boards in seven
  are not solutions -- and `trials` measures it: 200 boards, the share
  solved and the moves taken. Sideways moves cross the plateaus, random
  restarts keep trying, simulated annealing wanders downhill early and
  settles late, local beam passes its effort to wherever climbing is
  good, and the genetic algorithm breeds a solution from strings of
  rows. On the tour, a move reverses a stretch of it, which uncrosses two
  edges; the chart below plots the value as it goes."
  (:require [allgo.random :as random]
            [allgo.search.local :as local]
            [allgo.search.problems :as problems]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:problem   "8-queens"
       :algorithm "steepest ascent"
       :queens    8
       :cities    40
       :seed      1
       :speed     1})

(def ^:private names
  ["steepest ascent" "sideways moves" "stochastic" "first-choice" "random restart"
   "simulated annealing" "local beam" "genetic algorithm"])

(defn- tsp? [] (= "traveling salesman" (.-problem controls)))

(defn- problem [seed]
  (let [rng (random/rng seed)]
    (if (tsp?)
      (assoc (problems/tsp (.-cities controls) rng) :rng rng :goal? (constantly false))
      (assoc (problems/queens (.-queens controls) rng) :rng rng))))

(defn- states
  "The seq of search states of the chosen algorithm on `p`."
  [p]
  (let [n (.-queens controls)]
    (case (.-algorithm controls)
      "steepest ascent" (local/hill-climbing-seq p)
      "sideways moves" (local/hill-climbing-seq (assoc p :sideways 100))
      "stochastic" (local/hill-climbing-seq (assoc p :variant :stochastic))
      "first-choice" (local/hill-climbing-seq (assoc p :variant :first-choice))
      "random restart" (local/random-restart-seq local/hill-climbing-seq (assoc p :restarts (if (tsp?) 20 1000)))
      "simulated annealing" (local/simulated-annealing-seq
                             (cond-> p (tsp?) (assoc :schedule (local/exponential-schedule 0.3 0.9998 1e-4))))
      "local beam" (local/local-beam-seq (assoc p :k (if (tsp?) 4 8)))
      "genetic algorithm" (let [g (problems/queens-genetic n)]
                            (local/genetic-seq (merge g {:rng (:rng p)
                                                         :population (local/random-population (:rng p) 100 n (:alphabet g))}))))))

(defn- solved? [st]
  (and (not (tsp?))
       (zero? (problems/attacking-pairs (:current st)))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        trials (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (when (and (tsp?) (= "genetic algorithm" (.-algorithm controls)))
                ;; One-point crossover of two tours is not a tour.
                (set! (.-algorithm controls) "simulated annealing"))
              ;; A tour takes tens of thousands of moves; queens a handful.
              (when (and (tsp?) (< (.-speed controls) 100)) (set! (.-speed controls) 200))
              (let [p (problem (.-seed controls))
                    ss (states p)]
                (reset! trials nil)
                (reset! model {:problem p :rest ss :st (first ss) :history [(:current-value (first ss))]
                               :steps 0 :done? false})))
            (advance! []
              (let [{:keys [rest done?]} @model]
                (when-not done?
                  (loop [k 0 rest rest hist (:history @model)]
                    (let [st (first rest)
                          hist (conj hist (:current-value st))
                          finished? (or (nil? (next rest)) (solved? st))]
                      (if (or finished? (>= k (dec (.-speed controls))))
                        (swap! model assoc :rest (if finished? rest (next rest)) :st st :done? finished?
                               ;; Thinned as it grows, so the chart stays cheap.
                               :history (if (> (count hist) 4000) (vec (take-nth 2 hist)) hist))
                        (recur (inc k) (next rest) hist)))))))
            (run-trials! []
              (when-not (tsp?)
                (reset! trials {:done 0 :solved 0 :moves-solved 0 :moves-stuck 0})
                (letfn [(batch []
                          (when-let [{:keys [done]} @trials]
                            (when (< done 200)
                              (doseq [k (range done (min 200 (+ done 5)))]
                                (let [p (problem (+ 1000 k))
                                      end (last (states p))
                                      ok (solved? end)
                                      moves (or (:steps end) 0)]
                                  (swap! trials (fn [t] (cond-> (update t :done inc)
                                                          ok (update :solved inc)
                                                          true (update (if ok :moves-solved :moves-stuck) + moves))))))
                              (js/setTimeout batch 0))))]
                  (batch))))
            (draw-queens! [x0 y0 side rows]
              (let [n (count rows) cell (/ side n)
                    center (fn [c r] [(+ x0 (* cell (+ c 0.5))) (+ y0 (* cell (- n r 0.5)))])]
                (dotimes [c n]
                  (dotimes [r n]
                    (set! (.-fillStyle ctx) (if (even? (+ c r)) "#2a3140" "#202633"))
                    (.fillRect ctx (+ x0 (* cell c)) (+ y0 (* cell (- n 1 r))) cell cell)))
                (set! (.-strokeStyle ctx) "rgba(255,90,90,0.55)")
                (set! (.-lineWidth ctx) (max 1 (* 0.05 cell)))
                (doseq [i (range n) j (range (inc i) n)
                        :let [a (rows i) b (rows j)]
                        :when (or (= a b) (= (abs (- a b)) (- j i)))
                        :let [[ax ay] (center i a) [bx by] (center j b)]]
                  (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
                (set! (.-fillStyle ctx) "#ffd479")
                (dotimes [c n]
                  (let [[x y] (center c (rows c))]
                    (.beginPath ctx) (.arc ctx x y (* 0.32 cell) 0 (* 2 math/PI)) (.fill ctx)))))
            (draw-tour! [x0 y0 side cities tour color width]
              (set! (.-strokeStyle ctx) color)
              (set! (.-lineWidth ctx) width)
              (.beginPath ctx)
              (doseq [[k c] (map-indexed vector (conj tour (first tour)))
                      :let [[x y] (cities c)]]
                ((if (zero? k) #(.moveTo ctx %1 %2) #(.lineTo ctx %1 %2)) (+ x0 (* side x)) (+ y0 (* side (- 1 y)))))
              (.stroke ctx)
              (set! (.-fillStyle ctx) "#ffd479")
              (doseq [[x y] cities]
                (.beginPath ctx) (.arc ctx (+ x0 (* side x)) (+ y0 (* side (- 1 y))) 3 0 (* 2 math/PI)) (.fill ctx)))
            (draw-chart! [x0 y0 w h hist]
              (set! (.-fillStyle ctx) "#151a23")
              (.fillRect ctx x0 y0 w h)
              (when (> (count hist) 1)
                (let [lo (reduce min hist) hi (reduce max hist)
                      span (if (= lo hi) 1 (- hi lo))
                      n (count hist)]
                  (set! (.-strokeStyle ctx) "#7dd3ff")
                  (set! (.-lineWidth ctx) 1.5)
                  (.beginPath ctx)
                  (doseq [[k v] (map-indexed vector hist)
                          :let [x (+ x0 (* w (/ k (dec n)))) y (+ y0 h (- (* (- h 6) (/ (- v lo) span))) -3)]]
                    ((if (zero? k) #(.moveTo ctx %1 %2) #(.lineTo ctx %1 %2)) x y))
                  (.stroke ctx))))
            (draw! []
              (let [{:keys [problem st history done?]} @model
                    cw (.-width canvas) ch (.-height canvas)
                    chart-h 90
                    side (min (- cw 40) (- ch chart-h 90))
                    x0 (* 0.5 (- cw side)) y0 20]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (if (tsp?)
                  (do (doseq [[s] (:beam st)]
                        (draw-tour! x0 y0 side (:cities problem) s "rgba(125,211,255,0.25)" 1))
                      ;; Annealed, the best tour seen, not where it stopped.
                      (draw-tour! x0 y0 side (:cities problem) (if (and done? (:best st)) (:best st) (:current st))
                                  "#7dd3ff" 2))
                  (draw-queens! x0 y0 side (:current st)))
                (draw-chart! x0 (+ y0 side 10) side chart-h history)
                (set! (.-textContent caption)
                      (str (.-algorithm controls) "  |  "
                           (if (:generation st) (str "generation " (:generation st)) (str "move " (:steps st)))
                           (when-let [r (:restart st)] (str "  |  restart " r))
                           (when-let [t (:temperature st)] (str "  |  T " (.toFixed t 4)))
                           (if (tsp?)
                             (str "  |  tour length " (.toFixed (- (:current-value st)) 3)
                                  (when (:best st) (str " (best " (.toFixed (- (:best-value st)) 3) ")")))
                             (str "  |  attacking pairs " (problems/attacking-pairs (:current st))))
                           (cond (solved? st) "  |  solved" (and done? (tsp?)) "  |  done" done? "  |  stuck" :else "")
                           (when-let [{:keys [done solved moves-solved moves-stuck]} @trials]
                             (str "\ntrials: " done "/200 boards, " (.toFixed (* 100 (/ solved (max 1 done))) 1) "% solved"
                                  (when (pos? solved) (str ", " (.toFixed (/ moves-solved solved) 1) " moves when solved"))
                                  (when (< solved done) (str ", " (.toFixed (/ moves-stuck (- done solved)) 1) " when stuck"))))
                           "\nchart: value per move (higher is better)"))))
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
        (.onChange (.add gui controls "problem" #js ["8-queens" "traveling salesman"]) restart!)
        (.onChange (.listen (.add gui controls "algorithm" (clj->js names))) restart!)
        (.onChange (.add gui controls "queens" 4 40 1) restart!)
        (.onChange (.add gui controls "cities" 8 120 1) restart!)
        (.onChange (.add gui controls "seed" 1 99 1) restart!)
        (.listen (.add gui controls "speed" 1 400 1))
        (.add gui #js {:restart restart!} "restart")
        (.add gui #js {:trials run-trials!} "trials"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false) (reset! trials nil))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "local-search")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
