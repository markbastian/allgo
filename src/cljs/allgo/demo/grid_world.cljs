(ns allgo.demo.grid-world
  "Canvas demo for `allgo.decision.mdp` and `allgo.decision.reinforcement`:
  the 4 x 3 world of Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapters 17 and 22.

  Each square shows its utility and an arrow for the action the current
  utilities favor. Value iteration sweeps Bellman updates, a step a sweep,
  and the values spread out from the two exits; policy iteration
  alternates evaluating its policy and improving it, and is done in a
  handful of steps. The learning agents know nothing of the 80/10/10
  slips: the passive ones follow a fixed policy and learn what it is
  worth -- TD from each move, ADP by learning the model -- while
  Q-learning, SARSA and REINFORCE learn what to do, a trial at a time,
  the white dot the agent on its walk."
  (:require [allgo.decision.mdp :as mdp]
            [allgo.decision.reinforcement :as rl]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:method "value iteration" :speed 2})

(def ^:private w mdp/four-by-three)

(def ^:private fixed-policy
  ;; The optimal policy, for the passive learners to evaluate.
  (mdp/policy-iteration w))

(defn- utilities-of
  "What to draw: `{:U :pi :agent :label}` for the method's state."
  [method st]
  (case method
    "value iteration" {:U (:U st) :pi (mdp/best-policy w (:U st)) :label (str "sweep " (:sweeps st))}
    "policy iteration" {:U (:U st) :pi (:pi st) :label (str "iteration " (:iterations st) ", " (:changed st) " actions changed")}
    ("TD agent" "ADP agent") {:U (:U (:memory st)) :pi fixed-policy :agent (:state st) :label (str "trial " (:trial st))}
    ("Q-learning" "SARSA") (let [Q (:Q (:memory st))
                                 U (into {} (for [s (:states w)]
                                              [s (if (seq ((:actions w) s))
                                                   (apply max (map #(get Q [s %] 0.0) ((:actions w) s)))
                                                   (get Q [s nil] 0.0))]))]
                             {:U U :pi (rl/greedy-policy w Q) :agent (:state st) :label (str "trial " (:trial st))})
    "REINFORCE" {:pi (into {} (for [s (:states w)] [s (let [d (rl/softmax-policy w (:theta st) s)] (when (seq d) (key (apply max-key val d))))]))
                 :U {} :label (str "trial " (:trial st))}))

(defn- learner-states
  "A seq of states, one per move, across trials: each `{:state :memory
  :trial}`."
  [[program memory] rng]
  (letfn [(trials [m k]
            (lazy-seq
             (let [t (vec (rl/trial-seq w [program m] [0 0] rng 2000))]
               (concat (map #(assoc % :trial k) t) (trials (:memory (peek t)) (inc k))))))]
    (trials memory 1)))

(defn- states [method]
  (let [rng (random/rng 7)]
    (case method
      "value iteration" (mdp/value-iteration-seq w 0.0001)
      "policy iteration" (mdp/policy-iteration-seq w)
      "TD agent" (learner-states (rl/passive-td-agent w fixed-policy) rng)
      "ADP agent" (learner-states (rl/passive-adp-agent w fixed-policy) rng)
      "Q-learning" (learner-states (rl/q-learning-agent w {:ne 30 :rng rng}) rng)
      "SARSA" (learner-states (rl/sarsa-agent w {:epsilon 0.1 :rng rng}) rng)
      "REINFORCE" (map-indexed (fn [k theta] {:theta theta :trial k}) (rl/reinforce-seq w [0 0] 0.1 rng)))))

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
              (let [ss (states (.-method controls))]
                (reset! model {:rest ss :st (first ss) :clock 0.0 :last-time nil})))
            (advance! []
              (let [now (js/performance.now)
                    {:keys [clock last-time]} @model
                    ;; Planners step a few times a second; learners move faster.
                    rate (* (.-speed controls) (if (#{"value iteration" "policy iteration"} (.-method controls)) 1.0 15.0))
                    clock (+ clock (* rate 0.001 (- now (or last-time now))))]
                (swap! model assoc :last-time now)
                (loop [clock clock]
                  (let [{:keys [rest]} @model]
                    (if (and (>= clock 1.0) (next rest))
                      (do (swap! model assoc :rest (next rest) :st (fnext rest)) (recur (- clock 1.0)))
                      (swap! model assoc :clock (min clock 1.0)))))))
            (draw! []
              (let [{:keys [st]} @model
                    {:keys [U pi agent label]} (utilities-of (.-method controls) st)
                    cw (.-width canvas) ch (.-height canvas)
                    cell (min (/ (- cw 40) 4) (/ (- ch 120) 3))
                    x0 (* 0.5 (- cw (* 4 cell))) y0 50
                    px (fn [[x y]] [(+ x0 (* cell x)) (+ y0 (* cell (- 2 y)))])
                    text (fn [s x y size color]
                           (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str size "px ui-monospace,Menlo,monospace"))
                           (set! (.-textAlign ctx) "center") (set! (.-textBaseline ctx) "middle") (.fillText ctx s x y))]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (doseq [x (range 4) y (range 3)
                        :let [s [x y] [sx sy] (px s) u (get U s)]]
                  (set! (.-fillStyle ctx)
                        (cond (= s [1 1]) "#2b3140"
                              (= s [3 2]) "#1f5a33"
                              (= s [3 1]) "#6a2020"
                              (nil? u) "#161b25"
                              :else (let [t (max -1.0 (min 1.0 u))]
                                      (if (pos? t) (str "rgba(70,190,110," (* 0.6 t) ")") (str "rgba(220,70,70," (* -0.6 t) ")")))))
                  (.fillRect ctx (+ sx 2) (+ sy 2) (- cell 4) (- cell 4))
                  (when (not= s [1 1])
                    (when u (text (.toFixed u 3) (+ sx (* 0.5 cell)) (+ sy (* 0.8 cell)) (int (* 0.12 cell)) "#e6edf7"))
                    (when-let [[dx dy] (pi s)]
                      (let [cx (+ sx (* 0.5 cell)) cy (+ sy (* 0.42 cell)) r (* 0.22 cell)]
                        (set! (.-strokeStyle ctx) "#ffd479") (set! (.-lineWidth ctx) 3)
                        (.beginPath ctx) (.moveTo ctx (- cx (* r dx)) (+ cy (* r dy))) (.lineTo ctx (+ cx (* r dx)) (- cy (* r dy))) (.stroke ctx)
                        (.beginPath ctx) (.arc ctx (+ cx (* r dx)) (- cy (* r dy)) 5 0 (* 2 math/PI)) (set! (.-fillStyle ctx) "#ffd479") (.fill ctx)))))
                (when agent
                  (let [[sx sy] (px agent)]
                    (set! (.-fillStyle ctx) "#ffffff")
                    (.beginPath ctx) (.arc ctx (+ sx (* 0.2 cell)) (+ sy (* 0.2 cell)) (* 0.07 cell) 0 (* 2 math/PI)) (.fill ctx)))
                (set! (.-textContent caption)
                      (str (.-method controls) "  |  " label
                           "\nnumbers: utilities  arrows: the action they favor (moves go that way 80% of the time, 10% to each side)  white: the agent"))))
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
        (.onChange (.add gui controls "method" #js ["value iteration" "policy iteration" "TD agent" "ADP agent"
                                                    "Q-learning" "SARSA" "REINFORCE"])
                   restart!)
        (.add gui controls "speed" 0.25 20 0.25)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "grid-world")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
