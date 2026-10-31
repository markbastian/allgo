(ns allgo.demo.game-search
  "Canvas demo for `allgo.search.adversarial`: Connect Four or
  tic-tac-toe, each side played by a person -- click a column, or a
  square -- or by one of Russell and Norvig's chapter 5 searches.

  Every search is a lazy seq of steps, and the computer thinks `speed`
  steps a frame, so the page never stops to think: the score of each move
  fills in as the search finishes with it, and the caption counts the
  nodes it visited and the moves it never had to look at. Minimax looks
  at everything to `depth`; alpha-beta finds the same value with far
  less; with move ordering and a transposition table, less again. Under
  alpha-beta a score marked ≤ is a bound, not a value: the move was cut
  short once it was clear it could not beat the best. Deepening searches
  1, 2, 3 ... plies, each starting from the last one's choice, until the
  depth or `budget` steps run out. Monte Carlo has no evaluation function
  at all: it plays `playouts` random games through a growing tree, and
  shows each move's playouts and how often they went well.

  It opens with the computer playing itself, Monte Carlo against
  alpha-beta; make either side `human` to play. A depth of 0 means to the end of the game -- quick for tic-tac-toe, and
  for Connect Four, never."
  (:require [allgo.random :as random]
            [allgo.search.adversarial :as adv]
            [allgo.search.games :as games]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:game     "Connect Four"
       :x        "Monte Carlo"
       :o        "alpha-beta + ordering + table"
       :depth    6
       :budget   60000
       :playouts 3000
       :speed    1500})

(def ^:private players
  ["human" "minimax" "alpha-beta" "alpha-beta + ordering + table" "iterative deepening" "Monte Carlo"])

(defn- game [] (if (= "Connect Four" (.-game controls)) games/connect-four games/tic-tac-toe))

(defn- c4? [] (= "Connect Four" (.-game controls)))

(defn- depth []
  (let [d (.-depth controls)] (when (pos? d) d)))

(defonce ^:private searches (atom 0))

(defn- start-search
  "A new search for the player to move in `state`: its kind and its seq."
  [kind state]
  (let [g (game)]
    (case kind
      "minimax" {:kind :tree :rest (adv/minimax-seq g state {:depth (depth)})}
      "alpha-beta" {:kind :tree :rest (adv/alpha-beta-seq g state {:depth (depth)})}
      "alpha-beta + ordering + table" {:kind :tree :rest (adv/alpha-beta-seq g state {:depth (depth) :order? true :table? true})}
      "iterative deepening" {:kind :deepening :rest (adv/deepening-seq g state {:max-depth (or (depth) 42) :table? true})}
      ;; A stream of its own each search, so games differ, and repeat.
      "Monte Carlo" {:kind :mcts :rest (adv/monte-carlo-seq g state {:rng (random/rng (swap! searches inc))})}
      nil)))

(defn- fmt [v bound]
  (str (when (= bound :upper) "≤")
       (cond (>= v 1.0) "win" (<= v -1.0) "loss" :else (.toFixed v 2))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        layout (atom nil)]
    (set! (.-style canvas) "display:block;cursor:pointer")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(new-game! []
              (reset! model {:state (:initial (game)) :search nil :last-search nil :pause 0}))
            (player-of [side] (if (= side :x) (.-x controls) (.-o controls)))
            (play! [move]
              (swap! model (fn [m] (assoc m :state ((:result (game)) (:state m) move) :search nil :pause 20))))
            (finish-search
              ;; The move a search, run as far as it has been, chooses.
              [{:keys [kind st done-pass]}]
              (case kind
                :mcts (:best (adv/monte-carlo-choice st))
                :deepening (or (:best (:context done-pass)) (:best (:context st)))
                (:best (:context st))))
            (think! []
              (let [{:keys [state search pause]} @model
                    g (game)]
                (cond
                  ((:terminal? g) state) nil
                  (pos? pause) (swap! model update :pause dec)
                  (= "human" (player-of (:to-move state))) nil
                  (nil? search) (swap! model assoc :search (assoc (start-search (player-of (:to-move state)) state)
                                                                  :steps 0))
                  :else
                  (let [{:keys [kind rest steps]} search
                        budget (case kind :mcts (.-playouts controls) :deepening (.-budget controls) nil)
                        n (if (= kind :mcts) (max 1 (quot (.-speed controls) 40)) (.-speed controls))
                        [rest st done-pass steps]
                        (loop [k 0 rest rest st (:st search) done-pass (:done-pass search) steps steps]
                          (if (or (= k n) (nil? (seq rest)) (and budget (>= steps budget)))
                            [rest st done-pass steps]
                            (let [s (first rest)]
                              (recur (inc k) (next rest) s
                                     ;; A deepening pass ends on a state with a value.
                                     (if (and (= kind :deepening) (contains? s :value)) s done-pass)
                                     (inc steps)))))
                        search (assoc search :rest rest :st st :done-pass done-pass :steps steps)
                        over? (or (nil? (seq rest)) (and budget (>= steps budget)))]
                    (if over?
                      (do (swap! model assoc :last-search search)
                          (play! (finish-search search)))
                      (swap! model assoc :search search))))))
            (scores
              ;; What to write on each move: from the search thinking now,
              ;; or the one that chose the last move.
              []
              (let [{:keys [search last-search]} @model
                    {:keys [kind st done-pass]} (or search last-search)]
                (cond
                  ;; Nothing to show before a search's first step.
                  (or (nil? kind) (nil? st)) {}
                  (= kind :mcts) (into {} (for [[m {:keys [n q]}] (:moves (adv/monte-carlo-choice st))
                                                :when (pos? n)]
                                            [m (str n " / " (.toFixed q 2))]))
                  :else
                  (let [c (if (and (= kind :deepening) done-pass (nil? (:value st))) (:context done-pass) (:context st))]
                    (into {} (for [[m v] (:scores c)] [m (fmt v (get-in c [:bounds m]))]))))))
            (draw! []
              (let [{:keys [state search last-search]} @model
                    g (game)
                    cw (.-width canvas) ch (.-height canvas)
                    [cols rws] (if (c4?) [games/columns games/rows] [3 3])
                    cell (min (/ (- cw 40) cols) (/ (- ch 140) (inc rws)))
                    x0 (* 0.5 (- cw (* cell cols))) y0 (+ 20 cell)
                    board (:board state)
                    marks (scores)]
                (reset! layout {:x0 x0 :y0 y0 :cell cell :cols cols :rows rws})
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (set! (.-fillStyle ctx) (if (c4?) "#1d3a8a" "#1a202b"))
                (.fillRect ctx x0 y0 (* cell cols) (* cell rws))
                (dotimes [c cols]
                  (dotimes [r rws]
                    (let [k (if (c4?) (+ c (* games/columns r)) (+ c (* 3 r)))
                          ;; Connect Four's row 0 is the bottom.
                          ry (if (c4?) (- rws 1 r) r)
                          cx (+ x0 (* cell (+ c 0.5))) cy (+ y0 (* cell (+ ry 0.5)))
                          p (board k)]
                      (if (c4?)
                        (do (set! (.-fillStyle ctx) (case p :x "#ffcc33" :o "#e8453c" "#0d1118"))
                            (.beginPath ctx) (.arc ctx cx cy (* 0.4 cell) 0 (* 2 math/PI)) (.fill ctx))
                        (do (set! (.-strokeStyle ctx) "#3a4356")
                            (.strokeRect ctx (+ x0 (* cell c)) (+ y0 (* cell r)) cell cell)
                            (when p
                              (set! (.-fillStyle ctx) (if (= p :x) "#ffcc33" "#e8453c"))
                              (set! (.-font ctx) (str (int (* 0.6 cell)) "px sans-serif"))
                              (set! (.-textAlign ctx) "center") (set! (.-textBaseline ctx) "middle")
                              (.fillText ctx (if (= p :x) "X" "O") cx cy)))))))
                ;; The scores: above each column, or in each empty square.
                (set! (.-font ctx) "12px ui-monospace,Menlo,monospace")
                (set! (.-textAlign ctx) "center") (set! (.-textBaseline ctx) "middle")
                (doseq [[m text] marks]
                  (let [[cx cy] (if (c4?)
                                  [(+ x0 (* cell (+ m 0.5))) (- y0 (* 0.5 cell))]
                                  [(+ x0 (* cell (+ (mod m 3) 0.5))) (+ y0 (* cell (+ (quot m 3) 0.5)))])
                        {:keys [kind st]} (or search last-search)
                        best? (= m (:best (if (= :mcts kind) (adv/monte-carlo-choice st) (:context st))))]
                    (set! (.-fillStyle ctx) (if best? "#7dffa0" "#c9d2e0"))
                    (.fillText ctx text cx cy)))
                (set! (.-textContent caption)
                      (let [to-move (:to-move state)
                            {:keys [kind st steps]} (or search last-search)
                            c (:context st)]
                        (str (cond
                               (:winner state) (str (name (:winner state)) " wins")
                               ((:terminal? g) state) "a draw"
                               (= "human" (player-of to-move)) (str (name to-move) " to play: click "
                                                                    (if (c4?) "a column" "a square"))
                               search (str (name to-move) " (" (player-of to-move) ") is thinking")
                               :else (str (name to-move) " to play"))
                             (when kind
                               (str "\n" (if search "this search" "last search") ": "
                                    (case kind
                                      :mcts (str (:playouts st) " playouts")
                                      (str steps " steps, " (:nodes c) " nodes"
                                           (when (pos? (:pruned c)) (str ", " (:pruned c) " moves pruned"))
                                           (when-let [d (:depth st)] (str ", depth " d))))))
                             "\nscores: the value of each move to the player choosing"
                             (when (= kind :mcts) " -- Monte Carlo: playouts / average result"))))))
            (click! [^js e]
              (let [{:keys [x0 y0 cell cols rows]} @layout
                    r (.getBoundingClientRect canvas)
                    x (- (.-clientX e) (.-left r)) y (- (.-clientY e) (.-top r))
                    c (long (math/floor (/ (- x x0) cell)))
                    rw (long (math/floor (/ (- y y0) cell)))
                    {:keys [state]} @model
                    g (game)
                    move (if (c4?) c (+ c (* 3 rw)))]
                (when (and (= "human" (player-of (:to-move state)))
                           (not ((:terminal? g) state))
                           (< -1 c cols)
                           (or (c4?) (< -1 rw rows))
                           (some #{move} ((:actions g) state)))
                  (play! move))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (think!)
                (draw!)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h))))]
      (.addEventListener canvas "click" click!)
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "game" #js ["Connect Four" "tic-tac-toe"])
                   (fn [v] (set! (.-depth controls) (if (= v "tic-tac-toe") 0 6)) (new-game!)))
        (.onChange (.add gui controls "x" (clj->js players)) new-game!)
        (.onChange (.add gui controls "o" (clj->js players)) new-game!)
        (.listen (.add gui controls "depth" 0 10 1))
        (.add gui controls "budget" 1000 400000 1000)
        (.add gui controls "playouts" 100 30000 100)
        (.add gui controls "speed" 100 20000 100)
        (.add gui #js {:new-game new-game!} "new-game"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (new-game!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "game-search")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
