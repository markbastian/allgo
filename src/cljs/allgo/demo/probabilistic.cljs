(ns allgo.demo.probabilistic
  "Canvas demo for `allgo.learning.probabilistic`: the learning of
  probabilistic models of Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 20, a step at a time.

  Candy: candies unwrapped from a bag of one of the five kinds, the
  posterior over the kinds after each, and the predictions -- Bayesian,
  averaged over every kind, and MAP, from the likeliest alone. Beta: the
  beta posterior over a bag's fraction of cherries narrowing as candies
  come, from the uniform prior. EM: the chapter's two bags mixed
  together, their parameters found by EM an iteration at a time, the log
  likelihood passing the true model's. Density: points from three blobs,
  and the density estimated from them by Gaussian kernels or by nearest
  neighbors, beside the truth. Structure: samples of the sprinkler
  network, and hill climbing finding its links, a move a step."
  (:require [allgo.learning.probabilistic :as pl]
            [allgo.probability.bayes :as bayes]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]
            [clojure.string :as str]))

(def ^:private ^js controls
  #js {:method "candy" :bag 5 :estimator "kernels" :width 0.07 :k 10 :seed 1 :speed 3})

(def ^:private kinds [1.0 0.75 0.5 0.25 0.0])
(def ^:private bags (mapv (fn [cherry] (fn [d] (if (= d :cherry) cherry (- 1.0 cherry)))) kinds))
(def ^:private bag-prior [0.1 0.2 0.4 0.2 0.1])
(def ^:private colors ["#ff8f70" "#ffd479" "#7dffa0" "#7dc4ff" "#d28fff"])

(def ^:private candy-counts
  "Section 20.3.2's thousand candies from two bags mixed: [[cherry? red?
  hole?] count]."
  [[[true true true] 273] [[true true false] 93] [[true false true] 104] [[true false false] 90]
   [[false true true] 79] [[false true false] 100] [[false false true] 94] [[false false false] 167]])

(def ^:private blobs [[[0.25 0.3] 0.1] [[0.7 0.7] 0.12] [[0.75 0.25] 0.09]])

(defn- truth-density [[x y]]
  (/ (reduce + (for [[[mx my] s] blobs]
                 (/ (math/exp (- (/ (+ (* (- x mx) (- x mx)) (* (- y my) (- y my))) (* 2 s s)))) (* 2 math/PI s s))))
     (count blobs)))

(def ^:private node-at
  {:Cloudy [0.5 0.12] :Sprinkler [0.2 0.5] :Rain [0.8 0.5] :WetGrass [0.5 0.88]})

(defn- init [method]
  (let [rng (random/rng (.-seed controls))]
    (case method
      "candy" (let [bag (nth kinds (dec (.-bag controls)))
                    data (vec (repeatedly 30 #(if (< (rng) bag) :cherry :lime)))]
                {:data data :posts (vec (pl/posterior-seq bags bag-prior data)) :shown 0})
      "beta" (let [data (vec (repeatedly 200 #(< (rng) 0.75)))]
               {:data data :betas (vec (pl/beta-seq [1 1] data)) :shown 0})
      "EM" {:states (vec (take 31 (pl/naive-bayes-em-seq {:prior [0.6 0.4] :features [[0.6 0.6 0.6] [0.4 0.4 0.4]]} candy-counts)))
            :truth (pl/mixture-log-likelihood {:prior [0.5 0.5] :features [[0.8 0.8 0.8] [0.3 0.3 0.3]]} candy-counts)
            :shown 0}
      "density" {:points (vec (repeatedly 128 #(let [[[mx my] s] (random/pick rng blobs)]
                                                 [(random/gaussian rng mx s) (random/gaussian rng my s)])))}
      "structure" (let [data (vec (repeatedly 3000 #(bayes/prior-sample bayes/sprinkler rng)))]
                    ;; walked a move a tick, the search being slow to run all at once
                    {:states (pl/structure-search-seq (zipmap (:variables bayes/sprinkler) (repeat [true false])) data
                                                      :restarts 10 :rng (random/rng (+ 100 (.-seed controls))))
                     :step 0 :shown 0}))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        scratch (js/document.createElement "canvas")
        running? (atom false)
        model (atom nil)
        clock (atom 0.0)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! [] (reset! model (init (.-method controls))) (reset! clock 0.0))
            (advance! []
              ;; a step every so often, faster with the speed
              (swap! clock + (* 0.016 (.-speed controls)))
              (when (>= @clock 1.0)
                (reset! clock 0.0)
                (swap! model (fn [m] (if (contains? m :shown) (update m :shown inc) m)))))
            (text! [s x y color & [size align]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (set! (.-textAlign ctx) (or align "left")) (.fillText ctx s x y) (set! (.-textAlign ctx) "left"))
            (axes! [x y w h]
              (set! (.-strokeStyle ctx) "#2a3242") (set! (.-lineWidth ctx) 1) (.strokeRect ctx x y w h))
            (plot! [x y w h pts lo hi color]
              ;; a polyline of [t v] with t in [0 1], v in [lo hi]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 2)
              (.beginPath ctx)
              (doseq [[i [t v]] (map-indexed vector pts)]
                (let [px (+ x (* w t)) py (- (+ y h) (* h (/ (- v lo) (- hi lo))))]
                  (if (zero? i) (.moveTo ctx px py) (.lineTo ctx px py))))
              (.stroke ctx))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas) method (.-method controls)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "candy" (draw-candy! cw ch)
                             "beta" (draw-beta! cw ch)
                             "EM" (draw-em! cw ch)
                             "density" (draw-density! cw ch)
                             "structure" (draw-structure! cw ch))))))
            (draw-candy! [cw ch]
              (let [{:keys [data posts shown]} @model
                    n (min shown 30)
                    w (- (* 0.5 cw) 70) h (- ch 260) y 110]
                (axes! 40 y w h) (axes! (+ 60 w) y w h)
                (text! "P(h | d): the posterior over the five kinds of bag" 40 (- y 8) "#9fb4d8" 11)
                (text! "P(next is lime | d): Bayesian (green), MAP (red)" (+ 60 w) (- y 8) "#9fb4d8" 11)
                (doseq [i (range 5)]
                  (plot! 40 y w h (map (fn [k] [(/ k 30.0) ((posts k) i)]) (range (inc n))) 0.0 1.0 (colors i))
                  (text! (str "h" (inc i)) (+ 46 (* w (/ n 30.0))) (- (+ y h) (* h ((posts n) i)) 4) (colors i) 11))
                (plot! (+ 60 w) y w h (map (fn [k] [(/ k 30.0) (pl/predict bags (posts k) :lime)]) (range (inc n))) 0.0 1.0 "#7dffa0")
                (plot! (+ 60 w) y w h (map (fn [k] [(/ k 30.0) ((bags (pl/map-hypothesis (posts k))) :lime)]) (range (inc n))) 0.0 1.0 "#ff8f70")
                (text! (str/join " " (map #(if (= % :cherry) "C" "L") (take n data))) 40 (+ y h 28) "#e6edf7" 13)
                (when (>= shown 40) (swap! model assoc :shown 0))
                (str n " candies from an h" (.-bag controls) " bag  |  most probable: h" (inc (pl/map-hypothesis (posts n)))
                     "  |  P(next lime) " (.toFixed (pl/predict bags (posts n) :lime) 3)
                     "\nthe prior 0.1 0.2 0.4 0.2 0.1, as the manufacturer advertises; C cherry, L lime")))
            (draw-beta! [cw ch]
              (let [{:keys [data betas shown]} @model
                    ;; one candy a step, then faster
                    n (min 200 (if (< shown 20) shown (+ 20 (* 6 (- shown 20)))))
                    [a b] (betas n)
                    w (- cw 80) h (- ch 260) y 110
                    peak (max 1.0 (apply max (map #(pl/beta-density a b (/ (+ % 0.5) 200.0)) (range 200))))]
                (axes! 40 y w h)
                (plot! 40 y w h (map (fn [i] (let [t (/ (+ i 0.5) 200.0)] [t (pl/beta-density a b t)])) (range 200)) 0.0 peak "#ffd479")
                (set! (.-strokeStyle ctx) "rgba(125,255,160,0.6)")
                (.beginPath ctx) (.moveTo ctx (+ 40 (* w 0.75)) y) (.lineTo ctx (+ 40 (* w 0.75)) (+ y h)) (.stroke ctx)
                (text! "true fraction 0.75" (+ 44 (* w 0.75)) (+ y 14) "#7dffa0" 11)
                (text! "0" 40 (+ y h 16) "#9fb4d8" 11) (text! "1" (+ 30 w) (+ y h 16) "#9fb4d8" 11)
                (when (>= shown 60) (swap! model assoc :shown 0))
                (str n " candies, " (count (filter true? (take n data))) " cherry  |  posterior beta[" a ", " b "], mean "
                     (.toFixed (pl/beta-mean [a b]) 3)
                     "\nthe beta family is closed under update: each cherry adds one to a, each lime one to b -- virtual counts")))
            (draw-em! [cw ch]
              (let [{:keys [states truth shown]} @model
                    i (min shown 30)
                    {:keys [params log-likelihood]} (states i)
                    w (- (* 0.5 cw) 70) h (- ch 260) y 110
                    lls (mapv :log-likelihood states)
                    lo (apply min lls) hi (+ 5 (max truth (apply max lls)))]
                (axes! 40 y w h)
                (text! "log likelihood by iteration" 40 (- y 8) "#9fb4d8" 11)
                (plot! 40 y w h (map (fn [k] [(/ k 30.0) (lls k)]) (range (inc i))) lo hi "#ffd479")
                (plot! 40 y w h [[0.0 truth] [1.0 truth]] lo hi "rgba(125,255,160,0.5)")
                (text! "the model that made the data" 46 (- (+ y h) (* h (/ (- truth lo) (- hi lo))) 6) "#7dffa0" 11)
                (let [x0 (+ 80 w)]
                  (text! "parameter          learned   true" x0 (- y 8) "#9fb4d8" 11)
                  (doseq [[r [label v t]] (map-indexed vector
                                                       (concat [["P(bag 1)" ((:prior params) 0) 0.5]]
                                                               (for [b [0 1] [f nm] (map vector (range 3) ["cherry" "red" "hole"])]
                                                                 [(str "bag " (inc b) " " nm) (get-in params [:features b f]) (if (zero? b) 0.8 0.3)])))]
                    (let [yy (+ y 18 (* 24 r))]
                      (text! label x0 yy "#e6edf7" 12)
                      (set! (.-fillStyle ctx) "#7dffa0") (.fillRect ctx (+ x0 130) (- yy 10) (* 160 v) 9)
                      (set! (.-strokeStyle ctx) "#ffd479") (.strokeRect ctx (+ x0 130) (- yy 10) (* 160 t) 9)
                      (text! (.toFixed v 4) (+ x0 300) yy "#7dffa0" 12))))
                (when (>= shown 45) (swap! model assoc :shown 0))
                (str "iteration " i "  |  log likelihood " (.toFixed log-likelihood 2) " against the true model's " (.toFixed truth 3)
                     "\nthe chapter's two bags mixed, the bag hidden: EM from 0.6 and 0.4; bars learned, outlines true")))
            (draw-density! [cw ch]
              (let [{:keys [points cache]} @model
                    kernels? (= "kernels" (.-estimator controls))
                    key [kernels? (.-width controls) (.-k controls)]
                    n 72
                    grid (if (= key (first cache))
                           (second cache)
                           (let [g (vec (for [j (range n) i (range n) :let [p [(/ (+ i 0.5) n) (- 1.0 (/ (+ j 0.5) n))]]]
                                          [(truth-density p)
                                           (if kernels? (pl/kernel-density points (.-width controls) p) (pl/knn-density points (.-k controls) p))]))]
                             (swap! model assoc :cache [key g])
                             g))
                    top (apply max (map first grid))
                    side (min (- (* 0.5 cw) 60) (- ch 220))
                    paint! (fn [x0 which]
                             (set! (.-width scratch) n) (set! (.-height scratch) n)
                             (let [c2 (.getContext scratch "2d") data (.createImageData c2 n n) d (.-data data)]
                               (dotimes [i (* n n)]
                                 (let [v (min 1.0 (/ (which (grid i)) top)) o (* 4 i)]
                                   (aset d o (* 255 (min 1 (* 1.4 v)))) (aset d (+ o 1) (* 255 (* v v))) (aset d (+ o 2) (* 255 (* 0.5 (- 1 v))))
                                   (aset d (+ o 3) 255)))
                               (.putImageData c2 data 0 0))
                             (set! (.-imageSmoothingEnabled ctx) false)
                             (.drawImage ctx scratch x0 110 side side))]
                (paint! 40 second)
                (set! (.-fillStyle ctx) "#e6edf7")
                (doseq [[x y] points] (.fillRect ctx (+ 40 (* side x) -1) (+ 110 (* side (- 1 y)) -1) 3 3))
                (paint! (+ 80 side) first)
                (text! (if kernels? (str "Gaussian kernels, width " (.-width controls)) (str (.-k controls) " nearest neighbors")) 40 100 "#9fb4d8" 11)
                (text! "the true density" (+ 80 side) 100 "#9fb4d8" 11)
                (str "128 points from three blobs"
                     "\nkernels: too narrow is spiky, too wide too smooth, 0.07 about right; neighbors: k = 10 between 3 and 40")))
            (draw-structure! [cw ch]
              (let [{:keys [states step shown]} @model
                    _ (when (and (> shown step) (next states)) (swap! model #(-> % (update :states rest) (update :step inc))))
                    {:keys [states step]} @model
                    done? (nil? (next states))
                    {:keys [parents score move restart best]} (first states)
                    side (min (- cw 80) (- ch 240))
                    at (fn [v] (let [[x y] (node-at v)] [(+ 40 (* side x)) (+ 110 (* side y 0.8))]))
                    arrow! (fn [a b color width]
                             (let [[ax ay] (at a) [bx by] (at b)
                                   dx (- bx ax) dy (- by ay) l (math/hypot dx dy)
                                   [ux uy] [(/ dx l) (/ dy l)]
                                   [ex ey] [(- bx (* 46 ux)) (- by (* 26 uy))]]
                               (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                               (.beginPath ctx) (.moveTo ctx (+ ax (* 46 ux)) (+ ay (* 26 uy))) (.lineTo ctx ex ey) (.stroke ctx)
                               (.beginPath ctx) (.moveTo ctx ex ey)
                               (.lineTo ctx (- ex (* 10 ux) (* -6 uy)) (- ey (* 10 uy) (* 6 ux)))
                               (.lineTo ctx (- ex (* 10 ux) (* 6 uy)) (- ey (* 10 uy) (* -6 ux)))
                               (.closePath ctx) (set! (.-fillStyle ctx) color) (.fill ctx)))]
                (doseq [[v ps] (:parents bayes/sprinkler) p ps] (arrow! p v "rgba(255,255,255,0.12)" 8))
                (doseq [[v ps] (if done? (:parents best) parents) p ps] (arrow! p v (if done? "#7dffa0" "#ffd479") 2.5))
                (doseq [v (keys node-at)]
                  (let [[x y] (at v)]
                    (set! (.-fillStyle ctx) "#1b2230") (.fillRect ctx (- x 46) (- y 16) 92 32)
                    (set! (.-strokeStyle ctx) "#9fb4d8") (set! (.-lineWidth ctx) 1) (.strokeRect ctx (- x 46) (- y 16) 92 32)
                    (text! (name v) x (+ y 4) "#e6edf7" 12 "center")))
                (when (and done? (> shown (+ step 12))) (restart!))
                (str (if (zero? restart) "the climb from no links" (str "restart " restart " of 10, from a random network"))
                     ", step " step (when move (str ": " (name (first move)) " " (name (second move)) " -> " (name (nth move 2))))
                     "  |  score " (.toFixed score 1) "  |  best so far " (.toFixed (:score best) 1)
                     (when done? "  |  done: the best network is shown")
                     "\nhill climbing on 3000 samples, scored by log likelihood less BIC's penalty, restarting from random networks; links found: yellow while climbing, green the best; faint: the true network's")))
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
        (.onChange (.add gui controls "method" #js ["candy" "beta" "EM" "density" "structure"]) restart!)
        (.onChange (.add gui controls "bag" 1 5 1) restart!)
        (.add gui controls "estimator" #js ["kernels" "neighbors"])
        (.add gui controls "width" #js [0.02 0.07 0.2])
        (.add gui controls "k" #js [3 10 40])
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 20 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "probabilistic")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
