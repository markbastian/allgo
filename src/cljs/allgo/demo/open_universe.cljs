(ns allgo.demo.open-universe
  "Canvas demo for chapter 15 of Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition: probabilistic models of
  worlds with unknown objects, a step at a time.

  Ratings: skill rated from team games alone (15.1.2), each player's
  rating narrowing round their true skill as the games come. Tracking:
  two aircraft crossing, blips without labels (15.3.1) -- the
  nearest-neighbor filter, the most probable assignment and a particle
  filter, against the true paths; turn on clutter to see the first two
  lose their way. Citations: which citations are of the same paper
  (15.2.3), the section's own two among them, by Gibbs sampling over
  groupings. Seismic: events found among detections nine in ten of which
  are false alarms, as NET-VISA does (15.2.3). Reading: a noisy image of
  a word read by MCMC over the generative program's traces (15.4), and
  the exact readings with independent letters and with a bigram model
  of English."
  (:require [allgo.numerics.linear :as lin]
            [allgo.probability.association :as association]
            [allgo.probability.identity :as identity]
            [allgo.probability.rating :as rating]
            [allgo.probability.reading :as reading]
            [allgo.probability.seismic :as seismic]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:method "tracking" :clutter true :word "uncertainty" :noise 0.3 :seed 1 :speed 4})

(def ^:private colors ["#ff8f70" "#ffd479" "#7dffa0" "#7dc4ff" "#d28fff" "#ff7ab8" "#9fe3d6" "#c9d2e0"])

;; ---------------------------------------------------------------- ratings

(def ^:private skills {:ana 12.0 :ben 16.0 :cy 20.0 :dee 24.0 :eli 28.0 :fay 32.0 :gus 36.0 :hal 40.0})
(def ^:private beta 4.0)

(defn- ratings-init [rng]
  (let [players (vec (keys skills))
        perform (fn [team] (reduce + (map #(random/gaussian rng (skills %) beta) team)))
        games (vec (for [_ (range 400)
                         :let [[a b c d] (take 4 (random/shuffle rng players))]]
                     (if (> (perform [a b]) (perform [c d])) {:winner [a b] :loser [c d]} {:winner [c d] :loser [a b]})))]
    {:games games
     :states (vec (rating/rating-seq (zipmap players (repeat (rating/rating 25 8))) games {:beta beta :tau 0.05}))
     :shown 0}))

;; --------------------------------------------------------------- tracking

(def ^:private track-model (association/constant-velocity 0.01 0.5))
(def ^:private starts [[0.0 0.0 1.0 0.3] [0.0 6.0 1.0 -0.3]])
(def ^:private area [[-2.0 -3.0] [22.0 9.0]])

(defn- tracking-init [rng clutter?]
  (let [m (if clutter? (assoc track-model :p-detect 0.9 :clutter (/ 2.0 (* 24.0 12.0))) track-model)
        step (fn [x] (mapv + (lin/mat-vec (:F m) x) [0.0 0.0 (random/gaussian rng 0.0 0.05) (random/gaussian rng 0.0 0.05)]))
        paths (mapv (fn [s] (vec (take 20 (rest (iterate step s))))) starts)
        [[x0 y0] [x1 y1]] area
        scans (vec (for [t (range 20)]
                     (random/shuffle rng (concat (for [p paths :when (< (rng) (:p-detect m))]
                                                   (let [[x y] (p t)] [(random/gaussian rng x 0.5) (random/gaussian rng y 0.5)]))
                                                 (for [_ (range 10) :when (and clutter? (< (rng) 0.2))]
                                                   [(+ x0 (* (rng) (- x1 x0))) (+ y0 (* (rng) (- y1 y0)))])))))
        prior (mapv (fn [s] [s (lin/mat-scale (lin/eye 4) 0.25)]) starts)]
    {:model m :paths paths :scans scans
     :nearest (vec (association/nearest-neighbor-seq m prior scans))
     :probable (vec (association/most-probable-seq m prior scans))
     :particles (vec (association/particle-seq m prior scans rng 150))
     :shown 0}))

;; -------------------------------------------------------------- citations

(def ^:private citation-texts
  [["Lashkari, Metral, Maes" "Yezdi Lashkari, Max Metral, and Pattie Maes" "Collaborative Interface Agents"]
   ["Metral, Lashkari, Maes" "Metral M. Lashkari, Y. and P. Maes" "Collaborative interface agents"]
   ["Lashkari et al." "Lashkari, Y., Metral, M., Maes, P." "Collaborative agents for interfaces"]
   ["Maes 1994" "Maes, P." "Agents that reduce work and information overload"]
   ["Maes, CACM" "P. Maes" "Agents that reduce work and information overload"]
   ["Russell & Norvig" "Russell, S. and Norvig, P." "Artificial Intelligence: A Modern Approach"]
   ["Russell, Norvig" "S. Russell, P. Norvig" "Artificial Intelligence a modern approach"]])

(defn- citations-init [rng]
  (let [words (mapv (fn [[_ authors title]] (identity/words (str authors " " title))) citation-texts)
        m (identity/model words {})]
    {:samples (vec (take 400 (identity/gibbs-seq m words rng))) :shown 0}))

;; ---------------------------------------------------------------- seismic

(def ^:private stations (mapv #(* 2.0 math/PI (/ % 12.0)) (range 12)))
(def ^:private seismic-model (seismic/model stations :false-rate 9.0 :event-rate 1.0))

(defn- seismic-init [rng]
  (let [events [{:place (* 2.0 math/PI (rng)) :time (math/round (+ 15 (* 30 (rng)))) :mag 5.0}
                {:place (* 2.0 math/PI (rng)) :time (math/round (+ 55 (* 30 (rng)))) :mag 5.0}]
        detections (into {} (for [[s station] (map-indexed vector stations)]
                              [s (vec (concat (for [_ (range 100) :when (< (rng) 0.09)] (* 100.0 (rng)))
                                              (for [{:keys [place time mag]} events
                                                    :when (< (rng) (seismic/p-detect seismic-model station place mag))]
                                                (random/gaussian rng (seismic/arrival seismic-model station place time) 0.25))))]))]
    {:events events :detections detections
     :bulletin (seismic/bulletin seismic-model detections (vec (seismic/grid seismic-model 90 100)))
     :shown 0}))

;; ---------------------------------------------------------------- reading

(def ^:private bigrams (delay (reading/english-bigrams reading/common-words)))

(defn- reading-init [rng word noise]
  (let [image (reading/degrade rng (reading/render word 14 false) noise)
        m (reading/model)]
    {:word word :image image
     :trace (vec (take 60 (reading/mcmc-seq m image rng)))
     :independent (:text (reading/exact-posterior m image))
     :english (:text (reading/exact-posterior (reading/model :bigrams @bigrams) image))
     :shown 0}))

(defn- init [method]
  (let [rng (random/rng (.-seed controls))]
    (assoc (case method
             "ratings" (ratings-init rng)
             "tracking" (tracking-init rng (.-clutter controls))
             "citations" (citations-init rng)
             "seismic" (seismic-init rng)
             "reading" (reading-init rng (.-word controls) (.-noise controls)))
           :method method)))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
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
              (swap! clock + (* 0.016 (.-speed controls)))
              (when (>= @clock 1.0)
                (reset! clock 0.0)
                (swap! model update :shown inc)))
            (text! [s x y color & [size align]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (set! (.-textAlign ctx) (or align "left")) (.fillText ctx s x y) (set! (.-textAlign ctx) "left"))
            (line! [ax ay bx by color width]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
              (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
            (dot! [x y r color]
              (set! (.-fillStyle ctx) color) (.beginPath ctx) (.arc ctx x y r 0 (* 2 math/PI)) (.fill ctx))
            (polyline! [pts color width]
              (when (seq pts)
                (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                (.beginPath ctx)
                (doseq [[i [x y]] (map-indexed vector pts)] (if (zero? i) (.moveTo ctx x y) (.lineTo ctx x y)))
                (.stroke ctx)))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas) method (:method @model)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "ratings" (draw-ratings! cw ch)
                             "tracking" (draw-tracking! cw ch)
                             "citations" (draw-citations! cw ch)
                             "seismic" (draw-seismic! cw ch)
                             "reading" (draw-reading! cw ch))))))
            (draw-ratings! [cw ch]
              (let [{:keys [states games shown]} @model
                    k (min shown (dec (count states)))
                    rs (states k)
                    players (sort-by skills (keys skills))
                    x0 140 x1 (- cw 60) y0 110 row (/ (- ch 230) (count players))
                    px (fn [s] (+ x0 (* (- x1 x0) (/ (- s 0) 52.0))))]
                (doseq [s (range 0 53 4)]
                  (line! (px s) (- y0 20) (px s) (+ y0 (* row (count players))) "#1b2230" 1)
                  (text! (str s) (px s) (- y0 26) "#56627a" 10 "center"))
                (doseq [[i p] (map-indexed vector players)
                        :let [{:keys [mu sigma]} (rs p) y (+ y0 (* row i) (/ row 2))]]
                  (text! (name p) 40 (+ y 4) (colors i) 12)
                  (line! (px (- mu (* 2 sigma))) y (px (+ mu (* 2 sigma))) y (colors i) 4)
                  (dot! (px mu) y 5 "#ffffff")
                  (line! (px (skills p)) (- y 9) (px (skills p)) (+ y 9) "#e6edf7" 2))
                (when (pos? k)
                  (let [{:keys [winner loser]} (games (dec k))]
                    (text! (str "game " k ": " (apply str (interpose " + " (map name winner))) " beat " (apply str (interpose " + " (map name loser))))
                           40 80 "#9fb4d8" 12)))
                (when (>= shown (+ (count states) 20)) (restart!))
                (str "2-on-2 games, the teams reshuffled every game; each bar is a rating, mean plus or minus two sigma; the tick is the true skill"
                     "\nonly team results are seen, yet the ratings sort the players and narrow round their skills (they agree up to a shift: only differences show)")))
            (draw-tracking! [cw ch]
              (let [{:keys [model paths scans nearest probable particles shown]} @model
                    k (min shown 19)
                    [[ax ay] [bx by]] area
                    s (min (/ (- cw 80) (- bx ax)) (/ (- ch 200) (- by ay)))
                    at (fn [[x y]] [(+ 40 (* s (- x ax))) (+ 90 (* s (- by y)))])
                    pos (fn [track] (at (association/position model track)))]
                (set! (.-fillStyle ctx) "#141a24") (let [[x y] (at [ax by])] (.fillRect ctx x y (* s (- bx ax)) (* s (- by ay))))
                (doseq [p paths] (polyline! (map #(at (subvec % 0 2)) (take (inc k) p)) "rgba(255,255,255,0.18)" 6))
                (doseq [t (range (inc k)) z (scans t)]
                  (let [[x y] (at z)] (dot! x y 2.5 (if (= t k) "#e6edf7" "rgba(200,210,230,0.35)"))))
                (doseq [p (particles k)]
                  (doseq [tr (:tracks p)] (let [[x y] (pos tr)] (dot! x y 1.2 "rgba(125,255,160,0.25)"))))
                (doseq [[seq-of color] [[nearest "#ff8f70"] [probable "#ffd479"]]
                        i (range 2)]
                  (polyline! (map #(pos ((seq-of %) i)) (range (inc k))) color 1.5))
                (doseq [i (range 2)]
                  (polyline! (map (fn [t] (at (lin/mean (map #(association/position model (get-in % [:tracks i])) (particles t))))) (range (inc k)))
                             "#7dffa0" 2.5))
                (when (>= shown 40) (restart!))
                (str "scan " (inc k) " of 20" (when (.-clutter controls) ", with false alarms and missed detections")
                     "\nfaint: the true paths; dots: blips, unlabeled; red: nearest neighbor; gold: most probable assignment; green: particles (their cloud and mean)")))
            (draw-citations! [_cw ch]
              (let [{:keys [samples shown]} @model
                    k (min shown (dec (count samples)))
                    seen (subvec samples 0 (inc k))
                    n (count citation-texts)
                    same (fn [i j] (/ (count (filter (fn [g] (some #(and (% i) (% j)) g)) seen)) (double (count seen))))
                    y0 110
                    cell (max 12 (min 36 (/ (- ch y0 (* n 34) 170) n)))
                    x0 60 m0 (+ y0 (* n 34) 40)
                    groups (samples k)]
                (doseq [[i [label authors title]] (map-indexed vector citation-texts)
                        :let [g (first (keep-indexed (fn [gi grp] (when (grp i) gi)) groups))]]
                  (text! (str (inc i) ". " authors " -- " title) 40 (+ y0 (* i 34) 14) (colors (mod g (count colors))) 11)
                  (text! label 40 (+ y0 (* i 34) 28) "#56627a" 10))
                (doseq [i (range n) j (range n)]
                  (let [p (same i j)]
                    (set! (.-fillStyle ctx) (str "rgba(125,255,160," (.toFixed p 2) ")"))
                    (.fillRect ctx (+ x0 (* j cell)) (+ m0 (* i cell)) (dec cell) (dec cell))))
                (doseq [i (range n)]
                  (text! (str (inc i)) (+ x0 (* i cell) (/ cell 2)) (- m0 4) "#9fb4d8" 10 "center")
                  (text! (str (inc i)) (- x0 6) (+ m0 (* i cell) (/ cell 2) 4) "#9fb4d8" 10 "right"))
                (text! "P(same paper)" (+ x0 (* n cell) 16) (+ m0 12) "#9fb4d8" 11)
                (when (>= shown (+ (count samples) 20)) (restart!))
                (str "Gibbs sweep " k "; colors: the current grouping into papers; the grid: P(same paper) over the sweeps so far"
                     "\n1 and 2 are the section's own two citations -- the same paper, it says -- though written differently")))
            (draw-seismic! [cw ch]
              (let [{:keys [events detections bulletin shown]} @model
                    r (* 0.32 (min cw (- ch 160)))
                    cx (* 0.3 cw) cy (+ 80 (* 0.5 (- ch 160)))
                    on (fn [a rr] [(+ cx (* rr (math/cos a))) (- cy (* rr (math/sin a)))])
                    k (min shown (count bulletin))
                    tx0 (* 0.58 cw) tx1 (- cw 40) ty0 100 ty1 (- ch 120)
                    tpx (fn [t] (+ tx0 (* (- tx1 tx0) (/ t 100.0))))
                    spy (fn [s] (+ ty0 (* (- ty1 ty0) (/ (+ s 0.5) 12.0))))]
                (set! (.-strokeStyle ctx) "#3a4458") (set! (.-lineWidth ctx) 2)
                (.beginPath ctx) (.arc ctx cx cy r 0 (* 2 math/PI)) (.stroke ctx)
                (doseq [[s a] (map-indexed vector stations)]
                  (let [[x y] (on a r)] (dot! x y 4 "#9fb4d8") (text! (str s) (first (on a (+ r 14))) (+ 4 (second (on a (+ r 14)))) "#56627a" 10 "center")))
                (doseq [{:keys [place]} events] (let [[x y] (on place r)] (dot! x y 9 "rgba(255,255,255,0.25)")))
                (doseq [[i e] (map-indexed vector (take k bulletin))]
                  (let [[x y] (on (:place e) r)] (dot! x y 10 "#ffffff") (dot! x y 8 (colors (inc i)))))
                (doseq [s (range 12) t (detections s)]
                  (let [made (some (fn [[i e]] (when (= t (get-in e [:detections s])) i)) (map-indexed vector (take k bulletin)))]
                    (line! (tpx t) (- (spy s) 6) (tpx t) (+ (spy s) 6) (if made (colors (inc made)) "#56627a") (if made 3 1))))
                (doseq [s (range 12)] (text! (str s) (- tx0 8) (+ (spy s) 4) "#56627a" 10 "right"))
                (text! "detections at each station over the window" tx0 (- ty0 14) "#9fb4d8" 11)
                (when (>= shown (+ (count bulletin) 12)) (restart!))
                (str (reduce + (map count (vals detections))) " detections, about nine in ten false alarms; the bulletin's events, found one at a time: "
                     (count (take k bulletin)) " of " (count bulletin)
                     "\nfaint circles: where the events really were; colored: the bulletin's, and the detections each explains")))
            (draw-reading! [cw _ch]
              (let [{:keys [word image trace independent english shown]} @model
                    k (min shown (dec (count trace)))
                    w (count (first image))
                    px (min 10 (/ (- cw 80) w))]
                (doseq [r (range reading/rows) c (range w)]
                  (let [v (max 0.0 (min 1.0 (get-in image [r c])))
                        g (math/round (* 255 (- 1.0 v)))]
                    (set! (.-fillStyle ctx) (str "rgb(" g "," g "," g ")"))
                    (.fillRect ctx (+ 40 (* c px)) (+ 100 (* r px)) px px)))
                (let [y (+ 120 (* reading/rows px))
                      {:keys [text noise]} (trace k)]
                  (text! (str "MCMC, sweep " k ": " (if (seq text) text "(empty)")) 40 (+ y 20) "#ffd479" 16)
                  (text! (str "its noise variance: " (.toFixed noise 2)) 40 (+ y 42) "#9fb4d8" 12)
                  (text! (str "exact, independent letters: " independent) 40 (+ y 76) "#ff8f70" 14)
                  (text! (str "exact, bigrams of English:  " english) 40 (+ y 100) "#7dffa0" 14))
                (when (>= shown (+ (count trace) 20)) (restart!))
                (str "the word \"" word "\" rendered and degraded by noise of variance " (.-noise controls)
                     "\nthe trace starts empty; each sweep resamples every letter, proposes a longer or shorter text, and resamples the noise")))
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
        (.onChange (.add gui controls "method" #js ["ratings" "tracking" "citations" "seismic" "reading"]) restart!)
        (.onChange (.add gui controls "clutter") restart!)
        (.onChange (.add gui controls "word" #js ["uncertainty" "reading" "letters" "program" "planet" "window"]) restart!)
        (.onChange (.add gui controls "noise" #js [0.1 0.3 0.6 1.0]) restart!)
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
                   (when-let [container (js/document.getElementById "open-universe")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
