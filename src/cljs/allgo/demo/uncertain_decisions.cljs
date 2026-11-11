(ns allgo.demo.uncertain-decisions
  "Canvas demo for `allgo.decision.nature` and `allgo.decision.assistance`:
  deciding when something you cannot see decides the outcome too -- nature,
  in LaValle's games against nature (*Planning Algorithms*, Cambridge
  University Press, 2006, section 9.2), or a person whose wishes you serve,
  in Russell and Norvig's assistance games (*Artificial Intelligence: A
  Modern Approach*, 4th edition, sections 16.7.2 and 18.2.5).

  Games against nature: a cost for each action against each of nature's
  choices, and the choice of action under four views of nature -- the
  worst case (minimax), an expected cost under a distribution you drag
  about the triangle (each point colored by the action it favors), the
  worst regret, and Hurwicz's blend of best and worst. Click a column's
  head to rule that choice of nature's out, as an observation would.

  The off-switch game: Robbie's prior over what his plan is worth to
  Harriet, and the worth of acting, switching off, and deferring to her.
  Deferring cuts off the downside, so it is never worth less than acting
  -- until her time costs something or she sometimes decides wrong.

  The paperclip game: Harriet makes a few paperclips or staples and Robbie
  then makes many. Round by round of myopic best response, her choices
  become a signal and his replies its reading. Click to pick a Harriet,
  by how much she values a paperclip."
  (:require [allgo.decision.assistance :as assist]
            [allgo.decision.nature :as nature]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:scene "games against nature" :seed 1 :alpha 0.5
       :lo -40 :hi 60 :cost 0 :error 0})

(def ^:private colors {:u1 "#7fd1ff" :u2 "#ffd479" :u3 "#ff8f70"})
(def ^:private actions [:u1 :u2 :u3])
(def ^:private natures [:t1 :t2 :t3])

(def ^:private lavalle
  ;; LaValle's Examples 9.8 to 9.10
  {:u1 {:t1 1 :t2 -1 :t3 0}
   :u2 {:t1 -1 :t2 2 :t3 -2}
   :u3 {:t1 2 :t2 -1 :t3 1}})

(defn- cost-table
  "LaValle's table for seed 1, otherwise costs from -3 to 3 at random."
  [seed]
  (if (= seed 1)
    lavalle
    (let [rng (random/rng seed)]
      (into {} (for [a actions] [a (into {} (for [t natures] [t (- (random/below rng 7) 3)]))])))))

(defn- best-expected
  "The action of least expected cost at barycentric `p` over the `possible`
  choices of nature's, renormalized."
  [costs possible p]
  (let [z (reduce + (map p possible))]
    (when (pos? z)
      (nature/expected (nature/observe costs possible) (into {} (map (fn [t] [t (/ (p t) z)]) possible))))))

(def ^:private choice-colors {[2 0] "#7fd1ff" [1 1] "#ffd479" [0 2] "#ff8f70"})
(defn- clips [[p s]]
  (cond (zero? s) (str p " paperclips") (zero? p) (str s " staples") :else (str p " + " s)))

(def ^:private paperclip-prior (assist/uniform-prior 0.0 1.0 1001))
(def ^:private paperclip-rounds
  (delay (vec (take 3 (assist/best-response-seq assist/paperclip-game paperclip-prior)))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        state (atom {:p {:t1 0.2 :t2 0.2 :t3 0.6} :ruled-out #{} :theta 0.45 :t 0.0 :dragging false :costs lavalle})]
    (set! (.-style canvas) "display:block;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(zoom []
              ;; drawn in a space at least 1100 by 720, scaled up to fill a larger canvas
              (max 1.0 (min (/ (.-width canvas) 1100) (/ (.-height canvas) 720))))
            (vw [] (/ (.-width canvas) (zoom)))
            (vh [] (/ (.-height canvas) (zoom)))
            (text! [s x y color & [font align]]
              (set! (.-fillStyle ctx) color)
              (set! (.-font ctx) (or font "12px ui-monospace,Menlo,monospace"))
              (set! (.-textAlign ctx) (or align "left"))
              (.fillText ctx s x y)
              (set! (.-textAlign ctx) "left"))
            (line! [ax ay bx by color width]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
              (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
            (fmt [x] (if (integer? x) (str x) (.toFixed x 2)))
            ;; ---- games against nature --------------------------------------
            (nature-layout []
              (let [cw (vw) ch (vh)
                    side (min (* 0.36 cw) (- ch 260))]
                {:table [40 70] :col 74 :row 30
                 ;; the simplex: t1 at the bottom left, t2 bottom right, t3 the top
                 :tri [[(- cw side 60) (+ 80 (* side 0.87))] [(- cw 60) (+ 80 (* side 0.87))] [(- cw 60 (* 0.5 side)) 80]]
                 :plot [40 (+ 70 (* 30 6)) (min 420 (- cw side 160)) 150]}))
            (bary->xy [[a b c] p]
              (let [{:keys [t1 t2 t3]} p]
                [(+ (* t1 (a 0)) (* t2 (b 0)) (* t3 (c 0))) (+ (* t1 (a 1)) (* t2 (b 1)) (* t3 (c 1)))]))
            (xy->bary [[[ax ay] [bx by] [cx cy]] [x y]]
              (let [d (+ (* (- by cy) (- ax cx)) (* (- cx bx) (- ay cy)))
                    l1 (/ (+ (* (- by cy) (- x cx)) (* (- cx bx) (- y cy))) d)
                    l2 (/ (+ (* (- cy ay) (- x cx)) (* (- ax cx) (- y cy))) d)
                    l1 (max 0.0 l1) l2 (max 0.0 l2) l3 (max 0.0 (- 1 l1 l2))
                    z (+ l1 l2 l3)]
                {:t1 (/ l1 z) :t2 (/ l2 z) :t3 (/ l3 z)}))
            (draw-nature! []
              (let [{:keys [costs p ruled-out]} @state
                    possible (vec (remove ruled-out natures))
                    seen (nature/observe costs possible)
                    alpha (.-alpha controls)
                    {[tx ty] :table col :col row :row tri :tri [px py pw ph] :plot} (nature-layout)
                    worst (when (seq possible) (nature/worst-case seen))
                    expected (best-expected costs possible p)
                    regret (when (seq possible) (nature/minimax-regret seen))
                    hurwicz (when (seq possible) (nature/hurwicz seen alpha))
                    regrets (when (seq possible) (nature/regrets seen))
                    z (reduce + (map p possible))
                    expect-of (fn [a] (reduce + (map (fn [t] (* (get-in costs [a t]) (/ (p t) z))) possible)))]
                ;; the table
                (text! "action" tx ty "#9fb4d8")
                (doseq [[k t] (map-indexed vector natures)]
                  (text! (str (name t) (when (ruled-out t) " (ruled out)")) (+ tx (* col (inc k))) ty (if (ruled-out t) "#5b6884" "#9fb4d8")))
                (doseq [[k h] (map-indexed vector ["worst" "expected" "regret" "Hurwicz"])]
                  (text! h (+ tx (* col (+ 4 k))) ty "#9fb4d8"))
                (doseq [[r a] (map-indexed vector actions)
                        :let [y (+ ty (* row (inc r)))
                              pick (fn [[w] k] (when (= w a) (set! (.-fillStyle ctx) "rgba(255,212,121,0.18)")
                                                     (.fillRect ctx (- (+ tx (* col k)) 6) (- y 18) (- col 6) 26)))]]
                  (pick worst 4) (pick expected 5) (pick regret 6) (pick hurwicz 7)
                  (text! (name a) tx y (colors a) "bold 13px ui-monospace,Menlo,monospace")
                  (doseq [[k t] (map-indexed vector natures)]
                    (text! (str (get-in costs [a t])) (+ tx (* col (inc k))) y (if (ruled-out t) "#3a4458" "#c9d2e0")))
                  (when (seq possible)
                    (let [cs (map #(get-in costs [a %]) possible)]
                      (text! (str (reduce max cs)) (+ tx (* col 4)) y "#c9d2e0")
                      (text! (fmt (expect-of a)) (+ tx (* col 5)) y "#c9d2e0")
                      (text! (str (reduce max (vals (regrets a)))) (+ tx (* col 6)) y "#c9d2e0")
                      (text! (fmt (+ (* alpha (reduce max cs)) (* (- 1 alpha) (reduce min cs)))) (+ tx (* col 7)) y "#c9d2e0"))))
                ;; Hurwicz against alpha
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx px py pw ph)
                (when (seq possible)
                  (let [vals* (for [a actions :let [cs (map #(get-in costs [a %]) possible)]] [a (reduce min cs) (reduce max cs)])
                        lo (reduce min (map second vals*)) hi (reduce max (map #(nth % 2) vals*))
                        span (max 1 (- hi lo))
                        y-of (fn [v] (- (+ py ph -10) (* (- ph 20) (/ (- v lo) span))))]
                    (doseq [[a mn mx] vals*]
                      (line! px (y-of mn) (+ px pw) (y-of mx) (colors a) 2))
                    (line! (+ px (* pw alpha)) py (+ px (* pw alpha)) (+ py ph) "rgba(255,255,255,0.6)" 1)
                    (text! "Hurwicz value against alpha: 0, the optimist, to 1, minimax" px (- py 8) "#9fb4d8" "11px ui-monospace,Menlo,monospace")))
                ;; the simplex, each cell colored by the expected-cost choice
                (let [n 48]
                  (doseq [i (range n) j (range (- n i))
                          :let [b {:t1 (/ (+ i 0.33) n) :t2 (/ (+ j 0.33) n) :t3 (- 1 (/ (+ i j 0.66) n))}
                                [a] (best-expected costs possible b)
                                corners [[i j] [(inc i) j] [i (inc j)]]]]
                    (when a
                      (set! (.-fillStyle ctx) (str (colors a) "66"))
                      (.beginPath ctx)
                      (doseq [[k [ci cj]] (map-indexed vector corners)
                              :let [[x y] (bary->xy tri {:t1 (/ ci n) :t2 (/ cj n) :t3 (- 1 (/ (+ ci cj) n))})]]
                        (if (zero? k) (.moveTo ctx x y) (.lineTo ctx x y)))
                      (.fill ctx)
                      (when (< (+ i j 1) n)
                        (.beginPath ctx)
                        (doseq [[k [ci cj]] (map-indexed vector [[(inc i) j] [(inc i) (inc j)] [i (inc j)]])
                                :let [[x y] (bary->xy tri {:t1 (/ ci n) :t2 (/ cj n) :t3 (- 1 (/ (+ ci cj) n))})]]
                          (if (zero? k) (.moveTo ctx x y) (.lineTo ctx x y)))
                        (.fill ctx))))
                  (let [[[ax ay] [bx by] [cx cy]] tri]
                    (text! "t1" (- ax 20) (+ ay 4) "#9fb4d8") (text! "t2" (+ bx 6) (+ by 4) "#9fb4d8") (text! "t3" (- cx 6) (- cy 8) "#9fb4d8")
                    (text! "nature's distribution: drag the dot" ax (+ ay 26) "#9fb4d8" "11px ui-monospace,Menlo,monospace"))
                  (let [[x y] (bary->xy tri p)]
                    (set! (.-fillStyle ctx) "#ffffff")
                    (.beginPath ctx) (.arc ctx x y 6 0 (* 2 math/PI)) (.fill ctx)))
                (str (if (seq possible)
                       (str "minimax " (name (first worst)) "  |  expected " (name (first expected))
                            " (P = " (apply str (interpose ", " (map #(.toFixed (/ (p %) z) 2) possible))) ")"
                            "  |  minimax regret " (name (first regret)) "  |  Hurwicz at " (.toFixed alpha 2) " " (name (first hurwicz)))
                       "every choice of nature's ruled out")
                     "\ngold: each criterion's choice  triangle: the expected-cost choice at every distribution  click a column head: rule it out")))
            ;; ---- the off-switch game ----------------------------------------
            (draw-off-switch! []
              (let [cw (vw) ch (vh)
                    lo (min (.-lo controls) (.-hi controls)) hi (max (.-lo controls) (.-hi controls) (+ 1 (min (.-lo controls) (.-hi controls))))
                    prior (assist/uniform-prior lo hi 2001)
                    cost (.-cost controls) error (.-error controls)
                    {:keys [act off defer best]} (assist/off-switch prior cost error)
                    ;; the prior over u, from -100 to 100
                    [x0 y0 w h] [40 70 (- (* 0.55 cw) 60) 120]
                    x-of (fn [u] (+ x0 (* w (/ (+ u 100) 200))))]
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x0 y0 w h)
                (when (neg? lo)
                  (set! (.-fillStyle ctx) "rgba(255,143,112,0.5)") (.fillRect ctx (x-of lo) (+ y0 30) (- (x-of (min 0 hi)) (x-of lo)) (- h 30)))
                (when (pos? hi)
                  (set! (.-fillStyle ctx) "rgba(255,212,121,0.6)") (.fillRect ctx (x-of (max 0 lo)) (+ y0 30) (- (x-of hi) (x-of (max 0 lo))) (- h 30)))
                (line! (x-of 0) y0 (x-of 0) (+ y0 h) "#ffffff" 1)
                (text! "Robbie's prior over u, what his plan is worth to Harriet (-100 to 100)" x0 (- y0 8) "#9fb4d8" "11px ui-monospace,Menlo,monospace")
                (text! "coral: she would switch him off" x0 (+ y0 h 18) "#ff8f70" "11px ui-monospace,Menlo,monospace")
                (text! "gold: she lets him go ahead" (+ x0 260) (+ y0 h 18) "#ffd479" "11px ui-monospace,Menlo,monospace")
                ;; the three values
                (let [bx (+ (* 0.55 cw) 20) by 70 bw (- (* 0.45 cw) 80) bh 160
                      top (max 1.0 (abs act) (abs defer) 10.0)
                      mid (+ by (* 0.5 bh))
                      y-of (fn [v] (- mid (* 0.45 bh (/ v top))))]
                  (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx bx by bw bh)
                  (line! bx mid (+ bx bw) mid "rgba(255,255,255,0.3)" 1)
                  (doseq [[k [lbl v]] (map-indexed vector [["act" act] ["off" off] ["defer" defer]])
                          :let [x (+ bx 20 (* k (/ (- bw 40) 3))) wbar (- (/ (- bw 40) 3) 20)
                                kw (keyword lbl)]]
                    (set! (.-fillStyle ctx) (if (= kw best) "#ffd479" "#5b6884"))
                    (.fillRect ctx x (min mid (y-of v)) wbar (max 1 (abs (- (y-of v) mid))))
                    (text! (str lbl " " (.toFixed v 1)) x (+ by bh 18) (if (= kw best) "#ffd479" "#c9d2e0")))
                  (text! "what each choice is worth to Harriet" bx (- by 8) "#9fb4d8" "11px ui-monospace,Menlo,monospace"))
                ;; defer's worth as she errs more
                (let [[x0 y0 w h] [40 (+ 70 190) (- cw 80) (- ch 70 190 90)]
                      errs (map #(* 0.5 (/ % 50.0)) (range 51))
                      ds (map #(:defer (assist/off-switch prior cost %)) errs)
                      vs (concat ds [act 0.0])
                      mn (reduce min vs) mx (reduce max vs) span (max 1e-9 (- mx mn))
                      y-of (fn [v] (- (+ y0 h -10) (* (- h 20) (/ (- v mn) span))))
                      x-of (fn [e] (+ x0 (* w (/ e 0.5))))]
                  (when (> h 60)
                    (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x0 y0 w h)
                    (line! x0 (y-of act) (+ x0 w) (y-of act) "#7fd1ff" 2)
                    (line! x0 (y-of 0.0) (+ x0 w) (y-of 0.0) "rgba(255,255,255,0.3)" 1)
                    (set! (.-strokeStyle ctx) "#ffd479") (set! (.-lineWidth ctx) 2)
                    (.beginPath ctx)
                    (doseq [[k e d] (map vector (range) errs ds)] (if (zero? k) (.moveTo ctx (x-of e) (y-of d)) (.lineTo ctx (x-of e) (y-of d))))
                    (.stroke ctx)
                    (line! (x-of error) y0 (x-of error) (+ y0 h) "rgba(255,255,255,0.6)" 1)
                    (text! "defer (gold) against act (blue) as Harriet's chance of deciding wrong goes from 0 to 0.5" x0 (- y0 8) "#9fb4d8" "11px ui-monospace,Menlo,monospace")))
                (str "act " (.toFixed act 1) ", off 0, defer " (.toFixed defer 1) ": Robbie should " (name best)
                     (when (and (zero? cost) (zero? error)) " -- deferring keeps the gold and loses the coral, never worse than acting")
                     "\nsliders: the prior's ends, the cost of Harriet's time, her chance of deciding wrong")))
            ;; ---- the paperclip game ------------------------------------------
            (draw-paperclips! []
              (let [cw (vw) ch (vh)
                    {:keys [theta t]} @state
                    rounds @paperclip-rounds
                    shown (min (count rounds) (inc (long (/ t 1.5))))
                    [x0 y0 w] [40 80 (- cw 80)]
                    rh (min 90 (/ (- ch 220) 3))
                    x-of (fn [th] (+ x0 (* w th)))
                    labels ["greedy Harriet, on her own -- and Robbie's reply" "Harriet's reply, knowing his -- and his to hers" "and again: nothing changes, an equilibrium"]]
                (doseq [k (range shown)
                        :let [{:keys [human robot]} (rounds k)
                              y (+ y0 (* k (+ rh 30)))]]
                  (text! (labels k) x0 (- y 8) "#9fb4d8" "11px ui-monospace,Menlo,monospace")
                  (doseq [choice (:human assist/paperclip-game)
                          :let [[a b] (assist/signal-range human choice)]
                          :when a]
                    (set! (.-fillStyle ctx) (str (choice-colors choice) "aa"))
                    (.fillRect ctx (x-of a) y (max 1 (- (x-of b) (x-of a))) rh)
                    (when (> (- (x-of b) (x-of a)) 70)
                      (text! (clips choice) (+ (x-of a) 6) (+ y 18) "#0d1118" "bold 12px ui-monospace,Menlo,monospace")
                      (when-let [r (robot choice)] (text! (str "Robbie: " (clips r)) (+ (x-of a) 6) (+ y 36) "#0d1118" "12px ui-monospace,Menlo,monospace")))
                    (when (pos? a) (text! (.toFixed a 3) (x-of a) (+ y rh 14) "#c9d2e0" "10px ui-monospace,Menlo,monospace" "center"))))
                (line! (x-of theta) (- y0 20) (x-of theta) (+ y0 (* 3 (+ rh 30))) "#ffffff" 1.5)
                (text! (str "theta " (.toFixed theta 2)) (+ (x-of theta) 4) (- y0 24) "#ffffff")
                (let [eq (peek rounds)
                      alone ((assist/greedy-human assist/paperclip-game paperclip-prior) theta)
                      mine ((:human eq) theta)
                      his ((:robot eq) mine)]
                  (str "a Harriet who values a paperclip at " (.toFixed theta 2) " (a staple at " (.toFixed (- 1 theta) 2) "): alone she makes " (clips alone)
                       "; at the equilibrium " (clips mine) " and Robbie " (clips his) ", worth $"
                       (.toFixed (+ (assist/payoff theta mine) (assist/payoff theta his)) 2)
                       "\ntheta, a paperclip's worth to Harriet, from 0 (left) to 1 (right)  click: pick a Harriet"))))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas) k (zoom)]
                (.setTransform ctx 1 0 0 1 0 0)
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (.setTransform ctx k 0 0 k 0 0)
                (set! (.-textContent caption)
                      (str (.-scene controls) "  |  "
                           (case (.-scene controls)
                             "games against nature" (draw-nature!)
                             "off-switch game" (draw-off-switch!)
                             "paperclip game" (draw-paperclips!))))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (swap! state update :t + (/ 1.0 60))
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))
            (at [^js e]
              (let [r (.getBoundingClientRect canvas) k (zoom)]
                [(/ (- (.-clientX e) (.-left r)) k) (/ (- (.-clientY e) (.-top r)) k)]))
            (near-tri? [[x y]]
              (let [{[[ax ay] [bx _] [_ cy]] :tri} (nature-layout)]
                (and (<= (- ax 10) x (+ bx 10)) (<= (- cy 10) y (+ ay 10)))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (let [[x y :as xy] (at e)]
                             (case (.-scene controls)
                               "games against nature"
                               (let [{[tx ty] :table col :col tri :tri} (nature-layout)]
                                 (cond
                                   (near-tri? xy) (swap! state assoc :dragging true :p (xy->bary tri xy))
                                   (< (- ty 20) y (+ ty 6))
                                   (let [k (dec (long (math/floor (/ (- x tx) col))))]
                                     (when (< -1 k 3) (swap! state update :ruled-out #(let [t (natures k)] (if (% t) (disj % t) (conj % t))))))
                                   :else nil))
                               "paperclip game"
                               (swap! state assoc :theta (max 0.0 (min 1.0 (/ (- x 40) (- (vw) 80)))))
                               nil))))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (when (and (:dragging @state) (= "games against nature" (.-scene controls)))
                             (swap! state assoc :p (xy->bary (:tri (nature-layout)) (at e))))))
      (doseq [ev ["pointerup" "pointercancel" "pointerleave"]]
        (.addEventListener canvas ev (fn [_] (swap! state assoc :dragging false))))
      (let [gui (GUI. #js {:container container})
            scene (.add gui controls "scene" #js ["games against nature" "off-switch game" "paperclip game"])
            ;; each scene's own sliders, the rest hidden
            own {"games against nature" [(.onChange (.add gui controls "seed" 1 50 1)
                                                    #(swap! state assoc :costs (cost-table (.-seed controls)) :ruled-out #{}))
                                         (.add gui controls "alpha" 0 1 0.01)]
                 "off-switch game" [(.add gui controls "lo" -100 100 1) (.add gui controls "hi" -100 100 1)
                                    (.add gui controls "cost" 0 30 0.5) (.add gui controls "error" 0 0.5 0.01)]
                 "paperclip game" [(.add gui #js {:replay #(swap! state assoc :t 0.0)} "replay")]}
            show! (fn [] (doseq [[k cs] own ^js c cs] (.show c (= k (.-scene controls)))))]
        (.onChange scene (fn [_] (swap! state assoc :t 0.0) (show!)))
        (show!))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "uncertain-decisions")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
