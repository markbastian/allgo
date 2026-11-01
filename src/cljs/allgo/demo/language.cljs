(ns allgo.demo.language
  "Canvas demo for `allgo.language`: the deep learning for language of
  Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
  edition, chapter 24, each model training as you watch -- a little each
  frame, so the page never stalls.

  Embeddings: word vectors learned by skip-gram, drawn on their two
  principal components; as they settle, man -> woman runs parallel to
  king -> queen. Language model: a recurrent network reading grammar
  sentences, its perplexity falling and its own sentences growing
  grammatical. Attention: an encoder-decoder learning to reverse digits,
  its attention weights a grid that turns into the anti-diagonal. Tagging:
  a transformer learning parts of speech where a word can be two -- the
  duck will duck -- and one head's attention over the sentence."
  (:require [allgo.language.embeddings :as emb]
            [allgo.language.parsing :as parsing]
            [allgo.language.recurrent :as rec]
            [allgo.language.seq2seq :as s2s]
            [allgo.language.transformer :as tf]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.string :as str]))

(def ^:private ^js controls #js {:method "embeddings" :cell "lstm" :seed 1})

;; ---------------------------------------------------------------------------
;; Data

(defn- royal-corpus [rng n]
  (let [people {"king" [:m :royal :old] "queen" [:f :royal :old] "prince" [:m :royal :young] "princess" [:f :royal :young]
                "man" [:m :common :old] "woman" [:f :common :old] "boy" [:m :common :young] "girl" [:f :common :young]}
        verbs {:royal ["rules" "crowns" "reigns"] :common ["works" "plows" "sells"]}
        things {:royal ["throne" "crown" "kingdom"] :common ["field" "cart" "goods"]}
        pron {:m ["his" "he"] :f ["her" "she"]}
        age {:old ["grown" "elder"] :young ["little" "small"]}]
    (vec (repeatedly n (fn [] (let [[w [g r a]] (random/pick rng (vec people))]
                                [(random/pick rng (age a)) w (random/pick rng (verbs r)) (random/pick rng (pron g)) (random/pick rng (things r))]))))))

(def ^:private shown-words ["king" "queen" "man" "woman" "prince" "princess" "boy" "girl" "rules" "works" "throne" "field"])
(def ^:private probe-sentence ["the" "duck" "will" "duck" "a" "dog"])
(def ^:private probe-digits [0 1 2 3 4])

(defn- digits [rng] (vec (repeatedly (+ 3 (random/below rng 3)) #(random/below rng 5))))

(defn- init [method]
  (let [rng (random/rng (.-seed controls))]
    (case method
      "embeddings" {:models (emb/sgns-seq (royal-corpus rng 1500) rng :dim 12 :window 2 :chunk 100) :steps 0}
      "language model" (let [sents (vec (repeatedly 400 #(mapv first (parsing/generate parsing/ambiguous-grammar rng))))]
                         {:states (rec/train-seq (rec/model (keyword (.-cell controls)) (rec/vocabulary sents) 8 16 rng) (subvec sents 0 300))
                          :test (subvec sents 300) :rng rng :steps 0 :losses []})
      "attention" (let [pairs (vec (repeatedly 500 (fn [] (let [d (digits rng)] [d (vec (reverse d))]))))]
                    {:states (s2s/train-seq (s2s/model :rnn (s2s/vocabulary (range 5)) 8 16 rng) pairs) :steps 0 :losses []})
      "tagging" (let [tagged (vec (repeatedly 400 #(parsing/generate parsing/ambiguous-grammar rng)))
                      vocab (tf/vocabulary (map #(mapv first %) tagged))]
                  {:states (tf/train-tagger-seq (tf/tagger (tf/encoder vocab 16 2 32 1 rng) [:Det :N :V :Adj :Modal :Adv] rng) tagged)
                   :steps 0 :losses []}))))

;; ---------------------------------------------------------------------------

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
    (letfn [(restart! [] (reset! model (init (.-method controls))))
            (advance! []
              ;; train for about a dozen milliseconds a frame
              (let [t0 (js/performance.now) budget 12]
                (loop []
                  (swap! model (fn [{:keys [models states] :as m}]
                                 (if models
                                   (assoc m :models (rest models) :steps (inc (:steps m)))
                                   (let [st (second states)]
                                     (-> m (assoc :states (rest states) :steps (inc (:steps m)))
                                         (update :losses #(let [l (conj % (:loss st))] (if (> (count l) 400) (subvec l 1) l))))))))
                  (when (< (- (js/performance.now) t0) budget) (recur)))))
            (text! [s x y color & [size]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (.fillText ctx s x y))
            (loss-curve! [losses x y w h]
              (when (> (count losses) 1)
                (let [sm (mapv (fn [i] (let [ws (subvec losses (max 0 (- i 20)) (inc i))] (/ (reduce + ws) (count ws)))) (range (count losses)))
                      hi (apply max sm) lo (min 0 (apply min sm))]
                  (set! (.-strokeStyle ctx) "#ffd479") (set! (.-lineWidth ctx) 1.5)
                  (.beginPath ctx)
                  (doseq [[i v] (map-indexed vector sm)]
                    (let [px (+ x (* w (/ i (dec (count sm))))) py (- (+ y h) (* h (/ (- v lo) (max 1e-9 (- hi lo)))))]
                      (if (zero? i) (.moveTo ctx px py) (.lineTo ctx px py))))
                  (.stroke ctx)
                  (text! "training loss" x (- y 6) "#9fb4d8" 11))))
            (heatmap! [rows x y cell row-labels col-labels]
              (doseq [[i r] (map-indexed vector rows) [j v] (map-indexed vector r)]
                (set! (.-fillStyle ctx) (str "rgba(255,190,90," (max 0.04 v) ")"))
                (.fillRect ctx (+ x (* j cell)) (+ y (* i cell)) (- cell 2) (- cell 2)))
              (doseq [[i l] (map-indexed vector row-labels)] (text! (str l) (- x 70) (+ y (* i cell) (* 0.6 cell)) "#c9d2e0"))
              (doseq [[j l] (map-indexed vector col-labels)] (text! (str l) (+ x (* j cell) 4) (- y 8) "#c9d2e0")))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas)
                    method (.-method controls)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "embeddings" (draw-embeddings! cw ch)
                             "language model" (draw-lm! cw ch)
                             "attention" (draw-attention! cw ch)
                             "tagging" (draw-tagging! cw ch))))))
            (draw-embeddings! [cw ch]
              (let [{:keys [models steps]} @model
                    m (first models)
                    pts (emb/project-2d (mapv #(emb/word-vector m %) shown-words))
                    xs (map first pts) ys (map second pts)
                    [x0 x1 y0 y1] [(apply min xs) (apply max xs) (apply min ys) (apply max ys)]
                    pad 80 w (- cw (* 2 pad)) h (- ch 220)
                    at (zipmap shown-words (map (fn [[x y]] [(+ pad (* w (/ (- x x0) (max 1e-9 (- x1 x0)))))
                                                             (+ 110 (* h (/ (- y y0) (max 1e-9 (- y1 y0)))))])
                                                pts))
                    arrow! (fn [a b color]
                             (let [[ax ay] (at a) [bx by] (at b)]
                               (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) 2)
                               (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx)))]
                (arrow! "man" "woman" "rgba(125,196,255,0.7)") (arrow! "king" "queen" "rgba(125,196,255,0.7)")
                (arrow! "boy" "girl" "rgba(125,196,255,0.7)") (arrow! "prince" "princess" "rgba(125,196,255,0.7)")
                (doseq [[w [x y]] at]
                  (set! (.-fillStyle ctx) "#ffd479") (.fillRect ctx (- x 3) (- y 3) 6 6)
                  (text! w (+ x 6) (- y 6) "#e6edf7" 13))
                (str (* 100 steps) " word pairs seen  |  king - man + woman = " (emb/analogy m "man" "woman" "king")
                     "\nwords on their vectors' two principal components; blue: male to female, parallel when the vectors learn the factor")))
            (draw-lm! [cw _ch]
              (let [{:keys [states test steps losses rng]} @model
                    m (:model (first states))]
                (when (zero? (mod steps 40))
                  (swap! model assoc
                         :samples (vec (repeatedly 6 #(rec/sample m rng)))
                         :ppl (rec/perplexity m (take 40 test))))
                (let [{:keys [samples ppl]} @model]
                  (loss-curve! losses 40 140 (- cw 80) 140)
                  (text! "its own sentences (green: they parse under the grammar)" 40 330 "#9fb4d8" 11)
                  (doseq [[i s] (map-indexed vector samples)]
                    (text! (str/join " " s) 40 (+ 356 (* 22 i)) (if (parsing/parse-tree s parsing/ambiguous-grammar) "#7dffa0" "#ff8f70") 14))
                  (str (name (:cell m)) ", " steps " sentences read"
                       (when ppl (str "  |  perplexity on unseen sentences " (.toFixed ppl 2)))
                       "\na recurrent language model reading sentences from a grammar, writing its own"))))
            (draw-attention! [cw ch]
              (let [{:keys [states steps losses]} @model
                    m (:model (first states))]
                (when (zero? (mod steps 10))
                  (swap! model assoc :probe (s2s/greedy m probe-digits)))
                (let [{:keys [target attention]} (:probe @model)
                      cell (min 60 (/ (- ch 260) 7))]
                  (loss-curve! losses 40 140 (- cw 80) 80)
                  (when attention
                    (heatmap! (vec (take 6 attention)) 160 290 cell (concat target ["</s>"]) probe-digits))
                  (str steps " pairs seen  |  " (str/join " " probe-digits) " -> " (str/join " " target)
                       "\nrows: each output token; columns: the source; bright where the decoder attends -- reversing, the anti-diagonal"))))
            (draw-tagging! [cw ch]
              (let [{:keys [states steps losses]} @model
                    m (:model (first states))]
                (when (zero? (mod steps 10))
                  (swap! model assoc :probe (tf/tag m probe-sentence)))
                (let [{:keys [tags attention]} (:probe @model)
                      truth [:Det :N :Modal :V :Det :N]
                      cell (min 60 (/ (- ch 320) 7))]
                  (loss-curve! losses 40 140 (- cw 80) 80)
                  (doseq [[i [w t g]] (map-indexed vector (map vector probe-sentence truth tags))]
                    (text! w (+ 40 (* 110 i)) 280 "#e6edf7" 15)
                    (text! (name (or g :?)) (+ 40 (* 110 i)) 302 (if (= t g) "#7dffa0" "#ff8f70") 13))
                  (when attention
                    (heatmap! (get-in attention [0 0]) 160 360 cell probe-sentence probe-sentence))
                  (str steps " sentences seen"
                       "\ngreen: tags right -- the two ducks differ only by context; the grid: one head's attention, each row a word looking at the others"))))
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
        (.onChange (.add gui controls "method" #js ["embeddings" "language model" "attention" "tagging"]) restart!)
        (.onChange (.add gui controls "cell" #js ["lstm" "rnn"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "language")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
