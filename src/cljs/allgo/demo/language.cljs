(ns allgo.demo.language
  "Canvas demo for `allgo.language`: the deep learning for language of
  Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
  edition, chapter 24, each model training as you watch -- a little each
  frame, so the page never stalls.

  Embeddings: word vectors learned by skip-gram, drawn on their two
  principal components; as they settle, man -> woman runs parallel to
  king -> queen; or by GloVe, fitting the table of how often words occur
  together, an epoch a step. Language model: a recurrent network reading grammar
  sentences, its perplexity falling and its own sentences growing
  grammatical. Attention: an encoder-decoder learning to reverse digits,
  its attention weights a grid that turns into the anti-diagonal. Tagging:
  a transformer learning parts of speech where a word can be two -- the
  duck will duck -- and one head's attention over the sentence.
  Bidirectional: two recurrent taggers learning sentences whose first
  word is a noun or a verb by the word after it, \"duck swim today\" and
  \"duck stones today\": the one reading left to right can only guess,
  the one reading both ways learns it. N-grams: chapter 23's unigram,
  bigram and trigram models of grammar sentences, their samples and
  perplexities as EM fits the trigram's interpolation weights, and an HMM
  tagging both ducks. Inside-outside: EM learning the grammar's rule
  probabilities from unparsed sentences, the likelihood climbing and the
  probabilities settling on the true ones."
  (:require [allgo.language.embeddings :as emb]
            [allgo.language.ngrams :as ng]
            [allgo.language.parsing :as parsing]
            [allgo.language.recurrent :as rec]
            [allgo.language.seq2seq :as s2s]
            [allgo.language.tagging :as tagging]
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

(def ^:private right-grammar
  "The first word a noun before an intransitive verb, a verb before an
  object: only the word after it decides."
  {:rules [[:S :A :Adv 1.0] [:A :N :Vi 0.5] [:A :V :Obj 0.5]]
   :lexicon {"duck" [[:N 0.2] [:V 0.2]] "fish" [[:N 0.2] [:V 0.2]] "watch" [[:N 0.2] [:V 0.2]] "saw" [[:N 0.2] [:V 0.2]]
             "dogs" [[:N 0.2]] "feed" [[:V 0.2]]
             "swim" [[:Vi 0.34]] "sleep" [[:Vi 0.33]] "run" [[:Vi 0.33]]
             "stones" [[:Obj 0.34]] "boats" [[:Obj 0.33]] "bread" [[:Obj 0.33]]
             "today" [[:Adv 0.34]] "again" [[:Adv 0.33]] "here" [[:Adv 0.33]]}})

(def ^:private right-probes [["duck" "swim" "today"] ["duck" "stones" "today"] ["fish" "sleep" "here"] ["fish" "bread" "here"]])

(defn- held
  "The first `n` of `xs`, then the last of those for ever: a run of a
  converging algorithm that the demo shows and then holds."
  [n xs]
  (let [s (take n xs)] (concat s (lazy-seq (repeat (last s))))))

(defn- uniform-grammar
  "The ambiguous grammar's rules, each left-hand side's equally likely."
  []
  (let [g parsing/ambiguous-grammar
        weights (fn [ks] (let [totals (frequencies (map first ks))] (fn [x] (/ 1.0 (totals x)))))
        bin (map (fn [[x y z]] [x y z]) (:rules g))
        lex (for [[w xs] (:lexicon g) [x] xs] [x w])
        wb (weights bin) wl (weights lex)]
    {:rules (mapv (fn [[x y z]] [x y z (wb x)]) bin)
     :lexicon (reduce (fn [m [x w]] (update m w (fnil conj []) [x (wl x)])) {} lex)}))

(defn- first-ambiguous [w i] (and (zero? i) (#{"duck" "fish" "watch" "saw"} w)))

(defn- digits [rng] (vec (repeatedly (+ 3 (random/below rng 3)) #(random/below rng 5))))

(defn- init [method]
  (let [rng (random/rng (.-seed controls))]
    (case method
      "embeddings" {:models (emb/sgns-seq (royal-corpus rng 1500) rng :dim 12 :window 2 :chunk 100) :steps 0}
      "GloVe" {:models (emb/glove-seq (royal-corpus rng 2000) rng :dim 16 :window 4 :chunk 20) :steps 0 :glove true}
      "n-grams" (let [tagged (vec (repeatedly 2000 #(parsing/generate parsing/ambiguous-grammar rng)))
                      sents (mapv #(mapv first %) tagged)
                      [train held-out test] [(subvec sents 0 1500) (subvec sents 1500 1800) (subvec sents 1800)]
                      ms (mapv #(ng/model % train) [1 2 3])
                      fit (fn [m] (:lambdas (nth (ng/interpolation-seq m held-out) 20)))]
                  {:models (held 25 (mapcat #(repeat 15 %) (ng/interpolation-seq (ms 2) held-out)))
                   :ms ms :lambdas [(fit (ms 0)) (fit (ms 1))] :test test :rng rng
                   :tagger (tagging/train (subvec tagged 0 1500)) :steps 0 :paced true})
      "inside-outside" (let [sents (vec (repeatedly 150 #(parsing/leaves (parsing/generate-tree parsing/ambiguous-grammar rng))))]
                         {:models (held 360 (mapcat #(repeat 30 %) (parsing/inside-outside-seq (uniform-grammar) sents))) :steps 0 :paced true :history []})
      "bidirectional" (let [data (vec (repeatedly 300 #(parsing/generate right-grammar rng)))
                            vocab (rec/vocabulary (map #(mapv first %) data))
                            tagger (fn [bi] (rec/tagger :lstm vocab [:N :V :Vi :Obj :Adv] 8 12 1 bi (random/rng (inc (.-seed controls)))))]
                        {:pair [(rec/train-tagger-seq (tagger false) (subvec data 0 200))
                                (rec/train-tagger-seq (tagger true) (subvec data 0 200))]
                         :test (subvec data 200 230) :steps 0})
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
              ;; train for about a dozen milliseconds a frame -- or, paced,
              ;; a step a frame
              (if (:paced @model)
                (swap! model (fn [m] (-> m (update :models rest) (update :steps inc))))
                (let [t0 (js/performance.now) budget 12]
                  (loop []
                    (swap! model (fn [{:keys [models states pair] :as m}]
                                   (cond
                                     models (assoc m :models (rest models) :steps (inc (:steps m)))
                                   ;; the two taggers take turns
                                     pair (-> m (update-in [:pair (mod (:steps m) 2)] rest) (update :steps inc))
                                     :else
                                     (let [st (second states)]
                                       (-> m (assoc :states (rest states) :steps (inc (:steps m)))
                                           (update :losses #(let [l (conj % (:loss st))] (if (> (count l) 400) (subvec l 1) l))))))))
                    (when (< (- (js/performance.now) t0) budget) (recur))))))
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
                             ("embeddings" "GloVe") (draw-embeddings! cw ch)
                             "bidirectional" (draw-bidirectional! cw ch)
                             "n-grams" (draw-ngrams! cw ch)
                             "inside-outside" (draw-inside-outside! cw ch)
                             "language model" (draw-lm! cw ch)
                             "attention" (draw-attention! cw ch)
                             "tagging" (draw-tagging! cw ch))))))
            (draw-embeddings! [cw ch]
              (let [{:keys [models steps glove]} @model
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
                (str (if glove (str (* 20 steps) " entries of the co-occurrence table fitted") (str (* 100 steps) " word pairs seen"))
                     "  |  king - man + woman = " (emb/analogy m "man" "woman" "king")
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
            (draw-ngrams! [_cw _ch]
              (let [{:keys [models ms lambdas test rng steps tagger]} @model
                    ls (conj lambdas (:lambdas (first models)))
                    ppl (mapv #(ng/perplexity (ms %) [:interpolate (ls %)] test) [0 1 2])]
                (when (zero? (mod steps 90))
                  (swap! model assoc :samples (mapv #(ng/sample (ms %) [:interpolate (ls %)] rng :max-len 12) [0 1 2])))
                (doseq [[i label] (map-indexed vector ["unigram" "bigram" "trigram"])]
                  (let [y (+ 140 (* 70 i))]
                    (text! (str label "  perplexity " (.toFixed (ppl i) 2)) 40 y "#9fb4d8" 13)
                    (when-let [s (get (:samples @model) i)]
                      (text! (str/join " " s) 40 (+ y 24) (if (parsing/parse-tree s parsing/ambiguous-grammar) "#7dffa0" "#ff8f70") 14))))
                (text! (str "the trigram's weights, uniform to trigram: " (str/join "  " (map #(.toFixed % 3) (ls 2)))) 40 370 "#ffd479" 13)
                (text! "an HMM tagging, Viterbi over the whole sentence:" 40 420 "#9fb4d8" 13)
                (doseq [[i [w t]] (map-indexed vector (map vector probe-sentence (tagging/tag tagger probe-sentence)))]
                  (text! w (+ 40 (* 110 i)) 446 "#e6edf7" 15)
                  (text! (name t) (+ 40 (* 110 i)) 468 (if (= t ([:Det :N :Modal :V :Det :N] i)) "#7dffa0" "#ff8f70") 13))
                (str "EM iteration " (quot steps 15) " on held-out sentences"
                     "\nsamples from each model, green where they parse under the grammar; perplexity on unseen sentences, lower the better")))
            (draw-inside-outside! [_cw _ch]
              (let [{:keys [models steps history]} @model
                    {:keys [grammar log-likelihood]} (first models)
                    ;; an iteration every thirty frames, so that it can be watched
                    iteration (quot (dec steps) 30)
                    history (if (and (<= (count history) iteration) (< iteration 12)) (conj history log-likelihood) history)
                    _ (swap! model assoc :history history)
                    truth parsing/ambiguous-grammar
                    p-of (fn [g [x y z]] (or (some (fn [[a b c p]] (when (= [a b c] [x y z]) p)) (:rules g)) 0.0))]
                (text! "log likelihood of the sentences, iteration by iteration -- EM never lets it fall:" 40 140 "#9fb4d8" 13)
                (text! (str/join "   " (map #(.toFixed % 1) (take 8 history))) 40 166 "#ffd479" 14)
                (text! "rule                  learned   true" 40 280 "#9fb4d8" 12)
                (doseq [[i [x y z p]] (map-indexed vector (:rules truth))]
                  (let [y0 (+ 300 (* 22 i)) q (p-of grammar [x y z])]
                    (text! (str (name x) " -> " (name y) " " (name z)) 40 y0 "#e6edf7" 13)
                    (set! (.-fillStyle ctx) "#7dffa0") (.fillRect ctx 240 (- y0 10) (* 200 q) 8)
                    (set! (.-strokeStyle ctx) "#ffd479") (.strokeRect ctx 240 (- y0 10) (* 200 p) 8)
                    (text! (.toFixed q 2) 460 y0 "#7dffa0" 12)))
                (str "iteration " (min iteration 11) "  |  log likelihood of the sentences " (.toFixed log-likelihood 1)
                     "\nfrom every rule equally likely, EM on unparsed sentences -- the rules known, one iteration nearly finds them; bars: each rule's learned probability, green, against its true one, outlined")))
            (draw-bidirectional! [_cw _ch]
              (let [{:keys [pair test steps]} @model
                    [fwd bi] (map (comp :model first) pair)]
                ;; the held-out sentences and the four shown, tagged by both a
                ;; couple a frame, every fifty sentences trained
                (let [{:keys [checking]} @model]
                  (cond (and (nil? checking) (zero? (mod steps 100)))
                        (swap! model assoc :checking {:models [fwd bi] :left (concat right-probes (map #(mapv first %) test)) :done []})
                        checking
                        (let [{:keys [models left done]} checking
                              done (into done (for [s (take 2 left)] (mapv #(rec/tag % s) models)))]
                          (if (seq (drop 2 left))
                            (swap! model assoc-in [:checking] (assoc checking :left (drop 2 left) :done done))
                            (let [[shown held] (split-at (count right-probes) done)
                                  right (fn [k] (count (filter true? (for [[got s] (map vector held test)
                                                                           :when (first-ambiguous (ffirst s) 0)]
                                                                       (= (second (first s)) (first (get got k)))))))
                                  n (count (filter #(first-ambiguous (ffirst %) 0) test))]
                              (swap! model assoc :checking nil
                                     :scores (mapv #(/ (right %) (max 1 n)) [0 1])
                                     :probes (mapv (fn [k] (mapv #(get % k) shown)) [0 1])))))))
                (let [{:keys [scores probes]} @model]
                  (doseq [[col label ps] (map vector [0 1] ["left to right" "both ways"] (or probes [nil nil]))]
                    (let [x (+ 40 (* col 420))]
                      (text! label x 140 "#9fb4d8" 13)
                      (when scores (text! (str (.toFixed (* 100 (scores col)) 0) "% of first words right") x 162 "#e6edf7" 13))
                      (doseq [[i [s tags]] (map-indexed vector (map vector right-probes (or ps (repeat nil))))]
                        (let [y (+ 210 (* 56 i)) truth (if (even? i) [:N :Vi :Adv] [:V :Obj :Adv])]
                          (doseq [[j w] (map-indexed vector s)]
                            (text! w (+ x (* 110 j)) y "#e6edf7" 15)
                            (when tags (text! (name (tags j)) (+ x (* 110 j)) (+ y 20) (if (= (truth j) (tags j)) "#7dffa0" "#ff8f70") 13)))))))
                  (str (quot steps 2) " sentences each"
                       "\nthe first word is a noun before \"swim\", a verb before \"stones\": reading left to right, it is tagged before the word that decides"))))
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
        (.onChange (.add gui controls "method" #js ["embeddings" "GloVe" "language model" "attention" "tagging" "bidirectional" "n-grams" "inside-outside"]) restart!)
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
