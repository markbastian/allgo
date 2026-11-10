(ns allgo.demo.generative
  "Canvas demo for `allgo.learning.generative`: the unsupervised deep
  models of Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, section 21.7.1, learning two-dimensional data a
  step at a time.

  Probabilistic PCA fits a Gaussian, so its samples fill an ellipse
  whatever the data's shape. The linear autoencoder squeezes each point
  to one number and back, onto a line: the first principal component.
  The deep autoencoder squeezes to one number too, but its decoder bends
  that line into a curve through the data. The variational autoencoder
  and the GAN each learn to turn Gaussian noise into new points like the
  data -- the GAN's discriminator drawn behind, brighter where it calls a
  point real. Autoregressive: the logistic map, a chaotic series, and the
  next value predicted from the two before by the linear AR model the
  Yule-Walker equations fit, and by a deep one."
  (:require [allgo.learning.generative :as gen]
            [allgo.learning.neural :as neural]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:model "variational autoencoder" :data "ring" :points 400 :seed 3 :speed 4})

(defn- dataset [kind n rng]
  (let [jitter #(random/gaussian rng 0.0 %)]
    (vec (repeatedly n (fn []
                         (case kind
                           "ring" (let [a (* 2 math/PI (rng))] [(+ (* 1.5 (math/cos a)) (jitter 0.06)) (+ (* 1.5 (math/sin a)) (jitter 0.06))])
                           "moon" (let [a (* math/PI (rng))] [(+ (* 1.6 (math/cos a)) (jitter 0.08)) (+ -0.6 (* 1.6 (math/sin a)) (jitter 0.08))])
                           "two blobs" (let [c (if (< (rng) 0.5) -1 1)] [(+ c (jitter 0.3)) (+ (* 0.6 c) (jitter 0.3))])
                           "spiral" (let [t (+ 0.5 (* 3.5 (rng))) r (* 0.5 t)] [(+ (* r (math/cos (* 1.4 t))) (jitter 0.05)) (+ (* r (math/sin (* 1.4 t))) (jitter 0.05))])))))))

(def ^:private series (vec (take 300 (iterate #(* 3.9 % (- 1 %)) 0.3))))

(defn- vae-round
  "A step of the VAE on a minibatch of 32 drawn from `X`, a frame each."
  [X rng m]
  (gen/vae-step m (vec (repeatedly 32 (fn [] (random/pick rng X)))) rng 0.005))

(defn- trainer
  "The model's states, a seq, a step each."
  [model X rng]
  (case model
    "probabilistic PCA" (gen/ppca-em-seq X 2 rng)
    "linear autoencoder" (gen/linear-autoencoder-seq X 1 rng :lr 0.02)
    "deep autoencoder" (gen/autoencoder-seq (gen/autoencoder 2 24 1 rng) X {:alpha 0.02 :rng rng})
    "variational autoencoder" (iterate (partial vae-round X rng) (gen/vae 2 32 2 32 0.01 rng))
    "GAN" (gen/gan-seq (gen/gan 2 2 32 rng) X rng :lr 0.002)
    "autoregressive" (let [[Xs Ys] (gen/windows series 2)]
                       (map (fn [net] {:net net :linear (gen/yule-walker series 2)})
                            (neural/train-seq [(neural/dense 2 16 :tanh rng) (neural/dense 16 1 :identity rng)] Xs Ys
                                              {:alpha 0.05 :rng rng})))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        region (js/document.createElement "canvas")
        rctx (.getContext region "2d")
        running? (atom false)
        model (atom nil)
        res 48]
    (set! (.-width region) res) (set! (.-height region) res)
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls))
                    X (dataset (.-data controls) (.-points controls) rng)
                    states (trainer (.-model controls) X rng)]
                (reset! model {:X X :rest states :st (first states) :step 0})))
            (advance! []
              (let [n (.-speed controls)]
                (swap! model (fn [{:keys [rest step] :as m}]
                               (let [rest (nthrest rest n)]
                                 (assoc m :rest rest :st (first rest) :step (+ step n)))))))
            (dot! [[px py] r color]
              (set! (.-fillStyle ctx) color)
              (.beginPath ctx) (.arc ctx px py r 0 (* 2 math/PI)) (.fill ctx))
            (line! [[ax ay] [bx by] color]
              (set! (.-strokeStyle ctx) color)
              (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
            (draw-plane! [cw ch]
              (let [{:keys [X st step]} @model
                    kind (.-model controls)
                    side (min (- cw 40) (- ch 110))
                    x0 (* 0.5 (- cw side)) y0 40
                    to-px (fn [[x y]] [(+ x0 (* side (/ (+ x 2.5) 5))) (+ y0 (* side (/ (- 2.5 y) 5)))])
                    ;; the same noise every frame, so the samples move as the model learns
                    fixed (random/rng 99)]
                (set! (.-fillStyle ctx) "#151b26")
                (.fillRect ctx x0 y0 side side)
                (.save ctx)
                (.beginPath ctx) (.rect ctx x0 y0 side side) (.clip ctx)
                (when (= kind "GAN")
                  (let [grid (vec (for [j (range res) i (range res)] [(- (* 5 (/ (+ i 0.5) res)) 2.5) (- 2.5 (* 5 (/ (+ j 0.5) res)))]))
                        ps (gen/discriminate st grid)
                        img (.createImageData rctx res res)
                        data (.-data img)]
                    (dotimes [k (* res res)]
                      (let [p (nth ps k) o (* 4 k)]
                        (aset data o (int (+ 21 (* 90 p)))) (aset data (+ o 1) (int (+ 27 (* 60 p))))
                        (aset data (+ o 2) (int (+ 38 (* 30 p)))) (aset data (+ o 3) 255)))
                    (.putImageData rctx img 0 0)
                    (set! (.-imageSmoothingEnabled ctx) true)
                    (.drawImage ctx region x0 y0 side side)))
                (doseq [x X] (dot! (to-px x) 2.5 "#7dc4ff"))
                (set! (.-lineWidth ctx) 1)
                (case kind
                  "probabilistic PCA"
                  (let [C (lin/mat-add (lin/mat-mul (:W st) (lin/transpose (:W st))) (lin/mat-scale (lin/eye 2) (:var st)))]
                    (doseq [x (repeatedly 300 #(gen/ppca-sample st fixed))] (dot! (to-px x) 2 "#ff8f70"))
                    ;; the one and two standard deviation ellipses
                    (when-let [L (lin/cholesky C)]
                      (doseq [k [1 2]]
                        (set! (.-strokeStyle ctx) "#ffd479")
                        (.beginPath ctx)
                        (doseq [i (range 65) :let [a (* 2 math/PI (/ i 64))
                                                   [px py] (to-px (lin/add (:mean st) (lin/scale (lin/mat-vec L [(math/cos a) (math/sin a)]) k)))]]
                          (if (zero? i) (.moveTo ctx px py) (.lineTo ctx px py)))
                        (.stroke ctx))))
                  "linear autoencoder"
                  (do (doseq [x X] (line! (to-px x) (to-px (gen/linear-decode st (gen/linear-encode st x))) "rgba(255,143,112,0.35)"))
                      (line! (to-px (gen/linear-decode st [-3])) (to-px (gen/linear-decode st [3])) "#ffd479")
                      (doseq [x X] (dot! (to-px (gen/linear-decode st (gen/linear-encode st x))) 2 "#ff8f70")))
                  "deep autoencoder"
                  (let [zs (map #(first (gen/encode st %)) X)
                        lo (reduce min zs) hi (reduce max zs)]
                    (doseq [x X] (line! (to-px x) (to-px (gen/decode st (gen/encode st x))) "rgba(255,143,112,0.35)"))
                    (set! (.-strokeStyle ctx) "#ffd479") (set! (.-lineWidth ctx) 1.5)
                    (.beginPath ctx)
                    (doseq [i (range 121) :let [[px py] (to-px (gen/decode st [(+ lo (* (- hi lo) (/ i 120)))]))]]
                      (if (zero? i) (.moveTo ctx px py) (.lineTo ctx px py)))
                    (.stroke ctx)
                    (doseq [x X] (dot! (to-px (gen/decode st (gen/encode st x))) 2 "#ff8f70")))
                  "variational autoencoder"
                  (doseq [x (repeatedly 400 #(gen/vae-sample st fixed))] (dot! (to-px x) 2 "#ff8f70"))
                  "GAN"
                  (doseq [x (gen/gan-sample st fixed 400)] (dot! (to-px x) 2 "#ff8f70")))
                (.restore ctx)
                (set! (.-textContent caption)
                      (str kind " on " (.-data controls) "  |  step " step
                           (case kind
                             "probabilistic PCA" (str "  |  sigma^2 " (.toFixed (:var st) 3))
                             ("linear autoencoder" "deep autoencoder")
                             (str "  |  mean squared reconstruction error "
                                  (.toFixed (gen/reconstruction-error (if (= kind "deep autoencoder")
                                                                        #(gen/decode st (gen/encode st %))
                                                                        #(gen/linear-decode st (gen/linear-encode st %)))
                                                                      X)
                                            4))
                             "")
                           "\nblue: the data; orange: "
                           (case kind
                             "probabilistic PCA" "samples from the fitted model; yellow: its 1 and 2 sd ellipses"
                             "linear autoencoder" "each point encoded as one number and decoded; yellow: the decoder's line"
                             "deep autoencoder" "each point encoded as one number and decoded; yellow: every number decoded"
                             "variational autoencoder" "new points, z drawn from N(0, I) and decoded"
                             "GAN" "the generator's points; behind, brighter where the discriminator calls a point real")))))
            (draw-series! [cw ch]
              (let [{:keys [st step]} @model
                    {:keys [net linear]} st
                    n 60
                    x0 30 y0 40 w (- cw 60) h (- ch 150)
                    to-px (fn [t v] [(+ x0 (* w (/ t (dec n)))) (+ y0 (* h (- 1 v)))])
                    trace! (fn [vals color]
                             (set! (.-strokeStyle ctx) color)
                             (.beginPath ctx)
                             (doseq [[t v] (map-indexed vector vals) :let [[px py] (to-px t v)]]
                               (if (zero? t) (.moveTo ctx px py) (.lineTo ctx px py)))
                             (.stroke ctx))
                    predict (fn [f] (concat (take 2 series) (map #(f (subvec series % (+ % 2))) (range (- n 2)))))
                    rms (fn [f] (math/sqrt (/ (reduce + (map #(let [d (- (f (subvec series % (+ % 2))) (series (+ % 2)))] (* d d)) (range 200))) 200)))]
                (set! (.-fillStyle ctx) "#151b26")
                (.fillRect ctx x0 y0 w h)
                (set! (.-lineWidth ctx) 1.5)
                (trace! (take n series) "#7dc4ff")
                (trace! (predict (partial gen/ar-predict linear)) "#7dffa0")
                (trace! (predict (partial gen/deep-ar-predict net)) "#ff8f70")
                (set! (.-textContent caption)
                      (str "autoregressive, order 2, on the logistic map x' = 3.9 x (1 - x)  |  epoch " step
                           "  |  rms error: linear " (.toFixed (rms (partial gen/ar-predict linear)) 3)
                           ", deep " (.toFixed (rms (partial gen/deep-ar-predict net)) 3)
                           "\nblue: the series; green: the linear AR model's prediction of each value from the two before; orange: the deep one's"))))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (if (= "autoregressive" (.-model controls)) (draw-series! cw ch) (draw-plane! cw ch))))
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
        (.onChange (.add gui controls "model" #js ["probabilistic PCA" "linear autoencoder" "deep autoencoder"
                                                   "variational autoencoder" "GAN" "autoregressive"])
                   restart!)
        (.onChange (.add gui controls "data" #js ["ring" "moon" "two blobs" "spiral"]) restart!)
        (.onChange (.add gui controls "points" 100 800 50) restart!)
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
                   (when-let [container (js/document.getElementById "generative")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
