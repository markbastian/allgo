(ns allgo.demo.classifiers
  "Canvas demo for `allgo.learning`: the learners of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapters 19
  to 21, on two-dimensional data, each drawing the regions it assigns to
  each class behind the points.

  Logistic regression and the perceptron can only draw a line, and fail
  on the circles and on XOR; k nearest neighbors follows every bump of
  the data; naive Bayes draws the boundary two Gaussians make; the SVM's
  Gaussian kernel curves round the classes, its support vectors ringed;
  a network with a hidden layer learns curves an epoch at a time, and a
  wider one with two can follow the noise -- unless trained with dropout
  (section 21.5.4), half its hidden units deactivated each step. The
  trees split on the coordinates binned into ten, so their regions are
  rectangles: one tree, a random forest's vote, or AdaBoost's weighted
  stumps."
  (:require [allgo.learning.models :as models]
            [allgo.learning.neural :as neural]
            [allgo.learning.trees :as trees]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:data "moons" :learner "network" :points 160 :noise 0.15 :seed 3 :speed 4})

(defn- dataset
  "`n` points in [-2, 2]^2 and labels 0 and 1."
  [kind n noise rng]
  (let [jitter #(* noise (random/gaussian rng 0.0 1.0))
        pts (for [i (range n)]
              (let [c (mod i 2)]
                (case kind
                  "blobs" [[(+ (if (zero? c) -0.8 0.8) (* 2 (jitter))) (+ (if (zero? c) -0.5 0.5) (* 2 (jitter)))] c]
                  "circles" (let [a (* 2 math/PI (rng)) r (if (zero? c) 0.6 1.5)]
                              [[(+ (* r (math/cos a)) (jitter)) (+ (* r (math/sin a)) (jitter))] c])
                  "XOR" (let [x (- (* 3.6 (rng)) 1.8) y (- (* 3.6 (rng)) 1.8)]
                          [[x y] (if (= (pos? x) (pos? y)) 0 1)])
                  "moons" (let [a (* math/PI (rng))]
                            (if (zero? c)
                              [[(+ -0.5 (math/cos a) (jitter)) (+ -0.25 (math/sin a) (jitter))] 0]
                              [[(+ 0.5 (- (math/cos a)) (jitter)) (+ 0.25 (- (math/sin a)) (jitter))] 1]))
                  "spiral" (let [t (* 3.5 (rng)) a (+ t (* c math/PI)) r (* 0.5 t)]
                             [[(+ (* r (math/cos a)) (jitter)) (+ (* r (math/sin a)) (jitter))] c]))))]
    [(mapv first pts) (mapv second pts)]))

(defn- binned [[x y]] {:bx (int (min 9 (max 0 (math/floor (* 2.5 (+ x 2)))))) :by (int (min 9 (max 0 (math/floor (* 2.5 (+ y 2))))))})

(def ^:private bin-values {:bx (vec (range 10)) :by (vec (range 10))})

(defn- trainer
  "The learner's states, a seq, and how to classify with one: `[states
  predict]`, predict returning P(class 1) or the class."
  [learner X ys rng]
  (let [examples (mapv (fn [x y] (assoc (binned x) :y y)) X ys)]
    (case learner
      "logistic regression" [(models/logistic-regression-seq X ys 1.0) (fn [m x] (models/predict-logistic m x))]
      "perceptron" [(models/perceptron-seq X ys 0.05) (fn [m x] (models/predict-threshold m x))]
      "k nearest neighbors" [(repeat (models/knn X ys 7)) (fn [h x] (h x))]
      "naive Bayes" [(repeat (models/gaussian-naive-bayes X ys)) (fn [h x] (h x))]
      "SVM, linear" [(models/svm-seq X (mapv #(if (= 1 %) 1 -1) ys) {:C 5.0 :rng rng})
                     (fn [m x] (if (= 1 (models/svm-predict m x)) 1 0))]
      "SVM, Gaussian kernel" [(models/svm-seq X (mapv #(if (= 1 %) 1 -1) ys) {:C 5.0 :kernel (models/gaussian-kernel 2.0) :rng rng})
                              (fn [m x] (if (= 1 (models/svm-predict m x)) 1 0))]
      "network" [(neural/train-seq [(neural/dense 2 12 :tanh rng) (neural/dense 12 1 :sigmoid rng)]
                                   X (mapv vector ys) {:alpha 0.05 :loss :cross-entropy :rng rng})
                 (fn [net x] (first (neural/predict net x)))]
      ("wide network" "wide network, dropout")
      [(neural/train-seq (cond-> [(neural/dense 2 32 :tanh rng) (neural/dense 32 32 :tanh rng) (neural/dense 32 1 :sigmoid rng)]
                           (= learner "wide network, dropout") (neural/with-dropout 0.5))
                         X (mapv vector ys) {:alpha 0.05 :loss :cross-entropy :rng rng})
       (fn [net x] (first (neural/predict net x)))]
      "decision tree" [(repeat (trees/decision-tree-learning examples [:bx :by] :y bin-values))
                       (fn [t x] (trees/classify t (binned x)))]
      "random forest" [(repeat (trees/random-forest examples [:bx :by] :y bin-values 25 rng :m 1))
                       (fn [f x] (f (binned x)))]
      "AdaBoost" [(trees/adaboost-seq examples [:bx :by] :y bin-values)
                  (fn [st x] (if (seq (:hypotheses st)) (trees/weighted-majority st (binned x)) 0))])))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        region (js/document.createElement "canvas")
        rctx (.getContext region "2d")
        running? (atom false)
        model (atom nil)
        res 64]
    (set! (.-width region) res) (set! (.-height region) res)
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls))
                    [X ys] (dataset (.-data controls) (.-points controls) (.-noise controls) rng)
                    ;; as many again from the same distribution, held out
                    [Xt yts] (dataset (.-data controls) (.-points controls) (.-noise controls)
                                      (random/rng (+ 1000 (.-seed controls))))
                    [states predict] (trainer (.-learner controls) X ys rng)]
                (reset! model {:X X :ys ys :Xt Xt :yts yts :rest states :st (first states) :predict predict :step 0})))
            (advance! []
              (let [n (.-speed controls)]
                (swap! model (fn [{:keys [rest step] :as m}]
                               (let [rest (nthrest rest n)]
                                 (assoc m :rest rest :st (first rest) :step (+ step n)))))))
            (draw! []
              (let [{:keys [X ys Xt yts st predict step]} @model
                    cw (.-width canvas) ch (.-height canvas)
                    side (min (- cw 40) (- ch 110))
                    x0 (* 0.5 (- cw side)) y0 40
                    to-px (fn [[x y]] [(+ x0 (* side (/ (+ x 2) 4))) (+ y0 (* side (/ (- 2 y) 4)))])
                    right (fn [X ys] (count (filter true? (map (fn [x y] (= y (if (> (predict st x) 0.5) 1 0))) X ys))))
                    correct (right X ys)]
                ;; The regions only when the model has changed: for the
                ;; learners that are done at once, never again.
                (when-not (identical? st (:drawn @model))
                  (let [img (.createImageData rctx res res)
                        data (.-data img)]
                    (dotimes [j res]
                      (dotimes [i res]
                        (let [x (- (* 4 (/ (+ i 0.5) res)) 2) y (- 2 (* 4 (/ (+ j 0.5) res)))
                              p (double (predict st [x y]))
                              o (* 4 (+ i (* res j)))]
                          (aset data o (int (+ 40 (* 180 p)))) (aset data (+ o 1) (int (+ 60 (* 40 (- 1 (abs (- p 0.5)))))))
                          (aset data (+ o 2) (int (+ 40 (* 180 (- 1 p))))) (aset data (+ o 3) 255))))
                    (.putImageData rctx img 0 0)
                    (swap! model assoc :drawn st)))
                (set! (.-fillStyle ctx) "#0d1118")
                (.fillRect ctx 0 0 cw ch)
                (set! (.-imageSmoothingEnabled ctx) true)
                (.drawImage ctx region x0 y0 side side)
                (doseq [[x y] (map vector X ys) :let [[px py] (to-px x)]]
                  (set! (.-fillStyle ctx) (if (= 1 y) "#ff8f70" "#7dc4ff"))
                  (set! (.-strokeStyle ctx) "#0d1118")
                  (.beginPath ctx) (.arc ctx px py 4 0 (* 2 math/PI)) (.fill ctx) (.stroke ctx))
                (when (:alphas st)
                  (doseq [sv (models/support-vectors st) :let [[px py] (to-px sv)]]
                    (set! (.-strokeStyle ctx) "#ffffff") (set! (.-lineWidth ctx) 1.5)
                    (.beginPath ctx) (.arc ctx px py 7 0 (* 2 math/PI)) (.stroke ctx)))
                (set! (.-textContent caption)
                      (str (.-learner controls) " on " (.-data controls) "  |  step " step
                           "  |  " correct " of " (count X) " training points right, " (right Xt yts) " of " (count Xt) " held out"
                           "\nred and blue: the classes, and the regions the learner gives each"
                           (when (:alphas st) "  |  ringed: support vectors")))))
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
        (.onChange (.add gui controls "data" #js ["moons" "circles" "XOR" "blobs" "spiral"]) restart!)
        (.onChange (.add gui controls "learner" #js ["logistic regression" "perceptron" "k nearest neighbors" "naive Bayes"
                                                     "SVM, linear" "SVM, Gaussian kernel" "network"
                                                     "wide network" "wide network, dropout"
                                                     "decision tree" "random forest" "AdaBoost"])
                   restart!)
        (.onChange (.add gui controls "points" 40 400 10) restart!)
        (.onChange (.add gui controls "noise" 0.0 0.5 0.01) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 50 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "classifiers")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
