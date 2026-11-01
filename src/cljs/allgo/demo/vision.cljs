(ns allgo.demo.vision
  "Canvas demo for `allgo.vision`: the computer vision of Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 25, each method a step at a time on images drawn so that the
  right answer is known.

  Edges: Canny's stages -- smoothing, gradient, thinning, and the edges
  growing out from the strong ones. Segmentation: the regions merging an
  edge at a time, cheapest first; or k-means on the colors. Optical flow:
  a textured disk moving over a drifting backdrop, and the flow Lucas and
  Kanade find, as arrows. Stereo: the left and right views of boxes at
  different depths and the disparity found a row at a time, brighter
  nearer. Detection: a convolutional network learning disks from
  squares, a few examples a frame; at the end of each epoch it is slid
  over a new scene, and after a few, the windows it wrongly fires on are
  mined and added to its negatives."
  (:require [allgo.graph :as graph]
            [allgo.learning.neural :as neural]
            [allgo.learning.unsupervised :as unsupervised]
            [allgo.random :as random]
            [allgo.vision.detection :as det]
            [allgo.vision.edges :as edges]
            [allgo.vision.flow :as flow]
            [allgo.vision.image :as img]
            [allgo.vision.segmentation :as seg]
            [allgo.vision.stereo :as stereo]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "edges" :seed 1 :speed 3})

;; ---------------------------------------------------------------------------
;; Scenes

(defn- scene
  "A color scene of a few shapes over a faint pattern, with a little noise."
  [rng]
  (let [palette [[0.85 0.35 0.25] [0.3 0.75 0.4] [0.95 0.8 0.3] [0.35 0.5 0.9] [0.8 0.4 0.8]]
        back (img/pattern rng 4 20.0)
        shapes (for [i (range 5)]
                 (let [c (palette i) x (+ 8 (* 80 (rng))) y (+ 8 (* 48 (rng)))]
                   (if (< (rng) 0.5)
                     {:kind :disk :center [x y] :radius (+ 6 (* 8 (rng))) :value c}
                     (let [a (+ 6 (* 10 (rng))) b (+ 6 (* 8 (rng)))]
                       {:kind :rect :box [(- x a) (- y b) (+ x a) (+ y b)] :value c}))))]
    (img/add-noise (img/render 96 64 (fn [x y] (let [v (+ 0.1 (* 0.12 (back x y)))] [v v (+ v 0.05)])) shapes) 0.02 rng)))

(def ^:private det-size 12)

(defn- det-scene [rng]
  (let [slots (vec (random/shuffle rng (for [x [8 20 32 44] y [8 22]] [x y])))
        disks (mapv (fn [[x y]] [(+ x (- (* 2 (rng)) 1.0)) (+ y (- (* 2 (rng)) 1.0))]) (take 3 slots))
        squares (take 3 (drop 3 slots))
        shapes (concat (map (fn [c] {:kind :disk :center c :radius (+ 3.5 (rng)) :value 0.9}) disks)
                       (map (fn [[x y]] (let [r (+ 3 (rng))] {:kind :rect :box [(- x r) (- y r) (+ x r) (+ y r)] :value 0.9})) squares))]
    {:image (img/add-noise (img/render 56 32 0.1 shapes) 0.05 rng)
     :boxes (mapv (fn [[x y]] (let [x0 (math/round (- x 6)) y0 (math/round (- y 6))] [x0 y0 (+ x0 det-size) (+ y0 det-size)])) disks)}))

;; ---------------------------------------------------------------------------
;; Each method as a seq of frames to show

(defn- region-colors
  "A color for each region root, by hashing it."
  [root]
  (let [h (hash root)]
    [(+ 0.25 (* 0.7 (/ (mod h 97) 97.0))) (+ 0.25 (* 0.7 (/ (mod (quot h 97) 89) 89.0))) (+ 0.25 (* 0.7 (/ (mod (quot h 8633) 83) 83.0)))]))

(defn- frames [method rng]
  (case method
    "edges"
    (let [color (scene rng)
          gray (img/gray color)]
      (map (fn [st] {:left color :stage st}) (edges/canny-seq gray :sigma 1.2 :low 0.03 :high 0.08)))

    "segmentation"
    (let [color (scene rng)
          smoothed (img/smooth color 0.8)
          states (vec (seg/graph-based-seq smoothed 1.5))
          ;; every 200th edge, then the specks tidied into their neighbors
          shown (concat (take-nth 200 states) [(peek states) (seg/merge-small (peek states) smoothed 20)])]
      (map (fn [st] {:left color :forest st :smoothed smoothed}) shown))

    "k-means"
    (let [color (scene rng)
          points (vec (:px color))]
      (map (fn [st] {:left color :kmeans st}) (take 15 (rest (unsupervised/k-means-seq points 6 rng)))))

    "optical flow"
    (let [back (img/pattern rng 6 10.0) disk (img/pattern rng 6 5.0)
          frame (fn [t] (img/image 72 48 (fn [x y]
                                           (let [dx (- x (+ 20 (* 1.2 t))) dy (- y (+ 18 (* 0.6 t)))]
                                             (if (< (+ (* dx dx) (* dy dy)) 100)
                                               (disk (- x (* 1.2 t)) (- y (* 0.6 t)))
                                               (* 0.6 (back (+ x (* 0.5 t)) y)))))))]
      (map (fn [t] (let [a (frame t) b (frame (inc t))]
                     {:left a :flow (flow/flow-field a b 6 6 4 3 :lucas-kanade)}))
           (range 0 28)))

    "stereo"
    (let [layers [{:texture (img/pattern rng 8 6.0) :disparity 1}
                  {:texture (img/pattern rng 8 4.0) :disparity 5 :box [8 8 30 40]}
                  {:texture (img/pattern rng 8 4.0) :disparity 9 :box [40 14 62 34]}]
          [l r] (stereo/stereo-pair layers 72 48)]
      (map (fn [d] {:left l :right r :disparity d}) (stereo/disparity-seq l r 3 12)))

    "detection" nil))

;; ---------------------------------------------------------------------------
;; Detection, trained a few examples a frame

(defn- detection-init [rng]
  (let [train (vec (repeatedly 6 #(det-scene rng)))
        ex (map #(det/training-windows (:image %) (:boxes %) det-size 10 0.2 rng) train)
        pos (vec (mapcat :positives ex)) neg (vec (mapcat :negatives ex))
        [X Y] (det/examples pos neg det-size)]
    {:train train :pos pos :neg neg :X X :Y Y :order (random/shuffle rng (range (count X))) :i 0
     :net (det/classifier det-size 4 8 rng) :epoch 0 :round 0 :test (det-scene rng) :dets nil :rng rng}))

(defn- detection-advance [{:keys [X Y order i net epoch round train pos neg test rng] :as st} n]
  (let [k (min n (- (count order) i))
        net (reduce (fn [net j] (neural/sgd-step net (X j) (Y j) 0.05 :cross-entropy)) net (subvec order i (+ i k)))
        i (+ i k)]
    (if (< i (count order))
      (assoc st :net net :i i)
      ;; an epoch done: detect on the test scene; every few, mine
      (let [epoch (inc epoch)
            dets (det/detect net (:image test) det-size 2 0.5 0.3)
            mine? (and (zero? (mod epoch 6)) (< round 2))]
        (if mine?
          (let [hard (vec (mapcat #(det/hard-negatives net (:image %) (:boxes %) det-size 2 0.3 0.2) train))
                neg (into neg hard)
                [X Y] (det/examples pos neg det-size)]
            (assoc st :net net :epoch epoch :dets dets :neg neg :X X :Y Y :i 0 :round (inc round) :mined (count hard)
                   :order (random/shuffle rng (range (count X)))))
          (assoc st :net net :epoch epoch :dets dets :i 0 :order (random/shuffle rng order)))))))

;; ---------------------------------------------------------------------------
;; Drawing

(defn- to-rgb [v]
  (cond (vector? v) v
        (nil? v) [0.0 0.0 0.0]
        :else [v v v]))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        scratch (js/document.createElement "canvas")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls)) method (.-method controls)]
                (reset! model (if (= method "detection")
                                {:det (detection-init rng) :clock 0.0}
                                (let [fs (frames method rng)] {:rest fs :frame (first fs) :step 0 :clock 0.0})))))
            (blit! [{:keys [w h px]} f x y side-w side-h]
              ;; an image drawn to fit side-w by side-h at (x, y), pixels square
              (set! (.-width scratch) w) (set! (.-height scratch) h)
              (let [c2 (.getContext scratch "2d") data (.createImageData c2 w h) d (.-data data)]
                (dotimes [i (* w h)]
                  (let [[r g b] (to-rgb (f (px i))) o (* 4 i)]
                    (aset d o (* 255 (max 0 (min 1 r)))) (aset d (+ o 1) (* 255 (max 0 (min 1 g))))
                    (aset d (+ o 2) (* 255 (max 0 (min 1 b)))) (aset d (+ o 3) 255)))
                (.putImageData c2 data 0 0))
              (let [s (min (/ side-w w) (/ side-h h))]
                (set! (.-imageSmoothingEnabled ctx) false)
                (.drawImage ctx scratch x y (* s w) (* s h))
                s))
            (panels [n]
              ;; n side-by-side panels below the controls: [[x y w h] ...]
              (let [cw (.-width canvas) ch (.-height canvas)
                    gap 16 pw (/ (- cw (* gap (inc n))) n) ph (- ch 200)]
                (mapv (fn [i] [(+ gap (* i (+ pw gap))) 120 pw ph]) (range n))))
            (label! [text x y]
              (set! (.-fillStyle ctx) "#9fb4d8") (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
              (.fillText ctx text x (- y 8)))
            (advance! []
              (let [method (.-method controls) n (.-speed controls)]
                (if (= method "detection")
                  (swap! model update :det detection-advance (* 8 n))
                  (let [{:keys [clock]} @model
                        rate (case method "segmentation" 7.0 "stereo" 12.0 "optical flow" 4.0 2.0)
                        clock (+ clock (* rate n 0.016))]
                    (if (>= clock 1.0)
                      (swap! model (fn [{:keys [rest step] :as m}]
                                     (if (next rest)
                                       (assoc m :rest (next rest) :frame (fnext rest) :step (inc step) :clock 0.0)
                                       ;; flow loops; the rest hold their last frame
                                       (if (= method "optical flow")
                                         (let [fs (frames method (random/rng (.-seed controls)))] (assoc m :rest fs :frame (first fs) :step 0 :clock 0.0))
                                         (assoc m :clock 0.0)))))
                      (swap! model assoc :clock clock))))))
            (draw! []
              (let [method (.-method controls)
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "edges" (draw-edges!)
                             "segmentation" (draw-segmentation!)
                             "k-means" (draw-kmeans!)
                             "optical flow" (draw-flow!)
                             "stereo" (draw-stereo!)
                             "detection" (draw-detection!))))))
            (draw-edges! []
              (let [{{:keys [left stage]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (case (:stage stage)
                  :smoothed (do (blit! (:image stage) identity x2 y2 w2 h2) (label! "smoothed by a Gaussian" x2 y2))
                  :gradient (let [m (apply max 1e-9 (:px (:image stage)))]
                              (blit! (:image stage) #(/ % m) x2 y2 w2 h2) (label! "gradient magnitude" x2 y2))
                  :thinned (let [m (apply max 1e-9 (:px (:image stage)))]
                             (blit! (:image stage) #(/ % m) x2 y2 w2 h2) (label! "thinned: maxima across the edge" x2 y2))
                  :hysteresis (let [e (:edges stage) {:keys [w h]} left]
                                (blit! (img/image w h (fn [x y] (if (e [x y]) 1.0 0.08))) identity x2 y2 w2 h2)
                                (label! "edges, grown from the strong ones" x2 y2)))
                (str (name (:stage stage)) (when-let [e (:edges stage)] (str ", " (count e) " edge pixels"))
                     "\nleft: the scene  right: Canny's current stage -- of brightness only, so colors as bright as each other meet without an edge")))
            (draw-segmentation! []
              (let [{{:keys [left forest smoothed]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    {:keys [w h]} smoothed
                    ;; colored by each region's root, which stays put as
                    ;; regions merge, so the colors do not flicker
                    roots (img/image w h (fn [x y] (graph/canonical (:forest forest) [x y])))
                    colors (into {} (map (fn [r] [r (region-colors r)]) (distinct (:px roots))))]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (blit! roots colors x2 y2 w2 h2) (label! "regions, merged cheapest edge first" x2 y2)
                (str (count colors) " regions"
                     "\nright: each region a color; pixels merge across the weakest boundaries first")))
            (draw-kmeans! []
              (let [{{:keys [left kmeans]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    {:keys [centers assignment cost]} kmeans]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (blit! (assoc left :px assignment) centers x2 y2 w2 h2) (label! "each pixel painted its cluster's center" x2 y2)
                (str "cost " (.toFixed cost 2)
                     "\nright: 6 clusters of color, no regard for where pixels are")))
            (draw-flow! []
              (let [{{:keys [left flow]} :frame} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! left identity x1 y1 w1 h1)]
                (set! (.-strokeStyle ctx) "#ffe08a") (set! (.-fillStyle ctx) "#ffe08a") (set! (.-lineWidth ctx) 2.5)
                (doseq [[[x y] [u v]] flow
                        :let [px (+ x1 (* s (+ x 0.5))) py (+ y1 (* s (+ y 0.5)))]]
                  (.beginPath ctx) (.moveTo ctx px py) (.lineTo ctx (+ px (* 4 s u)) (+ py (* 4 s v))) (.stroke ctx)
                  (.fillRect ctx (- px 1) (- py 1) 2 2))
                (str (count flow) " windows"
                     "\narrows: the flow by Lucas-Kanade, four times its length -- the disk moving, the backdrop drifting")))
            (draw-stereo! []
              (let [{{:keys [left right disparity]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2] [x3 y3 w3 h3]] (panels 3)]
                (blit! left identity x1 y1 w1 h1) (label! "left eye" x1 y1)
                (blit! right identity x2 y2 w2 h2) (label! "right eye" x2 y2)
                (blit! disparity (fn [d] (let [t (/ d 12.0)] [(* 1.2 t) (* 0.8 t) (- 1.0 t)])) x3 y3 w3 h3)
                (label! "disparity: nearer is warmer" x3 y3)
                "a row at a time: each pixel's best match along its row in the other eye"))
            (draw-detection! []
              (let [{{:keys [test dets epoch round i X mined]} :det} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! (:image test) identity x1 y1 w1 h1)
                    rect! (fn [[a b c d] color width]
                            (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                            (.strokeRect ctx (+ x1 (* s a)) (+ y1 (* s b)) (* s (- c a)) (* s (- d b))))]
                (doseq [b (:boxes test)] (rect! b "rgba(255,255,255,0.35)" 1))
                (doseq [{:keys [box]} dets]
                  (rect! box (if (some #(> (det/iou box %) 0.5) (:boxes test)) "#7dffa0" "#ff8f70") 2.5))
                (str "epoch " epoch ", example " i " of " (count X)
                     (when (pos? round) (str "  |  " round " round" (when (> round 1) "s") " of hard negatives mined"
                                             (when mined (str " (" mined " last time)"))))
                     "\nfaint: the true disks  green: detections that hit one  red: false alarms among the squares")))
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
        (.onChange (.add gui controls "method" #js ["edges" "segmentation" "k-means" "optical flow" "stereo" "detection"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 10 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "vision")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
