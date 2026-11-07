(ns allgo.vision.detection
  "Classifying images and detecting objects: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25.

    classification   a convolutional network -- `allgo.learning.neural`'s
                     conv layer, then dense ones -- trained on image
                     patches, with the data augmented: each example also
                     shifted by a pixel, as the chapter recommends (it
                     names shifts, rotations, stretches and hue), and
                     mirrored, which teaches the network that those do
                     not change what a thing is; and, for a detector,
                     the examples taken from whole scenes: the objects'
                     windows, and random windows that barely touch them,
                     then the windows the detector wrongly fires on mined
                     and added as negatives
    detection        a classifier slid over the image, a window every
                     `stride` pixels, each window's score whether it holds
                     the object; then non-maximum suppression, since an
                     object scores well in many overlapping windows: the
                     best window kept, those overlapping it by more than a
                     threshold of intersection over union dropped, and
                     again, a window a step
    test time        a window's shifted and mirrored copies scored too,
                     the scores averaged
    ROI pooling      a box of any size reduced to a fixed grid of
                     maxima, the same number of features from every box
    evaluation       precision and recall, the detections taken best
                     first and each matched to an object not yet found;
                     average precision, the area under that curve

  A box is `[x0 y0 x1 y1]`, the corners exclusive at the far side."
  (:require [allgo.learning.neural :as neural]
            [allgo.vision.image :as img]))

(defn area [[x0 y0 x1 y1]] (* (max 0 (- x1 x0)) (max 0 (- y1 y0))))

(defn iou
  "Intersection over union: the area two boxes share over the area they
  cover between them."
  [[ax0 ay0 ax1 ay1 :as a] [bx0 by0 bx1 by1 :as b]]
  (let [inter (area [(max ax0 bx0) (max ay0 by0) (min ax1 bx1) (min ay1 by1)])
        union (- (+ (area a) (area b)) inter)]
    (if (pos? union) (/ (double inter) union) 0.0)))

(defn nms-seq
  "Non-maximum suppression of `detections`, each `{:box :score}`, a window
  a step: each state `{:kept :remaining}`, the best remaining kept and
  those overlapping it by more than `threshold` dropped."
  [detections threshold]
  (take-while some?
              (iterate (fn [{:keys [kept remaining]}]
                         (when-let [best (first remaining)]
                           {:kept (conj kept best)
                            :remaining (remove #(> (iou (:box best) (:box %)) threshold) (rest remaining))}))
                       {:kept [] :remaining (sort-by (comp - :score) detections)})))

(defn non-maximum-suppression [detections threshold]
  (:kept (last (nms-seq detections threshold))))

(defn windows
  "Every `size` by `size` box fully inside the image, every `stride`
  pixels."
  [{:keys [w h]} size stride]
  (for [y (range 0 (inc (- h size)) stride) x (range 0 (inc (- w size)) stride)]
    [x y (+ x size) (+ y size)]))

(defn patch
  "The pixels of `box` of a gray image, a flat vector row by row -- a
  network's input."
  [image [x0 y0 x1 y1]]
  (:px (img/crop image x0 y0 (- x1 x0) (- y1 y0))))

(defn training-windows
  "Examples for a detector from an `image` whose objects are known, by
  their `size` by `size` boxes: `{:positives :negatives}`, patches. The
  positives are the objects' own windows; the negatives `n` random
  windows overlapping every object by less than `max-iou` -- partial
  objects and clutter included, which is what a sliding window will meet,
  where patches of whole things alone would not teach it."
  [image boxes size n max-iou rng]
  (let [{:keys [w h]} image]
    {:positives (mapv #(patch image %) boxes)
     :negatives (->> (repeatedly (fn [] (let [x (long (* (rng) (inc (- w size)))) y (long (* (rng) (inc (- h size))))]
                                          [x y (+ x size) (+ y size)])))
                     (filter (fn [b] (every? #(< (iou b %) max-iou) boxes)))
                     (take n)
                     (mapv #(patch image %)))}))

(defn augment
  "A patch of `size` by `size` and its variants: shifted a pixel each way
  (the edge pixels repeated) and mirrored left to right."
  [p size]
  (let [im {:w size :h size :px p}
        shifted (fn [dx dy] (:px (img/image size size (fn [x y] (img/at im (- x dx) (- y dy))))))
        mirrored (:px (img/image size size (fn [x y] (img/at im (- size 1 x) y))))]
    [p (shifted 1 0) (shifted -1 0) (shifted 0 1) (shifted 0 -1) mirrored]))

(defn classifier
  "A small convolutional network for `size` by `size` patches: `kernels`
  3 by 3 kernels, a dense hidden layer of `hidden`, and a sigmoid
  output, the probability the patch holds the object."
  [size kernels hidden rng]
  (let [conv (neural/conv [size size] 3 kernels :relu rng)]
    [conv (neural/dense (neural/layer-out-size conv) hidden :relu rng) (neural/dense hidden 1 :sigmoid rng)]))

(defn examples
  "The training set from `positives` and `negatives`, each augmented:
  `[X Y]`, the targets 1 and 0."
  [positives negatives size]
  (let [pos (mapcat #(augment % size) positives)
        neg (mapcat #(augment % size) negatives)]
    [(vec (concat pos neg)) (vec (concat (repeat (count pos) [1.0]) (repeat (count neg) [0.0])))]))

(defn train-seq
  "The classifier trained on `positives` and `negatives`, patches as flat
  vectors, each augmented: the networks after each epoch."
  [network positives negatives size rng & {:keys [alpha] :or {alpha 0.05}}]
  (let [[X Y] (examples positives negatives size)]
    (neural/train-seq network X Y {:alpha alpha :loss :cross-entropy :rng rng})))

(defn score [network p] (first (neural/predict network p)))

(defn detect
  "The objects `network` finds in a gray image: every window of `size`
  every `stride` pixels scored, those over `min-score` kept, and
  non-maximum suppression at `threshold`. `[{:box :score} ...]`."
  [network image size stride min-score threshold]
  (non-maximum-suppression
   (for [b (windows image size stride)
         :let [s (score network (patch image b))]
         :when (> s min-score)]
     {:box b :score s})
   threshold))

(defn hard-negatives
  "The windows of a training `image` that `network` wrongly takes for the
  object -- scoring over `min-score` yet overlapping every true box by
  less than `max-iou` -- as patches: the negatives to train on next.
  Random negatives seldom land squarely on the clutter most like the
  object; these are exactly that, so a round or two of mining them makes
  the detector far pickier."
  [network image boxes size stride min-score max-iou]
  (vec (for [b (windows image size stride)
             :when (every? #(< (iou b %) max-iou) boxes)
             :let [p (patch image b)]
             :when (> (score network p) min-score)]
         p)))

(defn score-augmented
  "The network's score for a patch with augmentation at test time, as
  section 25.4.2 offers: the patch's variants -- `augment`'s shifts and
  mirror -- each scored and the scores averaged, a vote."
  [network p size]
  (let [ss (map #(score network %) (augment p size))]
    (/ (reduce + ss) (count ss))))

;; ---------------------------------------------------------------------------
;; Pooling a region to a fixed size

(defn roi-pool
  "A region of interest of any size and shape reduced to `rows` by
  `cols` numbers: `box` of the image divided into a grid of that many
  cells, as nearly equal as whole pixels allow, and each cell's largest
  value taken -- so a classifier made for one size of input can read
  boxes of every size. Row by row, a flat vector."
  [image [x0 y0 x1 y1] rows cols]
  (let [edges (fn [a b n] (mapv #(+ a (quot (* % (- b a)) n)) (range (inc n))))
        xs (edges x0 x1 cols) ys (edges y0 y1 rows)]
    (vec (for [i (range rows) j (range cols)]
           (let [[ya yb] [(ys i) (max (ys (inc i)) (inc (ys i)))]
                 [xa xb] [(xs j) (max (xs (inc j)) (inc (xs j)))]]
             (reduce max (for [y (range ya yb) x (range xa xb)] (img/at image x y))))))))

;; ---------------------------------------------------------------------------
;; Evaluation

(defn precision-recall
  "How a detector did against the truth, over images: `results` a seq of
  `[detections boxes]` per image, each detection `{:box :score}`. Taken
  best score first, a detection is a hit if it overlaps by at least
  `min-iou` a true box no better detection has claimed, else a false
  alarm. The curve, a point per detection: `[{:score :precision
  :recall}]`, precision the fraction of detections so far that hit,
  recall the fraction of true boxes found."
  [results min-iou]
  (let [total (reduce + (map (comp count second) results))
        ranked (sort-by (comp - :score first)
                        (for [[i [dets _]] (map-indexed vector results) d dets] [d i]))
        truth (mapv second results)]
    (:curve (reduce (fn [{:keys [claimed hits n] :as acc} [{:keys [box score]} i]]
                      (let [free (remove #(claimed [i %]) (range (count (truth i))))
                            best (when (seq free) (apply max-key #(iou box ((truth i) %)) free))
                            hit? (and best (>= (iou box ((truth i) best)) min-iou))
                            hits (if hit? (inc hits) hits) n (inc n)]
                        (-> acc
                            (assoc :hits hits :n n)
                            (cond-> hit? (update :claimed conj [i best]))
                            (update :curve conj {:score score :precision (/ (double hits) n)
                                                 :recall (if (pos? total) (/ (double hits) total) 0.0)}))))
                    {:claimed #{} :hits 0 :n 0 :curve []}
                    ranked))))

(defn average-precision
  "The area under the precision-recall `curve`, precision at each recall
  taken as the best at that recall or beyond -- one number balancing
  finding everything against finding only what is there."
  [curve]
  (let [pts (vec curve)
        best-after (vec (reductions max (rseq (mapv :precision pts))))
        envelope (fn [i] (best-after (- (count pts) 1 i)))]
    (reduce + (map (fn [i] (* (envelope i) (- (:recall (pts i)) (if (pos? i) (:recall (pts (dec i))) 0.0))))
                   (range (count pts))))))
