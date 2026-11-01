(ns allgo.vision.detection
  "Classifying images and detecting objects: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25.

    classification   a convolutional network -- `allgo.learning.neural`'s
                     conv layer, then dense ones -- trained on image
                     patches, with the data augmented as the chapter
                     recommends: each example also shifted by a pixel and
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
