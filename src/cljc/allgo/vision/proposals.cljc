(ns allgo.vision.proposals
  "Region proposals and bounding-box regression: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25
  -- rather than score every window at every size, propose a few hundred
  boxes likely to hold something, classify those, and nudge each box onto
  the object it holds.

    selective search  Uijlings, van de Sande, Gevers and Smeulders'
                      grouping (\"Selective Search for Object
                      Recognition\", International Journal of Computer
                      Vision 104, 2013), the proposals of R-CNN: start
                      from a graph-based segmentation's regions, and
                      merge the two most similar neighbors, a merge a
                      step, until one is left; every region met along the
                      way, at every scale, is a proposal. Similarity is
                      the sum of color (the intersection of the regions'
                      histograms), size (small regions merge first, so
                      that growth is even across the image) and fill (how
                      well the two fill their joint box, so that pieces
                      of one thing go together); the paper's texture
                      measure is left out, since these drawn scenes have
                      none.
    anchors           the region proposal network of Ren, He, Girshick
                      and Sun's Faster R-CNN (NeurIPS 2015): boxes of a
                      few sizes and shapes at every position, each scored
                      for objectness by a learned classifier and the best,
                      after non-maximum suppression, kept. Here the
                      classifier is logistic regression on the anchor's
                      pixels, where the paper's shares the detector's
                      convolutional features.
    box regression    Girshick, Donahue, Darrell and Malik's (R-CNN,
                      CVPR 2014): from a proposal P to the object's box G,
                      the offsets of the center in units of P's size and
                      the logarithms of the size ratios, (Gx - Px) / Pw,
                      (Gy - Py) / Ph, log(Gw / Pw), log(Gh / Ph); each a
                      ridge regression on the proposal's features, learned
                      from proposals that overlap an object well.

  Boxes are `allgo.vision.detection`'s `[x0 y0 x1 y1]`; images gray. A
  proposal's features are its pixels -- the box with some context round
  it, resampled to a fixed grid, as R-CNN warps every proposal to its
  network's input."
  (:require [allgo.learning.models :as models]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            [allgo.vision.detection :as detection]
            [allgo.vision.image :as img]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Selective search

(def ^:private bins 25)

(defn- histogram
  "The normalized histogram of `values` in [0, 1], `bins` bins."
  [values]
  (let [n (double (count values))]
    (reduce (fn [h x] (update h (min (dec bins) (max 0 (long (* bins (double x))))) + (/ 1.0 n)))
            (vec (repeat bins 0.0)) values)))

(defn- union-box [[ax0 ay0 ax1 ay1] [bx0 by0 bx1 by1]]
  [(min ax0 bx0) (min ay0 by0) (max ax1 bx1) (max ay1 by1)])

(defn- similarity
  "Color, size and fill, each in [0, 1], summed."
  [image-size a b]
  (let [image-size (double image-size)
        sa (double (:size a)) sb (double (:size b))]
    (+ (reduce + (map min (:hist a) (:hist b)))
       (- 1.0 (/ (+ sa sb) image-size))
       (- 1.0 (/ (- (double (detection/area (union-box (:box a) (:box b)))) sa sb) image-size)))))

(defn- initial-regions
  "The regions of `labels`, each `{:size :box :hist :neighbors}`, by
  label, and the image's size."
  [image {:keys [w h px]}]
  (let [w (long w) h (long h)
        pixels (group-by (fn [i] (px i)) (range (* w h)))
        neighbors (reduce (fn [m [a b]] (-> m (update a (fnil conj #{}) b) (update b (fnil conj #{}) a)))
                          {}
                          (for [y (range h) x (range w)
                                [dx dy] [[1 0] [0 1]]
                                :let [x2 (+ x dx) y2 (+ y dy)]
                                :when (and (< x2 w) (< y2 h))
                                :let [a (px (+ (* y w) x)) b (px (+ (* y2 w) x2))]
                                :when (not= a b)]
                            [a b]))]
    (into {} (for [[label is] pixels]
               [label {:size (count is)
                       :box [(apply min (map #(mod % w) is)) (apply min (map #(quot % w) is))
                             (inc (long (apply max (map #(mod % w) is)))) (inc (long (apply max (map #(quot % w) is))))]
                       :hist (histogram (map (:px image) is))
                       :neighbors (get neighbors label #{})}]))))

(defn- merge-step
  "The two most similar neighboring regions merged into a new one."
  [{:keys [regions size next-label] :as st}]
  (let [pairs (for [[a ra] regions b (:neighbors ra) :when (< (compare a b) 0)] [a b])]
    (when (seq pairs)
      (let [[a b] (apply max-key (fn [[a b]] (similarity size (regions a) (regions b))) pairs)
            ra (regions a) rb (regions b)
            sa (double (:size ra)) sb (double (:size rb))
            merged {:size (+ (long (:size ra)) (long (:size rb)))
                    :box (union-box (:box ra) (:box rb))
                    :hist (mapv (fn [x y] (/ (+ (* sa x) (* sb y)) (+ sa sb))) (:hist ra) (:hist rb))
                    :neighbors (disj (into (:neighbors ra) (:neighbors rb)) a b)}
            regions (reduce (fn [rs n] (update-in rs [n :neighbors] #(-> % (disj a b) (conj next-label))))
                            (-> regions (dissoc a b) (assoc next-label merged))
                            (:neighbors merged))]
        (-> st
            (assoc :regions regions :merged [a b next-label])
            (update :next-label inc)
            (update :proposals conj (:box merged)))))))

(defn selective-search-seq
  "Selective search over `image`, starting from the regions of `labels`
  (a segmentation, such as `allgo.vision.segmentation`'s), a merge a
  step. Each state `{:regions :proposals :merged}`: the regions left, by
  label, the boxes so far -- every starting region's, then each merged
  one's -- and the merge just made, `[a b new]`."
  [image labels]
  (let [regions (initial-regions image labels)
        {:keys [w h]} image]
    (take-while some?
                (iterate merge-step
                         {:regions regions
                          :size (* (long w) (long h))
                          :next-label (inc (long (apply max (keys regions))))
                          :proposals (vec (distinct (map :box (vals regions))))}))))

(defn selective-search
  "Every box selective search proposes for `image` from the regions of
  `labels`, without repeats."
  [image labels]
  (vec (distinct (:proposals (last (selective-search-seq image labels))))))

;; ---------------------------------------------------------------------------
;; Features

(defn features
  "What a regressor or classifier sees of `box` in `image`: the box grown
  by `context` (1.5 makes it half as big again, about its center),
  sampled bilinearly on a `k` by `k` grid."
  [image [x0 y0 x1 y1] k context]
  (let [k (long k)
        cx (* 0.5 (+ (double x0) (double x1))) cy (* 0.5 (+ (double y0) (double y1)))
        hw (* 0.5 (double context) (- (double x1) (double x0))) hh (* 0.5 (double context) (- (double y1) (double y0)))]
    (vec (for [j (range k) i (range k)]
           (img/sample image
                       (- (+ (- cx hw) (* (/ (+ i 0.5) k) 2.0 hw)) 0.5)
                       (- (+ (- cy hh) (* (/ (+ j 0.5) k) 2.0 hh)) 0.5))))))

;; ---------------------------------------------------------------------------
;; Box regression

(defn- center-size [[x0 y0 x1 y1]]
  [(* 0.5 (+ (double x0) (double x1))) (* 0.5 (+ (double y0) (double y1)))
   (- (double x1) (double x0)) (- (double y1) (double y0))])

(defn box-deltas
  "R-CNN's regression targets from proposal `p` to box `g`: `[tx ty tw
  th]`."
  [p g]
  (let [[px py pw ph] (center-size p) [gx gy gw gh] (center-size g)]
    [(/ (- gx px) pw) (/ (- gy py) ph) (math/log (/ gw pw)) (math/log (/ gh ph))]))

(defn apply-deltas
  "Proposal `p` moved by `[tx ty tw th]`: the inverse of `box-deltas`."
  [p [tx ty tw th]]
  (let [[px py pw ph] (center-size p)
        cx (+ px (* pw (double tx))) cy (+ py (* ph (double ty)))
        w (* pw (math/exp tw)) h (* ph (math/exp th))]
    [(- cx (* 0.5 w)) (- cy (* 0.5 h)) (+ cx (* 0.5 w)) (+ cy (* 0.5 h))]))

(defn- accumulate
  "The normal equations `{:xtx :xty}` with the example `x` (a bias in
  front) and its targets `ys` -- one for each regression -- added."
  [{:keys [xtx xty]} x ys]
  (let [xb (into [1.0] x)]
    {:xtx (mapv (fn [row a] (let [a (double a)] (mapv (fn [v b] (+ (double v) (* a (double b)))) row xb))) xtx xb)
     :xty (mapv (fn [col y] (let [y (double y)] (mapv (fn [v a] (+ (double v) (* y (double a)))) col xb))) xty ys)}))

(defn- empty-equations [n targets]
  {:xtx (vec (repeat (inc (long n)) (vec (repeat (inc (long n)) 0.0))))
   :xty (vec (repeat targets (vec (repeat (inc (long n)) 0.0))))})

(defn- solve-ridge
  "The weights of each regression from its normal equations, `lambda`
  added to the diagonal but for the bias."
  [{:keys [xtx xty]} lambda]
  (let [A (vec (map-indexed (fn [i row] (if (pos? (long i)) (update row i + lambda) row)) xtx))]
    (mapv #(lin/cholesky-solve A %) xty)))

(defn ridge
  "Ridge regression of `ys` on the rows of `X`: the weights, a bias first,
  minimizing the squared error plus `lambda` times the squared weights
  (the bias free) -- by the normal equations, (X^T X + lambda I) w = X^T
  y, through `allgo.numerics.linear`'s Cholesky."
  [X ys lambda]
  (first (solve-ridge (reduce (fn [eqs [x y]] (accumulate eqs x [y])) (empty-equations (count (first X)) 1) (map vector X ys))
                      lambda)))

(defn box-regression-seq
  "Learning a box regressor from `examples`, each `[image proposal
  truth]`, an example a step: the normal equations of the four ridge
  regressions -- one for each delta, on the proposal's `features` -- with
  each example added. `{:xtx :xty :n :k :context}`."
  [examples k context]
  (reductions (fn [eqs [image p g]]
                (-> eqs (merge (accumulate eqs (features image p k context) (box-deltas p g))) (update :n inc)))
              (assoc (empty-equations (* (long k) (long k)) 4) :n 0 :k k :context context)
              examples))

(defn regressor
  "The box regressor from the normal equations of `box-regression-seq`,
  with ridge penalty `lambda`: `{:weights :k :context}`."
  [{:keys [k context] :as eqs} lambda]
  {:weights (solve-ridge eqs lambda) :k k :context context})

(defn box-regressor
  "A box regressor learned from `examples`, each `[image proposal
  truth]`: four ridge regressions, one for each delta, on the proposal's
  `features`. `{:weights :k :context}`."
  [examples k context lambda]
  (regressor (last (box-regression-seq examples k context)) lambda))

(defn refine
  "Proposal `p` in `image` moved where `regressor` says its object is."
  [{:keys [weights k context]} image p]
  (let [x (features image p k context)]
    (apply-deltas p (mapv #(models/predict-linear % x) weights))))

;; ---------------------------------------------------------------------------
;; Anchors

(defn anchors
  "Boxes of each of `sizes` and `ratios` (height over width, the area
  kept) centered every `stride` pixels, those wholly inside the image."
  [{:keys [w h]} stride sizes ratios]
  (for [cy (range (* 0.5 (double stride)) h stride)
        cx (range (* 0.5 (double stride)) w stride)
        s sizes r ratios
        :let [bw (/ (double s) (math/sqrt r)) bh (* (double s) (math/sqrt r))
              box [(- cx (* 0.5 bw)) (- cy (* 0.5 bh)) (+ cx (* 0.5 bw)) (+ cy (* 0.5 bh))]]
        :when (and (>= (box 0) 0) (>= (box 1) 0) (<= (box 2) w) (<= (box 3) h))]
    box))

(defn objectness-examples
  "Training examples for objectness from `image` with objects at `truth`,
  sampled as Faster R-CNN samples them: anchors overlapping an object by
  at least `hi` are positive, those overlapping every object by under
  `lo` negative, the in-between left out; `batch` of them drawn at
  random, at most half positive, the rest negative. Only those drawn have
  their features computed. `[X ys]`."
  [image truth anchors k context hi lo batch rng]
  (let [best (fn [a] (reduce max 0.0 (map #(detection/iou a %) truth)))
        {pos true neg false} (group-by #(>= (double (best %)) (double hi))
                                       (filter #(let [b (double (best %))] (or (>= b (double hi)) (< b (double lo)))) anchors))
        pos (take (quot (long batch) 2) (random/shuffle rng pos))
        neg (take (- (long batch) (count pos)) (random/shuffle rng neg))]
    [(mapv #(features image % k context) (concat pos neg))
     (into (vec (repeat (count pos) 1.0)) (repeat (count neg) 0.0))]))

(defn objectness-seq
  "Logistic regression for objectness on `[X ys]`, an epoch a step: the
  weights, a bias first, after each."
  [[X ys] alpha]
  (models/logistic-regression-seq X ys alpha))

(defn propose
  "The anchors of `image` that `weights` scores over `min-score`, after
  non-maximum suppression at `threshold`: `[{:box :score} ...]`, best
  first."
  [weights image anchors k context min-score threshold]
  (detection/non-maximum-suppression
   (for [a anchors
         :let [s (models/predict-logistic weights (features image a k context))]
         :when (> s min-score)]
     {:box a :score s})
   threshold))
