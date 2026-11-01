(ns allgo.vision.segmentation
  "Segmentation: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 25 -- the image cut into regions that
  are each roughly one thing.

    graph-based   Felzenszwalb and Huttenlocher's method: the pixels a
                  graph, each joined to its neighbors by an edge weighted
                  by how different they are; the edges taken cheapest
                  first, two regions merged when the edge between them is
                  no bigger than the larger difference inside either,
                  plus k / size -- so small regions merge readily and
                  large ones only across a weak boundary. A step an edge.
    clustering    k-means on the pixels' values, the regions whatever
                  falls together: no notion of neighbors, so it is
                  what the chapter means by segmenting by color alone

  Regions are kept by `allgo.graph`'s disjoint sets, joined smaller into
  larger so their trees stay shallow; the clustering is
  `allgo.learning.unsupervised`'s k-means."
  (:require [allgo.graph :as graph]
            [allgo.learning.unsupervised :as unsupervised]
            [allgo.vision.image :as img]
            [clojure.math :as math]))

(defn- difference [a b]
  (if (vector? a)
    (math/sqrt (reduce + (map #(let [d (- %1 %2)] (* d d)) a b)))
    (abs (- a b))))

(defn pixel-edges
  "Every pair of 8-neighboring pixels, `[weight p q]`, weighted by how
  different they are, cheapest first."
  [{:keys [w h] :as image}]
  (sort-by first
           (for [y (range h) x (range w)
                 [dx dy] [[1 0] [0 1] [1 1] [1 -1]]
                 :let [x2 (+ x dx) y2 (+ y dy)]
                 :when (and (< x2 w) (< -1 y2 h))]
             [(difference (img/at image x y) (img/at image x2 y2)) [x y] [x2 y2]])))

(defn graph-based-seq
  "Felzenszwalb-Huttenlocher segmentation of `image` (best smoothed a
  little first), an edge a step: each state `{:forest :size :internal}`,
  the disjoint sets and each region's size and largest internal edge,
  keyed by its representative. `k` sets the scale: larger, larger
  regions."
  [image k]
  (let [pixels (for [y (range (:h image)) x (range (:w image))] [x y])]
    (reductions
     (fn [{:keys [forest size internal] :as st} [wt p q]]
       (let [a (graph/canonical forest p) b (graph/canonical forest q)]
         (if (and (not= a b)
                  (<= wt (min (+ (internal a) (/ k (size a))) (+ (internal b) (/ k (size b))))))
           ;; smaller into larger: `union` points the first's root at the
           ;; second's
           (let [[small big] (if (< (size a) (size b)) [a b] [b a])]
             {:forest (graph/union forest small big)
              :size (assoc size big (+ (size a) (size b)))
              :internal (assoc internal big (max wt (internal a) (internal b)))})
           st)))
     {:forest (graph/disjoint-sets pixels)
      :size (zipmap pixels (repeat 1))
      :internal (zipmap pixels (repeat 0.0))}
     (pixel-edges image))))

(defn labels
  "The region of every pixel as an image of region numbers, 0 up, in the
  order regions are first met reading the image."
  [{:keys [forest]} {:keys [w h]}]
  (let [roots (for [y (range h) x (range w)] (graph/canonical forest [x y]))
        number (reduce (fn [m r] (if (contains? m r) m (assoc m r (count m)))) {} roots)]
    {:w w :h h :px (mapv number roots)}))

(defn merge-small
  "The segmentation with regions smaller than `min-size` merged into a
  neighbor, cheapest boundary first -- the usual tidying after the main
  pass, which leaves specks of noise as regions of their own."
  [state image min-size]
  (reduce (fn [{:keys [forest size] :as st} [_ p q]]
            (let [a (graph/canonical forest p) b (graph/canonical forest q)]
              (if (and (not= a b) (or (< (size a) min-size) (< (size b) min-size)))
                (let [[small big] (if (< (size a) (size b)) [a b] [b a])]
                  (-> st
                      (assoc :forest (graph/union forest small big))
                      (assoc-in [:size big] (+ (size a) (size b)))))
                st)))
          state (pixel-edges image)))

(defn graph-based
  "The regions of `image` at scale `k`, specks under `min-size` merged:
  an image of region numbers."
  [image k min-size]
  (labels (merge-small (last (graph-based-seq image k)) image min-size) image))

(defn k-means
  "The pixels of `image` in `k` clusters of value: an image of cluster
  numbers. k-means settles in whichever local optimum its first centers
  lead to -- two of them in one large region, say -- so it is run
  `restarts` times (default 5) for `steps` iterations each (default 10)
  and the run of least cost kept."
  ([image k rng] (k-means image k rng 5 10))
  ([{:keys [w h px]} k rng restarts steps]
   (let [points (mapv #(if (vector? %) % [%]) px)
         best (apply min-key :cost (repeatedly restarts #(nth (unsupervised/k-means-seq points k rng) steps)))]
     {:w w :h h :px (:assignment best)})))

(defn region-count [{:keys [px]}] (count (distinct px)))
