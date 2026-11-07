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
            [allgo.numerics.eigen :as eigen]
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

;; ---------------------------------------------------------------------------
;; Normalized cuts

(defn- mean-value [vals]
  (if (vector? (first vals))
    (apply mapv (fn [& xs] (/ (reduce + xs) (count xs))) vals)
    (/ (reduce + vals) (count vals))))

(defn superpixel-graph
  "`image` over-segmented into superpixels -- graph-based segmentation at
  the small scale `k`, specks under `min-size` merged -- as a weighted
  graph for normalized cuts: the chapter's way of making a graph method
  affordable, a few dozen nodes in place of thousands of pixels. Two
  superpixels that touch are joined with weight the length of the
  boundary they share times exp(-|Fi - Fj|^2 / sigma^2), Fi a superpixel's
  mean value: alike and long in contact, strongly. `{:nodes :w}`, each
  node `{:pixels :mean}` and the weights a symmetric matrix."
  [image k min-size sigma]
  (let [{:keys [w h px]} (graph-based image k min-size)
        groups (group-by #(px (+ (* (second %) w) (first %))) (for [y (range h) x (range w)] [x y]))
        ids (vec (sort (keys groups)))
        index (zipmap ids (range))
        nodes (mapv (fn [id] (let [ps (groups id)] {:pixels ps :mean (mean-value (map (fn [[x y]] (img/at image x y)) ps))})) ids)
        contacts (frequencies (for [y (range h) x (range w) [dx dy] [[1 0] [0 1]]
                                    :let [x2 (+ x dx) y2 (+ y dy)]
                                    :when (and (< x2 w) (< y2 h))
                                    :let [a (index (px (+ (* y w) x))) b (index (px (+ (* y2 w) x2)))]
                                    :when (not= a b)]
                                (if (< a b) [a b] [b a])))
        n (count nodes)]
    {:nodes nodes
     :w (reduce (fn [m [[a b] len]]
                  (let [d (difference (:mean (nodes a)) (:mean (nodes b)))
                        wt (* len (math/exp (- (/ (* d d) (* sigma sigma)))))]
                    (-> m (assoc-in [a b] wt) (assoc-in [b a] wt))))
                (vec (repeat n (vec (repeat n 0.0))))
                contacts)}))

(defn- subgraph
  "The graph on just the nodes at `indices`."
  [{:keys [nodes w]} indices]
  (let [is (vec indices)]
    {:nodes (mapv nodes is) :w (mapv (fn [i] (mapv #(get-in w [i %]) is)) is)}))

(defn fiedler
  "Shi and Malik's relaxation of the best normalized cut of the graph of
  weights `w`: the generalized eigenvector y of (D - W) y = lambda D y
  with the second smallest eigenvalue -- the smallest, 0, has y constant.
  Found as D^-1/2 times the eigenvector of D^-1/2 W D^-1/2 with the
  second largest eigenvalue, by `allgo.numerics.eigen`'s Jacobi."
  [w]
  (let [d (mapv #(reduce + %) w)
        s (mapv #(math/sqrt (max % 1e-300)) d)
        m (vec (map-indexed (fn [i row] (vec (map-indexed (fn [j x] (/ x (* (s i) (s j)))) row))) w))
        {:keys [vectors]} (eigen/jacobi m)
        n (count w)]
    (mapv (fn [i] (/ (get-in vectors [i (- n 2)]) (s i))) (range n))))

(defn best-cut
  "The best normalized cut of the graph of weights `w` with the nodes
  taken in the order of `y`: every split of that order into a first part
  and the rest tried, Ncut(A, B) = cut(A, B) / assoc(A, V) + cut(A, B) /
  assoc(B, V) kept up to date as each node moves across. `{:a :b :ncut}`,
  the parts as vectors of indices."
  [w y]
  (let [d (mapv #(reduce + %) w)
        total (reduce + d)
        order (vec (sort-by y (range (count y))))
        [best] (reduce (fn [[best in-a cut assoc-a] [k i]]
                         (let [to-a (reduce + (map #(get-in w [i %]) in-a))
                               cut (+ cut (- (d i) to-a) (- to-a))
                               assoc-a (+ assoc-a (d i))
                               ncut (if (< 0 assoc-a total) (+ (/ cut assoc-a) (/ cut (- total assoc-a))) ##Inf)]
                           [(if (< ncut (:ncut best)) {:ncut ncut :split (inc k)} best) (conj in-a i) cut assoc-a]))
                       [{:ncut ##Inf :split 0} [] 0.0 0.0]
                       (map-indexed vector (pop order)))]
    {:ncut (:ncut best) :a (subvec order 0 (:split best)) :b (subvec order (:split best))}))

(defn normalized-cuts-seq
  "Segmentation of `image` by recursive normalized cuts, after Shi and
  Malik (\"Normalized Cuts and Image Segmentation\", IEEE PAMI 22, 2000),
  on its superpixels (`superpixel-graph` at `k`, `min-size` and
  `sigma`), a cut a step: a graph cut in two where the second eigenvector
  gives the least Ncut, and each part again, as long as the cut's Ncut is
  under `max-ncut`. Each state `{:pending :regions}`, the graphs still to
  try and the finished regions, each a seq of `[x y]`; until none is
  pending."
  [image {:keys [k min-size sigma max-ncut] :or {k 0.05 min-size 8 sigma 0.1 max-ncut 0.05}}]
  (let [pixels (fn [graph] (mapcat :pixels (:nodes graph)))
        step (fn [{:keys [pending regions]}]
               (when-let [graph (peek pending)]
                 (let [pending (pop pending)
                       done {:pending pending :regions (conj regions (pixels graph))}]
                   (if (< (count (:nodes graph)) 2)
                     done
                     (let [{:keys [ncut a b]} (best-cut (:w graph) (fiedler (:w graph)))]
                       (if (> ncut max-ncut)
                         done
                         {:pending (conj pending (subgraph graph a) (subgraph graph b)) :regions regions}))))))]
    (take-while some? (iterate step {:pending [(superpixel-graph image k min-size sigma)] :regions []}))))

(defn normalized-cuts
  "The regions `normalized-cuts-seq` finds: an image of region numbers."
  [image opts]
  (let [number (into {} (for [[k ps] (map-indexed vector (:regions (last (normalized-cuts-seq image opts)))) p ps] [p k]))]
    (img/image (:w image) (:h image) (fn [x y] (number [x y])))))
