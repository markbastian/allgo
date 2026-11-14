(ns allgo.spatial.kdtree
  "k-d trees: LaValle, *Planning Algorithms*, Cambridge University Press,
  2006, section 5.5.2 -- finding the nearest of many points without
  measuring them all, which a sampling planner does at every step.

  The points are split at the median of one coordinate, then each half at
  the median of the next, and so on round the coordinates (Bentley,
  \"Multidimensional Binary Search Trees Used for Associative
  Searching\", CACM 18, 1975). A query descends to the leaf its point
  would lie in, then backs out, visiting the far side of a split only if
  the splitting plane is nearer than the best found so far.

  Points are vectors of numbers; a tree is `{:point :axis :left :right}`
  or nil, each node holding the point it splits at."
  (:require [allgo.numerics.linear :as lin]))

(defn build
  "The k-d tree of `points`."
  ([points] (build (vec points) 0))
  ([points depth]
   (when (seq points)
     (let [k (count (first points))
           axis (mod depth k)
           sorted (vec (sort-by #(nth % axis) points))
           m (quot (count sorted) 2)]
       {:point (sorted m) :axis axis
        :left (build (subvec sorted 0 m) (inc depth))
        :right (build (subvec sorted (inc m)) (inc depth))}))))

(defn nearest
  "The point of `tree` nearest `q`, and its squared distance: `[point
  d2]`; in the metadata, the points the search measured, in order,
  `:measured`, and how many, `:visited`."
  [tree q]
  (let [measured (volatile! [])
        search (fn search [{:keys [point axis left right] :as node} [_ best-d2 :as best]]
                 (if (nil? node)
                   best
                   (let [_ (vswap! measured conj point)
                         d2 (lin/distance-squared point q)
                         best (if (< d2 best-d2) [point d2] best)
                         diff (- (nth q axis) (nth point axis))
                         [near far] (if (neg? diff) [left right] [right left])
                         best (search near best)]
                     (if (< (* diff diff) (second best)) (search far best) best))))]
    (with-meta (search tree [nil ##Inf]) {:visited (count @measured) :measured @measured})))

(defn k-nearest
  "The `k` points of `tree` nearest `q`, nearest first; the points
  measured, `:measured` in the metadata."
  [tree q k]
  (let [measured (volatile! [])
        search (fn search [{:keys [point axis left right] :as node} best]
                 (if (nil? node)
                   best
                   (let [_ (vswap! measured conj point)
                         best (vec (take k (sort-by second (conj best [point (lin/distance-squared point q)]))))
                         diff (- (nth q axis) (nth point axis))
                         [near far] (if (neg? diff) [left right] [right left])
                         best (search near best)
                         worst (if (< (count best) k) ##Inf (second (peek best)))]
                     (if (< (* diff diff) worst) (search far best) best))))]
    (with-meta (mapv first (search tree [])) {:measured @measured})))

(defn within
  "Every point of `tree` within `r` of `q`; the points measured,
  `:measured` in the metadata."
  [tree q r]
  (let [r2 (* r r)
        measured (volatile! [])
        search (fn search [{:keys [point axis left right] :as node}]
                 (when node
                   (vswap! measured conj point)
                   (let [diff (- (nth q axis) (nth point axis))
                         [near far] (if (neg? diff) [left right] [right left])]
                     (concat (when (<= (lin/distance-squared point q) r2) [point])
                             (search near)
                             (when (<= (* diff diff) r2) (search far))))))
        ;; realized before the measured points are read
        found (vec (search tree))]
    (with-meta found {:measured @measured})))
