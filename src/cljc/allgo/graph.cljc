(ns allgo.graph
  "Graph algorithms on undirected graphs.

  A graph here is whatever you already have: a collection of nodes and a
  collection of edges, where an edge is a two-element set `#{u v}`. Nodes
  can be any value that works as a map key. Weights are supplied as a
  function of an edge rather than stored on it, so the same edge set can be
  measured different ways without being rebuilt."
  (:require [allgo.search :as search]
            [clojure.set :as set]))

;; ---------------------------------------------------------------------------
;; Disjoint sets

(defn disjoint-sets
  "A disjoint-set forest with every node in its own set."
  [nodes]
  (zipmap nodes nodes))

(defn canonical
  "The representative of the set containing `x`."
  [forest x]
  (loop [x x]
    (let [parent (forest x)]
      (if (or (nil? parent) (= parent x)) x (recur parent)))))

(defn union
  "Merges the sets containing `x` and `y`."
  [forest x y]
  (let [rx (canonical forest x)
        ry (canonical forest y)]
    (cond-> forest (not= rx ry) (assoc rx ry))))

(defn connected?
  "Whether `x` and `y` are in the same set."
  [forest x y]
  (= (canonical forest x) (canonical forest y)))

;; ---------------------------------------------------------------------------

(defn adjacency
  "`edges` as a map from each node to the set of its neighbors. Nodes with
  no edges are present only if listed in `nodes`."
  ([edges] (adjacency nil edges))
  ([nodes edges]
   (reduce (fn [adj e]
             (let [[u v] (vec e)]
               (-> adj
                   (update u (fnil conj #{}) v)
                   (update v (fnil conj #{}) u))))
           (zipmap nodes (repeat #{}))
           edges)))

(defn components
  "The connected components of the graph, as a set of sets of nodes."
  [nodes edges]
  (let [forest (reduce (fn [f e] (let [[u v] (vec e)] (union f u v)))
                       (disjoint-sets nodes)
                       edges)]
    (->> nodes
         (group-by #(canonical forest %))
         vals
         (map set)
         set)))

(defn components-by
  "The connected components of `nodes` under `neighbors`, a function from
  a node to its neighbors: a vector of sets, in the order of their first
  node in `nodes`, each filled breadth first. For graphs held as a
  neighbor function -- a roadmap, a lattice -- where gathering the edges
  for `components` would cost more than the search."
  [nodes neighbors]
  (let [fill (fn [n] (-> (search/breadth-first-seq {:start n :goal? (constantly false) :neighbors neighbors})
                         last :came-from keys set))]
    (second (reduce (fn [[seen comps :as acc] n]
                      (if (seen n) acc (let [c (fill n)] [(into seen c) (conj comps c)])))
                    [#{} []]
                    nodes))))

(defn connected-graph?
  "Whether every node is reachable from every other."
  [nodes edges]
  (<= (count (components nodes edges)) 1))

(defn minimum-spanning-tree
  "Kruskal's minimum spanning tree: the cheapest subset of `edges` that
  keeps every node in `nodes` as connected as the full graph did.

  Sorts the edges by `weight` and takes each one that joins two pieces not
  already linked. On a disconnected graph this yields a spanning forest --
  one tree per component -- rather than failing."
  [nodes edges weight]
  (loop [[e & more] (sort-by weight edges)
         forest     (disjoint-sets nodes)
         tree       #{}]
    (if (nil? e)
      tree
      (let [[u v] (vec e)]
        (if (connected? forest u v)
          (recur more forest tree)
          (recur more (union forest u v) (conj tree e)))))))

(defn complete-graph
  "Every unordered pair of `nodes`, as edges."
  [nodes]
  (let [nodes (vec nodes)]
    (set (for [i (range (count nodes))
               j (range (inc i) (count nodes))]
           #{(nodes i) (nodes j)}))))

(defn edges-not-in
  "The edges of a graph that its spanning `tree` left out."
  [edges tree]
  (set/difference (set edges) (set tree)))

;; ---------------------------------------------------------------------------
;; Assignment

(defn- relax-step
  "One step of the shortest augmenting path's search from row `i`, as in
  Kuhn's Hungarian method with potentials: mark column j0 reached, lower
  each unreached column's reduced cost through it, move to the cheapest,
  and shift the potentials by that cost. nil once a free column is
  reached."
  [cost {:keys [u v p way minv used j0] :as st}]
  (when-not (and (pos? j0) (zero? (p j0)))
    (let [used (assoc used j0 true)
          i0 (p j0)
          m (dec (count v))
          [minv way delta j1]
          (reduce (fn [[minv way delta j1 :as acc] j]
                    (if (used j)
                      acc
                      (let [cur (- (get-in cost [(dec i0) (dec j)]) (u i0) (v j))
                            [minv way] (if (< cur (minv j)) [(assoc minv j cur) (assoc way j j0)] [minv way])]
                        (if (< (minv j) delta) [minv way (minv j) j] [minv way delta j1]))))
                  [minv way ##Inf 0]
                  (range 1 (inc m)))
          [u v minv] (reduce (fn [[u v minv] j]
                               (if (used j)
                                 [(update u (p j) + delta) (update v j - delta) minv]
                                 [u v (update minv j - delta)]))
                             [u v minv]
                             (range (inc m)))]
      (assoc st :u u :v v :minv minv :way way :used used :j0 j1))))

(defn- augment-step
  "Shift the assignments back along the augmenting path by one column."
  [{:keys [p way j0] :as st}]
  (when (pos? j0)
    (let [j1 (way j0)]
      (assoc st :p (assoc p j0 (p j1)) :j0 j1))))

(defn- assign-row
  "Row `i` brought into the assignment: a shortest augmenting path to a
  free column, then the assignments shifted along it."
  [cost {:keys [v] :as st} i]
  (let [m (dec (count v))
        st (assoc st :p (assoc (:p st) 0 i) :j0 0
                  :minv (vec (repeat (inc m) ##Inf)) :used (vec (repeat (inc m) false)))
        st (last (search/step-seq (partial relax-step cost) st))]
    (last (search/step-seq augment-step st))))

(defn min-cost-assignment
  "The assignment of each row of `cost` (n rows of m >= n costs) to a
  distinct column of least total cost, by Kuhn's Hungarian method (Naval
  Research Logistics Quarterly 2, 1955) with potentials, O(n^2 m): the
  column of each row, as a vector."
  [cost]
  (let [n (count cost) m (count (first cost))
        {:keys [p]} (reduce (partial assign-row cost)
                            {:u (vec (repeat (inc n) 0.0)) :v (vec (repeat (inc m) 0.0))
                             :p (vec (repeat (inc m) 0)) :way (vec (repeat (inc m) 0))}
                            (range 1 (inc n)))]
    (reduce (fn [out j] (if (pos? (p j)) (assoc out (dec (p j)) (dec j)) out))
            (vec (repeat n nil))
            (range 1 (inc m)))))
