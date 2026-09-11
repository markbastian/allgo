(ns allgo.graph
  "Graph algorithms on undirected graphs.

  A graph here is whatever you already have: a collection of nodes and a
  collection of edges, where an edge is a two-element set `#{u v}`. Nodes
  can be any value that works as a map key. Weights are supplied as a
  function of an edge rather than stored on it, so the same edge set can be
  measured different ways without being rebuilt."
  (:require [clojure.set :as set]))

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
  "`edges` as a map from each node to the set of its neighbours. Nodes with
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
