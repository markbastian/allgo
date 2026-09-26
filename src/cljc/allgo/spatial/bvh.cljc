(ns allgo.spatial.bvh
  "A bounding volume hierarchy built by sorting. After Ten Minute
  Physics 24.

  A BVH is a binary tree of boxes: each node's box contains its children's,
  so a query that misses a node misses everything under it and a whole
  subtree is skipped on one test. The question is how to decide which
  objects go together, and the answer here costs almost nothing: sort the
  objects by `allgo.spatial.morton` code and split the sorted list down
  the middle, over and over.

  That works because Morton order keeps things that are near each other
  near each other in the list, so any contiguous run is already a
  reasonably compact cluster. No searching for a good split, no surface
  area heuristic, no comparisons at all beyond the sort -- which is itself
  a radix sort. It is how large trees are built on a GPU, and the reason
  this one is worth having over the alternatives is that it can be rebuilt
  from scratch every frame rather than updated.

  The tree is flat arrays rather than nodes: `left`, `right` and `box`
  indices, and the bounds six doubles at a time. A tree of `n` leaves has
  exactly `2n - 1` nodes when every split halves, which is what makes the
  arrays a fixed size known up front.

  Where it is weaker than `allgo.spatial.sweep` or `allgo.spatial.hash`:
  the Z curve jumps at every power-of-two boundary, so two objects either
  side of one are far apart in the list however close they are in space.
  That shows up as the occasional stretched node, and is the price of
  building a tree by sorting rather than by searching."
  (:require [allgo.array :as a]
            [allgo.spatial.morton :as morton]))

(def ^:private ^:const no-box -1)

(defn- centers
  "The center of every box, flat, 3 per box -- what the Morton code is
  taken of."
  [^doubles mins ^doubles maxs ^long n]
  (let [out (double-array (* 3 n))]
    (dotimes [i n]
      (let [b (* 3 i)]
        (dotimes [k 3]
          (aset out (+ b k) (* 0.5 (+ (aget mins (+ b k)) (aget maxs (+ b k))))))))
    out))

(defn build
  "A tree over `n` boxes, given as flat `mins` and `maxs`, 3 per box.

  Returns `{:n :root :left :right :box :lo :hi}`. `box` is the object a
  leaf holds, or -1 for an internal node."
  [mins maxs ^long n]
  (let [^doubles mins mins
        ^doubles maxs maxs]
    (if (zero? n)
      ;; The reference divides an empty range and recurses forever. An
      ;; empty tree is a legitimate thing to ask for and a legitimate
      ;; thing to query -- it just never hits.
      {:n 0 :root -1 :left (a/i32 0) :right (a/i32 0) :box (a/i32 0)
       :lo (double-array 0) :hi (double-array 0)}
      (let [nodes (dec (* 2 n))
            ^ints left  (a/i32 nodes)
            ^ints right (a/i32 nodes)
            ^ints box   (a/i32 nodes)
            ^doubles lo (double-array (* 3 nodes))
            ^doubles hi (double-array (* 3 nodes))
            ^ints order (morton/sorted-points (centers mins maxs n) n)
            next-node (volatile! 0)
            take-node! (fn [] (let [k @next-node] (vreset! next-node (inc k)) k))]
        (letfn [(leaf! [id]
                  (let [k (take-node!)
                        b (* 3 id)
                        c (* 3 k)]
                    (aset left k (int -1))
                    (aset right k (int -1))
                    (aset box k (int id))
                    (dotimes [t 3]
                      (aset lo (+ c t) (aget mins (+ b t)))
                      (aset hi (+ c t) (aget maxs (+ b t))))
                    k))
                (subtree! [begin end]
                  ;; Inclusive range, as the reference has it.
                  (if (= begin end)
                    (leaf! (aget order begin))
                    (let [mid (quot (+ begin end) 2)
                          l (subtree! begin mid)
                          r (subtree! (inc mid) end)
                          k (take-node!)
                          c (* 3 k) cl (* 3 l) cr (* 3 r)]
                      (aset left k (int l))
                      (aset right k (int r))
                      (aset box k (int no-box))
                      (dotimes [t 3]
                        (aset lo (+ c t) (min (aget lo (+ cl t)) (aget lo (+ cr t))))
                        (aset hi (+ c t) (max (aget hi (+ cl t)) (aget hi (+ cr t)))))
                      k)))]
          ;; Children are taken before their parent, so a node's index is
          ;; always above its children's and the root is the last one.
          (let [root (subtree! 0 (dec n))]
            {:n n :root root :left left :right right :box box :lo lo :hi hi}))))))

(defn node-count [{:keys [n]}] (if (zero? (long n)) 0 (dec (* 2 (long n)))))

(defn leaf? [{:keys [^ints box]} k] (not= no-box (aget box k)))

(defn depth
  "The deepest path in the tree. A halving split gives ceil(log2 n) + 1,
  and anything much worse means the split is not halving."
  [{:keys [root ^ints left ^ints right] :as t}]
  (if (neg? (long root))
    0
    (letfn [(go [k] (if (leaf? t k)
                      1
                      (inc (max (go (aget left k)) (go (aget right k))))))]
      (go root))))

(defn- node-hits?
  "Does node `k`'s box overlap the query box?

  The query arrives as one array of six -- min x y z then max x y z --
  rather than six arguments, because a function taking any primitive
  argument may have at most four of them, and this is the test every
  traversal runs at every node."
  [^doubles lo ^doubles hi ^long k ^doubles q]
  (let [c (* 3 k)]
    (and (<= (aget lo c) (aget q 3)) (>= (aget hi c) (aget q 0))
         (<= (aget lo (+ c 1)) (aget q 4)) (>= (aget hi (+ c 1)) (aget q 1))
         (<= (aget lo (+ c 2)) (aget q 5)) (>= (aget hi (+ c 2)) (aget q 2)))))

(defn query
  "Every box whose bounds overlap the given region, as a vector of object
  indices.

  The region is `[x0 y0 z0]` to `[x1 y1 z1]`."
  [{:keys [root ^ints left ^ints right ^ints box ^doubles lo ^doubles hi] :as t}
   [x0 y0 z0] [x1 y1 z1]]
  (if (neg? (long root))
    []
    (let [q (double-array [x0 y0 z0 x1 y1 z1])]
      (persistent!
       (loop [stack [root] acc (transient [])]
         (if (empty? stack)
           acc
           (let [k (peek stack)
                 stack (pop stack)]
             (if-not (node-hits? lo hi k q)
               (recur stack acc)
               (if (leaf? t k)
                 (recur stack (conj! acc (aget box k)))
                 (recur (conj stack (aget left k) (aget right k)) acc))))))))))

(defn query-point [t [x y z]] (query t [x y z] [x y z]))

(defn- overlap?
  [^doubles mins ^doubles maxs ^long i ^long j]
  (let [a (* 3 i) b (* 3 j)]
    (and (<= (aget mins a) (aget maxs b)) (>= (aget maxs a) (aget mins b))
         (<= (aget mins (+ a 1)) (aget maxs (+ b 1)))
         (>= (aget maxs (+ a 1)) (aget mins (+ b 1)))
         (<= (aget mins (+ a 2)) (aget maxs (+ b 2)))
         (>= (aget maxs (+ a 2)) (aget mins (+ b 2))))))

(defn overlapping-pairs
  "Every pair of overlapping boxes, each pair once, as `[i j]` with
  `i < j`.

  The same shape `allgo.spatial.hash/overlapping-pairs` and
  `allgo.spatial.sweep/overlapping-pairs` return.

  Each object is queried against the tree, and only partners numbered
  above it are kept. The reference keeps both, so every pair comes back
  twice -- while the brute force loop it is compared against reports each
  once, which makes the two counts it prints side by side differ by a
  factor of two for reasons that have nothing to do with the tree."
  [{:keys [root ^ints left ^ints right ^ints box ^doubles lo ^doubles hi] :as t}
   mins maxs ^long n]
  (let [^doubles mins mins
        ^doubles maxs maxs]
    (if (neg? (long root))
      []
      (persistent!
       (reduce
        (fn [acc i]
          (let [b (* 3 i)
                q (double-array [(aget mins b) (aget mins (+ b 1)) (aget mins (+ b 2))
                                 (aget maxs b) (aget maxs (+ b 1)) (aget maxs (+ b 2))])]
            (loop [stack [root] acc acc]
              (if (empty? stack)
                acc
                (let [k (peek stack)
                      stack (pop stack)]
                  (if-not (node-hits? lo hi k q)
                    (recur stack acc)
                    (if (leaf? t k)
                      (let [j (aget box k)]
                        (recur stack
                               (if (and (> j i) (overlap? mins maxs i j))
                                 (conj! acc [i j])
                                 acc)))
                      (recur (conj stack (aget left k) (aget right k)) acc))))))))
        (transient [])
        (range n))))))

(defn adjacency
  "Every overlapping pair as a compressed adjacency list, laid out as
  `allgo.spatial.hash/adjacency` lays it out."
  [t mins maxs ^long n]
  (let [pairs (overlapping-pairs t mins maxs n)
        by-hi (reduce (fn [m [i j]] (update m j (fnil conj []) i)) {} pairs)
        ^ints starts (a/i32 (inc n))]
    (loop [i 0 acc (transient []) total 0]
      (if (= i n)
        (do (aset starts n (int total))
            {:starts starts
             :ids (let [v (persistent! acc)
                        ^ints out (a/i32 (count v))]
                    (dotimes [k (count v)] (aset out k (int (v k))))
                    out)
             :pairs total})
        (let [_ (aset starts i (int total))
              partners (sort (get by-hi i))]
          (recur (inc i) (reduce conj! acc partners) (+ total (count partners))))))))
