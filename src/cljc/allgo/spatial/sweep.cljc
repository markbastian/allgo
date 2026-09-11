(ns allgo.spatial.sweep
  "Sweep and prune: sort the boxes along one axis, then only compare the
  ones that overlap on it. After Ten Minute Physics 23.

  Two boxes can only touch if they overlap on every axis, so overlapping
  on one axis is a cheap necessary condition. Sort the boxes by their
  lower edge and walk the list; for each box, the only candidates are the
  ones that start before it ends, and they are the next few entries. The
  moment one starts later, so does every box after it, and the scan stops.

  How this compares with `allgo.spatial.hash`:

    hash    bins by position. Best when objects are all about one size and
            spread out -- the query cost does not care how they are
            arranged, only how many share a cell.
    sweep   sorts by extent. Handles wildly different sizes without
            having to pick a cell size for them, and costs nothing to
            maintain when things move a little, because a nearly-sorted
            list is nearly free to re-sort.

  Both return the same thing, so they can be swapped: `overlapping-pairs`
  gives `[i j]` pairs with `i < j`, once each.

  Two things the reference does not do. It sorts the whole array every
  frame, where the point of sweep and prune is that objects barely move
  between frames, so the order is nearly right already -- keeping it and
  running an insertion pass costs about the number of things that actually
  changed places, not n log n. And it always sweeps x, where sweeping the
  axis the boxes are most spread along is what decides how many candidates
  survive the prune; on a tall thin scene, x is the worst of the three."
  (:require [allgo.array :as a]))

(defn sweep
  "Room for `max-objects` boxes.

  The order is kept between calls; that is the whole point."
  [max-objects]
  {:max-objects (long max-objects)
   :order       (a/i32 max-objects)
   ;; Which axis was swept last, and over how many objects, so that a
   ;; change in either is known to invalidate the kept order.
   :axis        (volatile! -1)
   :count       (volatile! 0)})

(defn axis-of
  "The axis the boxes are most spread along, 0, 1 or 2.

  Spread is measured over the box centres, by variance. Sweeping the
  widest axis is what makes the prune bite: on a scene twice as tall as it
  is wide, sweeping x leaves roughly twice as many candidates per box as
  sweeping y does."
  [^doubles mins ^doubles maxs ^long n]
  (if (zero? n)
    0
    (let [sums (double-array 3)
          sq   (double-array 3)]
      (dotimes [i n]
        (let [b (* 3 i)]
          (dotimes [k 3]
            (let [c (* 0.5 (+ (aget mins (+ b k)) (aget maxs (+ b k))))]
              (aset sums k (+ (aget sums k) c))
              (aset sq k (+ (aget sq k) (* c c)))))))
      (let [variance (fn [k]
                       (let [mean (/ (aget sums k) n)]
                         (- (/ (aget sq k) n) (* mean mean))))]
        (long (first (apply max-key second
                            (map-indexed (fn [k _] [k (variance k)]) (range 3)))))))))

(defn sort!
  "Puts the kept order back in order, by insertion.

  Insertion sort is quadratic on a shuffled list and linear on one that is
  nearly right, which is exactly the case here: between frames a box moves
  a little and swaps with a neighbour or two. The first call on a fresh
  structure is the shuffled case, so it starts from the identity order,
  which for boxes built in any spatial order is already close.

  Returns the axis swept."
  [{:keys [^ints order axis count]} ^doubles mins ^doubles maxs ^long n]
  (let [a (axis-of mins maxs n)]
    ;; A different axis, or a different set of objects, makes the kept
    ;; order meaningless -- start from identity rather than insertion
    ;; sorting something arbitrary.
    (when (or (not= a @axis) (not= n @count))
      (dotimes [i n] (aset order i (int i)))
      (vreset! axis a)
      (vreset! count n))
    (let [key-of (fn ^double [^long id] (aget mins (+ (* 3 id) a)))]
      (loop [i 1]
        (when (< i n)
          (let [id (aget order i)
                k  (key-of id)]
            (loop [j (dec i)]
              (if (and (>= j 0) (> (key-of (aget order j)) k))
                (do (aset order (inc j) (aget order j))
                    (recur (dec j)))
                (aset order (inc j) (int id)))))
          (recur (inc i)))))
    a))

(defn- overlap?
  "Do two boxes overlap on all three axes?"
  [^doubles mins ^doubles maxs ^long i ^long j]
  (let [a (* 3 i) b (* 3 j)]
    (and (<= (aget mins a) (aget maxs b))
         (>= (aget maxs a) (aget mins b))
         (<= (aget mins (+ a 1)) (aget maxs (+ b 1)))
         (>= (aget maxs (+ a 1)) (aget mins (+ b 1)))
         (<= (aget mins (+ a 2)) (aget maxs (+ b 2)))
         (>= (aget maxs (+ a 2)) (aget mins (+ b 2))))))

(defn overlapping-pairs
  "Every pair of overlapping boxes, each pair once, as `[i j]` with
  `i < j`.

  The same shape `allgo.spatial.hash/overlapping-pairs` returns, so the
  two are interchangeable."
  [{:keys [^ints order] :as s} mins maxs ^long n]
  (let [^doubles mins mins
        ^doubles maxs maxs
        a (sort! s mins maxs n)]
    (persistent!
     (loop [i 0 acc (transient [])]
       (if (>= i n)
         acc
         (let [id (aget order i)
               hi (aget maxs (+ (* 3 id) a))]
           (recur
            (inc i)
            (loop [j (inc i) acc acc]
              (if (>= j n)
                acc
                (let [jd (aget order j)]
                  ;; Sorted by lower edge, so once one box starts after
                  ;; this one ends, so does every box after it.
                  (if (> (aget mins (+ (* 3 jd) a)) hi)
                    acc
                    (recur (inc j)
                           (if (overlap? mins maxs id jd)
                             (conj! acc [(min id jd) (max id jd)])
                             acc)))))))))))))

(defn adjacency
  "Every overlapping pair as a compressed adjacency list.

  `{:starts :ids :pairs}`, laid out the same way
  `allgo.spatial.hash/adjacency` lays it out: object `i`'s partners are
  `ids[starts[i] .. starts[i+1])`, each pair recorded once under the
  higher-numbered object."
  [s mins maxs ^long n]
  (let [pairs (overlapping-pairs s mins maxs n)
        by-hi (reduce (fn [m [i j]] (update m j (fnil conj []) i)) {} pairs)
        ^ints starts (a/i32 (inc n))]
    (loop [i 0 acc (transient []) total 0]
      (if (= i n)
        (do (aset starts n (int total))
            {:starts starts
             :ids    (let [v (persistent! acc)
                           ^ints out (a/i32 (count v))]
                       (dotimes [k (count v)] (aset out k (int (v k))))
                       out)
             :pairs  total})
        (let [_ (aset starts i (int total))
              partners (sort (get by-hi i))]
          (recur (inc i)
                 (reduce conj! acc partners)
                 (+ total (count partners))))))))

;; ---------------------------------------------------------------------------
;; Building boxes

(defn spheres->boxes!
  "Fills `mins` and `maxs` from sphere centres and a radius.

  Lets a scene held as points be swept without being rewritten as boxes,
  which is how `allgo.spatial.hash` is usually fed."
  ([positions n radius] (spheres->boxes! positions n radius nil nil))
  ;; No primitive hint on `n`: five arguments is one past what a
  ;; primitive-taking function may have.
  ([^doubles positions n radius mins maxs]
   (let [n (long n)
         ^doubles mins (or mins (double-array (* 3 n)))
         ^doubles maxs (or maxs (double-array (* 3 n)))
         r (double radius)]
     (dotimes [i n]
       (let [b (* 3 i)]
         (dotimes [k 3]
           (aset mins (+ b k) (- (aget positions (+ b k)) r))
           (aset maxs (+ b k) (+ (aget positions (+ b k)) r)))))
     [mins maxs])))

(defn boxes-of
  "`[mins maxs]` for `n` boxes given per-object half-extents, flat."
  [^doubles positions ^doubles half-extents ^long n]
  (let [mins (double-array (* 3 n))
        maxs (double-array (* 3 n))]
    (dotimes [i n]
      (let [b (* 3 i)]
        (dotimes [k 3]
          (aset mins (+ b k) (- (aget positions (+ b k)) (aget half-extents (+ b k))))
          (aset maxs (+ b k) (+ (aget positions (+ b k)) (aget half-extents (+ b k)))))))
    [mins maxs]))

(defn inversions
  "How far out of order the kept list is, for checking that coherence is
  doing what it is supposed to.

  Counted as adjacent swaps, which is what an insertion pass costs."
  [{:keys [^ints order axis]} ^doubles mins ^long n]
  (let [a (long (max 0 @axis))]
    (loop [i 1 total 0]
      (if (>= i n)
        total
        (let [k (aget mins (+ (* 3 (aget order i)) a))]
          (recur (inc i)
                 (+ total
                    (loop [j (dec i) c 0]
                      (if (and (>= j 0) (> (aget mins (+ (* 3 (aget order j)) a)) k))
                        (recur (dec j) (inc c))
                        c)))))))))
