(ns allgo.spatial.hash
  "Uniform spatial hashing: finding which of thousands of moving objects
  are near each other, rebuilt from scratch every frame.

  After Ten Minute Physics 11,
  https://matthias-research.github.io/pages/tenMinutePhysics/

  Space is cut into cubes of side `spacing` and each object is filed under
  the cube it falls in. Two things make that fast, and neither is the
  hashing:

  The cell index is *hashed into a fixed table* rather than used to index
  a grid. A grid needs memory proportional to the volume it covers and has
  to know its bounds in advance; a hash table needs memory proportional to
  the number of objects and covers infinite space. Two cells occasionally
  collide and hand back each other's objects, which costs a few extra
  distance checks and never a wrong answer -- the caller filters by real
  distance regardless, so the hash is free to be imperfect.

  The table is built by *counting sort into two flat arrays*, which is
  what makes this worth writing out rather than reaching for a map of
  vectors. Count how many objects fall in each bucket, prefix-sum those
  counts into starting offsets, then walk the objects once more writing
  each into its slot. Everything lands in one contiguous array in bucket
  order, so a query reads a slice rather than chasing pointers, and a
  rebuild allocates nothing at all. A map of vectors allocates per cell
  per frame and scatters the contents across the heap.

  Positions are flat `[x y z x y z ...]` arrays -- the layout
  `allgo.physics.xpbd` already keeps particles in, so a soft body can be
  queried without copying.

      (def h (spatial-hash 0.05 5000))
      (rebuild! h positions n)
      (let [found (query! h positions i 0.05)]
        (dotimes [k found]
          (let [j (neighbor h k)] ...)))"
  (:require [allgo.array :as a]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------

(defn spatial-hash
  "A hash sized for up to `max-objects`, on a grid of side `spacing`.

  `spacing` should be about the distance you intend to query: smaller and
  a query sweeps many cells, larger and each cell holds objects too far
  away to matter. For equal spheres in contact it is one diameter.

  The table is at least twice `max-objects`, which keeps the load factor
  near a half -- enough that most buckets hold the one cell they were
  given -- and rounded up to a power of two so that reducing a hash to a
  bucket is a mask rather than a division. That is not a micro-optimization
  here: the modulo runs once per cell visited, which for a few thousand
  objects is a few hundred thousand times a frame."
  [spacing max-objects]
  (let [table-size (loop [p 1] (if (>= p (* 2 max-objects)) p (recur (* 2 p))))]
    {:spacing      (double spacing)
     ;; Reciprocal once, so placing a coordinate is a multiply.
     :inv-spacing  (/ 1.0 (double spacing))
     :table-size   table-size
     :mask         (dec table-size)
     ;; One longer than the table: the extra entry is a guard holding the
     ;; total, so the end of the last bucket can be read as `start[h+1]`
     ;; without a special case.
     :cell-start   (a/i32 (inc table-size))
     :cell-entries (a/i32 max-objects)
     :query-ids    (a/i32 max-objects)
     :query-size   (a/i32 1)
     ;; Two cells in one query range can hash to the same bucket, and
     ;; reading that bucket twice would report its objects twice. Stamping
     ;; each object with the query that last saw it costs one comparison
     ;; and makes the result a set.
     :stamp        (a/i32 max-objects)
     :generation   (a/i32 1)}))

(defn- cell-of ^long [coord inv-spacing]
  (long (math/floor (* (double coord) (double inv-spacing)))))

(defn- hash-cell
  "A bucket for an integer cell coordinate.

  Three large odd multipliers combined with xor -- Müller calls it a
  fantasy function, and the name is fair: it has no derivation, it just
  scatters well enough. Masked rather than reduced modulo, which needs the
  table to be a power of two and makes the sign irrelevant."
  ^long [^long xi ^long yi ^long zi ^long mask]
  (bit-and (bit-xor (* xi 92837111) (* yi 689287499) (* zi 283923481)) mask))

(defn- hash-of ^long [^doubles positions ^long i inv-spacing ^long mask]
  (let [b (* 3 i)]
    (hash-cell (cell-of (aget positions b) inv-spacing)
               (cell-of (aget positions (+ b 1)) inv-spacing)
               (cell-of (aget positions (+ b 2)) inv-spacing)
               mask)))

(defn rebuild!
  "Files `n` objects from the flat `positions` array into the table.

  Three linear passes and no allocation, which is why this can be thrown
  away and rebuilt every frame rather than updated -- for objects that all
  move at once, maintaining an index incrementally costs more than
  rebuilding one."
  ([h positions] (rebuild! h positions (quot (alength ^doubles positions) 3)))
  ([{:keys [inv-spacing table-size mask ^ints cell-start ^ints cell-entries]} positions n]
   (let [^doubles positions positions
         table-size (long table-size)
         mask (long mask)
         n (min n (alength cell-entries))]
     (a/ifill! cell-start 0)

     ;; 1. How many objects fall in each bucket.
     (dotimes [i n]
       (let [b (hash-of positions i inv-spacing mask)]
         (aset cell-start b (inc (aget cell-start b)))))

     ;; 2. Turn the counts into starting offsets, running total from the
     ;;    left, with the guard entry picking up the grand total.
     (loop [i 0 start 0]
       (if (= i table-size)
         (aset cell-start table-size start)
         (let [start (+ start (aget cell-start i))]
           (aset cell-start i start)
           (recur (inc i) start))))

     ;; 3. Walk the objects again, each one stepping its bucket's offset
     ;;    back and writing itself there. Filling backward is what lets
     ;;    the offsets double as cursors: when the pass ends every bucket's
     ;;    offset has walked down to exactly where that bucket begins.
     (dotimes [i n]
       (let [b (hash-of positions i inv-spacing mask)
             slot (dec (aget cell-start b))]
         (aset cell-start b slot)
         (aset cell-entries slot i)))
     n)))

(defn query-point!
  "Every object in the cells covering the box of `max-dist` around a point.

  Returns how many were found; `neighbor` reads them out. These are
  *candidates*: everything in the cells the box touches, which includes
  objects up to a cell beyond `max-dist` and, rarely, objects from an
  unrelated cell that hashed to the same bucket. Filter by actual distance.

  Each object is reported once, however many of the cells it was reached
  through."
  [{:keys [inv-spacing mask ^ints cell-start ^ints cell-entries
           ^ints query-ids ^ints query-size ^ints stamp ^ints generation]}
   x y z max-dist]
  (let [mask (long mask)
        x0 (cell-of (- x max-dist) inv-spacing) x1 (cell-of (+ x max-dist) inv-spacing)
        y0 (cell-of (- y max-dist) inv-spacing) y1 (cell-of (+ y max-dist) inv-spacing)
        z0 (cell-of (- z max-dist) inv-spacing) z1 (cell-of (+ z max-dist) inv-spacing)
        cap (alength query-ids)
        ;; A stamp of zero means "never seen", so generations start at 1.
        ;; On the (unreachable in practice) wrap, clear the stamps rather
        ;; than let a stale one read as current.
        gen (let [g (inc (aget generation 0))]
              (if (< g 0)
                (do (a/ifill! stamp 0) 1)
                g))]
    (aset generation 0 gen)
    (aset query-size 0 0)
    (loop [xi x0]
      (when (<= xi x1)
        (loop [yi y0]
          (when (<= yi y1)
            (loop [zi z0]
              (when (<= zi z1)
                (let [b     (hash-cell xi yi zi mask)
                      start (aget cell-start b)
                      end   (aget cell-start (inc b))]
                  (loop [k start]
                    (when (< k end)
                      (let [id (aget cell-entries k)
                            found (aget query-size 0)]
                        (when (and (not= (aget stamp id) gen) (< found cap))
                          (aset stamp id gen)
                          (aset query-ids found id)
                          (aset query-size 0 (inc found))))
                      (recur (inc k)))))
                (recur (inc zi))))
            (recur (inc yi))))
        (recur (inc xi))))
    (aget query-size 0)))

(defn query!
  "`query-point!` for object `i` of a flat positions array."
  [h positions i max-dist]
  (let [^doubles positions positions
        b (* 3 i)]
    (query-point! h (aget positions b) (aget positions (+ b 1)) (aget positions (+ b 2))
                  max-dist)))

(defn neighbor
  "The `k`th candidate from the last query."
  ^long [{:keys [^ints query-ids]} k]
  (aget query-ids k))

(defn neighbors
  "The candidates from the last query, as a vector.

  Convenient and allocating; the hot path is `query!` plus `neighbor`."
  [h found]
  (into [] (map #(neighbor h %)) (range found)))

;; ---------------------------------------------------------------------------
;; Conveniences over the raw query

(defn within
  "Every object whose center is within `max-dist` of object `i`, excluding
  `i` itself. The candidates filtered down to real hits."
  [h positions i max-dist]
  (let [^doubles positions positions
        found (query! h positions i max-dist)
        b     (* 3 i)
        x (aget positions b) y (aget positions (+ b 1)) z (aget positions (+ b 2))
        r2 (* max-dist max-dist)]
    (into []
          (comp (map #(neighbor h %))
                (remove #(= % i))
                (filter (fn [j]
                          (let [c  (* 3 j)
                                dx (- (aget positions c) x)
                                dy (- (aget positions (+ c 1)) y)
                                dz (- (aget positions (+ c 2)) z)]
                            (<= (+ (* dx dx) (* dy dy) (* dz dz)) r2)))))
          (range found))))

(defn adjacency
  "Every pair within `max-dist`, as a compressed adjacency list.

  Returns `{:starts :ids :pairs}`: `ids` holds the partners of every
  object end to end, and `starts` says where each object's run begins, so
  object `i`'s partners are `ids[starts[i] .. starts[i+1])`. One flat
  array and one index, which is how a neighbor list is kept when it is
  read far more often than it is built.

  Each pair is recorded once, under the higher-numbered object, so walking
  every object's run visits every pair exactly once.

  The point of building it at all, rather than querying as you go, is that
  a solver taking many small steps between frames can detect once and
  resolve many times -- provided the radius is widened to cover how far
  anything could travel in between."
  [h positions n max-dist]
  (let [^doubles positions positions
        ^ints starts (a/i32 (inc n))
        r2     (* max-dist max-dist)]
    (rebuild! h positions n)
    (loop [i 0 acc (transient []) total 0]
      (if (= i n)
        (do (aset starts n (int total))
            {:starts starts
             :ids    (a/i32 (persistent! acc))
             :pairs  total})
        (let [_     (aset starts i (int total))
              found (query! h positions i max-dist)
              a     (* 3 i)
              [acc total]
              (loop [k 0 acc acc total total]
                (if (= k found)
                  [acc total]
                  (let [j (neighbor h k)]
                    (if (>= j i)
                      (recur (inc k) acc total)
                      (let [b  (* 3 j)
                            dx (- (aget positions a) (aget positions b))
                            dy (- (aget positions (+ a 1)) (aget positions (+ b 1)))
                            dz (- (aget positions (+ a 2)) (aget positions (+ b 2)))]
                        (if (<= (+ (* dx dx) (* dy dy) (* dz dz)) r2)
                          (recur (inc k) (conj! acc j) (inc total))
                          (recur (inc k) acc total)))))))]
          (recur (inc i) acc total))))))

(defn adjacent
  "The partners of object `i` in an `adjacency`, as a vector."
  [{:keys [^ints starts ^ints ids]} i]
  (into [] (map #(aget ids %)) (range (aget starts i) (aget starts (inc i)))))

(defn overlapping-pairs
  "Every pair of objects closer than `max-dist`, each pair once.

  The whole point of the structure: this is the O(n^2) loop everyone
  writes first, made linear in the number of objects for any distribution
  where the cells are not all piled into one."
  [h positions n max-dist]
  (let [^doubles positions positions]
    (rebuild! h positions n)
    (persistent!
     (reduce (fn [acc i]
               (reduce (fn [acc j]
                         ;; Each pair is reported from both ends; keeping
                         ;; only the ordered one halves the work and makes
                         ;; the result a set of pairs rather than of
                         ;; directed adjacencies.
                         (if (< i j) (conj! acc [i j]) acc))
                       acc
                       (within h positions i max-dist)))
             (transient [])
             (range n)))))
