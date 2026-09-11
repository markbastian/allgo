(ns allgo.spatial.morton
  "Morton codes: one number that says roughly where a point is.

  Interleave the bits of a point's coordinates -- x0 y0 z0 x1 y1 z1 ... --
  and you get a single integer with a useful property: points near each
  other in space usually have codes near each other, because the high bits
  of the code are the high bits of every coordinate at once. Sorting by
  that code walks the points along a Z-shaped curve that stays within a
  region before moving on.

  That is worth having whenever a spatial structure wants an order:

    a BVH     sort the leaves by code and any contiguous run is a
              compact cluster, so splitting the sorted list anywhere gives
              two children that do not overlap much. `allgo.spatial.bvh`
              is built exactly that way.
    a cache   objects visited together end up adjacent in memory.
    a GPU     the sort is a radix sort, which parallelises, and the tree
              build afterwards needs no comparisons at all.

  The curve is not perfect: two points either side of a high-bit boundary
  are far apart in code however close they are in space. That shows up as
  the occasional stretched node in a Morton BVH, and is the price of
  building a tree by sorting instead of by searching.

  Codes here are 30 bits -- 10 per axis in 3D, 15 per axis in 2D -- which
  keeps them non-negative on both platforms, ClojureScript's bit
  operations being 32-bit and signed."
  (:require [clojure.math :as math]))

(def ^:const bits-3d
  "Bits per axis in 3D. Ten of them, interleaved three ways, is 30."
  10)

(def ^:const bits-2d
  "Bits per axis in 2D. Fifteen, interleaved two ways, is 30."
  15)

(def ^:const max-3d "Largest per-axis value in 3D." (dec (bit-shift-left 1 bits-3d)))
(def ^:const max-2d "Largest per-axis value in 2D." (dec (bit-shift-left 1 bits-2d)))

;; ---------------------------------------------------------------------------
;; Bit interleaving

(defn expand-3
  "Spreads the low 10 bits of `v` out, leaving two zeros between each.

  Done by shift-and-mask rather than the multiplication trick the
  reference uses: the same thing in the same number of steps, without
  relying on a float multiply landing exactly on an integer."
  ^long [^long v]
  (let [v (bit-and v 0x3FF)
        v (bit-and (bit-or v (bit-shift-left v 16)) 0x030000FF)
        v (bit-and (bit-or v (bit-shift-left v 8)) 0x0300F00F)
        v (bit-and (bit-or v (bit-shift-left v 4)) 0x030C30C3)
        v (bit-and (bit-or v (bit-shift-left v 2)) 0x09249249)]
    v))

(defn compact-3
  "The inverse of `expand-3`: takes every third bit."
  ^long [^long v]
  (let [v (bit-and v 0x09249249)
        v (bit-and (bit-or v (bit-shift-right v 2)) 0x030C30C3)
        v (bit-and (bit-or v (bit-shift-right v 4)) 0x0300F00F)
        v (bit-and (bit-or v (bit-shift-right v 8)) 0x030000FF)
        v (bit-and (bit-or v (bit-shift-right v 16)) 0x000003FF)]
    v))

(defn expand-2
  "Spreads the low 15 bits of `v` out, leaving one zero between each."
  ^long [^long v]
  (let [v (bit-and v 0x7FFF)
        v (bit-and (bit-or v (bit-shift-left v 8)) 0x00FF00FF)
        v (bit-and (bit-or v (bit-shift-left v 4)) 0x0F0F0F0F)
        v (bit-and (bit-or v (bit-shift-left v 2)) 0x33333333)
        v (bit-and (bit-or v (bit-shift-left v 1)) 0x55555555)]
    v))

(defn compact-2
  "The inverse of `expand-2`: takes every other bit."
  ^long [^long v]
  (let [v (bit-and v 0x55555555)
        v (bit-and (bit-or v (bit-shift-right v 1)) 0x33333333)
        v (bit-and (bit-or v (bit-shift-right v 2)) 0x0F0F0F0F)
        v (bit-and (bit-or v (bit-shift-right v 4)) 0x00FF00FF)
        v (bit-and (bit-or v (bit-shift-right v 8)) 0x0000FFFF)]
    v))

(defn encode-3
  "The Morton code of integer cell `x y z`, each 0 to 1023."
  ^long [^long x ^long y ^long z]
  (bit-or (expand-3 x)
          (bit-shift-left (expand-3 y) 1)
          (bit-shift-left (expand-3 z) 2)))

(defn decode-3
  "`[x y z]` back out of a 3D code."
  [^long code]
  [(compact-3 code)
   (compact-3 (bit-shift-right code 1))
   (compact-3 (bit-shift-right code 2))])

(defn encode-2
  "The Morton code of integer cell `x y`, each 0 to 32767."
  ^long [^long x ^long y]
  (bit-or (expand-2 x) (bit-shift-left (expand-2 y) 1)))

(defn decode-2
  "`[x y]` back out of a 2D code."
  [^long code]
  [(compact-2 code) (compact-2 (bit-shift-right code 1))])

;; ---------------------------------------------------------------------------
;; Points to codes

(defn quantize
  "Puts `v` on an integer grid of `2^bits` steps spanning `lo` to `hi`.

  Clamped at both ends, so a point outside the bounds lands on the edge
  rather than wrapping to the far side of the space -- which, in a code,
  would place it as far away as it is possible to be."
  ^long [v lo hi ^long bits]
  (let [steps (bit-shift-left 1 bits)
        span  (- (double hi) (double lo))
        t     (if (pos? span) (/ (- (double v) (double lo)) span) 0.0)]
    (min (dec steps) (max 0 (long (math/floor (* t steps)))))))

(defn bounds
  "`[[min-x min-y min-z] [max-x max-y max-z]]` over `n` points, flat, 3 per
  point."
  [^doubles positions ^long n]
  (loop [i 0 lo [##Inf ##Inf ##Inf] hi [##-Inf ##-Inf ##-Inf]]
    (if (= i n)
      [lo hi]
      (let [b (* 3 i)
            x (aget positions b) y (aget positions (+ b 1)) z (aget positions (+ b 2))]
        (recur (inc i)
               [(min (lo 0) x) (min (lo 1) y) (min (lo 2) z)]
               [(max (hi 0) x) (max (hi 1) y) (max (hi 2) z)])))))

(defn code-of
  "The 3D code of a point, quantized into `[lo hi]`."
  ^long [x y z [lo hi]]
  (encode-3 (quantize x (lo 0) (hi 0) bits-3d)
            (quantize y (lo 1) (hi 1) bits-3d)
            (quantize z (lo 2) (hi 2) bits-3d)))

(defn codes
  "The code of every point, as an int array.

  `box` is the region to quantize into; by default the points' own bounds,
  which spreads them over the whole code range and is what a tree build
  wants."
  ([positions n] (codes positions n (bounds positions n)))
  ([^doubles positions ^long n box]
   (let [out #?(:clj (int-array n) :cljs (js/Int32Array. n))]
     (dotimes [i n]
       (let [b (* 3 i)]
         (aset out i (int (code-of (aget positions b)
                                   (aget positions (+ b 1))
                                   (aget positions (+ b 2))
                                   box)))))
     out)))

;; ---------------------------------------------------------------------------
;; Ordering

(defn- radix-pass!
  "One 10-bit counting-sort pass over `src`, keyed on `shift`."
  ;; No primitive hints on `n` and `shift`: a function taking any
  ;; primitive argument is capped at four arguments, and this needs six.
  [^ints codes ^ints src ^ints dst n shift ^ints counts]
  (let [n (long n) shift (long shift)
        buckets (bit-shift-left 1 10)
        mask (dec buckets)]
    (dotimes [i buckets] (aset counts i (int 0)))
    (dotimes [i n]
      (let [k (bit-and (bit-shift-right (aget codes (aget src i)) shift) mask)]
        (aset counts k (inc (aget counts k)))))
    ;; Exclusive prefix sum: where each bucket starts.
    (loop [i 0 total 0]
      (when (< i buckets)
        (let [c (aget counts i)]
          (aset counts i (int total))
          (recur (inc i) (+ total c)))))
    (dotimes [i n]
      (let [id (aget src i)
            k  (bit-and (bit-shift-right (aget codes id) shift) mask)]
        (aset dst (aget counts k) (int id))
        (aset counts k (inc (aget counts k)))))))

(defn order
  "The object indices sorted by Morton code, as an int array.

  Radix sorted: three passes of ten bits, each a counting sort, which is
  linear in the number of objects rather than n log n and needs no
  comparisons. That matters less on one core than it does on a GPU, where
  it is the reason Morton order is how large trees get built at all.

  The sort is stable, so objects sharing a code keep their original order
  and the result is the same every run."
  [^ints codes ^long n]
  (let [a #?(:clj (int-array n) :cljs (js/Int32Array. n))
        b #?(:clj (int-array n) :cljs (js/Int32Array. n))
        counts #?(:clj (int-array 1024) :cljs (js/Int32Array. 1024))]
    (dotimes [i n] (aset a i (int i)))
    ;; 30 bits, low ten first, so the last pass leaves the highest bits
    ;; dominant -- which is what makes a stable sort come out ordered.
    (radix-pass! codes a b n 0 counts)
    (radix-pass! codes b a n 10 counts)
    (radix-pass! codes a b n 20 counts)
    b))

(defn sorted-points
  "Convenience: `order` straight from positions."
  ([positions n] (order (codes positions n) n))
  ([positions n box] (order (codes positions n box) n)))
