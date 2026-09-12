(ns allgo.procedural.terrain
  "Fractal terrain by diamond-square.

  Start with the four corners of a square and repeatedly do two things.
  Give every square's centre the average of its four corners, which is the
  *diamond* step because the points it uses form one. Then give every
  diamond's centre the average of its four points, which is the *square*
  step for the same reason. Each round doubles the resolution, and each
  round adds a smaller random displacement than the last -- that is the
  fractal part, and the reason the result looks like landscape rather than
  like noise. Mountains have small bumps on them; the bumps are not as big
  as the mountains.

  ## Shape of the thing

  A grid is a value: a flat array of heights, how many cells on a side,
  and how large a displacement the next round may add. `step` takes one
  and returns another, so the whole algorithm is

      (nth (iterate step (init-grid opts)) iterations)

  and every intermediate is a real grid you can hold, draw, or keep.

  The array inside is the only concession to speed, and it is a large one.
  This was a map from `[i j]` vectors to heights, which is the obvious
  Clojure representation and costs an allocated key and a hash on every
  single access -- about two microseconds a cell, which is two and a half
  seconds for a grid a thousand on a side. A flat array indexed `i*dim + j`
  needs neither. `step` is still a pure function of a value to a value;
  what changed is what the value is made of, not how it is used.

  There is no need to record which cells have been filled, either, and
  that is worth seeing. After doubling the resolution, a cell's parity says
  exactly what it is: both coordinates even and it came from the previous
  round, both odd and it is a diamond centre, one of each and it is a
  square centre. The old code tested a map for absence; parity is free."
  (:require [allgo.array :as a])
  #?(:clj (:require [clojure.java.io :as io]))
  #?(:clj (:import (java.awt Color)
                   (java.awt.image BufferedImage)
                   (javax.imageio ImageIO))))

;; After https://danielbeard.wordpress.com/2010/08/07/terrain-generation-and-smoothing/

(defn height
  "The height at row `i`, column `j`."
  ^double [{:keys [^doubles heights dim]} i j]
  (aget heights (+ (* (long i) (long dim)) (long j))))

(defn init-grid
  "Seeds the 2x2 corner grid the algorithm grows from.

  Corners are `[top-left top-right bottom-right bottom-left]`. `:rng` is a
  function of no arguments returning a number in [0, 1) -- `rand` unless
  given, and the hook that makes a terrain reproducible."
  [{:keys [width corners rng] :or {corners [0.0 0.0 0.0 0.0] rng rand}}]
  (let [[c00 c01 c11 c10] corners]
    {:heights (a/f64 [(double c00) (double c01) (double c10) (double c11)])
     :dim 2
     :width (double width)
     :rng rng}))

(defn- mean-diagonal
  "Average of the four corners around a diamond centre."
  ^double [^doubles out ^long dim ^long i ^long j]
  (loop [k 0 sum 0.0 n 0.0]
    (if (= k 4)
      (if (pos? n) (/ sum n) 0.0)
      (let [a (if (or (= k 0) (= k 3)) (inc i) (dec i))
            b (if (or (= k 0) (= k 1)) (inc j) (dec j))]
        (if (and (<= 0 a) (< a dim) (<= 0 b) (< b dim))
          (recur (inc k) (+ sum (aget out (+ (* a dim) b))) (inc n))
          (recur (inc k) sum n))))))

(defn- mean-orthogonal
  "Average of the four points around a square centre."
  ^double [^doubles out ^long dim ^long i ^long j]
  (loop [k 0 sum 0.0 n 0.0]
    (if (= k 4)
      (if (pos? n) (/ sum n) 0.0)
      (let [a (case k 0 (inc i) 2 (dec i) i)
            b (case k 1 (inc j) 3 (dec j) j)]
        (if (and (<= 0 a) (< a dim) (<= 0 b) (< b dim))
          (recur (inc k) (+ sum (aget out (+ (* a dim) b))) (inc n))
          (recur (inc k) sum n))))))

(defn step
  "One round: double the resolution, then fill the two kinds of new cell.

  A pure function of a grid to a grid -- the array it returns is freshly
  its own, so intermediate grids from `iterate` stay independent of one
  another."
  [{:keys [^doubles heights dim width rng]}]
  (let [n (long dim)
        n' (dec (* 2 n))
        w (* 0.5 (double width))
        ^doubles out (a/f64 (* n' n'))
        jitter (fn ^double [] (* w (- (* 2.0 (rng)) 1.0)))]
    ;; The previous round's cells, spread out to make room.
    (dotimes [i n]
      (dotimes [j n]
        (aset out (+ (* (* 2 i) n') (* 2 j))
              (aget heights (+ (* i n) j)))))
    ;; Diamond: both coordinates odd, from the four corners around it.
    (loop [i 1]
      (when (< i n')
        (loop [j 1]
          (when (< j n')
            (aset out (+ (* i n') j)
                  (+ (mean-diagonal out n' i j) (jitter)))
            (recur (+ j 2))))
        (recur (+ i 2))))
    ;; Square: one coordinate odd, from the four points around it. Every
    ;; neighbour it wants is either an old cell or a diamond centre, so no
    ;; square centre ever waits on another.
    (loop [i 0]
      (when (< i n')
        (loop [j 0]
          (when (< j n')
            (when (not= (odd? i) (odd? j))
              (aset out (+ (* i n') j)
                    (+ (mean-orthogonal out n' i j) (jitter))))
            (recur (inc j))))
        (recur (inc i))))
    {:heights out :dim n' :width w :rng rng}))

(defn generate
  "Runs diamond-square for `iterations` rounds, giving a grid
  `2^iterations + 1` on a side."
  [{:keys [iterations] :as opts}]
  (nth (iterate step (init-grid opts)) iterations))

(defn cells->grid
  "The grid as a `dim` by `dim` vector of row vectors, indexed
  `[row][col]`."
  [{:keys [^doubles heights dim]}]
  (let [dim (long dim)]
    (mapv (fn [i] (mapv (fn [j] (aget heights (+ (* i dim) j))) (range dim)))
          (range dim))))

(defn heights-seq
  "Every height, row by row."
  [{:keys [^doubles heights]}]
  (map #(aget heights %) (range (alength heights))))

(defn cells
  "The grid as the sparse map of `[i j]` to height it used to be.

  Kept for callers that want to look a cell up by coordinate without
  knowing how the grid is stored."
  [{:keys [dim] :as grid}]
  (let [dim (long dim)]
    (persistent!
     (reduce (fn [m [i j]] (assoc! m [i j] (height grid i j)))
             (transient {})
             (for [i (range dim) j (range dim)] [i j])))))

(defn bounds
  "`[lowest highest]` over the whole grid."
  [{:keys [^doubles heights]}]
  (let [n (alength heights)]
    (loop [i 0 lo ##Inf hi ##-Inf]
      (if (= i n)
        [lo hi]
        (let [v (aget heights i)]
          (recur (inc i) (min lo v) (max hi v)))))))

#?(:clj
   (defn create-image-map
     "Writes the grid to img.png as a greyscale heightmap."
     [{:keys [dim] :as grid}]
     (let [dim (long dim)
           img (BufferedImage. dim dim BufferedImage/TYPE_INT_RGB)
           [lo hi] (bounds grid)
           span (max 1e-12 (- hi lo))]
       (dotimes [i dim]
         (dotimes [j dim]
           (let [c (float (/ (- (height grid i j) lo) span))]
             (.setRGB img i j (.getRGB (Color. c c c))))))
       (ImageIO/write img "png" (io/file "img.png")))))
