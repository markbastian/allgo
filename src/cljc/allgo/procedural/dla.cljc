(ns allgo.procedural.dla
  "Diffusion-limited aggregation, and mountain ridges made out of it.

  Release a particle far from a seed, let it wander at random, and stick
  it where it first touches what is already there. Repeat. The shape that
  grows is dendritic -- a trunk with branches, branches with twigs -- and
  it is not put there by any rule that mentions branching. It happens
  because the cluster shadows itself: an arriving walker is far more
  likely to meet a tip that sticks out than to find its way down into a
  gap between two, so tips grow faster than gaps fill, and every tip is
  soon a branch with tips of its own.

  Witten and Sander's model of soot and electrodeposition, in other
  words, borrowed here because that is also what a mountain range looks
  like from above. Ridges branch, and the branching is not self-similar
  the way fractal noise is -- it has a *direction*, out from the trunk to
  the tips, which is what `allgo.procedural.fractal` cannot produce at
  any setting.

  ## The cluster is a tree, and that is what makes it terrain

  Every particle remembers what it stuck to, so the cluster is not a set
  of cells but a tree rooted at the seed. That is the whole reason this
  is worth having for terrain rather than for textures: the tree can be
  measured. A node carrying half the cluster behind it is a trunk and
  should be high; a node with nothing behind it is a twig and should be
  low. `subtree-sizes` counts that, and `heightmap` reads a ridge line
  off it -- high in the middle, tapering to every tip, with the height
  ordering following the branching structure rather than the geometry.

  ## Getting it to finish

  Walkers released from anywhere and wandering until they arrive is the
  definition and not an algorithm: a walk on a plane is recurrent but
  takes forever, and most of the plane is nowhere near the cluster. Two
  standard bounds fix it, and they are what the whole namespace's speed
  rests on:

    birth   a walker starts on a circle just outside the cluster's
            current reach, not at the edge of the grid. Everything
            further out is empty, and crossing it is a random walk with
            nothing to hit.
    kill    a walker that wanders well past that circle is abandoned and
            a new one released. It has not failed -- by symmetry a walker
            that far out is as likely to come back from anywhere else, so
            restarting it loses nothing and stops the tail.

  With both, the cost is roughly the cluster's own size rather than the
  grid's area."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

(def defaults
  {:dim 129
   ;; As a fraction of the grid's cells. Past a few percent the cluster
   ;; fills the disk it lives in and the branches stop reading as
   ;; branches.
   :density 0.06
   :rng rand
   ;; How far outside the cluster a walker starts, and how far past that
   ;; it is allowed to stray, both in cells.
   :birth-margin 4
   :kill-margin 24})

;; ---------------------------------------------------------------------------
;; The aggregate
;;
;; Cells are flat indices, `y * dim + x`, and occupancy and parentage
;; live in arrays rather than in a set and a map of `[x y]` vectors. That
;; is not premature: a walker asks "is this cell taken?" several times
;; per step and takes thousands of steps, so the question is asked
;; millions of times for a cluster of a few thousand. Hashing a
;; two-element vector to answer it cost six times the whole rest of the
;; algorithm -- 3.6 seconds for a 257-square against 600ms for a
;; 129-square, where the arrays do both in a fraction of that.

(defn xy
  "The `[x y]` of a flat cell index."
  [{:keys [dim]} cell]
  [(rem (long cell) (long dim)) (quot (long cell) (long dim))])

(defn aggregate
  "Grows a cluster from a seed in the middle of a `:dim` square grid.

  Returns `{:dim :seed :parent :order}`. Cells are flat indices, `:parent`
  is an int array giving the cell each one stuck to (-1 for the seed), and
  `:order` is the cells in the order they arrived -- which is also an
  order in which every cell follows its parent."
  [opts]
  (let [{:keys [dim density rng birth-margin kill-margin]} (merge defaults opts)
        dim (long dim)
        center (quot dim 2)
        seed (+ (* center dim) center)
        target (max 1 (long (* (double density) dim dim)))
        half (double (dec (quot dim 2)))
        ;; Ints rather than booleans: ClojureScript has no boolean array,
        ;; and `allgo.array` exists so that a flat array means the same
        ;; thing on both platforms.
        occ (a/i32 (* dim dim))
        parent (a/i32 (* dim dim) -1)
        taken? (fn [^long x ^long y]
                 (and (<= 0 x) (< x dim) (<= 0 y) (< y dim)
                      (pos? (aget ^ints occ (+ (* y dim) x)))))]
    (aset ^ints occ seed 1)
    (loop [order (transient [seed])
           grown 1
           reach 1.0
           misses 0]
      (if (or (>= grown target) (> misses 60))
        ;; Either the cluster is the size asked for, or sixty walkers in a
        ;; row failed to arrive and it is as grown as it is going to get.
        {:dim dim :seed seed :parent parent :order (persistent! order)}
        (let [birth (min half (+ reach (double birth-margin)))
              kill (+ birth (double kill-margin))
              ang (* 2.0 math/PI (rng))
              sx (+ center (long (math/round (* birth (math/cos ang)))))
              sy (+ center (long (math/round (* birth (math/sin ang)))))
              ;; Walk until it touches something, strays past the kill
              ;; circle, or has plainly got lost.
              hit (loop [x sx y sy steps 0]
                    (cond
                      (or (neg? x) (neg? y) (>= x dim) (>= y dim)) nil
                      (> steps 20000) nil
                      (> (math/sqrt (+ (* (- x center) (- x center))
                                       (* (- y center) (- y center))))
                         kill)
                      nil

                      (taken? (inc x) y) [x y (+ (* y dim) (inc x))]
                      (taken? (dec x) y) [x y (+ (* y dim) (dec x))]
                      (taken? x (inc y)) [x y (+ (* (inc y) dim) x)]
                      (taken? x (dec y)) [x y (+ (* (dec y) dim) x)]

                      :else
                      (let [d (long (* (rng) 4.0))]
                        (recur (case d 0 (inc x) 1 (dec x) x)
                               (case d 2 (inc y) 3 (dec y) y)
                               (inc steps)))))]
          (if (nil? hit)
            (recur order grown reach (inc misses))
            (let [[hx hy stuck-to] hit
                  cell (+ (* hy dim) hx)]
              (aset ^ints occ cell 1)
              (aset ^ints parent cell (int stuck-to))
              (recur (conj! order cell)
                     (inc grown)
                     (max reach (math/sqrt (+ (* (- hx center) (- hx center))
                                              (* (- hy center) (- hy center)))))
                     0))))))))

(defn subtree-sizes
  "How many cells hang off each cell, itself included, as an int array.

  The trunk's measure. `:order` has every cell after its parent, so one
  pass backward over it accumulates each node into its parent and every
  child is finished before the parent is read."
  ^ints [{:keys [dim order parent]}]
  (let [^ints parent parent
        sizes (a/i32 (* (long dim) (long dim)))
        n (count order)]
    (dotimes [i n] (aset ^ints sizes (long (nth order i)) 1))
    (loop [i (dec n)]
      (when (>= i 0)
        (let [cell (long (nth order i))
              p (aget ^ints parent cell)]
          (when (>= p 0)
            (aset ^ints sizes p (+ (aget ^ints sizes p) (aget ^ints sizes cell)))))
        (recur (dec i))))
    sizes))

;; ---------------------------------------------------------------------------
;; Turning it into ground

(defn- blur!
  "One pass of a 3x3 average over a heightmap, in place."
  [^doubles h ^long dim ^long passes]
  (dotimes [_ passes]
    (let [^doubles src (a/f64 (* dim dim))]
      (dotimes [i (* dim dim)] (aset src i (aget h i)))
      (dotimes [y dim]
        (dotimes [x dim]
          (let [sum (loop [k 0 s 0.0 n 0]
                      (if (= k 9)
                        (/ s n)
                        (let [ax (+ x (- (rem k 3) 1))
                              ay (+ y (- (quot k 3) 1))]
                          (if (and (<= 0 ax) (< ax dim) (<= 0 ay) (< ay dim))
                            (recur (inc k) (+ s (aget src (+ (* ay dim) ax))) (inc n))
                            (recur (inc k) s n)))))]
            (aset h (+ (* y dim) x) (double sum)))))))
  h)

(defn- upsample
  "Doubles a heightmap's resolution, bilinearly."
  [^doubles src ^long dim]
  (let [dim' (dec (* 2 dim))
        ^doubles out (a/f64 (* dim' dim'))]
    (dotimes [y dim']
      (dotimes [x dim']
        (let [fx (/ (double x) 2.0) fy (/ (double y) 2.0)
              x0 (min (- dim 2) (long fx)) y0 (min (- dim 2) (long fy))
              tx (- fx x0) ty (- fy y0)
              at (fn [i j] (aget src (+ (* j dim) i)))]
          (aset out (+ (* y dim') x)
                (+ (* (- 1.0 tx) (- 1.0 ty) (at x0 y0))
                   (* tx (- 1.0 ty) (at (inc x0) y0))
                   (* (- 1.0 tx) ty (at x0 (inc y0)))
                   (* tx ty (at (inc x0) (inc y0))))))))
    [out dim']))

(defn heightmap
  "A heightmap with the cluster's ridges standing up out of it.

  Height comes from the tree rather than from the picture: a cell is as
  high as the log of what hangs off it, so the trunk is a ridge line, a
  major branch is a spur, and a twig is barely a rise. Logarithmic
  because subtree sizes are wildly skewed -- the seed carries everything
  and a tip carries one -- and linearly that is a spike in the middle of
  a plain.

  Then blurred, and optionally doubled and blurred again, which is what
  turns a one-cell-wide skeleton into something with sides. The
  refinement is the point rather than a finish: blurring at the coarse
  size and then interpolating gives smooth flanks that still remember
  where the ridge was, which is not what blurring at the final size
  would give.

  Returns `{:heights :dim}`, the shape `allgo.procedural.erosion` and
  `allgo.demo.terrain-webgl` already take."
  ([cluster] (heightmap cluster {}))
  ([{:keys [dim order] :as cluster} {:keys [blur refine] :or {blur 2 refine 1}}]
   (let [dim (long dim)
         ^ints sizes (subtree-sizes cluster)
         top (math/log (+ 1.0 (double (reduce max 1 (map #(aget sizes (long %)) order)))))
         ^doubles h (a/f64 (* dim dim))]
     (doseq [cell order]
       (aset h (long cell)
             (/ (math/log (+ 1.0 (double (aget sizes (long cell))))) (double top))))
     (blur! h dim blur)
     (loop [^doubles hs h d dim n (long refine)]
       (if (zero? n)
         {:heights hs :dim d}
         (let [[up d'] (upsample hs d)]
           (recur (blur! up d' blur) d' (dec n))))))))
