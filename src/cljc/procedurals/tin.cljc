(ns procedurals.tin
  (:require [procedurals.delaunay :as delaunay]
            [procedurals.perlin :as perlin]))

;; Triangulated Irregular Network terrain, per Paul Bourke's "An Algorithm
;; for Interpolating Irregularly-Spaced Data with Applications in Terrain
;; Modelling" (https://paulbourke.net/papers/triangulate/): Delaunay-
;; triangulate a set of scattered (x,y) sample points, then estimate the
;; height at any other (x,y) by intersecting it with the plane of whichever
;; triangular facet contains it. An alternative to `procedurals.terrain`'s
;; regular-grid diamond-square approach, useful when samples are naturally
;; irregular (e.g. survey data) rather than a lattice.

(defn scatter-points
  "`n` random 3D sample points [x y z] scattered over a `size` x `size`
  square, z assigned via multi-octave Perlin noise so nearby samples have
  correlated elevations rather than pure noise. `noise-scale` is the
  world-space wavelength of the lowest (dominant) octave -- it must be
  well under `size`, or that octave never completes a cycle across the
  domain and the whole terrain reads as one smooth tilt instead of hills."
  [n size noise-scale]
  (vec (repeatedly n (fn []
                       (let [x (rand size) y (rand size)]
                         [x y (perlin/operlin (/ x noise-scale) (/ y noise-scale) 0.0 0.5 4)])))))

(defn triangulate-points
  "Delaunay-triangulates the (x,y) projection of 3D `points` [x y z],
  returning triangles as `[p1 p2 p3]` with each vertex's original z intact."
  [points]
  (let [height-of (into {} (map (fn [[x y z]] [[x y] z])) points)
        xy        (mapv (fn [[x y _]] [x y]) points)]
    (->> (delaunay/triangulate xy)
         (map (fn [{[[ax ay] [bx by] [cx cy]] :points}]
                [[ax ay (height-of [ax ay])]
                 [bx by (height-of [bx by])]
                 [cx cy (height-of [cx cy])]])))))

(defn- sign [[x1 y1] [x2 y2] [x3 y3]]
  (- (* (- x1 x3) (- y2 y3)) (* (- x2 x3) (- y1 y3))))

(defn point-in-triangle?
  "True if 2D point `p` lies inside (or on the boundary of) triangle abc
  (each a 2D or 3D point; only x,y are used)."
  [p a b c]
  (let [d1 (sign p a b) d2 (sign p b c) d3 (sign p c a)]
    (not (and (or (neg? d1) (neg? d2) (neg? d3))
              (or (pos? d1) (pos? d2) (pos? d3))))))

(defn plane-height
  "Height of the plane through 3D triangle vertices a b c at 2D position
  [x y] -- Bourke's facet-intersection method: the query column through
  [x y] meets the triangle's plane at exactly one z."
  [[ax ay az] [bx by bz] [cx cy cz] [x y]]
  (let [ux (- bx ax) uy (- by ay) uz (- bz az)
        vx (- cx ax) vy (- cy ay) vz (- cz az)
        nx (- (* uy vz) (* uz vy))
        ny (- (* uz vx) (* ux vz))
        nz (- (* ux vy) (* uy vx))]
    (if (zero? nz)
      az
      (- az (/ (+ (* nx (- x ax)) (* ny (- y ay))) nz)))))

(defn height-at
  "Interpolated terrain height at 2D position `p`, via Bourke's method:
  finds the `triangles` facet enclosing `p` and intersects it with that
  facet's plane. nil if `p` falls outside the convex hull of the samples."
  [triangles p]
  (some (fn [[a b c]] (when (point-in-triangle? p a b c) (plane-height a b c p))) triangles))

(defn- sq [x] (* x x))

(defn- nearest-height
  "Height of whichever sample `points` is closest to `p`, for filling in
  grid cells outside the triangulated convex hull."
  [points [x y]]
  (let [[_ _ z] (apply min-key (fn [[px py _]] (+ (sq (- px x)) (sq (- py y)))) points)]
    z))

(defn generate
  "Scatters `n` random 3D sample points over a `size` x `size` domain and
  Delaunay-triangulates them into a TIN. `noise-scale` defaults to a
  quarter of `size`, giving roughly four hills across the domain."
  [{:keys [n size noise-scale] :or {size 1.0}}]
  (let [noise-scale (or noise-scale (/ size 4.0))
        points      (scatter-points n size noise-scale)]
    {:points points :triangles (triangulate-points points) :size size}))

(defn sample-grid
  "Resamples a TIN (as returned by `generate`) onto a `dim` x `dim` regular
  grid of heights, matching `procedurals.terrain/cells->grid`'s shape so
  the same renderers can consume either terrain source. Cells outside the
  samples' convex hull take the nearest sample's height."
  [{:keys [points triangles size]} dim]
  (let [step (/ size (dec dim))]
    (vec (for [i (range dim)]
           (vec (for [j (range dim)]
                  (let [p [(* j step) (* i step)]]
                    (or (height-at triangles p) (nearest-height points p)))))))))

(defn smooth-grid
  "Averages each cell with its orthogonal neighbors, `passes` times. Linear
  per-facet interpolation is only C0 continuous, so slopes jump abruptly at
  triangle edges and independently-sampled neighboring points can pinch
  into sharp cone-like peaks -- exactly the 'cone shaped peaks about local
  minima and maxima' Bourke's paper names as Method 1's weakness. A real
  elevation survey varies gently enough between neighbors to avoid this;
  our randomly-scattered, independently-heighted samples don't, so this
  rounds off the resulting kinks after the fact."
  [grid passes]
  (let [dim   (count grid)
        clamp (fn [v] (max 0 (min (dec dim) v)))
        step  (fn [g]
                (let [h (fn [ii jj] (get-in g [(clamp ii) (clamp jj)]))]
                  (vec (for [i (range dim)]
                         (vec (for [j (range dim)]
                                (/ (+ (* 4.0 (h i j)) (h (dec i) j) (h (inc i) j) (h i (dec j)) (h i (inc j))) 8.0)))))))]
    (nth (iterate step grid) passes)))

(comment
  (def tin (generate {:n 400 :size 1.0}))
  (count (:triangles tin))
  (-> tin (sample-grid 65) (smooth-grid 2)))
