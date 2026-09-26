(ns allgo.geometry.delaunay
  (:require #?(:clj [clojure.math :as math] :cljs [cljs.math :as math])))

;; Bowyer-Watson incremental Delaunay triangulation, plus the Voronoi
;; diagram as its dual (Voronoi vertices are Delaunay triangle circumcenters;
;; a site's Voronoi cell is the fan of circumcenters of triangles incident
;; to it). https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm

(defn- sq [x] (* x x))

(defn- distance-sq [[x1 y1] [x2 y2]] (+ (sq (- x2 x1)) (sq (- y2 y1))))

(defn circumcircle
  "Center and squared radius of the circle through three points. nil if the
  points are collinear (no finite circumcircle)."
  [[ax ay] [bx by] [cx cy]]
  (let [d (* 2.0 (+ (* ax (- by cy)) (* bx (- cy ay)) (* cx (- ay by))))]
    (when-not (zero? d)
      (let [ax2ay2 (+ (sq ax) (sq ay))
            bx2by2 (+ (sq bx) (sq by))
            cx2cy2 (+ (sq cx) (sq cy))
            ux     (/ (+ (* ax2ay2 (- by cy)) (* bx2by2 (- cy ay)) (* cx2cy2 (- ay by))) d)
            uy     (/ (+ (* ax2ay2 (- cx bx)) (* bx2by2 (- ax cx)) (* cx2cy2 (- bx ax))) d)]
        {:center [ux uy] :radius-sq (distance-sq [ux uy] [ax ay])}))))

(defn- make-triangle [p1 p2 p3]
  {:points [p1 p2 p3] :circle (circumcircle p1 p2 p3)})

(defn- in-circumcircle? [{:keys [circle]} p]
  (and circle (<= (distance-sq (:center circle) p) (:radius-sq circle))))

(defn- super-triangle
  "A triangle guaranteed to strictly contain every point (and their
  circumcircles), for the algorithm to start from and discard at the end."
  [points]
  (let [xs        (map first points)
        ys        (map second points)
        mid-x     (/ (+ (apply min xs) (apply max xs)) 2.0)
        mid-y     (/ (+ (apply min ys) (apply max ys)) 2.0)
        delta-max (* 20 (max 1.0 (- (apply max xs) (apply min xs)) (- (apply max ys) (apply min ys))))]
    (make-triangle [(- mid-x delta-max) (- mid-y delta-max)]
                   [mid-x (+ mid-y delta-max)]
                   [(+ mid-x delta-max) (- mid-y delta-max)])))

(defn- edge-key
  "A Delaunay edge, in an order that does not depend on which triangle
  offered it."
  [a b]
  (if (neg? (compare a b)) [a b] [b a]))

(defn- tri-edges [{[a b c] :points}]
  [(edge-key a b) (edge-key b c) (edge-key c a)])

;; ---------------------------------------------------------------------------
;; The index
;;
;; Bowyer-Watson asks one question per inserted point: which triangles
;; have `p` inside their circumcircle? Asked of every triangle there is,
;; that is what makes the plain algorithm quadratic -- and it is nearly
;; all wasted, because the answer is a handful of triangles around `p`
;; and the rest are nowhere near it.
;;
;; The way out is a necessary condition that is cheap to index on: a
;; circumcircle containing `p` has a bounding box containing `p`. So keep
;; a uniform grid of triangles by circumcircle bounding box and look only
;; in the cell `p` falls in. Nothing is approximated -- a triangle that
;; would have been found by the scan is registered in that cell, because
;; its box covers the point -- so the same triangles come out.
;;
;; This wants no adjacency between triangles, no point location and no
;; argument about the cavity being connected, which is what the textbook
;; fix needs and where its bugs live.
;;
;; The one awkward case is a circumcircle far larger than a grid cell.
;; Triangles touching the super-triangle have them, and registering one
;; in every cell it covers would cost more than the scan it replaces. So
;; a triangle covering more than `max-index-cells` is held aside in a
;; list that is searched every time. That list stays short -- those are
;; the triangles on the hull, and there are far fewer of them than there
;; are points.

(def ^:private max-index-cells
  "How many cells a triangle may occupy before it is held aside instead.

  Both ends of this are expensive and the middle was found by measuring.
  Set it low and ordinary triangles fall into the held-aside list, which
  is scanned on every remaining insertion -- at twelve, a fifth of all
  live triangles ended up there and the quadratic term came straight
  back. Set it high and the enormous circumcircles around the
  super-triangle get filed into thousands of cells apiece, and the
  bookkeeping costs more than the search it saves.

  Measured on uniform random points, milliseconds to triangulate:

      threshold    n=1000   n=2000   n=4000   n=8000
             12        54      165      426     1306
             32        46      103      236      582
             64        51      114      240      563
            256        62      136      294      627
        no limit       86      203      492     1125

  Flat between about 32 and 128, so 48 rather than either edge."
  48)

(defn- make-grid
  "A uniform grid over the points, about one cell per point.

  One per point is the useful density: many more and the bookkeeping
  outweighs the search, many fewer and each cell holds a crowd."
  [points]
  (let [xs (map first points)
        ys (map second points)
        x0 (double (apply min xs))
        x1 (double (apply max xs))
        y0 (double (apply min ys))
        y1 (double (apply max ys))
        side (max 1 (long (math/ceil (math/sqrt (double (count points))))))]
    {:ox x0 :oy y0
     :cw (/ (max 1e-9 (- x1 x0)) side)
     :ch (/ (max 1e-9 (- y1 y0)) side)
     :n side}))

(defn- cell-index [{:keys [ox oy cw ch n]} [x y]]
  (let [n (long n)
        i (min (dec n) (max 0 (long (math/floor (/ (- (double x) ox) cw)))))
        j (min (dec n) (max 0 (long (math/floor (/ (- (double y) oy) ch)))))]
    (+ (* j n) i)))

(defn- cells-covering
  "Flat indices of the cells a circumcircle reaches, clipped to the grid,
  or nil when it reaches none.

  Clipping is safe because the only points ever looked up are the input
  points, and those are inside the grid by construction: a cell outside
  it can never be queried, so leaving a triangle out of one loses
  nothing."
  [{:keys [ox oy cw ch n]} {:keys [center radius-sq]}]
  (let [[cx cy] center
        r (math/sqrt (double radius-sq))
        n (long n)
        last-cell (dec n)
        i0 (max 0 (long (math/floor (/ (- (- (double cx) r) ox) cw))))
        i1 (min last-cell (long (math/floor (/ (- (+ (double cx) r) ox) cw))))
        j0 (max 0 (long (math/floor (/ (- (- (double cy) r) oy) ch))))
        j1 (min last-cell (long (math/floor (/ (- (+ (double cy) r) oy) ch))))]
    (when (and (<= i0 i1) (<= j0 j1))
      (for [j (range j0 (inc j1)) i (range i0 (inc i1))] (+ (* j n) i)))))

;; The grid is a flat array of vectors of triangle ids, and ids are never
;; taken back out of it. A dead triangle is left where it is and skipped
;; when the cell is read, which is worth the stale entries: removal means
;; touching every cell the triangle covered, and the first version of
;; this did it with nested persistent maps. On the JVM that was still a
;; large win; in a browser, where those operations cost several times
;; more, it made a five hundred point triangulation *slower* than the
;; linear scan it replaced, and only paid off past a few thousand.
;;
;; Leaving the dead behind is safe as long as they do not pile up, and
;; they do not: a triangulation creates O(n) triangles over its life and
;; each is filed in a handful of cells, so a cell sees a couple of dozen
;; ids in total however large the input is.
;;
;; `big` is the exception and is pruned properly. It is read in full on
;; every insertion, so anything dead in it is paid for again and again.

(defn- place
  "Files triangle `id` in the grid, or returns it for the held-aside list."
  [^objects grid spec id {:keys [circle]}]
  (let [cells (when circle (cells-covering spec circle))]
    (if (or (nil? circle) (> (count cells) max-index-cells))
      ;; A triangle with no circumcircle is three collinear points; it can
      ;; never contain anything, but it still has to be somewhere.
      false
      (do (doseq [c cells]
            (aset grid (long c) (conj (or (aget grid (long c)) []) id)))
          true))))

(defn- add-point
  "Bowyer-Watson's step: discard every triangle whose circumcircle
  contains `p`, and re-fan the hole from `p`.

  The hole's boundary is the edges of the discarded triangles that were
  offered exactly once -- an edge offered twice was interior to the hole
  and is gone with it."
  [{:keys [tris big spec] :as state} p]
  (let [^objects grid (:grid state)
        live? (fn [id] (some? (nth tris id)))
        bad (persistent!
             (reduce (fn [acc id]
                       (if (and (live? id) (in-circumcircle? (nth tris id) p))
                         (conj! acc id)
                         acc))
                     (transient [])
                     (concat (aget grid (cell-index spec p)) big)))]
    (if (empty? bad)
      state
      (let [counts (reduce (fn [m id]
                             (reduce (fn [m e] (assoc m e (inc (long (get m e 0)))))
                                     m
                                     (tri-edges (nth tris id))))
                           {}
                           bad)
            tris (reduce #(assoc! %1 %2 nil) tris bad)
            big (reduce disj big bad)]
        (reduce-kv (fn [st [p1 p2] c]
                     (if (= 1 c)
                       (let [t (make-triangle p1 p2 p)
                             id (count (:tris st))
                             st (update st :tris conj! t)]
                         (if (place grid spec id t)
                           st
                           (update st :big conj id)))
                       st))
                   (assoc state :tris tris :big big)
                   counts)))))

(defn triangulate
  "Delaunay triangulation of `points` (a coll of [x y]). Returns a seq of
  triangles, each `{:points [p1 p2 p3] :circle {:center [x y] :radius-sq r}}`."
  [points]
  (let [pts (vec points)]
    (if (empty? pts)
      ()
      (let [st (super-triangle pts)
            super-pts (set (:points st))
            spec (make-grid pts)
            grid (object-array (* (long (:n spec)) (long (:n spec))))
            indexed? (place grid spec 0 st)
            state {:tris (transient [st])
                   :grid grid
                   :spec spec
                   :big (if indexed? #{} #{0})}]
        (->> (reduce add-point state pts)
             :tris
             persistent!
             (remove nil?)
             (remove #(some super-pts (:points %))))))))

(defn- angle-around [[ox oy] [x y]] (math/atan2 (- y oy) (- x ox)))

(defn voronoi-cells
  "The Voronoi diagram dual to Delaunay `triangles` (as returned by
  `triangulate`): a map from each site point to its cell, a seq of vertices
  (triangle circumcenters) ordered counter-clockwise around the site. Cells
  for sites on the convex hull of the input are open -- their true Voronoi
  region is unbounded -- so this gives their finite triangle fan rather
  than a closed polygon."
  [triangles]
  (->> triangles
       (mapcat (fn [{:keys [points circle]}] (map (fn [p] [p (:center circle)]) points)))
       (group-by first)
       (reduce-kv (fn [cells site vs] (assoc cells site (sort-by (partial angle-around site) (map second vs))))
                  {})))

(defn voronoi-diagram
  "Voronoi diagram of `points`; see `voronoi-cells`."
  [points]
  (voronoi-cells (triangulate points)))

(defn- sentinel-ring
  "Eight sites ringing `[[x0 y0] [x1 y1]]`, far enough out that none of their
  bisectors can reach inside it.

  A bisector lies half-way to its site, so a sentinel at distance r from a
  site puts a cut r/2 away. To leave the box untouched that has to exceed
  the box diagonal -- not merely its half-width, or the cut clears the edges
  but still shaves the corners. Three diagonals from the center clears it
  for any site in the box."
  [[[x0 y0] [x1 y1]]]
  (let [cx (/ (+ x0 x1) 2.0)
        cy (/ (+ y0 y1) 2.0)
        w  (- x1 x0)
        h  (- y1 y0)
        r  (* 3.0 (math/sqrt (+ (* w w) (* h h))))]
    (for [k (range 8)
          :let [a (* k (/ math/PI 4.0))]]
      [(+ cx (* r (math/cos a))) (+ cy (* r (math/sin a)))])))

(defn- intersect-x [[ax ay] [bx by] k]
  (let [t (/ (- k ax) (- bx ax))] [k (+ ay (* t (- by ay)))]))

(defn- intersect-y [[ax ay] [bx by] k]
  (let [t (/ (- k ay) (- by ay))] [(+ ax (* t (- bx ax))) k]))

(defn- clip-halfplane [ring inside? cut]
  (->> (partition 2 1 (conj (vec ring) (first ring)))
       (into [] (mapcat (fn [[a b]]
                          (case [(boolean (inside? a)) (boolean (inside? b))]
                            [true true]   [b]
                            [true false]  [(cut a b)]
                            [false true]  [(cut a b) b]
                            []))))))

(defn clip-polygon
  "Sutherland-Hodgman clip of convex `poly` to the axis-aligned rectangle
  `[[x0 y0] [x1 y1]]`, returning the clipped ring (empty when the polygon
  lies wholly outside). Voronoi cells are always convex, so the result is
  exact for them."
  [poly [[x0 y0] [x1 y1]]]
  (reduce (fn [ring [inside? cut]]
            (if (seq ring) (clip-halfplane ring inside? cut) ring))
          (vec poly)
          [[#(>= (first %) x0)  #(intersect-x %1 %2 x0)]
           [#(<= (first %) x1)  #(intersect-x %1 %2 x1)]
           [#(>= (second %) y0) #(intersect-y %1 %2 y0)]
           [#(<= (second %) y1) #(intersect-y %1 %2 y1)]]))

(defn bounded-cells
  "Voronoi cells for `points`, each a closed polygon clipped to `bounds`
  (`[[x0 y0] [x1 y1]]`). Together they tile the rectangle exactly.

  `voronoi-cells` leaves hull sites open, since their true regions really are
  unbounded. Ringing the input with sentinel sites outside `bounds` makes
  every real site interior, so each gets a finite cell; the sentinels are
  dropped and what remains is clipped back to `bounds`."
  [points bounds]
  (let [ring (sentinel-ring bounds)]
    (reduce-kv (fn [cells site cell]
                 (let [clipped (clip-polygon cell bounds)]
                   (cond-> cells (>= (count clipped) 3) (assoc site clipped))))
               {}
               (apply dissoc (voronoi-diagram (into (vec points) ring)) ring))))

(comment
  (def pts (repeatedly 12 #(vector (rand-int 100) (rand-int 100))))
  (def tris (triangulate pts))
  (count tris)
  (voronoi-cells tris)
  (voronoi-diagram pts))
