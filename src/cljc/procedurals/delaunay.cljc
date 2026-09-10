(ns procedurals.delaunay
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

(defn- edges [{[a b c] :points}]
  [#{a b} #{b c} #{c a}])

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

(defn- add-point [triangles p]
  (let [{bad true ok false} (group-by #(in-circumcircle? % p) triangles)
        boundary             (->> bad (mapcat edges) frequencies
                                  (keep (fn [[edge n]] (when (= n 1) edge))))]
    (into (vec ok) (map (fn [edge] (let [[p1 p2] (vec edge)] (make-triangle p1 p2 p)))) boundary)))

(defn triangulate
  "Delaunay triangulation of `points` (a coll of [x y]). Returns a seq of
  triangles, each `{:points [p1 p2 p3] :circle {:center [x y] :radius-sq r}}`."
  [points]
  (let [st         (super-triangle points)
        super-pts  (set (:points st))
        triangles  (reduce add-point [st] points)]
    (remove #(some super-pts (:points %)) triangles)))

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
  "Eight sites ringing `[[x0 y0] [x1 y1]]` at one box-width's remove -- the
  neighbours of the centre cell in a 3x3 tiling of the box."
  [[[x0 y0] [x1 y1]]]
  (let [w  (- x1 x0)
        h  (- y1 y0)
        xs [(- x0 w) (+ x0 (/ w 2.0)) (+ x1 w)]
        ys [(- y0 h) (+ y0 (/ h 2.0)) (+ y1 h)]]
    (for [x xs y ys
          :when (not (and (== x (second xs)) (== y (second ys))))]
      [x y])))

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
