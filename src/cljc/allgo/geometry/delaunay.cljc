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

(defn- add-point
  "Bowyer-Watson's step: discard every triangle whose circumcircle
  contains `p`, and re-fan the hole from `p`.

  The hole's boundary is the edges of the discarded triangles that were
  offered exactly once -- an edge offered twice was interior to the hole
  and is gone with it.

  Written as one pass with transients rather than `group-by` over
  `mapcat` over `frequencies`. Those read better and allocate a map, two
  vectors and three sets for every triangle examined, which this does
  once per point for every triangle there is. It is the same algorithm
  and the same result; it is about twice as fast, and this is the whole
  cost of a Voronoi diagram."
  [triangles p]
  (let [n (count triangles)]
    (loop [i 0 ok (transient []) counts (transient {})]
      (if (= i n)
        (let [counts (persistent! counts)]
          (reduce-kv (fn [acc [p1 p2] c]
                       (if (= 1 c) (conj acc (make-triangle p1 p2 p)) acc))
                     (persistent! ok)
                     counts))
        (let [t (nth triangles i)]
          (if (in-circumcircle? t p)
            (let [[a b c] (:points t)
                  bump (fn [m u v] (let [k (edge-key u v)]
                                     (assoc! m k (inc (long (get m k 0))))))]
              (recur (inc i) ok (-> counts (bump a b) (bump b c) (bump c a))))
            (recur (inc i) (conj! ok t) counts)))))))

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
  "Eight sites ringing `[[x0 y0] [x1 y1]]`, far enough out that none of their
  bisectors can reach inside it.

  A bisector lies half-way to its site, so a sentinel at distance r from a
  site puts a cut r/2 away. To leave the box untouched that has to exceed
  the box diagonal -- not merely its half-width, or the cut clears the edges
  but still shaves the corners. Three diagonals from the centre clears it
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
