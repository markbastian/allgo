(ns allgo.geometry.polygon
  "Plane polygons and segments: the arithmetic a city plan is made of.

  A point is `[x y]` and a polygon a vector of points, counterclockwise
  for a positive area, with the closing edge implied: edge `i` runs from
  point `i` to point `i + 1`, wrapping. Some of the functions take a
  parallel vector of *tags*, one per edge, and carry each through to the
  pieces of that edge -- which is how a lot remembers which of its sides
  face a street after the block it came from has been cut up."
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Points and segments

(defn sub [[ax ay] [bx by]] [(- (double ax) (double bx)) (- (double ay) (double by))])
(defn add [[ax ay] [bx by]] [(+ (double ax) (double bx)) (+ (double ay) (double by))])
(defn scale [[x y] k] [(* (double x) (double k)) (* (double y) (double k))])
(defn dot ^double [[ax ay] [bx by]] (+ (* (double ax) (double bx)) (* (double ay) (double by))))
(defn cross ^double [[ax ay] [bx by]] (- (* (double ax) (double by)) (* (double ay) (double bx))))
(defn norm ^double [[x y]] (math/sqrt (+ (* (double x) (double x)) (* (double y) (double y)))))
(defn distance ^double [a b] (norm (sub a b)))

(defn lerp [[ax ay] [bx by] t]
  (let [t (double t)]
    [(+ (double ax) (* t (- (double bx) (double ax))))
     (+ (double ay) (* t (- (double by) (double ay))))]))

(defn unit [v]
  (let [l (norm v)] (if (zero? l) [0.0 0.0] (scale v (/ 1.0 l)))))

(defn heading
  "The unit vector at `angle` radians counterclockwise from +x."
  [angle]
  [(math/cos angle) (math/sin angle)])

(defn angle-of ^double [[x y]] (math/atan2 y x))

(defn segment-intersection
  "Where segments `a`-`b` and `c`-`d` cross, as `[t u]`: the fractions
  along each, both in [0, 1] if they do cross, or nil if they are
  parallel."
  [a b c d]
  (let [r (sub b a) s (sub d c)
        den (cross r s)]
    (when-not (< (abs den) 1e-12)
      (let [ca (sub c a)]
        [(/ (cross ca s) den) (/ (cross ca r) den)]))))

(defn closest-on-segment
  "The fraction along `a`-`b`, clamped to [0, 1], of the point nearest
  `p`."
  ^double [p a b]
  (let [ab (sub b a) l2 (dot ab ab)]
    (if (zero? l2) 0.0 (max 0.0 (min 1.0 (/ (dot (sub p a) ab) l2))))))

(defn segment-distance
  "How far `p` is from the segment `a`-`b`."
  ^double [p a b]
  (distance p (lerp a b (closest-on-segment p a b))))

;; ---------------------------------------------------------------------------
;; Polygons

(defn signed-area
  "The shoelace formula: positive when the points run counterclockwise."
  ^double [poly]
  (let [n (count poly)]
    (* 0.5 (loop [i 0 acc 0.0]
             (if (< i n)
               (recur (inc i) (+ acc (cross (nth poly i) (nth poly (mod (inc i) n)))))
               acc)))))

(defn area ^double [poly] (abs (signed-area poly)))

(defn centroid
  "The center of mass of the polygon's area."
  [poly]
  (let [n (count poly)
        a (signed-area poly)]
    (if (< (abs a) 1e-12)
      (scale (reduce add poly) (/ 1.0 n))
      (let [[sx sy] (reduce (fn [[sx sy] i]
                              (let [[x0 y0 :as p] (nth poly i)
                                    [x1 y1 :as q] (nth poly (mod (inc i) n))
                                    c (cross p q)]
                                [(+ sx (* (+ x0 x1) c)) (+ sy (* (+ y0 y1) c))]))
                            [0.0 0.0] (range n))]
        [(/ sx (* 6.0 a)) (/ sy (* 6.0 a))]))))

(defn edges
  "The polygon's edges, as `[a b]` pairs."
  [poly]
  (let [n (count poly)]
    (mapv (fn [i] [(nth poly i) (nth poly (mod (inc i) n))]) (range n))))

(defn convex?
  "Whether every turn goes the same way: counterclockwise and convex, with
  straight runs allowed."
  [poly]
  (let [n (count poly)]
    (and (>= n 3)
         (pos? (signed-area poly))
         (every? (fn [i]
                   (let [a (nth poly i) b (nth poly (mod (inc i) n)) c (nth poly (mod (+ i 2) n))]
                     (>= (cross (sub b a) (sub c b)) (* -1e-9 (distance a b) (distance b c)))))
                 (range n)))))

(defn contains-point?
  "Whether `p` is inside the polygon, by counting crossings of a ray."
  [poly [px py]]
  (let [n (count poly) px (double px) py (double py)]
    (loop [i 0 j (dec n) inside? false]
      (if (< i n)
        (let [[xi yi] (nth poly i) [xj yj] (nth poly j)
              xi (double xi) yi (double yi) xj (double xj) yj (double yj)]
          (recur (inc i) i
                 (if (and (not= (> yi py) (> yj py))
                          (< px (+ xi (/ (* (- xj xi) (- py yi)) (- yj yi)))))
                   (not inside?)
                   inside?)))
        inside?))))

(defn clip
  "The part of the polygon on the left of the directed line through `p`
  along `d` -- Sutherland-Hodgman against one half-plane -- with each
  edge's tag carried to what is left of it and `cut-tag` on the new edge
  along the line. Returns `[points tags]`, empty when nothing is left.

  Exact for a convex polygon. A concave one can come back with a
  zero-width bridge along the line, which `convex?` will refuse."
  ([poly p d] (first (clip poly (vec (repeat (count poly) nil)) p d nil)))
  ([poly tags p d cut-tag]
   (let [n (count poly)
         side (fn [q] (cross d (sub q p)))
         out (loop [i 0 pts [] tgs []]
               (if (< i n)
                 (let [a (nth poly i) b (nth poly (mod (inc i) n))
                       sa (side a) sb (side b)
                       tag (nth tags i)
                       cut (fn [] (lerp a b (/ sa (- sa sb))))]
                   (cond
                     (and (>= sa 0.0) (>= sb 0.0)) (recur (inc i) (conj pts a) (conj tgs tag))
                     ;; Leaving: this edge is kept as far as the line,
                     ;; and the line itself runs on to where we re-enter.
                     (>= sa 0.0) (recur (inc i) (conj pts a (cut)) (conj tgs tag cut-tag))
                     ;; Entering: start from the crossing.
                     (>= sb 0.0) (recur (inc i) (conj pts (cut)) (conj tgs tag))
                     :else (recur (inc i) pts tgs)))
                 [pts tgs]))
         [pts tgs] out]
     (if (< (count pts) 3) [[] []] [pts tgs]))))

(defn- line-intersection
  "Where the line through `p` along `d` meets the one through `q` along
  `e`, or nil for parallel lines."
  [p d q e]
  (let [den (cross d e)]
    (when-not (< (abs den) 1e-12)
      (add p (scale d (/ (cross (sub q p) e) den))))))

(defn inset
  "The polygon with each edge moved inward by its own distance, from
  `distances` (one per edge, or one number for all), the corners mitered.

  Returns nil when the result has collapsed or turned inside out -- an
  edge shorter than the corners eat, or a reversed edge -- rather than a
  polygon that is not one."
  [poly distances]
  (let [n (count poly)
        ds (if (number? distances) (vec (repeat n distances)) (vec distances))
        lines (mapv (fn [i]
                      (let [a (nth poly i) b (nth poly (mod (inc i) n))
                            d (unit (sub b a))
                            inward [(- (double (second d))) (double (first d))]]
                        [(add a (scale inward (nth ds i))) d]))
                    (range n))
        pts (mapv (fn [i]
                    (let [[p d] (nth lines (mod (dec (+ i n)) n))
                          [q e] (nth lines i)]
                      ;; Straight on: no corner to miter.
                      (or (when (> (abs (cross d e)) 1e-9) (line-intersection p d q e)) q)))
                  (range n))]
    (when (and (pos? (signed-area pts))
               ;; Every edge still runs the way it did.
               (every? (fn [i]
                         (pos? (dot (sub (nth pts (mod (inc i) n)) (nth pts i))
                                    (sub (nth poly (mod (inc i) n)) (nth poly i)))))
                       (range n)))
      pts)))

(defn simplify
  "Drops points closer than `eps` to the one before and points in the
  middle of a straight run, keeping tags with the edges that survive."
  ([poly eps] (first (simplify poly (vec (repeat (count poly) nil)) eps)))
  ([poly tags eps]
   (let [eps (double eps)
         ;; Near-duplicates first: the edge from a dropped point is merged
         ;; into the one before it.
         [pts tgs] (reduce (fn [[pts tgs] [p t]]
                             (if (and (seq pts) (< (distance (peek pts) p) eps))
                               [pts tgs]
                               [(conj pts p) (conj tgs t)]))
                           [[] []] (map vector poly tags))
         [pts tgs] (if (and (> (count pts) 1) (< (distance (peek pts) (first pts)) eps))
                     [(pop pts) (pop tgs)]
                     [pts tgs])
         n (count pts)
         straight? (fn [i]
                     (let [a (nth pts (mod (dec (+ i n)) n)) b (nth pts i) c (nth pts (mod (inc i) n))]
                       (< (abs (cross (unit (sub b a)) (unit (sub c b)))) 1e-6)))
         keep (filterv #(not (straight? %)) (range n))]
     (if (or (< n 3) (< (count keep) 3))
       [pts tgs]
       ;; A point dropped from a straight run takes its edge with it; the
       ;; edge before it carries on, and keeps its own tag.
       [(mapv pts keep) (mapv tgs keep)]))))

(defn oriented-extent
  "The polygon's extent in the frame of unit vector `u` (and its left
  normal): `[min-u max-u min-v max-v]`, taken about `origin`."
  [poly origin u]
  (let [v [(- (double (second u))) (double (first u))]
        us (map #(dot (sub % origin) u) poly)
        vs (map #(dot (sub % origin) v) poly)]
    [(reduce min us) (reduce max us) (reduce min vs) (reduce max vs)]))
