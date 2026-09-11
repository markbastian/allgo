(ns allgo.geometry.hex
  "Hexagonal grids, following Amit Patel's guide and its companion
  implementation notes:

    https://www.redblobgames.com/grids/hexagons/
    https://www.redblobgames.com/grids/hexagons/implementation.html

  A hex is an axial pair `[q r]` -- a plain vector, so it works as a map
  key and prints readably. The third cube coordinate is always `s = -q-r`,
  so it carries no information and is derived on demand by `cube`;
  functions that are clearer in cube coordinates convert internally and
  hand back axial. Anywhere a hex is taken, a 3-element cube vector is
  accepted too.

  The payoff for cube coordinates is that a hex grid becomes a plane
  through a cubic lattice, and the usual vector operations are suddenly
  meaningful: distance is the Chebyshev distance, rotation is a cyclic
  shift of the three coordinates, and reflection is a swap of two of them.

  `layout` connects hexes to the screen: it holds the orientation, the
  size of a hex and where the origin sits, and `->pixel` / `pixel->hex`
  move between the two. Nothing above that layer knows about pixels."
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Coordinates

(defn cube
  "`[q r s]` for a hex, with `s = -q-r`. Idempotent on cube input."
  [h]
  (if (= 3 (count h))
    (vec h)
    (let [[q r] h] [q r (- (+ q r))])))

(defn axial
  "`[q r]` for a hex. Idempotent on axial input."
  [h]
  (subvec (vec h) 0 2))

(defn cube?
  "Whether `[q r s]` satisfies the constraint that makes it a hex."
  [h]
  (and (= 3 (count h)) (zero? (reduce + h))))

;; ---------------------------------------------------------------------------
;; Arithmetic

(defn add [a b] (mapv + (axial a) (axial b)))
(defn subtract [a b] (mapv - (axial a) (axial b)))
(defn scale [h k] (mapv #(* % k) (axial h)))

(defn hex-length
  "Steps from the origin to `h`.

  Half the sum of the absolute cube coordinates: any move changes exactly
  two of the three by one, so the sum counts every step twice."
  [h]
  (long (/ (reduce + (map abs (cube h))) 2)))

(defn distance
  "Steps between two hexes."
  [a b]
  (hex-length (subtract a b)))

;; ---------------------------------------------------------------------------
;; Directions

(def directions
  "The six neighbouring directions, as axial offsets. Index 0 is `[1 0]`
  and they proceed counter-clockwise in a pointy-top layout."
  [[1 0] [1 -1] [0 -1] [-1 0] [-1 1] [0 1]])

(def diagonals
  "The six hexes two steps away that are *not* along an axis -- the ones
  sharing only a corner. Not the directions doubled: `[2 0]` is two steps
  along one axis, whereas `[2 -1]` is the diagonal between two of them."
  [[2 -1] [1 -2] [-1 -1] [-2 1] [-1 2] [1 1]])

(defn direction [i] (directions (mod i 6)))
(defn diagonal [i] (diagonals (mod i 6)))

(defn neighbor [h i] (add h (direction i)))
(defn neighbors [h] (mapv #(add h %) directions))

(defn diagonal-neighbor [h i] (add h (diagonal i)))
(defn diagonal-neighbors [h] (mapv #(add h %) diagonals))

(defn walk
  "An infinite seq of hexes stepping from `h` in direction `i`."
  [h i]
  (iterate #(neighbor % i) h))

;; ---------------------------------------------------------------------------
;; Rotation and reflection

(defn rotate-left
  "60 degrees counter-clockwise about the origin, or about `center`.

  In cube coordinates a rotation is just a cyclic shift with a sign flip,
  which is the clearest dividend of the three-coordinate view."
  ([h] (let [[q r s] (cube h)] (axial [(- s) (- q) (- r)])))
  ([h center] (add center (rotate-left (subtract h center)))))

(defn rotate-right
  "60 degrees clockwise about the origin, or about `center`."
  ([h] (let [[q r s] (cube h)] (axial [(- r) (- s) (- q)])))
  ([h center] (add center (rotate-right (subtract h center)))))

(defn reflect-q
  "Mirror across the q axis, about the origin or about `center`. The axis
  coordinate is kept and the other two swap."
  ([h] (let [[q r s] (cube h)] (axial [q s r])))
  ([h center] (add center (reflect-q (subtract h center)))))

(defn reflect-r
  "Mirror across the r axis."
  ([h] (let [[q r s] (cube h)] (axial [s r q])))
  ([h center] (add center (reflect-r (subtract h center)))))

(defn reflect-s
  "Mirror across the s axis."
  ([h] (let [[q r s] (cube h)] (axial [r q s])))
  ([h center] (add center (reflect-s (subtract h center)))))

;; ---------------------------------------------------------------------------
;; Rounding and lines

(defn round
  "The hex containing a fractional cube or axial coordinate.

  Rounding each coordinate separately can break `q + r + s = 0`, so the
  one that moved furthest is recomputed from the other two -- discarding
  the least trustworthy of the three rather than an arbitrary one."
  [h]
  (let [[q r s]    (cube h)
        [rq rr rs] (map #(math/round (double %)) [q r s])
        [dq dr ds] (map #(abs (- %1 %2)) [rq rr rs] [q r s])]
    (cond
      (and (> dq dr) (> dq ds)) (axial [(- (- rr) rs) rr])
      (> dr ds)                 [rq (- (- rq) rs)]
      :else                     [rq rr])))

(defn lerp
  "Fractional cube coordinate `t` of the way from `a` to `b`."
  [a b t]
  (mapv (fn [x y] (+ x (* (- y x) t))) (cube a) (cube b)))

(defn line
  "Every hex a straight line from `a` to `b` passes through, inclusive.

  The endpoints are nudged off the lattice first. A line that runs exactly
  along the boundary between two hexes is genuinely ambiguous, and without
  the nudge the rounding picks one side or the other depending on
  floating-point noise, which makes the line wander. The nudge commits to
  one side consistently."
  [a b]
  (let [n    (distance a b)
        [aq ar as] (cube a)
        [bq br bs] (cube b)
        a*   [(+ aq 1e-6) (+ ar 1e-6) (- as 2e-6)]
        b*   [(+ bq 1e-6) (+ br 1e-6) (- bs 2e-6)]
        step (/ 1.0 (max n 1))]
    (mapv #(round (lerp a* b* (* step %))) (range (inc n)))))

;; ---------------------------------------------------------------------------
;; Regions

(defn hexes-within
  "Every hex within `n` steps of `center`, itself included. 1 + 3n(n+1) of
  them."
  [center n]
  (vec (for [q (range (- n) (inc n))
             r (range (max (- n) (- (- q) n)) (inc (min n (+ (- q) n))))]
         (add center [q r]))))

(defn intersecting-ranges
  "The hexes within `n1` of `c1` and also within `n2` of `c2`.

  Each range is a box in cube coordinates, so the intersection is found by
  intersecting the three coordinate intervals rather than by generating
  both sets and comparing them."
  [c1 n1 c2 n2]
  (let [[q1 r1 s1] (cube c1)
        [q2 r2 s2] (cube c2)
        qmin (max (- q1 n1) (- q2 n2)) qmax (min (+ q1 n1) (+ q2 n2))
        rmin (max (- r1 n1) (- r2 n2)) rmax (min (+ r1 n1) (+ r2 n2))
        smin (max (- s1 n1) (- s2 n2)) smax (min (+ s1 n1) (+ s2 n2))]
    (vec (for [q (range qmin (inc qmax))
               r (range (max rmin (- (- q) smax)) (inc (min rmax (- (- q) smin))))]
           [q r]))))

(defn ring
  "The `n` hexes-per-side ring exactly `radius` steps from `center`, in
  order around it. Radius 0 is the centre alone."
  [center radius]
  (if (zero? radius)
    [center]
    ;; Start on the corner in direction 4 so that stepping in direction i
    ;; walks the side rather than immediately leaving the ring.
    (loop [hex (add center (scale (direction 4) radius))
           i   0
           out []]
      (if (= i 6)
        out
        (let [[hex side] (reduce (fn [[hex acc] _]
                                   [(neighbor hex i) (conj acc hex)])
                                 [hex []]
                                 (range radius))]
          (recur hex (inc i) (into out side)))))))

(defn spiral
  "Every hex within `radius`, ordered from the centre outward ring by
  ring."
  [center radius]
  (into [center] (mapcat #(ring center %)) (range 1 (inc radius))))

;; ---------------------------------------------------------------------------
;; Map shapes

(defn parallelogram [q1 q2 r1 r2]
  (vec (for [q (range q1 (inc q2)) r (range r1 (inc r2))] [q r])))

(defn triangle [size]
  (vec (for [q (range 0 (inc size)) r (range 0 (inc (- size q)))] [q r])))

(defn hexagon [radius] (hexes-within [0 0] radius))

(defn rectangle
  "A rectangular block of `width` by `height` hexes. The rows (pointy) or
  columns (flat) are sheared back as they advance so the block comes out
  square on screen rather than as a parallelogram."
  [orientation width height]
  (if (= :flat orientation)
    (vec (for [q (range width)
               r (range (- (long (math/floor (/ q 2.0)))) (- height (long (math/floor (/ q 2.0)))))]
           [q r]))
    (vec (for [r (range height)
               q (range (- (long (math/floor (/ r 2.0)))) (- width (long (math/floor (/ r 2.0)))))]
           [q r]))))

;; ---------------------------------------------------------------------------
;; Movement over obstacles

(defn reachable
  "Every hex within `n` steps of `start` that can actually be walked to,
  given `blocked?`. Breadth-first, so it flows around walls instead of
  through them the way `hexes-within` would."
  [start n blocked?]
  (loop [frontier [start] visited #{start} k 0]
    (if (or (= k n) (empty? frontier))
      visited
      (let [next-frontier (for [hex frontier
                                nb  (neighbors hex)
                                :when (and (not (visited nb)) (not (blocked? nb)))]
                            nb)
            next-frontier (distinct next-frontier)]
        (recur next-frontier (into visited next-frontier) (inc k))))))

(defn visible
  "The hexes within `n` steps of `from` with an unobstructed line to it.

  Field of view the simple way the guide describes: trace a line to each
  candidate and keep it only if nothing blocks the way. A blocking hex is
  itself visible -- you can see the wall, just not past it."
  [from n blocked?]
  (into #{}
        (filter (fn [target]
                  (every? (fn [hex] (or (= hex target) (not (blocked? hex))))
                          (line from target))))
        (hexes-within from n)))

;; ---------------------------------------------------------------------------
;; Offset coordinates
;;
;; Rectangular storage for a grid that is not rectangular: every other row
;; (or column) is shifted half a hex, which is how a hex map gets stored in
;; a plain 2D array. Four variants, depending on whether the odd or the
;; even line is the one pushed out.

(def ^:private EVEN 1)
(def ^:private ODD -1)

(defn- offset-parity [kind]
  (case kind (:odd-r :odd-q) ODD (:even-r :even-q) EVEN))

(defn ->offset
  "Hex to `[col row]` in one of `:odd-r`, `:even-r` (shifted rows, for
  pointy-top) or `:odd-q`, `:even-q` (shifted columns, for flat-top)."
  [kind h]
  (let [[q r] (axial h)
        off   (offset-parity kind)]
    (case kind
      (:odd-q :even-q) [q (+ r (quot (+ q (* off (bit-and q 1))) 2))]
      (:odd-r :even-r) [(+ q (quot (+ r (* off (bit-and r 1))) 2)) r])))

(defn offset->
  "`[col row]` back to a hex."
  [kind [col row]]
  (let [off (offset-parity kind)]
    (case kind
      (:odd-q :even-q) [col (- row (quot (+ col (* off (bit-and col 1))) 2))]
      (:odd-r :even-r) [(- col (quot (+ row (* off (bit-and row 1))) 2)) row])))

;; ---------------------------------------------------------------------------
;; Doubled coordinates
;;
;; The other way to store a hex grid in a rectangle: instead of shifting
;; alternate lines, double one axis so every hex lands on its own integer
;; pair. Unlike offset coordinates these have a distance formula of their
;; own, so you can work in them directly.

(defn ->doublewidth [h] (let [[q r] (axial h)] [(+ (* 2 q) r) r]))
(defn doublewidth-> [[col row]] [(/ (- col row) 2) row])

(defn ->doubleheight [h] (let [[q r] (axial h)] [q (+ (* 2 r) q)]))
(defn doubleheight-> [[col row]] [col (/ (- row col) 2)])

(defn doublewidth-distance [[c1 r1] [c2 r2]]
  (let [dcol (abs (- c1 c2)) drow (abs (- r1 r2))]
    (+ drow (max 0 (/ (- dcol drow) 2)))))

(defn doubleheight-distance [[c1 r1] [c2 r2]]
  (let [dcol (abs (- c1 c2)) drow (abs (- r1 r2))]
    (+ dcol (max 0 (/ (- drow dcol) 2)))))

;; ---------------------------------------------------------------------------
;; Layout: hexes to pixels
;;
;; The forward matrix takes axial coordinates to the plane, the inverse
;; brings a point back. They differ between orientations only by a 30
;; degree turn, which `start-angle` also applies to the corners.

(def orientations
  {:pointy {:f [(math/sqrt 3.0) (/ (math/sqrt 3.0) 2.0) 0.0 (/ 3.0 2.0)]
            :b [(/ (math/sqrt 3.0) 3.0) (/ -1.0 3.0) 0.0 (/ 2.0 3.0)]
            :start-angle 0.5}
   :flat   {:f [(/ 3.0 2.0) 0.0 (/ (math/sqrt 3.0) 2.0) (math/sqrt 3.0)]
            :b [(/ 2.0 3.0) 0.0 (/ -1.0 3.0) (/ (math/sqrt 3.0) 3.0)]
            :start-angle 0.0}})

(defn layout
  "How hexes sit on the screen: `orientation` is `:pointy` or `:flat`,
  `size` the `[x y]` radius of one hex, `origin` where `[0 0]` lands."
  ([orientation size] (layout orientation size [0.0 0.0]))
  ([orientation size origin]
   {:orientation orientation
    :size        (if (number? size) [size size] (vec size))
    :origin      (vec origin)}))

(defn ->pixel
  "The centre of hex `h` on screen."
  [{:keys [orientation size origin]} h]
  (let [[f0 f1 f2 f3] (:f (orientations orientation))
        [q r]  (axial h)
        [sx sy] size
        [ox oy] origin]
    [(+ (* (+ (* f0 q) (* f1 r)) sx) ox)
     (+ (* (+ (* f2 q) (* f3 r)) sy) oy)]))

(defn pixel->fractional
  "The fractional hex coordinate at point `p` -- cube, and generally not
  integral."
  [{:keys [orientation size origin]} p]
  (let [[b0 b1 b2 b3] (:b (orientations orientation))
        [sx sy] size
        [ox oy] origin
        x (/ (- (nth p 0) ox) sx)
        y (/ (- (nth p 1) oy) sy)
        q (+ (* b0 x) (* b1 y))
        r (+ (* b2 x) (* b3 y))]
    [q r (- (+ q r))]))

(defn pixel->hex
  "The hex containing point `p`."
  [layout p]
  (round (pixel->fractional layout p)))

(defn corner-offset
  "Offset from a hex's centre to its `i`th corner."
  [{:keys [orientation size]} i]
  (let [{:keys [start-angle]} (orientations orientation)
        [sx sy] size
        angle   (* 2.0 math/PI (/ (+ start-angle i) 6.0))]
    [(* sx (math/cos angle)) (* sy (math/sin angle))]))

(defn corners
  "The six corners of hex `h`, as a polygon."
  [layout h]
  (let [[cx cy] (->pixel layout h)]
    (mapv (fn [i]
            (let [[dx dy] (corner-offset layout i)]
              [(+ cx dx) (+ cy dy)]))
          (range 6))))

;; ---------------------------------------------------------------------------
;; Wraparound

(defn mirror-centers
  "The six centres a hexagonal map of `radius` repeats around, plus the
  origin.

  A hexagonal map tiles the plane, so a map that wraps can be built by
  translating it to each mirror centre; a hex that walks off one edge is
  the one at the matching offset on the other side."
  [radius]
  (let [seed [(+ (* 2 radius) 1) (- radius)]]
    (into [[0 0]] (take 6 (iterate #(rotate-right %) seed)))))

(defn wrap
  "The hex inside a hexagonal map of `radius` that `h` corresponds to,
  wrapping around the edges. `h` itself when it is already inside.

  Applied repeatedly, because one subtraction only brings back a hex that
  has strayed into an adjacent copy of the map. Something several map
  widths out -- anything walking in a straight line for a while -- needs
  one subtraction per copy crossed, and each one strictly reduces the
  distance from the origin, so this terminates."
  [h radius]
  (let [centers (rest (mirror-centers radius))]      ; the origin moves nothing
    (loop [h (axial h)]
      (if (<= (hex-length h) radius)
        h
        (let [closer (->> centers
                          (map #(subtract h %))
                          (apply min-key hex-length))]
          (if (< (hex-length closer) (hex-length h))
            (recur closer)
            ;; Cannot happen for a well-formed map, but a guard beats a
            ;; silent infinite loop if one is ever passed a bad radius.
            h))))))
