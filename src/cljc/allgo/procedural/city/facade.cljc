(ns allgo.procedural.city.facade
  "Facades as layered grids: Parish and Mueller's section 4.3.

  Their observations: a facade is one or several grids laid over each
  other or nested, most of whose cells do the same thing (a window, a
  door); some cells change the cells around them; and irregularities
  take whole rows and columns, not single cells. So the unit is the
  *interval group* -- an ordered set of non-overlapping intervals along
  one axis -- and a change to a row or a column is a change to one
  interval.

  Two interval groups and an evaluation function make a *layer*: a point
  `(s, t)` is active where the function of whether `s` falls in an
  interval of the first group and `t` in one of the second is true --
  AND for a grid of windows, OR for crossing bands -- and an active point
  gets the layer's color function's answer. Layers are *stacked*: a point
  not active in one is asked of the next, down to the wall. And a layer's
  evaluation function may ask another layer, which is the paper's
  functions between layers: brick lintels that exist only above windows,
  a door that takes the place of one window.

  The paper's `Grid(3).scale(0.95)` is `(scale (grid 3 length) 0.95)`
  here, and its PULSE is the test for falling inside an interval. `s`
  runs along the wall and `t` up it, in meters, so a style fits any size
  of wall: the number of bays and floors follows from the size."
  (:require [allgo.procedural.noise :as noise]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Interval groups

(defn grid
  "`n` equal intervals filling `[0, length]`."
  [n length]
  (let [n (max 1 (long n)) w (/ (double length) n)]
    (mapv (fn [i] [(* i w) (* (inc i) w)]) (range n))))

(defn rows
  "Intervals of height `h` from `from` up to `to`: floors."
  [h from to]
  (let [h (double h) from (double from)
        n (max 0 (long (math/floor (/ (- (double to) from) h))))]
    (mapv (fn [i] [(+ from (* i h)) (+ from (* (inc i) h))]) (range n))))

(defn scale
  "Every interval shrunk about its middle to `k` of its width, or, with
  `[below above]`, to that fraction of it from each end -- a window sits
  above a sill, not in the middle of a floor."
  [group k]
  (mapv (fn [[a b]]
          (let [w (- (double b) (double a))]
            (if (vector? k)
              (let [[lo hi] k] [(+ a (* w (double lo))) (- b (* w (double hi)))])
              (let [m (* 0.5 (+ a b)) k (double k)] [(- m (* 0.5 k w)) (+ m (* 0.5 k w))]))))
        group))

(defn locate
  "Which interval of `group` holds `x`, and how far through it, as
  `[index fraction]`, or nil: the pulse."
  [group x]
  (let [x (double x)]
    (loop [lo 0 hi (dec (count group))]
      (when (<= lo hi)
        (let [mid (quot (+ lo hi) 2)
              [a b] (nth group mid)]
          (cond (< x (double a)) (recur lo (dec mid))
                (>= x (double b)) (recur (inc mid) hi)
                :else [mid (/ (- x (double a)) (- (double b) (double a)))]))))))

;; ---------------------------------------------------------------------------
;; Layers and stacks

(defn layer
  "A layer of the interval groups `s` and `t`. `eval` is `:and`, `:or`,
  or a function of the two pulses (each `[index fraction]` or nil) and
  the point; `color` is a function of a map of the cell -- `:i :j :u :v`,
  where `u` and `v` are the fractions through it, or nil along a band --
  and the point, `:s :t`."
  [s t eval color]
  {:s s :t t :eval eval :color color})

(defn active
  "The cell of `layer` at `(s, t)` if the point is active in it, or nil."
  [{s-group :s t-group :t eval :eval} s t]
  (let [ps (locate s-group s)
        pt (locate t-group t)
        on (case eval
             :and (and ps pt)
             :or (or ps pt)
             (eval ps pt s t))]
    (when on
      {:i (first ps) :u (second ps) :j (first pt) :v (second pt) :s s :t t})))

(defn evaluate
  "The color of the point `(s, t)`: from the first layer of `stack` it is
  active in, or `base`."
  [stack base s t]
  (or (some (fn [{:keys [color] :as l}]
              (when-let [cell (active l s t)] (color cell)))
            stack)
      (base s t)))

;; ---------------------------------------------------------------------------
;; Styles

(def ^:private table (memoize noise/permutation))

(defn- cell-random
  "A number in [0, 1) belonging to one cell of one layer: the same every
  time it is asked, as a procedural texture's has to be. Hashed through
  Perlin's permutation table, as his lattice is: a seeded stream folded
  over neighboring cells comes out a constant step apart, and lights
  whole runs of windows at once."
  ^double [seed tag i j]
  (let [^ints p (table (long seed))
        a (aget p (+ (aget p (+ (aget p (bit-and (long (or i 0)) 255)) (bit-and (long (or j 0)) 255)))
                     (bit-and (long tag) 255)))]
    (/ (+ a (/ (aget p (inc a)) 256.0)) 256.0)))

(defn- shade [c k] (mapv #(min 1.0 (* (double %) (double k))) c))

(defn- mix [a b t] (mapv #(+ (double %1) (* (double t) (- (double %2) (double %1)))) a b))

(defn- wall-tone
  "A wall that is not one flat color: a little variation from course to
  course."
  [c seed course]
  (fn [_ t] (shade c (+ 0.94 (* 0.08 (cell-random seed 9 0 (long (math/floor (/ t course)))))))))

(defn- window-glass
  "A pane: night-dark or lit, by cell, with a blind drawn some way down
  some of them."
  [seed tag lit dark]
  (fn [{:keys [i j v]}]
    (let [r (cell-random seed tag i j)
          blind (cell-random seed (inc tag) i j)
          pane (if (< r lit) [0.98 0.86 0.55] (shade dark (+ 0.8 (* 0.4 r))))]
      (if (and (< blind 0.35) v (> v (- 1.0 (* 0.8 blind)))) [0.86 0.84 0.78] pane))))

(def styles
  "The kinds of facade, each `(fn [width height floor seed])` returning
  `[stack base]`."
  {:tower
   (fn [w h f seed]
     ;; A curtain wall: glass between thin mullions and spandrels.
     (let [bays (grid (math/round (/ w 1.6)) w)
           floors (rows f 0.0 h)
           glass (window-glass seed 1 0.18 [0.25 0.36 0.45])]
       [[(layer (scale bays 0.9) (scale floors [0.12 0.04]) :and glass)]
        (fn [_ _] [0.62 0.66 0.70])]))

   :base
   (fn [w h f seed]
     ;; A tower's podium: tall shop glass, stone above it.
     (let [bays (grid (math/round (/ w 4.5)) w)]
       [[(layer (scale bays 0.86) [[0.3 (* 0.85 (min h (* 1.6 f)))]] :and
                (window-glass seed 3 0.5 [0.2 0.24 0.26]))
         (layer (scale bays 0.55) (scale (rows f (* 1.6 f) h) [0.25 0.15]) :and
                (window-glass seed 5 0.25 [0.2 0.24 0.3]))]
        (wall-tone [0.72 0.68 0.6] seed 0.6)]))

   :office
   (fn [w h f seed]
     ;; Ribbon windows: rows of glass between bands of concrete.
     (let [bays (grid (math/round (/ w 3.2)) w)
           floors (scale (rows f f h) [0.3 0.15])]
       [[(layer (scale bays 0.92) floors :and (window-glass seed 7 0.22 [0.18 0.24 0.3]))
         ;; The ground floor: an entrance in one bay, glass in the rest.
         (layer bays [[0.0 (* 0.85 f)]] :and
                (fn [{:keys [i]}]
                  (if (= i (quot (count bays) 2)) [0.12 0.12 0.13] [0.3 0.36 0.38])))]
        (wall-tone [0.78 0.77 0.74] seed 0.5)]))

   :shops
   (fn [w h f seed]
     (let [units (grid (max 1 (math/round (/ w 6.0))) w)
           awning (fn [{:keys [i]}]
                    (nth [[0.7 0.2 0.18] [0.2 0.45 0.3] [0.85 0.65 0.25] [0.25 0.3 0.6]]
                         (long (* 4 (cell-random seed 11 i 0)))))]
       [[;; A canvas awning over each shopfront, then the glass.
         (layer (scale units 0.94) [[(* 0.78 f) (* 0.95 f)]] :and awning)
         (layer (scale units [0.08 0.08]) [[0.2 (* 0.76 f)]] :and
                (window-glass seed 13 0.6 [0.22 0.26 0.28]))
         (layer (scale (grid (math/round (/ w 2.6)) w) 0.5) (scale (rows f f h) [0.3 0.2]) :and
                (window-glass seed 15 0.3 [0.2 0.22 0.28]))]
        (wall-tone [0.68 0.56 0.46] seed 0.4)]))

   :house
   (fn [w h f seed]
     ;; Brick, two windows a floor, a door; and the paper's example of a
     ;; function between layers -- a lintel of darker brick, only over
     ;; each opening.
     (let [bays (grid (max 2 (math/round (/ w 3.4))) w)
           floors (rows f 0.0 h)
           door-bay (long (* (count bays) (cell-random seed 17 0 0)))
           openings (layer (scale bays [0.3 0.3]) (scale floors [0.3 0.25]) :and
                           (window-glass seed 19 0.2 [0.2 0.23 0.28]))
           door (layer (scale bays [0.36 0.36]) [[0.0 (* 0.78 f)]]
                       (fn [ps pt _ _] (and ps pt (= (first ps) door-bay)))
                       (fn [_] [0.32 0.2 0.14]))
           lintel (layer bays floors
                         (fn [ps pt s t]
                           (and ps pt
                                (or (active openings s (- t 0.3)) (active door s (- t 0.3)))
                                (not (active openings s t)) (not (active door s t))))
                         (fn [_] [0.42 0.2 0.15]))
           brick (cell-random seed 23 0 0)]
       [[door openings lintel]
        (let [c (mix [0.66 0.36 0.27] [0.82 0.74 0.62] brick)]
          (fn [s t]
            ;; Courses of brick, a half brick offset each course.
            (let [course (long (math/floor (/ t 0.075)))
                  b (long (math/floor (/ (+ s (* 0.11 (mod course 2))) 0.22)))]
              (shade c (+ 0.88 (* 0.18 (cell-random seed 29 b course)))))))]))

   :plant
   (fn [_ _ _ _]
     [[(layer (grid 40 40.0) [[0.0 40.0]] :and (fn [{:keys [i]}] (shade [0.55 0.56 0.58] (+ 0.9 (* 0.1 (mod i 2))))))]
      (fn [_ _] [0.5 0.5 0.52])])})

(defn facade
  "The color of a facade of `style`, `width` by `height` meters with
  floors `floor` high, at `(s, t)` meters along and up it."
  [style width height floor seed s t]
  (let [[stack base] ((get styles style (:office styles)) (double width) (double height) (double floor) seed)]
    (evaluate stack base s t)))

(defn painter
  "A function of `(s, t)` for one facade, with its layers built once."
  [style width height floor seed]
  (let [[stack base] ((get styles style (:office styles)) (double width) (double height) (double floor) seed)]
    (fn [s t] (evaluate stack base s t))))
