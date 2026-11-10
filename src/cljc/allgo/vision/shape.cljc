(ns allgo.vision.shape
  "Shape from one image: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 25 -- the cues a single picture
  gives to the 3D shape of what it shows.

    shading   a matte surface of known albedo, lit from a known direction,
              is brighter the more its normal faces the light; but a
              brightness is one number and a normal two, so the surface
              is also assumed smooth. Ikeuchi and Horn's relaxation
              (\"Numerical Shape from Shading and Occluding Boundaries\",
              Artificial Intelligence 17, 1981): each normal moved to the
              mean of its neighbors', then along the brightness error --
              an iteration a step -- with the normals at the occluding
              contour, where the surface turns away and the normal lies
              in the image plane, known from the silhouette. Normals
              are kept in stereographic coordinates (f, g), finite at
              the contour where the gradient (p, q) is infinite.
    texture   a uniform texture of round texels images as ellipses:
              foreshortened by the cosine of the slant (the angle between
              the normal and the line of sight) across the tilt (the
              direction the surface falls away), and smaller by 1 / the
              distance. Each texel's ellipse gives its slant and tilt;
              their sizes give their relative distances, and so, along
              their rays, points of the surface up to scale, and a plane
              through them.
    contour   an object standing on the ground plane, seen by a camera at
              height h with its axis level, has its foot f h / Z below
              the horizon: the lower, the nearer, and its height follows
              from its top. Where one region occludes another the
              occluder's outline runs on and the other's stops against
              it -- a T-junction, the occluder the bar of the T, covering
              about half the ground round it; the junctions, voted, give
              which region is in front of which. Objects of one known
              height place the horizon by how their size goes with
              how low their feet are, and the horizon rules out
              detections no such object could make.

  Images are `allgo.vision.image`'s; a mask is an image of booleans. In
  the shading the image's x runs right and y down, and a normal `[nx ny
  nz]` has nz toward the viewer."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.numerics.eigen :as eigen]
            [clojure.math :as math]
            [clojure.set :as set]))

;; ---------------------------------------------------------------------------
;; Shading

(defn normal
  "The unit normal with stereographic coordinates `(f, g)`: f = 2 nx / (1
  + nz), g likewise; f^2 + g^2 = 4 at the occluding contour."
  [f g]
  (let [f (double f) g (double g) d (+ 4.0 (* f f) (* g g))]
    [(/ (* 4.0 f) d) (/ (* 4.0 g) d) (/ (- 4.0 (* f f) (* g g)) d)]))

(defn- reflectance
  "The brightness of a normal at `(f, g)` with `albedo`, lit from `light`,
  and its partial derivatives by f and g: `[R Rf Rg]`. Unclipped: the
  relaxation needs the slope where the surface is in shadow too."
  [albedo [sx sy sz] f g]
  (let [f (double f) g (double g) albedo (double albedo)
        d (+ 4.0 (* f f) (* g g)) d2 (* d d)
        [nx ny nz] (normal f g)
        dot (fn [a b c] (* albedo (+ (* (double sx) a) (* (double sy) b) (* (double sz) c))))]
    [(dot nx ny nz)
     (dot (/ (* 4.0 (- d (* 2.0 f f))) d2) (/ (* -8.0 f g) d2) (/ (* -16.0 f) d2))
     (dot (/ (* -8.0 f g) d2) (/ (* 4.0 (- d (* 2.0 g g))) d2) (/ (* -16.0 g) d2))]))

(defn occluding-boundary
  "The pixels of `mask` on its edge -- a 4-neighbor outside it -- each with
  the stereographic coordinates of the normal there: in the image plane,
  pointing out of the silhouette, the way the pixels outside round it
  lie. `{index [f g]}`."
  [{:keys [w h px]}]
  (let [w (long w) h (long h)
        in? (fn [x y] (and (< -1 x w) (< -1 y h) (px (+ (* y w) x))))]
    (into {}
          (for [y (range h) x (range w)
                :when (and (in? x y) (some (fn [[dx dy]] (not (in? (+ x dx) (+ y dy)))) [[1 0] [-1 0] [0 1] [0 -1]]))
                :let [[ox oy] (reduce (fn [[a b] [dx dy]] (if (in? (+ x dx) (+ y dy)) [a b] [(+ a dx) (+ b dy)]))
                                      [0.0 0.0] (for [dy (range -2 3) dx (range -2 3)] [dx dy]))
                      l (math/hypot ox oy)]
                :when (pos? l)]
            [(+ (* y w) x) [(/ (* 2.0 ox) l) (/ (* 2.0 oy) l)]]))))

(defn- relax
  "One iteration of Ikeuchi and Horn's: every inner pixel's (f, g) the mean
  of its four neighbors', moved along the gradient of the reflectance by
  `lambda` times the brightness error there."
  [{:keys [w px]} inner albedo light lambda {:keys [f g] :as st}]
  (let [w (long w) lambda (double lambda)
        mean (fn [v ^long i] (* 0.25 (+ (double (v (inc i))) (double (v (dec i))) (double (v (+ i w))) (double (v (- i w))))))
        [f' g'] (reduce (fn [[f' g'] i]
                          (let [fb (mean f i) gb (mean g i)
                                [r rf rg] (reflectance albedo light fb gb)
                                e (* lambda (- (double (px i)) (double r)))]
                            [(assoc! f' i (+ fb (* e (double rf)))) (assoc! g' i (+ gb (* e (double rg))))]))
                        [(transient f) (transient g)] inner)]
    (-> st (assoc :f (persistent! f') :g (persistent! g')) (update :iteration inc))))

(defn shading-seq
  "Shape from shading: the normals of the surface in `image` inside
  `mask`, of known `albedo` and lit from the unit direction `light`, an
  iteration a step. Each state `{:f :g :iteration}`, `:f` and `:g` vectors
  over the pixels: the occluding contour's fixed, the inside's starting
  facing the viewer. `:lambda` (default 2) weighs brightness against
  smoothness; more is worse, and by 8 the relaxation diverges."
  [image mask albedo light & {:keys [lambda] :or {lambda 2.0}}]
  (let [{:keys [w h]} image
        n (* (long w) (long h))
        boundary (occluding-boundary mask)
        inner (vec (for [i (range n) :when (and ((:px mask) i) (not (boundary i)))] i))
        start {:f (reduce (fn [v [i [f _]]] (assoc v i f)) (vec (repeat n 0.0)) boundary)
               :g (reduce (fn [v [i [_ g]]] (assoc v i g)) (vec (repeat n 0.0)) boundary)
               :iteration 0}]
    (iterate (partial relax image inner albedo light lambda) start)))

(defn normals
  "The normals of a shading state, an image of `[nx ny nz]` inside `mask`
  and nil outside."
  [{:keys [f g]} {:keys [w h px]}]
  {:w w :h h :px (mapv (fn [i] (when (px i) (normal (f i) (g i)))) (range (* (long w) (long h))))})

;; ---------------------------------------------------------------------------
;; Texture

(defn blobs
  "The 4-connected regions of pixels at least `threshold` bright, each the
  vector of its pixels' centers `[x y]`."
  [{:keys [w h px]} threshold]
  (let [w (long w) h (long h) n (* w h)
        bright? (fn [^long i] (>= (double (px i)) (double threshold)))]
    (loop [i 0 seen (transient #{}) out []]
      (cond (= i n) out
            (or (seen i) (not (bright? i))) (recur (inc i) seen out)
            :else (let [[region seen]
                        (loop [stack [i] seen (conj! seen i) region []]
                          (if-let [j (peek stack)]
                            (let [x (mod j w) y (quot j w)
                                  next (for [[dx dy] [[1 0] [-1 0] [0 1] [0 -1]]
                                             :let [a (+ x dx) b (+ y dy)]
                                             :when (and (< -1 a w) (< -1 b h))
                                             :let [k (+ (* b w) a)]
                                             :when (and (not (seen k)) (bright? k))]
                                         k)]
                              (recur (into (pop stack) next) (reduce conj! seen next) (conj region j)))
                            [region seen]))]
                    (recur (inc i) seen (conj out (mapv (fn [j] [(+ 0.5 (mod j w)) (+ 0.5 (quot j w))]) region))))))))

(defn- ellipse
  "The ellipse with the second moments of `pts`: `[center major minor
  angle]`, the axes full lengths (four standard deviations, a uniform
  ellipse's), `angle` the minor axis's direction."
  [pts]
  (let [n (double (count pts))
        [mx my] (mapv #(/ (double %) n) (reduce (fn [[a b] [x y]] [(+ a x) (+ b y)]) [0.0 0.0] pts))
        [sxx sxy syy] (mapv #(/ (double %) n)
                            (reduce (fn [[a b c] [x y]] (let [dx (- x mx) dy (- y my)] [(+ a (* dx dx)) (+ b (* dx dy)) (+ c (* dy dy))]))
                                    [0.0 0.0 0.0] pts))
        half-trace (* 0.5 (+ sxx syy))
        disc (math/sqrt (max 0.0 (- (* half-trace half-trace) (- (* sxx syy) (* sxy sxy)))))]
    [[mx my]
     (* 4.0 (math/sqrt (+ half-trace disc)))
     (* 4.0 (math/sqrt (max 0.0 (- half-trace disc))))
     (+ (* 0.5 (math/atan2 (* 2.0 sxy) (- sxx syy))) (/ math/PI 2))]))

(defn- direction
  "The unit ray, in camera coordinates, through pixel `[u v]` of a camera
  with focal length `f` and center `[cx cy]`."
  [f [cx cy] [u v]]
  (v/normalize [(/ (- u cx) f) (/ (- v cy) f) 1.0]))

(defn texel
  "What the texel whose pixels are `pts` says, seen by a camera of focal
  length `f` and center `center`. It is measured as a camera turned to
  look straight at it would see it -- each pixel's ray projected onto the
  plane square to the texel's -- since off the axis the image plane is
  itself oblique to the line of sight. `{:pixel :ray :major :minor :slant
  :tilt}`: its center and ray, its ellipse's axes (in units of that
  plane, so the major axis is the texel's size over its distance), the
  slant acos(minor / major), and the tilt, the direction in that plane
  -- x right, y down -- the surface falls away, up to a half turn;
  `:frame` that plane's two axes, in camera coordinates."
  [f center pts]
  (let [n (double (count pts))
        pixel (mapv #(/ (double %) n) (reduce (fn [[a b] [x y]] [(+ a x) (+ b y)]) [0.0 0.0] pts))
        ray (direction f center pixel)
        e1 (v/normalize (v/cross [0.0 1.0 0.0] ray))
        e2 (v/cross ray e1)
        local (mapv (fn [p] (let [r (direction f center p) s (v/dot r ray)] [(/ (v/dot r e1) s) (/ (v/dot r e2) s)])) pts)
        [_ major minor angle] (ellipse local)]
    {:pixel pixel :ray ray :major major :minor minor
     :slant (math/acos (min 1.0 (/ minor major)))
     :tilt (mod angle math/PI)
     :frame [e1 e2]}))

(defn texels
  "The texels of `image` -- regions brighter than `threshold` of at least
  `min-size` pixels, not touching its edge -- as `texel` measures them."
  [{:keys [w h] :as image} f center threshold min-size]
  (for [pts (blobs image threshold)
        :when (and (>= (count pts) (long min-size))
                   (every? (fn [[x y]] (and (< 1 x (dec (long w))) (< 1 y (dec (long h))))) pts))]
    (texel f center pts)))

(defn plane-from-texels
  "The unit normal, in camera coordinates and facing the camera, of the
  plane through the texels placed along their rays at distances
  proportional to 1 / their size -- the texture's scale cue, needing no
  slant at all. The least-squares plane: the direction their scatter is
  least along."
  [texels]
  (let [ps (mapv (fn [{:keys [ray major]}] (v/scale ray (/ 1.0 (double major)))) texels)
        m (v/scale (reduce v/add ps) (/ 1.0 (count ps)))
        scatter (reduce (fn [C p] (let [d (v/sub p m)] (mapv (fn [i] (mapv (fn [j] (+ (double (get-in C [i j])) (* (double (d i)) (double (d j))))) (range 3))) (range 3))))
                        (vec (repeat 3 [0.0 0.0 0.0])) ps)
        n (mapv #(% 0) (:vectors (eigen/jacobi scatter)))]
    (if (pos? (v/dot n m)) (v/negate n) n)))

;; ---------------------------------------------------------------------------
;; Contour

(defn ground-distance
  "How far ahead, along a level axis, a point on the ground plane is when
  it lands at image row `v`, for a camera of focal length `f` and center
  row `cy` at height `camera-height`: f h / (v - cy). Nil at or above the
  horizon."
  [f cy camera-height v]
  (let [below (- (double v) (double cy))]
    (when (pos? below) (/ (* (double f) (double camera-height)) below))))

(defn height-on-ground
  "The height of an object standing on the ground plane whose foot lands
  at row `v-foot` and top at row `v-top`, for a level camera at
  `camera-height` with center row `cy`: h (v-foot - v-top) / (v-foot -
  cy). The focal length cancels."
  [cy camera-height v-foot v-top]
  (/ (* (double camera-height) (- (double v-foot) (double v-top))) (- (double v-foot) (double cy))))

(defn- line-fit
  "The least-squares line h = a v + b through `points` `[v h]`: `[a b]`."
  [points]
  (let [n (count points)
        sv (reduce + (map first points)) sh (reduce + (map second points))
        svv (reduce + (map #(* (first %) (first %)) points)) svh (reduce + (map #(* (first %) (second %)) points))
        a (/ (- (* n svh) (* sv sh)) (- (* n svv) (* sv sv)))]
    [a (/ (- sh (* a sv)) n)]))

(defn- horizon-fit [kept]
  (let [[a b] (line-fit (map (juxt :foot :height) kept))]
    {:kept kept :scale a :horizon (/ (- b) a)}))

(defn- off-line
  "How far `d`'s height is from the line's, as a fraction of it."
  [{:keys [horizon scale]} {:keys [foot height]}]
  (/ (abs (- height (* scale (- foot horizon)))) height))

(defn horizon-seq
  "The horizon from `detections` of objects all of one height standing on
  the ground -- pedestrians -- each `{:foot :height}`, the row of its foot
  and its height in the image, for a level camera: an object H tall at
  depth Z has its foot f c / Z below the horizon, c the camera's height,
  and is f H / Z tall, so its height is (H / c) (v-foot - horizon) -- a
  line, zero at the horizon. Section 25.6.4's two ways round: the
  detections place the horizon, and the horizon rules out detections no
  pedestrian could make, too big or too small for where their feet are.

  One wrong detection far from the rest can drag a least-squares line
  anywhere, so the start is the consensus: of the lines through every
  pair of detections, the one most of the others lie within `tolerance`
  of (as a fraction of their heights). A step fits the line by least
  squares to the detections kept and keeps those within the tolerance of
  it; the states `{:kept :dropped :horizon :scale}`, until the kept stop
  changing."
  [detections tolerance]
  (let [ds (vec detections)
        within (fn [fit] (filterv #(<= (off-line fit %) tolerance) ds))
        pairs (for [i (range (count ds)) j (range (inc i) (count ds))
                    :when (not= (:foot (ds i)) (:foot (ds j)))]
                (horizon-fit [(ds i) (ds j)]))
        start (apply max-key (comp count within) pairs)
        step (fn [{:keys [kept]}]
               (let [fit (horizon-fit kept)]
                 (assoc fit :kept (within fit) :dropped (filterv #(> (off-line fit %) tolerance) ds))))
        states (iterate step (step {:kept (within start)}))]
    (into [(first states)] (map second (take-while (fn [[a b]] (not= (:kept a) (:kept b))) (partition 2 1 states))))))

(defn horizon
  "The last of `horizon-seq`: `{:kept :dropped :horizon :scale}`."
  [detections tolerance]
  (peek (horizon-seq detections tolerance)))

(defn t-junctions
  "Where three regions of `labels` meet: the corners shared by a 2 by 2
  block of pixels holding three different labels, `{:at [x y] :regions
  #{a b c}}`, the corner in pixel-edge coordinates."
  [{:keys [w h px]}]
  (let [w (long w)]
    (for [y (range 1 h) x (range 1 w)
          :let [block (hash-set (px (+ (* (dec y) w) (dec x))) (px (+ (* (dec y) w) x))
                                (px (+ (* y w) (dec x))) (px (+ (* y w) x)))]
          :when (= 3 (count block))]
      {:at [x y] :regions block})))

(defn- branches
  "The directions from the junction at `[jx jy]` along the boundary
  between each pair of its `regions`, within `radius`: `{#{a b} [dx
  dy]}`, the sum of the unit directions to the midpoints of the pixel
  edges between them."
  [{:keys [w h px]} [jx jy] regions radius]
  (let [w (long w) h (long h) r (long (math/ceil radius))]
    (reduce (fn [m [pair [dx dy]]] (update m pair (fn [[a b]] [(+ (double (or a 0.0)) dx) (+ (double (or b 0.0)) dy)])))
            {}
            (for [y (range (- (long jy) r) (+ (long jy) r)) x (range (- (long jx) r) (+ (long jx) r))
                  [sx sy] [[1 0] [0 1]]
                  :let [x2 (+ x sx) y2 (+ y sy)]
                  :when (and (< -1 x w) (< -1 y h) (< x2 w) (< y2 h))
                  :let [a (px (+ (* y w) x)) b (px (+ (* y2 w) x2))]
                  :when (and (not= a b) (regions a) (regions b))
                  :let [ox (- (+ x 0.5 (* 0.5 sx)) (double jx)) oy (- (+ y 0.5 (* 0.5 sy)) (double jy))
                        d (math/hypot ox oy)]
                  :when (<= 1.0 d (double radius))]
              [#{a b} [(/ ox d) (/ oy d)]]))))

(defn occlusions
  "What each T-junction of `labels` says. Of the three boundaries meeting
  there, the two that run on most nearly straight through it are the
  occluder's outline -- the bar of the T -- and the third, the stem, is
  the boundary between what it hides and what is behind that. The region
  the two straight ones share is in front of the other two. A junction
  where no pair runs clearly straight through -- at less than about 145
  degrees, or hardly straighter than the next pair, as where an outline
  turns a corner or curves tightly within `radius` of it -- says nothing.
  `[{:at :front :behind} ...]`, `:behind` the set of the other two;
  `radius` how far along the boundaries to look."
  [labels radius]
  (for [{:keys [at regions]} (t-junctions labels)
        :let [bs (branches labels at regions radius)]
        :when (= 3 (count bs))
        :let [cosine (fn [[a b] [c d]] (/ (+ (* a c) (* b d)) (* (math/hypot a b) (math/hypot c d))))
              ks (vec (keys bs))
              [[c1 p q] [c2]] (sort-by first (for [[i j] [[0 1] [0 2] [1 2]]]
                                               [(cosine (bs (ks i)) (bs (ks j))) (ks i) (ks j)]))]
        :when (and (< (double c1) -0.8) (> (- (double c2) (double c1)) 0.4))
        :let [front (first (set/intersection p q))]]
    {:at at :front front :behind (disj regions front)}))

(defn in-front
  "What the T-junctions of `labels` decide: the pairs `[a b]`, region a in
  front of region b, where more junctions say so than the reverse. A
  partial order at best: regions no junction joins stay unrelated."
  [labels radius]
  (let [votes (frequencies (for [{:keys [front behind]} (occlusions labels radius) b behind] [front b]))]
    (set (for [[[a b] n] votes :when (> (long n) (long (votes [b a] 0)))] [a b]))))

(defn depth-order
  "The regions of `labels` nearest first, by `in-front`: each comes next
  once nothing left is in front of it -- a topological order, so regions
  the junctions leave unrelated come in order of label, and on a cycle
  the region fewest are in front of comes first. The `background`, if
  given, last."
  [labels radius & [background]]
  (let [pairs (in-front labels radius)]
    (loop [left (disj (set (:px labels)) background) order []]
      (if (empty? left)
        (cond-> order (some? background) (conj background))
        (let [hidden-by (fn [a] (count (filter #(pairs [% a]) left)))
              next (first (sort-by (fn [a] [(hidden-by a) a]) left))]
          (recur (disj left next) (conj order next)))))))
