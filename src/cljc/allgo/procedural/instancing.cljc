(ns allgo.procedural.instancing
  "Procedural geometric instancing: a model as a few rules that place
  copies of each other, expanded into as many primitives as the viewer
  can see.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, Hart's
  \"Procedural Synthesis of Geometry\". Ordinary instancing is a scene
  graph whose nodes say \"draw that, moved here\". Two changes make it a
  way of *generating* geometry rather than storing it:

  - **The graph may have cycles.** A rule that places shrunken copies of
    itself is recursion, and its limit is a fractal: four half-size
    copies of a tetrahedron is Sierpinski's, twenty third-size copies of a
    cube is Menger's sponge. This is an iterated function system written
    as a scene graph -- the attractor is what is left after infinitely
    many expansions, and any finite expansion is an approximation of it.
  - **Instances carry parameters,** evaluated when the instance is
    expanded. A branch that is told its own random seed can choose how
    many twigs to sprout and at what angles, so every branch of a tree is
    different though there is one rule for all of them.

  A rule is a map:

      {:parts  (fn [params] -> parts)   ; what one copy is made of
       :proxy  (fn [params] -> parts)   ; what to draw in its place, when
                                        ; it is not worth expanding
       :radius r :center [x y z]}       ; a sphere, in the rule's own
                                        ; frame, holding all of it

  and a part is either a primitive, `{:shape :cylinder :transform m
  :color [r g b]}`, or another instance, `{:instance :branch :transform m
  :params p}`. A model is `{:root :rules :params}`.

  ## Level of detail

  The expansion is the model's only cost, so it is also where detail is
  decided. Each instance is measured before it is expanded: the radius of
  its bounding sphere, carried through its transform, and -- given an
  `:eye` -- divided by its distance, which is the angle it fills. Too
  small, and it is drawn as its proxy instead: a branch too far away to
  resolve becomes a blob of foliage, a sponge too fine to see becomes a
  cube. Nothing smaller than the cutoff is ever generated, so the cost
  follows what is on the screen rather than what is in the model.

  For that to hold, a parameter must not depend on the order in which the
  graph is walked. Every random choice here is drawn from a stream keyed
  by the instance's own seed, and children get seeds from their parent's
  stream, so a branch is the same branch whether the tree around it is
  expanded to two levels or ten. Detail adds twigs; it does not move the
  boughs.

  ## Transforms

  Affine, as a row-major 3x4 vector `[a b c x, d e f y, g h i z]`: the
  linear part's rows, each followed by that row of the translation."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Affine transforms

(def identity-transform [1.0 0.0 0.0 0.0, 0.0 1.0 0.0 0.0, 0.0 0.0 1.0 0.0])

(defn translation [[x y z]]
  [1.0 0.0 0.0 (double x), 0.0 1.0 0.0 (double y), 0.0 0.0 1.0 (double z)])

(defn scaling
  "Uniform by `s`, or by `[sx sy sz]` along the axes."
  [s]
  (let [[sx sy sz] (if (number? s) [s s s] s)]
    [(double sx) 0.0 0.0 0.0, 0.0 (double sy) 0.0 0.0, 0.0 0.0 (double sz) 0.0]))

(defn rotation
  "About the unit `axis` by `angle` radians, counterclockwise looking down
  the axis. Rodrigues' formula."
  [[x y z] angle]
  (let [c (math/cos angle) s (math/sin angle) t (- 1.0 c)]
    [(+ c (* t x x)) (- (* t x y) (* s z)) (+ (* t x z) (* s y)) 0.0
     (+ (* t x y) (* s z)) (+ c (* t y y)) (- (* t y z) (* s x)) 0.0
     (- (* t x z) (* s y)) (+ (* t y z) (* s x)) (+ c (* t z z)) 0.0]))

(defn frame
  "The transform taking the unit axes to `x`, `y` and `z` and the origin
  to `origin`: its columns."
  [[xx xy xz] [yx yy yz] [zx zy zz] [ox oy oz]]
  [(double xx) (double yx) (double zx) (double ox)
   (double xy) (double yy) (double zy) (double oy)
   (double xz) (double yz) (double zz) (double oz)])

(defn compose
  "`a` after `b`: the transform that applies `b` and then `a`."
  ([a b]
   (let [[a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11] a
         [b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11] b]
     [(+ (* a0 b0) (* a1 b4) (* a2 b8))
      (+ (* a0 b1) (* a1 b5) (* a2 b9))
      (+ (* a0 b2) (* a1 b6) (* a2 b10))
      (+ (* a0 b3) (* a1 b7) (* a2 b11) a3)
      (+ (* a4 b0) (* a5 b4) (* a6 b8))
      (+ (* a4 b1) (* a5 b5) (* a6 b9))
      (+ (* a4 b2) (* a5 b6) (* a6 b10))
      (+ (* a4 b3) (* a5 b7) (* a6 b11) a7)
      (+ (* a8 b0) (* a9 b4) (* a10 b8))
      (+ (* a8 b1) (* a9 b5) (* a10 b9))
      (+ (* a8 b2) (* a9 b6) (* a10 b10))
      (+ (* a8 b3) (* a9 b7) (* a10 b11) a11)]))
  ([a b & more] (reduce compose (compose a b) more)))

(defn transform-point [m [x y z]]
  (let [[a b c tx d e f ty g h i tz] m]
    [(+ (* a x) (* b y) (* c z) tx)
     (+ (* d x) (* e y) (* f z) ty)
     (+ (* g x) (* h y) (* i z) tz)]))

(defn stretch
  "The most `m` lengthens any of the axes: the length of its longest
  column. For the rotations and uniform scalings models are built from it
  is the exact factor by which a sphere's radius grows; for a shear it is
  an underestimate, which a bound can absorb."
  ^double [m]
  (let [[a b c _ d e f _ g h i _] m]
    (max (math/sqrt (+ (* a a) (* d d) (* g g)))
         (math/sqrt (+ (* b b) (* e e) (* h h)))
         (math/sqrt (+ (* c c) (* f f) (* i i))))))

;; ---------------------------------------------------------------------------
;; Expansion

(defn- size-class
  "Which queue an instance of apparent size `s` waits in: quarter steps of
  the natural log, so that siblings of one size are served together. An
  eye inside the bound measures it negative, which is as big as it gets."
  ^long [^double s]
  (cond (neg? s) 100000
        (zero? s) -100000
        :else (long (math/floor (* 4.0 (math/log s))))))

(defn expand
  "The primitives of `model`, each `{:shape :transform :color :depth}`,
  and what it took to find them:

      {:primitives [...] :instances n :proxies n :depth deepest}

  Options:
  - `:min-size` -- below this an instance is drawn as its proxy. In world
    units, or with `:eye` in radians: the angle the bounding sphere fills.
  - `:eye` -- a point to measure that angle from.
  - `:max-depth` -- a hard cap on recursion, past which every instance is
    a proxy (default 24).
  - `:limit` -- a budget of primitives (default 500000). Each instance
    still waiting is counted as one, so a proxy of several parts can
    carry the total past it by that much.
  - `:transform` -- where to put the whole model.

  The largest instances are expanded first, so that when the budget runs
  out it is the finest detail that is left as proxies, evenly, all over
  the model. Expanding depth first would spend all of it on the first
  corner and leave the rest of the model as a few big stand-ins.

  A proxy's parts are drawn as they are: they should be primitives."
  [{:keys [root rules params]}
   {:keys [min-size eye max-depth limit transform]
    :or {min-size 0.0 max-depth 24 limit 500000 transform identity-transform}}]
  (let [min-size (double min-size)
        size (fn [m {:keys [radius center] :or {radius 1.0 center v/zero}}]
               (let [r (* (double radius) (stretch m))]
                 (if eye
                   (/ r (max 1e-9 (- (v/distance eye (transform-point m center)) r)))
                   r)))
        item (fn [kind params m depth]
               (let [s (size m (get rules kind))]
                 {:kind kind :params params :m m :depth depth :size s}))
        enqueue (fn [q {:keys [size] :as it}]
                  (update q (size-class size) (fnil conj []) it))
        primitive (fn [m depth part]
                    {:shape (:shape part)
                     :transform (if-let [t (:transform part)] (compose m t) m)
                     :color (:color part) :depth depth})]
    ;; Classes in descending order, each a stack: the biggest first.
    (loop [q (enqueue (sorted-map-by >) (item root params transform 0))
           waiting 1
           out (transient [])
           instances 0 proxies 0 deepest 0]
      (if-let [[k stack] (first q)]
        (let [{:keys [kind params m depth size]} (peek stack)
              q (if (= 1 (count stack)) (dissoc q k) (assoc q k (pop stack)))
              {:keys [parts proxy]} (get rules kind)
              s (double size)
              ;; Everything still waiting will cost at least a primitive.
              small? (or (>= depth max-depth)
                         (>= (+ (count out) waiting) limit)
                         (and (pos? s) (< s min-size)))
              deepest (max deepest depth)]
          (if small?
            (recur q (dec waiting)
                   (reduce #(conj! %1 (primitive m depth %2)) out (when proxy (proxy params)))
                   (inc instances) (inc proxies) deepest)
            (let [[q waiting out]
                  (reduce (fn [[q waiting out] part]
                            (if-let [child (:instance part)]
                              [(enqueue q (item child (:params part)
                                                (compose m (or (:transform part) identity-transform))
                                                (inc depth)))
                               (inc waiting) out]
                              [q waiting (conj! out (primitive m depth part))]))
                          [q (dec waiting) out]
                          (parts params))]
              (recur q waiting out (inc instances) proxies deepest))))
        {:primitives (persistent! out)
         :instances instances
         :proxies proxies
         :depth deepest}))))

;; ---------------------------------------------------------------------------
;; Models

(def ^:private tetrahedron-corners
  ;; Alternate corners of a cube, on the unit sphere: what three's
  ;; TetrahedronGeometry of radius 1 is built on.
  (mapv #(v/scale % (/ 1.0 (math/sqrt 3.0)))
        [[1 1 1] [1 -1 -1] [-1 1 -1] [-1 -1 1]]))

(defn- corner-color [p]
  (let [[x y z] p] [(+ 0.55 (* 0.4 x)) (+ 0.55 (* 0.4 y)) (+ 0.55 (* 0.4 z))]))

(defn sierpinski
  "Four half-size copies of a tetrahedron, one at each corner, forever. The
  copies touch at their corners and the limit has no volume at all, which
  is why each level of expansion leaves the same silhouette with half as
  much behind it."
  []
  {:root :tetra
   :params {:color [0.9 0.7 0.4]}
   :rules
   {:tetra {:radius 1.0
            :parts (fn [{:keys [color]}]
                     (for [c tetrahedron-corners]
                       {:instance :tetra
                        :transform (compose (translation (v/scale c 0.5)) (scaling 0.5))
                        :params {:color (mapv #(+ (* 0.6 %1) (* 0.4 %2)) color (corner-color c))}}))
            :proxy (fn [{:keys [color]}] [{:shape :tetrahedron :color color}])}}})

(def ^:private menger-offsets
  ;; The 27 subcubes less the six face centers and the middle: those with
  ;; at most one zero coordinate.
  (vec (for [i [-1 0 1] j [-1 0 1] k [-1 0 1]
             :when (<= (count (filter zero? [i j k])) 1)]
         [i j k])))

(defn menger
  "Menger's sponge: a cube cut into 27 and the middle cross taken out,
  leaving 20 that are each a sponge. Its volume goes to nothing and its
  surface to infinity, and a cube at the cutoff is exactly right for it."
  []
  {:root :sponge
   :rules
   {:sponge {:radius (math/sqrt 3.0)
             :parts (fn [_]
                      (for [o menger-offsets]
                        {:instance :sponge
                         :transform (compose (translation (v/scale o (/ 2.0 3.0)))
                                             (scaling (/ 1.0 3.0)))}))
             :proxy (fn [_] [{:shape :box :color [0.78 0.80 0.84]}])}}})

(def ^:private up [0.0 1.0 0.0])

(defn- bark [gen]
  (let [t (min 1.0 (* 0.12 gen))]
    [(+ 0.30 (* 0.12 t)) (+ 0.22 (* 0.12 t)) (+ 0.15 (* 0.06 t))]))

(defn- leaf-color [rng]
  (let [t (double (rng))]
    [(+ 0.20 (* 0.18 t)) (+ 0.42 (* 0.18 t)) (+ 0.12 (* 0.05 t))]))

(defn- tree-rules
  "A branch is a unit length of limb along +y, from the origin, and the
  branches it carries at its tip, each a branch in its own smaller frame."
  [{:keys [generations ratio spread branching]
    :or {generations 9 ratio 0.74 spread 0.62 branching 3}}]
  (let [ratio (double ratio) spread (double spread)
        ;; The longest a child can be, and so the farthest the limbs of a
        ;; subtree reach past their root: a geometric series.
        longest (* 1.08 ratio)
        reach (/ 1.0 (- 1.0 longest))
        leaf-size 1.6]
    {:branch
     {:radius (+ reach leaf-size) :center v/zero
      :parts
      (fn [{:keys [seed gen width]}]
        (let [rng (random/rng seed)
              limb {:shape :cylinder :transform (scaling [width 1.0 width]) :color (bark gen)}]
          (if (>= gen generations)
            (cons limb
                  (for [k (range 4)]
                    {:shape :leaf :color (leaf-color rng)
                     :transform (compose (translation [0.0 (- 1.0 (* 0.18 k)) 0.0])
                                         (rotation up (+ (* 1.9 k) (random/uniform rng 0.6)))
                                         (rotation [1.0 0.0 0.0] (+ 0.5 (random/uniform rng 0.6)))
                                         (scaling leaf-size))}))
            (let [n (+ (dec (long branching)) (random/below rng 2))
                  phase (random/uniform rng (* 2.0 math/PI))]
              (cons limb
                    (for [i (range n)
                          :let [;; A fork of two splits evenly; a fork of three keeps
                                ;; a leader, nearly straight, between two sides.
                                leader? (and (zero? i) (>= n 3))
                                lean (if leader?
                                       (random/uniform rng (* 0.25 spread))
                                       (* spread (random/uniform rng 0.7 1.25)))
                                turn (+ phase (/ (* 2.0 math/PI i) n) (random/uniform rng -0.4 0.4))
                                k (* ratio (random/uniform rng 0.85 1.08) (if leader? 1.0 0.88))]]
                      {:instance :branch
                       :transform (compose (translation [0.0 1.0 0.0])
                                           (rotation up turn)
                                           (rotation [1.0 0.0 0.0] lean)
                                           (scaling k))
                       ;; The limb thins faster than it shortens: a child's
                       ;; frame is already shrunk by k.
                       :params {:seed (random/below rng 2147483646)
                                :gen (inc gen)
                                :width (* width (if leader? 0.95 0.8))}}))))))
      :proxy
      (fn [{:keys [seed gen width]}]
        ;; A subtree not worth expanding: its limb and a crown of foliage
        ;; about as big as the twigs it would have grown.
        (let [rng (random/rng seed)
              crown (* 0.3 reach)]
          [{:shape :cylinder :transform (scaling [width 1.0 width]) :color (bark gen)}
           {:shape :sphere :color (leaf-color rng)
            :transform (compose (translation [0.0 (+ 1.0 (* 0.6 crown)) 0.0])
                                (scaling [crown (* 0.8 crown) crown]))}]))}}))

(defn tree
  "A tree from one rule: a limb, and two to four smaller limbs at its tip,
  each with its own seed. Options `:seed`, `:generations`, `:ratio` (a
  child's length to its parent's), `:spread` (radians a side branch leans)
  and `:branching`."
  ([] (tree {}))
  ([{:keys [seed] :or {seed 1} :as opts}]
   {:root :branch
    :params {:seed seed :gen 0 :width 0.07}
    :rules (tree-rules opts)}))

(defn forest
  "Trees scattered over a disk of radius `:extent`, each a tree instance
  with a seed, a size and a turn of its own. Seen from inside it, this is
  where measuring by angle earns its keep: the near trees are expanded
  to the leaf, the far ones stop at a few crowns."
  ([] (forest {}))
  ([{:keys [seed trees extent] :or {seed 1 trees 60 extent 14.0} :as opts}]
   (let [rules (tree-rules opts)
         rng (random/rng seed 7 0 0)
         placed (vec (for [_ (range trees)]
                       (let [r (* extent (math/sqrt (double (rng))))
                             a (random/uniform rng (* 2.0 math/PI))]
                         {:at [(* r (math/cos a)) 0.0 (* r (math/sin a))]
                          :size (random/uniform rng 0.8 1.4)
                          :turn (random/uniform rng (* 2.0 math/PI))
                          :seed (random/below rng 2147483646)})))]
     {:root :forest
      :rules (assoc rules :forest
                    {:radius (+ extent (* 1.4 (:radius (:branch rules))))
                     :parts (fn [_]
                              (for [{:keys [at size turn seed]} placed]
                                {:instance :branch
                                 :transform (compose (translation at) (rotation up turn) (scaling size))
                                 :params {:seed seed :gen 0 :width 0.07}}))})})))

(def catalog
  "The models the demo offers, with a camera distance that frames each."
  [{:name "tree" :make tree :distance 8.0 :target [0.0 2.4 0.0]}
   {:name "forest" :make forest :distance 26.0 :target [0.0 2.0 0.0]}
   {:name "sierpinski" :make (fn [_] (sierpinski)) :distance 2.8 :target [0.0 0.0 0.0]}
   {:name "menger" :make (fn [_] (menger)) :distance 5.0 :target [0.0 0.0 0.0]}])
