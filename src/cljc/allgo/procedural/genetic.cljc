(ns allgo.procedural.genetic
  "Genetic textures: textures nobody designed, bred by someone choosing
  which ones they like.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, chapter 19
  (Musgrave, \"Genetic Textures\"), after Karl Sims' \"Artificial
  Evolution for Computer Graphics\" (SIGGRAPH 1991). A texture here is an
  expression -- a tree of arithmetic and noise over the coordinates of
  the point being shaded -- and a palette its value is looked up in. Its
  genome *is* the expression: there is no fitness function, because
  nobody can write down what makes a texture good. A person looks at a
  screenful, picks the ones they like, and the next screenful is bred
  from those: mutants of a single parent, children of two.

  What evolution buys is a way to search a space too large and too
  strange to search by hand. Nobody would write
  `(cells (sin (* y 2.6)) (turb x z y) x)`; but nobody needs to,
  because it can be found by picking the best of twelve a dozen times.

  ## Expressions

  Plain data, so that they can be printed, read back, mutated and
  crossed:

  - a number, or one of `:x` `:y` `:z`;
  - `[op & args]`, with `op` from `ops`: arithmetic, `sin`, `abs`, a
    blend, and noise, fBm, turbulence and two cellular bases, each of
    which takes its *point* as three expressions. That is where the
    variety comes from: a noise whose point is another noise is domain
    distortion, discovered rather than designed.

  The value is squashed into the palette by `0.5 + 0.5 tanh v`, which
  keeps every expression drawable however large it gets. The noise bases
  have fixed seeds, so a genome alone decides its texture.

  ## Breeding

  Every random choice is drawn from a stream the caller passes in, so a
  generation is repeatable from its seed.

  - `mutate`: each node, with probability `rate`, is replaced by a new
    random subtree, nudged (a constant), swapped for another operator of
    the same arity, wrapped in a new operator, or collapsed into one of
    its arguments; each palette color may drift.
  - `cross`: a random subtree of one parent is replaced by a random
    subtree of the other, and the palette is taken color by color from
    either.

  Trees are kept under `max-nodes`: past that, breeding mostly produces
  slower textures, not better ones."
  (:require [allgo.math :as am]
            [allgo.procedural.cellular :as cellular]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.texture :as tx]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.string :as str]))

(def max-nodes 40)

;; ---------------------------------------------------------------------------
;; The operators

(def ^:private gradient (noise/gradient-basis {:seed 101}))
(def ^:private fbm (fractal/fbm (noise/gradient-basis {:seed 102}) {:octaves 5}))
(def ^:private turb (fractal/turbulence (noise/gradient-basis {:seed 103}) {:octaves 5}))
(def ^:private f1 (cellular/f1 {:seed 104}))
(def ^:private crack (cellular/f2-f1 {:seed 105}))

(defn- point-op
  "An operator whose three arguments are a point to sample `basis` at."
  [basis]
  (fn [fx fy fz]
    (fn ^double [x y z]
      (basis (double (fx x y z)) (double (fy x y z)) (double (fz x y z))))))

(def ops
  "Operator -> `{:arity n :f (fn [& arg-fns] -> fn of x y z) :weight w}`.
  `weight` is how often a random tree picks it: noise is weighted up,
  being what makes a texture rather than a gradient."
  {:add {:arity 2 :weight 1.0 :f (fn [a b] (fn ^double [x y z] (+ (double (a x y z)) (double (b x y z)))))}
   :sub {:arity 2 :weight 0.8 :f (fn [a b] (fn ^double [x y z] (- (double (a x y z)) (double (b x y z)))))}
   :mul {:arity 2 :weight 1.0 :f (fn [a b] (fn ^double [x y z] (* (double (a x y z)) (double (b x y z)))))}
   :max {:arity 2 :weight 0.5 :f (fn [a b] (fn ^double [x y z] (max (double (a x y z)) (double (b x y z)))))}
   :min {:arity 2 :weight 0.5 :f (fn [a b] (fn ^double [x y z] (min (double (a x y z)) (double (b x y z)))))}
   :sin {:arity 1 :weight 1.0 :f (fn [a] (fn ^double [x y z] (math/sin (* 2.0 (double (a x y z))))))}
   :abs {:arity 1 :weight 0.6 :f (fn [a] (fn ^double [x y z] (abs (double (a x y z)))))}
   :blend {:arity 3 :weight 0.8
           :f (fn [a b t]
                (fn ^double [x y z]
                  (let [s (am/clamp (+ 0.5 (* 0.5 (double (t x y z)))) 0.0 1.0)]
                    (+ (* (- 1.0 s) (double (a x y z))) (* s (double (b x y z)))))))}
   :noise {:arity 3 :weight 1.6 :f (fn [a b c] ((point-op gradient) a b c))}
   :fbm {:arity 3 :weight 1.4 :f (fn [a b c] ((point-op fbm) a b c))}
   :turb {:arity 3 :weight 1.2 :f (fn [a b c] ((point-op turb) a b c))}
   :cells {:arity 3 :weight 1.0 :f (fn [a b c] ((point-op f1) a b c))}
   :crack {:arity 3 :weight 0.8 :f (fn [a b c] ((point-op crack) a b c))}})

(def ^:private op-names (vec (sort (keys ops))))

(defn- arity [op] (:arity (ops op)))

;; ---------------------------------------------------------------------------
;; Expressions as data

(defn compile-expr
  "`expr` as a function of `x y z`."
  [expr]
  (cond
    (number? expr) (let [c (double expr)] (fn ^double [_ _ _] c))
    (= :x expr) (fn ^double [x _ _] (double x))
    (= :y expr) (fn ^double [_ y _] (double y))
    (= :z expr) (fn ^double [_ _ z] (double z))
    :else (apply (:f (ops (first expr))) (map compile-expr (rest expr)))))

(defn node-count ^long [expr]
  (if (vector? expr) (reduce + 1 (map node-count (rest expr))) 1))

(defn- paths
  "The path to every node: the indices to `get-in` by, `[]` the root."
  [expr]
  (if (vector? expr)
    (into [[]] (mapcat (fn [i] (map #(into [i] %) (paths (nth expr i))))
                       (range 1 (count expr))))
    [[]]))

(defn- put [expr path sub] (if (empty? path) sub (assoc-in expr path sub)))

(defn- fmt-num [x]
  (let [r (/ (math/round (* 100.0 (double x))) 100.0)]
    (str r)))

(defn expr->string
  "As a Lisp would print it, with constants to two places."
  [expr]
  (cond
    (number? expr) (fmt-num expr)
    (keyword? expr) (name expr)
    :else (str "(" (str/join " " (cons (name (first expr)) (map expr->string (rest expr)))) ")")))

;; ---------------------------------------------------------------------------
;; Random genomes

(defn- weighted-op [rng]
  (let [total (reduce + (map #(:weight (ops %)) op-names))
        x (* total (double (rng)))]
    (loop [[op & more] op-names acc 0.0]
      (let [acc (+ acc (double (:weight (ops op))))]
        (if (or (< x acc) (empty? more)) op (recur more acc))))))

(defn- random-terminal [rng]
  (if (< (double (rng)) 0.7)
    (random/pick rng [:x :y :z])
    (/ (math/round (* 100.0 (random/uniform rng -3.0 3.0))) 100.0)))

(defn random-expr
  "A random tree no deeper than `depth`."
  [rng depth]
  (if (or (<= depth 0) (< (double (rng)) 0.2))
    (random-terminal rng)
    (let [op (weighted-op rng)]
      (into [op] (repeatedly (arity op) #(random-expr rng (dec depth)))))))

(defn- random-color [rng]
  (vec (repeatedly 3 #(double (rng)))))

(defn random-genome [rng]
  {:expr (random-expr rng 4)
   :palette (vec (repeatedly 5 #(random-color rng)))})

(defn- squash ^double [^double v] (+ 0.5 (* 0.5 (math/tanh v))))

(defn texture
  "`(fn [x y z] -> [r g b])` for `genome`."
  [{:keys [expr palette]}]
  (let [f (compile-expr expr)
        ramp (tx/color-map palette)]
    (fn [x y z] (ramp (squash (double (f x y z)))))))

(defn interesting?
  "Whether a genome draws something other than a flat color: its value
  varies over a coarse grid of the tile. A good part of random trees do
  not -- `(max 2.4 (abs x))` is nearly flat in a palette -- and a screen
  of them is a screen wasted."
  [{:keys [expr]}]
  (let [f (compile-expr expr)
        vs (for [i (range 5) j (range 5)]
             (squash (double (f (- (* 0.5 i) 1.0) (- (* 0.5 j) 1.0) 0.3))))
        mean (/ (reduce + vs) (count vs))
        var (/ (reduce + (map #(let [d (- % mean)] (* d d)) vs)) (count vs))]
    (> var 0.004)))

(defn- until-interesting [make]
  (loop [tries 0]
    (let [g (make)]
      (if (or (and (interesting? g) (<= (node-count (:expr g)) max-nodes)) (> tries 30))
        g
        (recur (inc tries))))))

(defn population
  "`n` random, drawable genomes."
  [rng n]
  (vec (repeatedly n #(until-interesting
                       (fn []
                          ;; A bare coordinate is a ramp: drawable, and dull.
                         (loop [] (let [g (random-genome rng)]
                                    (if (< (node-count (:expr g)) 5) (recur) g))))))))

;; ---------------------------------------------------------------------------
;; Breeding

(defn- mutate-node [rng node]
  (let [choice (random/below rng 5)]
    (cond
      (number? node)
      (if (< choice 3)
        (+ node (* 0.25 (- (* 2.0 (double (rng))) 1.0) (max 1.0 (abs node))))
        (random-expr rng 2))

      (keyword? node)
      (if (< choice 2) (random-expr rng 2) (random-terminal rng))

      :else
      (case (int choice)
        0 (random-expr rng 3)
        1 (let [same (filterv #(= (arity %) (arity (first node))) op-names)]
            (assoc node 0 (random/pick rng same)))
        2 (let [op (weighted-op rng)
                at (random/below rng (arity op))]
            (into [op] (for [i (range (arity op))]
                         (if (= i at) node (random-expr rng 1)))))
        (nth node (inc (random/below rng (dec (count node)))))))))

(defn- mutate-expr [rng expr rate]
  (let [expr (if (< (double (rng)) rate) (mutate-node rng expr) expr)]
    (if (vector? expr)
      (into [(first expr)] (map #(mutate-expr rng % rate) (rest expr)))
      expr)))

(defn- drift [rng color amount]
  (mapv #(am/clamp (+ % (* amount (- (* 2.0 (double (rng))) 1.0))) 0.0 1.0) color))

(defn mutate
  "A child of one parent. `rate` is the chance that any one node, or any
  one palette color, changes."
  [rng {:keys [expr palette]} rate]
  (until-interesting
   (fn []
     (let [e (mutate-expr rng expr rate)
           ;; One change at least: a child identical to its parent is a
           ;; wasted tile.
           e (if (= e expr)
               (let [p (random/pick rng (paths expr))]
                 (put expr p (mutate-node rng (if (empty? p) expr (get-in expr p)))))
               e)]
       {:expr e
        :palette (mapv #(if (< (double (rng)) rate) (drift rng % 0.35) %) palette)}))))

(defn cross
  "A child of two parents: `a` with one of its subtrees replaced by one
  of `b`'s, and each palette color from one or the other."
  [rng a b]
  (until-interesting
   (fn []
     (let [pa (random/pick rng (paths (:expr a)))
           pb (random/pick rng (paths (:expr b)))
           sub (if (empty? pb) (:expr b) (get-in (:expr b) pb))]
       {:expr (put (:expr a) pa sub)
        :palette (mapv #(if (< (double (rng)) 0.5) %1 %2) (:palette a) (:palette b))}))))

(defn breed
  "The next generation of `n`, from the genomes `parents` chosen out of
  the last. None chosen is a fresh start; one, its mutants; more, their
  children, mutated a little. The parents themselves are kept, first, so
  that a choice is never lost by breeding from it."
  [rng parents n rate]
  (let [parents (vec parents)]
    (case (count parents)
      0 (population rng n)
      1 (into parents (repeatedly (dec n) #(mutate rng (first parents) rate)))
      (into parents
            (repeatedly (- n (count parents))
                        #(let [[a b] (take 2 (random/shuffle rng parents))]
                           (mutate rng (cross rng a b) (* 0.5 rate))))))))
