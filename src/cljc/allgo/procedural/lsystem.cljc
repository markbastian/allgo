(ns allgo.procedural.lsystem
  "Lindenmayer systems: plants grown by rewriting a string, and drawn by a
  turtle that reads it.

  The other half of Hart's \"Procedural Synthesis of Geometry\" in Ebert
  et al., *Texturing & Modeling*; the notation is Prusinkiewicz and
  Lindenmayer's, *The Algorithmic Beauty of Plants*. Where
  `allgo.procedural.instancing` describes a plant as a rule that places
  copies of itself, an L-system describes it as a rule that *rewrites*
  itself: every symbol of the string is replaced, all at once, by its
  production, and after `n` rounds the string is the plant at age `n`.
  The two are the same idea from opposite ends -- a production with a
  bracketed copy of its own symbol is a recursive instance -- and the
  difference is when the work is done. An L-system writes the whole
  string out before any of it is drawn; instancing expands only what
  it is about to draw.

  ## Productions

  A map from a character to what replaces it: a string, or for a
  stochastic system a vector of `[weight string]` choices, one drawn at
  random for each occurrence. Characters without a production are copied.

  ## The turtle

  A position and an orthonormal frame -- heading `H`, left `L`, up `U` --
  read one symbol at a time:

  | Symbol  | Turtle                                        |
  |---------|-----------------------------------------------|
  | `F`     | move forward a step, drawing a limb           |
  | `f`     | move forward a step without drawing           |
  | `+` `-` | turn left, right, about `U`                   |
  | `&` `^` | pitch down, up, about `L`                     |
  | `\\` `/`| roll left, right, about `H`                   |
  | `\\|`   | turn around                                   |
  | `[` `]` | push and pop the turtle: a branch             |
  | `!`     | thin the limbs from here on                   |
  | `'`     | the next color in the palette                 |
  | `L`     | a leaf, here                                  |

  and any symbol listed in `:leaves` draws a leaf too, which is how an
  apex that was never rewritten -- a bud at the tip of the derivation --
  becomes foliage without a production of its own.

  After every step forward the heading is bent toward `:tropism`, by
  `:elasticity` times the sine of the angle between them: gravity for a
  weeping tree, the sun for one that leans into the light.

  Limbs and leaves come out as `allgo.procedural.instancing` primitives,
  so the same renderer draws both."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.procedural.instancing :as in]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.string :as str]))

(defn- choose [rng production]
  (if (string? production)
    production
    (let [total (reduce + (map first production))
          x (* total (double (rng)))]
      (loop [[[w s] & more] production acc 0.0]
        (let [acc (+ acc (double w))]
          (if (or (< x acc) (empty? more)) s (recur more acc)))))))

(defn derive-string
  "The string after `n` rounds of rewriting `axiom` by `rules`. Random
  choices are drawn from one stream of `seed`, in order along the string,
  so the same seed always grows the same plant. Stops early, at the last
  round that kept the string under `:max-length` symbols."
  ([system n] (derive-string system n {}))
  ([{:keys [axiom rules]} n {:keys [seed max-length] :or {seed 1 max-length 400000}}]
   (let [rng (random/rng seed 3 0 0)]
     (loop [s axiom i 0]
       (if (>= i n)
         s
         (let [next (str/join (map (fn [c] (if-let [p (get rules c)] (choose rng p) c)) s))]
           (if (> (count next) max-length) s (recur next (inc i)))))))))

(defn- rotate
  "`x` turned about the unit `axis` by `angle`, Rodrigues."
  [x axis angle]
  (let [c (math/cos angle) s (math/sin angle)]
    (v/add (v/scale x c)
           (v/scale (v/cross axis x) s)
           (v/scale axis (* (v/dot axis x) (- 1.0 c))))))

(defn- turn
  "The turtle's frame, turned about one of its own axes."
  [t axis angle]
  (let [k (get t axis)]
    (reduce (fn [t a] (if (= a axis) t (update t a rotate k angle))) t [:H :L :U])))

(defn- bend
  "The frame bent toward the tropism vector: turned about `H x T` by
  `e |H x T|`, so the pull is strongest across the vector and vanishes
  along it."
  [t tropism e]
  (let [a (v/cross (:H t) tropism)
        m (v/length a)]
    (if (< m 1e-9)
      t
      (let [axis (v/scale a (/ 1.0 m))
            angle (* e m)]
        (reduce #(update %1 %2 rotate axis angle) t [:H :L :U])))))

(defn interpret
  "The limbs and leaves the turtle draws reading `s`. Options:
  `:angle` (degrees, default 25), `:step` (1), `:width` (0.08),
  `:taper` (what `!` multiplies the width by, 0.75), `:leaf-size` (0.4),
  `:leaves` (symbols that draw a leaf, default `#{\\L}`), `:tropism`
  and `:elasticity`, `:jitter` (radians of random error in every turn),
  `:seed`, and `:palette`, the limb colors `'` steps through."
  [s {:keys [angle step width taper leaf-size leaves tropism elasticity jitter seed palette]
      :or {angle 25.0 step 1.0 width 0.08 taper 0.75 leaf-size 0.4 leaves #{\L}
           elasticity 0.0 jitter 0.0 seed 1
           palette [[0.36 0.25 0.16] [0.42 0.32 0.20] [0.30 0.45 0.20]]}}]
  (let [d (math/to-radians angle)
        rng (random/rng seed 5 0 0)
        wobble (fn [] (if (pos? jitter) (* jitter (- (* 2.0 (double (rng))) 1.0)) 0.0))
        start {:p v/zero :H [0.0 1.0 0.0] :L [-1.0 0.0 0.0] :U [0.0 0.0 1.0]
               :w (double width) :c 0}
        limb (fn [{:keys [p H L U w c]}]
               {:shape :cylinder :color (nth palette (mod c (count palette)))
                :transform (in/frame (v/scale L (- w)) (v/scale H step) (v/scale U w) p)})
        leaf (fn [{:keys [p H L U]}]
               (let [g (double (rng))]
                 {:shape :leaf :color [(+ 0.22 (* 0.15 g)) (+ 0.45 (* 0.15 g)) 0.14]
                  :transform (in/frame (v/scale L (- leaf-size)) (v/scale H leaf-size)
                                       (v/scale U leaf-size) p)}))]
    (loop [i 0 t start stack () out (transient [])]
      (if (>= i (count s))
        (persistent! out)
        (let [ch (nth s i)]
          (case ch
            \F (let [t' (cond-> (assoc t :p (v/add-scaled (:p t) (:H t) step))
                          tropism (bend tropism elasticity))]
                 (recur (inc i) t' stack (conj! out (limb t))))
            \f (recur (inc i) (assoc t :p (v/add-scaled (:p t) (:H t) step)) stack out)
            \+ (recur (inc i) (turn t :U (+ d (wobble))) stack out)
            \- (recur (inc i) (turn t :U (- (wobble) d)) stack out)
            \& (recur (inc i) (turn t :L (+ d (wobble))) stack out)
            \^ (recur (inc i) (turn t :L (- (wobble) d)) stack out)
            \\ (recur (inc i) (turn t :H (+ d (wobble))) stack out)
            \/ (recur (inc i) (turn t :H (- (wobble) d)) stack out)
            \| (recur (inc i) (turn t :U math/PI) stack out)
            \[ (recur (inc i) t (conj stack t) out)
            \] (recur (inc i) (first stack) (rest stack) out)
            \! (recur (inc i) (update t :w * taper) stack out)
            \' (recur (inc i) (update t :c inc) stack out)
            (recur (inc i) t stack (if (contains? leaves ch) (conj! out (leaf t)) out))))))))

;; ---------------------------------------------------------------------------
;; Plants

(def plants
  "Some systems, each with the rounds and turtle that suit it.

  - `shrub`: stochastic, so every seed is another bush of the same
    species. An apex `A` either forks three ways, forks two ways, or
    just grows, and the buds left at the end are leaves.
  - `conifer`: a leader that lays down a whorl of four branches and a
    length of trunk every round. The lowest whorl is the oldest and has
    grown longest, so the cone shape is not drawn anywhere: it is age.
  - `weeping`: forks that fan out and then sag, by tropism toward the
    ground alone.
  - `weed`: a plant in one plane, the plainest bracketed system: a stem
    that sprouts a smaller copy of itself to either side."
  [{:name "shrub"
    :system {:axiom "FA"
             :rules {\A [[0.45 "!F[&+A]/[&A]//[&-A]A"]
                         [0.35 "!F[&+A]///[&-A]"]
                         [0.20 "FA"]]
                     \F "F"}}
    :rounds 9
    :turtle {:angle 30 :step 0.42 :width 0.06 :taper 0.8 :leaves #{\A} :leaf-size 0.3
             :jitter 0.15 :tropism [0.0 1.0 0.0] :elasticity 0.12}
    :distance 8.0 :target [0.0 2.6 0.0]}
   {:name "conifer"
    ;; The whorl's `!!!!!!` starts each branch at half the trunk's width
    ;; there; every `!` after that thins it a little more.
    :system {:axiom "T"
             :rules {\T "F!/[&&&!!!!!!B]////[&&&!!!!!!B]////[&&&!!!!!!B]////[&&&!!!!!!B]//T"
                     \B "!FC[-L][+L][&L][^L]B"
                     \C "C"}}
    :rounds 10
    :turtle {:angle 24 :step 0.42 :width 0.12 :taper 0.89 :leaves #{\L \C} :leaf-size 0.26
             :jitter 0.05 :tropism [0.0 -1.0 0.0] :elasticity 0.08}
    :distance 10.0 :target [0.0 2.4 0.0]}
   {:name "weeping"
    :system {:axiom "FFFA"
             :rules {\A "!F[&&B]////[&&B]////[&&B]A"
                     \B "!F[-B]F[+B]FL"}}
    :rounds 6
    :turtle {:angle 28 :step 0.4 :width 0.1 :taper 0.82 :leaf-size 0.3
             :jitter 0.1 :tropism [0.0 -1.0 0.0] :elasticity 0.35}
    :distance 9.0 :target [0.0 2.5 0.0]}
   {:name "weed"
    :system {:axiom "X"
             :rules {\X "F[+X][-X]FXL"
                     \F "FF"}}
    :rounds 6
    :turtle {:angle 24 :step 0.05 :width 0.016 :taper 0.9 :leaf-size 0.1}
    :distance 8.0 :target [0.0 3.0 0.0]}])

(defn grow
  "The primitives of one of `plants`, at `rounds` (default its own)."
  [{:keys [system rounds turtle]} {:keys [seed] :as opts}]
  (let [n (get opts :rounds rounds)
        s (derive-string system n {:seed (or seed 1)})]
    {:string-length (count s)
     :primitives (interpret s (assoc turtle :seed (or seed 1)))}))
