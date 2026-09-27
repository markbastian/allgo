(ns allgo.random
  "Repeatable randomness, the same on the JVM and in a browser.

  `rand` cannot be seeded, and `java.util.Random` is the JVM's alone, so
  nothing that drew from either could be reproduced: not a dungeon, not a
  cave, not an arm's restarts. This is the one generator everything draws
  from instead. Two layers:

  * **Stateless** -- `advance`, `unit`, `cell-seed`. A state is a double,
    and the next one is a function of it. This is what the noise bases
    use, because a lattice cell's randomness has to come from its
    coordinates and not from how many cells were asked for before it.
  * **Streams** -- `(rng seed)` is a function of no arguments returning
    the next number in [0, 1), which is exactly the shape `rand` has. So
    anything that takes an `:rng` takes `rand` for a different result
    every time or `(rng 42)` for the same one; and `pick`, `shuffle`,
    `sample` and `gaussian` are `rand-nth`, `shuffle`, `random-sample`
    and a normal draw that take one."
  (:refer-clojure :exclude [shuffle])
  (:require [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; The generator
;;
;; Park-Miller, done in doubles rather than in bit twiddling, so that the
;; JVM's 64-bit longs and JavaScript's 32-bit bitwise operators cannot
;; disagree about what a seed means. Every product below stays under 2^53
;; and so is exact on both.

(def ^:private prng-modulus 2147483647.0)
(def ^:private prng-multiplier 16807.0)

(defn advance
  "The next state of the random stream. Park-Miller's minimal standard.

  The modulus is taken by hand rather than with `rem`, and that is worth a
  paragraph because it is the innermost loop of every cellular and
  sparse-convolution sample. A remainder of two large doubles misses the
  integer fast path in a JavaScript engine and costs on the order of a
  hundred nanoseconds; subtracting `m * floor(x / m)` costs a divide and a
  floor, both of which are single instructions, and measured eight times
  faster here.

  It is also exact, which is the part that has to be argued rather than
  measured. `16807 * s` never exceeds 2^53, so the product is an exact
  integer. Its true quotient by `m` is either an integer -- in which case
  the division is exact -- or at least `1/m` away from one, which is four
  hundred times the largest error the division can introduce. So the floor
  is never off by one, and neither platform can disagree with the other
  about what a seed means."
  ^double [^double s]
  (let [x (* prng-multiplier s)
        v (- x (* prng-modulus (math/floor (/ x prng-modulus))))]
    (if (pos? v) v 1.0)))

(defn unit
  "A random stream state as a number in [0, 1)."
  ^double [^double s]
  (/ s prng-modulus))

(defn signed
  "A random stream state as a number in [-1, 1)."
  ^double [^double s]
  (- (* 2.0 (unit s)) 1.0))

(defn- fold
  "Mixes one integer into a stream state.

  Coordinates are taken modulo 65536, so the lattice repeats every 65536
  cells along each axis. That is no worse than the 256-entry permutation
  table underneath the classic noise functions and far past anything a
  scene will reach."
  ^double [^double s ^long v]
  ;; `bit-and` with a mask is `mod` by a power of two, on both platforms
  ;; and for negative coordinates too, and it is a single instruction
  ;; where `mod` is a call.
  (advance (advance (+ s (double (bit-and v 0xFFFF)) 1.0))))

(defn cell-seed
  "A repeatable random stream for the lattice cell `[i j k]` of world `seed`.

  This is what lets a procedural texture be unbounded and still be a
  function: no cell is stored anywhere, and asking twice gives the same
  answer because the answer is derived from the coordinates."
  ^double [^long seed ^long i ^long j ^long k]
  (-> (+ 1.0 (double (mod seed 2147483646)))
      (fold i) (fold j) (fold k)))

;; ---------------------------------------------------------------------------
;; Streams

(defn rng
  "A seeded stream: a function of no arguments returning the next number
  in [0, 1), a drop-in for `rand`. The same keys give the same numbers, on
  either platform. Extra integer keys pick one of many independent streams
  under one seed."
  ([seed] (rng seed 0 0 0))
  ([seed i j k]
   (let [s (volatile! (cell-seed (long seed) (long i) (long j) (long k)))]
     (fn [] (let [x @s] (vswap! s advance) (unit x))))))

(defn uniform
  "A number drawn uniformly from [`lo`, `hi`), or [0, `hi`)."
  (^double [rng hi] (* (double (rng)) (double hi)))
  (^double [rng lo hi] (+ (double lo) (* (double (rng)) (- (double hi) (double lo))))))

(defn below
  "An integer drawn uniformly from [0, `n`): `rand-int`."
  ^long [rng n]
  (min (dec (long n)) (long (* (double (rng)) (double n)))))

(defn pick
  "An element of `coll` chosen uniformly: `rand-nth`."
  [rng coll]
  (nth coll (below rng (count coll))))

(defn shuffle
  "The elements of `coll` in a uniformly random order, as a vector.
  Fisher-Yates, which draws each position's element from those not yet
  placed."
  [rng coll]
  (let [v (transient (vec coll))]
    (loop [i (dec (count v)) v v]
      (if (pos? i)
        (let [j (below rng (inc i))
              vi (nth v i)]
          (recur (dec i) (-> v (assoc! i (nth v j)) (assoc! j vi))))
        (persistent! v)))))

(defn sample
  "Each element of `coll` kept with probability `p`: `random-sample`."
  [rng p coll]
  (filter (fn [_] (< (double (rng)) (double p))) coll))

(defn gaussian
  "One draw from a normal distribution, by the Box-Muller transform."
  ^double [rng ^double mean ^double sd]
  (let [u1 (max 1e-12 (double (rng)))
        u2 (double (rng))]
    (+ mean (* sd (math/sqrt (* -2.0 (math/log u1))) (math/cos (* 2.0 math/PI u2))))))
