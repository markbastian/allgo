(ns allgo.simulation.flock
  "One interface over the two flocking implementations, so a caller can
  hold either and a benchmark can run both.

  `allgo.simulation.boids` is the reference: a vector of `{:pos :vel}`
  maps, dimension-generic, and the only one that steers around obstacles.
  `allgo.simulation.boids-flat` is the same three rules on flat arrays,
  several times faster and limited to flocking.

  They agree to floating-point noise -- `boids-flat-test` runs them side
  by side and compares -- so which one a demo holds is a performance
  decision, not a behavioural one."
  (:require [allgo.simulation.boids :as boids]
            [allgo.simulation.boids-flat :as flat #?@(:cljs [:refer [FlatFlock]])])
  #?(:clj (:import [allgo.simulation.boids_flat FlatFlock])))

(defprotocol Flock
  (advance [flock bounds params] "One tick, returning a flock of the same kind.")
  (as-boids [flock] "The flock as `[{:pos :vel} ...]`, whatever it is inside.")
  (flock-size [flock]))

(extend-protocol Flock
  #?(:clj clojure.lang.PersistentVector :cljs cljs.core/PersistentVector)
  (advance [flock bounds params] (boids/step flock bounds params))
  (as-boids [flock] flock)
  (flock-size [flock] (count flock)))

(extend-type FlatFlock
  Flock
  (advance [flock bounds params] (flat/step flock bounds params))
  (as-boids [flock] (flat/to-boids flock))
  (flock-size [flock] (:n flock)))

(defn reference
  "A reference flock of `n` boids, randomly placed."
  ([n bounds] (reference n bounds boids/defaults))
  ([n bounds params] (boids/flock n bounds params)))

(defn fast
  "A flat-array flock holding the same state as `flock`."
  [flock]
  (flat/from-boids (as-boids flock)))

(defn simulate
  "Run `ticks` steps of whichever implementation is passed."
  [flock bounds params ticks]
  (loop [flock flock i 0]
    (if (= i ticks) flock (recur (advance flock bounds params) (inc i)))))

(defn divergence
  "The largest difference in any coordinate between two flocks.

  What makes the two implementations comparable rather than merely
  similar: run both from the same start and this should stay at the level
  of floating-point noise."
  [a b]
  (let [pairs (map vector (as-boids a) (as-boids b))]
    (reduce max 0.0
            (for [[x y] pairs
                  k     [:pos :vel]
                  [u v] (map vector (k x) (k y))]
              (abs (- (double u) (double v)))))))
