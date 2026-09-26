(ns allgo.simulation.flock
  "One interface over the two flocking implementations, so a caller can
  hold either and a benchmark can run both.

  `allgo.simulation.boids` is the reference: a vector of `{:pos :vel}`
  maps, and the clearer statement of the rules.
  `allgo.simulation.boids-flat` is the same rules on flat arrays, several
  times faster.

  They agree to floating-point noise -- `boids-flat-test` runs them side
  by side and compares -- so which one a demo holds is a performance
  decision, not a behavioral one.

  This namespace knows only about the reference. The flat implementation
  implements `Flock` where its record is defined, which is what keeps the
  dependency pointing one way and means neither type has to be imported."
  (:require [allgo.simulation.boids :as boids])
  #?(:clj (:import [clojure.lang PersistentVector])))

(defprotocol Flock
  (advance [flock bounds params] "One tick, returning a flock of the same kind.")
  (as-boids [flock] "The flock as `[{:pos :vel} ...]`, whatever it is inside.")
  (flock-size [flock]))

(extend-protocol Flock
  PersistentVector
  (advance [flock bounds params] (boids/step flock bounds params))
  (as-boids [flock] flock)
  (flock-size [flock] (count flock)))

(defn reference
  "A reference flock of `n` boids, randomly placed."
  ([n bounds] (reference n bounds boids/defaults))
  ([n bounds params] (boids/flock n bounds params)))

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
