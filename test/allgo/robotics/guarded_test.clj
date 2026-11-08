(ns allgo.robotics.guarded-test
  (:require [allgo.random :as random]
            [allgo.robotics.guarded :as guarded]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Guarded moves, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, section 26.6: a peg into a hole when
;; every motion may stray 20 degrees from the command.

(def ^:private world {:width 0.4 :depth 0.6})
(def ^:private cone (math/to-radians 20))
(def ^:private start [0.0 2.0])

(deftest guarded-moves-test
  (let [rng (random/rng 1)
        down (guarded/direction (- (/ math/PI 2)))
        ;; 30 degrees left of straight down: the whole cone misses the hole
        left (guarded/direction (- (+ (/ math/PI 2) (math/to-radians 30))))
        ;; 45 degrees below the horizontal, to the right: every direction in
        ;; the cone moves right and down
        right (guarded/direction (math/to-radians -45))
        bottom (fn [p _] (guarded/at-bottom? world p))
        trials 200]
    (testing "straight down, stopping at contact, the peg lands in the hole only now and then"
      (let [ends (repeatedly trials #(peek (guarded/guarded-moves world start [[down guarded/on-contact]] cone 0.01 rng)))]
        (is (< 20 (count (filter #(guarded/in-hole? world %) ends)) 100))))
    (testing "the first guarded move always ends on the floor left of the hole, wherever exactly"
      (let [ends (repeatedly trials #(peek (guarded/guarded-moves world start [[left guarded/on-contact]] cone 0.01 rng)))]
        (is (every? (fn [[x y]] (and (< x (- (/ (:width world) 2))) (< (abs y) 0.02))) ends))))
    (testing "and the second slides it right along the floor, into the hole and down the far wall to the bottom -- every time"
      (let [ends (repeatedly trials #(peek (guarded/guarded-moves world start [[left guarded/on-contact] [right bottom]] cone 0.01 rng)))]
        (is (every? #(guarded/at-bottom? world %) ends))))))
