(ns allgo.robotics.guarded
  "Guarded moves: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, section 26.6 -- reaching a goal despite
  uncertainty in the motion itself, with no map to plan in more finely.

  A peg is to go into a hole, seen in its two-dimensional configuration
  space: a floor at y = 0 with a hole `width` wide and `depth` deep
  centered on x = 0. A motion command is a direction to move in, but the
  robot's actual direction is anywhere within `cone` radians of it. A
  guarded move is a command and a termination condition, a test on what
  the robot senses -- here, contact. Against a surface the robot slides
  along it (compliant motion), so long as its direction has any part
  along the surface, and stops where it has none.

  Moving straight down from above the hole, the robot lands somewhere in
  the cone's footprint, more likely beside the hole than in it, and does
  not know which side. The section's strategy: first move down and to
  the left, steeply enough that the whole cone misses the hole, and stop
  on contact -- now the robot is on the floor left of the hole, wherever
  exactly; then move down and to the right, every direction of the cone
  sliding it right along the floor until it drops into the hole, slides
  down the far wall, and stops at the bottom."
  (:require [clojure.math :as math]))

(defn blocked?
  "Whether point `p` is inside the solid: below the floor, except in the
  hole."
  [{:keys [width depth]} [x y]]
  (and (< y 0.0) (or (> (abs x) (/ width 2.0)) (< y (- depth)))))

(defn at-bottom?
  "Whether `p` rests on the hole's bottom, to within `tolerance` (default
  0.02) -- a step stops short of a surface by less than its length."
  ([world p] (at-bottom? world p 0.02))
  ([{:keys [width depth]} [x y] tolerance]
   (and (<= (abs x) (/ width 2.0)) (<= y (+ (- depth) tolerance)))))

(defn in-hole?
  "Whether `p` is down in the hole."
  [{:keys [width]} [x y]]
  (and (<= (abs x) (/ width 2.0)) (< y 0.0)))

(defn move-seq
  "The robot's positions under the motion command `command` (a unit
  direction) from `p`, its actual direction `command` turned by `error`
  radians, a step of `ds` a step: a step that would enter the solid slides
  along the surface -- along x on a floor, along y on a wall -- or, if
  neither is possible, stays put. `stop?` (fn [p contact?]) ends the
  motion: the positions up to and including where it stops."
  [world p [cx cy] error ds stop?]
  (let [c (math/cos error) s (math/sin error)
        [dx dy] [(- (* c cx) (* s cy)) (+ (* s cx) (* c cy))]
        step (fn [[x y]]
               (let [full [(+ x (* ds dx)) (+ y (* ds dy))]
                     along-x [(+ x (* ds dx)) y]
                     along-y [x (+ y (* ds dy))]]
                 (cond (not (blocked? world full)) [full false]
                       (and (> (abs dx) 1e-9) (not (blocked? world along-x))) [along-x true]
                       (and (> (abs dy) 1e-9) (not (blocked? world along-y))) [along-y true]
                       :else [[x y] true])))]
    (loop [p p out [p] n 0]
      (let [[q contact?] (step p)]
        (if (or (stop? q contact?) (= q p) (> n 100000))
          (conj out q)
          (recur q (conj out q) (inc n)))))))

(defn on-contact
  "The termination condition of contact with any surface."
  [_ contact?]
  contact?)

(defn guarded-moves
  "A sequence of guarded moves, each `[command stop?]`, run from `p` with
  actual directions `command` turned by an error drawn uniformly within
  the `cone` for each move: the positions, all moves in turn."
  [world p moves cone ds rng]
  (reduce (fn [path [command stop?]]
            (let [error (* cone (- (* 2.0 (rng)) 1.0))]
              (into path (rest (move-seq world (peek path) command error ds stop?)))))
          [p] moves))

(defn direction
  "The unit direction at `angle` radians counterclockwise from +x."
  [angle]
  [(math/cos angle) (math/sin angle)])
