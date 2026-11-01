(ns allgo.robotics.localization-test
  (:require [allgo.random :as random]
            [allgo.robotics.localization :as mcl]
            [clojure.test :refer [deftest is testing]]))

;; Monte Carlo localization, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 26.

(deftest raycast-test
  (testing "rays meet the walls where geometry says"
    (let [{:keys [walls max-range]} (mcl/room)]
      (is (< (abs (- 5.5 (mcl/raycast walls [1.5 4.5] 0.0 max-range))) 1e-9) "east to the pillar's west face, at x = 7")
      (is (< (abs (- 8.5 (mcl/raycast walls [1.5 5.0] 0.0 max-range))) 1e-9) "and above it, on to the wall at x = 10")
      (is (< (abs (- 1.5 (mcl/raycast walls [1.5 4.5] Math/PI max-range))) 1e-9) "west to x = 0")
      (is (< (abs (- 0.5 (mcl/raycast walls [7.5 3.0] (/ Math/PI 2) max-range))) 1e-9) "north to the pillar"))))

(deftest localization-test
  (testing "from not knowing where it is, the robot's particles find it: within a meter at the end of the route, every time, and usually within half"
    (let [errors (for [seed [5 6 7 8 9 10]]
                   (let [world (assoc (mcl/room) :rng (random/rng seed))
                         controls (concat (repeat 40 [1.0 0.0]) (repeat 16 [0.5 -1.0]) (repeat 15 [1.0 0.0]))
                         {:keys [truth particles]} (last (mcl/mcl-seq world [1.5 4.5 0.0] controls 1500))
                         [[x y]] (mcl/estimate particles)]
                     (Math/hypot (- x (truth 0)) (- y (truth 1)))))]
      (is (every? #(< % 1.0) errors))
      (is (< (nth (sort errors) 3) 0.5)))))
