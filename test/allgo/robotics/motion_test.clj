(ns allgo.robotics.motion-test
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.random :as random]
            [allgo.robotics.motion :as motion]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Motion planning, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 26. The checks are exact answers
;; worked out here: the shortest way round a square.

(def ^:private square [[1.0 1.0] [3.0 1.0] [3.0 3.0] [1.0 3.0]])
(def ^:private room {:bounds [[-1.0 -1.0] [5.0 5.0]] :obstacles [square] :radius 0.0})
(def ^:private start [0.0 2.0])
(def ^:private goal [4.0 2.0])

(defn- shortest
  "Round the square grown by r, corners mitered: to a corner, along a
  side, to the goal -- 2 sqrt((1 - r)^2 + (1 + r)^2) + 2 + 2r."
  [r]
  (+ (* 2 (math/sqrt (+ (math/pow (- 1 r) 2) (math/pow (+ 1 r) 2)))) 2 (* 2 r)))

(deftest collision-test
  (let [w (assoc room :radius 0.2)]
    (testing "a disk is free where it clears the obstacles by its radius"
      (is (motion/free? w [0.0 2.0]))
      (is (not (motion/free? w [0.9 2.0])))
      (is (not (motion/free? w [2.0 2.0])))
      (is (not (motion/free? w [-0.9 2.0])) "nor too near the wall"))
    (testing "a straight move is free unless it passes within the radius"
      (is (motion/segment-free? w [0.0 0.5] [4.0 0.5]))
      (is (not (motion/segment-free? w [0.0 0.9] [4.0 0.9])))
      (is (not (motion/segment-free? w start goal))))))

(deftest visibility-test
  (testing "the visibility graph gives the shortest path, for a point and for a disk"
    (doseq [r [0.0 0.2 0.5]]
      (let [p (motion/visibility-path (assoc room :radius r) start goal)]
        (is (< (abs (- (shortest r) (motion/path-length p))) 1e-9) r)
        ;; it grazes the grown square, so it clears a slightly smaller robot
        (when (pos? r) (is (motion/path-free? (assoc room :radius (* 0.999 r)) p) r))))))

;; The walls of `room`, each its own obstacle to the Voronoi graph.
(def ^:private walls [[[-1.0 -1.0] [5.0 -1.0]] [[5.0 -1.0] [5.0 5.0]] [[5.0 5.0] [-1.0 5.0]] [[-1.0 5.0] [-1.0 -1.0]]])

(defn- distances
  "The distance from `p` to the square and to each wall, nearest first."
  [p]
  (sort (cons (apply min (map (fn [[a b]] (polygon/segment-distance p a b)) (polygon/edges square)))
              (map (fn [[a b]] (polygon/segment-distance p a b)) walls))))

(defn- least-clearance [w path]
  (apply min (for [[a b] (partition 2 1 path) t (range 0.0 1.0001 0.05)] (motion/clearance w (polygon/lerp a b t)))))

(deftest voronoi-test
  (let [w (assoc room :radius 0.2)
        h 0.2
        g (motion/voronoi-graph w h)]
    (testing "every node of the graph is about equally far from its two nearest obstacles"
      (is (seq (:nodes g)))
      (doseq [n (:nodes g)]
        (let [[d1 d2] (distances n)]
          (is (< (- d2 d1) h) n))))
    (testing "the chapter's plan: onto the graph, along it, off it -- free, and well clear"
      (let [p (motion/voronoi-path w g start goal)
            v (motion/visibility-path w start goal)]
        (is (motion/path-free? w p))
        (is (= [start goal] [(first p) (peek p)]))
        (is (> (least-clearance w (subvec p 1 (dec (count p)))) 0.9) "the graph's part a meter clear")
        (is (< (least-clearance w v) 0.21) "the shortest path grazes the square")
        (is (> (motion/path-length p) (motion/path-length v)) "and is the longer for it")))
    (testing "no path when the goal cannot reach the graph"
      (is (nil? (motion/voronoi-path w g start [2.0 2.0]))))))

(deftest planners-test
  (let [w (assoc room :radius 0.2)
        best (shortest 0.2)]
    (testing "every planner's path is collision free, and none is shorter than the shortest"
      (let [cell (motion/cell-path w start goal 0.1)
            prm (motion/prm-path w (nth (motion/prm-seq w (random/rng 2)) 150) start goal)
            rrt (motion/tree-path (first (filter :goal (motion/rrt-seq w start goal (random/rng 1)))))]
        (doseq [[name p] {:cell cell :prm prm :rrt rrt}]
          (is (motion/path-free? w p) name)
          (is (>= (motion/path-length p) (- best 1e-9)) name)
          (is (= [start goal] [(first p) (peek p)]) name))))
    (testing "RRT* keeps shortening its path toward the shortest as it grows"
      (let [lengths (map #(motion/path-length (motion/tree-path (nth (motion/rrt-seq w start goal (random/rng 1) :star? true) %)))
                         [200 500 1000 2000])]
        (is (apply >= lengths))
        (is (< (last lengths) (* 1.01 best)))))
    (testing "shortcutting never lengthens a path or makes it collide"
      (let [p (motion/tree-path (first (filter :goal (motion/rrt-seq w start goal (random/rng 4)))))
            s (nth (motion/shortcut-seq w p (random/rng 5)) 100)]
        (is (motion/path-free? w s))
        (is (<= (motion/path-length s) (motion/path-length p)))))
    (testing "a passage narrower than a cell is lost to the grid but found by sampling"
      (let [gap {:bounds [[0.0 0.0] [6.0 4.0]] :radius 0.1
                 :obstacles [[[2.9 0.0] [3.1 0.0] [3.1 1.75] [2.9 1.75]] [[2.9 2.25] [3.1 2.25] [3.1 4.0] [2.9 4.0]]]}]
        (is (nil? (motion/cell-path gap [1.0 2.0] [5.0 2.0] 0.5)))
        (is (some? (motion/visibility-path gap [1.0 2.0] [5.0 2.0])))))))

(deftest potential-field-test
  (testing "the potential field reaches the goal past a post, and is trapped in a cup -- a local minimum"
    (let [post {:bounds [[-1.0 -1.0] [6.0 6.0]] :obstacles [[[2.0 1.0] [2.5 1.0] [2.5 1.5] [2.0 1.5]]] :radius 0.0}
          cup (assoc post :obstacles [[[2.0 1.0] [3.0 1.0] [3.0 4.0] [2.0 4.0] [2.0 3.5] [2.5 3.5] [2.5 1.5] [2.0 1.5]]])
          run (fn [w] (nth (motion/potential-seq w [0.5 2.5] [5.0 2.5] :influence 0.5 :repulse 0.05) 1000))]
      (is (= [5.0 2.5] (run post)))
      (let [[x y] (run cup)]
        (is (< 2.0 x 2.5))
        (is (< 1.5 y 3.5))))))

(deftest arm-test
  (let [arm (motion/arm [0.0 0.0] [2.0 1.5])
        w {:obstacles [[[2.0 1.0] [3.0 1.0] [3.0 2.0] [2.0 2.0]]]}]
    (testing "forward kinematics through the chain: shoulder up, elbow bent right"
      (let [[[s e] [_ h]] (motion/arm-links arm [(/ math/PI 2) (- (/ math/PI 2))])]
        (is (every? #(< (abs %) 1e-12) (map - (concat s e h) [0.0 0.0 0.0 2.0 1.5 2.0])))))
    (testing "configuration space marks exactly the configurations that collide"
      (let [n 32 grid (motion/arm-cspace arm w n)]
        (doseq [i (range n) j (range n)]
          (is (= (get-in grid [i j]) (motion/arm-collides? arm w (motion/arm-angles n [i j])))))
        (is (some true? (flatten grid)))
        (is (some false? (flatten grid)))))
    (testing "a plan through configuration space never collides"
      (let [grid (motion/arm-cspace arm w 64)
            plan (motion/arm-plan grid [-0.5 0.0] [1.4 0.3])]
        (is (seq plan))
        (is (not-any? #(motion/arm-collides? arm w %) plan))))))
