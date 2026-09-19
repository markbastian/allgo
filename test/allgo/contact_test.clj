(ns allgo.contact-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.contact :as ct]
            [allgo.physics.rigid :as rigid]
            [clojure.test :refer [deftest is testing]]))

(defn- box
  ([pos] (box pos [1.0 1.0 1.0] q/identity-q))
  ([pos size] (box pos size q/identity-q))
  ([pos size rot] (rigid/box {:pos pos :size size :rot rot :density 1.0})))

(defn- ball [pos r] (rigid/ball {:pos pos :radius r :density 1.0}))

(defn- unit? [n] (< (abs (- 1.0 (v/length n))) 1e-9))

(deftest face-contact-test
  (testing "a box resting squarely on a box gives a whole face, not a point"
    ;; The reason this namespace exists. GJK with EPA answers the same
    ;; question with one point, and a brick supported at one point rocks.
    (let [cs (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.0 0.95 0.0]))]
      (is (= 4 (count cs)))
      (is (every? #(< (v/distance (:normal %) [0.0 1.0 0.0]) 1e-9) cs))
      (is (every? #(< (abs (- 0.05 (:depth %))) 1e-9) cs))
      (testing "and the four points are spread out, not stacked"
        (let [pts (map :point cs)]
          (is (> (count (distinct pts)) 3))))))

  (testing "a box hanging over the edge is clipped to the overlap"
    ;; Sutherland-Hodgman against the reference face's sides. Without the
    ;; clip the contact points sit out in mid-air past the supporting
    ;; brick, and a wall leans on nothing.
    (let [cs (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.6 0.95 0.0]))]
      (is (seq cs))
      (is (every? (fn [{[px _ pz] :point}]
                    (and (<= -0.501 px 0.501) (<= -0.501 pz 0.501)))
                  cs)))))

(deftest separation-test
  (testing "boxes with daylight between them do not touch"
    (is (empty? (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.0 3.0 0.0]))))
    (is (empty? (ct/between 0 1 (box [0.0 0.0 0.0]) (box [2.5 0.0 0.0])))))

  (testing "boxes exactly abutting do not overlap either"
    ;; Touching is not penetrating. Bodies laid out edge to edge report
    ;; nothing until gravity presses them together, which is correct and
    ;; is worth knowing when a scene looks inert on its first frame.
    (is (empty? (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.0 1.0 0.0])))))

  (testing "a rotated box that clears the corner is separated"
    ;; Diagonally offset far enough that an edge-edge axis separates
    ;; them, which the face axes alone would not catch.
    (is (empty? (ct/between 0 1 (box [0.0 0.0 0.0])
                            (box [0.9 0.9 0.0] [1.0 1.0 1.0] (q/from-euler 0.0 0.0 0.7)))))))

(deftest sphere-test
  (testing "two spheres meet at one point on the line between them"
    (let [cs (ct/between 0 1 (ball [0.0 0.0 0.0] 1.0) (ball [1.5 0.0 0.0] 1.0))]
      (is (= 1 (count cs)))
      (is (< (abs (- 0.5 (:depth (first cs)))) 1e-9))
      (is (< (v/distance (:normal (first cs)) [1.0 0.0 0.0]) 1e-9))))

  (testing "a sphere on a box top is pushed straight up"
    (let [cs (ct/between 0 1 (ball [0.0 0.9 0.0] 0.5) (box [0.0 0.0 0.0]))]
      (is (= 1 (count cs)))
      ;; Normal runs from the sphere towards the box, so downwards.
      (is (< (v/distance (:normal (first cs)) [0.0 -1.0 0.0]) 1e-9))
      (is (< (abs (- 0.1 (:depth (first cs)))) 1e-9))))

  (testing "a sphere swallowed by a box comes out of the nearest face"
    ;; No nearest surface point exists in the usual sense once the centre
    ;; is inside, and a solver handed a zero normal does nothing at all.
    (let [cs (ct/between 0 1 (ball [0.1 0.0 0.0] 0.2) (box [0.0 0.0 0.0] [4.0 1.0 1.0]))]
      (is (= 1 (count cs)))
      (is (unit? (:normal (first cs))))
      (is (pos? (:depth (first cs))))))

  (testing "a sphere clear of a box does not touch it"
    (is (empty? (ct/between 0 1 (ball [0.0 2.0 0.0] 0.5) (box [0.0 0.0 0.0]))))))

(deftest orientation-test
  (let [a (box [0.0 0.0 0.0])
        b (box [0.0 0.95 0.0])]
    (testing "the normal runs from a to b, whichever way the pair is given"
      (let [fwd (ct/between 0 1 a b)
            rev (ct/between 1 0 b a)]
        (is (every? #(pos? (v/dot (:normal %) [0.0 1.0 0.0])) fwd))
        (is (every? #(neg? (v/dot (:normal %) [0.0 1.0 0.0])) rev))))

    (testing "a box against a sphere reverses the sphere-box normal"
      (let [s (ball [0.0 0.9 0.0] 0.5)
            sb (ct/between 0 1 s a)
            bs (ct/between 0 1 a s)]
        (is (seq sb))
        (is (seq bs))
        (is (< (v/distance (:normal (first sb))
                           (v/negate (:normal (first bs))))
               1e-9))))))

(deftest invariants-test
  (let [bodies (vec (cons (rigid/box {:pos [3.0 -0.25 0.0] :size [40.0 0.5 5.0]})
                          (for [row (range 6) col (range 5)]
                            (box [(+ (* col 0.99) (if (odd? row) 0.5 0.0))
                                  (+ 0.245 (* row 0.49)) 0.0]
                                 [1.0 0.5 0.5]))))
        cs (ct/all bodies)]

    (testing "a wall of bricks produces contacts everywhere it touches"
      (is (> (count cs) 100)))

    (testing "every normal is a unit vector and every depth is positive"
      (is (every? #(unit? (:normal %)) cs))
      (is (every? #(pos? (:depth %)) cs)))

    (testing "nothing is reported as deeply buried in a resting wall"
      (is (every? #(< (:depth %) 0.2) cs)))

    (testing "contacts name real bodies, a before b"
      (is (every? #(< -1 (:a %) (count bodies)) cs))
      (is (every? #(< -1 (:b %) (count bodies)) cs))
      (is (every? #(< (:a %) (:b %)) cs)))

    (testing "two static bodies are never tested against each other"
      ;; Neither can move, so whatever they would say cannot matter, and
      ;; a scene with a lot of scenery would spend its frame on it.
      (let [statics [(rigid/box {:pos [0.0 0.0 0.0] :size [2.0 2.0 2.0]})
                     (rigid/box {:pos [0.5 0.0 0.0] :size [2.0 2.0 2.0]})]]
        (is (empty? (ct/all statics)))))))
