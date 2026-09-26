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

  (testing "boxes exactly abutting touch, at zero depth"
    ;; They are not penetrating, and they still get a contact: the solver
    ;; is told about a touch it can hold rather than a penetration it has
    ;; to undo. A stack laid out brick on brick is in contact on its
    ;; first frame, before gravity has pressed anything together, which
    ;; is the difference between a column that settles and one that free
    ;; falls a frame and then has to be caught.
    (let [cs (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.0 1.0 0.0]))]
      (is (= 4 (count cs)))
      (is (every? #(zero? (:depth %)) cs))))

  (testing "and a gap inside the speculative margin is reported as a gap"
    ;; Negative depth is the signal: these are apart by that much, and
    ;; the solver may let them close it but no faster.
    (let [gap (* 0.5 ct/speculative)
          cs  (ct/between 0 1 (box [0.0 0.0 0.0]) (box [0.0 (+ 1.0 gap) 0.0]))]
      (is (seq cs))
      (is (every? #(< (abs (- (- gap) (:depth %))) 1e-9) cs)))
    (is (empty? (ct/between 0 1 (box [0.0 0.0 0.0])
                            (box [0.0 (+ 1.0 (* 2.0 ct/speculative)) 0.0])))))

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
      ;; Normal runs from the sphere toward the box, so downward.
      (is (< (v/distance (:normal (first cs)) [0.0 -1.0 0.0]) 1e-9))
      (is (< (abs (- 0.1 (:depth (first cs)))) 1e-9))))

  (testing "a sphere swallowed by a box comes out of the nearest face"
    ;; No nearest surface point exists in the usual sense once the center
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

    (testing "every normal is a unit vector, and no depth is a tunnel"
      (is (every? #(unit? (:normal %)) cs))
      ;; Depth is signed now: an overlap is positive, and a gap inside
      ;; the speculative margin is negative but never further away than
      ;; the margin itself.
      (is (every? #(> (:depth %) (- ct/speculative)) cs)))

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

(deftest margin-test
  (testing "a still pair is given the fixed margin and no more"
    (let [a (box [0.0 0.0 0.0])
          b (box [0.0 3.0 0.0])]
      (is (= ct/speculative (ct/margin a b (/ 1.0 60.0))))))

  (testing "and a moving one is given what it can cross before the next look"
    ;; This is the whole of continuous detection's cheap half. A contact
    ;; the pair is not near enough to be offered is a contact the solver
    ;; never sees, and a body crossing more ground in a step than the
    ;; margin is wide steps straight over the window.
    (let [still (box [0.0 3.0 0.0])
          fast  (rigid/ball {:pos [0.0 0.0 0.0] :radius 0.5 :density 1.0
                             :vel [0.0 60.0 0.0]})
          dt    (/ 1.0 60.0)]
      (is (< (- (ct/margin still fast dt) (+ ct/speculative 1.0)) 1e-9))
      ;; A body nothing can move contributes nothing, however its
      ;; velocity field happens to read.
      (is (= ct/speculative (ct/margin still (assoc fast :inv-mass 0.0) dt)))))

  (testing "a body is offered a contact with what it is about to reach"
    ;; A meter short of the floor and closing at sixty meters a second,
    ;; which is a meter of travel in the step. Standing still it is a
    ;; meter of daylight and nothing to report; moving, it is a contact
    ;; with a negative depth -- a gap the solver is allowed to see
    ;; coming and stop at.
    (let [floor (rigid/box {:pos [0.0 -0.5 0.0] :size [20.0 1.0 20.0]})
          ball  (rigid/ball {:pos [0.0 1.5 0.0] :radius 0.5 :density 1.0
                             :vel [0.0 -60.0 0.0]})]
      (is (empty? (ct/all [floor ball] nil 0.0)))
      (let [cs (ct/all [floor ball] nil (/ 1.0 60.0))]
        (is (= 1 (count cs)))
        (is (neg? (:depth (first cs))))))))

(deftest torus-box-test
  (let [floor (rigid/box {:pos [0.0 -0.5 0.0] :size [20.0 1.0 20.0]})
        big 0.3 r 0.06
        wheel (fn [pos rot] (rigid/torus {:pos pos :rot rot :major big :minor r :density 1.0}))]
    (testing "an upright wheel touches the floor once, straight under its hub"
      (let [cs (ct/between 0 1 (wheel [0.2 (- (+ big r) 0.01) 0.1] q/identity-q) floor)
            c (first cs)]
        (is (= 1 (count cs)))
        (is (< (abs (- 0.01 (double (:depth c)))) 1e-6))
        (is (< (v/distance (:point c) [0.2 0.0 0.1]) 1e-4))
        (is (< (v/distance (:normal c) [0.0 -1.0 0.0]) 1e-9) "from the wheel toward the floor")))
    (testing "a leaning wheel touches at the bottom of its own rim, not under its hub"
      ;; This is the whole reason for a torus. Leaned 40 degrees, a ball
      ;; would still touch directly beneath its center; a tire touches
      ;; where its rim is lowest, which is out in the plane of the wheel.
      (let [lean 0.7
            rot (q/from-axis-angle [1.0 0.0 0.0] lean)
            axis (q/rotate rot [0.0 0.0 1.0])
            hub [0.0 (- (+ (* big (Math/cos lean)) r) 0.005) 0.0]
            c (first (ct/between 0 1 (wheel hub rot) floor))
            ;; The point is on the floor's surface, so the tube's center is
            ;; the tube's radius less the overlap above it.
            tube-center (v/add (:point c) [0.0 (- r (double (:depth c))) 0.0])]
        (is (some? c))
        (is (< (abs (- 0.005 (double (:depth c)))) 1e-4))
        (is (< (abs (v/dot (v/sub tube-center hub) axis)) 1e-4) "in the plane of the wheel")
        (is (< (abs (- big (v/distance tube-center hub))) 1e-4) "on the rim")
        (is (> (abs (double (nth (:point c) 2))) 0.15) "well out from under the hub")))
    (testing "a wheel clear of the floor by more than the margin does not touch it"
      (is (empty? (ct/between 0 1 (wheel [0.0 (+ big r 0.1) 0.0] q/identity-q) floor))))
    (testing "the pair works in either order"
      (let [c (first (ct/between 0 1 floor (wheel [0.0 (- (+ big r) 0.01) 0.0] q/identity-q)))]
        (is (< (v/distance (:normal c) [0.0 1.0 0.0]) 1e-9))))))
