(ns allgo.joint-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.joint :as j]
            [allgo.physics.rigid :as r]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 60.0))

(defn- base [] (r/box {:size [0.3 0.3 0.3] :pos [0.0 0.0 0.0]}))
(defn- arm [pos] (r/box {:size [1.0 0.2 0.2] :density 800.0 :pos pos}))

(defn- world [bodies cs & [opts]]
  (merge {:bodies (vec bodies) :constraints (vec cs)
          :gravity [0.0 -10.0 0.0] :dt dt}
         opts))

(defn- advance [w n] (nth (iterate j/step w) n))

(defn- anchor-gap [bodies jt]
  (let [[p0 _ p1 _] (j/frames bodies jt)] (v/distance p0 p1)))

(defn- hinge-angle
  "How far the second frame has turned about the joint axis."
  [bodies jt]
  (let [[_ r0 _ r1] (j/frames bodies jt)]
    (j/signed-angle (q/rotate r0 [1 0 0]) (q/rotate r0 [0 1 0]) (q/rotate r1 [0 1 0]))))

;; The joint axis laid along world z, so gravity actually swings an arm
;; that reaches out along x.
(def ^:private z-axis (q/from-axis-angle [0 1 0] (- (/ Math/PI 2))))

(deftest signed-angle-test
  (testing "the angle from one vector to another about an axis"
    (is (< (abs (j/signed-angle [0 0 1] [1 0 0] [1 0 0])) 1e-12))
    (is (< (abs (- (/ Math/PI 2) (j/signed-angle [0 0 1] [1 0 0] [0 1 0]))) 1e-9))
    (is (< (abs (+ (/ Math/PI 2) (j/signed-angle [0 0 1] [0 1 0] [1 0 0]))) 1e-9)))

  (testing "it works past a quarter turn, where the sine alone would not"
    ;; asin cannot tell 3pi/4 from pi/4; the cosine settles it.
    (doseq [phi [0.3 1.2 2.0 3.0 -0.3 -1.2 -2.0 -3.0]]
      (let [b (q/rotate (q/from-axis-angle [0 0 1] phi) [1 0 0])]
        (is (< (abs (- phi (j/signed-angle [0 0 1] [1 0 0] b))) 1e-9)
            (str "at " phi)))))

  (testing "and wraps rather than running off"
    (is (<= (- Math/PI) (j/signed-angle [0 0 1] [1 0 0]
                                        (q/rotate (q/from-axis-angle [0 0 1] 3.1) [1 0 0]))
            Math/PI))))

(deftest frames-test
  (testing "a joint's frames start where it was put"
    (let [bodies [(base) (arm [0.5 0.0 0.0])]
          jt (j/ball bodies {:a 0 :b 1 :at [0.0 0.0 0.0]})]
      (is (< (anchor-gap bodies jt) 1e-12))))

  (testing "and are carried by their bodies"
    (let [bodies [(base) (arm [0.5 0.0 0.0])]
          jt (j/ball bodies {:a 0 :b 1 :at [0.0 0.0 0.0]})
          moved (assoc bodies 1 (r/set-pose (bodies 1) [2.0 1.0 0.0] q/identity-q))
          [_ _ p1 _] (j/frames moved jt)]
      (is (< (v/distance p1 [1.5 1.0 0.0]) 1e-12)
          "the frame moved with the body it is fixed to"))))

(deftest ball-test
  (let [bodies [(base) (assoc (r/box {:size [0.3 0.3 0.3] :density 800.0
                                      :pos [0.0 -1.0 0.0]})
                              :vel [3.0 0.0 0.0])]
        jt (j/ball bodies {:a 0 :b 1 :at [0.0 0.0 0.0]})
        w (world bodies [jt])]

    (testing "the anchors stay together while it swings"
      (doseq [n [1 60 200 400]]
        (is (< (anchor-gap (:bodies (advance w n)) jt) 0.01)
            (str "at step " n))))

    (testing "and it really is swinging"
      (let [tilts (map (fn [n]
                         (let [p (:pos (nth (:bodies (advance w n)) 1))]
                           (Math/acos (max -1.0 (min 1.0 (/ (- (second p)) (v/length p)))))))
                       (range 0 400 20))]
        (is (> (apply max tilts) 0.3))))

    (testing "the joint holds tighter the more substeps it is given"
      ;; What a hard constraint should do, and the check that caught a
      ;; scaled angular correction that made it diverge instead.
      (let [worst (fn [substeps]
                    (reduce max 0.0 (map #(anchor-gap (:bodies (advance (assoc w :substeps substeps) %)) jt)
                                         (range 0 200 4))))]
        (is (< (worst 40) (worst 20) (worst 10)))))

    (testing "it is the same constraint as a rigid distance of zero"
      (let [rc (r/distance-constraint bodies {:a 0 :b 1 :at [0.0 0.0 0.0]
                                              :other-at [0.0 0.0 0.0] :distance 0.0})
            wr (world bodies [rc])]
        (is (< (v/distance (:pos (nth (:bodies (advance w 120)) 1))
                           (:pos (nth (:bodies (advance wr 120)) 1)))
               1e-9)))))

  (testing "a swing limit keeps it inside a cone"
    (let [bodies [(base) (assoc (r/box {:size [0.3 0.3 0.3] :density 800.0
                                        :pos [0.0 -1.0 0.0]})
                                :vel [4.0 0.0 0.0])]
          limited (j/ball bodies {:a 0 :b 1 :at [0.0 0.0 0.0]
                                  :swing-min -0.5 :swing-max 0.5})
          free (j/ball bodies {:a 0 :b 1 :at [0.0 0.0 0.0]})
          tilt-of (fn [jt n]
                    (let [p (:pos (nth (:bodies (advance (world bodies [jt]) n)) 1))]
                      (Math/acos (max -1.0 (min 1.0 (/ (- (second p)) (v/length p)))))))
          worst (fn [jt] (reduce max 0.0 (map #(tilt-of jt %) (range 0 400 8))))]
      (is (< (worst limited) 0.65) (str "cone held to " (worst limited)))
      (is (> (worst free) (worst limited)) "and it is the limit doing it"))))

(deftest hinge-test
  (let [bodies [(base) (arm [0.5 0.0 0.0])]
        jt (j/hinge bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis})
        w (world bodies [jt])]

    (testing "the arm swings about the hinge axis"
      (is (> (abs (hinge-angle (:bodies (advance w 60)) jt)) 0.5)))

    (testing "and does not leave the plane the axis allows"
      ;; A hinge is the constraint that everything but one rotation is
      ;; taken away.
      (doseq [n [30 120 400]]
        (let [p (:pos (nth (:bodies (advance w n)) 1))]
          (is (< (abs (nth p 2)) 0.02) (str "drifted off-plane at step " n)))))

    (testing "the anchors stay together"
      (is (< (anchor-gap (:bodies (advance w 200)) jt) 0.02))))

  (testing "limits stop it where they are set"
    (let [bodies [(base) (arm [0.5 0.0 0.0])]
          jt (j/hinge bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis
                              :swing-min -0.4 :swing-max 0.4})
          w (world bodies [jt])]
      (doseq [n [60 200 400]]
        (let [phi (hinge-angle (:bodies (advance w n)) jt)]
          (is (<= -0.45 phi 0.45) (str "at step " n " the angle was " phi))))
      (testing "and it comes to rest against the limit it fell onto"
        (is (< (abs (+ 0.4 (hinge-angle (:bodies (advance w 400)) jt))) 0.02))))))

(deftest servo-test
  (testing "a servo drives to the angle it is told"
    (doseq [target [0.8 -0.6 1.5]]
      (let [bodies [(base) (arm [0.5 0.0 0.0])]
            jt (j/servo bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis
                                :target-angle target})
            w (world bodies [jt])]
        (is (< (abs (- target (hinge-angle (:bodies (advance w 400)) jt))) 0.05)
            (str "target " target))))))

(deftest motor-test
  (testing "a motor turns at the rate it is given"
    (let [bodies [(base) (arm [0.5 0.0 0.0])]
          jt (j/motor bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis :velocity 2.0})
          w (world bodies [jt] {:gravity v/zero})]
      (doseq [[n expected] [[30 1.0] [60 2.0] [90 3.0]]]
        (is (< (abs (- expected (hinge-angle (:bodies (advance w n)) jt))) 0.02)
            (str "after " n " steps")))))

  (testing "and at that rate whatever the solver is set to"
    ;; The target advances once a frame, not once a substep, or how fast a
    ;; motor turned would depend on how the solver was configured.
    (let [turn (fn [substeps]
                 (let [bodies [(base) (arm [0.5 0.0 0.0])]
                       jt (j/motor bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis
                                           :velocity 2.0})]
                   (hinge-angle (:bodies (advance (world bodies [jt]
                                                         {:gravity v/zero
                                                          :substeps substeps})
                                                  60))
                                jt)))]
      (is (< (abs (- (turn 5) (turn 20))) 0.02))))

  (testing "reversing the velocity reverses the turn"
    (let [spin (fn [velocity]
                 (let [bodies [(base) (arm [0.5 0.0 0.0])]
                       jt (j/motor bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis
                                           :velocity velocity})]
                   (hinge-angle (:bodies (advance (world bodies [jt] {:gravity v/zero}) 30)) jt)))]
      (is (pos? (spin 2.0)))
      (is (neg? (spin -2.0))))))

(deftest fixed-test
  (testing "a fixed joint welds two bodies into one"
    (let [a (r/box {:size [0.4 0.4 0.4] :density 800.0 :pos [0.0 2.0 0.0]})
          b (r/box {:size [0.4 0.4 0.4] :density 800.0 :pos [0.5 2.0 0.0]})
          bodies [(assoc a :omega [0.0 0.0 2.0]) b]
          jt (j/fixed bodies {:a 0 :b 1 :at [0.25 2.0 0.0]})
          w (world bodies [jt] {:gravity v/zero})
          after (:bodies (advance w 200))]
      (is (< (abs (- 0.5 (v/distance (:pos (after 0)) (:pos (after 1))))) 0.02)
          "they stay the distance apart they started")
      (is (< (q/angle (q/between (:rot (after 0)) (:rot (after 1)))) 0.05)
          "and keep the same orientation as one another"))))

(deftest prismatic-test
  (testing "a prismatic joint slides along its axis and nothing else"
    (let [bodies [(base) (r/box {:size [0.3 0.3 0.3] :density 800.0 :pos [0.0 1.0 0.0]})]
          ;; Joint axis pointing up, so gravity drives the slide.
          jt (j/prismatic bodies {:a 0 :b 1 :at [0.0 1.0 0.0]
                                  :rot (q/from-axis-angle [0 0 1] (/ Math/PI 2))
                                  :distance-min -0.5 :distance-max 0.5})
          after (:bodies (advance (world bodies [jt]) 400))
          [x y z] (:pos (after 1))]
      (is (< (abs x) 0.01) "no sideways movement")
      (is (< (abs z) 0.01))
      (is (< (abs (- 0.5 y)) 0.02) "and it came to rest against the lower stop")))

  (testing "it does not turn"
    (let [bodies [(base) (assoc (r/box {:size [0.3 0.3 0.3] :density 800.0
                                        :pos [0.0 1.0 0.0]})
                                :omega [0.0 3.0 0.0])]
          jt (j/prismatic bodies {:a 0 :b 1 :at [0.0 1.0 0.0]
                                  :rot (q/from-axis-angle [0 0 1] (/ Math/PI 2))})
          after (:bodies (advance (world bodies [jt] {:gravity v/zero}) 200))]
      (is (< (q/angle (:rot (after 1))) 0.15)))))

(deftest mixing-test
  (testing "joints and plain distance constraints live in one world"
    ;; Joints are `allgo.physics.rigid` constraints and nothing more, so
    ;; the step does not need to know which is which.
    (let [bodies [(base)
                  (r/box {:size [0.4 0.4 0.4] :density 800.0 :pos [1.0 0.0 0.0]})
                  (r/box {:size [0.4 0.4 0.4] :density 800.0 :pos [2.0 0.0 0.0]})]
          hinge (j/hinge bodies {:a 0 :b 1 :at [0.0 0.0 0.0] :rot z-axis})
          rod (r/distance-constraint bodies {:a 1 :b 2 :at [1.0 0.0 0.0]
                                             :other-at [2.0 0.0 0.0] :distance 1.0})
          after (:bodies (advance (world bodies [hinge rod]) 200))]
      (is (every? #(v/finite? (:pos %)) after))
      (is (< (abs (- 1.0 (v/distance (:pos (after 1)) (:pos (after 2))))) 0.02)
          "the rod holds its length")
      (is (< (anchor-gap after hinge) 0.05) "and the hinge holds its anchor"))))

(deftest relaxation-test
  (testing "applying an angular correction in full beats scaling it down"
    ;; The reference halves every position-level angular correction. It is
    ;; not a relaxation: the correction was sized against the full angular
    ;; response, so scaling only that half delivers something that no
    ;; longer opposes the error.
    (let [worst (fn [relax]
                  (let [bodies [(base) (assoc (r/box {:size [0.3 0.3 0.3] :density 800.0
                                                      :pos [0.0 -1.0 0.0]})
                                              :vel [3.0 0.0 0.0])]
                        jt (j/ball bodies (cond-> {:a 0 :b 1 :at [0.0 0.0 0.0]}
                                            relax (assoc :angular-relaxation relax)))]
                    (reduce max 0.0
                            (map #(anchor-gap (:bodies (advance (world bodies [jt]) %)) jt)
                                 (range 0 200 4)))))]
      (is (< (worst nil) 0.05) "in full, the anchors stay together")
      (is (> (worst 0.5) (* 5.0 (worst nil))) "halved, the joint comes apart")
      (is (< (worst nil) (worst 0.25))
          "and a quarter is better than a half, which no real damping is"))))
