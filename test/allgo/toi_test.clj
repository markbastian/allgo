(ns allgo.toi-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.toi :as toi]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 60.0))

(def ^:private slab
  (rigid/box {:pos [0.0 0.0 0.0] :size [20.0 20.0 0.4]}))

(defn- ball [z speed]
  (rigid/ball {:pos [0.0 0.0 z] :radius 0.2 :density 7.8
               :vel [0.0 0.0 (- (double speed))]}))

(defn- close? [^double a ^double b ^double tol] (< (abs (- a b)) tol))

(deftest time-of-impact-test
  (testing "it finds the moment a ball reaches a slab"
    ;; The slab's front face is at z = 0.2 and the ball's surface leads
    ;; its center by 0.2, so the two meet after (z0 - 0.4) / speed. There
    ;; is an exact answer here and conservative advancement should walk
    ;; down onto it.
    (doseq [[z0 speed] [[5.0 1000.0] [1.0 5000.0] [0.9 60.0] [2.0 300.0]]]
      (let [exact (/ (- (double z0) 0.4) (double speed))
            t (toi/time-of-impact (ball z0 speed) slab dt)]
        (is (some? t) (str "no impact found for z0=" z0 " v=" speed))
        (is (close? (double t) exact 1e-4)
            (str "z0=" z0 " v=" speed ": got " t " want " exact)))))

  (testing "and says nothing when nothing is going to happen"
    ;; Too slow to arrive inside the step, and moving the other way.
    (is (nil? (toi/time-of-impact (ball 5.0 10.0) slab dt)))
    (is (nil? (toi/time-of-impact (ball 5.0 -300.0) slab dt))))

  (testing "bodies already touching answer zero"
    ;; Which is the ordinary case, and is why a caller has to tell the
    ;; difference between `touching now` and `about to touch`.
    (is (zero? (double (toi/time-of-impact (ball 0.4 300.0) slab dt)))))

  (testing "it answers early rather than late when a body is spinning"
    ;; Rotation is bounded by `|omega| r` rather than simulated, so a
    ;; spinning body's answer is conservative -- sooner than the truth,
    ;; never later. A spinning ball is the case where the bound is pure
    ;; overstatement, since turning a ball moves none of its surface
    ;; anywhere new.
    (let [still (toi/time-of-impact (ball 2.0 300.0) slab dt)
          spun (toi/time-of-impact (assoc (ball 2.0 300.0) :omega [0.0 40.0 0.0])
                                   slab dt)]
      (is (< (double spun) (double still)))))

  (testing "closing speed comes back with the time"
    ;; What a caller needs to ask for a little past the impact rather
    ;; than exactly onto it.
    (let [{:keys [t closing]} (toi/impact (ball 2.0 300.0) slab dt)]
      (is (close? (double t) (/ 1.6 300.0) 1e-4))
      (is (close? (double closing) 300.0 1e-6)))))

(deftest support-test
  (testing "a box's support mapping turns with the box"
    (let [b (rigid/box {:pos [0.0 0.0 0.0] :size [2.0 0.4 0.4]
                        :rot (q/from-axis-angle [0.0 0.0 1.0] (/ Math/PI 2))})
          s ((toi/support b (:pos b)) [0.0 1.0 0.0])]
      ;; Stood on end, its furthest point along +y is a meter up rather
      ;; than a fifth of one.
      (is (close? (double (second s)) 1.0 1e-9) (str s))))

  (testing "and a ball's does not, because there is nothing to turn"
    (let [b (rigid/ball {:pos [1.0 2.0 3.0] :radius 0.5})
          s ((toi/support b (:pos b)) [1.0 0.0 0.0])]
      (is (close? (double (first s)) 1.5 1e-9))))

  (testing "reach is the furthest any part of a body is from its center"
    (is (close? (toi/reach (rigid/ball {:pos [0.0 0.0 0.0] :radius 0.7})) 0.7 1e-9))
    (is (close? (toi/reach (rigid/box {:pos [0.0 0.0 0.0] :size [2.0 2.0 2.0]}))
                (* 0.5 (Math/sqrt 12.0)) 1e-9))))
