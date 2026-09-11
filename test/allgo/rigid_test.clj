(ns allgo.rigid-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as r]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 60.0))

(defn- world [bodies constraints & [opts]]
  (merge {:bodies (vec bodies) :constraints (vec constraints)
          :gravity [0.0 -10.0 0.0] :dt dt}
         opts))

(defn- advance [w n] (nth (iterate r/step w) n))

(deftest construction-test
  (testing "a box gets the mass and inertia of a box"
    (let [b (r/box {:size [2 3 4] :density 5.0})
          mass (* 5.0 2 3 4)]
      (is (< (abs (- mass (:mass b))) 1e-9))
      (is (< (abs (- (/ 1.0 mass) (:inv-mass b))) 1e-12))
      (is (< (abs (- (* (/ 1.0 12.0) mass (+ 9 16)) (first (:inertia b)))) 1e-9))
      (testing "and it is heavier to spin about its long axis than its short one"
        (let [thin (r/box {:size [0.1 3 0.1] :density 1.0})
              [ix iy _] (:inertia thin)]
          (is (> ix iy) "harder about x, which the length is across")))))

  (testing "a ball gets 2/5 m r^2 about every axis"
    (let [b (r/ball {:radius 2.0 :density 3.0})
          mass (* 3.0 (/ 4.0 3.0) Math/PI 8.0)]
      (is (< (abs (- mass (:mass b))) 1e-9))
      (is (apply = (:inertia b)))
      (is (< (abs (- (* 0.4 mass 4.0) (first (:inertia b)))) 1e-9))))

  (testing "no density makes it static"
    (let [b (r/box {:size [10 1 10]})]
      (is (r/static? b))
      (is (zero? (:inv-mass b)))
      (is (zero? (r/inverse-mass b [0 1 0] [1 2 3])))))

  (testing "a body keeps its inverse orientation in step with its orientation"
    (let [b (r/box {:size [1 1 1] :density 1.0 :rot (q/from-euler 0.3 0.4 0.5)})]
      (is (< (v/distance [1.0 2.0 3.0]
                         (q/rotate (:inv-rot b) (q/rotate (:rot b) [1 2 3])))
             1e-12)))))

(deftest frames-test
  (testing "local and world are inverses of one another"
    (let [b (r/box {:size [1 1 1] :density 1.0 :pos [3 -1 2]
                    :rot (q/from-euler 0.2 -0.9 1.4)})
          p [1.5 0.25 -3.0]]
      (is (< (v/distance p (r/local->world b (r/world->local b p))) 1e-12))))

  (testing "a body's own centre is the origin of its frame"
    (let [b (r/box {:size [1 1 1] :density 1.0 :pos [3 -1 2]})]
      (is (< (v/length (r/world->local b [3 -1 2])) 1e-12)))))

(deftest integration-test
  (testing "free fall matches the closed form"
    (let [w (world [(r/box {:size [1 1 1] :density 1.0 :pos [0 10 0]})] [])
          [_ y _] (:pos (first (:bodies (advance w 60))))]
      ;; Half g t squared, one second in.
      (is (< (abs (- y 5.0)) 0.05) (str "fell to " y))))

  (testing "a static body is not moved by anything"
    (let [w (world [(r/box {:size [10 1 10] :pos [0 0 0]})] [])
          b (first (:bodies (advance w 120)))]
      (is (= [0.0 0.0 0.0] (mapv double (:pos b))))
      (is (= v/zero (mapv double (:vel b))))))

  (testing "a free body keeps spinning at the rate it was given"
    ;; Nothing acts on it, so angular momentum is conserved, and a cube's
    ;; inertia is the same about every axis so the rate is too.
    (let [b (assoc (r/box {:size [1 1 1] :density 1.0}) :omega [0.0 3.0 0.0])
          w (world [b] [] {:gravity v/zero})
          after (first (:bodies (advance w 180)))]
      (is (< (v/distance [0.0 3.0 0.0] (:omega after)) 0.05)
          (str "omega drifted to " (:omega after)))))

  (testing "orientation stays a rotation as it is integrated"
    ;; A first-order step walks off the unit sphere; renormalizing is what
    ;; keeps it a rotation rather than a rotation plus a scale.
    (let [b (assoc (r/box {:size [1 2 3] :density 1.0}) :omega [1.5 -2.0 0.7])
          w (world [b] [] {:gravity v/zero})
          after (first (:bodies (advance w 600)))]
      (is (< (abs (- 1.0 (q/length (:rot after)))) 1e-9))))

  (testing "damping takes energy out and nothing puts it back"
    (let [spin (fn [damping]
                 (let [b (assoc (r/box {:size [1 1 1] :density 1.0
                                        :angular-damping damping})
                                :omega [0.0 3.0 0.0])]
                   (v/length (:omega (first (:bodies (advance (world [b] [] {:gravity v/zero}) 60)))))))]
      (is (< (spin 1.0) (spin 0.0)))
      (is (< (spin 3.0) (spin 1.0))))))

(deftest inverse-mass-test
  (let [b (r/box {:size [1 1 1] :density 1.0})]
    (testing "a correction through the centre sees only the mass"
      (is (< (abs (- (:inv-mass b) (r/inverse-mass b [0 1 0] (:pos b)))) 1e-12)))

    (testing "one at arm's length sees more, because it also spins the body"
      (is (> (r/inverse-mass b [0 1 0] [1.0 0.0 0.0])
             (r/inverse-mass b [0 1 0] (:pos b)))))

    (testing "and the further out, the more"
      (is (< (r/inverse-mass b [0 1 0] [0.5 0 0])
             (r/inverse-mass b [0 1 0] [1.0 0 0])
             (r/inverse-mass b [0 1 0] [2.0 0 0]))))

    (testing "a correction along the arm does not spin it at all"
      ;; Pushing a body straight at its own centre cannot turn it.
      (is (< (abs (- (:inv-mass b) (r/inverse-mass b [1 0 0] [1.0 0.0 0.0]))) 1e-12)))

    (testing "with no point at all the question is purely angular"
      (is (< (abs (- (first (:inv-inertia b)) (r/inverse-mass b [1 0 0] nil))) 1e-12))
      (is (not= (r/inverse-mass b [1 0 0] nil) (r/inverse-mass b [1 0 0] (:pos b)))))))

(deftest distance-constraint-test
  (testing "a pendulum keeps its length, exactly, indefinitely"
    (let [anchor [0.0 0.0 0.0]
          b (r/box {:size [0.4 0.4 0.4] :density 1.0 :pos [1.0 0.0 0.0]})
          bodies [b]
          c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                           :other-at anchor :distance 1.0})
          w (world bodies [c])]
      (doseq [n [1 30 120 600]]
        (let [p (:pos (first (:bodies (advance w n))))]
          (is (< (abs (- 1.0 (v/distance p anchor))) 1e-6)
              (str "at step " n " the length was " (v/distance p anchor)))))))

  (testing "and it swings rather than hanging still"
    (let [anchor [0.0 0.0 0.0]
          b (r/box {:size [0.4 0.4 0.4] :density 1.0 :pos [1.0 0.0 0.0]})
          bodies [b]
          c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                           :other-at anchor :distance 1.0})
          w (world bodies [c])
          ys (map (fn [n] (second (:pos (first (:bodies (advance w n))))))
                  (range 0 90 5))]
      (is (< (apply min ys) -0.9) "it reaches the bottom of its swing")
      (is (> (apply max ys) -0.2) "and comes back up")))

  (testing "a pendulum does not gain energy"
    ;; Released from horizontal, it should never rise above where it
    ;; started. A solver that injects energy shows it here first.
    (let [anchor [0.0 0.0 0.0]
          b (r/box {:size [0.4 0.4 0.4] :density 1.0 :pos [1.0 0.0 0.0]})
          bodies [b]
          c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                           :other-at anchor :distance 1.0})
          w (world bodies [c])
          highest (apply max (map (fn [n] (second (:pos (first (:bodies (advance w n))))))
                                  (range 0 600 7)))]
      (is (< highest 0.02) (str "rose to " highest))))

  (testing "two bodies joined by a rod stay a rod apart"
    (let [a (r/box {:size [0.3 0.3 0.3] :density 1.0 :pos [-1.0 0.0 0.0]})
          b (r/box {:size [0.3 0.3 0.3] :density 1.0 :pos [1.0 0.0 0.0]})
          bodies [a b]
          c (r/distance-constraint bodies {:a 0 :b 1 :at (:pos a) :other-at (:pos b)
                                           :distance 2.0})
          w (world bodies [c] {:gravity [0.0 0.0 0.0]})
          after (:bodies (advance w 120))]
      (is (< (abs (- 2.0 (v/distance (:pos (after 0)) (:pos (after 1))))) 1e-6))))

  (testing "a rope pulls but does not push"
    (let [anchor [0.0 0.0 0.0]
          make (fn [unilateral?]
                 (let [b (r/box {:size [0.2 0.2 0.2] :density 1.0 :pos [0.0 -0.5 0.0]})
                       bodies [b]
                       c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                                        :other-at anchor :distance 1.0
                                                        :unilateral? unilateral?})]
                   (world bodies [c])))
          rope (:pos (first (:bodies (advance (make true) 1))))
          rod  (:pos (first (:bodies (advance (make false) 30))))]
      ;; Starting slack, a rope lets the body keep falling; a rod yanks it
      ;; out to length immediately.
      (is (< (v/distance rope anchor) 1.0) "the rope stays slack")
      (is (< (abs (- 1.0 (v/distance rod anchor))) 1e-6) "the rod does not")))

  (testing "compliance makes it springy"
    (let [hang (fn [compliance]
                 (let [b (r/box {:size [0.2 0.2 0.2] :density 1.0 :pos [0.0 -1.0 0.0]})
                       bodies [b]
                       c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                                        :other-at [0.0 0.0 0.0]
                                                        :distance 1.0
                                                        :compliance compliance})]
                   (v/distance (:pos (first (:bodies (advance (world bodies [c]) 120))))
                               [0.0 0.0 0.0])))]
      (is (< (abs (- 1.0 (hang 0.0))) 1e-6) "rigid holds its length")
      (is (> (hang 0.01) 1.0) "and a compliant one stretches under the weight")))

  (testing "an unknown constraint is an error rather than a silent nothing"
    (is (thrown? clojure.lang.ExceptionInfo
                 (r/solve [(r/box {:size [1 1 1] :density 1.0})] {:kind :nonsense} dt)))))

(deftest forces-test
  (testing "the step reports what each constraint carried"
    (let [b (r/box {:size [0.2 0.2 0.2] :density 1.0 :pos [0.0 -1.0 0.0]})
          bodies [b]
          c (r/distance-constraint bodies {:a 0 :b nil :at (:pos b)
                                           :other-at [0.0 0.0 0.0] :distance 1.0})
          w (advance (world bodies [c]) 120)]
      (is (= 1 (count (:forces w))))
      ;; Hanging still, the rod carries the body's weight.
      (is (< (abs (- (abs (double (first (:forces w)))) (* (:mass b) 10.0)))
             (* 0.25 (:mass b) 10.0))
          (str "carried " (first (:forces w)) " for a weight of " (* (:mass b) 10.0))))))
