(ns allgo.solver-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.solver :as s]
            [clojure.test :refer [deftest is testing]]))

(def ^:private solvers [:sequential-impulse :tgs :xpbd])
(def ^:private dt (/ 1.0 60.0))

(defn- floor [] (rigid/box {:pos [0.0 -0.5 0.0] :size [60.0 1.0 60.0]}))

(defn- brick
  ([pos] (brick pos [1.0 0.5 0.5]))
  ([pos size] (rigid/box {:pos pos :size size :density 1.0})))

(defn- run [solver bodies steps & [opts]]
  (loop [w (s/world bodies (merge {:solver solver} opts)) i 0]
    (if (= i steps) w (recur (s/step w dt) (inc i)))))

(defn- dynamics [w] (rest (:bodies w)))
(defn- y-of [b] (nth (:pos b) 1))
(defn- speed [b] (v/length (:vel b)))
(defn- finite? [b]
  (every? #(and (not (Double/isNaN %)) (not (Double/isInfinite %)))
          (concat (:pos b) (:vel b) (:omega b) (:rot b))))

(deftest resting-test
  (testing "a box dropped on the floor comes to rest on it"
    ;; Not above it and not through it, and -- the part that was wrong --
    ;; not launched off it. XPBD resolves the whole overlap in a substep
    ;; and reads velocity back from the position change, so without a
    ;; pass to take that energy out again a dropped box left the floor at
    ;; fifteen metres a second.
    (doseq [solver solvers]
      (let [w (run solver [(floor) (brick [0.0 3.0 0.0])] 240)
            b (first (dynamics w))]
        (is (< (abs (- 0.25 (y-of b))) 0.02) (str (name solver) " y=" (y-of b)))
        (is (< (speed b) 0.1) (str (name solver) " speed=" (speed b)))
        (is (finite? b) (name solver)))))

  (testing "and stays there rather than sinking through"
    (doseq [solver solvers]
      (let [w (run solver [(floor) (brick [0.0 0.25 0.0])] 300)]
        (is (> (y-of (first (dynamics w))) 0.2) (name solver))))))

(deftest stack-test
  (testing "a stack of bricks stands up"
    ;; The thing contact manifolds exist for. One contact point per pair
    ;; would let every brick rock, and the tower would walk itself apart.
    (doseq [solver solvers]
      (let [bodies (vec (cons (floor)
                              (for [i (range 6)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            w (run solver bodies 180)
            ys (sort (map y-of (dynamics w)))]
        (is (every? finite? (dynamics w)) (name solver))
        (is (every? (fn [[a b]] (< 0.42 (- b a) 0.56)) (partition 2 1 ys))
            (str (name solver) " courses out of order or crushed: " (vec ys)))
        (is (every? #(< (abs (nth (:pos %) 0)) 0.2) (dynamics w))
            (str (name solver) " drifted sideways"))
        (is (< (apply max (map speed (dynamics w))) 0.3)
            (str (name solver) " never settled")))))

  (testing "nothing climbs while resting"
    ;; A settled stack should not be rising. Baumgarte stabilisation does
    ;; put energy in, which is why the bias is a fraction of the overlap
    ;; and not all of it.
    (doseq [solver solvers]
      (let [bodies (vec (cons (floor)
                              (for [i (range 4)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            settled (run solver bodies 120)
            later (loop [w settled i 0] (if (= i 120) w (recur (s/step w dt) (inc i))))]
        (is (< (apply max (map y-of (dynamics later)))
               (+ 0.05 (apply max (map y-of (dynamics settled)))))
            (name solver))))))

(deftest left-alone-test
  (testing "a settled stack is still settled a long time later"
    ;; The failure this exists for does not show up in the first second.
    ;; XPBD reads velocity back off the position correction, at the
    ;; substep rate, so anything the position solve leaves behind is
    ;; multiplied by 240 -- and with nothing to take it out again it
    ;; compounds. A six brick column sat at six millimetres a second for
    ;; four seconds and was doing twenty nine metres a second at five.
    ;; Any test that stops at three seconds calls that stable.
    (doseq [solver solvers]
      (let [bodies (vec (cons (floor)
                              (for [i (range 5)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            w (run solver bodies 600)]
        (is (every? finite? (dynamics w)) (name solver))
        (is (< (apply max (map speed (dynamics w))) 0.3)
            (str (name solver) " wound itself up: "
                 (vec (map speed (dynamics w)))))
        (is (every? #(< (abs (nth (:pos %) 0)) 0.3) (dynamics w))
            (str (name solver) " walked sideways"))))))

(deftest stays-put-test
  (testing "a settled column does not merely settle, it stops"
    ;; The difference is the whole point. A stack under a Baumgarte bias
    ;; comes to rest in the sense that it stops going anywhere, and never
    ;; stops *moving*: it jitters at a centimetre a second forever,
    ;; because the velocity invented to push the overlap out is left in
    ;; the bodies and gravity puts it back. Over a minute that walks a
    ;; wall sideways and over it goes. A soft contact has no such
    ;; velocity to leave behind, and the column below is still to every
    ;; digit printed.
    ;; Sequential impulse only, and the other two say why. TGS divides
    ;; the iteration budget by the substeps, so each substep gets two
    ;; sweeps to cancel the gravity it just added and a five course stack
    ;; keeps about one substep of it -- 0.05 m/s that no contact model
    ;; fixes. XPBD is still quietly pushing overlap out through the
    ;; position solve at this point and creeps a centimetre over the five
    ;; seconds. Both are in TODOs.
    (doseq [solver [:sequential-impulse]]
      (let [bodies (vec (cons (floor)
                              (for [i (range 6)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            settled (run solver bodies 300)
            where (mapv :pos (dynamics settled))
            later (loop [w settled i 0] (if (= i 300) w (recur (s/step w dt) (inc i))))]
        (is (< (apply max (map speed (dynamics later))) 0.005)
            (str (name solver) " still moving: "
                 (vec (map speed (dynamics later)))))
        (is (every? #(< % 0.002)
                    (map v/distance where (map :pos (dynamics later))))
            (str (name solver) " wandered: "
                 (vec (map v/distance where (map :pos (dynamics later))))))))))

(deftest sleeping-test
  (testing "a settled stack goes to sleep, and is then exactly still"
    ;; Exactly is the word. Everything else in this file measures how
    ;; little a stack moves; a sleeping one does not move, because it is
    ;; not being integrated at all. That is the only thing that takes the
    ;; solver's own residual to zero, and the residual is what an
    ;; unstable arrangement amplifies until it falls over.
    (doseq [solver solvers]
      (let [bodies (vec (cons (floor)
                              (for [i (range 5)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            w (run solver bodies 300)
            bs (dynamics w)]
        (is (every? rigid/sleeping? bs)
            (str (name solver) " still awake after five seconds"))
        (is (every? #(zero? (double (speed %))) bs) (name solver))
        (let [where (mapv :pos bs)
              later (loop [w w i 0] (if (= i 300) w (recur (s/step w dt) (inc i))))]
          (is (= where (mapv :pos (dynamics later)))
              (str (name solver) " moved in its sleep"))))))

  (testing "and wakes when something runs into it"
    (doseq [solver solvers]
      (let [bodies (vec (cons (floor)
                              (for [i (range 5)]
                                (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            settled (run solver bodies 300)
            ;; A ball, fired along the ground at the foot of the stack.
            fired (update settled :bodies conj
                          (rigid/ball {:pos [4.0 0.3 0.0] :radius 0.3
                                       :density 4.0 :vel [-18.0 0.0 0.0]}))
            hit (loop [w fired i 0] (if (= i 40) w (recur (s/step w dt) (inc i))))]
        (is (not-any? rigid/sleeping? (dynamics hit))
            (str (name solver) " slept through being hit"))
        (is (every? finite? (dynamics hit)) (name solver)))))

  (testing "and does not sleep when it is told not to"
    (let [bodies (vec (cons (floor) [(brick [0.0 0.25 0.0])]))
          w (run :sequential-impulse bodies 300 {:allow-sleep? false})]
      (is (not-any? rigid/sleeping? (dynamics w))))))

(deftest pushout-bound-test
  (testing "however deep the overlap, nothing is fired out of it"
    ;; A box dropped half inside another is not a scene anyone builds on
    ;; purpose; it is what a stack too tall for the solver turns into,
    ;; one course at a time, and what the solver does about it decides
    ;; whether the wall falls over or leaves the building. XPBD reads
    ;; velocity back off the position correction, over the substep, so
    ;; an unbounded correction is an unbounded velocity -- a sixteen
    ;; course wall reached sixty metres a second that way. The bound is
    ;; `max-push-speed`, which the impulse solvers have always had.
    (doseq [solver solvers]
      (let [limit 3.0
            bodies [(floor)
                    (brick [0.0 0.25 0.0])
                    ;; Half a brick down into the one below it.
                    (brick [0.0 0.5 0.0])]
            w (run solver bodies 1 {:max-push-speed limit})
            fastest (apply max 0.0 (map speed (dynamics w)))]
        ;; The bound, and one step of gravity on top of it.
        (is (< fastest (+ limit 0.2))
            (str (name solver) " left the overlap at " fastest " m/s"))
        (is (every? finite? (dynamics w)) (name solver))))))

(deftest restitution-test
  (testing "a bouncy ball bounces and a dead one does not"
    (doseq [solver solvers]
      (let [drop-it (fn [e]
                      (let [w (run solver [(floor) (rigid/ball {:pos [0.0 2.0 0.0] :radius 0.3
                                                                :density 1.0})]
                                   60 {:restitution e})]
                        (apply max (map y-of (dynamics w)))))]
        ;; Measured while still in flight: the elastic one is higher up
        ;; at the same moment because it has already come back.
        (is (> (drop-it 0.8) (drop-it 0.0)) (name solver))))))

(deftest friction-test
  (testing "friction decides whether a box slides down a ramp"
    ;; A slope steep enough to slide on when it is slippery and not when
    ;; it grips. The comparison is the test -- the absolute distance
    ;; depends on the solver.
    (doseq [solver solvers]
      (let [ramp (rigid/box {:pos [0.0 -0.5 0.0] :size [40.0 1.0 40.0]
                             :rot (q/from-euler 0.0 0.0 -0.35)})
            block (rigid/box {:pos [0.0 0.75 0.0] :size [0.6 0.4 0.6] :density 1.0})
            slide (fn [mu]
                    (let [w (run solver [ramp block] 120 {:friction mu})]
                      (nth (:pos (first (dynamics w))) 0)))]
        (is (> (slide 0.02) (slide 1.2))
            (str (name solver) " slippery=" (slide 0.02) " grippy=" (slide 1.2)))))))

(deftest determinism-test
  (testing "the same scene steps the same way every time"
    (doseq [solver solvers]
      (let [scene #(vec (cons (floor) (for [i (range 4)] (brick [0.01 (+ 0.3 (* i 0.55)) 0.0]))))
            a (run solver (scene) 60)
            b (run solver (scene) 60)]
        (is (= (map :pos (:bodies a)) (map :pos (:bodies b))) (name solver))))))

(deftest warm-start-test
  (testing "warm starting holds a stack that cold starting lets sag"
    ;; The impulse a contact needed last step is a good guess at what it
    ;; needs now, and starting from it is most of why eight iterations is
    ;; enough. Only the impulse solvers have it to turn off.
    (doseq [solver [:sequential-impulse :tgs]]
      (let [bodies (vec (cons (floor)
                              (for [i (range 8)] (brick [0.0 (+ 0.25 (* i 0.5)) 0.0]))))
            top (fn [warm?]
                  (apply max (map y-of (dynamics (run solver bodies 120
                                                      {:iterations 4 :warm-start? warm?})))))]
        (is (> (top true) (top false))
            (str (name solver) " warm=" (top true) " cold=" (top false)))))))
