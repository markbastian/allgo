(ns allgo.world-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.world :as w]
            [clojure.test :refer [deftest is testing]]))

(def ^:private idq [0.0 0.0 0.0 1.0])
(def ^:private g [0.0 -9.81 0.0])

(defn- box-inertia [[sx sy sz] m]
  (let [f (/ (double m) 12.0)
        sx (double sx) sy (double sy) sz (double sz)]
    [[(* f (+ (* sy sy) (* sz sz))) 0.0 0.0]
     [0.0 (* f (+ (* sx sx) (* sz sz))) 0.0]
     [0.0 0.0 (* f (+ (* sx sx) (* sy sy)))]]))

(def ^:private floor
  (rigid/box {:pos [0.0 -0.5 0.0] :size [40.0 1.0 40.0]}))

(def ^:private arm
  "A meter of limb on a ball joint, bolted at the origin and hanging."
  [{:parent -1 :joint :spherical :origin {:rot nil :pos [0.0 0.0 0.0]}
    :mass 3.0 :com [0.0 -0.5 0.0] :inertia (box-inertia [0.12 1.0 0.12] 3.0)
    :shape :box :size [0.12 1.0 0.12]}])

(defn- swinging [wz] {:q [idq] :qd [0.0 0.0 (double wz)]})

(defn- brick
  ([] (brick [0.30 -0.5 0.0]))
  ([pos] (rigid/box {:pos pos :size [0.25 0.25 0.25] :density 300.0})))

(defn- run [world n dt] (reduce (fn [s _] (w/step s dt)) world (range n)))

(defn- tip [pose]
  (let [f (first (ab/poses arm (:q pose)))]
    (v/add (:pos f) (q/rotate (:rot f) [0.0 -1.0 0.0]))))

(defn- close? [^double a ^double b ^double tol] (< (abs (- a b)) tol))

(deftest rigid-alone-test
  (testing "a loose body falls, and lands on a static floor"
    (let [wd (w/world [floor (brick [0.0 2.0 0.0])] [] {:gravity g})
          landed (run wd 480 (/ 1.0 240.0))
          b (second (:bodies landed))]
      ;; Half the brick's height, less the slop a soft contact keeps.
      (is (close? (double (second (:pos b))) 0.125 0.02)
          (str "came to rest at " (:pos b)))
      (is (< (v/length (:vel b)) 0.05)))))

(deftest model-alone-test
  (testing "a model with nothing to touch moves exactly as it would alone"
    ;; The world is supposed to add nothing when there is nothing to add.
    ;; If this drifts, the world is integrating the model differently
    ;; from the namespace that owns it, and every other result here is
    ;; suspect.
    (let [dt (/ 1.0 480.0)
          alone (reduce (fn [p _] (ab/step p arm dt {:gravity g}))
                        (swinging 4.0) (range 120))
          in-world (:pose (first (:models (run (w/world [] [{:model arm :pose (swinging 4.0)}]
                                                        {:gravity g})
                                               120 dt))))]
      (is (< (v/distance (tip alone) (tip in-world)) 1e-12))
      (is (every? #(< (abs (double %)) 1e-12) (map - (:qd alone) (:qd in-world)))))))

(deftest model-driven-test
  (testing "joint forces carried on a model move it exactly as they would alone"
    ;; A hinge driven by a torque, with nothing to touch: the world has to
    ;; hand `:tau` to the dynamics the way `ab/step` does, or a motor on a
    ;; model in a scene does nothing at all.
    (let [dt (/ 1.0 480.0)
          hinge [{:parent -1 :joint :revolute :axis [0.0 0.0 1.0]
                  :origin {:rot nil :pos [0.0 0.0 0.0]}
                  :mass 2.0 :com [0.5 0.0 0.0] :inertia (box-inertia [1.0 0.1 0.1] 2.0)}]
          tau [3.0]
          start {:q [0.0] :qd [0.0]}
          alone (reduce (fn [p _] (ab/step p hinge dt {:gravity [0.0 0.0 0.0] :tau tau}))
                        start (range 120))
          driven (:pose (first (:models (run (w/world [] [{:model hinge :pose start :tau tau}]
                                                      {:gravity [0.0 0.0 0.0]})
                                             120 dt))))]
      (is (pos? (double (first (:qd driven)))) "the torque turned it")
      (is (< (abs (- (double (first (:qd alone))) (double (first (:qd driven))))) 1e-12))
      (is (< (abs (- (double (first (:q alone))) (double (first (:q driven))))) 1e-12)))))

(deftest model-pushes-rigid-test
  (testing "a swung limb knocks a loose brick, and the brick pushes back"
    ;; The thing this namespace exists for. No gravity and no floor, so
    ;; the only thing that can move either of them is the other.
    (let [opts {:gravity [0.0 0.0 0.0] :iterations 12 :friction 0.4}
          dt (/ 1.0 480.0)
          free (run (w/world [] [{:model arm :pose (swinging 4.0)}] opts) 120 dt)
          hit (run (w/world [(brick)] [{:model arm :pose (swinging 4.0)}] opts) 120 dt)
          b (first (:bodies hit))]
      ;; Nothing in the way: the limb keeps the rate it started with.
      (is (close? (double (nth (:qd (:pose (first (:models free)))) 2)) 4.0 1e-9))
      ;; Something in the way: it is slowed, and by a lot.
      (is (< (double (nth (:qd (:pose (first (:models hit)))) 2)) 2.0)
          "the limb should have been slowed by what it hit")
      ;; And the brick has been sent off in the direction of the swing.
      (is (> (double (first (:pos b))) 0.4) (str "brick at " (:pos b)))
      (is (> (double (first (:vel b))) 0.5) (str "brick velocity " (:vel b))))))

(def ^:private free-arm
  "The same limb on a root that is free to move.

  Momentum only means anything for this one. A bolted model is nailed to
  the world and the nail takes whatever momentum it likes."
  {:base {:mass 5.0 :com [0.0 0.0 0.0]
          :inertia (box-inertia [0.3 0.3 0.3] 5.0)
          :shape :box :size [0.3 0.3 0.3]
          :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
   :links [{:parent -1 :joint :spherical :origin {:rot nil :pos [0.0 -0.15 0.0]}
            :mass 3.0 :com [0.0 -0.5 0.0] :inertia (box-inertia [0.12 1.0 0.12] 3.0)
            :shape :box :size [0.12 1.0 0.12]}]})

(defn- free-swinging [wz]
  {:q [idq] :qd [0.0 0.0 (double wz)]
   :base {:rot idq :pos [0.0 0.0 0.0] :vel [0.0 0.0 0.0 0.0 0.0 0.0]}})

(deftest exchange-conserves-momentum-test
  (testing "the exchange loses momentum only as fast as the integrator does"
    ;; A contact between a rigid body and an articulated one is an equal
    ;; and opposite pair of impulses, so it cannot move the total. What
    ;; can, and does, is the first-order integration underneath -- so the
    ;; test is that halving the step halves the loss, not that there is
    ;; none.
    (let [total (fn [s]
                  (let [m (first (:models s))
                        b (first (:bodies s))
                        bm (/ 1.0 (double (:inv-mass b)))]
                    (v/add (:linear (ab/momentum (:model m) (:pose m)))
                           (v/scale (:vel b) bm))))
          drift (fn [dt]
                  (let [wd (w/world [(brick [0.45 -0.8 0.0])]
                                    [{:model free-arm :pose (free-swinging 5.0)}]
                                    {:gravity [0.0 0.0 0.0] :iterations 20
                                     :restitution 0.0})
                        before (total wd)]
                    (v/distance before (total (run wd (long (/ 0.5 dt)) dt)))))
          coarse (drift (/ 1.0 240.0))
          fine (drift (/ 1.0 960.0))]
      (is (pos? coarse))
      ;; Four times the step, within a factor of two of four times the
      ;; drift.
      (is (< 2.0 (/ coarse fine) 8.0) (str "coarse " coarse " fine " fine)))))

(deftest contact-kinds-test
  (testing "all three kinds of contact are found"
    (let [model {:base {:mass 4.0 :com [0.0 0.0 0.0]
                        :inertia (box-inertia [0.3 0.3 0.3] 4.0)
                        :shape :box :size [0.3 0.3 0.3]
                        :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
                 :links [{:parent -1 :joint :spherical
                          :origin {:rot nil :pos [0.0 -0.1 0.0]}
                          :mass 2.0 :com [0.0 -0.3 0.0]
                          :inertia (box-inertia [0.12 0.6 0.12] 2.0)
                          :shape :box :size [0.12 0.6 0.12]}]}
          ;; Resting on the floor, with a brick against the limb, and
          ;; the limb folded back against its own root.
          pose {:q [(q/from-axis-angle [0.0 0.0 1.0] 2.9)] :qd [0.0 0.0 0.0]
                :base {:rot idq :pos [0.0 0.4 0.0]
                       :vel [0.0 0.0 0.0 0.0 0.0 0.0]}}
          wd (w/world [floor (brick [0.22 0.55 0.0])] [{:model model :pose pose}]
                      {:gravity g})
          cfgs [(ab/configuration model (:q pose) (:base pose))]
          cs (w/contacts wd cfgs)
          kind (fn [c] [(first (:a c)) (first (:b c))])
          kinds (set (map kind cs))]
      (is (seq cs))
      (is (contains? kinds [:link :rigid]) (str "kinds found: " kinds)))))

(deftest pin-to-the-world-test
  (testing "a loose box pinned to a fixed point swings from it and stays on it"
    (let [bx (rigid/box {:pos [0.6 2.0 0.0] :size [0.2 0.2 0.2] :density 500.0})
          wd (w/world [bx] [] {:gravity g
                               :pins [{:a [:rigid 0] :pa [-0.1 0.0 0.0]
                                       :b [:static] :pb [0.5 2.0 0.0]}]})
          states (reductions (fn [s _] (w/step s (/ 1.0 240.0))) wd (range 480))
          worst (apply max (map #(v/distance (rigid/local->world (first (:bodies %)) [-0.1 0.0 0.0])
                                             [0.5 2.0 0.0])
                                states))
          lowest (apply min (map #(double (second (:pos (first (:bodies %))))) states))]
      (is (< worst 0.01) (str "the pinned point strayed " worst))
      ;; A tenth of a meter from the pin, so it swings down to 1.9.
      (is (< lowest 1.91) "and it swung down"))))

(def ^:private hanging
  "A meter of limb on a hinge about z, bolted at `x` and hanging straight down."
  (fn [x]
    [{:parent -1 :joint :revolute :axis [0.0 0.0 1.0] :origin {:rot nil :pos [x 0.0 0.0]}
      :mass 2.0 :com [0.0 -0.5 0.0] :inertia (box-inertia [0.1 1.0 0.1] 2.0)}]))

(deftest pin-between-models-test
  (testing "two models pinned tip to tip close a loop"
    ;; Two hinged limbs half a meter apart, their tips pinned to a point
    ;; between them: a closed loop, which no single tree of reduced
    ;; coordinates can describe.
    (let [left (hanging -0.25)
          right (hanging 0.25)
          tip [0.0 -1.0 0.0]
          ;; Each tilted a quarter of a radian in, so the tips meet.
          wd (w/world [] [{:model left :pose {:q [0.2527] :qd [3.0]}}
                          {:model right :pose {:q [-0.2527] :qd [0.0]}}]
                      {:gravity g
                       :pins [{:a [:link 0 0] :pa tip :b [:link 1 0] :pb tip}]})
          ends (fn [world]
                 (let [configs (mapv (fn [{:keys [model pose]}] (ab/configuration model (:q pose) nil))
                                     (:models world))]
                   (apply v/distance (w/pin-points world configs (first (:pins world))))))
          start-gap (ends wd)
          later (run wd 480 (/ 1.0 240.0))]
      (is (< start-gap 0.01) "they start joined")
      (is (< (ends later) 0.01) (str "and stay joined, " (ends later) " apart"))
      ;; Two bars hinged apart and pinned at their tips are a triangle,
      ;; which cannot move: the swing given to the left one has to have
      ;; been taken out of it by the right one, through the pin.
      (is (< (abs (double (first (:qd (:pose (first (:models later))))))) 0.05)
          "the swing was stopped by the other model"))))

(def ^:private four-bar
  "A parallelogram as one model: two meter bars hinged half a meter
  apart, and a coupler hinged to the left one's tip. Its far end pinned
  to the right bar's tip closes the loop."
  [{:parent -1 :joint :revolute :axis [0.0 0.0 1.0] :origin {:rot nil :pos [-0.25 0.0 0.0]}
    :mass 2.0 :com [0.0 -0.5 0.0] :inertia (box-inertia [0.1 1.0 0.1] 2.0)}
   {:parent 0 :joint :revolute :axis [0.0 0.0 1.0] :origin {:rot nil :pos [0.0 -1.0 0.0]}
    :mass 1.0 :com [0.25 0.0 0.0] :inertia (box-inertia [0.5 0.1 0.1] 1.0)}
   {:parent -1 :joint :revolute :axis [0.0 0.0 1.0] :origin {:rot nil :pos [0.25 0.0 0.0]}
    :mass 2.0 :com [0.0 -0.5 0.0] :inertia (box-inertia [0.1 1.0 0.1] 2.0)}])

(deftest loop-in-one-model-test
  (testing "a pin between two links of one model closes a loop in it, and leaves it the freedom it has"
    ;; A tree cannot describe a four-bar linkage; the pin is the fourth
    ;; bar. What is left is one degree of freedom: the bars swing
    ;; together, the coupler stays level, and the swing keeps its height.
    (doseq [solver [:convex :sequential-impulse]]
      (let [wd (w/world [] [{:model four-bar :pose {:q [0.0 0.0 0.0] :qd [2.0 -2.0 2.0]}}]
                        {:gravity g :solver solver :self-collide? false
                         :pins [{:a [:link 0 1] :pa [0.5 0.0 0.0] :b [:link 0 2] :pb [0.0 -1.0 0.0]}]})
            gap (fn [s] (let [configs (mapv (fn [{:keys [model pose]}] (ab/configuration model (:q pose) nil))
                                            (:models s))]
                          (apply v/distance (w/pin-points s configs (first (:pins s))))))
            states (take 960 (iterate #(w/step % (/ 1.0 240.0)) wd))
            qs (map #(:q (:pose (first (:models %)))) states)
            worst (fn [f] (apply max (map f qs)))]
        (is (< (apply max (map gap states)) 1e-5) (str solver " stays joined"))
        (is (< (worst (fn [[a b _]] (abs (+ (double a) (double b))))) 1e-5) (str solver " the coupler level"))
        (is (< (worst (fn [[a _ c]] (abs (- (double a) (double c))))) 1e-5) (str solver " the bars together"))
        ;; From 2 rad/s at the bottom a meter bar's swing tops out near
        ;; 0.57 rad, on both sides, every time.
        (is (< 0.56 (worst (fn [[a]] (double a))) 0.58) (str solver " and swinging as high as it started"))
        (is (< 0.56 (worst (fn [[a]] (- (double a)))) 0.58))))))

(deftest damping-test
  (testing "a joint's passive damping brings a swinging limb to rest, and lets it sleep"
    ;; Nothing else touches the limb, so without damping it swings for
    ;; ever: a ragdoll with a shin hanging over an edge never slept.
    ;; Damping is the joint's own friction, not a drive, so it does not
    ;; keep the model awake the way a `:tau` does.
    (let [free (hanging 0.0)
          damped [(assoc (first free) :damping 2.0 :armature (* 2.0 (/ 1.0 60.0)))]
          after (fn [model] (first (:models (run (w/world [] [{:model model :pose {:q [0.0] :qd [3.0]}}] {:gravity g})
                                                 600 (/ 1.0 60.0)))))
          swinging (after free)
          stilled (after damped)]
      (is (> (abs (double (first (:qd (:pose swinging))))) 0.5) "undamped, it is still going")
      (is (not (:sleeping? swinging)))
      (is (:sleeping? stilled) "damped, it has stopped and gone to sleep"))))

(deftest pin-breaks-test
  (testing "a pin asked to hold more than its break force lets go"
    (let [bx (rigid/box {:pos [0.0 2.0 0.0] :size [0.3 0.3 0.3] :density 1000.0})
          pin {:a [:rigid 0] :pa [0.0 0.15 0.0] :b [:static] :pb [0.0 2.15 0.0]}
          weak (run (w/world [bx] [] {:gravity g :pins [(assoc pin :break-force 50.0)]})
                    120 (/ 1.0 240.0))
          strong (run (w/world [bx] [] {:gravity g :pins [(assoc pin :break-force 1000.0)]})
                      120 (/ 1.0 240.0))]
      ;; 27 kilograms hanging is 265 newtons.
      (is (empty? (:pins weak)))
      (is (= 1 (count (:broken-pins weak))))
      (is (< (double (second (:pos (first (:bodies weak))))) 1.9) "and it falls")
      (is (= 1 (count (:pins strong))))
      (is (< (abs (- 2.0 (double (second (:pos (first (:bodies strong))))))) 0.01)))))

(deftest limits-in-a-world-test
  (testing "a hinge's limit holds inside a world as it does alone"
    ;; `world` used to skip limits entirely, so a limited model in a scene
    ;; -- a ragdoll's knees -- went wherever gravity took it.
    (let [limb (assoc-in (hanging 0.0) [0 :limit] [-0.3 0.3])
          end (run (w/world [] [{:model limb :pose {:q [0.0] :qd [4.0]}}] {:gravity g})
                   480 (/ 1.0 240.0))
          angle (double (first (:q (:pose (first (:models end))))))]
      (is (< -0.32 angle 0.32) (str "at " angle)))))

(deftest wedge-test
  (testing "bricks started deep inside a tilted slab are pushed out, not thrown"
    ;; Each brick is wedged between the floor and a slab tilted over it,
    ;; whose normals are 162 degrees apart: asked to separate from both
    ;; at the push speed, the only answer held rigid is sideways at six
    ;; times it -- sixteen meters a second, and deeper, NaN.
    (let [slab (rigid/box {:pos [0.0 0.25 0.0] :size [1.6 0.5 1.6]
                           :rot (q/from-axis-angle [0.0 0.0 1.0] (Math/toRadians 18.0))})
          ;; four courses of four, the wall the ragdoll demo once stood
          ;; where its slab was
          bricks (for [r (range 4) c (range 4) :let [o (if (odd? r) 0.11 0.0)]]
                   (rigid/box {:pos [(+ (* (- c 1.5) 0.22) o 0.55) (+ 0.055 (* 0.11 r)) 0.0]
                               :size [0.22 0.11 0.14] :density 500.0}))
          states (take 120 (iterate #(w/step % (/ 1.0 60.0))
                                    (w/world (into [floor slab] bricks) [] {:gravity g :solver :convex})))
          fastest (apply max (for [s states b (drop 2 (:bodies s))] (v/length (:vel b))))]
      (is (< fastest 5.0) (str "a brick left at " fastest " m/s")))))

(def ^:private dummy
  "The ragdoll demo's figure: a free pelvis, ten box limbs on ball joints
  and hinges with the demo's cones, twists and ranges, and each joint
  damped a little, with the armature that damping wants."
  (let [deg #(Math/toRadians (double %))
        damp (fn [l] (assoc l :damping 0.2 :armature (/ 0.2 60.0)))
        bone (fn [dir parent joint axis origin size m]
               (let [[_ sy _] size]
                 (damp (cond-> {:parent parent :joint joint :origin {:rot nil :pos origin}
                                :mass m :com [0.0 (* 0.5 dir (double sy)) 0.0]
                                :inertia (box-inertia size m) :shape :box :size size}
                         axis (assoc :axis axis)))))
        ball (fn [l cone twist] (assoc l :cone (deg cone) :twist (deg twist)))
        hinge (fn [l lo hi] (assoc l :limit [(deg lo) (deg hi)]))]
    {:base {:mass 12.0 :com [0.0 0.0 0.0] :inertia (box-inertia [0.28 0.20 0.18] 12.0)
            :shape :box :size [0.28 0.20 0.18] :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
     :links
     [(ball (bone 1.0 -1 :spherical nil [0.0 0.10 0.0] [0.34 0.42 0.20] 25.0) 25 35)
      (ball (bone 1.0 0 :spherical nil [0.0 0.42 0.0] [0.20 0.22 0.20] 5.0) 40 60)
      (ball (bone -1.0 0 :spherical nil [0.23 0.38 0.0] [0.10 0.28 0.10] 2.5) 85 70)
      (hinge (bone -1.0 2 :revolute [1.0 0.0 0.0] [0.0 -0.28 0.0] [0.09 0.26 0.09] 1.8) 0 150)
      (ball (bone -1.0 0 :spherical nil [-0.23 0.38 0.0] [0.10 0.28 0.10] 2.5) 85 70)
      (hinge (bone -1.0 4 :revolute [1.0 0.0 0.0] [0.0 -0.28 0.0] [0.09 0.26 0.09] 1.8) 0 150)
      (ball (bone -1.0 -1 :spherical nil [0.09 -0.10 0.0] [0.14 0.40 0.14] 8.0) 60 35)
      (hinge (bone -1.0 6 :revolute [1.0 0.0 0.0] [0.0 -0.40 0.0] [0.12 0.38 0.12] 3.5) -150 0)
      (ball (bone -1.0 -1 :spherical nil [-0.09 -0.10 0.0] [0.14 0.40 0.14] 8.0) 60 35)
      (hinge (bone -1.0 8 :revolute [1.0 0.0 0.0] [0.0 -0.40 0.0] [0.12 0.38 0.12] 3.5) -150 0)]}))

(defn- heap
  "Ten dummies dropped in rows of five, as the demo drops them, under the
  sequential impulse."
  [opts]
  (w/world [floor]
           (vec (for [i (range 10)
                      :let [col (mod i 5) row (quot i 5)]]
                  {:model dummy
                   :pose {:q (mapv (fn [l] (if (= :spherical (:joint l)) idq 0.0)) (:links dummy))
                          :qd (vec (repeat (ab/dof dummy) 0.0))
                          :base {:rot (q/mul (q/from-axis-angle [0.0 1.0 0.0] (* 0.7 i))
                                             (q/from-axis-angle [0.0 0.0 1.0] 0.15))
                                 :pos [(* 0.75 (- col 2)) (+ 1.8 (* 0.35 row)) (* 0.75 (- row 1))]
                                 :vel [0.0 2.0 0.0 0.0 0.0 0.0]}}}))
           (merge {:gravity g :friction 0.7 :solver :sequential-impulse} opts)))

(deftest heap-sleeps-test
  (testing "a heap of ragdolls goes to sleep, given its limits' slop"
    ;; A limb resting on its limit with none is pushed just inside, falls
    ;; back past and is pushed again; one pressed past it against its own
    ;; chest by the heap fights the chest for good. Three degrees of slop
    ;; and ten of them sleep in eleven seconds; with a hundredth of a
    ;; radian, in thirty.
    (let [after (run (heap {:limit-slop 0.05}) 1200 (/ 1.0 60.0))]
      (is (every? :sleeping? (:models after))
          (str (count (filter :sleeping? (:models after))) " of 10 asleep")))))

(deftest limit-warm-start-test
  (testing "a limit pushed against remembers what it pushed with, for the next step"
    ;; Held at the bottom of [0.3, 0.6] by its own weight: two kilograms
    ;; at half a meter, so 2 g 0.5 sin 0.3 newton meters, a sixtieth of a
    ;; second of it each step.
    (doseq [solver [:convex :sequential-impulse]]
      (let [limb (assoc-in (hanging 0.0) [0 :limit] [0.3 0.6])
            end (run (w/world [] [{:model limb :pose {:q [0.35] :qd [0.0]}}]
                              {:gravity g :solver solver :allow-sleep? false})
                     240 (/ 1.0 60.0))
            held (get (:limit-impulses end) [0 0 :lo])
            weight (* 2.0 9.81 0.5 (Math/sin 0.3) (/ 1.0 60.0))]
        (is (some? held) (str solver " remembered " (:limit-impulses end)))
        (is (< (abs (- (double held) weight)) (* 0.1 weight)) (str solver " " held " " weight))))))

(defn- free-box
  "A loose box as an articulated model: a free root with nothing hung off
  it, so it collides as a model rather than as a rigid body."
  []
  {:base {:mass 2.0 :com [0.0 0.0 0.0] :inertia (box-inertia [0.3 0.3 0.3] 2.0)
          :shape :box :size [0.3 0.3 0.3]
          :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
   :links []})

(defn- box-pose [x vx]
  {:q [] :qd [] :base {:rot idq :pos [(double x) 0.0 0.0] :vel [0.0 0.0 0.0 (double vx) 0.0 0.0]}})

(defn- collide [dt opts]
  (run (w/world [] [{:model (free-box) :pose (box-pose -0.6 2.0)}
                    {:model (free-box) :pose (box-pose 0.0 0.0)}]
                (merge {:gravity [0.0 0.0 0.0]} opts))
       (long (/ 1.0 dt)) dt))

(defn- model-momentum [world]
  (reduce v/add (map (fn [{:keys [model pose]}] (:linear (ab/momentum model pose)))
                     (:models world))))

(deftest models-touch-test
  (testing "one model runs into another and they part"
    ;; No gravity and no floor: a box model flying at a resting one. If
    ;; models did not see each other the first would pass straight
    ;; through.
    (let [end (collide (/ 1.0 240.0) nil)
          x (fn [i] (double (first (:pos (:base (:pose (nth (:models end) i)))))))
          vx (fn [i] (let [{:keys [rot vel]} (:base (:pose (nth (:models end) i)))]
                       (double (first (q/rotate rot (subvec (vec vel) 3 6))))))]
      (is (pos? (vx 1)) "the resting one was pushed")
      (is (< (x 0) (x 1)) "and the flying one did not pass through it")))
  (testing "the push between them loses momentum only as fast as the integrator does"
    ;; Equal and opposite impulses cannot move the pair's total; the
    ;; first-order step underneath can, a little -- so the test, as for a
    ;; model against a rigid body, is that a finer step loses less.
    (let [start (v/scale [2.0 0.0 0.0] 2.0)
          drift (fn [dt] (v/distance start (model-momentum (collide dt nil))))
          coarse (drift (/ 1.0 240.0))
          fine (drift (/ 1.0 960.0))]
      (is (< coarse 1e-3) (str "coarse " coarse))
      ;; Or none at all: the convex solve's one exact Newton step on a
      ;; single contact can leave the pair's total exact to the bit.
      (is (<= fine coarse) (str "coarse " coarse " fine " fine))))
  (testing "told not to, they pass through each other"
    (let [end (collide (/ 1.0 240.0) {:collide-models? false})]
      (is (> (double (first (:pos (:base (:pose (first (:models end))))))) 1.0)))))

(deftest warm-start-test
  (let [column (fn [warm?]
                 (w/world (into [floor] (for [r (range 4)] (brick [0.0 (+ 0.125 (* 0.25 r)) 0.0])))
                          [] {:gravity g :solver :sequential-impulse :warm-start? warm?}))
        top (fn [s] (double (second (:pos (last (:bodies s))))))
        sink (fn [warm?] (- 0.875 (top (run (column warm?) 300 (/ 1.0 60.0)))))]
    (testing "each contact remembers what it pushed with"
      ;; A brick on the floor, settled: the normal impulses it is held
      ;; up with add to its weight over the step.
      (let [s (run (w/world [floor (brick [0.0 0.125 0.0])] [] {:gravity g}) 120 (/ 1.0 60.0))
            m (double (:mass (second (:bodies s))))
            held (reduce + (map first (vals (:contact-impulses s))))]
        (is (= 4 (count (:contact-impulses s))))
        (is (close? held (* m 9.81 (/ 1.0 60.0)) (* 0.01 held)) (str held))))
    (testing "and starting from it, a column of four sags less than starting cold"
      ;; Cold, eight iterations have to rediscover the column's weight
      ;; every step and never quite do.
      (is (< (sink true) (* 0.5 (sink false))) (str (sink true) " " (sink false))))))

(deftest soft-contacts-test
  (let [;; three bricks started 12 cm into the floor
        sunk (fn [soft?]
               (w/world (into [floor] (for [r (range 3)] (brick [(* 0.3 r) (- 0.125 0.12) 0.0])))
                        [] {:gravity g :solver :sequential-impulse :soft? soft?}))
        fastest (fn [s] (reduce (fn [[m s] _]
                                  (let [s (w/step s (/ 1.0 60.0))]
                                    [(apply max m (map #(double (second (:vel %))) (rest (:bodies s)))) s]))
                                [0.0 s] (range 120)))
        [hard-speed hard] (fastest (sunk false))
        [soft-speed soft] (fastest (sunk true))
        rest-at (fn [s] (double (second (:pos (second (:bodies s))))))]
    (testing "a soft contact pushes a deep overlap out gently, and neither overshoots"
      (is (< soft-speed (* 0.5 hard-speed)) (str soft-speed " " hard-speed))
      (is (< (rest-at hard) 0.125))
      (is (< (rest-at soft) 0.125)))
    (testing "and comes to rest a little deeper, as a spring holding the weight"
      (is (< 0.0 (- (rest-at hard) (rest-at soft)) 0.003) (str (rest-at hard) " " (rest-at soft))))))

(deftest recall-test
  (let [column (fn [reach]
                 (w/world (into [floor] (for [r (range 6)] (brick [0.0 (+ 0.125 (* 0.25 r)) 0.0])))
                          [] {:gravity g :solver :sequential-impulse :recall-reach reach}))
        drift (fn [s] (apply max (map (fn [b] (let [[x _ z] (:pos b)] (Math/hypot x z))) (rest (:bodies s)))))]
    (testing "a column of six stands when a renamed contact keeps its impulse, and falls when it loses it"
      ;; The bricks sway by fractions of a millimeter, their corners cross
      ;; the edges of the bricks below, and the box test renames them.
      (is (< (drift (run (column 0.01) 600 (/ 1.0 60.0))) 0.05))
      (is (> (drift (run (column nil) 600 (/ 1.0 60.0))) 0.5)))))

(def ^:private hinged
  "Two boxes on a hinge, free to fall."
  {:base {:mass 5.0 :com [0.0 0.0 0.0] :inertia (box-inertia [0.4 0.1 0.2] 5.0)
          :shape :box :size [0.4 0.1 0.2] :shape-pose {:rot idq :pos [0.0 0.0 0.0]}}
   :links [{:parent -1 :joint :revolute :axis [0.0 0.0 1.0] :origin {:rot nil :pos [0.2 0.0 0.0]}
            :mass 2.0 :com [0.2 0.0 0.0] :inertia (box-inertia [0.4 0.1 0.2] 2.0)
            :shape :box :size [0.4 0.1 0.2] :shape-pose {:rot idq :pos [0.2 0.0 0.0]}}]})

(defn- dropped [bodies]
  (w/world (into [floor] bodies)
           [{:model hinged
             :pose {:q [0.3] :qd [0.0]
                    :base {:rot idq :pos [0.0 0.5 0.0] :vel [0.0 0.0 0.0 0.0 0.0 0.0]}}}]
           {:gravity g}))

(defn- frames [s n] (vec (take n (iterate #(w/step % (/ 1.0 60.0)) s))))

(deftest sleeping-test
  (testing "a column at rest falls asleep, all of it at once, and remembers what held it up"
    (let [s (run (w/world (into [floor] (for [r (range 4)] (brick [0.0 (+ 0.125 (* 0.25 r)) 0.0])))
                          [] {:gravity g})
                 180 (/ 1.0 60.0))
          bricks (rest (:bodies s))]
      (is (every? rigid/sleeping? bricks))
      (is (= 1 (count (set (map :island bricks)))))
      (is (= 16 (count (:contact-impulses s))))
      (testing "and the step a ball reaches it, the whole column wakes"
        (let [ball (assoc (rigid/ball {:pos [-1.0 0.875 0.0] :radius 0.1 :density 1000.0})
                          :vel [6.0 0.0 0.0])
              fs (frames (update s :bodies conj ball) 30)
              touched (first (keep-indexed (fn [i s] (when (seq (filter #(= [:rigid 5] (second %))
                                                                        (keys (:contact-impulses s))))
                                                       i))
                                           fs))]
          (is (some? touched))
          (is (not-any? rigid/sleeping? (rest (:bodies (nth fs touched)))))))))
  (testing "a model at rest falls asleep, a brick dropped on it wakes it, and so does `wake`"
    (let [fs (frames (dropped []) 200)
          asleep (first (keep-indexed (fn [i s] (when (:sleeping? (first (:models s))) i)) fs))
          s (nth fs (inc asleep))]
      (is (< 30 asleep 120) (str asleep))
      (is (not (:sleeping? (first (:models (w/wake s [:link 0 0]))))))
      (let [brick (rigid/box {:pos [0.0 0.5 0.0] :size [0.1 0.1 0.1] :density 2000.0})
            fs (frames (update s :bodies conj brick) 40)]
        (is (some #(not (:sleeping? (first (:models %)))) fs)))))
  (testing "a model with a torque on it never sleeps"
    (is (not-any? #(:sleeping? (first (:models %)))
                  (frames (assoc-in (dropped []) [:models 0 :tau] [0.01]) 300))))
  (testing "a brick resting on a model sleeps with it, in one island"
    (let [s (run (dropped [(rigid/box {:pos [-0.05 0.62 0.0] :size [0.1 0.1 0.1] :density 500.0})])
                 500 (/ 1.0 60.0))
          m (first (:models s))
          b (second (:bodies s))]
      (is (and (:sleeping? m) (rigid/sleeping? b)))
      (is (= (:island m) (:island b))))))

(deftest convex-test
  (testing "every constraint at once: a column of thirty cubes stands, a millimeter a cube shorter,
            where the sequential impulse holds six"
    (let [s (run (w/world (into [floor] (for [r (range 30)] (brick [0.0 (+ 0.125 (* 0.25 r)) 0.0])))
                          [] {:gravity g :solver :convex :allow-sleep? false})
                 600 (/ 1.0 60.0))
          top (last (:bodies s))
          [x y z] (:pos top)]
      (is (< (Math/hypot x z) 0.01))
      (is (< (- (+ 0.125 (* 0.25 29)) (double y)) 0.04))))
  (testing "and the pins, joint limits and models touching do under the sequential impulse what they do under it"
    (with-redefs [w/defaults (assoc w/defaults :solver :sequential-impulse)]
      (doseq [t [#'model-pushes-rigid-test #'pin-to-the-world-test #'pin-between-models-test #'pin-breaks-test
                 #'limits-in-a-world-test #'contact-kinds-test #'model-alone-test #'sleeping-test]]
        ((:test (meta t)))))))
