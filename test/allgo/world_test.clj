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
      (is (> (double (first (:pos b))) 0.45) (str "brick at " (:pos b)))
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
      (is (< fine coarse) (str "coarse " coarse " fine " fine))))
  (testing "told not to, they pass through each other"
    (let [end (collide (/ 1.0 240.0) {:collide-models? false})]
      (is (> (double (first (:pos (:base (:pose (first (:models end))))))) 1.0)))))
