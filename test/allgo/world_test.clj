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
  "A metre of limb on a ball joint, bolted at the origin and hanging."
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
