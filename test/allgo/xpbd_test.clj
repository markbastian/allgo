(ns allgo.xpbd-test
  (:require [allgo.geometry.tet-mesh :as tm]
            [allgo.physics.xpbd :as x]
            [clojure.test :refer [deftest is testing]]))

(def ^:private no-gravity {:gravity [0.0 0.0 0.0] :floor nil})

(defn- finite? [v] (and (not (Double/isNaN v)) (not (Double/isInfinite v))))

(defn- all-finite? [body] (every? finite? (x/positions body)))

(defn- ys [body] (map second (partition 3 (x/positions body))))

(defn- centre-of-mass
  "Mass-weighted, so it is the quantity momentum actually conserves."
  [body]
  (let [n  (:n body)
        ms (for [i (range n)
                 :let [w (x/inv-mass body i)]]
             (if (pos? w) (/ 1.0 w) 0.0))
        total (reduce + ms)]
    (mapv (fn [axis]
            (/ (reduce + (map (fn [i m] (* m (nth (x/particle body i) axis)))
                              (range n) ms))
               total))
          [0 1 2])))

(defn- separation
  "The distance between the two particles. Distance, not signed
  difference: a distance constraint says nothing about which side of each
  other they are on, and they can end up swapped."
  [body]
  (abs (- (first (x/particle body 1)) (first (x/particle body 0)))))

(defn- two-particle-spring
  "Two particles a unit apart at rest, stretched to `stretched`."
  [stretched compliance]
  (let [mesh {:verts [0.0 0.0 0.0 1.0 0.0 0.0] :edge-ids [0 1]}
        body (-> (x/body [0.0 0.0 0.0 stretched 0.0 0.0])
                 (x/add-constraint (x/distance-constraint mesh compliance)))]
    (x/set-inv-mass! body 0 1.0)
    (x/set-inv-mass! body 1 1.0)
    body))

;; ---------------------------------------------------------------------------

(deftest free-fall-test
  (testing "with nothing constraining it, a particle falls under gravity"
    (let [body (x/body [0.0 100.0 0.0])]
      (x/set-inv-mass! body 0 1.0)
      (dotimes [_ 60] (x/step! body {:floor nil}))
      (let [[_ y _] (x/particle body 0)]
        ;; After one second of g = 10, analytically 100 - 5 = 95. Symplectic
        ;; Euler lands half a substep out, which at 600 substeps is small.
        (is (< (abs (- 95.0 y)) 0.02)))))

  (testing "velocity is read back from the motion, not integrated separately"
    (let [body (x/body [0.0 100.0 0.0])]
      (x/set-inv-mass! body 0 1.0)
      (dotimes [_ 60] (x/step! body {:floor nil}))
      (let [[_ vy _] (let [v (:vel body)] [(aget v 0) (aget v 1) (aget v 2)])]
        (is (< (abs (- -10.0 vy)) 0.05)))))

  (testing "a particle with zero inverse mass is immovable"
    (let [body (x/body [0.0 5.0 0.0])]
      (dotimes [_ 60] (x/step! body {}))
      (is (= [0.0 5.0 0.0] (x/particle body 0))))))

(deftest distance-constraint-test
  (testing "one projection satisfies a rigid constraint exactly"
    ;; Measured on a single substep, which isolates the projection from
    ;; what happens to the velocity it produces.
    (doseq [stretched [3.0 1.5 0.25 0.01]]
      (let [body (two-particle-spring stretched 0.0)]
        (x/substep! body no-gravity (/ 1.0 600))
        (is (< (abs (- 1.0 (separation body))) 1e-9) (str "from " stretched)))))

  (testing "and it settles there once the energy is allowed to leave"
    ;; Undamped, the correction is handed to the velocity and the pair
    ;; oscillates about the rest length forever -- from far enough out
    ;; they overshoot hard enough to pass through each other and settle
    ;; the same distance apart on opposite sides. Damping is what makes
    ;; "settles" meaningful.
    (doseq [stretched [3.0 0.25]]
      (let [body (two-particle-spring stretched 0.0)]
        (dotimes [_ 30] (x/step! body (assoc no-gravity :damping 4.0)))
        (is (< (abs (- 1.0 (separation body))) 1e-6) (str "from " stretched)))))

  (testing "equal masses are corrected symmetrically"
    (let [body (two-particle-spring 3.0 0.0)]
      (x/step! body no-gravity)
      (let [[x0] (x/particle body 0) [x1] (x/particle body 1)]
        ;; They started at 0 and 3, so the midpoint is 1.5 and both must
        ;; have moved the same distance toward it.
        (is (< (abs (- (- 1.5 x0) (- x1 1.5))) 1e-9)))))

  (testing "a heavier particle moves less"
    (let [body (two-particle-spring 3.0 0.0)]
      (x/set-inv-mass! body 0 0.25)          ; four times the mass of the other
      (x/set-inv-mass! body 1 1.0)
      (x/step! body no-gravity)
      (let [moved0 (abs (- (first (x/particle body 0)) 0.0))
            moved1 (abs (- (first (x/particle body 1)) 3.0))]
        (is (< (abs (- (* 4 moved0) moved1)) 1e-9)
            "displacement is in proportion to inverse mass"))))

  (testing "compliance is softness: more of it corrects less"
    (let [errors (for [c [0.0 1e-5 1e-4 1e-3]]
                   (let [body (two-particle-spring 2.0 c)]
                     (x/substep! body no-gravity (/ 1.0 600))
                     (abs (- 1.0 (separation body)))))]
      (is (apply < errors) "each softer material ends further from rest")
      (is (< (first errors) 1e-9) "and compliance 0 is rigid"))))

(deftest volume-constraint-test
  (let [mesh (tm/lattice-box 1 1 1 1.0)]
    (testing "a squashed tetrahedral block recovers its volume"
      (let [body (x/soft-body mesh)
            v0   (x/mesh-volume body (:tet-ids mesh))]
        ;; Flatten it: every particle onto one plane, a total collapse.
        (dotimes [i (:n body)]
          (let [[px _ pz] (x/particle body i)]
            (x/set-particle! body i [px 0.5 pz])))
        (is (< (x/mesh-volume body (:tet-ids mesh)) (* 0.01 v0)) "it really was flattened")
        (dotimes [_ 120] (x/step! body no-gravity))
        (is (all-finite? body))
        (is (> (x/mesh-volume body (:tet-ids mesh)) (* 0.5 v0))
            "and it inflates again rather than staying flat")))

    (testing "volume is held through a long fall and landing"
      (let [mesh (tm/lattice-box 2 2 2 0.2)
            body (x/soft-body (tm/translate mesh [0.0 2.0 0.0]))
            v0   (x/mesh-volume body (:tet-ids mesh))]
        (dotimes [_ 240] (x/step! body {}))
        (let [v1 (x/mesh-volume body (:tet-ids mesh))]
          (is (< (abs (/ (- v1 v0) v0)) 0.02) "within 2% after four seconds"))))))

(deftest unbreakable-test
  ;; The claim the tutorial is named for. Position based dynamics never
  ;; integrates a force, so there is no configuration that can diverge:
  ;; the worst case is a bad position, which the next projection corrects.
  (let [mesh (tm/lattice-box 2 2 2 0.2)]
    (testing "surviving a total collapse to a point"
      (let [body (x/soft-body mesh)]
        (dotimes [i (:n body)] (x/set-particle! body i [0.0 0.5 0.0]))
        (dotimes [_ 180] (x/step! body {}))
        (is (all-finite? body))))

    (testing "surviving absurd velocities"
      (let [body (x/soft-body (tm/translate mesh [0.0 1.0 0.0]))]
        (dotimes [i (:n body)] (x/set-velocity! body i [1e4 -1e4 1e4]))
        (dotimes [_ 120] (x/step! body {}))
        (is (all-finite? body))))

    (testing "surviving an inside-out start"
      (let [body (x/soft-body mesh)]
        (dotimes [i (:n body)]
          (let [[px py pz] (x/particle body i)]
            (x/set-particle! body i [(- px) (+ 0.4 (- py)) pz])))
        (dotimes [_ 180] (x/step! body {}))
        (is (all-finite? body))))

    (testing "surviving a huge timestep, where a force-based solver would explode"
      (let [body (x/soft-body (tm/translate mesh [0.0 3.0 0.0]))]
        (dotimes [_ 60] (x/step! body {:dt 1.0 :substeps 5}))
        (is (all-finite? body))))))

(deftest floor-test
  (let [mesh (tm/lattice-box 2 2 2 0.2)]
    (testing "a body dropped on the floor stays above it and settles"
      (let [body (x/soft-body (tm/translate mesh [0.0 1.5 0.0]))]
        (dotimes [_ 300] (x/step! body {}))
        (is (all-finite? body))
        ;; The floor is applied to the predicted position, and the
        ;; constraint solve runs after it, so a particle can be pushed a
        ;; little way back through. On a 0.2 cell this is a thousandth of
        ;; an element -- invisible, but not zero, and worth stating.
        (is (>= (apply min (ys body)) -1e-3) "nothing ends meaningfully below the floor")
        (is (< (apply max (ys body)) 1.5) "and it did fall")))

    (testing "the floor can be moved, or removed"
      (let [body (x/soft-body (tm/translate mesh [0.0 3.0 0.0]))]
        (dotimes [_ 180] (x/step! body {:floor 2.0}))
        (is (>= (apply min (ys body)) (- 2.0 1e-3))))
      (let [body (x/soft-body (tm/translate mesh [0.0 1.0 0.0]))]
        (dotimes [_ 60] (x/step! body {:floor nil}))
        (is (< (apply min (ys body)) 0.0) "without a floor it keeps going")))))

(deftest momentum-test
  (testing "with no gravity and no floor, the centre of mass drifts in a straight line"
    (let [mesh (tm/lattice-box 2 2 2 0.2)
          body (x/soft-body mesh)]
      (dotimes [i (:n body)] (x/set-velocity! body i [1.0 0.0 0.0]))
      (let [c0 (centre-of-mass body)
            _  (dotimes [_ 60] (x/step! body no-gravity))
            c1 (centre-of-mass body)]
        ;; One second at one unit per second.
        (is (< (abs (- 1.0 (- (c1 0) (c0 0)))) 0.02))
        (is (< (abs (- (c1 1) (c0 1))) 1e-6) "and does not wander off-axis")
        (is (< (abs (- (c1 2) (c0 2))) 1e-6))))))

(deftest substep-test
  (testing "substeps buy accuracy: more of them hold the volume better"
    (let [mesh  (tm/lattice-box 2 2 2 0.2)
          drift (fn [substeps]
                  (let [body (x/soft-body (tm/translate mesh [0.0 1.0 0.0])
                                          {:edge-compliance 1e-4})
                        v0   (x/mesh-volume body (:tet-ids mesh))]
                    (dotimes [_ 120] (x/step! body {:substeps substeps}))
                    (abs (/ (- (x/mesh-volume body (:tet-ids mesh)) v0) v0))))]
      (is (< (drift 20) (drift 1))))))

(deftest grab-test
  (let [mesh (tm/lattice-box 2 2 2 0.2)]
    (testing "grabbing pins the nearest particle and dragging carries the body"
      (let [body       (x/soft-body mesh)
            [body id]  (x/grab! body [0.2 0.4 0.2])]
        (is (zero? (x/inv-mass body id)) "the held particle is immovable")
        (x/move-grabbed! body [0.5 2.0 0.5])
        (is (= [0.5 2.0 0.5] (x/particle body id)))
        (dotimes [_ 60] (x/step! body {}))
        (is (= [0.5 2.0 0.5] (x/particle body id)) "and stays where it is put")
        (is (all-finite? body))))

    (testing "releasing gives the mass back"
      (let [body      (x/soft-body mesh)
            before    (x/inv-mass body 0)
            [body id] (x/grab! body (x/particle body 0))
            body      (x/release! body [1.0 0.0 0.0])]
        (is (= 0 id))
        (is (= before (x/inv-mass body id)))
        (is (nil? (:grabbed body)))
        (dotimes [_ 30] (x/step! body {}))
        (is (all-finite? body))))))

(deftest soft-body-test
  (let [mesh (tm/lattice-box 2 2 2 0.2)
        body (x/soft-body mesh)]
    (testing "it is built with both constraint kinds"
      (is (= 2 (count (:constraints body)))))

    (testing "mass follows the mesh, and every particle has some"
      (is (every? #(pos? (x/inv-mass body %)) (range (:n body))))
      (testing "corner particles are lighter than interior ones"
        ;; A corner belongs to fewer tetrahedra, so it carries less mass
        ;; and has the larger inverse mass.
        (let [corner (apply max (map #(x/inv-mass body %) (range (:n body))))
              middle (apply min (map #(x/inv-mass body %) (range (:n body))))]
          (is (> corner middle)))))

    (testing "density scales mass, not geometry"
      (let [heavy (x/soft-body mesh {:density 4000.0})]
        (is (< (abs (- (* 4 (x/inv-mass heavy 0)) (x/inv-mass body 0))) 1e-9))
        (is (= (x/positions body) (x/positions heavy)))))))
