(ns allgo.articulated-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.numerics.linear :as lin]
            [allgo.physics.articulated :as ab]
            [allgo.physics.rigid :as rigid]
            [clojure.test :refer [deftest is testing]]))

(def ^:private g [0.0 -9.81 0.0])
(def ^:private opts {:gravity g})

(defn- rod
  "A uniform rod of length `L` and mass `m`, lying along +x from its own
  joint frame and hinged about z. `q` is measured up from +x."
  [parent L m]
  (let [i (* (/ 1.0 12.0) m L L)]
    {:parent parent :joint :revolute :axis [0.0 0.0 1.0]
     :origin {:rot nil :pos (if (neg? parent) [0.0 0.0 0.0] [L 0.0 0.0])}
     :mass m :com [(* 0.5 L) 0.0 0.0]
     :inertia [[0.0 0.0 0.0] [0.0 i 0.0] [0.0 0.0 i]]}))

(defn- chain [n L m] (mapv #(rod (dec (long %)) L m) (range n)))

(defn- close? [^double a ^double b ^double tol] (< (abs (- a b)) tol))

(deftest spatial-algebra-test
  (testing "the force cross product is the negative transpose of the motion one"
    ;; Not a definition to be taken on trust: it is what makes the inward
    ;; pass of both algorithms a transpose rather than a second recursion.
    (let [v [0.3 -1.1 0.7 2.0 0.5 -0.9]]
      (is (= (ab/crf v) (lin/mat-scale (lin/transpose (ab/crm v)) -1.0)))))

  (testing "a transform that moves nothing is the identity"
    (is (= (lin/eye 6) (ab/transform {:rot nil :pos [0.0 0.0 0.0]}))))

  (testing "skew is the cross product written as a matrix"
    (let [a [0.2 -0.6 1.4] b [1.0 0.3 -0.8]]
      (is (every? #(< (abs %) 1e-15)
                  (map - (lin/mat-vec (ab/skew a) b)
                       [(- (* -0.6 -0.8) (* 1.4 0.3))
                        (- (* 1.4 1.0) (* 0.2 -0.8))
                        (- (* 0.2 0.3) (* -0.6 1.0))])))))

  (testing "a spatial inertia carries the mass in its lower-right block"
    (let [m 2.5
          i (ab/spatial-inertia m [0.4 0.0 0.0] (lin/eye 3))]
      (is (= m (get-in i [3 3])))
      (is (= m (get-in i [4 4])))
      (is (= m (get-in i [5 5]))))))

(deftest closed-form-test
  (testing "a hinged rod accelerates as the textbook says it does"
    ;; The one place there is an answer to check against that shares no
    ;; code with the thing being checked. A rod pivoted at one end:
    ;; qdd = -(m g d / I_pivot) cos q, with d the distance to the center
    ;; of mass and I_pivot = m L^2 / 3 by the parallel axis theorem.
    (let [L 1.3 m 2.4
          model [(rod -1 L m)]
          i-pivot (* (/ 1.0 3.0) m L L)]
      (doseq [q [0.0 0.3 1.0 -0.7 1.5707963267948966 3.0]]
        (let [qdd (first (:qdd (ab/forward-dynamics model [q] [0.0] [0.0] opts)))
              exact (/ (* (- (* m 9.81 (* 0.5 L))) (Math/cos q)) i-pivot)]
          (is (close? qdd exact 1e-12) (str "at q=" q))))))

  (testing "a mass on a frictionless slider takes gravity's component along it"
    (let [axis (mapv #(/ (double %) (Math/sqrt 2.0)) [1.0 -1.0 0.0])
          model [{:parent -1 :joint :prismatic :axis axis
                  :origin {:rot nil :pos [0.0 0.0 0.0]}
                  :mass 3.1 :com [0.0 0.0 0.0] :inertia (lin/eye 3)}]]
      (is (close? (first (:qdd (ab/forward-dynamics model [0.0] [0.0] [0.0] opts)))
                  (* 9.81 (- (double (second axis))))
                  1e-12))))

  (testing "a fixed joint has nothing to move"
    ;; Not a joint with one coordinate that stays put -- no coordinate.
    ;; A weld is there to hold two shapes together as one body, and
    ;; giving it a degree of freedom to ignore means every pass carries
    ;; a column of zeroes around.
    (let [model [(assoc (rod -1 1.0 1.0) :joint :fixed)]]
      (is (zero? (ab/dof model)))
      (is (= [] (:qdd (ab/forward-dynamics model [0.0] [] [] opts)))))))

(deftest forward-against-inverse-test
  (testing "the fast path agrees with M inverse times what is left"
    ;; The articulated body algorithm never forms the mass matrix. The
    ;; recursive Newton-Euler algorithm can be made to hand one over a
    ;; column at a time. Solving with it must give what the fast path
    ;; gave, and the two share only the spatial algebra underneath.
    (doseq [n [1 2 3 5 8]]
      (let [model (chain n 0.8 1.7)
            q (mapv #(* 0.37 (inc (long %))) (range n))
            qd (mapv #(* -0.23 (inc (long %))) (range n))
            tau (mapv #(* 0.9 (- 2 (long %))) (range n))
            fast (:qdd (ab/forward-dynamics model q qd tau opts))
            slow (lin/mat-vec (lin/inverse (ab/mass-matrix model q))
                              (mapv - tau (ab/bias-forces model q qd opts)))]
        (is (every? #(< (abs (double %)) 1e-9) (map - fast slow))
            (str "n=" n " " fast " vs " slow)))))

  (testing "the mass matrix is symmetric"
    (let [m (ab/mass-matrix (chain 4 0.8 1.7) [0.2 -0.5 1.1 0.3])]
      (is (every? #(< (abs (double %)) 1e-12)
                  (flatten (lin/mat-sub m (lin/transpose m)))))))

  (testing "and positive definite, so the chain has a unique acceleration"
    ;; Cholesky exists exactly when it is, which is the cheapest way to
    ;; ask and the one that fails loudly if a link is given no inertia.
    (is (some? (lin/cholesky (ab/mass-matrix (chain 4 0.8 1.7) [0.2 -0.5 1.1 0.3]))))))

(deftest kinematics-test
  (testing "joint angles turn back into places"
    (let [model (chain 2 1.0 1.0)
          origin-of (fn [q i] (mapv #(Math/round (* 1000.0 (double %)))
                                    (:pos (nth (ab/poses model q) i))))]
      (is (= [1000 0 0] (origin-of [0.0 0.0] 1)))
      (is (= [0 1000 0] (origin-of [(/ Math/PI 2) 0.0] 1)))
      ;; The second joint sits at the first link's tip whatever it does
      ;; itself, so bending it does not move its own frame.
      (is (= [1000 0 0] (origin-of [0.0 (/ Math/PI 2)] 1))))))

(deftest energy-test
  (testing "a chain under no torque loses energy only as fast as the integrator does"
    ;; Semi-implicit Euler is first order, so the drift over a fixed
    ;; stretch of time must fall off in step with the step size. That it
    ;; does is a statement about the dynamics; that it is not zero is a
    ;; statement about the integrator, and the two are worth keeping
    ;; apart.
    (doseq [n [1 2]]
      (let [model (chain n 0.8 1.7)
            q0 (mapv #(+ 0.4 (* 0.3 (double %))) (range n))
            e0 (ab/energy model {:q q0 :qd (vec (repeat n 0.0))} opts)
            drift (fn [h]
                    (let [final (reduce (fn [st _] (ab/step st model h opts))
                                        {:q q0 :qd (vec (repeat n 0.0))}
                                        (range (long (/ 2.0 h))))]
                      (abs (- (ab/energy model final opts) e0))))
            coarse (drift 1e-3)
            fine (drift 1e-4)]
        (is (pos? coarse))
        ;; Ten times the step, within a factor of two of ten times the
        ;; drift. Looser than first order exactly, because two links are
        ;; chaotic and the trajectories part company.
        (is (< 5.0 (/ coarse fine) 20.0)
            (str "n=" n " coarse " coarse " fine " fine))))))

;; ---------------------------------------------------------------------------
;; A base that is free to move

(defn- floating
  "A root body with `links` hanging off it."
  [links]
  {:base {:mass 3.0 :com [0.1 0.0 0.0]
          :inertia [[0.4 0.0 0.0] [0.0 0.6 0.0] [0.0 0.0 0.9]]}
   :links (vec links)})

(def ^:private at-rest
  {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0] :vel [0.0 0.0 0.0 0.0 0.0 0.0]})

(deftest floating-base-test
  (testing "a free body with nothing attached falls at g and nothing else"
    (let [{:keys [base-acc]} (ab/forward-dynamics (floating []) [] [] []
                                                  (assoc opts :base at-rest))]
      (is (every? #(< (abs (double %)) 1e-12)
                  (map - base-acc [0.0 0.0 0.0 0.0 -9.81 0.0])))))

  (testing "and falls the distance it should"
    (let [model (floating [])
          final (reduce (fn [st _] (ab/step st model 1e-4 opts))
                        {:q [] :qd [] :base at-rest}
                        (range 10000))]
      ;; Half g t squared, to the accuracy a first-order step has after
      ;; ten thousand of them.
      (is (close? (double (second (:pos (:base final)))) (* -0.5 9.81) 1e-3))))

  (testing "a bolted model reports no base acceleration at all"
    (is (= [0.0 0.0 0.0 0.0 0.0 0.0]
           (:base-acc (ab/forward-dynamics (chain 2 0.8 1.7) [0.1 0.2] [0.0 0.0]
                                           [0.0 0.0] opts))))))

(deftest momentum-test
  (testing "a floating chain with no gravity conserves both momenta"
    ;; The strongest thing available to say about a free multibody, and
    ;; the reason it is worth saying: momentum is a property of the whole
    ;; system rather than of any link, so a sign lost anywhere in the
    ;; recursion shows up here even when every individual link looks
    ;; plausible.
    (let [model (floating [(rod -1 0.7 1.2) (rod 0 0.6 0.9) (rod 1 0.5 0.6)])
          st0 {:q [0.3 -0.5 0.8] :qd [1.4 -2.1 0.7]
               :base {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6)
                      :pos [0.3 -0.2 0.5] :vel [0.4 -0.3 0.9 1.1 0.2 -0.6]}}
          free {:gravity [0.0 0.0 0.0]}
          m0 (ab/momentum model st0)
          m1 (ab/momentum model (reduce (fn [st _] (ab/step st model 1e-4 free))
                                        st0 (range 10000)))]
      (is (< (v/distance (:linear m0) (:linear m1)) 1e-2)
          (str (:linear m0) " -> " (:linear m1)))
      (is (< (v/distance (:angular m0) (:angular m1)) 1e-2)
          (str (:angular m0) " -> " (:angular m1)))))

  (testing "and with gravity on, linear momentum grows at exactly M g"
    ;; The check that catches gravity being applied in the wrong frame,
    ;; which is otherwise nearly invisible: sideways momentum is
    ;; conserved whatever the joints do, and a base that both falls and
    ;; spins stops conserving it the moment world-frame gravity is added
    ;; to a body-frame acceleration.
    (let [model (floating [(rod -1 0.7 1.2) (rod 0 0.6 0.9)])
          total (+ 3.0 1.2 0.9)
          st0 {:q [0.3 -0.5] :qd [1.4 -2.1]
               :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0]
                      :vel [0.2 0.1 -0.3 0.5 0.0 0.1]}}
          m0 (ab/momentum model st0)
          m1 (ab/momentum model (reduce (fn [st _] (ab/step st model 1e-4 opts))
                                        st0 (range 10000)))]
      (is (< (v/distance (mapv + (:linear m0) (mapv #(* (double %) total) g))
                         (:linear m1))
             1e-2)
          (str (:linear m0) " -> " (:linear m1))))))

;; ---------------------------------------------------------------------------
;; Being pushed

(defn- sphere-ish
  "A free body whose center of mass is at its frame origin and whose
  inertia is the same about every axis, so the textbook formulae for
  being hit apply without any further care."
  [m i]
  {:base {:mass m :com [0.0 0.0 0.0]
          :inertia [[i 0.0 0.0] [0.0 i 0.0] [0.0 0.0 i]]}
   :links []})

(deftest impulse-test
  (testing "through the center of mass, a free body weighs what it weighs"
    (let [m 2.5]
      (is (close? (:effective-mass (ab/impulse-at (sphere-ish m 0.7) [] at-rest -1
                                                  [0.0 0.0 0.0] [1.0 0.0 0.0]))
                  m 1e-12))))

  (testing "off the center of mass it weighs less, by the textbook amount"
    ;; 1/m_eff = 1/m + (r x d) . I^-1 (r x d). An impulse that can spin
    ;; the body as well as move it meets less resistance.
    (let [m 2.5 i 0.7
          r [0.0 0.4 0.0] dir [1.0 0.0 0.0]
          rn (v/cross r dir)
          exact (/ 1.0 (+ (/ 1.0 m) (/ (v/dot rn rn) i)))]
      (is (close? (:effective-mass (ab/impulse-at (sphere-ish m i) [] at-rest -1 r dir))
                  exact 1e-12))))

  (testing "at the tip of a hinged rod it weighs what the hinge leaves"
    ;; A transverse impulse P at the tip gives qd = P L / I_pivot, so the
    ;; mass felt there is I_pivot / L^2 -- less than the rod's own mass,
    ;; because the far end is what moves.
    (let [L 1.3 m 2.4
          i-pivot (* (/ 1.0 3.0) m L L)
          {:keys [effective-mass delta-u]}
          (ab/impulse-at [(rod -1 L m)] [0.0] nil 0 [L 0.0 0.0] [0.0 1.0 0.0])]
      (is (close? effective-mass (/ i-pivot (* L L)) 1e-12))
      (is (close? (first delta-u) (/ L i-pivot) 1e-12))))

  (testing "and along the rod it cannot be pushed at all"
    ;; The hinge takes the whole of it. An infinite effective mass is the
    ;; right answer and a contact solver has to be able to receive it.
    (let [{:keys [effective-mass delta-u]}
          (ab/impulse-at [(rod -1 1.3 2.4)] [0.0] nil 0 [1.3 0.0 0.0] [1.0 0.0 0.0])]
      (is (zero? (double (first delta-u))))
      (is (infinite? effective-mass))))

  (testing "the inverse mass matrix is the inverse of the mass matrix"
    ;; One comes from the articulated body algorithm without ever forming
    ;; a matrix, the other from inverse dynamics a column at a time.
    (doseq [n [1 2 4 6]]
      (let [model (chain n 0.8 1.7)
            q (mapv #(* 0.41 (inc (long %))) (range n))
            prod (lin/mat-mul (ab/inverse-mass-matrix model q) (ab/mass-matrix model q))]
        (is (every? #(< (abs (double %)) 1e-12)
                    (flatten (lin/mat-sub prod (lin/eye n))))
            (str "n=" n)))))

  (testing "an impulse on a free model changes its momentum by exactly that impulse"
    ;; The one that ties everything together. The Jacobian, the inverse
    ;; mass matrix and the momentum sum are three separate pieces of
    ;; arithmetic, and this holds only if all three agree -- for a
    ;; rotated base, an impulse on a middle link, and an oblique
    ;; direction, none of which are special-cased anywhere.
    (let [model (floating [(rod -1 0.7 1.1) (rod 0 0.7 1.1)])
          st {:q [0.3 -0.5] :qd [0.0 0.0]
              :base {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6)
                     :pos [0.3 -0.2 0.5] :vel [0.0 0.0 0.0 0.0 0.0 0.0]}}
          p [1.0 0.5 -0.2]
          dir (v/normalize [0.3 1.0 -0.4])
          mag 2.7
          m0 (ab/momentum model st)
          m1 (ab/momentum model (ab/apply-impulse model st 1 p dir mag))]
      (is (< (v/distance (v/sub (:linear m1) (:linear m0)) (v/scale dir mag)) 1e-12))
      (is (< (v/distance (v/sub (:angular m1) (:angular m0))
                         (v/cross p (v/scale dir mag)))
             1e-12)))))

;; ---------------------------------------------------------------------------
;; Shapes, and the ground

(def ^:private floor
  (rigid/box {:pos [0.0 -0.5 0.0] :size [40.0 1.0 40.0]}))

(defn- crate
  "A free box, its shape centered on its own frame rather than on a
  center of mass somewhere else."
  [[sx sy sz] m]
  (let [f (/ (double m) 12.0)]
    {:base {:mass m :com [0.0 0.0 0.0]
            :inertia [[(* f (+ (* sy sy) (* sz sz))) 0.0 0.0]
                      [0.0 (* f (+ (* sx sx) (* sz sz))) 0.0]
                      [0.0 0.0 (* f (+ (* sx sx) (* sy sy)))]]
            :shape :box :size [sx sy sz]
            :shape-pose {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0]}}
     :links []}))

(defn- shaped-rod
  "A rod that can also be collided with."
  [parent L m]
  (assoc (rod parent L m) :shape :box :size [L 0.12 0.12]))

(defn- fall
  "`n` steps of `dt`, under gravity, against `obstacles`."
  [model state n dt obstacles & [more]]
  (reduce (fn [s _] (ab/step s model dt (merge {:gravity g :obstacles obstacles} more)))
          state
          (range n)))

(defn- resting [y] {:q [] :qd [] :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 y 0.0]
                                        :vel [0.0 0.0 0.0 0.0 0.0 0.0]}})

(defn- base-y [s] (double (second (:pos (:base s)))))
(defn- base-speed [s] (v/length (subvec (vec (:vel (:base s))) 3 6)))

(deftest ground-contact-test
  (testing "a dropped box comes to rest on the floor and stays there"
    (let [model (crate [0.6 0.4 0.5] 2.0)
          landed (fall model (resting 2.0) 240 (/ 1.0 120.0) [floor])
          later (fall model landed 480 (/ 1.0 120.0) [floor])]
      ;; Half the box's height, less the slop a soft contact keeps.
      (is (close? (base-y landed) 0.2 0.01) (str "landed at " (base-y landed)))
      (is (< (base-speed landed) 0.01))
      (is (close? (base-y later) (base-y landed) 1e-3) "it did not creep")))

  (testing "and one that starts on the floor does not sink into it"
    (let [model (crate [0.6 0.4 0.5] 2.0)
          s (fall model (resting 0.2) 600 (/ 1.0 120.0) [floor])]
      (is (> (base-y s) 0.19) (str "sank to " (base-y s)))))

  (testing "nothing ends up buried"
    ;; Contacts are offered before the touch, so a settled body sits in
    ;; a gap rather than in a hole. Every depth here should be negative
    ;; or within the slop; a positive one means the solver let something
    ;; through and pushed it back out afterward.
    (let [model (crate [0.6 0.4 0.5] 2.0)
          s (fall model (resting 2.0) 360 (/ 1.0 120.0) [floor])
          depths (map :depth (ab/contacts-with model (:q s) (:base s) [floor]))]
      (is (seq depths) "it should be touching the floor at all")
      (is (every? #(< (double %) 0.01) depths) (str depths))))

  (testing "a link with no shape is not collided with"
    (let [bare (assoc-in (crate [0.6 0.4 0.5] 2.0) [:base :shape] nil)]
      (is (empty? (ab/collision-bodies bare [] (:base (resting 1.0)))))
      (is (empty? (ab/contacts-with bare [] (:base (resting 0.0)) [floor])))))

  (testing "the normal points the way the body has to be pushed"
    ;; Out of the floor, not into it. Stating it in the direction the
    ;; solver uses is what keeps the sign in one place.
    (let [model (crate [0.6 0.4 0.5] 2.0)
          cs (ab/contacts-with model [] (:base (resting 0.2)) [floor])]
      (is (seq cs))
      (is (every? #(> (double (second (:normal %))) 0.9) cs) (str (map :normal cs))))))

(deftest friction-test
  (testing "a box on a slope slides when it is slippery and grips when it is not"
    (let [tilt (q/from-euler 0.0 0.0 -0.35)
          ramp (rigid/box {:pos [0.0 -0.5 0.0] :size [40.0 1.0 40.0] :rot tilt})
          model (crate [0.6 0.4 0.5] 2.0)
          travel (fn [mu]
                   (let [st {:q [] :qd []
                             :base {:rot tilt :pos (q/rotate tilt [0.0 0.21 0.0])
                                    :vel [0.0 0.0 0.0 0.0 0.0 0.0]}}]
                     (abs (double (first (:pos (:base (fall model st 240 (/ 1.0 120.0)
                                                            [ramp] {:friction mu}))))))))]
      (is (> (travel 0.0) 2.0) "a frictionless box should be well down the slope")
      (is (< (travel 1.2) 0.2) "a gripping one should barely move")
      (is (> (travel 0.0) (travel 1.2))))))

(deftest bounce-test
  (testing "restitution returns a dropped box to the height it should"
    ;; Falls 0.8 m, so an ideal bounce comes back to 0.2 + 0.8 e^2.
    ;; Measured after the first descent bottoms out, because the drop
    ;; itself starts higher than any of these rebounds.
    (let [model (crate [0.4 0.4 0.4] 2.0)
          rebound (fn [e]
                    (let [ys (map base-y
                                  (take 600 (iterate #(ab/step % model (/ 1.0 240.0)
                                                               {:gravity g :obstacles [floor]
                                                                :restitution e})
                                                     (resting 1.0))))
                          low (apply min (take 130 ys))]
                      (apply max (take 400 (drop-while #(> (double %) (+ low 1e-9)) ys)))))]
      (doseq [e [0.0 0.5 0.9]]
        (is (close? (rebound e) (+ 0.2 (* 0.8 e e)) 0.02)
            (str "e=" e " came back to " (rebound e)))))))

(deftest chain-on-the-ground-test
  (testing "a chain on a free root falls over, settles, and stays settled"
    ;; The thing all of this was for. Nothing here is a special case:
    ;; the root is collided with like any link, the joints take whatever
    ;; the floor does to the links through the same inverse inertia, and
    ;; it comes to rest.
    (let [model {:base {:mass 2.0 :com [0.0 0.0 0.0]
                        :inertia [[0.03 0.0 0.0] [0.0 0.03 0.0] [0.0 0.0 0.03]]
                        :shape :box :size [0.3 0.3 0.3]}
                 :links [(shaped-rod -1 0.6 0.8) (shaped-rod 0 0.6 0.8)]}
          st0 {:q [-0.4 -0.7] :qd [0.0 0.0]
               :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 1.5 0.0]
                      :vel [0.0 0.0 0.0 0.0 0.0 0.0]}}
          landed (fall model st0 480 (/ 1.0 120.0) [floor])
          later (fall model landed 720 (/ 1.0 120.0) [floor])]
      (is (< (base-speed landed) 0.05) (str "still moving at " (base-speed landed)))
      (is (close? (base-y later) (base-y landed) 0.01) "it drifted after settling")
      (is (every? #(> (double (second (:pos %))) -0.05)
                  (ab/poses model (:q later) (:base later)))
          "a link ended up under the floor")
      ;; And the joints have stopped turning, which is the part a
      ;; constraint formulation has to be talked into.
      (is (every? #(< (abs (double %)) 0.05) (:qd later)) (str (:qd later))))))

(deftest impulse-response-matches-the-matrix-test
  (testing "the O(n) impulse response gives what the inverse inertia matrix gives"
    ;; Two entirely different ways to the same number, and the slow one
    ;; is the one with independent evidence behind it -- it is checked
    ;; against the mass matrix from inverse dynamics, which is checked
    ;; against the closed form. So this is what carries that evidence
    ;; over to the fast path the contact solver actually runs.
    (doseq [model [(chain 4 0.8 1.7)
                   (floating [(rod -1 0.7 1.1) (rod 0 0.6 0.9) (rod 1 0.5 0.7)])]]
      (let [n (ab/dof model)
            root? (some? (ab/base model))
            q (mapv #(* 0.37 (inc (long %))) (range n))
            root-state (when root?
                         {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6)
                          :pos [0.3 -0.2 0.5]
                          :vel [0.0 0.0 0.0 0.0 0.0 0.0]})
            hinv (ab/inverse-mass-matrix model q)]
        (doseq [i (if root? [-1 0 2] [0 2 3])
                dir [[1.0 0.0 0.0] [0.0 1.0 0.0] (v/normalize [0.3 -0.8 0.5])]]
          (let [p (v/add (:pos (ab/frame-of model q root-state i)) [0.11 -0.07 0.23])
                fast (:delta-u (ab/impulse-at model q root-state i p dir))
                ;; The same thing the slow way: J^T d through H^-1.
                jt (lin/transpose (ab/point-jacobian model q root-state i p))
                slow (lin/mat-vec hinv (lin/mat-vec jt (v/normalize dir)))]
            (is (every? #(< (abs (double %)) 1e-9) (map - fast slow))
                (str "body " i " along " dir "\n  fast " fast "\n  slow " slow))))))))

;; ---------------------------------------------------------------------------
;; Ball joints

(def ^:private limb-width 0.14)

(defn- limb
  "A rod with the inertia of a thin box rather than of a mathematical
  line. A ball joint needs the third one: a line has no inertia about
  its own length, and a joint free to turn about an axis nothing resists
  has no acceleration there to compute."
  [parent kind L m]
  (let [f (/ (double m) 12.0)
        w limb-width]
    (cond-> {:parent parent :joint kind
             :origin {:rot nil :pos (if (neg? (long parent)) [0.0 0.0 0.0] [L 0.0 0.0])}
             :mass m :com [(* 0.5 L) 0.0 0.0]
             :inertia [[(* f 2.0 w w) 0.0 0.0]
                       [0.0 (* f (+ (* L L) (* w w))) 0.0]
                       [0.0 0.0 (* f (+ (* L L) (* w w)))]]}
      (= kind :revolute) (assoc :axis [0.0 0.0 1.0]))))

(def ^:private idq [0.0 0.0 0.0 1.0])

(deftest joint-dof-test
  (testing "a joint is as wide as the motion it allows"
    (is (= 1 (ab/joint-dof {:joint :revolute :axis [0.0 0.0 1.0]})))
    (is (= 1 (ab/joint-dof {:joint :prismatic :axis [0.0 0.0 1.0]})))
    (is (= 3 (ab/joint-dof {:joint :spherical})))
    (is (= 0 (ab/joint-dof {:joint :fixed})))
    ;; And the model's is the sum, which is what `qd` and `tau` are long.
    (is (= 4 (ab/dof [(limb -1 :spherical 0.7 1.0) (limb 0 :revolute 0.6 1.0)])))))

(deftest ball-joint-test
  (testing "a ball joint swinging about one axis is a hinge about that axis"
    ;; The same rod, the same gravity, two different ways of saying
    ;; where it is. Nothing about the answer should depend on which.
    (let [L 1.3 m 2.4
          hinge [(limb -1 :revolute L m)]
          ball [(limb -1 :spherical L m)]
          run (fn [model st] (reduce (fn [s _] (ab/step s model 1e-4 opts)) st (range 20000)))
          hs (run hinge {:q [0.0] :qd [0.0]})
          bs (run ball {:q [idq] :qd [0.0 0.0 0.0]})
          ;; The ball joint's angle about z, read back off its quaternion.
          [bx by _] (q/rotate (first (:q bs)) [1.0 0.0 0.0])]
      (is (close? (double (first (:q hs))) (Math/atan2 (double by) (double bx)) 1e-6))
      (is (close? (double (first (:qd hs))) (double (nth (:qd bs) 2)) 1e-5))))

  (testing "and weighs the same as the hinge does, about the hinge's axis"
    (let [L 1.3 m 2.4]
      (is (close? (:effective-mass (ab/impulse-at [(limb -1 :spherical L m)] [idq] nil
                                                  0 [L 0.0 0.0] [0.0 1.0 0.0]))
                  (:effective-mass (ab/impulse-at [(limb -1 :revolute L m)] [0.0] nil
                                                  0 [L 0.0 0.0] [0.0 1.0 0.0]))
                  1e-9))))

  (testing "a ball joint adds rotations, not slides"
    ;; Pushing a rod along its own length through the joint still meets
    ;; an immovable object, because no rotation takes it that way.
    (is (infinite? (:effective-mass (ab/impulse-at [(limb -1 :spherical 1.3 2.4)] [idq] nil
                                                   0 [1.3 0.0 0.0] [1.0 0.0 0.0])))))

  (testing "there is no orientation it stops working in"
    ;; The reason the configuration is a quaternion and the velocity is
    ;; three numbers rather than both being three. Turn the whole
    ;; problem about the gravity axis and the answer must turn with it;
    ;; a three-angle parameterization has orientations where it does
    ;; not, and a ragdoll finds them by tumbling.
    (let [L 1.3 m 2.4
          model [(limb -1 :spherical L m)]
          swing (fn [r]
                  (let [st (reduce (fn [s _] (ab/step s model 1e-4 opts))
                                   {:q [r] :qd [0.0 0.0 0.0]} (range 3000))
                        {:keys [rot pos]} (first (ab/poses model (:q st)))]
                    (v/add pos (q/rotate rot [L 0.0 0.0]))))
          turn (q/from-axis-angle [0.0 1.0 0.0] (/ Math/PI 2))]
      (is (< (v/distance (swing turn) (q/rotate turn (swing idq))) 1e-9))))

  (testing "a ball joint on a body with no inertia about an axis welds itself"
    ;; Worth stating because the failure is quiet. A mathematically thin
    ;; rod has nothing to resist rotation about its own length, so the
    ;; three by three the joint has to invert is singular, and a joint
    ;; whose accelerations cannot be computed does not move at all.
    ;; Give limbs a real inertia; this is what happens if you do not.
    (let [thin (assoc (limb -1 :spherical 1.3 2.4)
                      :inertia [[0.0 0.0 0.0] [0.0 0.3 0.0] [0.0 0.0 0.3]])]
      (is (every? zero? (:qdd (ab/forward-dynamics [thin] [idq] [0.0 0.0 0.0]
                                                   [0.0 0.0 0.0] opts)))))))

(deftest ball-joint-invariants-test
  (testing "the fast impulse path still matches the matrix with ball joints in the chain"
    (let [model {:base {:mass 3.0 :com [0.1 0.0 0.0]
                        :inertia [[0.4 0.0 0.0] [0.0 0.6 0.0] [0.0 0.0 0.9]]}
                 :links [(limb -1 :spherical 0.7 1.2)
                         (limb 0 :revolute 0.6 0.9)
                         (limb 1 :spherical 0.5 0.6)]}
          qv [(q/from-axis-angle [0.3 0.5 0.8] 0.4) 0.37
              (q/from-axis-angle [0.9 0.1 0.2] -0.6)]
          rs {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6) :pos [0.3 -0.2 0.5]
              :vel [0.0 0.0 0.0 0.0 0.0 0.0]}
          hinv (ab/inverse-mass-matrix model qv)]
      (is (= 13 (ab/generalized-dof model)))
      (doseq [i [-1 0 2]
              dir [[1.0 0.0 0.0] (v/normalize [0.3 -0.8 0.5])]]
        (let [p (v/add (:pos (ab/frame-of model qv rs i)) [0.11 -0.07 0.23])
              fast (:delta-u (ab/impulse-at model qv rs i p dir))
              jt (lin/transpose (ab/point-jacobian model qv rs i p))
              slow (lin/mat-vec hinv (lin/mat-vec jt (v/normalize dir)))]
          (is (every? #(< (abs (double %)) 1e-9) (map - fast slow))
              (str "body " i " along " dir))))))

  (testing "and a floating chain of them still conserves both momenta"
    (let [model {:base {:mass 3.0 :com [0.1 0.0 0.0]
                        :inertia [[0.4 0.0 0.0] [0.0 0.6 0.0] [0.0 0.0 0.9]]}
                 :links [(limb -1 :spherical 0.7 1.2) (limb 0 :spherical 0.6 0.9)]}
          st0 {:q [(q/from-axis-angle [0.3 0.5 0.8] 0.4)
                   (q/from-axis-angle [0.9 0.1 0.2] -0.6)]
               :qd [1.1 -0.7 0.4 -0.3 0.9 1.5]
               :base {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6)
                      :pos [0.3 -0.2 0.5] :vel [0.4 -0.3 0.9 1.1 0.2 -0.6]}}
          free {:gravity [0.0 0.0 0.0]}
          m0 (ab/momentum model st0)
          m1 (ab/momentum model (reduce (fn [s _] (ab/step s model 1e-4 free))
                                        st0 (range 10000)))]
      (is (< (v/distance (:linear m0) (:linear m1)) 1e-2))
      (is (< (v/distance (:angular m0) (:angular m1)) 1e-2)))))

;; ---------------------------------------------------------------------------
;; Joint limits

(defn- swung
  "How far a ball joint's bone has swung from where it rests, in degrees."
  [link x]
  (let [rest-dir (v/normalize (:com link))
        bone (q/rotate (vec x) rest-dir)]
    (Math/toDegrees (Math/acos (max -1.0 (min 1.0 (v/dot bone rest-dir)))))))

(defn- twisted
  "How far it has turned about the bone itself, in degrees."
  [link x]
  (let [rest-dir (v/normalize (:com link))
        r (let [r (vec x)] (if (neg? (double (nth r 3))) (mapv - r) r))]
    (Math/toDegrees (* 2.0 (Math/atan2 (v/dot [(nth r 0) (nth r 1) (nth r 2)] rest-dir)
                                       (double (nth r 3)))))))

(defn- swinging-limb
  "A limb hanging off a ball joint, free to be kicked."
  [extra]
  [(merge (limb -1 :spherical 0.6 2.0) extra)])

(deftest hinge-limit-test
  (testing "a hinge stops at its limit instead of going round"
    (let [free [(limb -1 :revolute 0.6 2.0)]
          held [(assoc (limb -1 :revolute 0.6 2.0)
                       :limit [(Math/toRadians -40.0) (Math/toRadians 40.0)])]
          swing (fn [model]
                  (->> (iterate #(ab/step % model (/ 1.0 240.0) opts) {:q [0.0] :qd [12.0]})
                       (take 900)
                       (map #(Math/toDegrees (double (first (:q %)))))))]
      ;; Unlimited, it goes over the top and keeps going.
      (is (> (apply max (swing free)) 90.0))
      ;; Limited, it never leaves the range by more than the soft
      ;; constraint's overshoot.
      (is (< (apply max (swing held)) 42.0))
      (is (> (apply min (swing held)) -42.0)))))

(deftest cone-limit-test
  (testing "a ball joint's bone stays inside its cone"
    (let [kick {:q [[0.0 0.0 0.0 1.0]] :qd [0.0 0.0 6.0]}
          worst (fn [extra]
                  (let [model (swinging-limb extra)]
                    (->> (iterate #(ab/step % model (/ 1.0 240.0) opts) kick)
                         (take 1200)
                         (map #(swung (first model) (first (:q %))))
                         (apply max))))]
      ;; Free, the kick takes it most of the way over.
      (is (> (worst {}) 70.0))
      (is (< (worst {:cone (Math/toRadians 30.0)}) 32.0))
      (is (< (worst {:cone (Math/toRadians 60.0)}) 62.0))))

  (testing "and the cone says nothing about twist, which is its own limit"
    ;; The distinction is worth a test because it is the one that is easy
    ;; to miss by eye: a head can be within forty degrees of upright and
    ;; still be facing backward.
    ;;
    ;; Both cases keep the cone on. Swing and twist are read out of the
    ;; joint's quaternion by projecting onto the bone, and that
    ;; projection stops meaning anything as the swing approaches half a
    ;; turn -- so a twist measured on a limb that has swung right over
    ;; is noise, whatever the limit did.
    (let [spin {:q [[0.0 0.0 0.0 1.0]] :qd [5.0 0.0 0.0]}
          cone (Math/toRadians 30.0)
          worst (fn [extra]
                  (let [model (swinging-limb (assoc extra :cone cone))]
                    (->> (iterate #(ab/step % model (/ 1.0 240.0) opts) spin)
                         (take 600)
                         (map #(abs (twisted (first model) (first (:q %)))))
                         (apply max))))]
      (is (> (worst {}) 120.0) "a cone alone should leave twist free")
      (is (< (worst {:twist (Math/toRadians 45.0)}) 47.0)))))

(deftest limits-without-contacts-test
  (testing "limits hold with nothing to stand on"
    ;; They are solved in the same sweep as contacts, and for a while
    ;; that sweep only ran when there were contacts -- so a limb waving
    ;; in mid-air folded through itself freely and nothing said so.
    (let [model (swinging-limb {:cone (Math/toRadians 25.0)})]
      (is (ab/limited? model))
      (is (not (ab/limited? (swinging-limb {}))))
      (is (< (->> (iterate #(ab/step % model (/ 1.0 240.0) opts)
                           {:q [[0.0 0.0 0.0 1.0]] :qd [0.0 0.0 6.0]})
                  (take 600)
                  (map #(swung (first model) (first (:q %))))
                  (apply max))
             27.0)))))

;; ---------------------------------------------------------------------------
;; Self-collision

(defn- shaped-limb
  "A limb that can be collided with, hanging from `origin`."
  [parent origin]
  {:parent parent :joint :spherical :origin {:rot nil :pos origin}
   :mass 1.5 :com [0.0 -0.25 0.0]
   :inertia (let [f (/ 1.5 12.0)]
              [[(* f (+ (* 0.5 0.5) (* 0.1 0.1))) 0.0 0.0]
               [0.0 (* f (+ (* 0.1 0.1) (* 0.1 0.1))) 0.0]
               [0.0 0.0 (* f (+ (* 0.1 0.1) (* 0.5 0.5)))]])
   :shape :box :size [0.1 0.5 0.1]})

(defn- two-armed
  "A free root with two limbs hanging off it, `reach` apart."
  [reach]
  {:base {:mass 4.0 :com [0.0 0.0 0.0]
          :inertia [[0.05 0.0 0.0] [0.0 0.08 0.0] [0.0 0.0 0.08]]
          :shape :box :size [0.4 0.2 0.2]
          :shape-pose {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0]}}
   :links [(shaped-limb -1 [reach 0.0 0.0]) (shaped-limb -1 [(- reach) 0.0 0.0])]})

(def ^:private posed
  {:q [[0.0 0.0 0.0 1.0] [0.0 0.0 0.0 1.0]]
   :qd [0.0 0.0 0.0 0.0 0.0 0.0]
   :base {:rot (q/from-axis-angle [0.2 0.9 0.3] 0.6) :pos [0.3 -0.2 0.5]
          :vel [0.0 0.0 0.0 0.0 0.0 0.0]}})

(deftest internal-impulse-test
  (testing "a push between two of a model's own parts changes no momentum"
    ;; The sharpest check there is that the two halves of the push really
    ;; are equal and opposite. Nothing outside the model is touched, so
    ;; both momenta have to come out exactly as they went in -- while
    ;; something inside it does move.
    (let [model (two-armed 0.25)
          p [0.4 -0.3 0.5]
          dir (v/normalize [0.3 -0.8 0.5])
          m0 (ab/momentum model posed)
          after (ab/apply-pair-impulse model posed 0 1 p dir 3.0)
          m1 (ab/momentum model after)]
      (is (< (v/distance (:linear m0) (:linear m1)) 1e-12))
      (is (< (v/distance (:angular m0) (:angular m1)) 1e-12))
      (is (> (v/length (mapv - (:qd after) (:qd posed))) 1.0)
          "and it should actually have done something")))

  (testing "the mass it meets is a reduced mass"
    ;; Smaller than either body's alone, which sounds backward and is
    ;; not: both ends of this push give, where a contact with the floor
    ;; has one end that does not. With a root heavy enough to make the
    ;; two limbs independent it is the textbook identity exactly.
    (let [model (-> (two-armed 0.25)
                    (assoc-in [:base :mass] 100000.0)
                    (assoc-in [:base :inertia]
                              [[3000.0 0.0 0.0] [0.0 5000.0 0.0] [0.0 0.0 5000.0]]))
          p [0.4 -0.3 0.5]
          dir [1.0 0.0 0.0]
          mi (:effective-mass (ab/impulse-at model (:q posed) (:base posed) 0 p dir))
          mj (:effective-mass (ab/impulse-at model (:q posed) (:base posed) 1 p dir))
          mp (:effective-mass (ab/pair-impulse-at model (:q posed) (:base posed) 0 1 p dir))]
      (is (< mp (min mi mj)))
      (is (close? (/ 1.0 mp) (+ (/ 1.0 mi) (/ 1.0 mj)) 1e-4)))))

(deftest self-contact-pairs-test
  (testing "parts that are nowhere near each other make no contacts"
    (let [model (two-armed 0.25)]
      (is (empty? (ab/self-contacts model (:q posed) (:base posed)
                                    (ab/poses model (:q posed) (:base posed)))))))

  (testing "a link and its parent never collide, however much they overlap"
    ;; They meet at the joint and overlap there by construction. A
    ;; contact between them would be permanent and would push the model
    ;; apart from the inside.
    (let [chainy {:base (:base (two-armed 0.25))
                  :links [(shaped-limb -1 [0.0 0.0 0.0])
                          (shaped-limb 0 [0.0 -0.25 0.0])]}
          st {:q [[0.0 0.0 0.0 1.0] [0.0 0.0 0.0 1.0]] :qd (vec (repeat 6 0.0))
              :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0]
                     :vel (vec (repeat 6 0.0))}}]
      (is (empty? (ab/self-contacts chainy (:q st) (:base st)
                                    (ab/poses chainy (:q st) (:base st)))))))

  (testing "and `:no-collide` keeps a named pair apart too"
    (let [overlapping (two-armed 0.04)
          st {:q [[0.0 0.0 0.0 1.0] [0.0 0.0 0.0 1.0]] :qd (vec (repeat 6 0.0))
              :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 0.0 0.0]
                     :vel (vec (repeat 6 0.0))}}
          muted (assoc-in overlapping [:links 0 :no-collide] #{1})]
      (is (seq (ab/self-contacts overlapping (:q st) (:base st)
                                 (ab/poses overlapping (:q st) (:base st)))))
      (is (empty? (ab/self-contacts muted (:q st) (:base st)
                                    (ab/poses muted (:q st) (:base st))))))))

(deftest self-collision-test
  (testing "two limbs swung together stop instead of passing through"
    (let [model (two-armed 0.06)
          st0 {:q [(q/from-axis-angle [0.0 0.0 1.0] -0.9)
                   (q/from-axis-angle [0.0 0.0 1.0] 0.9)]
               :qd [0.0 0.0 -4.0 0.0 0.0 4.0]
               :base {:rot [0.0 0.0 0.0 1.0] :pos [0.0 2.0 0.0]
                      :vel (vec (repeat 6 0.0))}}
          tips-closest
          (fn [self?]
            (let [o {:gravity [0.0 0.0 0.0] :self-collide? self? :iterations 10}]
              (->> (iterate #(ab/step % model (/ 1.0 240.0) o) st0)
                   (take 200)
                   (map (fn [s]
                          (let [fs (ab/poses model (:q s) (:base s))
                                tip #(let [{:keys [rot pos]} (nth fs %)]
                                       (v/add pos (q/rotate rot [0.0 -0.5 0.0])))]
                            (v/distance (tip 0) (tip 1)))))
                   (apply min))))]
      ;; Without it they end up in the same place.
      (is (< (tips-closest false) 0.02))
      ;; With it they never get closer than their own width.
      (is (> (tips-closest true) 0.1)))))

(deftest flat-algebra-agrees-test
  (testing "the flat spatial algebra says what the readable one says"
    ;; The section above is the statement of what a Plucker transform and
    ;; a spatial inertia are, and the tests elsewhere check it. The flat
    ;; builders are what actually runs, for speed, and this is the only
    ;; thing tying the two together -- without it they are free to drift.
    (let [flat36 (fn [m] (vec (mapcat identity m)))
          near (fn [a b] (every? #(< (abs (double %)) 1e-12) (map - a b)))]
      (doseq [pose [{:rot nil :pos nil}
                    {:rot nil :pos [0.3 -0.7 1.1]}
                    {:rot (q/from-axis-angle [0.2 0.9 -0.3] 1.1) :pos [0.0 0.0 0.0]}
                    {:rot (q/from-axis-angle [0.7 -0.1 0.4] -2.3) :pos [-1.2 0.4 0.8]}]]
        (is (near (flat36 (ab/transform pose)) (seq (#'ab/transform-flat pose)))
            (str pose)))
      (doseq [[m com inertia]
              [[2.5 [0.0 0.0 0.0] (lin/eye 3)]
               [3.1 [0.4 -0.2 0.9] [[0.3 0.0 0.0] [0.0 0.5 0.0] [0.0 0.0 0.7]]]
               [0.7 [-1.1 0.6 0.2] [[0.2 0.03 -0.01] [0.03 0.4 0.05] [-0.01 0.05 0.6]]]]]
        (is (near (flat36 (ab/spatial-inertia m com inertia))
                  (seq (#'ab/spatial-inertia-flat m com inertia)))
            (str m " " com))))))
