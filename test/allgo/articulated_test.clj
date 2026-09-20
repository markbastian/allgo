(ns allgo.articulated-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.numerics.linear :as lin]
            [allgo.physics.articulated :as ab]
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
    ;; qdd = -(m g d / I_pivot) cos q, with d the distance to the centre
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

  (testing "a fixed joint does not move"
    (let [model [(assoc (rod -1 1.0 1.0) :joint :fixed)]]
      (is (= [0.0] (:qdd (ab/forward-dynamics model [0.0] [0.0] [5.0] opts)))))))

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
  "A free body whose centre of mass is at its frame origin and whose
  inertia is the same about every axis, so the textbook formulae for
  being hit apply without any further care."
  [m i]
  {:base {:mass m :com [0.0 0.0 0.0]
          :inertia [[i 0.0 0.0] [0.0 i 0.0] [0.0 0.0 i]]}
   :links []})

(deftest impulse-test
  (testing "through the centre of mass, a free body weighs what it weighs"
    (let [m 2.5]
      (is (close? (:effective-mass (ab/impulse-at (sphere-ish m 0.7) [] at-rest -1
                                                  [0.0 0.0 0.0] [1.0 0.0 0.0]))
                  m 1e-12))))

  (testing "off the centre of mass it weighs less, by the textbook amount"
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
