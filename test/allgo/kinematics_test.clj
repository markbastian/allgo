(ns allgo.kinematics-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.kinematics.analytic :as ik]
            [allgo.kinematics.chain :as k]
            [allgo.kinematics.numeric :as num]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private half-pi (/ math/PI 2))

(def puma
  "The PUMA 560, in Craig's modified DH parameters and metres. The
  canonical six-axis arm, and the one his chapter 4 solves."
  (k/chain [{:alpha 0.0         :a 0.0    :d 0.0    :limits [-2.8 2.8]}
            {:alpha (- half-pi) :a 0.0    :d 0.0    :limits [-3.9 0.8]}
            {:alpha 0.0         :a 0.4318 :d 0.1245 :limits [-0.8 3.9]}
            {:alpha (- half-pi) :a 0.0203 :d 0.4318 :limits [-3.5 3.5]}
            {:alpha half-pi     :a 0.0    :d 0.0    :limits [-1.8 1.8]}
            {:alpha (- half-pi) :a 0.0    :d 0.0    :limits [-5.2 5.2]}]))

(defn- rng-angles [^java.util.Random r n]
  (vec (repeatedly n #(* math/PI (- (.nextDouble r) 0.5) 2))))

(defn- pose-error [a b]
  (max (v/distance (:pos a) (:pos b))
       ;; A quaternion and its negation are the same orientation.
       (min (q/angle (q/between (:rot a) (:rot b)))
            (q/angle (q/between (:rot a) (mapv - (:rot b)))))))

(deftest dh-transform-test
  (testing "the quaternion form is Craig's matrix"
    ;; The transform is built by composing rotations rather than writing
    ;; the matrix out, so it is worth pinning against the matrix itself.
    (let [r (java.util.Random. 5)]
      (dotimes [_ 200]
        (let [alpha (* math/PI (- (.nextDouble r) 0.5) 2)
              a (* 3 (.nextDouble r))
              d (* 3 (.nextDouble r))
              theta (* math/PI (- (.nextDouble r) 0.5) 2)
              ct (math/cos theta) st (math/sin theta)
              ca (math/cos alpha) sa (math/sin alpha)
              expected-rot [[ct (- st) 0.0]
                            [(* st ca) (* ct ca) (- sa)]
                            [(* st sa) (* ct sa) ca]]
              expected-pos [a (* (- sa) d) (* ca d)]
              got (k/link-transform (k/link {:alpha alpha :a a :d d}) theta)]
          (is (< (reduce max (map (fn [r1 r2] (reduce max (map (comp abs -) r1 r2)))
                                  expected-rot (q/to-matrix (:rot got))))
                 1e-12))
          (is (< (v/distance expected-pos (:pos got)) 1e-12)))))))

(deftest forward-kinematics-test
  (testing "the PUMA at rest reaches where its table says it should"
    (let [[x y z] (:pos (k/pose puma (k/home puma)))]
      (is (< (abs (- x (+ 0.4318 0.0203))) 1e-12) "out by the two link lengths")
      (is (< (abs (- y 0.1245)) 1e-12) "across by the shoulder offset")
      (is (< (abs (+ z 0.4318)) 1e-12) "down by the forearm")))

  (testing "there is a frame per joint, plus the base"
    (is (= 7 (count (k/frames puma (k/home puma))))))

  (testing "a tool offset moves the tool and not the wrist"
    (let [with-tool (k/chain (:links puma) {:rot q/identity-q :pos [0.0 0.0 0.1]})
          bare (:pos (k/pose puma (k/home puma)))
          tooled (:pos (k/pose with-tool (k/home with-tool)))]
      (is (< (abs (- 0.1 (v/distance bare tooled))) 1e-12))))

  (testing "composing a transform with its inverse is the identity"
    (let [t (k/link-transform (k/link {:alpha 0.7 :a 1.3 :d 0.4}) 1.1)
          i (k/compose t (k/invert t))]
      (is (< (v/length (:pos i)) 1e-12))
      (is (< (q/angle (:rot i)) 1e-12)))))

(deftest jacobian-test
  (testing "the Jacobian is the derivative of forward kinematics"
    ;; Checked against finite differences rather than against another
    ;; derivation, so an error in either one shows up.
    (let [r (java.util.Random. 11)
          eps 1e-6]
      (dotimes [_ 40]
        (let [qs (rng-angles r 6)
              jac (k/jacobian puma qs)]
          (dotimes [i 6]
            (let [pp (k/pose puma (update qs i + eps))
                  pm (k/pose puma (update qs i - eps))
                  dv (v/scale (v/sub (:pos pp) (:pos pm)) (/ 1.0 (* 2 eps)))
                  dw (v/scale (vec (take 3 (q/mul (:rot pp) (q/inverse (:rot pm)))))
                              (/ 1.0 eps))
                  col (nth jac i)]
              (is (< (v/distance dv (vec (take 3 col))) 1e-6))
              (is (< (v/distance dw (vec (drop 3 col))) 1e-6))))))))

  (testing "a prismatic joint slides and does not turn"
    (let [c (k/chain [{:alpha 0.0 :a 0.0 :d 0.0 :joint :prismatic}])
          col (first (k/jacobian c [0.5]))]
      (is (< (v/length (vec (drop 3 col))) 1e-12) "no angular part")
      (is (< (abs (- 1.0 (v/length (vec (take 3 col))))) 1e-12)))))

(deftest analytic-solutions-test
  (testing "the PUMA has a spherical wrist and a chain of five does not"
    (is (ik/spherical-wrist? puma))
    (is (not (ik/spherical-wrist? (k/chain (take 5 (:links puma))))))
    (is (thrown? clojure.lang.ExceptionInfo
                 (ik/solutions (k/chain (take 5 (:links puma)))
                               {:rot q/identity-q :pos [0.5 0.0 0.0]}))))

  (testing "every reachable pose has exactly eight configurations"
    ;; Two shoulder, two elbow, two wrist. The multiplicity is the
    ;; structure, not an artefact.
    (let [r (java.util.Random. 17)]
      (dotimes [_ 300]
        (let [qs (rng-angles r 6)
              target (k/pose puma qs)]
          (is (= 8 (count (ik/solutions puma target))))))))

  (testing "and every one of them reaches the pose"
    (let [r (java.util.Random. 19)]
      (dotimes [_ 200]
        (let [qs (rng-angles r 6)
              target (k/pose puma qs)]
          (doseq [s (ik/solutions puma target)]
            (is (< (pose-error (k/pose puma s) target) 1e-9)))))))

  (testing "including the configuration it was generated from"
    (let [r (java.util.Random. 23)]
      (dotimes [_ 200]
        (let [qs (rng-angles r 6)
              target (k/pose puma qs)]
          (is (some (fn [s]
                      (< (reduce max (map (fn [a b] (abs (k/wrap-angle (- a b)))) s qs))
                         1e-6))
                    (ik/solutions puma target))
              "the arm cannot reach its own pose")))))

  (testing "a target out of reach has no solutions at all"
    ;; Not an error and not an approximation: the honest answer is none.
    (is (empty? (ik/solutions puma {:rot q/identity-q :pos [2.0 0.0 0.0]})))
    (is (empty? (ik/solutions puma {:rot q/identity-q :pos [0.0 0.0 5.0]})))
    (is (seq (ik/solutions puma {:rot q/identity-q :pos [0.5 0.0 0.0]}))))

  (testing "joint limits cut the eight down"
    (let [r (java.util.Random. 29)
          rows (for [_ (range 200)
                     :let [qs (rng-angles r 6)]
                     :when (k/within-limits? puma qs)
                     :let [target (k/pose puma qs)]]
                 [(count (ik/solutions puma target))
                  (count (ik/solutions puma target {:limits? true}))])]
      (is (seq rows))
      (is (every? (fn [[all limited]] (and (= 8 all) (<= limited all))) rows))
      (is (some (fn [[_ limited]] (< limited 8)) rows)
          "limits should rule some out")
      (is (every? (fn [[_ limited]] (pos? limited)) rows)
          "the configuration it came from is always one of them")))

  (testing "a limited solution really is within the limits"
    (let [r (java.util.Random. 31)]
      (dotimes [_ 100]
        (let [qs (rng-angles r 6)]
          (when (k/within-limits? puma qs)
            (doseq [s (ik/solutions puma (k/pose puma qs) {:limits? true})]
              (is (k/within-limits? puma s))))))))

  (testing "a wrist singularity collapses solutions rather than failing"
    ;; Joint 5 straight: joints 4 and 6 turn about the same line, so only
    ;; their sum is determined and a whole family of configurations works.
    (let [target (k/pose puma [0.3 -0.6 0.9 0.4 0.0 0.2])
          sols (ik/solutions puma target)]
      (is (seq sols))
      (is (< (count sols) 8))
      (doseq [s sols]
        (is (< (pose-error (k/pose puma s) target) 1e-9))))))

(deftest nearest-test
  (testing "it picks the configuration that moves the joints least"
    (let [target (k/pose puma [0.3 -0.6 0.9 0.4 0.7 0.2])
          from (k/home puma)
          chosen (ik/nearest puma target from)
          travel (fn [s] (reduce + (map (fn [a b] (abs (k/wrap-angle (- a b)))) s from)))]
      (is (some? chosen))
      (is (= (travel chosen) (reduce min (map travel (ik/solutions puma target)))))
      (is (< (pose-error (k/pose puma chosen) target) 1e-9))))

  (testing "and gives nothing when nothing reaches"
    (is (nil? (ik/nearest puma {:rot q/identity-q :pos [3.0 0.0 0.0]} (k/home puma))))))

(deftest limits-test
  (testing "clamping brings a configuration into range"
    (let [c (k/chain [{:limits [-1.0 1.0]} {:limits [0.0 2.0]} {}])]
      (is (= [1.0 0.0 5.0] (k/clamp c [3.0 -1.0 5.0])))
      (is (k/within-limits? c (k/clamp c [3.0 -1.0 5.0])))
      (is (not (k/within-limits? c [3.0 0.5 0.0])))))

  (testing "angles wrap to the half turn either side of zero"
    (is (< (abs (- 0.5 (k/wrap-angle (+ 0.5 (* 4 math/PI))))) 1e-12))
    (is (< (abs (- -0.5 (k/wrap-angle (- -0.5 (* 4 math/PI))))) 1e-12))
    (is (<= (- math/PI) (k/wrap-angle 100.0) math/PI))
    (testing "and a wrapped angle is the same rotation it started as"
      (let [r (java.util.Random. 37)]
        (dotimes [_ 200]
          (let [a (* 40 (- (.nextDouble r) 0.5))]
            (is (< (q/angle (q/between (q/from-axis-angle [0 0 1] a)
                                       (q/from-axis-angle [0 0 1] (k/wrap-angle a))))
                   1e-9))))))))

(def redundant
  "Seven joints and no spherical wrist, so no closed form exists."
  (k/chain [{:alpha 0.0         :a 0.0  :d 0.3}
            {:alpha (- half-pi) :a 0.0  :d 0.0}
            {:alpha half-pi     :a 0.0  :d 0.4}
            {:alpha (- half-pi) :a 0.05 :d 0.0}
            {:alpha half-pi     :a 0.0  :d 0.4}
            {:alpha (- half-pi) :a 0.0  :d 0.0}
            {:alpha half-pi     :a 0.0  :d 0.1}]))

(deftest numeric-solver-test
  (testing "pose error is zero exactly when the poses match"
    (let [p (k/pose puma [0.2 0.3 -0.4 0.5 0.6 0.7])]
      (is (< (v/length (vec (take 3 (num/pose-error p p)))) 1e-12))
      (is (< (v/length (vec (drop 3 (num/pose-error p p)))) 1e-12))))

  (testing "it reaches poses on a chain that has a closed form"
    (let [r (java.util.Random. 41)]
      (dotimes [_ 60]
        (let [qs (rng-angles r 6)
              target (k/pose puma qs)
              {:keys [values converged?]} (num/solve puma target (k/home puma)
                                                     {:iterations 300 :limits? false})]
          (is converged?)
          (is (< (pose-error (k/pose puma values) target) 1e-5))))))

  (testing "and on one that does not"
    ;; Seven joints, so the solution is not even unique. Descent does not
    ;; care; it walks to whichever one it finds.
    (is (not (ik/spherical-wrist? redundant)))
    (let [r (java.util.Random. 43)]
      (dotimes [_ 60]
        (let [qs (rng-angles r 7)
              target (k/pose redundant qs)
              {:keys [values converged?]} (num/solve redundant target
                                                     (k/home redundant)
                                                     {:iterations 300 :limits? false})]
          (is converged?)
          (is (< (pose-error (k/pose redundant values) target) 1e-5))))))

  (testing "an unreachable target does not converge, and says so"
    ;; And still reports where it got to, which is the nearest approach.
    (let [{:keys [converged? values error]}
          (num/solve puma {:rot q/identity-q :pos [5.0 0.0 0.0]} (k/home puma))]
      (is (not converged?))
      (is (pos? error))
      (is (every? #(Double/isFinite (double %)) values) "and has not flailed")))

  (testing "it respects joint limits when asked"
    ;; Aimed at poses the limits actually allow, so that respecting them
    ;; and reaching the target are not in conflict.
    (let [r (java.util.Random. 47)
          reachable (->> (repeatedly #(rng-angles r 6))
                         (filter #(k/within-limits? puma %))
                         (take 20))]
      (is (seq reachable))
      (doseq [qs reachable]
        (let [{:keys [values converged?]} (num/solve puma (k/pose puma qs) (k/home puma)
                                                     {:limits? true :iterations 300})]
          (is (k/within-limits? puma values) "stayed inside the limits")
          (when converged?
            (is (< (pose-error (k/pose puma values) (k/pose puma qs)) 1e-5)))))))

  (testing "it agrees with the exact solver where both apply"
    (let [r (java.util.Random. 53)]
      (dotimes [_ 30]
        (let [qs (rng-angles r 6)
              target (k/pose puma qs)
              {:keys [values converged?]} (num/solve puma target (k/home puma)
                                                     {:iterations 300 :limits? false})]
          (when converged?
            ;; Not the same joint angles -- it found whichever configuration
            ;; it walked to -- but the same place.
            (is (< (pose-error (k/pose puma values) target) 1e-5))
            (is (seq (ik/solutions puma target))))))))

  (testing "descent cannot tell you how many configurations there are"
    ;; The reason the closed form is worth having. Several seeds turn up
    ;; several answers and never say whether that is all of them.
    (let [target (k/pose puma [0.3 -0.6 0.9 0.4 0.7 0.2])
          exact (ik/solutions puma target)
          found (num/solve-from-many puma target 8 {:iterations 300 :limits? false})]
      (is (= 8 (count exact)))
      (is (pos? (count found)))
      (is (<= (count found) (count exact))
          "descent can only ever find configurations that exist")
      (doseq [s found]
        (is (< (pose-error (k/pose puma s) target) 1e-5))))))
