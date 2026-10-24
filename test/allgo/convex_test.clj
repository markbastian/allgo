(ns allgo.convex-test
  "The convex contact solve on its own: the cone projection, its
  derivative, and the solve, on problems small enough to answer by hand."
  (:require [allgo.physics.convex :as cvx]
            [clojure.test :refer [deftest is testing]]))

(def ^:private project! @#'cvx/project!)

(defn- impulse [mu rn rt y]
  (let [g (double-array 3) gm (double-array 9)]
    (project! mu rn rt (double-array y) g gm)
    [(vec g) (vec gm)]))

(defn- close? [a b tol] (every? #(< (abs (double %)) tol) (map - a b)))

(deftest projection-test
  (testing "inside the cone the impulse is what the contact asks for"
    (is (close? [2.0 0.1 0.0] (first (impulse 0.5 1.0 1.0 [2.0 0.1 0.0])) 1e-12)))
  (testing "in the polar cone it is nothing"
    (is (close? [0.0 0.0 0.0] (first (impulse 0.5 1.0 1.0 [-1.0 0.2 0.0])) 1e-12)))
  (testing "between them it lands on the cone, the friction at its limit and opposing the slip"
    (let [[[gn g1 g2]] (impulse 0.5 1.0 1.0 [1.0 3.0 0.0])]
      (is (pos? gn))
      (is (< (abs (- (abs g1) (* 0.5 gn))) 1e-12))
      (is (pos? g1))
      (is (zero? g2))))
  (testing "G is the derivative of the impulse: -dg/dvc, checked by differences, in every region"
    ;; y = -R^-1 (vc - vhat), so dg/dvc = -dg/dy R^-1, and G = dg/dy R^-1
    (doseq [[mu rn rt y] [[0.5 0.3 0.01 [2.0 0.1 0.05]]
                          [0.5 0.3 0.01 [1.0 3.0 -2.0]]
                          [0.8 1.0 0.2 [0.4 -0.7 0.9]]]]
      (let [[_ gm] (impulse mu rn rt y)
            h 1e-7
            r [rn rt rt]]
        (doseq [j (range 3)]
          (let [bump (fn [s] (first (impulse mu rn rt (update y j + s))))
                col (mapv #(/ (- %1 %2) (* 2.0 h)) (bump h) (bump (- h)))]
            (doseq [i (range 3)]
              (is (< (abs (- (nth gm (+ (* 3 i) j)) (/ (nth col i) (nth r j)))) 1e-5)
                  (str [mu rn rt y] " G[" i "," j "]")))))))))

(deftest solve-test
  (testing "a unit mass falling onto a contact below it is stopped, and the impulse is what stopping it took"
    ;; one body, a block of six; one contact along +y with the ground, through its center
    (let [eye (vec (for [r (range 6) c (range 6)] (if (= r c) 1.0 0.0)))
          row (fn [d] [nil (into d [0.0 0.0 0.0])])
          [b i] (cvx/add-block (cvx/builder) 6 eye eye [0.3 -2.0 0.0 0.0 0.0 0.0])
          problem (cvx/build (cvx/add-constraint b cvx/contact -1 i
                                                 [(row [0.0 1.0 0.0]) (row [1.0 0.0 0.0]) (row [0.0 0.0 1.0])]
                                                 1e-6 1e-6 0.0 0.5))
          {:keys [v impulses iterations]} (cvx/solve problem)]
      (is (< (abs (aget ^doubles v 1)) 1e-4) "stopped")
      (is (< (abs (- (aget ^doubles impulses 0) 2.0)) 1e-3) "by an impulse of its momentum")
      (testing "and its slide is stopped too, the friction well inside the cone"
        (is (< (abs (aget ^doubles v 0)) 1e-4))
        (is (< (abs (- (aget ^doubles impulses 1) -0.3)) 1e-3)))
      (is (< iterations 10))))
  (testing "a pin holds two unit masses together, either way; a limit pushes and never pulls"
    (let [eye (vec (for [r (range 6) c (range 6)] (if (= r c) 1.0 0.0)))
          [b0 p] (cvx/add-block (cvx/builder) 6 eye eye [1.0 0.0 0.0 0.0 0.0 0.0])
          [b0 q] (cvx/add-block b0 6 eye eye [-1.0 0.0 0.0 0.0 0.0 0.0])
          x [1.0 0.0 0.0 0.0 0.0 0.0] -x (mapv - x)
          pinned (cvx/build (cvx/add-constraint b0 cvx/bilateral p q [[x -x]] 1e-6 0.0 0.0 0.0))
          {v :v} (cvx/solve pinned)]
      (is (< (abs (- (aget ^doubles v 0) (aget ^doubles v 6))) 1e-4) "moving together")
      (is (< (abs (+ (aget ^doubles v 0) (aget ^doubles v 6))) 1e-6) "with the momentum they had")
      (let [limited (fn [vx] (let [[b i] (cvx/add-block (cvx/builder) 6 eye eye [vx 0.0 0.0 0.0 0.0 0.0])]
                               (aget ^doubles (:v (cvx/solve (cvx/build (cvx/add-constraint b cvx/unilateral i -1 [[x nil]]
                                                                                            1e-6 0.0 0.0 0.0)))) 0)))]
        (is (< (abs (limited -1.0)) 1e-4) "pushed back to the limit")
        (is (= 1.0 (limited 1.0)) "and left alone moving away from it"))))
  (testing "islands are solved apart, and the answer is each one's alone"
    ;; a body falling onto a contact, a pinned pair, and one nothing
    ;; touches -- together, and each by itself
    (let [eye (vec (for [r (range 6) c (range 6)] (if (= r c) 1.0 0.0)))
          heavy (mapv #(* 3.0 %) eye) light (mapv #(/ % 3.0) eye)
          x [1.0 0.0 0.0 0.0 0.0 0.0] -x (mapv - x)
          row (fn [d] [nil (into d [0.0 0.0 0.0])])
          falling (fn [b] (let [[b i] (cvx/add-block b 6 eye eye [0.3 -2.0 0.0 0.0 0.0 0.0])]
                            (cvx/add-constraint b cvx/contact -1 i
                                                [(row [0.0 1.0 0.0]) (row [1.0 0.0 0.0]) (row [0.0 0.0 1.0])]
                                                1e-6 1e-6 0.0 0.5)))
          pinned (fn [b] (let [[b p] (cvx/add-block b 6 heavy light [1.0 0.5 0.0 0.0 0.0 0.0])
                               [b q] (cvx/add-block b 6 eye eye [-1.0 0.0 0.0 0.0 0.2 0.0])]
                           (cvx/add-constraint b cvx/bilateral p q [[x -x]] 1e-6 0.0 0.0 0.0)))
          alone (fn [b] (first (cvx/add-block b 6 eye eye [0.0 -1.0 0.0 0.0 0.0 0.0])))
          whole (cvx/solve (cvx/build (-> (cvx/builder) falling alone pinned)))
          a (cvx/solve (cvx/build (falling (cvx/builder))))
          b (cvx/solve (cvx/build (pinned (cvx/builder))))]
      (is (close? (vec (:v whole)) (concat (:v a) [0.0 -1.0 0.0 0.0 0.0 0.0] (:v b)) 1e-12))
      (is (close? (vec (:impulses whole)) (concat (:impulses a) (:impulses b)) 1e-12))
      (is (= (:iterations whole) (max (:iterations a) (:iterations b))))))
  (testing "a factor kept from another step preconditions the Newton steps, and the answer is the same"
    ;; the same blocks and constraints each time, going different
    ;; speeds, and once weighing ten times as much
    (let [eye (vec (for [r (range 6) c (range 6)] (if (= r c) 1.0 0.0)))
          x [1.0 0.0 0.0 0.0 0.0 0.0] -x (mapv - x)
          row (fn [d] [nil (into d [0.0 0.0 0.0])])
          problem (fn [m speed]
                    (let [mass (mapv #(* m %) eye) minv (mapv #(/ % m) eye)
                          [b i] (cvx/add-block (cvx/builder) 6 mass minv [0.3 (- speed) 0.0 0.0 0.0 0.0])
                          b (cvx/add-constraint b cvx/contact -1 i
                                                [(row [0.0 1.0 0.0]) (row [1.0 0.0 0.0]) (row [0.0 0.0 1.0])]
                                                1e-6 1e-6 0.0 0.5)
                          [b p] (cvx/add-block b 6 mass minv [speed 0.5 0.0 0.0 0.0 0.0])
                          [b q] (cvx/add-block b 6 eye eye [-1.0 0.0 0.0 0.0 0.2 0.0])]
                      (cvx/build (cvx/add-constraint b cvx/bilateral p q [[x -x]] 1e-6 0.0 0.0 0.0))))
          kept (:factors (cvx/solve (problem 1.0 2.0)))
          heavy (:factors (cvx/solve (problem 10.0 2.0)))]
      (is (= 2 (count kept)) "one an island")
      (doseq [speed [2.0 2.5 4.0]]
        (let [fresh (cvx/solve (problem 1.0 speed))]
          (doseq [factors [kept heavy]]
            (is (close? (vec (:v fresh)) (vec (:v (cvx/solve (problem 1.0 speed) nil {:factors factors})))
                        1e-9))))))))

(deftest fallback-test
  (testing "a Newton solve that runs out of iterations gives way to projected Gauss-Seidel, and says so"
    ;; Stopped short, a Newton step is lower in cost than where it began
    ;; and can be nothing like the minimum; the fallback is slower to it
    ;; but bounded every step. Given no iterations at all, it is what
    ;; answers, and on problems this small it answers as Newton does.
    (let [eye (vec (for [r (range 6) c (range 6)] (if (= r c) 1.0 0.0)))
          row (fn [d] [nil (into d [0.0 0.0 0.0])])
          [b i] (cvx/add-block (cvx/builder) 6 eye eye [0.3 -2.0 0.0 0.0 0.0 0.0])
          falling (cvx/build (cvx/add-constraint b cvx/contact -1 i
                                                 [(row [0.0 1.0 0.0]) (row [1.0 0.0 0.0]) (row [0.0 0.0 1.0])]
                                                 1e-6 1e-6 0.0 0.5))
          x [1.0 0.0 0.0 0.0 0.0 0.0] -x (mapv - x)
          [b0 p] (cvx/add-block (cvx/builder) 6 eye eye [1.0 0.0 0.0 0.0 0.0 0.0])
          [b0 q] (cvx/add-block b0 6 eye eye [-1.0 0.0 0.0 0.0 0.0 0.0])
          pinned (cvx/build (cvx/add-constraint b0 cvx/bilateral p q [[x -x]] 1e-6 0.0 0.0 0.0))]
      (doseq [problem [falling pinned]]
        (let [newton (cvx/solve problem nil {})
              fallback (cvx/solve problem nil {:max-iterations 0})]
          (is (nil? (:fallback newton)))
          (is (= :not-converged (:fallback fallback)))
          (is (close? (vec (:v newton)) (vec (:v fallback)) 1e-4))
          (is (close? (vec (:impulses newton)) (vec (:impulses fallback)) 1e-3)))))))
