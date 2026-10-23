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
    ;; one body; one contact along +y with the ground, through its center
    (let [jac (double-array 36)
          _ (doseq [[k d] [[0 [0.0 1.0 0.0]] [1 [1.0 0.0 0.0]] [2 [0.0 0.0 1.0]]]
                    j (range 3)]
              (aset jac (+ (* 12 k) 6 j) (double (nth d j))))
          problem {:nb 1 :mass (double-array [1.0])
                   :inertia (double-array [1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0])
                   :vstar (double-array [0.3 -2.0 0.0 0.0 0.0 0.0])
                   :nc 1 :oa (int-array [-1]) :ob (int-array [0]) :jac jac
                   :rn (double-array [1e-6]) :rt (double-array [1e-6])
                   :vhat (double-array [0.0]) :mu (double-array [0.5])}
          {:keys [v impulses iterations]} (cvx/solve problem)]
      (is (< (abs (aget ^doubles v 1)) 1e-4) "stopped")
      (is (< (abs (- (aget ^doubles impulses 0) 2.0)) 1e-3) "by an impulse of its momentum")
      (testing "and its slide is stopped too, the friction well inside the cone"
        (is (< (abs (aget ^doubles v 0)) 1e-4))
        (is (< (abs (- (aget ^doubles impulses 1) -0.3)) 1e-3)))
      (is (< iterations 10)))))
