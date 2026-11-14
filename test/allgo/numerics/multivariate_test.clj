(ns allgo.numerics.multivariate-test
  "Polynomials in several variables against their expansions by hand and
  their derivatives against central differences."
  (:require [allgo.numerics.multivariate :as mv]
            [clojure.test :refer [deftest is testing]]))

(def ^:private x (mv/variable 3 0))
(def ^:private y (mv/variable 3 1))
(def ^:private z (mv/variable 3 2))

(defn- torus [R r]
  (let [s (mv/add (mv/power x 2) (mv/power y 2) (mv/power z 2) (mv/constant 3 (- (* R R) (* r r))))]
    (mv/sub (mv/power s 2) (mv/scale (mv/add (mv/power x 2) (mv/power y 2)) (* 4 R R)))))

(deftest arithmetic-test
  (testing "(x + y)^2 = x^2 + 2xy + y^2, and (x + y)(x - y) = x^2 - y^2: the xy terms cancel and are dropped"
    (is (= {[2 0 0] 1.0 [1 1 0] 2.0 [0 2 0] 1.0} (mv/power (mv/add x y) 2)))
    (is (= {[2 0 0] 1.0 [0 2 0] -1.0} (mv/mul (mv/add x y) (mv/sub x y)))))
  (testing "anything to the power 0 is 1"
    (is (= {[0 0 0] 1.0} (mv/power z 0)))))

(deftest evaluation-test
  (let [p (torus 1.0 0.4)]
    (testing "the torus is zero on its outer equator and its top circle"
      (is (< (abs (mv/evaluate p [1.4 0.0 0.0])) 1e-12))
      (is (< (abs (mv/evaluate p [0.0 1.0 0.4])) 1e-12)))
    (testing "the evaluator agrees with evaluate"
      (doseq [q [[0.3 -0.7 0.2] [1.1 0.5 -0.9] [0.0 0.0 0.0]]]
        (is (< (abs (- (mv/evaluate p q) ((mv/evaluator p) q))) 1e-12))))))

(deftest derivative-test
  (let [p (torus 1.0 0.4)
        q [0.3 -0.7 0.2]
        h 1e-5
        nudge (fn [v i d] (update v i + d))
        central (fn [f i v] (/ (- (f (nudge v i h)) (f (nudge v i (- h)))) (* 2 h)))]
    (testing "d/dx of x^3 y is 3 x^2 y"
      (is (= {[2 1 0] 3.0} (mv/derivative (mv/mul (mv/power x 3) y) 0))))
    (testing "the gradient and Hessian against central differences"
      (doseq [i (range 3)]
        (is (< (abs (- (mv/evaluate ((mv/gradient p 3) i) q) (central #(mv/evaluate p %) i q))) 1e-7))
        (doseq [j (range 3)]
          (is (< (abs (- (mv/evaluate (get-in (mv/hessian p 3) [i j]) q)
                         (central #(mv/evaluate ((mv/gradient p 3) j) %) i q)))
                 1e-6)))))))
